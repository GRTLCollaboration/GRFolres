/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 *
 * Port of GRChombo/Tests/FourDerivScalarTensorTest.
 *
 * GRChombo's version fed hand-built d1/d2/advec structs into the RHS functions
 * and compared against Mathematica. GRTeclyn computes every derivative from the
 * grid, so instead this test evaluates the full modified-CCZ4 + 4dST right hand
 * side on a fixed polynomial grid (the standard random_ccz4_initial_data plus a
 * polynomial scalar field) and compares one interior cell against reference
 * numbers in values1.hpp.
 *
 * The reference numbers are produced by
 *   GRChombo/Tests/FourDerivScalarTensorGridTest
 * which runs GRChombo's ModifiedCCZ4RHS::compute on the identical grid with the
 * identical parameters and constant coupling. Regenerate values1.hpp from that
 * program's stdout whenever the 4dST equations change.
 */

// Doctest header
#include "doctest.h"

// Common test headers
#include "InitialData.hpp"
#include "doctestCLIArgs.hpp"

// Test headers
#include "FourDerivScalarTensorTest.hpp"
#include "TestCouplingAndPotential.hpp"

// GRTeclyn headers
#include "CCZ4RHS.hpp"
#include "FourDerivScalarTensor.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRParmParse.hpp"
#include "ModifiedCCZ4RHS.hpp"
#include "ModifiedPunctureGauge.hpp"
#include "StateVariables.hpp"

// AMReX headers
#include "AMReX.H"
#include "AMReX_FArrayBox.H"

// System headers
#include <array>
#include <cmath>
#include <iomanip>

namespace
{
using TheoryType =
    FourDerivScalarTensor<TestCouplingAndPotential, FourthOrderDerivatives>;
using ModifiedRHSType = ModifiedCCZ4RHS<TheoryType, FourthOrderDerivatives>;
using GaugeType       = ModifiedPunctureGauge<FourthOrderDerivatives>;

void set_test_parameters()
{
    GRParmParse pp;

    // CCZ4 evolution
    pp.add("ccz4.formulation", CCZ4RHS<>::USE_CCZ4);
    pp.add("ccz4.kappa1", 0.1);
    pp.add("ccz4.kappa2", 0.0);
    pp.add("ccz4.kappa3", 1.0);
    pp.add("ccz4.covariantZ4", false);
    pp.add("evolution.sigma", 1.0);

    // Moving-puncture gauge
    pp.add("gauge.lapse_advec_coeff", 1.0);
    pp.add("gauge.lapse_power", 1.0);
    pp.add("gauge.lapse_coeff", 2.0);
    pp.add("gauge.shift_advec_coeff", 1.0);
    pp.add("gauge.shift_Gamma_coeff", 0.75);
    pp.add("gauge.eta", 1.0);

    // Modified gauge (constant a(x), b(x))
    pp.add("mod_gauge.mod_a", 0.35);
    pp.add("mod_gauge.mod_b", 0.55);

    // 4dST theory
    pp.add("four_deriv_scalar_tensor.G_Newton", 1.0);
}
} // namespace

void run_four_deriv_scalar_tensor_test()
{
    int amrex_argc    = doctest::cli_args.argc();
    char **amrex_argv = doctest::cli_args.argv();
    // NOLINTNEXTLINE(bugprone-casting-through-void) // Open MPI triggers this
    amrex::Initialize(amrex_argc, amrex_argv);
    {
        constexpr int num_cells  = 16;
        constexpr int num_ghosts = 3;
        constexpr amrex::Real dx = 0.5 / num_cells;

        // The single interior cell we compare against the reference
        const amrex::IntVect probe(num_cells / 2, num_cells / 2,
                                   num_cells / 2);

        amrex::Box box(amrex::IntVect(0, 0, 0),
                       amrex::IntVect(num_cells - 1, num_cells - 1,
                                      num_cells - 1));
        amrex::Box ghosted_box = box;
        ghosted_box.grow(num_ghosts);

        amrex::FArrayBox in_fab{ghosted_box, NUM_VARS,
                                amrex::The_Managed_Arena()};
        amrex::FArrayBox out_fab{box, NUM_VARS, amrex::The_Managed_Arena()};
        out_fab.setVal(0.0);

        const auto &in_array   = in_fab.array();
        const auto &in_c_array = in_fab.const_array();
        const auto &out_array  = out_fab.array();

        // Polynomial initial data (identical on the GRChombo side)
        amrex::ParallelFor(
            ghosted_box,
            [=] AMREX_GPU_DEVICE(int ix, int iy, int iz)
            {
                const amrex::IntVect iv{ix, iy, iz};
                const amrex::RealVect coords = amrex::RealVect{iv} * dx;
                random_ccz4_initial_data(iv, in_array, coords);
                fdst_scalar_initial_data(iv, in_array, coords);
            });
        amrex::Gpu::streamSynchronize();

        set_test_parameters();

        const ModifiedRHSType modified_ccz4(dx);
        const GaugeType modified_puncture_gauge(dx);

        // 1. chi and h_ij
        amrex::ParallelFor(box,
                           [=] AMREX_GPU_DEVICE(int ix, int iy, int iz)
                           {
                               modified_ccz4.compute_chi_and_h_ij(
                                   ix, iy, iz, out_array, in_c_array);
                           });

        // 2. A_ij, Theta and Gamma
        amrex::ParallelFor(box,
                           [=] AMREX_GPU_DEVICE(int ix, int iy, int iz)
                           {
                               modified_ccz4.compute_A_ij_and_Theta_and_Gamma(
                                   ix, iy, iz, out_array, in_c_array);
                           });

        // 3. modified gauge, b(x) gauge terms, effective EM tensor, scalar
        //    field evolution, principal-part solve and dissipation
        amrex::ParallelFor(
            box,
            [=] AMREX_GPU_DEVICE(int ix, int iy, int iz)
            {
                modified_puncture_gauge.calculate_rhs(ix, iy, iz, out_array,
                                                      in_c_array);
                modified_ccz4.add_b_rhs(ix, iy, iz, out_array, in_c_array);
                modified_ccz4.add_emtensor_rhs(ix, iy, iz, out_array,
                                               in_c_array);
                modified_ccz4.add_theory_rhs(ix, iy, iz, out_array, in_c_array);
                modified_ccz4.solve_lhs(ix, iy, iz, out_array, in_c_array);
                modified_ccz4.apply_dissipation(ix, iy, iz, out_array,
                                                in_c_array);
            });

        amrex::Gpu::streamSynchronize();

        // Reference values for the probe cell, indexed by component. See the
        // file header for how to regenerate this.
        std::array<amrex::Real, NUM_VARS> known{};
        bool have_reference = true;
#include "values1.hpp"

        constexpr amrex::Real tol    = 1.0e-9;
        constexpr int cout_precision = 17;

        if (!have_reference)
        {
            MESSAGE("values1.hpp has no reference numbers yet - run "
                    "GRChombo/Tests/FourDerivScalarTensorGridTest and paste "
                    "its output into values1.hpp");
        }

        for (int comp = 0; comp < NUM_VARS; ++comp)
        {
            const amrex::Real computed = out_fab.array()(probe, comp);
            const amrex::Real reference = known[comp];
            const amrex::Real diff      = std::abs(computed - reference);

            INFO("component " << StateVariables::names[comp] << " (" << comp
                              << "): computed "
                              << std::setprecision(cout_precision) << computed
                              << ", reference " << reference << ", diff "
                              << diff);
            if (have_reference)
            {
                CHECK(diff <= tol);
            }
        }

        CHECK(!out_fab.contains_nan(box, 0, NUM_VARS));
    }
    amrex::Finalize();
}
