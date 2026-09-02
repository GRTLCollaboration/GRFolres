/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "BinaryBH4dSTLevel.hpp"

// GRTeclyn infrastructure
#include "AlgebraicConstraintsEnforcer.hpp"
#include "BinaryBHInitialData.hpp" // also pulls in BoostedBHInitialData.hpp
#include "ChiTagger.hpp"
#include "ExtractionTagger.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRParmParse.hpp"
#include "PositiveChiAndLapse.hpp"
#include "PunctureTagger.hpp"
#include "PunctureTracker.hpp"
#include "SixthOrderDerivatives.hpp"
#include "SphericalExtractionParameters.hpp"
#include "WeylExtraction.hpp"
#ifdef USE_TWOPUNCTURES
#include "TwoPuncturesInitialData.hpp"
#endif

// Scalar field
#include "CouplingAndPotential.hpp"
#include "InitialScalarData.hpp"

// 4dST / modified-gravity source classes.
// NOTE: the include files below do not exist yet and Llibert will port them in
// Source and presumably they will be in
// GRTeclyn/Source/{FourDerivScalarTensor,ModifiedCCZ4}/
// For this example I assume the following: the same class names as the GRFolres
// ones, but expressed in GRTeclyn's split-kernel / derived-variable style, and
// with NO gauge template parameter (GRTeclyn's CCZ4RHS is not templated on the
// gauge -- the gauge is a separate object the Level constructs and calls):
//
//   FourDerivScalarTensor<coupling_and_potential_t, deriv_t>
//       - default ctor + explicit ctor(coupling_and_potential_t)
//       - nested Vars : public CCZ4Vars  (adds phi(), Pi())
//       - params_t::{check_params, fill_params}  (reads G_Newton)
//
//   ModifiedCCZ4RHS<theory_t, deriv_t> : public CCZ4RHS<deriv_t>
//       - ctor(amrex::Real dx); owns a theory_t and reads its params plus the
//         modified-gauge constants a0, b0 (modified_gauge.* scope)
//       - the RHS is built by the Level from these per-cell device methods,
//         following ScalarFieldLevel and the order of GRChombo's
//         ModifiedCCZ4RHS::compute(Cell):
//           compute_chi_and_h_ij              (inherited from CCZ4RHS)
//           compute_A_ij_and_Theta_and_Gamma  (inherited from CCZ4RHS)
//           add_a_and_b_rhs      -- a(x)/b(x) modified-gauge terms, called
//                                   just before add_emtensor_rhs
//           add_emtensor_rhs     -- kappa * T sources
//           add_matter_rhs       -- theory (phi, Pi) field evolution
//           solve_lhs            -- per-cell principal-part linear solve (4dST)
//           apply_dissipation    -- Kreiss-Oliger dissipation
//         each (int ix,int iy,int iz, Array4<Real> &rhs,
//               Array4<const Real> &state) const
//
//   ModifiedPunctureGauge<deriv_t>
//       - ctor(amrex::Real dx); standalone gauge object (like GRTeclyn's
//         MovingPunctureGauge / IntegratedMovingPunctureGauge)
//       - calculate_rhs(int,int,int, Array4<Real>&, Array4<const Real>&) const
//         sets the base moving-puncture lapse/shift/B RHS
//       - params_t with a0, b0 (+ standard gauge params) + check/fill_params
//
//   ModifiedGravityConstraints<theory_t>   -> derived record "constraints"
//   ModifiedGravityWeyl4<theory_t, deriv_t> -> derived record "Weyl4"
//   RhoDiagnostics<theory_t>               -> rho_phi/rho_g2/rho_g3/rho_GB
//       each with a static set_up(int state_index) + compute_mf(...)
//
//   ScalarExtraction : public SphericalExtraction<1>
//       - ctor(spherical_extraction_params_t, dt, time, first_step, restart)
//       - execute_query(ParticleInterpolator<1>*)   (adds state var c_phi)
//
#include "FourDerivScalarTensor.hpp"
#include "ModifiedCCZ4RHS.hpp"
#include "ModifiedGravityConstraints.hpp"
#include "ModifiedGravityWeyl4.hpp"
#include "ModifiedPunctureGauge.hpp"
#include "RhoDiagnostics.hpp"
#include "ScalarExtraction.hpp"

#include <array>
#include <string>
#include <type_traits>

namespace
{
// The full modified-CCZ4 + 4dST right hand side for one derivative order.
//
// GRTeclyn splits the CCZ4 RHS into several device kernels so that not all the
// first/second derivatives have to live in GPU registers at once (see
// ScalarFieldLevel and BinaryBHLevel). We follow the same pattern and add the
// modified-gravity pieces explicitly, in the order of GRChombo's
// ModifiedCCZ4RHS::compute(Cell):
//   vacuum CCZ4  ->  base moving-puncture gauge  ->  a(x)/b(x) gauge terms
//   ->  kappa * T sources  ->  theory (phi, Pi) evolution
//   ->  principal-part solve  ->  Kreiss-Oliger dissipation
template <class deriv_t>
void eval_rhs_impl(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                   amrex::Real a_dx)
{
    using theory_t = FourDerivScalarTensorWithCouplingAndPotential<deriv_t>;

    const auto &const_soln_arrays = a_soln.const_arrays();
    const auto &rhs_arrays        = a_rhs.arrays();

    const ModifiedCCZ4RHS<theory_t, deriv_t> modified_ccz4(a_dx);
    const ModifiedPunctureGauge<deriv_t> modified_puncture_gauge(a_dx);

    // 1. chi and h_ij
    amrex::ParallelFor(a_rhs,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       {
                           modified_ccz4.compute_chi_and_h_ij(
                               ix, iy, iz, rhs_arrays[box_no],
                               const_soln_arrays[box_no]);
                       });

    // 2. A_ij, Theta and Gamma
    amrex::ParallelFor(a_rhs,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       {
                           modified_ccz4.compute_A_ij_and_Theta_and_Gamma(
                               ix, iy, iz, rhs_arrays[box_no],
                               const_soln_arrays[box_no]);
                       });

    // 3. base gauge, modified-gauge a(x)/b(x) terms, matter sources, theory
    //    field evolution, principal-part solve and dissipation
    amrex::ParallelFor(
        a_rhs,
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
        {
            // base moving-puncture lapse/shift/B RHS (sets, does not add)
            modified_puncture_gauge.calculate_rhs(ix, iy, iz,
                                                  rhs_arrays[box_no],
                                                  const_soln_arrays[box_no]);
            // a(x)/b(x) modified-gauge terms, added just before the EM tensor
            modified_ccz4.add_a_and_b_rhs(ix, iy, iz, rhs_arrays[box_no],
                                          const_soln_arrays[box_no]);
            modified_ccz4.add_emtensor_rhs(ix, iy, iz, rhs_arrays[box_no],
                                           const_soln_arrays[box_no]);
            modified_ccz4.add_matter_rhs(ix, iy, iz, rhs_arrays[box_no],
                                         const_soln_arrays[box_no]);
            // solve the linear system for the fields that need it (4dST)
            modified_ccz4.solve_lhs(ix, iy, iz, rhs_arrays[box_no],
                                    const_soln_arrays[box_no]);
            modified_ccz4.apply_dissipation(ix, iy, iz, rhs_arrays[box_no],
                                            const_soln_arrays[box_no]);
        });
}
} // namespace

BinaryBH4dSTAmr *BinaryBH4dSTLevel::get_bh_amr_ptr()
{
    return dynamic_cast<BinaryBH4dSTAmr *>(get_gr_amr_ptr());
}

PunctureTracker<BinaryBH4dSTLevel::num_punctures> &
BinaryBH4dSTLevel::get_puncture_tracker()
{
    return get_bh_amr_ptr()->get_puncture_tracker();
}

void BinaryBH4dSTLevel::variableSetUp()
{
    BL_PROFILE("BinaryBH4dSTLevel::variableSetUp()");

    // Set up the state variables
    state_variable_set_up();

    using theory_t =
        FourDerivScalarTensorWithCouplingAndPotential<FourthOrderDerivatives>;

    // Register the modified-gravity diagnostics as AMReX derived records.
    ModifiedGravityConstraints<theory_t>::set_up(state_index);
    ModifiedGravityWeyl4<theory_t, FourthOrderDerivatives>::set_up(state_index);
    RhoDiagnostics<theory_t>::set_up(state_index);
}

// Things to do during the advance step after RK4 steps
void BinaryBH4dSTLevel::specific_advance()
{
    BL_PROFILE("BinaryBH4dSTLevel::specific_advance()");

    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto &state_arrays   = state_new.arrays();

    const AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;
    const PositiveChiAndLapse positive_chi_lapse;

    // Enforce det(h)=1, the trace-free A_ij condition and positive chi/lapse
    amrex::ParallelFor(state_new,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       {
                           algebraic_constraints_enforcer(
                               ix, iy, iz, state_arrays[box_no]);
                           positive_chi_lapse(ix, iy, iz,
                                              state_arrays[box_no]);
                       });
    amrex::Gpu::streamSynchronize();
}

// This initial data uses an approximation for the metric which is valid for
// small boosts
void BinaryBH4dSTLevel::initData()
{
    BL_PROFILE("BinaryBH4dSTLevel::initData");
    if (get_gr_amr_ptr()->Verbose() > 0)
    {
        amrex::Print() << "BinaryBH4dSTLevel::initData " << Level() << "\n";
    }

    const amrex::Real dx       = Geom().CellSize(0);
    amrex::MultiFab &state_new = get_new_data(state_index);

    const InitialScalarData initial_scalar_data(dx);
    static_assert(std::is_trivially_copyable_v<InitialScalarData>,
                  "InitialScalarData must be device copyable");

#ifdef USE_TWOPUNCTURES
    TwoPuncturesInitialData two_punctures_initial_data(dx);
    two_punctures_initial_data.solve(); // only solves the first time

#ifdef AMREX_USE_GPU
    amrex::MFInfo mf_info;
    mf_info.SetArena(amrex::The_Cpu_Arena());
    amrex::MultiFab host_state(state_new.boxArray(),
                               state_new.DistributionMap(), state_new.nComp(),
                               state_new.nGrowVect(), mf_info);
#else
    amrex::MultiFab &host_state = state_new;
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(state_new, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi)
    {
        const amrex::Box &grown_tile_box = mfi.growntilebox();
        const auto &state_array          = host_state.array(mfi);

        amrex::LoopOnCpu(grown_tile_box,
                         [=](int ix, int iy, int iz) {
                             two_punctures_initial_data(ix, iy, iz,
                                                        state_array);
                         });
#ifdef AMREX_USE_GPU
        amrex::Gpu::htod_memcpy_async(
            state_new[mfi].dataPtr(), host_state[mfi].dataPtr(),
            host_state[mfi].size() * sizeof(amrex::Real));
#endif
    }

    {
        // scalar field initial data (on device)
        const auto &state_arrays = state_new.arrays();
        amrex::ParallelFor(
            state_new, state_new.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            { initial_scalar_data(ix, iy, iz, state_arrays[box_no]); });
    }
#else
    const BinaryBHInitialData binary_initial_data(dx);
    static_assert(std::is_trivially_copyable_v<BinaryBHInitialData>,
                  "BinaryBHInitialData must be device copyable");

    // First set everything to zero (to avoid undefined values in constraints)
    // then calculate the initial data
    const auto &state_arrays = state_new.arrays();
    amrex::ParallelFor(
        state_new, state_new.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
        {
            amrex::CellData<amrex::Real> cell =
                state_arrays[box_no].cellData(ix, iy, iz);
            for (int n = 0; n < cell.nComp(); ++n)
            {
                cell[n] = 0.;
            }
            binary_initial_data(ix, iy, iz, state_arrays[box_no]);
            initial_scalar_data(ix, iy, iz, state_arrays[box_no]);
        });
#endif
    amrex::Gpu::streamSynchronize();

    if (get_bh_amr_ptr()->puncture_tracking_enabled() && Level() == 0)
    {
        // set the puncture coordinates: they are used for puncture tagging
        BoostedBHInitialData::params_t bh1_params(1);
        BoostedBHInitialData::params_t bh2_params(2);
#ifdef USE_TWOPUNCTURES
        two_punctures_initial_data.set_bh_params(bh1_params, bh2_params);
#else
        bh1_params.fill_params();
        bh2_params.fill_params();
#endif
        get_puncture_tracker().set_puncture_coords(
            {bh1_params.center[0], bh1_params.center[1], bh1_params.center[2],
             bh2_params.center[0], bh2_params.center[1], bh2_params.center[2]});
    }
}

// Calculate RHS during RK4 substeps
// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void BinaryBH4dSTLevel::specific_eval_rhs(amrex::MultiFab &a_soln,
                                          amrex::MultiFab &a_rhs,
                                          const amrex::Real /*a_time*/)
{
    BL_PROFILE("BinaryBH4dSTLevel::specific_eval_rhs()");

    const auto &soln_arrays = a_soln.arrays();
    const auto soln_ghosts  = a_soln.nGrowVect();
    const amrex::Real dx    = Geom().CellSize(0);

    // The classes to be used
    const AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;
    const PositiveChiAndLapse positive_chi_lapse;

    // Enforce positive chi and lapse, det(h)=1 and trace-free A on the solution
    amrex::ParallelFor(a_soln, soln_ghosts,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       {
                           algebraic_constraints_enforcer(
                               ix, iy, iz, soln_arrays[box_no]);
                           positive_chi_lapse(ix, iy, iz, soln_arrays[box_no]);
                       });

    if (m_evolution_spatial_derivative_order == 4)
    {
        eval_rhs_impl<FourthOrderDerivatives>(a_soln, a_rhs, dx);
    }
    else if (m_evolution_spatial_derivative_order == 6)
    {
        eval_rhs_impl<SixthOrderDerivatives>(a_soln, a_rhs, dx);
    }
    else
    {
        amrex::Abort("evolution.spatial_derivative_order must be 4 or 6");
    }

    amrex::Gpu::streamSynchronize();
}

// enforce algebraic constraints during RK4 substeps
// I think GRFolres doesn't do this?
void BinaryBH4dSTLevel::specific_update_ode(amrex::MultiFab &a_soln)
{
    BL_PROFILE("BinaryBH4dSTLevel::specific_update_ode()");

    const AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;
    const auto &soln_arrays = a_soln.arrays();

    amrex::ParallelFor(
        a_soln, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
        { algebraic_constraints_enforcer(ix, iy, iz, soln_arrays[box_no]); });
    amrex::Gpu::streamSynchronize();
}

void BinaryBH4dSTLevel::pre_tag_cells()
{
    BL_PROFILE("BinaryBH4dSTLevel::pre_tag_cells()");

    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto current_time    = get_state_data(state_index).curTime();

    // Only chi is used in the tagging criterion; 4th-order d2 needs 2 ghosts
    const int num_ghosts = 2;
    const int num_comps  = 1;
    FillPatch(*this, state_new, num_ghosts, current_time, state_index, c_chi,
              num_comps);
}

void BinaryBH4dSTLevel::tag_cells(amrex::TagBoxArray &a_tag_box_array,
                                  amrex::Real a_regrid_threshold)
{
    BL_PROFILE("BinaryBH4dSTLevel::tag_cells()");

    amrex::MultiFab &state_new = get_new_data(state_index);

    const auto &tag_arrays         = a_tag_box_array.arrays();
    const auto &state_const_arrays = state_new.const_arrays();

    const ChiTagger chi_tagger(Geom().CellSize(0), a_regrid_threshold);

    spherical_extraction_params_t extraction_params("weyl_extraction");
    extraction_params.fill_params();
    const ExtractionTagger extraction_tagger(Geom().CellSize(0), Level(),
                                             extraction_params);

    constexpr auto num_puncture_coords =
        static_cast<std::size_t>(AMREX_SPACEDIM * num_punctures);
    std::array<amrex::Real, num_puncture_coords> puncture_coords{};
    const bool puncture_tracking_enabled =
        get_bh_amr_ptr()->puncture_tracking_enabled();
    if (puncture_tracking_enabled)
    {
        puncture_coords = get_puncture_tracker().get_puncture_coords();
    }

    GRParmParse pp;
    amrex::Real bh1_mass{};
    amrex::Real bh2_mass{};
    pp.get("bh1.mass", bh1_mass);
    pp.get("bh2.mass", bh2_mass);

    const PunctureTagger<num_punctures> puncture_tagger(
        Geom().CellSize(0), Level(), get_gr_amr_ptr()->maxLevel(),
        puncture_coords, {bh1_mass, bh2_mass});

    amrex::ParallelFor(state_new, amrex::IntVect(0),
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       {
                           chi_tagger(ix, iy, iz, tag_arrays[box_no],
                                      state_const_arrays[box_no]);
                           extraction_tagger(ix, iy, iz, tag_arrays[box_no]);
                           if (puncture_tracking_enabled)
                           {
                               puncture_tagger(ix, iy, iz,
                                               tag_arrays[box_no]);
                           }
                       });
    amrex::Gpu::streamSynchronize();
}

void BinaryBH4dSTLevel::specific_post_init()
{
    BL_PROFILE("BinaryBH4dSTLevel::specific_post_init()");

    if (get_bh_amr_ptr()->puncture_tracking_enabled() && Level() == 0)
    {
        get_puncture_tracker().start_from_initial_punctures();
    }
}

void BinaryBH4dSTLevel::specific_post_restart()
{
    BL_PROFILE("BinaryBH4dSTLevel::specific_post_restart()");

    if (get_bh_amr_ptr()->puncture_tracking_enabled() && Level() == 0)
    {
        std::string restart_checkpoint{};
        GRParmParse pp("amr");
        pp.get("restart", restart_checkpoint);
        get_puncture_tracker().restart(restart_checkpoint);
    }
}

void BinaryBH4dSTLevel::specific_post_plotfile(const std::string &a_dir,
                                               std::ostream & /*a_os*/)
{
    if (get_bh_amr_ptr()->puncture_tracking_enabled() && Level() == 0)
    {
        get_puncture_tracker().write_plotfile(a_dir);
    }
}

void BinaryBH4dSTLevel::specific_post_checkpoint(const std::string &a_chk_dir,
                                                 std::ostream & /*a_os*/)
{
    if (get_bh_amr_ptr()->puncture_tracking_enabled() && Level() == 0)
    {
        get_puncture_tracker().checkpoint(a_chk_dir);
    }
}

void BinaryBH4dSTLevel::specific_post_timestep()
{
    BL_PROFILE("BinaryBH4dSTLevel::specific_post_timestep");

    // puncture tracking 
    if (get_bh_amr_ptr()->puncture_tracking_enabled())
    {
        GRParmParse puncture_tracking_pp("puncture_tracking");
        int puncture_tracking_level{};
        puncture_tracking_pp.get("level", puncture_tracking_level);
        int puncture_tracking_writeout_level{};
        puncture_tracking_pp.get("writeout_level",
                                 puncture_tracking_writeout_level);

        if (Level() == puncture_tracking_level)
        {
            BL_PROFILE("PunctureTracking");
            const bool write_punctures =
                at_level_timestep_multiple(puncture_tracking_writeout_level);
            const amrex::Real current_time =
                get_state_data(state_index).curTime();
            const amrex::Real dt = get_gr_amr_ptr()->dtLevel(Level());
            get_puncture_tracker().track(current_time, dt, write_punctures);
        }
    }

    // Weyl4 extraction
    spherical_extraction_params_t weyl_params("weyl_extraction");
    weyl_params.fill_params();
    if (weyl_params.enabled)
    {
        const int min_level = weyl_params.min_extraction_level();
        if (at_level_timestep_multiple(min_level) && Level() == min_level)
        {
            const amrex::Real m_time = get_state_data(state_index).curTime();
            const amrex::Real m_dt   = get_gr_amr_ptr()->dtLevel(Level());
            const amrex::Real restart_time =
                get_gr_amr_ptr()->get_restart_time();
            const bool first_step = (m_time <= m_dt);

            WeylExtraction my_extraction(weyl_params, m_dt, m_time, first_step,
                                        restart_time);
            my_extraction.execute_query(
                &get_bh_amr_ptr()->m_weyl_interpolator);
        }
    }

    // scalar-field extraction 
    spherical_extraction_params_t scalar_params("scalar_extraction");
    scalar_params.fill_params();
    if (scalar_params.enabled)
    {
        const int min_level = scalar_params.min_extraction_level();
        if (at_level_timestep_multiple(min_level) && Level() == min_level)
        {
            const amrex::Real m_time = get_state_data(state_index).curTime();
            const amrex::Real m_dt   = get_gr_amr_ptr()->dtLevel(Level());
            const amrex::Real restart_time =
                get_gr_amr_ptr()->get_restart_time();
            const bool first_step = (m_time <= m_dt);

            ScalarExtraction phi_extraction(scalar_params, m_dt, m_time,
                                            first_step, restart_time);
            phi_extraction.execute_query(
                &get_bh_amr_ptr()->m_scalar_interpolator);
        }
    }

    // NOTE: GRChombo's calculate_constraint_norms (needs AMRReductions) and the
    // apparent-horizon finder are not yet available in GRTeclyn, so I omit them
    // for now. The "constraints" derived record is still registered for
    // plotfile output.
}
