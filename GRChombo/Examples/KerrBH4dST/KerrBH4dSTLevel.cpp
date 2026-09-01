/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "KerrBH4dSTLevel.hpp"
#include "AMRReductions.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "FixedGridsTaggingCriterion.hpp"
#include "InitialScalarData.hpp"
#include "ModifiedCCZ4RHS.hpp"
#include "ModifiedGravityConstraints.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "RhoDiagnostics.hpp"
#include "SetValue.hpp"
#include "SixthOrderDerivatives.hpp"
#include "SmallDataIO.hpp"
#include "TraceARemoval.hpp"

// Initial data
#include "GammaCalculator.hpp"
#include "KerrBH.hpp"

void KerrBH4dSTLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck(), m_state_new, m_state_new,
                       EXCLUDE_GHOST_CELLS, disable_simd());
}

void KerrBH4dSTLevel::initialData()
{
    CH_TIME("KerrBH4dSTLevel::initialData");
    if (m_verbosity)
        pout() << "KerrBH4dSTLevel::initialData " << m_level << endl;

    // First set everything to zero then calculate initial data  Get the Kerr
    // solution in the variables, then calculate the \tilde\Gamma^i numerically
    // as these are non zero and not calculated in the Kerr ICs
    BoxLoops::loop(
        make_compute_pack(SetValue(0.), KerrBH(m_p.kerr_params, m_dx),
                          InitialScalarData(m_p.initial_params, m_dx)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    fillAllGhosts();
    BoxLoops::loop(GammaCalculator(m_dx), m_state_new, m_state_new,
                   EXCLUDE_GHOST_CELLS);

#ifdef USE_AHFINDER
    // Diagnostics needed for AHFinder
    CouplingAndPotential coupling_and_potential(
        m_p.coupling_and_potential_params);
    FourDerivScalarTensorWithCouplingAndPotential fdst(coupling_and_potential,
                                                       m_p.G_Newton);
    ModifiedGravityConstraints<FourDerivScalarTensorWithCouplingAndPotential>
        constraints(fdst, m_dx, m_p.center, m_p.G_Newton, c_Ham,
                    Interval(c_Mom1, c_Mom3));
    BoxLoops::loop(constraints, m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);
#endif
}

#ifdef CH_USE_HDF5
void KerrBH4dSTLevel::prePlotLevel()
{
#ifdef USE_AHFINDER
    // already calculated in 'specificPostTimeStep'
    if (m_bh_amr.m_ah_finder.need_diagnostics(m_dt, m_time))
        return;
#endif

    CouplingAndPotential coupling_and_potential(
        m_p.coupling_and_potential_params);
    FourDerivScalarTensorWithCouplingAndPotential fdst(coupling_and_potential,
                                                       m_p.G_Newton);
    fillAllGhosts();
    ModifiedGravityConstraints<FourDerivScalarTensorWithCouplingAndPotential>
        constraints(fdst, m_dx, m_p.center, m_p.G_Newton, c_Ham,
                    Interval(c_Mom1, c_Mom3));
    BoxLoops::loop(constraints, m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);
}
#endif /* CH_USE_HDF5 */

void KerrBH4dSTLevel::specificEvalRHS(GRLevelData &a_soln, GRLevelData &a_rhs,
                                      const double a_time)
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate ModifiedCCZ4 right hand side with theory_t =
    // FourDerivScalarTensor
    CouplingAndPotential coupling_and_potential(
        m_p.coupling_and_potential_params);
    FourDerivScalarTensorWithCouplingAndPotential fdst(coupling_and_potential,
                                                       m_p.G_Newton);
    ModifiedPunctureGauge modified_puncture_gauge(m_p.modified_ccz4_params);
    if (m_p.max_spatial_derivative_order == 4)
    {
        ModifiedCCZ4RHS<FourDerivScalarTensorWithCouplingAndPotential,
                        ModifiedPunctureGauge, FourthOrderDerivatives>
            my_modified_ccz4(fdst, m_p.modified_ccz4_params,
                             modified_puncture_gauge, m_dx, m_p.sigma,
                             m_p.center, m_p.G_Newton);
        BoxLoops::loop(my_modified_ccz4, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        ModifiedCCZ4RHS<FourDerivScalarTensorWithCouplingAndPotential,
                        ModifiedPunctureGauge, SixthOrderDerivatives>
            my_modified_ccz4(fdst, m_p.modified_ccz4_params,
                             modified_puncture_gauge, m_dx, m_p.sigma,
                             m_p.center, m_p.G_Newton);
        BoxLoops::loop(my_modified_ccz4, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

void KerrBH4dSTLevel::specificUpdateODE(GRLevelData &a_soln,
                                        const GRLevelData &a_rhs, Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void KerrBH4dSTLevel::preTagCells()
{
    // we don't need any ghosts filled for the fixed grids tagging criterion
    // used here so don't fill any
}

void KerrBH4dSTLevel::computeTaggingCriterion(FArrayBox &tagging_criterion,
                                              const FArrayBox &current_state)
{
    BoxLoops::loop(FixedGridsTaggingCriterion(m_dx, m_level, m_p.L, m_p.center),
                   current_state, tagging_criterion);
}

void KerrBH4dSTLevel::specificPostTimeStep()
{
    CH_TIME("KerrBHLevel::specificPostTimeStep");

    bool first_step =
        (m_time == 0.); // this form is used when 'specificPostTimeStep' was
                        // called during setup at t=0 from Main
                        //     // bool first_step = (m_time == m_dt); // if not
                        //     called in Main

    fillAllGhosts();
    CouplingAndPotential coupling_and_potential(
        m_p.coupling_and_potential_params);
    FourDerivScalarTensorWithCouplingAndPotential fdst(coupling_and_potential,
                                                       m_p.G_Newton);
    ModifiedPunctureGauge modified_puncture_gauge(m_p.modified_ccz4_params);
    RhoDiagnostics<FourDerivScalarTensorWithCouplingAndPotential>
        rho_diagnostics(fdst, m_dx, m_p.center);
    BoxLoops::loop(rho_diagnostics, m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);

    if (m_p.calculate_diagnostic_norms)
    {
        CouplingAndPotential coupling_and_potential(
            m_p.coupling_and_potential_params);
        FourDerivScalarTensorWithCouplingAndPotential fdst(
            coupling_and_potential, m_p.G_Newton);
        fillAllGhosts();
        BoxLoops::loop(ModifiedGravityConstraints<
                           FourDerivScalarTensorWithCouplingAndPotential>(
                           fdst, m_dx, m_p.center, m_p.G_Newton, c_Ham,
                           Interval(c_Mom1, c_Mom3)),
                       m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
        if (m_level == 0)
        {
            AMRReductions<VariableType::diagnostic> amr_reductions(m_gr_amr);
            bool normalise_by_volume = true;
            double L2_Ham = amr_reductions.norm(c_Ham, 2, normalise_by_volume);
            double L2_Mom = amr_reductions.norm(Interval(c_Mom1, c_Mom3), 2,
                                                normalise_by_volume);
            double L2_rho_phi =
                amr_reductions.norm(c_rho_phi, 2, normalise_by_volume);
            SmallDataIO diagnostics_file(m_p.data_path + "diagnostic_norms",
                                         m_dt, m_time, m_restart_time,
                                         SmallDataIO::APPEND, first_step);
            diagnostics_file.remove_duplicate_time_data();
            if (first_step)
            {
                diagnostics_file.write_header_line(
                    {"L^2_Ham", "L^2_Mom", "L^2_rho_phi"});
            }
            diagnostics_file.write_time_data_line({L2_Ham, L2_Mom, L2_rho_phi});
        }
    }

#ifdef USE_AHFINDER
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
#endif
}
