/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "BinaryBHCubicHorndeskiLevel.hpp"
#include "AMRReductions.hpp"
#include "BinaryBH.hpp"
#include "BinarySFTaggingCriterion.hpp"
#include "BoxLoops.hpp"
#include "ComputePack.hpp"
#include "CouplingAndPotential.hpp"
#include "FourthOrderDerivatives.hpp"
#include "InitialScalarData.hpp"
#include "ModifiedCCZ4RHS.hpp"
#include "ModifiedGravityConstraints.hpp"
#include "ModifiedGravityWeyl4.hpp"
#include "NanCheck.hpp"
#include "PositiveChiAndAlpha.hpp"
#include "RhoDiagnostics.hpp"
#include "ScalarExtraction.hpp"
#include "SetValue.hpp"
#include "SixthOrderDerivatives.hpp"
#include "TraceARemoval.hpp"
#include "WeylExtraction.hpp"

void BinaryBHCubicHorndeskiLevel::postRestart()
{
    if (m_p.make_SF_zero_at_restart)
    {
        pout() << "Restarting at time " << m_time << " on level " << m_level
               << "." << std::endl;
        pout() << "Setting Scalar Field to 0 everywhere." << std::endl;
        BoxLoops::loop(SetValue(0., Interval(c_phi, c_Pi)), m_state_new,
                       m_state_new, INCLUDE_GHOST_CELLS);
    }
}

// Things to do at each advance step, after the RK4 is calculated
void BinaryBHCubicHorndeskiLevel::specificAdvance()
{
    // Enforce the trace free A_ij condition and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   m_state_new, m_state_new, INCLUDE_GHOST_CELLS);

    // Check for nan's
    if (m_p.nan_check)
        BoxLoops::loop(NanCheck("NaNCheck in specific Advance: "), m_state_new,
                       m_state_new, EXCLUDE_GHOST_CELLS, disable_simd());
}

// Initial data for field and metric variables
void BinaryBHCubicHorndeskiLevel::initialData()
{
    CH_TIME("BinaryBHCubicHorndeskiLevel::initialData");

    if (m_verbosity)
        pout() << "BinaryBHCubicHorndeskiLevel::initialData " << m_level
               << endl;

    // Set up the compute class for the BinaryBH initial data
    BinaryBH binary(m_p.bh1_params, m_p.bh2_params, m_dx);

    // First set everything to zero (to avoid undefinded values in constraints)
    // then calculate initial data
    BoxLoops::loop(
        make_compute_pack(SetValue(0.), binary,
                          InitialScalarData(m_p.initial_params, m_dx)),
        m_state_new, m_state_new, INCLUDE_GHOST_CELLS, disable_simd());
}

// Things to do before outputting a plot file
void BinaryBHCubicHorndeskiLevel::prePlotLevel()
{
    fillAllGhosts();
    CouplingAndPotential coupling_and_potential(
        m_p.coupling_and_potential_params);
    CubicHorndeskiWithCouplingAndPotential cubic_horndeski(
        coupling_and_potential);

    // constraints Ham and Mom were calculated in specificPostTimeStep already
    ModifiedGravityConstraints<CubicHorndeskiWithCouplingAndPotential>
        constraints(cubic_horndeski, m_dx, m_p.center, m_p.G_Newton, c_Ham,
                    Interval(c_Mom1, c_Mom3));
    ModifiedPunctureGauge modified_puncture_gauge(m_p.modified_ccz4_params);
    ModifiedGravityWeyl4<CubicHorndeskiWithCouplingAndPotential,
                         ModifiedPunctureGauge, FourthOrderDerivatives>
        weyl4(cubic_horndeski, m_p.modified_ccz4_params,
              modified_puncture_gauge, m_p.extraction_params.extraction_center,
              m_dx, m_p.sigma, CCZ4RHS<>::USE_CCZ4);
    // CCZ4 is required since this code only works in this formulation
    RhoDiagnostics<CubicHorndeskiWithCouplingAndPotential> rho_diagnostics(
        cubic_horndeski, m_dx, m_p.center);
    auto compute_pack = make_compute_pack(weyl4, constraints, rho_diagnostics);
    BoxLoops::loop(compute_pack, m_state_new, m_state_diagnostics,
                   EXCLUDE_GHOST_CELLS);
}

// Things to do in RHS update, at each RK4 step
void BinaryBHCubicHorndeskiLevel::specificEvalRHS(GRLevelData &a_soln,
                                                  GRLevelData &a_rhs,
                                                  const double a_time)
{
    // Enforce trace free A_ij and positive chi and alpha
    BoxLoops::loop(make_compute_pack(TraceARemoval(), PositiveChiAndAlpha()),
                   a_soln, a_soln, INCLUDE_GHOST_CELLS);

    // Calculate ModifiedCCZ4 right hand side with theory_t = CubicHorndeski
    CouplingAndPotential coupling_and_potential(
        m_p.coupling_and_potential_params);
    CubicHorndeskiWithCouplingAndPotential cubic_horndeski(
        coupling_and_potential);
    ModifiedPunctureGauge modified_puncture_gauge(m_p.modified_ccz4_params);

    if (m_p.max_spatial_derivative_order == 4)
    {
        ModifiedCCZ4RHS<CubicHorndeskiWithCouplingAndPotential,
                        ModifiedPunctureGauge, FourthOrderDerivatives>
            my_ccz4_theory(cubic_horndeski, m_p.modified_ccz4_params,
                           modified_puncture_gauge, m_dx, m_p.sigma, m_p.center,
                           m_p.G_Newton);

        BoxLoops::loop(my_ccz4_theory, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
    else if (m_p.max_spatial_derivative_order == 6)
    {
        ModifiedCCZ4RHS<CubicHorndeskiWithCouplingAndPotential,
                        ModifiedPunctureGauge, SixthOrderDerivatives>
            my_ccz4_theory(cubic_horndeski, m_p.modified_ccz4_params,
                           modified_puncture_gauge, m_dx, m_p.sigma, m_p.center,
                           m_p.G_Newton);

        BoxLoops::loop(my_ccz4_theory, a_soln, a_rhs, EXCLUDE_GHOST_CELLS);
    }
}

// Things to do at ODE update, after soln + rhs
void BinaryBHCubicHorndeskiLevel::specificUpdateODE(GRLevelData &a_soln,
                                                    const GRLevelData &a_rhs,
                                                    Real a_dt)
{
    // Enforce the trace free A_ij condition
    BoxLoops::loop(TraceARemoval(), a_soln, a_soln, INCLUDE_GHOST_CELLS);
}

void BinaryBHCubicHorndeskiLevel::preTagCells()
{
    // We only use chi, phi, Pi in the tagging criterion so only fill the ghosts
    // of those
    fillAllGhosts(VariableType::evolution, Interval(c_chi, c_chi));
    fillAllGhosts(VariableType::evolution, Interval(c_phi, c_Pi));
}

void BinaryBHCubicHorndeskiLevel::computeTaggingCriterion(
    FArrayBox &tagging_criterion, const FArrayBox &current_state)
{
    BoxLoops::loop(BinarySFTaggingCriterion(m_dx, m_level, m_p.max_level,
                                            m_time, m_p.tag_params),
                   current_state, tagging_criterion);
}

void BinaryBHCubicHorndeskiLevel::specificPostTimeStep()
{
    CH_TIME("BinaryBHCubicHorndeskiLevel::specificPostTimeStep");

    bool first_step = (m_time == 0.); // called at t=0 from Main
    // bool first_step = (m_time == m_dt); // not called in Main

    bool ghosts_filled = false;
    bool constraints_calculated = false;
    bool interpolator_refreshed = false;

    double finest_dt = m_dt * pow(2., m_level - m_p.max_level);

    if (m_p.activate_extraction == 1)
    {
        CH_TIME(
            "BinaryBHCubicHorndeskiLevel::specificPostTimeStep::extraction");
        int weyl_level = m_p.extraction_params.min_extraction_level();
        bool calculate_weyl = at_level_timestep_multiple(weyl_level);
        if (calculate_weyl)
        {
            // Populate the Weyl Scalar values on the grid
            fillAllGhosts();
            ghosts_filled = true;

            CouplingAndPotential coupling_and_potential(
                m_p.coupling_and_potential_params);
            CubicHorndeskiWithCouplingAndPotential cubic_horndeski(
                coupling_and_potential);
            ModifiedPunctureGauge modified_puncture_gauge(
                m_p.modified_ccz4_params);
            ModifiedGravityWeyl4<CubicHorndeskiWithCouplingAndPotential,
                                 ModifiedPunctureGauge, FourthOrderDerivatives>
                weyl4(cubic_horndeski, m_p.modified_ccz4_params,
                      modified_puncture_gauge,
                      m_p.extraction_params.extraction_center, m_dx, m_p.sigma,
                      CCZ4RHS<>::USE_CCZ4);
            // CCZ4 is required since this code only works in this
            // formulation
            BoxLoops::loop(weyl4, m_state_new, m_state_diagnostics,
                           EXCLUDE_GHOST_CELLS);
            // Do the extraction on the min extraction level
            if (m_level == weyl_level)
            {
                // Now refresh the interpolator and do the interpolation
                if (!interpolator_refreshed)
                {
                    // Now refresh the interpolator and do the interpolation
                    // fill ghosts manually to minimise communication
                    bool fill_ghosts = false;
                    m_bh_amr.m_interpolator->refresh(fill_ghosts);
                    m_bh_amr.fill_multilevel_ghosts(
                        VariableType::diagnostic,
                        Interval(c_Weyl4_Re, c_Weyl4_Im), weyl_level);
                    interpolator_refreshed = true;
                }

                WeylExtraction weyl_extraction(m_p.extraction_params, m_dt,
                                               m_time, first_step,
                                               m_restart_time);
                weyl_extraction.execute_query(m_bh_amr.m_interpolator);

                ScalarExtraction phi_extraction(m_p.scalar_extraction_params,
                                                m_dt, m_time, first_step,
                                                m_restart_time);
                phi_extraction.execute_query(m_gr_amr.m_interpolator);
            }
        }
    }

    if (m_p.calculate_constraint_norms)
    {
        CH_TIME("BinaryBHCubicHorndeskiLevel::specificPostTimeStep::constraint_"
                "norms");
        int Ham_level = 0;
        double Ham_dt = m_p.calculate_constraint_norm_interval * m_dt *
                        pow(2., m_level - Ham_level);
        bool calculate_Ham = fabs(remainder(m_time, Ham_dt)) < finest_dt / 2.;
        // bool calculate_Ham = at_level_timestep_multiple(Ham_level);

        if (calculate_Ham)
        {
            if (!constraints_calculated)
            {
                if (!ghosts_filled)
                {
                    fillAllGhosts();
                    ghosts_filled = true;
                }

                CouplingAndPotential coupling_and_potential(
                    m_p.coupling_and_potential_params);
                CubicHorndeskiWithCouplingAndPotential cubic_horndeski(
                    coupling_and_potential);
                fillAllGhosts();
                BoxLoops::loop(
                    ModifiedGravityConstraints<
                        CubicHorndeskiWithCouplingAndPotential>(
                        cubic_horndeski, m_dx, m_p.center, m_p.G_Newton, c_Ham,
                        Interval(c_Mom1, c_Mom3), c_Ham_abs,
                        Interval(c_Mom_abs1, c_Mom_abs3)),
                    m_state_new, m_state_diagnostics, EXCLUDE_GHOST_CELLS);
            }

            if (m_level == Ham_level)
            {
                AMRReductions<VariableType::diagnostic> amr_reductions(
                    m_bh_amr);
                double L2_Ham = amr_reductions.norm(c_Ham, 2, true);
                double L2_Ham_abs = amr_reductions.norm(c_Ham_abs, 2, true);
                double L2_Mom =
                    amr_reductions.norm(Interval(c_Mom1, c_Mom3), 2, true);
                double L2_Mom_abs = amr_reductions.norm(
                    Interval(c_Mom_abs1, c_Mom_abs3), 2, true);
                SmallDataIO constraints_file("constraint_norms", m_dt, m_time,
                                             m_restart_time,
                                             SmallDataIO::APPEND, first_step);
                constraints_file.remove_duplicate_time_data();
                if (first_step)
                    constraints_file.write_header_line(
                        {"L^2_Ham", "L^2_Mom", "L^2_Ham_abs", "L^2_Mom_abs"});
                constraints_file.write_time_data_line(
                    {L2_Ham, L2_Mom, L2_Ham_abs, L2_Mom_abs});
            }
        }
    }

#ifdef USE_AHFINDER
    if (m_p.AH_activate && m_level == m_p.AH_params.level_to_run)
    {
        if (m_p.AH_set_origins_to_punctures && m_p.track_punctures)
        {
            m_bh_amr.m_ah_finder.set_origins(
                m_bh_amr.m_puncture_tracker.get_puncture_coords());
        }
        m_bh_amr.m_ah_finder.solve(m_dt, m_time, m_restart_time);
    }
#endif
}
