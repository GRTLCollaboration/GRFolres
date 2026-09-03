/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "KerrBH4dSTLevel.hpp"

#include "AlgebraicConstraintsEnforcer.hpp"
#include "EMTensor.hpp"
#include "FixedGridsTagger.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GammaCalculator.hpp"
#include "IntegratedMovingPunctureGauge.hpp"
#include "KerrBHInitialData.hpp"
#include "LineExtraction.hpp"
#include "ModifiedCCZ4RHS.hpp"
#include "ModifiedGravityConstraints.hpp"
#include "ModifiedGravityWeyl4.hpp"
#include "PositiveChiAndLapse.hpp"
#include "ScalarFieldInitialData.hpp"
#include "SixthOrderDerivatives.hpp"
#include "StateTypes.hpp"

#include <type_traits>

using theory_t = KerrBH4dSTLevel::KerrBH4dSTWithCouplingAndPotential<>;

using KerrBH4dSTEnergyDensity =
    EMTensor<theory_t, EMTensorOptions::justEnergyDensity>;
using KerrBH4dSTConstraints = ModifiedGravityConstraints<theory_t>;
using KerrBH4dSTWeyl4 = ModifiedGravityWeyl4<theory_t>;

KerrBH4dSTAmr *KerrBH4dSTLevel::get_kerrbh4dst_amr_ptr()
{
    return dynamic_cast<KerrBH4dSTAmr *>(get_gr_amr_ptr());
}

void KerrBH4dSTLevel::variableSetUp()
{
    BL_PROFILE("KerrBH4dSTLevel::variableSetUp()");
    state_variable_set_up();

    KerrBH4dSTConstraints::set_up(state_index);
    KerrBH4dSTWeyl4::set_up(state_index);
    KerrBH4dSTEnergyDensity::set_up(state_index);
}

void KerrBH4dSTLevel::specific_advance()
{
    BL_PROFILE("KerrBH4dSTLevel::specific_advance()");

    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto &state_arrays = state_new.arrays();

    const AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;
    const PositiveChiAndLapse positive_chi_and_lapse;

    amrex::ParallelFor(
        state_new,
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
        {
            algebraic_constraints_enforcer(ix, iy, iz, state_arrays[box_no]);
            positive_chi_and_lapse(ix, iy, iz, state_arrays[box_no]);
        });
    amrex::Gpu::streamSynchronize();
}

void KerrBH4dSTLevel::specific_post_timestep()
{
    BL_PROFILE("KerrBH4dSTLevel::specific_post_timestep()");

    if (Level() == 0)
    {
        const amrex::Real time = get_state_data(state_index).curTime();
        const amrex::Real dt = get_gr_amr_ptr()->dtLevel(0);
        const amrex::Real restart_time = get_gr_amr_ptr()->get_restart_time();
        const bool first_step = (time <= dt);

        const LineExtraction<1> phi_extraction("phi_line_extraction", c_phi, dt,
                                               time, restart_time, first_step);
        phi_extraction.execute_query(
            &get_kerrbh4dst_amr_ptr()->phi_interpolator);

        const LineExtraction<1> rho_extraction("rho_line_extraction", 0, dt,
                                               time, restart_time, first_step);
        const LineExtraction<1>::derived_vars_t rho_vars{
            KerrBH4dSTEnergyDensity::name, {"rho"}, {BCParity::even}};

        rho_extraction.execute_query(
            &get_kerrbh4dst_amr_ptr()->rho_interpolator, rho_vars);
    }
}

void KerrBH4dSTLevel::initData()
{
    BL_PROFILE("KerrBH4dSTLevel::initData()");

    if (get_gr_amr_ptr()->Verbose() > 0)
    {
        amrex::Print() << "KerrBH4dSTLevel::initData " << Level() << "\n";
    }

    const amrex::Real dx = Geom().CellSize(0);
    const KerrBHInitialData kerr_bh_initial_data(dx);
    static_assert(std::is_trivially_copyable_v<KerrBHInitialData>,
                  "KerrBHInitialData needs to be device copyable");

    const ScalarFieldInitialData scalar_field_initial_data(Geom().CellSize(0));
    static_assert(std::is_trivially_copyable_v<ScalarFieldInitialData>,
                  "ScalarFieldInitialData must be device copyable");

    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto &state_arrays = state_new.arrays();

    amrex::ParallelFor(
        state_new, state_new.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
        {
            const amrex::CellData<amrex::Real> cell =
                state_arrays[box_no].cellData(ix, iy, iz);
            for (int component = 0; component < cell.nComp(); ++component)
            {
                cell[component] = 0.0;
            }
            scalar_field_initial_data(ix, iy, iz, state_arrays[box_no]);
            kerr_bh_initial_data(ix, iy, iz, state_arrays[box_no]);
        });

    if (m_evolution_spatial_derivative_order == 4)
    {
        const GammaCalculator<FourthOrderDerivatives> gamma_calculator(
            Geom().CellSize(0));
        const IntegratedMovingPunctureGauge<FourthOrderDerivatives>
            integrated_moving_puncture_gauge(Geom().CellSize(0));
        amrex::ParallelFor(
            state_new,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                gamma_calculator(ix, iy, iz, state_arrays[box_no]);
                integrated_moving_puncture_gauge.set_initial_B_to_Gamma(
                    ix, iy, iz, state_arrays[box_no]);
            });
    }
    else if (m_evolution_spatial_derivative_order == 6)
    {
        const GammaCalculator<SixthOrderDerivatives> gamma_calculator(
            Geom().CellSize(0));
        const IntegratedMovingPunctureGauge<SixthOrderDerivatives>
            integrated_moving_puncture_gauge(Geom().CellSize(0));
        amrex::ParallelFor(
            state_new,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                gamma_calculator(ix, iy, iz, state_arrays[box_no]);
                integrated_moving_puncture_gauge.set_initial_B_to_Gamma(
                    ix, iy, iz, state_arrays[box_no]);
            });
    }

    amrex::Gpu::streamSynchronize();
}

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
void KerrBH4dSTLevel::specific_eval_rhs(amrex::MultiFab &a_soln,
                                        amrex::MultiFab &a_rhs,
                                        const amrex::Real /*a_time*/)
{
    BL_PROFILE("KerrBH4dSTLevel::specific_eval_rhs()");

    const auto &soln_arrays = a_soln.arrays();
    const auto &const_soln_arrays = a_soln.const_arrays();
    const auto &rhs_arrays = a_rhs.arrays();

    const AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;
    const PositiveChiAndLapse positive_chi_and_lapse;

    amrex::ParallelFor(
        a_soln, a_soln.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
        {
            algebraic_constraints_enforcer(ix, iy, iz, soln_arrays[box_no]);
            positive_chi_and_lapse(ix, iy, iz, soln_arrays[box_no]);
        });

    if (m_evolution_spatial_derivative_order != 4 &&
        m_evolution_spatial_derivative_order != 6)
    {
        amrex::Abort("spatial_derivative_order must be 4 or 6");
    }

    // NOLINTBEGIN(bugprone-branch-clone)
    if (m_evolution_spatial_derivative_order == 4)
    {
        const ModifiedCCZ4RHS<
            KerrBH4dSTWithCouplingAndPotential<FourthOrderDerivatives>,
            FourthOrderDerivatives>
            modified_ccz4_rhs(Geom().CellSize(0));
        const IntegratedMovingPunctureGauge<FourthOrderDerivatives>
            integrated_moving_puncture_gauge(Geom().CellSize(0));

        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                modified_ccz4_rhs.compute_chi_and_h_ij(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
            });

        // NOLINTBEGIN(bugprone-easily-swappable-parameters)
        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                modified_ccz4_rhs.compute_A_ij_and_Theta_and_Gamma(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
            });
        // NOLINTEND(bugprone-easily-swappable-parameters)

        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                modified_ccz4_rhs.add_emtensor_rhs(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
                integrated_moving_puncture_gauge.calculate_rhs(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
                modified_ccz4_rhs.add_theory_rhs(ix, iy, iz, rhs_arrays[box_no],
                                                 const_soln_arrays[box_no]);
                modified_ccz4_rhs.apply_dissipation(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
                modified_ccz4_rhs.add_a_and_b_rhs(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
            });
    }
    else if (m_evolution_spatial_derivative_order == 6)
    {
        const ModifiedCCZ4RHS<
            KerrBH4dSTWithCouplingAndPotential<SixthOrderDerivatives>,
            SixthOrderDerivatives>
            modified_ccz4_rhs(Geom().CellSize(0));
        const IntegratedMovingPunctureGauge<SixthOrderDerivatives>
            integrated_moving_puncture_gauge(Geom().CellSize(0));

        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                modified_ccz4_rhs.compute_chi_and_h_ij(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
            });

        // NOLINTBEGIN(bugprone-easily-swappable-parameters)
        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                modified_ccz4_rhs.compute_A_ij_and_Theta_and_Gamma(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
            });
        // NOLINTEND(bugprone-easily-swappable-parameters)

        amrex::ParallelFor(
            a_rhs,
            [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
            {
                modified_ccz4_rhs.add_emtensor_rhs(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
                integrated_moving_puncture_gauge.calculate_rhs(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
                modified_ccz4_rhs.add_theory_rhs(ix, iy, iz, rhs_arrays[box_no],
                                                 const_soln_arrays[box_no]);
                modified_ccz4_rhs.apply_dissipation(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
                modified_ccz4_rhs.add_a_and_b_rhs(
                    ix, iy, iz, rhs_arrays[box_no], const_soln_arrays[box_no]);
            });
    }
    // NOLINTEND(bugprone-branch-clone)

    amrex::Gpu::streamSynchronize();
}

void KerrBH4dSTLevel::specific_update_ode(amrex::MultiFab &a_soln)
{
    BL_PROFILE("KerrBH4dSTLevel::specific_update_ode()");

    const auto &soln_arrays = a_soln.arrays();
    const AlgebraicConstraintsEnforcer algebraic_constraints_enforcer;

    amrex::ParallelFor(
        a_soln, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
        { algebraic_constraints_enforcer(ix, iy, iz, soln_arrays[box_no]); });
    amrex::Gpu::streamSynchronize();
}

void KerrBH4dSTLevel::tag_cells(amrex::TagBoxArray &a_tag_box_array,
                                const amrex::Real /*a_regrid_threshold*/)
{
    BL_PROFILE("KerrBH4dSTLevel::tag_cells()");

    const auto &tag_arrays = a_tag_box_array.arrays();

    const FixedGridsTagger tagger(Geom().CellSize(0), Level());

    amrex::ParallelFor(a_tag_box_array,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       { tagger(ix, iy, iz, tag_arrays[box_no]); });
    amrex::Gpu::streamSynchronize();
}
