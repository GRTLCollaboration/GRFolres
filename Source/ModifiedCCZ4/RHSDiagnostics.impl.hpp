/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(RHSDIAGNOSTICS_HPP_)
#error "This file should only be included through RHSDiagnostics.hpp"
#endif

#ifndef RHSDIAGNOSTICS_IMPL_HPP_
#define RHSDIAGNOSTICS_IMPL_HPP_

template <class theory_t, class gauge_t, class deriv_t>
RHSDiagnostics<theory_t, gauge_t, deriv_t>::RHSDiagnostics(
    theory_t a_theory, modified_params_t a_params, gauge_t a_gauge, double a_dx,
    double a_sigma, const std::array<double, CH_SPACEDIM> a_center,
    double a_G_Newton)
    : ModifiedCCZ4RHS<theory_t, gauge_t, deriv_t>(
          a_theory, a_params, a_gauge, a_dx, a_sigma, a_center, a_G_Newton),
      my_theory(a_theory), /*m_params(a_params)*/ my_gauge(a_gauge),
      m_center(a_center), m_G_Newton(a_G_Newton)
{
}

template <class theory_t, class gauge_t, class deriv_t>
template <class data_t>
void RHSDiagnostics<theory_t, gauge_t, deriv_t>::compute(
    Cell<data_t> current_cell) const
{
    // copy data from chombo gridpoint into local variables
    const auto theory_vars = current_cell.template load_vars<Vars>();
    const auto d1 = this->m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = this->m_deriv.template diff2<Vars>(current_cell);
    const auto advec =
        this->m_deriv.template advection<Vars>(current_cell, theory_vars.shift);

    // Call CCZ4 RHS - work out GR RHS, no dissipation
    Vars<data_t> theory_rhs;
    this->rhs_equation(theory_rhs, theory_vars, d1, d2, advec);

    Coordinates<data_t> coords{current_cell, this->m_deriv.m_dx, m_center};

    // add functions a(x) and b(x) of the modified gauge
    // this->my_modified_class.add_a_and_b_rhs(theory_rhs, theory_vars, d1, d2,
    // advec, coords);

    // add RHS theory terms from EM Tensor
    // this->add_emtensor_rhs(theory_rhs, theory_vars, d1, d2, advec, coords);

    // add evolution of theory fields themselves
    this->my_theory.add_theory_rhs(theory_rhs, theory_vars, d1, d2, advec,
                                   coords);

    // solve linear system for the theory fields that require it (e.g. 4dST)
    this->my_theory.solve_lhs(theory_rhs, theory_vars, d1, d2, advec, coords);

    // Compute diagnostics
    WeakCouplingConditions<data_t> weak_coupling_conditions =
        my_theory.compute_weak_coupling_conditions(theory_rhs, theory_vars, d1,
                                                   d2, advec, coords);

    // Write the constraints into the output FArrayBox
    current_cell.store_vars(weak_coupling_conditions.g2,
                            c_weak_coupling_condition_g2);
    current_cell.store_vars(weak_coupling_conditions.g3,
                            c_weak_coupling_condition_g3);
    current_cell.store_vars(weak_coupling_conditions.GB,
                            c_weak_coupling_condition_GB);
}

#endif /* RHSDIAGNOSTICS_IMPL_HPP_ */
