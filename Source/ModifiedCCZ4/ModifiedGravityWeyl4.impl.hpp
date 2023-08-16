/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(MODIFIEDGRAVITYWEYL4_HPP_)
#error "This file should only be included through ModifiedGravityWeyl4.hpp"
#endif

#ifndef MODIFIEDGRAVITYWEYL4_IMPL_HPP_
#define MODIFIEDGRAVITYWEYL4_IMPL_HPP_

template <class theory_t, class gauge_t, class deriv_t>
template <class data_t>
void ModifiedGravityWeyl4<theory_t, gauge_t, deriv_t>::compute(
    Cell<data_t> current_cell) const
{

    // copy data from chombo gridpoint into local variables
    const auto vars = current_cell.template load_vars<Vars>();
    const auto d1 = m_deriv.template diff1<Vars>(current_cell);
    const auto d2 = m_deriv.template diff2<Diff2Vars>(current_cell);
    const auto advec =
        m_deriv.template advection<Vars>(current_cell, vars.shift);

    // Get the coordinates
    const Coordinates<data_t> coords(current_cell, m_dx, m_center);

    // Compute the inverse metric
    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // Compute the spatial volume element
    const auto epsilon3_LUU = compute_epsilon3_LUU(vars, h_UU);

    // Compute the E and B fields
    EBFields_t<data_t> ebfields =
        compute_EB_fields(vars, d1, d2, epsilon3_LUU, h_UU, chris);

    // Add modified gravity terms to E and B fields
    add_theory_EB(ebfields, vars, d1, d2, advec, coords);

    // work out the Newman Penrose scalar
    NPScalar_t<data_t> out =
        compute_Weyl4(ebfields, vars, d1, d2, h_UU, coords);

    // Write the rhs into the output FArrayBox
    current_cell.store_vars(out.Real, c_Weyl4_Re);
    current_cell.store_vars(out.Im, c_Weyl4_Im);
}

template <class theory_t, class gauge_t, class deriv_t>
template <class data_t>
void ModifiedGravityWeyl4<theory_t, gauge_t, deriv_t>::add_theory_EB(
    EBFields_t<data_t> &ebfields, const Vars<data_t> &vars,
    const Vars<Tensor<1, data_t>> &d1, const Diff2Vars<Tensor<2, data_t>> &d2,
    const Vars<data_t> &advec, const Coordinates<data_t> &coords) const
{
    Tensor<2, data_t> kappa_times_Sij_TF =
        this->get_full_kappa_times_Sij_TF(vars, d1, d2, advec, coords);
    // as we made the vacuum expression of Bij explictly symmetric and Eij
    // explictly trace-free, only Eij has matter terms
    FOR(i, j) { ebfields.E[i][j] += -0.5 * kappa_times_Sij_TF[i][j]; }
}

#endif /* MODIFIEDGRAVITYWEYL4_IMPL_HPP_ */
