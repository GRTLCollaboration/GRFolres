/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(MODIFIEDGRAVITYCONSTRAINTS_HPP_)
#error                                                                         \
    "This file should only be included through ModifiedGravityConstraints.hpp"
#endif

#ifndef MODIFIEDGRAVITYCONSTRAINTS_IMPL_HPP_
#define MODIFIEDGRAVITYCONSTRAINTS_IMPL_HPP_
#include "DimensionDefinitions.hpp"

template <class matter_t>
ModifiedGravityConstraints<matter_t>::ModifiedGravityConstraints(
    const matter_t a_matter, double dx,
    const std::array<double, CH_SPACEDIM> a_center, double G_Newton,
    int a_c_Ham, const Interval &a_c_Moms, int a_c_Ham_abs_terms /* defaulted*/,
    const Interval &a_c_Moms_abs_terms /*defaulted*/)
    : Constraints(dx, a_c_Ham, a_c_Moms, a_c_Ham_abs_terms, a_c_Moms_abs_terms,
                  0.0 /*No cosmological constant*/),
      my_matter(a_matter), m_center(a_center), m_G_Newton(G_Newton) {}

template <class matter_t>
template <class data_t>
void ModifiedGravityConstraints<matter_t>::compute(
    Cell<data_t> current_cell) const {
  // Load local vars and calculate derivs
  const auto vars = current_cell.template load_vars<BSSNMatterVars>();
  const auto d1 = m_deriv.template diff1<BSSNMatterVars>(current_cell);
  const auto d2 = m_deriv.template diff2<BSSNMatterVars>(current_cell);

  // Inverse metric and Christoffel symbol
  const auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);
  const auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

  // Coordinates
  Coordinates<data_t> coords{current_cell, this->m_deriv.m_dx, m_center};

  // Get the non matter terms for the constraints
  Vars<data_t> out = constraint_equations(vars, d1, d2, h_UU, chris);

  // Energy Momentum Tensor
  rho_and_Si_t<data_t> rho_and_Si =
      my_matter.compute_rho_and_Si(vars, d1, d2, coords);

  // Hamiltonian constraint
  if (m_c_Ham >= 0 || m_c_Ham_abs_terms >= 0) {
    out.Ham += -16. * M_PI * m_G_Newton * rho_and_Si.rho;
    out.Ham_abs_terms += 16. * M_PI * m_G_Newton * abs(rho_and_Si.rho);
  }

  // Momentum constraints
  if (m_c_Moms.size() > 0 || m_c_Moms_abs_terms.size() > 0) {
    FOR(i) {
      out.Mom[i] += -8. * M_PI * m_G_Newton * rho_and_Si.Si[i];
      out.Mom_abs_terms[i] += 8. * M_PI * m_G_Newton * abs(rho_and_Si.Si[i]);
    }
  }
  // Write the constraints into the output FArrayBox
  store_vars(out, current_cell);
}

#endif /* MODIFIEDGRAVITYCONSTRAINTS_IMPL_HPP_ */
