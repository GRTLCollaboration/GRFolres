/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(MODIFIEDCCZ4RHS_HPP_)
#error "This file should only be included through ModifiedCCZ4RHS.hpp"
#endif

#ifndef MODIFIEDCCZ4RHS_IMPL_HPP_
#define MODIFIEDCCZ4RHS_IMPL_HPP_
#include "DimensionDefinitions.hpp"

template <class theory_t, class gauge_t, class deriv_t>
ModifiedCCZ4RHS<theory_t, gauge_t, deriv_t>::ModifiedCCZ4RHS(
    theory_t a_theory, modified_params_t a_params, gauge_t a_gauge, double a_dx,
    double a_sigma, const std::array<double, CH_SPACEDIM> a_center,
    int a_formulation, double a_G_Newton)
    : CCZ4RHS<gauge_t, deriv_t>(a_params, a_dx, a_sigma, a_formulation,
                                0.0 /*No cosmological constant*/),
      my_theory(a_theory), my_gauge(a_gauge), m_center(a_center),
      m_G_Newton(a_G_Newton) {}

template <class theory_t, class gauge_t, class deriv_t>
template <class data_t>
void ModifiedCCZ4RHS<theory_t, gauge_t, deriv_t>::compute(
    Cell<data_t> current_cell) const {
  // copy data from chombo gridpoint into local variables
  const auto theory_vars = current_cell.template load_vars<Vars>();
  const auto d1 = this->m_deriv.template diff1<Vars>(current_cell);
  const auto d2 = this->m_deriv.template diff2<Diff2Vars>(current_cell);
  const auto advec =
      this->m_deriv.template advection<Vars>(current_cell, theory_vars.shift);

  // Call CCZ4 RHS - work out GR RHS, no dissipation
  Vars<data_t> theory_rhs;
  this->rhs_equation(theory_rhs, theory_vars, d1, d2, advec);

  Coordinates<data_t> coords{current_cell, this->m_deriv.m_dx, m_center};

  // add functions a(x) and b(x) of the modified gauge
  add_a_and_b_rhs(theory_rhs, theory_vars, d1, d2, advec, coords);

  // add RHS theory terms from EM Tensor
  add_emtensor_rhs(theory_rhs, theory_vars, d1, d2, advec, coords);

  // add evolution of theory fields themselves
  my_theory.add_theory_rhs(theory_rhs, theory_vars, d1, d2, advec, coords);

  // solve linear system for the theory fields that require it (e.g. 4dST)
  my_theory.solve_lhs(theory_rhs, theory_vars, d1, d2, advec, coords);

  // Add dissipation to all terms
  this->m_deriv.add_dissipation(theory_rhs, current_cell, this->m_sigma);

  // Write the rhs into the output FArrayBox
  current_cell.store_vars(theory_rhs);
}

// Function to add a(x) and b(x) from the modified gauge to CCZ4 rhs
template <class theory_t, class gauge_t, class deriv_t>
template <class data_t>
void ModifiedCCZ4RHS<theory_t, gauge_t, deriv_t>::add_a_and_b_rhs(
    Vars<data_t> &theory_rhs, const Vars<data_t> &theory_vars,
    const Vars<Tensor<1, data_t>> &d1, const Diff2Vars<Tensor<2, data_t>> &d2,
    const Vars<data_t> &advec, const Coordinates<data_t> &coords) const {
  data_t a_of_x = 0.;
  data_t b_of_x = 0.;

  my_gauge.compute_a_and_b(a_of_x, b_of_x, coords);

  const data_t chi_regularised = simd_max(1e-6, theory_vars.chi);
  using namespace TensorAlgebra;
  auto h_UU = compute_inverse_sym(theory_vars.h);
  auto chris = compute_christoffel(d1.h, h_UU);
  Tensor<1, data_t> Z_over_chi;
  Tensor<1, data_t> Z;
  if (this->m_formulation == CCZ4RHS<>::USE_BSSN) {
    FOR(i) Z_over_chi[i] = 0.0;
  } else {
    FOR(i)
    Z_over_chi[i] = 0.5 * (theory_vars.Gamma[i] - chris.contracted[i]);
  }
  auto ricci0 = CCZ4Geometry::compute_ricci_Z(theory_vars, d1, d2, h_UU, chris,
                                              {0., 0., 0.});
  // for some reason the function compute_ricci does not give the correct
  // answer

  Tensor<2, data_t> A_UU = raise_all(theory_vars.A, h_UU);
  // A^{ij} A_{ij}. - Note the abuse of the compute trace function.
  data_t tr_A2 = compute_trace(theory_vars.A, A_UU);

  Tensor<2, data_t> covdtilde_A[CH_SPACEDIM];

  FOR(i, j, k) {
    covdtilde_A[i][j][k] = d1.A[j][k][i];
    FOR(l) {
      covdtilde_A[i][j][k] += -chris.ULL[l][i][j] * theory_vars.A[l][k] -
                              chris.ULL[l][i][k] * theory_vars.A[l][j];
    }
  }

  Tensor<1, data_t> Ni;
  FOR(i) { Ni[i] = -(GR_SPACEDIM - 1.) * d1.K[i] / (double)GR_SPACEDIM; }
  FOR(i, j, k) {
    Ni[i] += h_UU[j][k] *
             (covdtilde_A[k][j][i] - GR_SPACEDIM * theory_vars.A[i][j] *
                                         d1.chi[k] / (2. * chi_regularised));
  }

  data_t kappa1_times_lapse;
  if (this->m_params.covariantZ4)
    kappa1_times_lapse = this->m_params.kappa1;
  else
    kappa1_times_lapse = this->m_params.kappa1 * theory_vars.lapse;

  if (this->m_formulation == CCZ4RHS<>::USE_BSSN) {
    theory_rhs.K += 0.;
  } else {
    theory_rhs.K +=
        b_of_x / (1. + b_of_x) *
        (theory_vars.lapse * (-0.5 * theory_vars.K * theory_vars.K +
                              0.5 * GR_SPACEDIM / (GR_SPACEDIM - 1.) *
                                  (tr_A2 - ricci0.scalar)) +
         kappa1_times_lapse * GR_SPACEDIM * theory_vars.Theta *
             (1. + this->m_params.kappa2 -
              0.5 * ((GR_SPACEDIM - 3.) / (GR_SPACEDIM - 1.) +
                     this->m_params.kappa2)));

    theory_rhs.Theta +=
        b_of_x / (1. + b_of_x) *
        (0.5 * theory_vars.lapse *
             (tr_A2 - ricci0.scalar -
              ((GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) * theory_vars.K *
                  theory_vars.K) +
         0.5 * theory_vars.Theta * kappa1_times_lapse *
             ((GR_SPACEDIM + 1.) + this->m_params.kappa2 * (GR_SPACEDIM - 1.)));
  }

  FOR(i) {
    theory_rhs.Gamma[i] +=
        b_of_x / (1. + b_of_x) *
        ((2.0 / (double)GR_SPACEDIM) *
             (theory_vars.lapse * theory_vars.K * Z_over_chi[i]) +
         2. * kappa1_times_lapse * Z_over_chi[i]);
    FOR(j) {
      theory_rhs.Gamma[i] += b_of_x / (1. + b_of_x) * 2. * h_UU[i][j] *
                             theory_vars.lapse * (-d1.Theta[j] - Ni[j]);
      FOR(k) {
        theory_rhs.Gamma[i] += b_of_x / (1. + b_of_x) * 2. * theory_vars.lapse *
                               A_UU[i][j] * theory_vars.h[k][j] * Z_over_chi[k];
      }
    }
  }

  theory_rhs.lapse += a_of_x / (1. + a_of_x) * this->m_params.lapse_coeff *
                      pow(theory_vars.lapse, this->m_params.lapse_power) *
                      (theory_vars.K - 2. * theory_vars.Theta);

  FOR(i) {
    theory_rhs.shift[i] -= a_of_x / (1. + a_of_x) *
                           this->m_params.shift_Gamma_coeff *
                           theory_vars.Gamma[i];
    FOR(j) {
      theory_rhs.shift[i] -= a_of_x / (1. + a_of_x) * theory_vars.lapse *
                             theory_vars.chi * h_UU[i][j] * d1.lapse[j];
    }
  }
}

// Function to add in EM Tensor theory terms to CCZ4 rhs
template <class theory_t, class gauge_t, class deriv_t>
template <class data_t>
void ModifiedCCZ4RHS<theory_t, gauge_t, deriv_t>::add_emtensor_rhs(
    Vars<data_t> &theory_rhs, const Vars<data_t> &theory_vars,
    const Vars<Tensor<1, data_t>> &d1, const Diff2Vars<Tensor<2, data_t>> &d2,
    const Vars<data_t> &advec, const Coordinates<data_t> &coords) const {
  using namespace TensorAlgebra;

  const auto h_UU = compute_inverse_sym(theory_vars.h);
  const auto chris = compute_christoffel(d1.h, h_UU);

  // Calculate elements of the decomposed stress energy tensor
  RhoAndSi<data_t> rho_and_Si =
      my_theory.compute_rho_and_Si(theory_vars, d1, d2, coords);
  SijTFAndS<data_t> Sij_TF_and_S =
      my_theory.compute_Sij_TF_and_S(theory_vars, d1, d2, advec, coords);

  data_t a_of_x = 0.;
  data_t b_of_x = 0.;

  my_gauge.compute_a_and_b(a_of_x, b_of_x, coords);

  // Update RHS for K and Theta depending on formulation
  if (this->m_formulation == CCZ4RHS<>::USE_BSSN) {
    theory_rhs.K += 4. * M_PI * m_G_Newton * theory_vars.lapse *
                    (Sij_TF_and_S.S + rho_and_Si.rho / (1. + b_of_x));
    theory_rhs.Theta += 0.0;
  } else {
    theory_rhs.K += 4.0 * M_PI * m_G_Newton * theory_vars.lapse *
                    (Sij_TF_and_S.S - 3 * rho_and_Si.rho / (1. + b_of_x));
    theory_rhs.Theta += -8. * M_PI * m_G_Newton * theory_vars.lapse *
                        rho_and_Si.rho / (1. + b_of_x);
  }

  // Update RHS for other variables

  FOR(i, j) {
    theory_rhs.A[i][j] += -8. * M_PI * m_G_Newton * theory_vars.lapse *
                          theory_vars.chi * Sij_TF_and_S.Sij_TF[i][j];
  }

  FOR(i, j) {
    theory_rhs.Gamma[i] += -16. * M_PI * m_G_Newton * theory_vars.lapse *
                           h_UU[i][j] * rho_and_Si.Si[j] / (1. + b_of_x);
  }
}

#endif /* MODIFIEDCCZ4RHS_IMPL_HPP_ */
