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

template <class matter_t, class gauge_t, class deriv_t, class modified_gauge_t>
ModifiedCCZ4RHS<matter_t, gauge_t, deriv_t, modified_gauge_t>::ModifiedCCZ4RHS(
    matter_t a_matter, CCZ4_params_t<typename gauge_t::params_t> a_params,
    modified_gauge_t a_modified_gauge, double a_dx, double a_sigma,
    const std::array<double, CH_SPACEDIM> a_center, int a_formulation,
    double a_G_Newton)
    : CCZ4RHS<gauge_t, deriv_t>(a_params, a_dx, a_sigma, a_formulation,
                                0.0 /*No cosmological constant*/),
      my_matter(a_matter), my_modified_gauge(a_modified_gauge),
      m_center(a_center), m_G_Newton(a_G_Newton) {}

template <class matter_t, class gauge_t, class deriv_t, class modified_gauge_t>
template <class data_t>
void ModifiedCCZ4RHS<matter_t, gauge_t, deriv_t, modified_gauge_t>::compute(
    Cell<data_t> current_cell) const {
  // copy data from chombo gridpoint into local variables
  const auto matter_vars = current_cell.template load_vars<Vars>();
  const auto d1 = this->m_deriv.template diff1<Vars>(current_cell);
  const auto d2 = this->m_deriv.template diff2<Diff2Vars>(current_cell);
  const auto advec =
      this->m_deriv.template advection<Vars>(current_cell, matter_vars.shift);

  // Call CCZ4 RHS - work out RHS without matter, no dissipation
  Vars<data_t> matter_rhs;
  this->rhs_equation(matter_rhs, matter_vars, d1, d2, advec);

  Coordinates<data_t> coords{current_cell, this->m_deriv.m_dx, m_center};

  // add functions a(x) and b(x) of the modified gauge
  add_a_b_rhs(matter_rhs, matter_vars, d1, d2, advec, coords);

  // add RHS matter terms from EM Tensor
  add_emtensor_rhs(matter_rhs, matter_vars, d1, d2, advec, coords);

  // add evolution of matter fields themselves
  my_matter.add_matter_rhs(matter_rhs, matter_vars, d1, d2, advec, coords);

  // solve linear system for the matter fields that require it (e.g. 4dST)
  my_matter.solve_lhs(matter_rhs, matter_vars, d1, d2, advec, coords);

  // Add dissipation to all terms
  this->m_deriv.add_dissipation(matter_rhs, current_cell, this->m_sigma);

  // Write the rhs into the output FArrayBox
  current_cell.store_vars(matter_rhs);
}

// Function to add in EM Tensor matter terms to CCZ4 rhs
template <class matter_t, class gauge_t, class deriv_t, class modified_gauge_t>
template <class data_t>
void ModifiedCCZ4RHS<matter_t, gauge_t, deriv_t, modified_gauge_t>::add_a_b_rhs(
    Vars<data_t> &matter_rhs, const Vars<data_t> &matter_vars,
    const Vars<Tensor<1, data_t>> &d1, const Diff2Vars<Tensor<2, data_t>> &d2,
    const Vars<data_t> &advec, const Coordinates<data_t> &coords) const {
  data_t a_of_x = 0.;
  data_t b_of_x = 0.;

  my_modified_gauge.compute_modified_gauge(a_of_x, b_of_x, coords);

  const data_t chi_regularised = simd_max(1e-6, matter_vars.chi);
  using namespace TensorAlgebra;
  auto h_UU = compute_inverse_sym(matter_vars.h);
  auto chris = compute_christoffel(d1.h, h_UU);
  Tensor<1, data_t> Z_over_chi;
  Tensor<1, data_t> Z;
  if (this->m_formulation == CCZ4RHS<>::USE_BSSN) {
    FOR(i) Z_over_chi[i] = 0.0;
  } else {
    FOR(i)
    Z_over_chi[i] = 0.5 * (matter_vars.Gamma[i] - chris.contracted[i]);
  }
  auto ricci0 = CCZ4Geometry::compute_ricci_Z(matter_vars, d1, d2, h_UU, chris,
                                              {0., 0., 0.});
  // for some reason the function compute_ricci does not give the correct
  // answer

  Tensor<2, data_t> A_UU = raise_all(matter_vars.A, h_UU);
  // A^{ij} A_{ij}. - Note the abuse of the compute trace function.
  data_t tr_A2 = compute_trace(matter_vars.A, A_UU);

  Tensor<2, data_t> covdtilde_A[CH_SPACEDIM];

  FOR(i, j, k) {
    covdtilde_A[i][j][k] = d1.A[j][k][i];
    FOR(l) {
      covdtilde_A[i][j][k] += -chris.ULL[l][i][j] * matter_vars.A[l][k] -
                              chris.ULL[l][i][k] * matter_vars.A[l][j];
    }
  }

  Tensor<1, data_t> Ni;
  FOR(i) { Ni[i] = -(GR_SPACEDIM - 1.) * d1.K[i] / (double)GR_SPACEDIM; }
  FOR(i, j, k) {
    Ni[i] += h_UU[j][k] *
             (covdtilde_A[k][j][i] - GR_SPACEDIM * matter_vars.A[i][j] *
                                         d1.chi[k] / (2. * chi_regularised));
  }

  data_t kappa1_times_lapse;
  if (this->m_params.covariantZ4)
    kappa1_times_lapse = this->m_params.kappa1;
  else
    kappa1_times_lapse = this->m_params.kappa1 * matter_vars.lapse;

  if (this->m_formulation == CCZ4RHS<>::USE_BSSN) {
    matter_rhs.K += 0.;
  } else {
    matter_rhs.K +=
        b_of_x / (1. + b_of_x) *
        (matter_vars.lapse * (-0.5 * matter_vars.K * matter_vars.K +
                              0.5 * GR_SPACEDIM / (GR_SPACEDIM - 1.) *
                                  (tr_A2 - ricci0.scalar)) +
         kappa1_times_lapse * GR_SPACEDIM * matter_vars.Theta *
             (1. + this->m_params.kappa2 -
              0.5 * ((GR_SPACEDIM - 3.) / (GR_SPACEDIM - 1.) +
                     this->m_params.kappa2)));

    matter_rhs.Theta +=
        b_of_x / (1. + b_of_x) *
        (0.5 * matter_vars.lapse *
             (tr_A2 - ricci0.scalar -
              ((GR_SPACEDIM - 1.0) / (double)GR_SPACEDIM) * matter_vars.K *
                  matter_vars.K) +
         0.5 * matter_vars.Theta * kappa1_times_lapse *
             ((GR_SPACEDIM + 1.) + this->m_params.kappa2 * (GR_SPACEDIM - 1.)));
  }

  FOR(i) {
    data_t b_term_gamma =
        b_of_x / (1. + b_of_x) *
        ((2.0 / (double)GR_SPACEDIM) *
             (matter_vars.lapse * matter_vars.K * Z_over_chi[i]) +
         2. * kappa1_times_lapse * Z_over_chi[i]);
    FOR(j) {
      b_term_gamma += b_of_x / (1. + b_of_x) * 2. * h_UU[i][j] *
                      matter_vars.lapse * (-d1.Theta[j] - Ni[j]);
      FOR(k) {
        b_term_gamma += b_of_x / (1. + b_of_x) * 2. * matter_vars.lapse *
                        A_UU[i][j] * matter_vars.h[k][j] * Z_over_chi[k];
      }
    }
    matter_rhs.Gamma[i] += b_term_gamma;
    matter_rhs.B[i] += b_term_gamma;
  }

  matter_rhs.lapse += a_of_x / (1. + a_of_x) * this->m_params.lapse_coeff *
                      pow(matter_vars.lapse, this->m_params.lapse_power) *
                      (matter_vars.K - 2. * matter_vars.Theta);

  FOR(i) {
    matter_rhs.shift[i] -= a_of_x / (1. + a_of_x) *
                           this->m_params.shift_Gamma_coeff * matter_vars.B[i];
    FOR(j) {
      matter_rhs.shift[i] -= a_of_x / (1. + a_of_x) * matter_vars.lapse *
                             matter_vars.chi * h_UU[i][j] * d1.lapse[j];
    }
  }
}

// Function to add in EM Tensor matter terms to CCZ4 rhs
template <class matter_t, class gauge_t, class deriv_t, class modified_gauge_t>
template <class data_t>
void ModifiedCCZ4RHS<matter_t, gauge_t, deriv_t, modified_gauge_t>::
    add_emtensor_rhs(Vars<data_t> &matter_rhs, const Vars<data_t> &matter_vars,
                     const Vars<Tensor<1, data_t>> &d1,
                     const Diff2Vars<Tensor<2, data_t>> &d2,
                     const Vars<data_t> &advec,
                     const Coordinates<data_t> &coords) const {
  using namespace TensorAlgebra;

  const auto h_UU = compute_inverse_sym(matter_vars.h);
  const auto chris = compute_christoffel(d1.h, h_UU);

  // Calculate elements of the decomposed stress energy tensor
  const auto emtensor =
      my_matter.compute_emtensor(matter_vars, d1, d2, advec, coords);

  data_t a_of_x = 0.;
  data_t b_of_x = 0.;

  my_modified_gauge.compute_modified_gauge(a_of_x, b_of_x, coords);

  // Update RHS for K and Theta depending on formulation
  if (this->m_formulation == CCZ4RHS<>::USE_BSSN) {
    matter_rhs.K += 4. * M_PI * m_G_Newton * matter_vars.lapse *
                    (emtensor.S + emtensor.rho / (1. + b_of_x));
    matter_rhs.Theta += 0.0;
  } else {
    matter_rhs.K += 4.0 * M_PI * m_G_Newton * matter_vars.lapse *
                    (emtensor.S - 3 * emtensor.rho / (1. + b_of_x));
    matter_rhs.Theta += -8. * M_PI * m_G_Newton * matter_vars.lapse *
                        emtensor.rho / (1. + b_of_x);
  }

  // Update RHS for other variables
  Tensor<2, data_t> Sij_TF = emtensor.Sij;
  make_trace_free(Sij_TF, matter_vars.h, h_UU);

  FOR(i, j) {
    // matter_rhs.A[i][j] += -matter_vars.chi *
    //                      matter_vars.lapse * Sij_TF[i][j];
    matter_rhs.A[i][j] +=
        -8. * M_PI * m_G_Newton * matter_vars.lapse *
        (matter_vars.chi * emtensor.Sij[i][j] -
         matter_vars.h[i][j] * emtensor.S / (double)GR_SPACEDIM);
  }

  FOR(i) {
    data_t matter_term_Gamma = 0.0;
    FOR(j) {
      matter_term_Gamma += -16. * M_PI * m_G_Newton * matter_vars.lapse *
                           h_UU[i][j] * emtensor.Si[j] / (1. + b_of_x);
    }

    matter_rhs.Gamma[i] += matter_term_Gamma;
    matter_rhs.B[i] += matter_term_Gamma;
  }
}

#endif /* MODIFIEDCCZ4RHS_IMPL_HPP_ */
