/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CUBICHORNDESKI_HPP_
#error "This file should only be included through CubicHorndeski.hpp"
#endif

#ifndef CUBICHORNDESKI_IMPL_HPP_
#define CUBICHORNDESKI_IMPL_HPP_
#include "DimensionDefinitions.hpp"

template <class coupling_and_potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
RhoAndSi<data_t> CubicHorndeski<coupling_and_potential_t>::compute_rho_and_Si(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const Coordinates<data_t> &coords) const
{
    UsefulQuantities<data_t> quantities;

    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    compute_useful_quantities(quantities, vars, d1, d2, h_UU, chris.ULL);

    RhoAndSi<data_t> out;

    // dphi_dot_dphi - not conformal
    data_t dphi_dot_dphi = 0.;
    FOR(i, j)
    {
        dphi_dot_dphi += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j];
    }

    // kinetic contribution to the energy density
    data_t Xplus = 0.5 * (vars.Pi * vars.Pi + dphi_dot_dphi);

    // rho = n^a n^b T_ab
    out.rho = quantities.dg3_dX * (quantities.tau * vars.Pi * vars.Pi -
                                   quantities.tau_ij_dot_dphi2 * vars.chi) +
              quantities.dg3_dphi * 2. * Xplus +
              quantities.dg2_dX * vars.Pi * vars.Pi - quantities.g2 + Xplus +
              quantities.V;

    // S_i (note lower index) = - n^a T_ai
    FOR(i)
    {
        out.Si[i] = quantities.dg3_dX *
                        (-quantities.tau * vars.Pi * d1.phi[i] -
                         vars.Pi * vars.Pi * quantities.tau_i[i] +
                         d1.phi[i] * quantities.tau_i_dot_dphi * vars.chi +
                         vars.Pi * quantities.tau_ij_dot_dphi[i]) -
                    vars.Pi * d1.phi[i] *
                        (1. + quantities.dg2_dX + 2. * quantities.dg3_dphi);
    }

    return out;
}

template <class coupling_and_potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
SijTFAndS<data_t>
CubicHorndeski<coupling_and_potential_t>::compute_Sij_TF_and_S(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec,
    const Coordinates<data_t> &coords) const
{
    UsefulQuantities<data_t> quantities;
    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    compute_useful_quantities(quantities, vars, d1, d2, h_UU, chris.ULL);

    SijTFAndS<data_t> out;

    // dphi_dot_dphi - not conformal
    data_t dphi_dot_dphi = 0.;
    FOR(i, j)
    {
        dphi_dot_dphi += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j];
    }

    // kinetic term in the Lagrangian
    data_t X = 0.5 * (vars.Pi * vars.Pi - dphi_dot_dphi);

    // compute potential and add contributions to EM Tensor
    const data_t chi_regularised = simd_max(vars.chi, 1e-6);

    FOR(i, j)
    {
        out.Sij_TF[i][j] =
            quantities.dg3_dX *
                (quantities.tau * d1.phi[i] * d1.phi[j] +
                 vars.Pi * (d1.phi[i] * quantities.tau_i[j] +
                            d1.phi[j] * quantities.tau_i[i]) -
                 (d1.phi[i] * quantities.tau_ij_dot_dphi[j] +
                  d1.phi[j] * quantities.tau_ij_dot_dphi[i]) -
                 vars.h[i][j] * (2. * vars.Pi * quantities.tau_i_dot_dphi -
                                 quantities.tau_ij_dot_dphi2) +
                 quantities.lie_deriv_Pi_no_lapse *
                     (-d1.phi[i] * d1.phi[j] +
                      vars.h[i][j] / chi_regularised * vars.Pi * vars.Pi)) +
            d1.phi[i] * d1.phi[j] *
                (1. + quantities.dg2_dX + 2. * quantities.dg3_dphi) +
            vars.h[i][j] / chi_regularised *
                (X + quantities.g2 - quantities.V +
                 2. * X * quantities.dg3_dphi);
    }

    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij_TF, h_UU);

    make_trace_free(out.Sij_TF, vars.h, h_UU); // make Sij tracefree

    return out;
}

// Adds in the RHS for the matter vars
template <class coupling_and_potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void CubicHorndeski<coupling_and_potential_t>::add_theory_rhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec,
    const Coordinates<data_t> &coords) const
{
    UsefulQuantities<data_t> quantities;
    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    compute_useful_quantities(quantities, vars, d1, d2, h_UU, chris.ULL);

    data_t lie_deriv_Pi_times_lapse =
        vars.lapse * quantities.lie_deriv_Pi_no_lapse;
    FOR(i, j)
    {
        lie_deriv_Pi_times_lapse +=
            h_UU[i][j] * vars.chi * d1.lapse[i] * d1.phi[j];
    }

    // adjust RHS for the potential term
    total_rhs.phi = advec.phi + vars.lapse * vars.Pi;
    total_rhs.Pi = advec.Pi + lie_deriv_Pi_times_lapse;
}

template <class coupling_and_potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void CubicHorndeski<coupling_and_potential_t>::solve_lhs(
    rhs_vars_t<data_t> &total_rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec,
    const Coordinates<data_t> &coords) const
{
}

template <class coupling_and_potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
void CubicHorndeski<coupling_and_potential_t>::compute_useful_quantities(
    UsefulQuantities<data_t> &quantities, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const Tensor<2, data_t> &h_UU,
    const Tensor<3, data_t> &chris_ULL) const
{
    using namespace TensorAlgebra;

    // lie_deriv_Pi = numerator/denominator + A^k D_k lapse (not conformal)
    data_t numerator = 0.;
    data_t denominator = 0.;

    // dphi_dot_dphi - not conformal
    data_t dphi_dot_dphi = 0.;
    FOR(i, j)
    {
        dphi_dot_dphi += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j];
    }

    // kinetic term in the Lagrangian
    data_t X = 0.5 * (vars.Pi * vars.Pi - dphi_dot_dphi);

    // kinetic contribution to the energy density
    data_t Xplus = 0.5 * (vars.Pi * vars.Pi + dphi_dot_dphi);

    ////////////////////////////////////////////////////////////////////////
    // Functions G2 & G3 and derivatives

    quantities.V = this->my_coupling_and_potential.V(vars.phi, X);
    quantities.g2 = this->my_coupling_and_potential.G2(vars.phi, X);

    quantities.dV_dphi = this->my_coupling_and_potential.dV_dphi(vars.phi, X);
    quantities.dg2_dphi = this->my_coupling_and_potential.dG2_dphi(vars.phi, X);
    quantities.dg3_dphi = this->my_coupling_and_potential.dG3_dphi(vars.phi, X);
    quantities.dg2_dX = this->my_coupling_and_potential.dG2_dX(vars.phi, X);
    quantities.dg3_dX = this->my_coupling_and_potential.dG3_dX(vars.phi, X);

    quantities.d2g2_dXX = this->my_coupling_and_potential.d2G2_dXX(vars.phi, X);
    quantities.d2g2_dXphi =
        this->my_coupling_and_potential.d2G2_dXphi(vars.phi, X);
    quantities.d2g3_dXX = this->my_coupling_and_potential.d2G3_dXX(vars.phi, X);
    quantities.d2g3_dXphi =
        this->my_coupling_and_potential.d2G3_dXphi(vars.phi, X);
    quantities.d2g3_dphiphi =
        this->my_coupling_and_potential.d2G3_dphiphi(vars.phi, X);

    ////////////////////////////////////////////////////////////////////////
    // tau auxiliary variables

    // covd2phi
    Tensor<2, data_t> covd2phi;
    FOR(i, j)
    {
        covd2phi[i][j] = d2.phi[i][j];
        FOR(k) { covd2phi[i][j] -= chris_ULL[k][i][j] * d1.phi[k]; }
    }

    // dphi_dot_dchi
    data_t dphi_dot_dchi = 0.;
    FOR(i, j) { dphi_dot_dchi += h_UU[i][j] * d1.phi[i] * d1.chi[j]; }

    // tau_ij
    FOR(i, j)
    {
        quantities.tau_ij[i][j] =
            vars.A[i][j] * vars.Pi + vars.K * vars.Pi * vars.h[i][j] / 3. +
            0.5 * (-vars.h[i][j] * dphi_dot_dchi + d1.phi[i] * d1.chi[j] +
                   d1.phi[j] * d1.chi[i] +
                   vars.chi * (covd2phi[i][j] + covd2phi[j][i]));
    }

    // tau
    quantities.tau = compute_trace(quantities.tau_ij, h_UU);

    // tau_i
    FOR(i)
    {
        quantities.tau_i[i] = vars.K * d1.phi[i] / 3. + d1.Pi[i];
        FOR(j, k)
        {
            quantities.tau_i[i] += h_UU[j][k] * vars.A[i][j] * d1.phi[k];
        }
    }

    ////////////////////////////////////////////////////////////////////////
    // common expressions

    FOR(i)
    {
        quantities.tau_ij_dot_dphi[i] = 0.;
        FOR(j, k)
        {
            quantities.tau_ij_dot_dphi[i] +=
                h_UU[j][k] * d1.phi[k] * quantities.tau_ij[i][j];
        }
    }

    quantities.tau_ij_dot_dphi2 = 0.;
    FOR(i, j)
    {
        quantities.tau_ij_dot_dphi2 +=
            h_UU[i][j] * d1.phi[i] * quantities.tau_ij_dot_dphi[j];
    }

    quantities.tau_i_dot_dphi = 0.;
    FOR(i, j)
    {
        quantities.tau_i_dot_dphi +=
            h_UU[i][j] * d1.phi[i] * quantities.tau_i[j];
    }

    data_t another_useful_quantity =
        2. * vars.Pi * quantities.tau_i_dot_dphi - quantities.tau_ij_dot_dphi2;

    ///////////////////////////////////////////////////////
    // denominator

    denominator =
        1. + quantities.dg2_dX + 2. * quantities.dg3_dphi +
        2. * quantities.tau * quantities.dg3_dX -
        X * X * quantities.dg3_dX * quantities.dg3_dX -
        quantities.tau_ij_dot_dphi2 * vars.chi * quantities.d2g3_dXX -
        2. * X * quantities.d2g3_dXphi +
        vars.Pi * vars.Pi *
            (2. * X * quantities.dg3_dX * quantities.dg3_dX +
             quantities.d2g2_dXX + quantities.tau * quantities.d2g3_dXX +
             2. * quantities.d2g3_dXphi);

    ///////////////////////////////////////////////////////
    // numerator

    numerator =
        quantities.tau * (1. + quantities.dg2_dX + 2. * quantities.dg3_dphi -
                          2. * X * quantities.d2g3_dXphi) -
        quantities.dg3_dX * quantities.dg3_dX *
            (quantities.tau * X - 2. * vars.chi * another_useful_quantity) * X +
        (quantities.d2g2_dXX + 2. * quantities.d2g3_dXphi) *
            another_useful_quantity * vars.chi -
        quantities.d2g3_dXX * vars.chi *
            (-quantities.tau * another_useful_quantity +
             vars.chi * quantities.tau_i_dot_dphi * quantities.tau_i_dot_dphi) +
        quantities.dg2_dphi - quantities.dV_dphi -
        2. * X * (quantities.d2g3_dphiphi + quantities.d2g2_dXphi) -
        quantities.dg3_dX *
            (-quantities.tau * quantities.tau + X * quantities.g2 +
             X * X * (2. + quantities.dg2_dX + 4. * quantities.dg3_dphi));

    FOR(i, j)
    {
        numerator += quantities.d2g3_dXX * h_UU[i][j] * vars.chi *
                         (vars.Pi * quantities.tau_i[i] -
                          quantities.tau_ij_dot_dphi[i]) *
                         (vars.Pi * quantities.tau_i[j] -
                          quantities.tau_ij_dot_dphi[j]) +
                     quantities.dg3_dX * h_UU[i][j] * quantities.tau_i[i] *
                         quantities.tau_i[j] * 2. * vars.chi;
        FOR(k, l)
        {
            numerator += -quantities.dg3_dX * h_UU[i][k] * h_UU[j][l] *
                         quantities.tau_ij[i][j] * quantities.tau_ij[k][l];
        }
    }

    ///////////////////////////////////////////////////////
    // Lie Derivative term: lie_deriv_Pi - A^k D_k ln(lapse)

    quantities.lie_deriv_Pi_no_lapse = numerator / denominator;
}

// Function to compute all the components of rho (used as diagnostics)
template <class coupling_and_potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
AllRhos<data_t> CubicHorndeski<coupling_and_potential_t>::compute_all_rhos(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const Coordinates<data_t> &coords) const
{
    AllRhos<data_t> out;

    UsefulQuantities<data_t> quantities;

    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    compute_useful_quantities(quantities, vars, d1, d2, h_UU, chris.ULL);

    // dphi_dot_dphi - not conformal
    data_t dphi_dot_dphi = 0.;
    FOR(i, j)
    {
        dphi_dot_dphi += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j];
    }

    // kinetic contribution to the energy density
    data_t Xplus = 0.5 * (vars.Pi * vars.Pi + dphi_dot_dphi);

    // rho = n^a n^b T_ab
    out.phi = Xplus + quantities.V;

    // Horndeski contribution
    out.g2 = quantities.dg2_dX * vars.Pi * vars.Pi - quantities.g2;
    out.g3 = quantities.dg3_dX * (quantities.tau * vars.Pi * vars.Pi -
                                  quantities.tau_ij_dot_dphi2 * vars.chi) +
             quantities.dg3_dphi * 2. * Xplus;

    out.GB = 0.;

    return out;
}

#endif /* CUBICHORNDESKI_IMPL_HPP_ */
