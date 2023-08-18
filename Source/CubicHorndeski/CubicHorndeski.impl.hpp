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
    expressions<data_t> E;

    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    pre_compute_no_gauge(E, vars, d1, d2, h_UU, chris.ULL);

    RhoAndSi<data_t> out;

    // rho = n^a n^b T_ab
    out.rho = E.dg3_dX * E.exprA + E.dg3_dphi * 2. * E.Xplus +
              E.dg2_dX * E.Pi2 - E.g2 + E.Xplus + E.V;

    // S_i (note lower index) = - n^a T_ai
    FOR(i)
    {
        out.Si[i] =
            E.dg3_dX * (-E.tau * vars.Pi * d1.phi[i] - E.Pi2 * E.tau_i[i] +
                        d1.phi[i] * E.tau_i_dot_dphi * vars.chi +
                        vars.Pi * E.tau_ij_dot_dphi[i]) -
            vars.Pi * d1.phi[i] * E.dcommon;
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
    expressions<data_t> E;
    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    pre_compute_no_gauge(E, vars, d1, d2, h_UU, chris.ULL);

    SijTFAndS<data_t> out;

    // compute potential and add contributions to EM Tensor
    const data_t chi_regularised = simd_max(vars.chi, 1e-6);

    FOR(i, j)
    {
        out.Sij_TF[i][j] =
            E.dg3_dX *
                (E.tau * d1.phi[i] * d1.phi[j] +
                 vars.Pi * (d1.phi[i] * E.tau_i[j] + d1.phi[j] * E.tau_i[i]) -
                 (d1.phi[i] * E.tau_ij_dot_dphi[j] +
                  d1.phi[j] * E.tau_ij_dot_dphi[i]) -
                 vars.h[i][j] * E.exprB +
                 E.lie_deriv_Pi_no_lapse *
                     (-d1.phi[i] * d1.phi[j] +
                      vars.h[i][j] / chi_regularised * E.Pi2)) +
            d1.phi[i] * d1.phi[j] * E.dcommon +
            vars.h[i][j] / vars.chi *
                (E.X + E.g2 - E.V + 2. * E.X * E.dg3_dphi);
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
    expressions<data_t> E;
    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    pre_compute_no_gauge(E, vars, d1, d2, h_UU, chris.ULL);

    data_t lie_deriv_Pi_times_lapse = vars.lapse * E.lie_deriv_Pi_no_lapse;
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
void CubicHorndeski<coupling_and_potential_t>::pre_compute_no_gauge(
    expressions<data_t> &E, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const Tensor<2, data_t> &h_UU,
    const Tensor<3, data_t> &chris_ULL) const
{
    using namespace TensorAlgebra;

    // lie_deriv_Pi = numerator/denominator + A^k D_k lapse (not conformal)
    data_t numerator = 0.;
    data_t denominator = 0.;

    ////////////////////////////////////////////////////////////////////////
    // DEFINE VARIABLES USED REPETEDLY

    data_t chi = vars.chi;
    data_t Pi = vars.Pi;

    E.Pi2 = Pi * Pi;

    // dphi_dot_dphi - not conformal
    data_t dphi_dot_dphi = 0.;
    FOR(i, j) { dphi_dot_dphi += h_UU[i][j] * d1.phi[i] * d1.phi[j]; }
    dphi_dot_dphi *= chi;

    // X
    E.X = 0.5 * (E.Pi2 - dphi_dot_dphi);
    E.Xplus = 0.5 * (E.Pi2 + dphi_dot_dphi);

    ////////////////////////////////////////////////////////////////////////
    // Functions G2 & G3 and derivatives

    E.V = this->my_coupling_and_potential.V(vars.phi, E.X);
    E.g2 = this->my_coupling_and_potential.G2(vars.phi, E.X);

    E.dV_dphi = this->my_coupling_and_potential.dV_dphi(vars.phi, E.X);
    E.dg2_dphi = this->my_coupling_and_potential.dG2_dphi(vars.phi, E.X);
    E.dg3_dphi = this->my_coupling_and_potential.dG3_dphi(vars.phi, E.X);
    E.dg2_dX = this->my_coupling_and_potential.dG2_dX(vars.phi, E.X);
    E.dg3_dX = this->my_coupling_and_potential.dG3_dX(vars.phi, E.X);

    E.d2g2_dXX = this->my_coupling_and_potential.d2G2_dXX(vars.phi, E.X);
    E.d2g2_dXphi = this->my_coupling_and_potential.d2G2_dXphi(vars.phi, E.X);
    E.d2g3_dXX = this->my_coupling_and_potential.d2G3_dXX(vars.phi, E.X);
    E.d2g3_dXphi = this->my_coupling_and_potential.d2G3_dXphi(vars.phi, E.X);
    E.d2g3_dphiphi =
        this->my_coupling_and_potential.d2G3_dphiphi(vars.phi, E.X);

    E.dcommon = 1. + E.dg2_dX + 2. * E.dg3_dphi;

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
        E.tau_ij[i][j] = vars.A[i][j] * Pi + vars.K * Pi * vars.h[i][j] / 3. +
                         0.5 * (-vars.h[i][j] * dphi_dot_dchi +
                                d1.phi[i] * d1.chi[j] + d1.phi[j] * d1.chi[i] +
                                chi * (covd2phi[i][j] + covd2phi[j][i]));
    }

    // tau
    E.tau = compute_trace(E.tau_ij, h_UU);

    // tau_i
    FOR(i)
    {
        E.tau_i[i] = vars.K * d1.phi[i] / 3. + d1.Pi[i];
        FOR(j, k) { E.tau_i[i] += h_UU[j][k] * vars.A[i][j] * d1.phi[k]; }
    }

    ////////////////////////////////////////////////////////////////////////
    // common expressions

    FOR(i)
    {
        E.tau_ij_dot_dphi[i] = 0.;
        FOR(j, k)
        {
            E.tau_ij_dot_dphi[i] += h_UU[j][k] * d1.phi[k] * E.tau_ij[i][j];
        }
    }

    E.tau_ij_dot_dphi2 = 0.;
    FOR(i, j)
    {
        E.tau_ij_dot_dphi2 += h_UU[i][j] * d1.phi[i] * E.tau_ij_dot_dphi[j];
    }

    E.exprA = E.tau * E.Pi2 - E.tau_ij_dot_dphi2 * chi;

    E.tau_i_dot_dphi = 0.;
    FOR(i, j) { E.tau_i_dot_dphi += h_UU[i][j] * d1.phi[i] * E.tau_i[j]; }

    E.exprB = 2. * Pi * E.tau_i_dot_dphi - E.tau_ij_dot_dphi2;

    ///////////////////////////////////////////////////////
    // denominator

    E.ders1 =
        E.dcommon + 2. * E.tau * E.dg3_dX - E.X * E.X * E.dg3_dX * E.dg3_dX -
        E.tau_ij_dot_dphi2 * vars.chi * E.d2g3_dXX - 2. * E.X * E.d2g3_dXphi;
    E.ders2 = 2. * E.X * E.dg3_dX * E.dg3_dX + E.d2g2_dXX + E.tau * E.d2g3_dXX +
              2. * E.d2g3_dXphi;

    denominator = E.ders1 + E.Pi2 * E.ders2;

    ///////////////////////////////////////////////////////
    // numerator

    numerator =
        E.tau * (E.dcommon - 2. * E.X * E.d2g3_dXphi) -
        E.dg3_dX * E.dg3_dX * (E.tau * E.X - 2. * chi * E.exprB) * E.X +
        (E.d2g2_dXX + 2. * E.d2g3_dXphi) * E.exprB * chi -
        E.d2g3_dXX * chi *
            (-E.tau * E.exprB + chi * E.tau_i_dot_dphi * E.tau_i_dot_dphi) +
        E.dg2_dphi - E.dV_dphi - 2. * E.X * (E.d2g3_dphiphi + E.d2g2_dXphi) -
        E.dg3_dX * (-E.tau * E.tau + E.X * E.g2 +
                    E.X * E.X * (2. + E.dg2_dX + 4. * E.dg3_dphi));

    FOR(i, j)
    {
        numerator += E.d2g3_dXX * h_UU[i][j] *
                         (Pi * E.tau_i[i] - E.tau_ij_dot_dphi[i]) *
                         (Pi * E.tau_i[j] - E.tau_ij_dot_dphi[j]) +
                     E.dg3_dX * h_UU[i][j] * E.tau_i[i] * E.tau_i[j] * 2. * chi;
        FOR(k, l)
        {
            numerator += -E.dg3_dX * h_UU[i][k] * h_UU[j][l] * E.tau_ij[i][j] *
                         E.tau_ij[k][l];
        }
    }

    ///////////////////////////////////////////////////////
    // Lie Derivative term: lie_deriv_Pi - A^k D_k ln(lapse)

    E.lie_deriv_Pi_no_lapse = numerator / denominator;
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

    expressions<data_t> E;

    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    pre_compute_no_gauge(E, vars, d1, d2, h_UU, chris.ULL);

    // rho = n^a n^b T_ab
    out.phi = E.Xplus + E.V;

    // Horndeski contribution
    out.g2 = E.dg2_dX * E.Pi2 - E.g2;
    out.g3 = E.dg3_dX * E.exprA + E.dg3_dphi * 2. * E.Xplus;

    out.GB = 0.;

    return out;
}

#endif /* CUBICHORNDESKI_IMPL_HPP_ */
