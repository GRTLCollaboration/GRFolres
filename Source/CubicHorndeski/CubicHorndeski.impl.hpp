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

    // Useful quantity Vt
    data_t Vt = -vars.Pi * vars.Pi;
    FOR(i, j) { Vt += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j]; }

    // S_i (note lower index) = - n^a T_ai
    FOR(i) { out.Si[i] = -d1.phi[i] * vars.Pi; }

    // rho = n^a n^b T_ab
    out.rho = vars.Pi * vars.Pi + 0.5 * Vt;

    // Horndeski contribution
    out.rho += E.dg3_dX * E.exprA + E.dg3_dphi * 2. * E.Xplus +
               E.dg2_dX * E.Pi2 - E.g2 + E.Xplus;

    FOR(i)
    {
        out.Si[i] +=
            E.dg3_dX * (-E.tau * vars.Pi * d1.phi[i] - E.Pi2 * E.taui[i] +
                        d1.phi[i] * E.taui_dot_dphi * vars.chi +
                        vars.Pi * E.tauij_dot_dphi[i]) -
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

    // compute potential and add constributions to EM Tensor

    const data_t chi_regularised = simd_max(vars.chi, 1e-6);

    // Useful quantity Vt
    data_t Vt = -vars.Pi * vars.Pi;
    FOR(i, j) { Vt += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j]; }

    // Calculate components of EM Tensor for the non-Horndeski sector
    // S_ij = T_ij
    FOR(i, j)
    {
        out.Sij_TF[i][j] =
            -0.5 * vars.h[i][j] * Vt / chi_regularised + d1.phi[i] * d1.phi[j];
    }

    FOR(i, j)
    {
        out.Sij_TF[i][j] =
            E.dg3_dX *
                (E.tau * d1.phi[i] * d1.phi[j] +
                 vars.Pi * (d1.phi[i] * E.taui[j] + d1.phi[j] * E.taui[i]) -
                 (d1.phi[i] * E.tauij_dot_dphi[j] +
                  d1.phi[j] * E.tauij_dot_dphi[i]) -
                 vars.h[i][j] * E.exprB +
                 E.lieD_Pi_no_lapse * (-d1.phi[i] * d1.phi[j] +
                                       vars.h[i][j] / vars.chi * E.Pi2)) +
            d1.phi[i] * d1.phi[j] * E.dcommon +
            vars.h[i][j] / vars.chi * (E.X + E.g2 + 2. * E.X * E.dg3_dphi);
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

    data_t lieD_Pi = E.lieD_Pi_no_lapse;
    FOR2(i, j)
    {
        lieD_Pi += h_UU[i][j] * vars.chi * d1.lapse[i] / vars.lapse * d1.phi[j];
    }

    // adjust RHS for the potential term
    total_rhs.phi = advec.phi + vars.lapse * vars.Pi;
    total_rhs.Pi = advec.Pi + vars.lapse * lieD_Pi;
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

    // lieD_Pi_potential = NUM/DEN + A^k D_k lapse (not conformal)
    data_t NUM = 0.;
    data_t DEN = 0.;

    ////////////////////////////////////////////////////////////////////////
    // DEFINE VARIABLES USED REPETEDLY

    data_t chi = vars.chi;
    data_t Pi = vars.Pi;

    E.Pi2 = Pi * Pi;

    // dphi_dot_dphi - not conformal
    data_t dphi_dot_dphi = 0.;
    FOR2(i, j) { dphi_dot_dphi += h_UU[i][j] * d1.phi[i] * d1.phi[j]; }
    dphi_dot_dphi *= chi;

    // X
    E.X = 0.5 * (E.Pi2 - dphi_dot_dphi);
    E.Xplus = 0.5 * (E.Pi2 + dphi_dot_dphi);

    ////////////////////////////////////////////////////////////////////////
    // Functions G2 & G3 and derivatives
    E.g2 = this->my_coupling_and_potential.G2(vars.phi, E.X);
    // E.g3             = this->my_potential.G3(vars.phi, X);

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
    // tau auxialiary variables

    // covd2phi
    Tensor<2, data_t> covd2phi;
    FOR2(i, j)
    {
        covd2phi[i][j] = d2.phi[i][j];
        FOR1(k) { covd2phi[i][j] -= chris_ULL[k][i][j] * d1.phi[k]; }
    }

    // dphi_dot_dchi
    data_t dphi_dot_dchi = 0.;
    FOR2(i, j) { dphi_dot_dchi += h_UU[i][j] * d1.phi[i] * d1.chi[j]; }

    // tau_ij
    // Tensor<2, data_t> tauij;
    FOR2(i, j)
    {
        E.tauij[i][j] = vars.A[i][j] * Pi + vars.K * Pi * vars.h[i][j] / 3. +
                        0.5 * (-vars.h[i][j] * dphi_dot_dchi +
                               d1.phi[i] * d1.chi[j] + d1.phi[j] * d1.chi[i] +
                               chi * (covd2phi[i][j] + covd2phi[j][i]));
    }

    // tau
    E.tau = compute_trace(E.tauij, h_UU);

    // taui
    FOR1(i)
    {
        E.taui[i] = vars.K * d1.phi[i] / 3. + d1.Pi[i];
        FOR2(j, k) { E.taui[i] += h_UU[j][k] * vars.A[i][j] * d1.phi[k]; }
    }

    ////////////////////////////////////////////////////////////////////////
    // common expressions

    FOR1(i)
    {
        E.tauij_dot_dphi[i] = 0.;
        FOR2(j, k)
        {
            E.tauij_dot_dphi[i] += h_UU[j][k] * d1.phi[k] * E.tauij[i][j];
        }
    }

    E.tauij_dot_dphi2 = 0.;
    FOR2(i, j)
    {
        E.tauij_dot_dphi2 += h_UU[i][j] * d1.phi[i] * E.tauij_dot_dphi[j];
    }

    E.exprA = E.tau * E.Pi2 - E.tauij_dot_dphi2 * chi;

    E.taui_dot_dphi = 0.;
    FOR2(i, j) { E.taui_dot_dphi += h_UU[i][j] * d1.phi[i] * E.taui[j]; }

    E.exprB = 2. * Pi * E.taui_dot_dphi - E.tauij_dot_dphi2;

    ///////////////////////////////////////////////////////
    // DEN

    E.ders1 =
        E.dcommon + 2. * E.tau * E.dg3_dX - E.X * E.X * E.dg3_dX * E.dg3_dX -
        E.tauij_dot_dphi2 * vars.chi * E.d2g3_dXX - 2. * E.X * E.d2g3_dXphi;
    E.ders2 = 2. * E.X * E.dg3_dX * E.dg3_dX + E.d2g2_dXX + E.tau * E.d2g3_dXX +
              2. * E.d2g3_dXphi;

    DEN = E.ders1 + E.Pi2 * E.ders2;

    ///////////////////////////////////////////////////////
    // NUM

    NUM = E.tau * (E.dcommon - 2. * E.X * E.d2g3_dXphi) -
          E.dg3_dX * E.dg3_dX * (E.tau * E.X - 2. * chi * E.exprB) * E.X +
          (E.d2g2_dXX + 2. * E.d2g3_dXphi) * E.exprB * chi -
          E.d2g3_dXX * chi *
              (-E.tau * E.exprB + chi * E.taui_dot_dphi * E.taui_dot_dphi) +
          E.dg2_dphi - 2. * E.X * (E.d2g3_dphiphi + E.d2g2_dXphi) -
          E.dg3_dX * (-E.tau * E.tau + E.X * E.g2 +
                      E.X * E.X * (2. + E.dg2_dX + 4. * E.dg3_dphi));

    FOR2(i, j)
    {
        NUM += E.d2g3_dXX * h_UU[i][j] *
                   (Pi * E.taui[i] - E.tauij_dot_dphi[i]) *
                   (Pi * E.taui[j] - E.tauij_dot_dphi[j]) +
               E.dg3_dX * h_UU[i][j] * E.taui[i] * E.taui[j] * 2. * chi;
        FOR2(k, l)
        {
            NUM -= E.dg3_dX * h_UU[i][k] * h_UU[j][l] * E.tauij[i][j] *
                   E.tauij[k][l];
        }
    }

    ///////////////////////////////////////////////////////
    // Lie Derivative term: lieD_Pi - A^k D_k ln(lapse)

    E.lieD_Pi_no_lapse = NUM / DEN;

    E.DEN = DEN; // for debugging
}

#endif /* CUBICHORNDESKI_IMPL_HPP_ */
