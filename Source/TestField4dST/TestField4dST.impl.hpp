/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(TestField4dST_HPP_)
#error "This file should only be included through TestField4dST.hpp"
#endif

#ifndef TESTFIELD4DST_IMPL_HPP_
#define TESTFIELD4DST_IMPL_HPP_

template <class coupling_and_potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
ScalarVectorTensor<data_t>
TestField4dST<coupling_and_potential_t>::compute_M_Ni_and_Mij(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2) const
{

    ScalarVectorTensor<data_t> out;

    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);
    const auto ricci0 =
        CCZ4Geometry::compute_ricci_Z(vars, d1, d2, h_UU, chris, {0., 0., 0.});
    const data_t chi_regularised = simd_max(1e-6, vars.chi);

    // M_{ij} = R_{ij} + KK_{ij} - K_{ik}K_j^{~k}
    FOR(i, j)
    {
        out.tensor[i][j] = ricci0.LL[i][j] +
                           vars.K / (3. * chi_regularised) *
                               (vars.A[i][j] + 2. / 3. * vars.K * vars.h[i][j]);
        FOR2(k, l)
        {
            out.tensor[i][j] -=
                vars.A[i][k] * vars.A[j][l] * h_UU[k][l] / chi_regularised;
        }
    }
    // M = \gamma^{ij}M_{ij} (GR Hamiltonian constraint)
    out.scalar = vars.chi * compute_trace(out.tensor, h_UU);

    Tensor<3, data_t> covdtilde_A;
    FOR(i, j, k)
    {
        covdtilde_A[j][k][i] = d1.A[j][k][i];
        FOR(l)
        {
            covdtilde_A[j][k][i] += -chris.ULL[l][i][j] * vars.A[l][k] -
                                    chris.ULL[l][i][k] * vars.A[l][j];
        }
    }
    // N_i = D^jK_{ij} - D_iK (GR momentum constraint)
    FOR(i) out.vector[i] = -(GR_SPACEDIM - 1.) * d1.K[i] / (double)GR_SPACEDIM;
    FOR(i, j, k)
    {
        out.vector[i] += h_UU[j][k] * (covdtilde_A[j][i][k] -
                                       GR_SPACEDIM * vars.A[i][j] * d1.chi[k] /
                                           (2. * chi_regularised));
    }
    return out;
}

// Calculate rho and Si
template <class coupling_and_potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
RhoAndSi<data_t> TestField4dST<coupling_and_potential_t>::compute_rho_and_Si(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const Coordinates<data_t> &coords) const
{
    RhoAndSi<data_t> out;
    FOR(i) { out.Si[i] = 0.; }
    out.rho = 0.;

    return out;
}

// Calculate the stress energy tensor elements
template <class coupling_and_potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
SijTFAndS<data_t> TestField4dST<coupling_and_potential_t>::compute_Sij_TF_and_S(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec,
    const Coordinates<data_t> &coords) const
{

    SijTFAndS<data_t> out;
    FOR(i, j) { out.Sij_TF[i][j] = 0.; }
    out.S = 0.;
    return out;
}

// Adds in the RHS for the theory vars
template <class coupling_and_potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void TestField4dST<coupling_and_potential_t>::add_theory_rhs(
    rhs_vars_t<data_t> &rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec,
    const Coordinates<data_t> &coords) const
{

    // set the coupling and potential values
    data_t dfdphi = 0.;
    data_t d2fdphi2 = 0.;
    data_t g2 = 0.;
    data_t dg2dphi = 0.;
    data_t V_of_phi = 0.;
    data_t dVdphi = 0.;

    // compute coupling and potential
    my_coupling_and_potential.compute_coupling_and_potential(
        dfdphi, d2fdphi2, g2, dg2dphi, V_of_phi, dVdphi, vars, coords);

    using namespace TensorAlgebra;

    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // Useful quantity Vt
    data_t Vt = -vars.Pi * vars.Pi;
    FOR(i, j) { Vt += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j]; }

    // evolution equations for scalar field and (minus) its conjugate momentum
    rhs.phi = vars.lapse * vars.Pi + advec.phi;
    rhs.Pi = vars.lapse * vars.K * vars.Pi + advec.Pi;

    FOR(i, j)
    {
        // includes non conformal parts of chris not included in chris_ULL
        rhs.Pi += -h_UU[i][j] * (0.5 * d1.chi[j] * vars.lapse * d1.phi[i] -
                                 vars.chi * vars.lapse * d2.phi[i][j] -
                                 vars.chi * d1.lapse[i] * d1.phi[j]);
        FOR(k)
        {
            rhs.Pi += -vars.chi * vars.lapse * h_UU[i][j] * chris.ULL[k][i][j] *
                      d1.phi[k];
        }
    }

    // Compute useful quantities for the Gauss-Bonnet sector

    const data_t chi_regularised = simd_max(1e-6, vars.chi);

    ScalarVectorTensor<data_t> SVT = compute_M_Ni_and_Mij(vars, d1, d2);
    data_t M = SVT.scalar;
    Tensor<1, data_t> Ni = SVT.vector;
    Tensor<2, data_t> Mij = SVT.tensor;

    data_t divshift = compute_trace(d1.shift);
    data_t dlapse_dot_dchi = compute_dot_product(d1.lapse, d1.chi, h_UU);

    Tensor<2, data_t> covdtilde2lapse;
    Tensor<2, data_t> covd2lapse_times_chi;
    FOR(k, l)
    {
        covdtilde2lapse[k][l] = d2.lapse[k][l];
        FOR(m) covdtilde2lapse[k][l] -= chris.ULL[m][k][l] * d1.lapse[m];
        covd2lapse_times_chi[k][l] =
            vars.chi * covdtilde2lapse[k][l] +
            0.5 * (d1.lapse[k] * d1.chi[l] + d1.chi[k] * d1.lapse[l] -
                   vars.h[k][l] * dlapse_dot_dchi);
    }
    data_t tr_covd2lapse = -(GR_SPACEDIM / 2.0) * dlapse_dot_dchi;
    FOR1(i)
    {
        tr_covd2lapse -= vars.chi * chris.contracted[i] * d1.lapse[i];
        FOR1(j)
        {
            tr_covd2lapse += h_UU[i][j] * (vars.chi * d2.lapse[i][j] +
                                           d1.lapse[i] * d1.chi[j]);
        }
    }

    Tensor<2, data_t> A_UU = raise_all(vars.A, h_UU);

    // A^{ij} A_{ij}. - Note the abuse of the compute trace function.
    data_t tr_A2 = compute_trace(vars.A, A_UU);

    // F_{ij} = \chi{\mathcal L}_nAphys_{ij} + \chi D_iD_j\alpha
    //+ \alphaA_{ik}A^k_{~j}) - \partial_tA_{ij}

    Tensor<2, data_t> Fij_times_lapse;
    FOR(i, j)
    {
        Fij_times_lapse[i][j] =
            -advec.A[i][j] + covd2lapse_times_chi[i][j] -
            2. / 3. * vars.A[i][j] * (vars.lapse * vars.K - divshift);
        FOR(k)
        {
            Fij_times_lapse[i][j] -=
                (vars.A[k][i] * d1.shift[k][j] + vars.A[k][j] * d1.shift[k][i]);
            FOR(l)
            Fij_times_lapse[i][j] +=
                vars.lapse * h_UU[k][l] * vars.A[i][k] * vars.A[l][j];
        }
    }

    // F_times_lapse = {\mathcal L}_nK + D^i_D_i\alpha - \alphaK_{ij}K^{ij} -
    // \partial_tK
    data_t F_times_lapse =
        -advec.K + tr_covd2lapse - vars.lapse * (tr_A2 + vars.K * vars.K / 3.);

    // other useful quantities
    Tensor<3, data_t> covdtilde_A;
    Tensor<3, data_t> covd_Aphys_times_chi;
    FOR(i, j, k)
    {
        covdtilde_A[j][k][i] = d1.A[j][k][i];
        FOR(l)
        {
            covdtilde_A[j][k][i] += -chris.ULL[l][i][j] * vars.A[l][k] -
                                    chris.ULL[l][i][k] * vars.A[l][j];
        }
        covd_Aphys_times_chi[j][k][i] =
            covdtilde_A[j][k][i] +
            (vars.A[i][k] * d1.chi[j] + vars.A[i][j] * d1.chi[k]) /
                (2. * chi_regularised);
        FOR(l, m)
        {
            covd_Aphys_times_chi[j][k][i] -=
                h_UU[l][m] * d1.chi[m] / (2. * chi_regularised) *
                (vars.h[i][j] * vars.A[k][l] + vars.h[i][k] * vars.A[j][l]);
        }
    }

    Tensor<2, data_t> Mij_TF = Mij;
    make_trace_free(Mij_TF, vars.h, h_UU);
    Tensor<2, data_t> Mij_TF_UU_over_chi =
        raise_all(Mij_TF, h_UU); // raise all indexs
    FOR(i, j) Mij_TF_UU_over_chi[i][j] *= vars.chi;

    // rhs of the Gauss-Bonnet curvature (multiplied by the lapse)
    data_t RGB_times_lapse = -4. / 3. * M * F_times_lapse;
    FOR(i, j)
    {
        RGB_times_lapse +=
            8 * Mij_TF_UU_over_chi[i][j] * Fij_times_lapse[i][j] +
            16. / 3. * vars.chi * vars.lapse * h_UU[i][j] * d1.K[i] *
                (Ni[j] + 1. / 3. * d1.K[j]) +
            8. * vars.chi * vars.lapse * h_UU[i][j] * Ni[i] * Ni[j];
        FOR(k, l, m, n)
        RGB_times_lapse -=
            8. * vars.chi * vars.lapse * h_UU[i][l] * h_UU[j][m] * h_UU[k][n] *
            covd_Aphys_times_chi[m][n][l] *
            (covd_Aphys_times_chi[j][k][i] - covd_Aphys_times_chi[i][j][k]);
    }
    rhs.Pi += dfdphi * RGB_times_lapse;
    rhs.Pi += -vars.lapse * dVdphi;

    // g2 contribution
    rhs.Pi += -3. / 4. * vars.lapse * dg2dphi * Vt * Vt -
              vars.lapse * vars.K * vars.Pi * g2 * Vt +
              advec.Pi * g2 * (2. * vars.Pi * vars.Pi - Vt);

    FOR(i, j)
    {
        // includes non conformal parts of chris not included in chris_ULL
        rhs.Pi += g2 * Vt * h_UU[i][j] *
                  (0.5 * d1.chi[j] * vars.lapse * d1.phi[i] -
                   vars.chi * vars.lapse * d2.phi[i][j]);
        FOR(k)
        {
            rhs.Pi += g2 * Vt * vars.chi * vars.lapse * h_UU[i][j] *
                      chris.ULL[k][i][j] * d1.phi[k];
        }
        rhs.Pi += g2 * h_UU[i][j] * vars.chi * d1.lapse[i] * d1.phi[j] *
                  (2. * vars.Pi * vars.Pi - Vt);
    }

    data_t dphi_dot_dchi = compute_dot_product(d1.phi, d1.chi, h_UU);
    Tensor<2, data_t> covdtilde2phi;
    Tensor<2, data_t> covd2phi_times_chi;
    FOR(k, l)
    {
        covdtilde2phi[k][l] = d2.phi[k][l];
        FOR(m) covdtilde2phi[k][l] -= chris.ULL[m][k][l] * d1.phi[m];
        covd2phi_times_chi[k][l] =
            vars.chi * covdtilde2phi[k][l] +
            0.5 * (d1.phi[k] * d1.chi[l] + d1.chi[k] * d1.phi[l] -
                   vars.h[k][l] * dphi_dot_dchi);
    }

    rhs.Pi +=
        2. / 3. * g2 * vars.lapse * vars.Pi * vars.K * (Vt + vars.Pi * vars.Pi);
    FOR(i, j)
    {
        rhs.Pi += 4. * g2 * h_UU[i][j] * vars.chi * d1.phi[j] * vars.lapse *
                  vars.Pi * d1.Pi[i];
        FOR(k, l)
        rhs.Pi += 2. * g2 * vars.lapse * h_UU[i][k] * h_UU[j][l] * vars.chi *
                  d1.phi[k] * d1.phi[l] *
                  (vars.Pi * vars.A[i][j] - covd2phi_times_chi[i][j]);
    }
}

// Function to solve LHS from the rhs
//! here we do not need to solve a linear system as gravity is not back-reacted
//! by GB terms one only needs to solve rhs of Pi from rhs of Aij and K
template <class coupling_and_potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void TestField4dST<coupling_and_potential_t>::solve_lhs(
    rhs_vars_t<data_t> &rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec,
    const Coordinates<data_t> &coords) const
{

    // set the coupling and potential values
    data_t dfdphi = 0.;
    data_t d2fdphi2 = 0.;
    data_t g2 = 0.;
    data_t dg2dphi = 0.;
    data_t V_of_phi = 0.;
    data_t dVdphi = 0.;

    // compute coupling and potential
    my_coupling_and_potential.compute_coupling_and_potential(
        dfdphi, d2fdphi2, g2, dg2dphi, V_of_phi, dVdphi, vars, coords);

    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);

    // Useful quantity Vt
    data_t Vt = -vars.Pi * vars.Pi;
    FOR(i, j) { Vt += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j]; }

    // Compute useful quantities for the Gauss-Bonnet sector

    ScalarVectorTensor<data_t> SVT = compute_M_Ni_and_Mij(vars, d1, d2);
    data_t M = SVT.scalar;
    Tensor<2, data_t> Mij = SVT.tensor;

    Tensor<2, data_t> Mij_TF = Mij;
    make_trace_free(Mij_TF, vars.h, h_UU);
    Tensor<2, data_t> Mij_TF_UU_over_chi =
        raise_all(Mij_TF, h_UU); // raise all indexs
    FOR(i, j) Mij_TF_UU_over_chi[i][j] *= vars.chi;

    // the lhs coefficient of  dtPi

    data_t lhs_dtPi = 1. + g2 * (2. * vars.Pi * vars.Pi - Vt);

    // the contribution of moving dtA_ij and dtK terms from lhs to rhs
    rhs.Pi -= rhs.K / 3. * dfdphi * M;

    FOR(i, j)
    {
        rhs.Pi -= -2. * dfdphi * Mij_TF_UU_over_chi[i][j] * rhs.A[i][j];
    }

    // solve the simple linear system
    rhs.Pi /= lhs_dtPi;
}

// Function to compute all the components of rho (used as diagnostics)
template <class coupling_and_potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
AllRhos<data_t> TestField4dST<coupling_and_potential_t>::compute_all_rhos(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const Coordinates<data_t> &coords) const
{
    AllRhos<data_t> out;

    // set the coupling and potential values
    data_t dfdphi = 0.;
    data_t d2fdphi2 = 0.;
    data_t g2 = 0.;
    data_t dg2dphi = 0.;
    data_t V_of_phi = 0.;
    data_t dVdphi = 0.;

    // compute coupling and potential
    my_coupling_and_potential.compute_coupling_and_potential(
        dfdphi, d2fdphi2, g2, dg2dphi, V_of_phi, dVdphi, vars, coords);

    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // Useful quantity Vt
    data_t Vt = -vars.Pi * vars.Pi;
    FOR(i, j) { Vt += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j]; }

    // rho = n^a n^b T_ab
    out.phi = vars.Pi * vars.Pi + 0.5 * Vt + V_of_phi;
    out.g2 = -g2 * Vt * (Vt / 4. + vars.Pi * vars.Pi);

    // Compute useful quantities for the Gauss-Bonnet sector
    const data_t chi_regularised = simd_max(1e-6, vars.chi);

    ScalarVectorTensor<data_t> SVT = compute_M_Ni_and_Mij(vars, d1, d2);
    data_t M = SVT.scalar;
    Tensor<2, data_t> Mij = SVT.tensor;

    // decomposition of Omega_{\mu\nu}
    Tensor<2, data_t> Omega_ij;

    Tensor<2, data_t> covdtilde2phi;
    Tensor<2, data_t> covd2phi;
    data_t dphi_dot_dchi = compute_dot_product(d1.phi, d1.chi, h_UU);
    FOR(k, l)
    {
        covdtilde2phi[k][l] = d2.phi[k][l];
        FOR1(m) { covdtilde2phi[k][l] -= chris.ULL[m][k][l] * d1.phi[m]; }
        covd2phi[k][l] = covdtilde2phi[k][l] +
                         0.5 *
                             (d1.phi[k] * d1.chi[l] + d1.chi[k] * d1.phi[l] -
                              vars.h[k][l] * dphi_dot_dchi) /
                             chi_regularised;
    }

    // Omega_{ij}=\gamma^{\mu}_{~i}\gamma^{\nu}_{~j}\Omega_{\mu\nu}
    FOR(i, j)
    {
        Omega_ij[i][j] = 4. * dfdphi *
                             (covd2phi[i][j] +
                              vars.Pi / chi_regularised *
                                  (vars.A[i][j] + vars.h[i][j] * vars.K /
                                                      (double)GR_SPACEDIM)) +
                         4. * d2fdphi2 * d1.phi[i] * d1.phi[j];
    }
    // trace of Omega_ij
    data_t Omega = vars.chi * compute_trace(Omega_ij, h_UU);

    Tensor<2, data_t> Omega_ij_UU =
        raise_all(Omega_ij, h_UU); // raise all indexs
    FOR(i, j) Omega_ij_UU[i][j] *= vars.chi * vars.chi;

    // Gauss-Bonnet contribution to rho
    out.GB = Omega * M;
    FOR(i, j) out.GB -= 2. * Mij[i][j] * Omega_ij_UU[i][j];

    out.g3 = 0.;

    return out;
}

#endif /* TESTFIELD4DST_IMPL_HPP_ */
