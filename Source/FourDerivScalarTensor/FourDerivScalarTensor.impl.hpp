/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(FOURDERIVSCALARTENSOR_HPP_)
#error "This file should only be included through FourDerivScalarTensor.hpp"
#endif

#ifndef FOURDERIVSCALARTENSOR_IMPL_HPP_
#define FOURDERIVSCALARTENSOR_IMPL_HPP_

template <class coupling_and_potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
ScalarVectorTensor<data_t>
FourDerivScalarTensor<coupling_and_potential_t>::compute_M_Ni_and_Mij(
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

template <class coupling_and_potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
ScalarVectorTensor<data_t>
FourDerivScalarTensor<coupling_and_potential_t>::compute_Omega_munu(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const Coordinates<data_t> &coords) const
{

    ScalarVectorTensor<data_t> out;

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

    // relevant quantities
    data_t chi_regularised = simd_max(1e-6, vars.chi);
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
        out.tensor[i][j] = 4. * dfdphi *
                               (covd2phi[i][j] +
                                vars.Pi / chi_regularised *
                                    (vars.A[i][j] + vars.h[i][j] * vars.K /
                                                        (double)GR_SPACEDIM)) +
                           4. * d2fdphi2 * d1.phi[i] * d1.phi[j];
    }
    // trace of Omega_ij
    out.scalar = vars.chi * compute_trace(out.tensor, h_UU);
    // Omega_i = -\gamma^{\mu}_{~i}n^{\nu}\Omega_{\mu\nu}
    FOR(i)
    {
        out.vector[i] =
            -4. * d2fdphi2 * vars.Pi * d1.phi[i] +
            4. * dfdphi *
                (-d1.Pi[i] - vars.K * d1.phi[i] / (double)GR_SPACEDIM);
        FOR(j, k)
        {
            out.vector[i] +=
                -4. * dfdphi * h_UU[j][k] * d1.phi[k] * vars.A[i][j];
        }
    }
    return out;
}

// Calculate rho and Si
template <class coupling_and_potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
RhoAndSi<data_t>
FourDerivScalarTensor<coupling_and_potential_t>::compute_rho_and_Si(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2,
    const Coordinates<data_t> &coords) const
{
    RhoAndSi<data_t> out;

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

    // S_i (note lower index) = - n^a T_ai
    FOR(i) { out.Si[i] = -d1.phi[i] * vars.Pi + g2 * Vt * vars.Pi * d1.phi[i]; }

    // rho = n^a n^b T_ab
    out.rho = vars.Pi * vars.Pi + 0.5 * Vt + V_of_phi -
              g2 * Vt * (Vt / 4. + vars.Pi * vars.Pi);

    // Compute useful quantities for the Gauss-Bonnet sector
    const data_t chi_regularised = simd_max(1e-6, vars.chi);

    ScalarVectorTensor<data_t> SVT = compute_M_Ni_and_Mij(vars, d1, d2);
    data_t M = SVT.scalar;
    Tensor<1, data_t> Ni = SVT.vector;
    Tensor<2, data_t> Mij = SVT.tensor;

    // decomposition of Omega_{\mu\nu}
    SVT = compute_Omega_munu(vars, d1, d2, coords);
    data_t Omega = SVT.scalar;
    Tensor<1, data_t> Omega_i = SVT.vector;
    Tensor<2, data_t> Omega_ij = SVT.tensor;

    Tensor<2, data_t> Omega_ij_UU =
        raise_all(Omega_ij, h_UU); // raise all indexs
    FOR(i, j) Omega_ij_UU[i][j] *= vars.chi * vars.chi;

    // Gauss-Bonnet contribution to rho
    data_t rhoGB = Omega * M;
    FOR(i, j) rhoGB -= 2. * Mij[i][j] * Omega_ij_UU[i][j];

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

    // Gauss-Bonnet contribution to Si
    Tensor<1, data_t> JGB;
    FOR(i)
    {
        JGB[i] = Omega_i[i] * M + 2. * Omega * (Ni[i] + d1.K[i] / 3.);
        FOR(j, k)
        {
            JGB[i] -= 2 * h_UU[j][k] * vars.chi *
                      (Mij[i][j] * Omega_i[k] +
                       Omega_ij[i][j] * (Ni[k] + d1.K[k] / 3.));
            FOR(l, m)
            {
                JGB[i] += 2 * h_UU[j][l] * h_UU[k][m] * Omega_ij[l][m] *
                          vars.chi *
                          (covd_Aphys_times_chi[j][k][i] -
                           covd_Aphys_times_chi[i][k][j]);
            }
        }
    }
    out.rho += rhoGB;
    FOR(i) out.Si[i] += JGB[i];

    return out;
}

// Calculate the stress energy tensor elements
template <class coupling_and_potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
SijTFAndS<data_t>
FourDerivScalarTensor<coupling_and_potential_t>::compute_Sij_TF_and_S(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec,
    const Coordinates<data_t> &coords) const
{

    SijTFAndS<data_t> out;

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

    const data_t chi_regularised = simd_max(vars.chi, 1e-6);

    // Useful quantity Vt
    data_t Vt = -vars.Pi * vars.Pi;
    FOR(i, j) { Vt += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j]; }

    // Calculate components of EM Tensor for the non-Gauss-Bonnet sector
    // S_ij = T_ij
    FOR(i, j)
    {
        out.Sij_TF[i][j] =
            -0.5 * vars.h[i][j] * Vt / chi_regularised + d1.phi[i] * d1.phi[j];
    }

    FOR(i, j)
    {
        out.Sij_TF[i][j] += -vars.h[i][j] * V_of_phi / chi_regularised;
    }

    FOR(i, j)
    {
        out.Sij_TF[i][j] +=
            g2 * Vt *
            (-d1.phi[i] * d1.phi[j] + vars.h[i][j] / chi_regularised * Vt / 4.);
    }

    // S = Tr_S_ij
    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij_TF, h_UU);

    make_trace_free(out.Sij_TF, vars.h, h_UU); // make Sij trace-free

    // Compute useful quantities for the Gauss-Bonnet sector
    ScalarVectorTensor<data_t> SVT = compute_M_Ni_and_Mij(vars, d1, d2);
    data_t M = SVT.scalar;
    Tensor<1, data_t> Ni = SVT.vector;
    Tensor<2, data_t> Mij = SVT.tensor;

    Tensor<2, data_t> covdtilde2phi;
    Tensor<2, data_t> covd2phi_times_chi;
    data_t dphi_dot_dchi = compute_dot_product(d1.phi, d1.chi, h_UU);
    FOR(k, l)
    {
        covdtilde2phi[k][l] = d2.phi[k][l];
        FOR1(m) { covdtilde2phi[k][l] -= chris.ULL[m][k][l] * d1.phi[m]; }
        covd2phi_times_chi[k][l] =
            vars.chi * covdtilde2phi[k][l] +
            0.5 * (d1.phi[k] * d1.chi[l] + d1.chi[k] * d1.phi[l] -
                   vars.h[k][l] * dphi_dot_dchi);
    }

    // decomposition of Omega_{\mu\nu}
    SVT = compute_Omega_munu(vars, d1, d2, coords);
    data_t Omega = SVT.scalar;
    Tensor<1, data_t> Omega_i = SVT.vector;
    Tensor<2, data_t> Omega_ij = SVT.tensor;

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

    // F_{ij} = \chi{\mathcal L}_nAphys_{ij}/\alpha + \chi D_iD_j\alpha/\alpha
    //+ A_{ik}A^k_{~j}) - \partial_tA_{ij}/\alpha

    data_t one_over_lapse = 1. / simd_max(1e-6, vars.lapse);
    Tensor<2, data_t> Fij;
    FOR(i, j)
    {
        Fij[i][j] =
            (-advec.A[i][j] + covd2lapse_times_chi[i][j]) * one_over_lapse -
            2. / 3. * vars.A[i][j] * (vars.K - divshift * one_over_lapse);
        FOR(k)
        {
            FOR(l)
            Fij[i][j] += h_UU[k][l] * vars.A[i][k] * vars.A[l][j];
            Fij[i][j] -= (vars.A[k][i] * d1.shift[k][j] +
                          vars.A[k][j] * d1.shift[k][i]) *
                         one_over_lapse;
        }
    }

    Tensor<2, data_t> Fij_TF = Fij;
    make_trace_free(Fij_TF, vars.h, h_UU);

    // F = {\mathcal L}_nK/\alpha + D^i_D_i\alpha/\alpha - K_{ij}K^{ij}) -
    // \partial_tK/\alpha
    data_t F = one_over_lapse * (-advec.K + tr_covd2lapse) - tr_A2 -
               vars.K * vars.K / 3.;

    Tensor<2, data_t> Mij_TF = Mij;
    make_trace_free(Mij_TF, vars.h, h_UU);
    Tensor<2, data_t> Mij_TF_UU_over_chi =
        raise_all(Mij_TF, h_UU); // raise all indexs
    FOR(i, j) Mij_TF_UU_over_chi[i][j] *= vars.chi;

    data_t RGB = -4. / 3. * M * F; // RHS terms of the Gauss-Bonnet curvature
    FOR(i, j)
    {
        RGB += 8. * Mij_TF_UU_over_chi[i][j] * Fij[i][j] +
               16. / 3. * vars.chi * h_UU[i][j] * d1.K[i] *
                   (Ni[j] + 1. / 3. * d1.K[j]);
        FOR(k, l, m, n)
        RGB -= 8. * vars.chi * h_UU[i][l] * h_UU[j][m] * h_UU[k][n] *
               covd_Aphys_times_chi[m][n][l] *
               (covd_Aphys_times_chi[j][k][i] - covd_Aphys_times_chi[i][j][k]);
    }

    Tensor<2, data_t> Omega_ij_UU =
        raise_all(Omega_ij, h_UU); // raise all indexs
    FOR(i, j) Omega_ij_UU[i][j] *= vars.chi * vars.chi;
    data_t rhoGB = Omega * M;
    FOR(i, j) rhoGB -= 2. * Mij[i][j] * Omega_ij_UU[i][j];

    // terms depending on g2 and V(phi) coming from having inserted the equation
    // for Pi in Sij and S
    data_t quadratic_terms = -dVdphi - 3. / 4. * dg2dphi * Vt * Vt -
                             2. * g2 * vars.Pi * vars.K * vars.Pi * vars.Pi;
    FOR(i, j)
    {
        // includes non conformal parts of chris not included in chris_ULL
        quadratic_terms +=
            2. * g2 * vars.Pi * vars.Pi * h_UU[i][j] *
            (0.5 * d1.chi[j] * d1.phi[i] - vars.chi * d2.phi[i][j]);
        FOR(k)
        {
            quadratic_terms += 2. * g2 * vars.Pi * vars.Pi * vars.chi *
                               h_UU[i][j] * chris.ULL[k][i][j] * d1.phi[k];
        }
    }
    quadratic_terms +=
        2. / 3. * g2 * vars.Pi * vars.K * (Vt + vars.Pi * vars.Pi);
    FOR(i, j)
    {
        quadratic_terms +=
            4. * g2 * h_UU[i][j] * vars.chi * d1.phi[j] * vars.Pi * d1.Pi[i];
        FOR(k, l)
        quadratic_terms += 2. * g2 * h_UU[i][k] * h_UU[j][l] * vars.chi *
                           d1.phi[k] * d1.phi[l] *
                           (vars.Pi * vars.A[i][j] - covd2phi_times_chi[i][j]);
    }
    quadratic_terms *= dfdphi / (1. + g2 * (-Vt + 2. * vars.Pi * vars.Pi));

    Tensor<2, data_t> Omega_ij_TF = Omega_ij;
    make_trace_free(Omega_ij_TF, vars.h, h_UU);
    Tensor<2, data_t> Omega_ij_TF_UU_over_chi2 =
        raise_all(Omega_ij_TF, h_UU); // raise all indexs

    // Gauss-Bonnet contribution for S
    data_t SGB = 4. / 3. * Omega * F +
                 4. * M * (-d2fdphi2 * Vt + quadratic_terms + Omega / 3.) -
                 rhoGB;
    FOR(i, j)
    SGB += -2. * Omega_ij_TF_UU_over_chi2[i][j] * vars.chi *
               (vars.chi * Mij_TF[i][j] + Fij[i][j]) -
           4. * h_UU[i][j] * vars.chi * Ni[i] * Omega_i[j];
    // add quadratic terms
    SGB += 4. * dfdphi * dfdphi * M * RGB /
           (1. + g2 * (-Vt + 2. * vars.Pi * vars.Pi));

    // Gauss-Bonnet contribution for SijTF
    Tensor<2, data_t> SijGB;
    FOR(i, j)
    {
        SijGB[i][j] = -2. / 3. * Omega_ij_TF[i][j] *
                          (F + 2. * (tr_covd2lapse * one_over_lapse - tr_A2)) -
                      2. * Mij_TF[i][j] *
                          (Omega - 4. * d2fdphi2 * Vt + 4. * quadratic_terms) -
                      2. * Omega / 3. * Fij_TF[i][j] / chi_regularised +
                      2. * ((Ni[i] + d1.K[i] / 3.) * Omega_i[j] +
                            (Ni[j] + d1.K[j] / 3.) * Omega_i[i]);
        FOR(k, l)
        {
            SijGB[i][j] +=
                2. * h_UU[k][l] *
                    (Omega_ij_TF[i][k] * Fij[l][j] +
                     Omega_ij_TF[j][k] * Fij[l][i] -
                     Omega_i[l] * (2. * covd_Aphys_times_chi[i][j][k] -
                                   covd_Aphys_times_chi[j][k][i] -
                                   covd_Aphys_times_chi[k][i][j])) -
                4. / 3. * vars.h[i][j] *
                    (Omega_ij_TF_UU_over_chi2[k][l] * Fij[k][l] +
                     h_UU[k][l] * Omega_i[k] * (2. * Ni[l] + d1.K[l]));
        }
        // add quadratic terms
        SijGB[i][j] += -8. * dfdphi * dfdphi * Mij_TF[i][j] * RGB /
                       (1. + g2 * (-Vt + 2. * vars.Pi * vars.Pi));
    }

    out.S += SGB;
    FOR(i, j) out.Sij_TF[i][j] += SijGB[i][j];

    return out;
}

// Adds in the RHS for the theory vars
template <class coupling_and_potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void FourDerivScalarTensor<coupling_and_potential_t>::add_theory_rhs(
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
                (Ni[j] + 1. / 3. * d1.K[j]);
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

// Computes the LHS matrix
template <class coupling_and_potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
void FourDerivScalarTensor<coupling_and_potential_t>::compute_lhs(
    const int N, data_t *LHS, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec,
    const Coordinates<data_t> &coords) const
{

    data_t LHS_mat[N][N];

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

    // Compute useful quantities for the Gauss-Bonnet sector

    ScalarVectorTensor<data_t> SVT = compute_M_Ni_and_Mij(vars, d1, d2);
    data_t M = SVT.scalar;
    Tensor<2, data_t> Mij = SVT.tensor;

    // decomposition of Omega_{\mu\nu}
    SVT = compute_Omega_munu(vars, d1, d2, coords);
    data_t Omega = SVT.scalar;
    Tensor<1, data_t> Omega_i = SVT.vector;
    Tensor<2, data_t> Omega_ij = SVT.tensor;

    Tensor<2, data_t> Mij_TF = Mij;
    make_trace_free(Mij_TF, vars.h, h_UU);
    Tensor<2, data_t> Mij_TF_UU_over_chi =
        raise_all(Mij_TF, h_UU); // raise all indexs
    FOR(i, j) Mij_TF_UU_over_chi[i][j] *= vars.chi;

    Tensor<2, data_t> Omega_ij_TF = Omega_ij;
    make_trace_free(Omega_ij_TF, vars.h, h_UU);
    Tensor<2, data_t> Omega_ij_TF_UU_over_chi =
        raise_all(Omega_ij_TF, h_UU); // raise all indexs
    FOR(i, j) Omega_ij_TF_UU_over_chi[i][j] *= vars.chi;

    data_t dfdphi2 =
        dfdphi * dfdphi / (1. + g2 * (-Vt + 2. * vars.Pi * vars.Pi));
    // comes from inserting the equation for Pi in Sij and S

    double G_factor = 16. * M_PI * m_G_Newton;
    int row = 0;
    FOR(i1, j1)
    {
        if (i1 > j1) // lower diagonal
            continue;

        int col = 0;
        FOR(i2, j2)
        {
            if (i2 > j2) // lower diagonal
                continue;
            // Lapack LHS argument has row-by-row layout
            // 'row' is the 'line' (multiplying variables 'x' in A.x = b), 'col'
            // is the column (which equation)
            LHS_mat[row][col] =
                -2. * G_factor *
                (1. / 3. * vars.h[i2][j2] * Omega_ij_TF_UU_over_chi[i1][j1] +
                 16. * vars.chi * dfdphi2 * Mij_TF[i2][j2] *
                     Mij_TF_UU_over_chi[i1][j1]);
            if (i1 == i2)
            {
                if (j1 == j2)
                    LHS_mat[row][col] += 1. - G_factor * Omega / 3.;
                FOR(k)
                LHS_mat[row][col] +=
                    G_factor * vars.chi * h_UU[j1][k] * Omega_ij_TF[j2][k];
            }
            if (j1 == j2)
                FOR(k)
            LHS_mat[row][col] +=
                G_factor * vars.chi * h_UU[i1][k] * Omega_ij_TF[i2][k];
            if (i1 != j1)
            {
                LHS_mat[row][col] +=
                    -2. * G_factor *
                    (1. / 3. * vars.h[i2][j2] *
                         Omega_ij_TF_UU_over_chi[j1][i1] +
                     16. * vars.chi * dfdphi2 * Mij_TF[i2][j2] *
                         Mij_TF_UU_over_chi[j1][i1]);
                if (i1 == j2)
                {
                    if (i2 == j1)
                        LHS_mat[row][col] += 1. - G_factor * Omega / 3.;
                    FOR(k)
                    LHS_mat[row][col] +=
                        G_factor * vars.chi * h_UU[j1][k] * Omega_ij_TF[i2][k];
                }
                if (i2 == j1)
                    FOR(k)
                LHS_mat[row][col] +=
                    G_factor * vars.chi * h_UU[i1][k] * Omega_ij_TF[j2][k];
            }

            ++col;
        }

        ++row;
    }

    int idx = 0;
    FOR(i, j)
    {
        if (i > j) // lower diagonal
            continue;

        LHS_mat[N - 2][idx] =
            G_factor * vars.chi / 3. *
            (-Omega_ij_TF[i][j] + 16. * dfdphi2 * M * Mij_TF[i][j]);

        LHS_mat[idx][N - 2] = 0.5 * G_factor *
                              (Omega_ij_TF_UU_over_chi[i][j] -
                               16. * dfdphi2 * M * Mij_TF_UU_over_chi[i][j]);
        if (i != j)
            LHS_mat[idx][N - 2] +=
                0.5 * G_factor *
                (Omega_ij_TF_UU_over_chi[j][i] -
                 16. * dfdphi2 * M * Mij_TF_UU_over_chi[j][i]);
        ++idx;
    }

    LHS_mat[N - 2][N - 2] =
        1. + G_factor / 3. * (-Omega + 4. * dfdphi2 * M * M);

    idx = 0;
    FOR(i, j)
    {
        if (i > j) // lower diagonal
            continue;

        LHS_mat[N - 1][idx] = 0.;
        LHS_mat[idx][N - 1] =
            -2. * G_factor * dfdphi * Mij_TF_UU_over_chi[i][j];
        if (i != j)
            LHS_mat[idx][N - 1] +=
                -2. * G_factor * dfdphi * Mij_TF_UU_over_chi[j][i];
        ++idx;
    }
    LHS_mat[N - 1][N - 2] = 0.;
    LHS_mat[N - 2][N - 1] = G_factor / 3. * dfdphi * M;
    LHS_mat[N - 1][N - 1] = 1. + G_factor * g2 * (2. * vars.Pi * vars.Pi - Vt);

    for (int row = 0; row < N; ++row)
    {
        for (int col = 0; col < N; ++col)
        {
            LHS[col * N + row] = LHS_mat[row][col];
        }
    }
}

// Function to solve LHS from the rhs and the LHS matrix
template <class coupling_and_potential_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void FourDerivScalarTensor<coupling_and_potential_t>::solve_lhs(
    rhs_vars_t<data_t> &rhs, const vars_t<data_t> &vars,
    const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec,
    const Coordinates<data_t> &coords) const
{

    const int N = GR_SPACEDIM * (GR_SPACEDIM + 1) / 2 + 2;
    data_t LHS[N][N];

    compute_lhs(N, (&LHS[0][0]), vars, d1, d2, advec, coords);

    data_t RHS[N];

    int row = 0;
    FOR(i1, j1)
    {
        if (i1 > j1) // lower diagonal
            continue;

        RHS[row] = rhs.A[i1][j1];
        ++row;
    }
    RHS[N - 2] = rhs.K;
    RHS[N - 1] = rhs.Pi;

    solve_linear_system(N, (&LHS[0][0]), RHS);

    row = 0;
    FOR(i1, j1)
    {
        if (i1 > j1) // lower diagonal
            continue;

        rhs.A[i1][j1] = RHS[row];
        if (i1 != j1)
            rhs.A[j1][i1] = RHS[row];

        ++row;
    }
    rhs.K = RHS[N - 2];
    rhs.Pi = RHS[N - 1];
}

#endif /* FOURDERIVSCALARTENSOR_IMPL_HPP_ */
