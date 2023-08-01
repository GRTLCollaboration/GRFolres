/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(ESGB_HPP_)
#error "This file should only be included through EsGB.hpp"
#endif

#ifndef ESGB_IMPL_HPP_
#define ESGB_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class coupling_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
rho_Si_t<data_t> EsGB<coupling_t>::compute_rho_Si(
    const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
    const diff2_vars_t<Tensor<2, data_t>> &d2) const
{
    rho_Si_t<data_t> out;

    // set the potential values
    data_t dfGB = 0.0;
    data_t df2GB = 0.0;

    // compute coupling
    my_coupling.compute_coupling(dfGB, df2GB, vars);

    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // Useful quantity Vt
    data_t Vt = -vars.Kphi * vars.Kphi;
    FOR(i, j) { Vt += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j]; }

    // S_i (note lower index) = - n^a T_ai
    FOR(i) { out.Si[i] = d1.phi[i] * vars.Kphi; }

    // rho = n^a n^b T_ab
    out.rho = vars.Kphi * vars.Kphi + 0.5 * Vt;

    // Compute useful quantities for the Gauss-Bonnet sector
    const auto ricci0 =
        CCZ4Geometry::compute_ricci_Z(vars, d1, d2, h_UU, chris, {0., 0., 0.});
    // for some reason the function compute_ricci does not give the correct
    // answer

    const data_t chi_regularised = simd_max(1e-6, vars.chi);

    Tensor<2, data_t> Mij; // M_{ij} = R_{ij} + KK_{ij} - K_{ik}K_j^{~k}
    FOR(i, j)
    {
        Mij[i][j] = ricci0.LL[i][j] +
                    vars.K / (3. * chi_regularised) *
                        (vars.A[i][j] + 2. / 3. * vars.K * vars.h[i][j]);
        FOR2(k, l)
        {
            Mij[i][j] -=
                vars.A[i][k] * vars.A[j][l] * h_UU[k][l] / chi_regularised;
        }
    }
    data_t M =
        vars.chi *
        compute_trace(Mij, h_UU); // trace of M_{ij} (which is the Einstein Ham)

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

    // decomposition of C_{ab} = \nabla_a\nabla_bf(\phi)
    Tensor<2, data_t> Cij; // C_{ij} = \gamma^a_{~i}\gamma^b_{~j}C_{ab}
    FOR(i, j)
    {
        Cij[i][j] = 4. * dfGB * (covd2phi[i][j] -
                            vars.Kphi / chi_regularised *
                                (vars.A[i][j] +
                                 vars.h[i][j] * vars.K / (double)GR_SPACEDIM)) +
                    4. * df2GB * d1.phi[i] * d1.phi[j];
    }
    data_t C = vars.chi * compute_trace(Cij, h_UU); // trace of Cij

    Tensor<1, data_t> Ci; // C_i = -\gamma^a_{~i}n^bC_{ab}
    FOR(i)
    {
        Ci[i] = 4. * df2GB * vars.Kphi * d1.phi[i] +
                4. * dfGB * (d1.Kphi[i] - vars.K * d1.phi[i] / (double)GR_SPACEDIM);
        FOR(j, k) { Ci[i] += - 4. * dfGB * h_UU[j][k] * d1.phi[k] * vars.A[i][j]; }
    }

    Tensor<2, data_t> Cij_UU = raise_all(Cij, h_UU); // raise all indexs
    // note that they have been raised with the conformal metric

    data_t rhoGB = C * M; // rho = n^a n^b T_ab
    FOR(i, j) rhoGB -= 2. * vars.chi * Mij[i][j] * vars.chi * Cij_UU[i][j];

    // other useful quantities
    Tensor<2, data_t> covdtilde_A[CH_SPACEDIM];
    Tensor<2, data_t> covd_Aphys_times_chi[CH_SPACEDIM];

    FOR(i, j, k)
    {
        covdtilde_A[i][j][k] = d1.A[j][k][i];
        FOR(l)
        {
            covdtilde_A[i][j][k] += -chris.ULL[l][i][j] * vars.A[l][k] -
                                    chris.ULL[l][i][k] * vars.A[l][j];
        }
        covd_Aphys_times_chi[i][j][k] =
            covdtilde_A[i][j][k] +
            (vars.A[i][k] * d1.chi[j] + vars.A[i][j] * d1.chi[k]) /
                (2. * chi_regularised);
        FOR(l, m)
        {
            covd_Aphys_times_chi[i][j][k] -=
                h_UU[l][m] * d1.chi[m] / (2. * chi_regularised) *
                (vars.h[i][j] * vars.A[k][l] + vars.h[i][k] * vars.A[j][l]);
        }
    }

    Tensor<1, data_t>
        Ni; // N_i = D^jK_{ij} - D_iK (which is the Einstein Momentum)
    FOR(i) Ni[i] = -(GR_SPACEDIM - 1.) * d1.K[i] / (double)GR_SPACEDIM;
    FOR(i, j, k)
    {
        Ni[i] += h_UU[j][k] * (covdtilde_A[k][j][i] -
                               GR_SPACEDIM * vars.A[i][j] * d1.chi[k] /
                                   (2. * chi_regularised));
    }

    Tensor<1, data_t> JGB; // S_i (note lower index) = - n^a T_ai
    FOR(i)
    {
        JGB[i] = Ci[i] * M + 2. * C * (Ni[i] + d1.K[i] / 3.);
        FOR(j, k)
        {
            JGB[i] -= 2 * h_UU[j][k] * vars.chi *
                      (Mij[i][j] * Ci[k] + Cij[i][j] * (Ni[k] + d1.K[k] / 3.));
            FOR(l, m)
            {
                JGB[i] += 2 * h_UU[j][l] * h_UU[k][m] * Cij[l][m] * vars.chi *
                          (covd_Aphys_times_chi[i][j][k] -
                           covd_Aphys_times_chi[j][i][k]);
            }
        }
    }
    out.rho += rhoGB;
    FOR(i) out.Si[i] += JGB[i];

    return out;
}

// Calculate the stress energy tensor elements
template <class coupling_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
emtensor_t<data_t>
EsGB<coupling_t>::compute_emtensor(const vars_t<data_t> &vars,
                                   const vars_t<Tensor<1, data_t>> &d1,
                                   const diff2_vars_t<Tensor<2, data_t>> &d2,
                                   const vars_t<data_t> &advec) const
{
    emtensor_t<data_t> out;

    rho_Si_t<data_t> rho_Si = compute_rho_Si(vars, d1, d2);
    out.rho = rho_Si.rho;
    FOR(i) out.Si[i] = rho_Si.Si[i];

    // set the potential values
    data_t dfGB = 0.0;
    data_t df2GB = 0.0;

    // compute coupling
    my_coupling.compute_coupling(dfGB, df2GB, vars);

    using namespace TensorAlgebra;
    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // Useful quantity Vt
    data_t Vt = -vars.Kphi * vars.Kphi;
    FOR(i, j) { Vt += vars.chi * h_UU[i][j] * d1.phi[i] * d1.phi[j]; }

    // Calculate components of EM Tensor for the non-Gauss-Bonnet sector
    // S_ij = T_ij
    FOR(i, j)
    {
        out.Sij[i][j] =
            -0.5 * vars.h[i][j] * Vt / vars.chi + d1.phi[i] * d1.phi[j];
    }

    // S = Tr_S_ij
    out.S = vars.chi * TensorAlgebra::compute_trace(out.Sij, h_UU);

    // Compute useful quantities for the Gauss-Bonnet sector

    const auto ricci0 =
        CCZ4Geometry::compute_ricci_Z(vars, d1, d2, h_UU, chris, {0., 0., 0.});
    // for some reason the function compute_ricci does not give the correct
    // answer

    const data_t chi_regularised = simd_max(1e-6, vars.chi);

    Tensor<2, data_t> Mij; // M_{ij} = R_{ij} + KK_{ij} - K_{ik}K_j^{~k}
    FOR(i, j)
    {
        Mij[i][j] = ricci0.LL[i][j] +
                    vars.K / (3. * chi_regularised) *
                        (vars.A[i][j] + 2. / 3. * vars.K * vars.h[i][j]);
        FOR2(k, l)
        {
            Mij[i][j] -=
                vars.A[i][k] * vars.A[j][l] * h_UU[k][l] / chi_regularised;
        }
    }
    data_t M =
        vars.chi *
        compute_trace(Mij, h_UU); // trace of M_{ij} (which is the Einstein Ham)

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

    // decomposition of C_{ab} = \nabla_a\nabla_bf(\phi)
    Tensor<2, data_t> Cij; // C_{ij} = \gamma^a_{~i}\gamma^b_{~j}C_{ab}
    FOR(i, j)
    {
        Cij[i][j] = 4. * dfGB * (covd2phi[i][j] -
                            vars.Kphi / chi_regularised *
                                (vars.A[i][j] +
                                 vars.h[i][j] * vars.K / (double)GR_SPACEDIM)) +
                    4. * df2GB * d1.phi[i] * d1.phi[j];
    }
    data_t C = vars.chi * compute_trace(Cij, h_UU); // trace of Cij

    Tensor<1, data_t> Ci; // C_i = -\gamma^a_{~i}n^bC_{ab}
    FOR(i)
    {
        Ci[i] = 4. * df2GB * vars.Kphi * d1.phi[i] +
                4. * dfGB * (d1.Kphi[i] - vars.K * d1.phi[i] / (double)GR_SPACEDIM);
        FOR(j, k) { Ci[i] += - 4. * dfGB * h_UU[j][k] * d1.phi[k] * vars.A[i][j]; }
    }

    Tensor<2, data_t> Cij_UU = raise_all(Cij, h_UU); // raise all indexs
    // note that they have been raised with the conformal metric

    Tensor<2, data_t> Mij_UU = raise_all(Mij, h_UU); // raise all indexs
    // note that they have been raised with the conformal metric

    Tensor<2, data_t> Mij_TF = Mij;
    make_trace_free(Mij_TF, vars.h, h_UU);
    Tensor<2, data_t> Mij_TF_UU = raise_all(Mij_TF, h_UU); // raise all indexs
    // note that they have been raised with the conformal metric

    Tensor<2, data_t> Cij_TF = Cij;
    make_trace_free(Cij_TF, vars.h, h_UU);
    Tensor<2, data_t> Cij_TF_UU = raise_all(Cij_TF, h_UU); // raise all indexs
    // note that they have been raised with the conformal metric

    // other useful quantities
    Tensor<2, data_t> covdtilde_A[CH_SPACEDIM];
    Tensor<2, data_t> covd_Aphys_times_chi[CH_SPACEDIM];

    FOR(i, j, k)
    {
        covdtilde_A[i][j][k] = d1.A[j][k][i];
        FOR(l)
        {
            covdtilde_A[i][j][k] += -chris.ULL[l][i][j] * vars.A[l][k] -
                                    chris.ULL[l][i][k] * vars.A[l][j];
        }
        covd_Aphys_times_chi[i][j][k] =
            covdtilde_A[i][j][k] +
            (vars.A[i][k] * d1.chi[j] + vars.A[i][j] * d1.chi[k]) /
                (2. * chi_regularised);
        FOR(l, m)
        {
            covd_Aphys_times_chi[i][j][k] -=
                h_UU[l][m] * d1.chi[m] / (2. * chi_regularised) *
                (vars.h[i][j] * vars.A[k][l] + vars.h[i][k] * vars.A[j][l]);
        }
    }

    Tensor<1, data_t>
        Ni; // N_i = D^jK_{ij} - D_iK (which is the Einstein Momentum)
    FOR(i) Ni[i] = -(GR_SPACEDIM - 1.) * d1.K[i] / (double)GR_SPACEDIM;
    FOR(i, j, k)
    {
        Ni[i] += h_UU[j][k] * (covdtilde_A[k][j][i] -
                               GR_SPACEDIM * vars.A[i][j] * d1.chi[k] /
                                   (2. * chi_regularised));
    }

    data_t divshift = compute_trace(d1.shift);
    data_t dlapse_dot_dchi = compute_dot_product(d1.lapse, d1.chi, h_UU);

    Tensor<2, data_t> covdtilde2lapse;
    Tensor<2, data_t> covd2lapse;
    FOR(k, l)
    {
        covdtilde2lapse[k][l] = d2.lapse[k][l];
        FOR(m) covdtilde2lapse[k][l] -= chris.ULL[m][k][l] * d1.lapse[m];
        covd2lapse[k][l] =
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

    data_t lapse_regularised = simd_max(1e-4, vars.lapse);
    Tensor<2, data_t> Fij;
    FOR(i, j)
    {
        Fij[i][j] = -advec.A[i][j] -
                    2. / 3. * vars.A[i][j] * (vars.lapse * vars.K - divshift) +
                    covd2lapse[i][j];
        FOR(k)
        {
            FOR(l)
            Fij[i][j] += vars.lapse * h_UU[k][l] * vars.A[i][k] * vars.A[l][j];
            Fij[i][j] -=
                (vars.A[k][i] * d1.shift[k][j] + vars.A[k][j] * d1.shift[k][i]);
        }
        Fij[i][j] /= lapse_regularised;
    }

    Tensor<2, data_t> Fij_TF = Fij;
    make_trace_free(Fij_TF, vars.h, h_UU);

    // F = {\mathcal L}_nK/\alpha + D^i_D_i\alpha/\alpha - K_{ij}K^{ij}) -
    // \partial_tK/\alpha
    data_t F =
        -advec.K + tr_covd2lapse - vars.lapse * (tr_A2 + vars.K * vars.K / 3.);
    F /= lapse_regularised;

    data_t RGB = -4. / 3. * M * F; // RHS terms of the Gauss-Bonnet curvature
    FOR(i, j)
    {
        RGB += 8. * Mij_TF_UU[i][j] * Fij[i][j] * vars.chi +
               16. / 3. * vars.chi * h_UU[i][j] * d1.K[i] *
                   (Ni[j] + 1. / 3. * d1.K[j]);
        FOR(k, l, m, n)
        RGB -= 8. * vars.chi * h_UU[i][l] * h_UU[j][m] * h_UU[k][n] *
               covd_Aphys_times_chi[l][m][n] *
               (covd_Aphys_times_chi[i][j][k] - covd_Aphys_times_chi[k][i][j]);
    }

    data_t SGB = 4. / 3. * C * F + 4. * M * (- 4. * df2GB / 4. * Vt + C / 3.) -
                 (out.rho - vars.Kphi * vars.Kphi - 0.5 * Vt);
    FOR(i, j)
    SGB += - 2. * Cij_TF_UU[i][j] * vars.chi * (vars.chi * Mij_TF[i][j] + Fij[i][j]) -
           4. * h_UU[i][j] * vars.chi * Ni[i] * Ci[j];
    SGB += 4. * dfGB * 4. * dfGB / 4. * M * RGB;

    Tensor<2, data_t> SijGB;
    FOR(i, j)
    {
        SijGB[i][j] =
            -2. / 3. * Cij_TF[i][j] * vars.chi *
                (F + 2. * (tr_covd2lapse / lapse_regularised - tr_A2)) -
            2. * vars.chi * Mij_TF[i][j] * (C - 4. * df2GB * Vt) - 2. * C / 3. * Fij_TF[i][j] +
            2. * vars.chi * ((Ni[i] + d1.K[i] / 3.) * Ci[j] +
                        (Ni[j] + d1.K[j] / 3.) * Ci[i]);
        FOR(k, l)
        {
            SijGB[i][j] +=
                2. * vars.chi * h_UU[k][l] *
                    (Cij_TF[i][k] * Fij[l][j] + Cij_TF[j][k] * Fij[l][i] -
                     Ci[l] * (2. * covd_Aphys_times_chi[k][i][j] -
                              covd_Aphys_times_chi[i][j][k] -
                              covd_Aphys_times_chi[j][k][i])) -
                4. / 3. * vars.h[i][j] * vars.chi *
                    (Cij_TF_UU[k][l] * Fij[k][l] +
                     h_UU[k][l] * Ci[k] * (2. * Ni[l] + d1.K[l]));
        }
        SijGB[i][j] += vars.h[i][j] * SGB / 3.;
        SijGB[i][j] /= chi_regularised;
        SijGB[i][j] += -4. * dfGB * 4. * dfGB / 2. * Mij_TF[i][j] * RGB;
    }

    out.S += SGB;
    FOR(i, j) out.Sij[i][j] += SijGB[i][j];

    return out;
}

// Adds in the RHS for the matter vars
template <class coupling_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void EsGB<coupling_t>::add_matter_rhs(rhs_vars_t<data_t> &rhs,
                                      const vars_t<data_t> &vars,
                                      const vars_t<Tensor<1, data_t>> &d1,
                                      const diff2_vars_t<Tensor<2, data_t>> &d2,
                                      const vars_t<data_t> &advec) const
{
    // first get the non potential part of the rhs
    // this may seem a bit long winded, but it makes the function
    // work for more multiple fields

    // set the potential values
    data_t dfGB = 0.0;
    data_t df2GB = 0.0;
    my_coupling.compute_coupling(dfGB, df2GB, vars);

    using namespace TensorAlgebra;

    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // evolution equations for scalar field and (minus) its conjugate momentum
    rhs.phi = -vars.lapse * vars.Kphi + advec.phi;
    rhs.Kphi = vars.lapse * vars.K * vars.Kphi + advec.Kphi;

    FOR(i, j)
    {
        // includes non conformal parts of chris not included in chris_ULL
        rhs.Kphi += h_UU[i][j] * (0.5 * d1.chi[j] * vars.lapse * d1.phi[i] -
                                  vars.chi * vars.lapse * d2.phi[i][j] -
                                  vars.chi * d1.lapse[i] * d1.phi[j]);
        FOR(k)
        {
            rhs.Kphi += vars.chi * vars.lapse * h_UU[i][j] *
                        chris.ULL[k][i][j] * d1.phi[k];
        }
    }

    // Compute useful quantities for the Gauss-Bonnet sector

    const auto ricci0 =
        CCZ4Geometry::compute_ricci_Z(vars, d1, d2, h_UU, chris, {0., 0., 0.});
    // for some reason the function compute_ricci does not give the correct
    // answer

    const data_t chi_regularised = simd_max(1e-6, vars.chi);

    Tensor<2, data_t> Mij; // M_{ij} = R_{ij} + KK_{ij} - K_{ik}K_j^{~k}
    FOR(i, j)
    {
        Mij[i][j] = ricci0.LL[i][j] +
                    vars.K / (3. * chi_regularised) *
                        (vars.A[i][j] + 2. / 3. * vars.K * vars.h[i][j]);
        FOR(k, l)
        {
            Mij[i][j] -=
                vars.A[i][k] * vars.A[j][l] * h_UU[k][l] / chi_regularised;
        }
    }
    data_t M =
        vars.chi *
        compute_trace(Mij, h_UU); // trace of M_{ij} (which is the Einstein Ham)

    Tensor<2, data_t> Mij_TF = Mij;
    make_trace_free(Mij_TF, vars.h, h_UU);
    Tensor<2, data_t> Mij_TF_UU = raise_all(Mij_TF, h_UU); // raise all indexs
    // note that they have been raised with the conformal metric

    data_t divshift = compute_trace(d1.shift);
    data_t dlapse_dot_dchi = compute_dot_product(d1.lapse, d1.chi, h_UU);

    Tensor<2, data_t> covdtilde2lapse;
    Tensor<2, data_t> covd2lapse;
    FOR(k, l)
    {
        covdtilde2lapse[k][l] = d2.lapse[k][l];
        FOR(m) covdtilde2lapse[k][l] -= chris.ULL[m][k][l] * d1.lapse[m];
        covd2lapse[k][l] =
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

    data_t lapse_regularised = simd_max(1e-4, vars.lapse);
    Tensor<2, data_t> Fij_times_lapse;
    FOR(i, j)
    {
        Fij_times_lapse[i][j] =
            -advec.A[i][j] -
            2. / 3. * vars.A[i][j] * (vars.lapse * vars.K - divshift) +
            covd2lapse[i][j];
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
    Tensor<2, data_t> covdtilde_A[CH_SPACEDIM];
    Tensor<2, data_t> covd_Aphys_times_chi[CH_SPACEDIM];

    FOR(i, j, k)
    {
        covdtilde_A[i][j][k] = d1.A[j][k][i];
        FOR(l)
        {
            covdtilde_A[i][j][k] += -chris.ULL[l][i][j] * vars.A[l][k] -
                                    chris.ULL[l][i][k] * vars.A[l][j];
        }
        covd_Aphys_times_chi[i][j][k] =
            covdtilde_A[i][j][k] +
            (vars.A[i][k] * d1.chi[j] + vars.A[i][j] * d1.chi[k]) /
                (2. * chi_regularised);
        FOR(l, m)
        {
            covd_Aphys_times_chi[i][j][k] -=
                h_UU[l][m] * d1.chi[m] / (2. * chi_regularised) *
                (vars.h[i][j] * vars.A[k][l] + vars.h[i][k] * vars.A[j][l]);
        }
    }

    Tensor<1, data_t>
        Ni; // N_i = D^jK_{ij} - D_iK (which is the Einstein Momentum)
    FOR(i) Ni[i] = -(GR_SPACEDIM - 1.) * d1.K[i] / (double)GR_SPACEDIM;
    FOR(i, j, k)
    {
        Ni[i] += h_UU[j][k] * (covdtilde_A[k][j][i] -
                               GR_SPACEDIM * vars.A[i][j] * d1.chi[k] /
                                   (2. * chi_regularised));
    }

    data_t RGB_times_lapse = -4. / 3. * M * F_times_lapse;
    FOR(i, j)
    {
        RGB_times_lapse +=
            8 * Mij_TF_UU[i][j] * Fij_times_lapse[i][j] * vars.chi +
            16. / 3. * vars.chi * vars.lapse * h_UU[i][j] * d1.K[i] *
                (Ni[j] + 1. / 3. * d1.K[j]);
        FOR(k, l, m, n)
        RGB_times_lapse -=
            8. * vars.chi * vars.lapse * h_UU[i][l] * h_UU[j][m] * h_UU[k][n] *
            covd_Aphys_times_chi[l][m][n] *
            (covd_Aphys_times_chi[i][j][k] - covd_Aphys_times_chi[k][i][j]);
    }
    rhs.Kphi += -4. * dfGB * RGB_times_lapse / 4.;
}

// Adds in the RHS for the matter vars
template <class coupling_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t>
void EsGB<coupling_t>::compute_lhs(const int N, data_t *LHS,
                                   const vars_t<data_t> &vars,
                                   const vars_t<Tensor<1, data_t>> &d1,
                                   const diff2_vars_t<Tensor<2, data_t>> &d2,
                                   const vars_t<data_t> &advec) const
{
    data_t LHS_mat[N][N];

    // first get the non potential part of the rhs
    // this may seem a bit long winded, but it makes the function
    // work for more multiple fields

    // set the potential values
    data_t dfGB = 0.0;
    data_t df2GB = 0.0;
    my_coupling.compute_coupling(dfGB, df2GB, vars);

    using namespace TensorAlgebra;

    const auto h_UU = compute_inverse_sym(vars.h);
    const auto chris = compute_christoffel(d1.h, h_UU);

    // Compute useful quantities for the Gauss-Bonnet sector

    const auto ricci0 =
        CCZ4Geometry::compute_ricci_Z(vars, d1, d2, h_UU, chris, {0., 0., 0.});
    // for some reason the function compute_ricci does not give the correct
    // answer

    const data_t chi_regularised = simd_max(1e-6, vars.chi);

    Tensor<2, data_t> Mij; // M_{ij} = R_{ij} + KK_{ij} - K_{ik}K_j^{~k}
    FOR(i, j)
    {
        Mij[i][j] = ricci0.LL[i][j] +
                    vars.K / (3. * chi_regularised) *
                        (vars.A[i][j] + 2. / 3. * vars.K * vars.h[i][j]);
        FOR2(k, l)
        {
            Mij[i][j] -=
                vars.A[i][k] * vars.A[j][l] * h_UU[k][l] / chi_regularised;
        }
    }
    data_t M =
        vars.chi *
        compute_trace(Mij, h_UU); // trace of M_{ij} (which is the Einstein Ham)

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

    // decomposition of C_{ab} = \nabla_a\nabla_bf(\phi)
    Tensor<2, data_t> Cij; // C_{ij} = \gamma^a_{~i}\gamma^b_{~j}C_{ab}
    FOR(i, j)
    {
        Cij[i][j] = 4. * dfGB * (covd2phi[i][j] -
                            vars.Kphi / chi_regularised *
                                (vars.A[i][j] +
                                 vars.h[i][j] * vars.K / (double)GR_SPACEDIM)) +
                    4. * df2GB * d1.phi[i] * d1.phi[j];
    }
    data_t C = vars.chi * compute_trace(Cij, h_UU); // trace of Cij

    Tensor<1, data_t> Ci; // C_i = -\gamma^a_{~i}n^bC_{ab}
    FOR(i)
    {
        Ci[i] = 4. * df2GB * vars.Kphi * d1.phi[i] +
                4. * dfGB * (d1.Kphi[i] - vars.K * d1.phi[i] / (double)GR_SPACEDIM);
        FOR(j, k) { Ci[i] += - 4. * dfGB * h_UU[j][k] * d1.phi[k] * vars.A[i][j]; }
    }

    Tensor<2, data_t> Mij_TF = Mij;
    make_trace_free(Mij_TF, vars.h, h_UU);
    Tensor<2, data_t> Mij_TF_UU = raise_all(Mij_TF, h_UU); // raise all indexs
    // note that they have been raised with the conformal metric

    Tensor<2, data_t> Cij_TF = Cij;
    make_trace_free(Cij_TF, vars.h, h_UU);
    Tensor<2, data_t> Cij_TF_UU = raise_all(Cij_TF, h_UU); // raise all indexs
    // note that they have been raised with the conformal metric

    int n1 = 0;
    FOR(a1, b1)
    {
        if (a1 > b1) // lower diagonal
            continue;

        int n2 = 0;
        FOR(a2, b2)
        {
            if (a2 > b2) // lower diagonal
                continue;
            // minus sign because the 2nd time derivative are moved from the RHS
            // to the LHS
            // Lapack LHS argument has row-by-row layout
            // n1 is the 'line' (multiplying variables 'x' in A.x = b), 'n2' is
            // the column (which equation)
            LHS_mat[n1][n2] =
                -2.0 * vars.chi *
                (1. / 3. * vars.h[a2][b2] * Cij_TF_UU[a1][b1] +
                 4. * dfGB * 4. * dfGB * Mij_TF[a2][b2] * Mij_TF_UU[a1][b1] * vars.chi);
            if (a1 == a2)
            {
                if (b1 == b2)
                    LHS_mat[n1][n2] += 1. - C / 3.;
                FOR(k)
                LHS_mat[n1][n2] += vars.chi * h_UU[b1][k] * Cij_TF[b2][k];
            }
            if (b1 == b2)
                FOR(k)
                LHS_mat[n1][n2] += vars.chi * h_UU[a1][k] * Cij_TF[a2][k];
            if (a1 != b1)
            {
                LHS_mat[n1][n2] +=
                    -2.0 * vars.chi *
                    (1. / 3. * vars.h[a2][b2] * Cij_TF_UU[b1][a1] +
                     4. * dfGB * 4. * dfGB * Mij_TF[a2][b2] * Mij_TF_UU[b1][a1] *
                         vars.chi);
                if (a1 == b2)
                {
                    if (a2 == b1)
                        LHS_mat[n1][n2] += 1. - C / 3.;
                    FOR(k)
                    LHS_mat[n1][n2] += vars.chi * h_UU[b1][k] * Cij_TF[a2][k];
                }
                if (a2 == b1)
                    FOR(k)
                    LHS_mat[n1][n2] += vars.chi * h_UU[a1][k] * Cij_TF[b2][k];
            }

            ++n2;
        }

        ++n1;
    }

    n1 = 0;
    FOR(a, b)
    {
        if (a > b) // lower diagonal
            continue;

        LHS_mat[N - 2][n1] =
            vars.chi / 3. * (-Cij_TF[a][b] + 4. * dfGB * 4. * dfGB * M * Mij_TF[a][b]);

        LHS_mat[n1][N - 2] =
            vars.chi / 2. *
            (Cij_TF_UU[a][b] - 4. * dfGB * 4. * dfGB * M * Mij_TF_UU[a][b]);
        if (a != b)
            LHS_mat[n1][N - 2] +=
                vars.chi / 2. *
                (Cij_TF_UU[b][a] - 4. * dfGB * 4. * dfGB * M * Mij_TF_UU[b][a]);
        ++n1;
    }

    LHS_mat[N - 2][N - 2] = 1. + 1. / 3. * (-C + dfGB * 4. * dfGB * M * M);

    n1 = 0;
    FOR(a, b)
    {
        if (a > b) // lower diagonal
            continue;

        LHS_mat[N - 1][n1] = 0.;
        LHS_mat[n1][N - 1] = 2. * dfGB * Mij_TF_UU[a][b] * vars.chi;
        if (a != b)
            LHS_mat[n1][N - 1] += 2. * dfGB * Mij_TF_UU[b][a] * vars.chi;
        ++n1;
    }
    LHS_mat[N - 1][N - 2] = 0.;
    LHS_mat[N - 2][N - 1] = -1. / 3. * dfGB * M;
    LHS_mat[N - 1][N - 1] = 1.;

    for (int n1 = 0; n1 < N; ++n1)
    {
        for (int n2 = 0; n2 < N; ++n2)
        {
            LHS[n1 * N + n2] = LHS_mat[n1][n2];
        }
    }
}

// Function to add in EM Tensor matter terms to CCZ4 rhs
template <class coupling_t>
template <class data_t, template <typename> class vars_t,
          template <typename> class diff2_vars_t,
          template <typename> class rhs_vars_t>
void EsGB<coupling_t>::solve_lhs(rhs_vars_t<data_t> &rhs,
                                 const vars_t<data_t> &vars,
                                 const vars_t<Tensor<1, data_t>> &d1,
                                 const diff2_vars_t<Tensor<2, data_t>> &d2,
                                 const vars_t<data_t> &advec) const
{
    const int N = GR_SPACEDIM * (GR_SPACEDIM + 1) / 2 + 2;
    data_t LHS[N][N];

    compute_lhs(N, (&LHS[0][0]), vars, d1, d2, advec);

    data_t RHS[N];

    int n1 = 0;
    FOR(a1, b1)
    {
        if (a1 > b1) // lower diagonal
            continue;

        RHS[n1] = rhs.A[a1][b1];
        if (a1 != b1)
            RHS[n1] = rhs.A[b1][a1];

        ++n1;
    }
    RHS[N - 2] = rhs.K;
    RHS[N - 1] = rhs.Kphi;

    EsGB_lapack::solve_system_lapack(N, (&LHS[0][0]), RHS);

    n1 = 0;
    FOR(a1, b1)
    {
        if (a1 > b1) // lower diagonal
            continue;

        rhs.A[a1][b1] = RHS[n1];
        if (a1 != b1)
            rhs.A[b1][a1] = RHS[n1];

        ++n1;
    }
    rhs.K = RHS[N - 2];
    rhs.Kphi = RHS[N - 1];
}

#endif /* ESGB_IMPL_HPP_ */
