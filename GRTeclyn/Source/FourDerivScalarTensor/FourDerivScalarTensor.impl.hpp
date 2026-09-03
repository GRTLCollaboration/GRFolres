/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(FOURDERIVSCALARTENSOR_HPP_)
#error "This file should only be included through FourDerivScalarTensor.hpp"
#endif

#ifndef FOURDERIVSCALARTENSOR_IMPL_HPP_
#define FOURDERIVSCALARTENSOR_IMPL_HPP_

template <class coupling_and_potential_t, class deriv_t>
AMREX_GPU_DEVICE ScalarVectorTensor FourDerivScalarTensor<coupling_and_potential_t, deriv_t>::compute_M_Ni_and_Mij(
    const int ix, const int iy, const int iz,
    const amrex::Array4<const amrex::Real> &state,
    const deriv_t &a_deriv) const
{
    ScalarVectorTensor out;

    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const typename theory_t::Vars vars(state_cell_data);

    // Construct derivatives
    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    auto d1_h = a_deriv.d1_sym_tensor(ix, iy, iz, state, c_h11);
    const auto chris = CCZ4Geometry::compute_christoffel(d1_h, h_UU);

    auto d1_A = a_deriv.d1_sym_tensor(ix, iy, iz, state, c_A11);
    auto d1_K     = a_deriv.d1_scalar(ix, iy, iz, state, c_K);
    auto d1_chi   = a_deriv.d1_scalar(ix, iy, iz, state, c_chi);
    auto d1_Gamma = a_deriv.d1_vector(ix, iy, iz, state, c_Gamma1);
    auto d2_chi   = a_deriv.d2_scalar(ix, iy, iz, state, c_chi);
    auto d2_h = a_deriv.d2_sym_tensor(ix, iy, iz, state, c_h11);

    // Compute Ricci
    auto ricci = CCZ4Geometry::compute_ricci(vars, d1_chi, d1_Gamma, d1_h,
                                                 d2_chi, d2_h, h_UU, chris);

    // M_{ij} = R_{ij} + KK_{ij} - K_{ik}K_j^{~k}
    FOR(i, j)
    {
	out.tensor(i, j) = ricci.LL(i, j) + vars.K() / (3.0 * chi) *
	                 (vars.A(i, j) + 2. / 3. * vars.K() * vars.h(i, j));
	FOR(k, l)
	{
            out.tensor(i, j) += -vars.A(i, k) * vars.A(j, l) * h_UU(k, l) / vars.chi();
	}
    }
    // M = \gamma^{ij}M_{ij} (vacuum GR Hamiltonian constraint)
    out.scalar = vars.chi() * TensorAlgebra::compute_trace(out.tensor, h_UU);

    // Covariant derivative of \bar A_ij
    Tensor::Rank3 covd_A{};
    FOR (i, j, k)
    {
        covd_A(i, j, k) = d1_A(j, k, i);
        FOR (l)
        {
            covd_A(i, j, k) += -chris.ULL(l, i, j) * vars.A(l, k) -
                               chris.ULL(l, i, k) * vars.A(l, j);
        }
    }
    // N_i = D^jK_{ij} - D_iK (vacuum GR momentum constraint)
    FOR (i)
    {
        out.vector(i) = -(GR_SPACEDIM - 1.) * d1_K(i) / GR_SPACEDIM;
	FOR (j, k)
	{
	    out.vector(i) += h_UU(j, k) * (covdA(i, j, k) - 
			    0.5 * GR_SPACEDIM * vars.A(i, j) * d1.chi(k) / vars.chi();
	}
    }
    return out;
}

template <class coupling_and_potential_t, class deriv_t>
AMREX_GPU_DEVICE ScalarVectorTensor FourDerivScalarTensor<coupling_and_potential_t, deriv_t>::compute_Omega_munu(
    const int ix, const int iy, const int iz,
    const amrex::Array4<const amrex::Real> &state,
    const deriv_t &a_deriv) const
{
    ScalarVectorTensor out;

    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const typename theory_t::Vars vars(state_cell_data);

    // Construct derivatives
    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    auto d1_h = a_deriv.d1_sym_tensor(ix, iy, iz, state, c_h11);
    const auto chris = CCZ4Geometry::compute_christoffel(d1_h, h_UU);

    auto d1_chi   = a_deriv.d1_scalar(ix, iy, iz, state, c_chi);
    auto d1_phi   = a_deriv.d1_scalar(ix, iy, iz, state, c_phi);
    auto d1_Pi    = a_deriv.d1_scalar(ix, iy, iz, state, c_Pi);
    auto d2_phi   = a_deriv.d2_scalar(ix, iy, iz, state, c_phi);

    // set the coupling and potential values
    amrex::Real dfdphi   = 0.0;
    amrex::Real d2fphi2  = 0.0;
    amrex::Real V_of_phi = 0.0;
    amrex::Real dVdphi   = 0.0;
    amrex::Real g2       = 0.0;
    amrex::Real dg2dphi  = 0.0;
    // compute coupling and potential and add constributions to EM Tensor
    m_coupling_and_potential.compute_coupling_and_potential(dfdphi, d2fdphi2, V_of_phi, dVdphi, g2, dg2dphi, vars);

    // relevant quantities
    amrex::Real dphi_dot_chi = 
	    TensorAlgebra::compute_dot_product(d1_phi, d1_chi, h_UU);
    Tensor::Rank2 covd2phi_phys_times_chi{};
    FOR(i, j)
    {
        covd2phi_phys_times_chi(i, j) = d2_phi(i, j)
	FOR (k) covd2phi_phys_times_chi(i, j) += -chris.ULL(k, i, j) * d1_phi(k);
	FOR(l, m)
	{
            covd2phi_phys_times_chi(i, j) += 0.5 * (d1_phi(i) * d1_chi(j) +
	        d1_phi(j) * d1_chi(i) - vars.h(i, j) * dphi_dot_dchi / vars.chi();
	}
    }

    // Omega_{ij}=\gamma^{\mu}_{~i}\gamma^{\nu}_{~j}\Omega_{\mu\nu}
    FOR(i, j)
    {
        out.tensor(i, j) = 4.0 * dfdphi * 
	                   (covd2phi_phys_times_chi(i, j) + vars.Pi() / vars.chi() *
			       (vars.A(i, j) + vars.h(i, j) * vars.K() / GR_SPACEDIM) +
                           4.0 * d2fdfphi2 * d1_phi(i) * d1_phi(j);
    }

    // trace of Omega_ij
    out.scalar = vars.chi() * TensorAlgebra::compute_trace(out.tensor, h_UU);

    // Omega_i = -\gamma^{\mu}_{~i}n^{\nu}\Omega_{\mu\nu}
    FOR (i)
    {
        out.vector(i) = -4.0 * d2fdphi2 * vars.Pi() * d1_phi(i) +
	                4.0 * dfdphi * (-d1_Pi(i) - vars.K() * d1_phi(i) / GR_SPACEDIM);
	FOR(j, k)
	{
	    out.vector(i) += -4.0 * dfdphi * h_UU(j, k) * d1_phi(k) * vars.A(i, j);
	}
    }
    return out;
}

// Calculate the rho and j components of the effective stress energy tensor
// (ONLY MINIMALLY COUPLED SCALAR FIELD TERMS FOR NOW)
template <class coupling_and_potential_t, class deriv_t>
AMREX_GPU_DEVICE emtensor_t FourDerivScalarTensor<coupling_and_potential_t, deriv_t>::compute_rho_and_j(
    const int ix, const int iy, const int iz,
    const amrex::Array4<const amrex::Real> &state, 
    const deriv_t &a_deriv,
    const Tensor::Rank2 &h_UU) const
{
    RhoAndJ out;

    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    const Vars vars(state_cell_data);

    auto d1_phi = a_deriv.d1_scalar(ix, iy, iz, state, c_phi);

    //    Useful quantity Vt
    amrex::Real Vt = -vars.Pi() * vars.Pi();
    FOR (i, j)
    {
        Vt += vars.chi() * h_UU(i, j) * d1_phi(i) * d1_phi(j);
    }

    // set the coupling and potential values
    amrex::Real dfdphi   = 0.0;
    amrex::Real d2fphi2  = 0.0;
    amrex::Real V_of_phi = 0.0;
    amrex::Real dVdphi   = 0.0;
    amrex::Real g2       = 0.0;
    amrex::Real dg2dphi  = 0.0;
    // compute coupling and potential and add constributions to EM Tensor
    m_coupling_and_potential.compute_coupling_and_potential(dfdphi, d2fdphi2, V_of_phi, dVdphi, g2, dg2dphi, vars);

    // Calculate components of EM Tensor
    // rho = n^a n^b T_ab
    out.rho = vars.Pi() * vars.Pi() + 0.5 * Vt + V_of_phi;

    //    j_i (note lower index) = - n^a T_ai
    FOR (i)
    {
        out.j(i) = -d1_phi(i) * vars.Pi();
    }

    return out;
}

// Calculate the S_TF and trS component of the effective stress energy tensor
// (ONLY MINIMALLY COUPLED SCALAR FIELD FOR NOW)
template <class coupling_and_potential_t, class deriv_t>
AMREX_GPU_DEVICE S_TFAndTrS FourDerivScalarTensor<coupling_and_potential_t, deriv_t>::compute_S_TF_and_trS(
    const int ix, const int iy, const int iz,
    const amrex::Array4<const amrex::Real> &state,
    const deriv_t &a_deriv,
    const Tensor::Rank2 &h_UU) const
{
    S_TFAndTrS out;

    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    const Vars vars(state_cell_data);

    auto d1_phi = a_deriv.d1_scalar(ix, iy, iz, state, c_phi);

    // Useful quantity Vt
    amrex::Real Vt = -vars.Pi() * vars.Pi();
    FOR (i, j)
    {
        Vt += vars.chi() * h_UU(i, j) * d1_phi(i) * d1_phi(j);
    }

    // set the coupling and potential values
    amrex::Real dfdphi   = 0.0;
    amrex::Real d2fphi2  = 0.0;
    amrex::Real V_of_phi = 0.0;
    amrex::Real dVdphi   = 0.0;
    amrex::Real g2       = 0.0;
    amrex::Real dg2dphi  = 0.0;
    // compute coupling and potential and add constributions to EM Tensor
    m_coupling_and_potential.compute_coupling_and_potential(dfdphi, d2fdphi2, V_of_phi, dVdphi, g2, dg2dphi, vars);

    // Calculate components of EM Tensor
    // S_TF = T_ij^{TF}
    FOR (i, j)
    {
        out.S_TF(i, j) = -0.5 * vars.h(i, j) * Vt / vars.chi() +
                      d1_phi(i) * d1_phi(j) -
                      vars.h(i, j) * V_of_phi / vars.chi();
    }

    // trS = Tr_S_ij
    out.trS = vars.chi() * TensorAlgebra::compute_trace(out.S, h_UU);
    
    TensorAlgebra::make_trace_free(out.S, vars.h(), h_UU);

    return out;
}


template <class coupling_and_potential_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE einstein_sources_t
FourDerivScalarTensor<coupling_and_potential_t, deriv_t>::compute_einstein_sources(
    int ix, int iy, int iz, const amrex::Array4<const amrex::Real> &state,
    const deriv_t &a_deriv, const Tensor::Rank2 &h_UU) const
{
    const RhoAndJ rho_and_j =
        compute_rho_and_j(ix, iy, iz, state, a_deriv, h_UU);
    const S_TFAndTrS S_TF_and_trS;
    const amrex::Real coupling = 8.0 * M_PI * m_G_Newton;

    einstein_sources_TF out;
    out.rho = coupling * rho_and_j.rho;
    out.trS = coupling * S_TF_and_trS.trS;
    FOR (i)
    {
        out.j(i) = coupling * rho_and_j.j(i);
    }
    FOR (i, j)
    {
        out.S_TF(i, j) = coupling * S_TF_and_trS.S_TF(i, j);
    }
    return out;
}

// Adds in the RHS for the theory vars (ONLY MINIMALLY COUPLED SCALAR FIELD TERMS FOR NOW)
template <class coupling_and_potential_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
FourDerivScalarTensor<coupling_and_potential_t, deriv_t>::add_theory_rhs(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<const amrex::Real> &state, const deriv_t &a_deriv) const
{
    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    const Vars vars(state_cell_data);

    // call the function for the rhs excluding the potential
    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    auto d1_h        = a_deriv.d1_sym_tensor(ix, iy, iz, state, c_h11);
    const auto chris = CCZ4Geometry::compute_christoffel(d1_h, h_UU);

    // calculate the derivatives
    auto d1_chi   = a_deriv.d1_scalar(ix, iy, iz, state, c_chi);
    auto d1_lapse = a_deriv.d1_scalar(ix, iy, iz, state, c_lapse);

    auto d1_phi = a_deriv.d1_scalar(ix, iy, iz, state, c_phi);
    auto d1_Pi  = a_deriv.d1_scalar(ix, iy, iz, state, c_Pi);

    auto d2_phi = a_deriv.d2_scalar(ix, iy, iz, state, c_phi);

    Tensor::Rank1 shift_vector{vars.shift(0), vars.shift(1), vars.shift(2)};

    auto advec_phi =
        a_deriv.advec_scalar(ix, iy, iz, state, shift_vector, c_phi);
    auto advec_Pi = a_deriv.advec_scalar(ix, iy, iz, state, shift_vector, c_Pi);

    // set the coupling and potential values
    amrex::Real dfdphi   = 0.0;
    amrex::Real d2fdphi2 = 0.0;
    amrex::Real V_of_phi = 0.0;
    amrex::Real dVdphi   = 0.0;
    amrex::Real g2       = 0.0;
    amrex::Real dg2dphi  = 0.0;
    m_coupling_and_potential.compute_coupling_and_potential(dfdphi, d2fdphi2, V_of_phi, dVdphi, g2, dg2dphi, vars);

    // evolution equations for scalar field and (minus) its conjugate momentum
    rhs_cell_data[c_phi] = vars.lapse() * vars.Pi() + advec_phi;

    rhs_cell_data[c_Pi] =
        vars.lapse() * (vars.K() * vars.Pi() - dVdphi) + advec_Pi;

    FOR (i, j)
    {
        // includes non conformal parts of chris not included in chris_ULL
        rhs_cell_data[c_Pi] +=
            h_UU(i, j) * (-0.5 * d1_chi(j) * vars.lapse() * d1_phi(i) +
                          vars.chi() * vars.lapse() * d2_phi(i, j) +
                          vars.chi() * d1_lapse(i) * d1_phi(j));
        FOR (k)
        {
            rhs_cell_data[c_Pi] += -vars.chi() * vars.lapse() * h_UU(i, j) *
                                   chris.ULL(k, i, j) * d1_phi(k);
        }
    }
}

// Computes the LHS matrix (IN PROGRESS)
template <class coupling_and_potential_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
FourDerivScalarTensor<coupling_and_potential_t, deriv_t>::compute_lhs(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<const amrex::Real> &state, const deriv_t &a_deriv,
    const int matrix_dim, amrex::Real *LHS) const
{
    amrex::Real LHS_mat[matrix_dim][matrix_dim];

    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    const Vars vars(state_cell_data);

    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    auto d1_h        = a_deriv.d1_sym_tensor(ix, iy, iz, state, c_h11);
    const auto chris = CCZ4Geometry::compute_christoffel(d1_h, h_UU);
    
    // set the coupling and potential values
    amrex::Real dfdphi   = 0.0;
    amrex::Real d2fdphi2 = 0.0;
    amrex::Real V_of_phi = 0.0;
    amrex::Real dVdphi   = 0.0;
    amrex::Real g2       = 0.0;
    amrex::Real dg2dphi  = 0.0;
    m_coupling_and_potential.compute_coupling_and_potential(dfdphi, d2fdphi2, V_of_phi, dVdphi, g2, dg2dphi, vars);

    // Compute useful quantities for the Gauss-Bonnet sector

    ScalarVectorTensor SVT = compute_M_Ni_and_Mij(ix, iy, iz, rhs_state, state, a_deriv);
    amrex::Real M = SVT.scalar;
    Tensor::Rank2 Mij = SVT.tensor;

    for (int row = 0; row < matrix_dim; ++row)
    {
        for (int col = 0; col < matrix_dim; ++col)
        {
            LHS[col * matrix_dim + row] = LHS_mat[row][matrix_dim];
        }
    }
}

template <class coupling_and_potential_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
FourDerivScalarTensor<coupling_and_potential_t, deriv_t>::solve_lhs(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<const amrex::Real> &state, const deriv_t &a_deriv) const
{
    const int matrix_dim = GR_SPACEDIM * (GR_SPACEDIM + 1) / 2 + 2;
    amrex::Real LHS[matrix_dim][matrix_dim];

    compute_lhs(ix, iy, iz, rhs_state, state, a_deriv, matrix_dim, (&LHS[0][0]));
    amrex::Real RHS[matrix_dim];

    int row = 0;
    FOR2_SYM(i, j)
    {
        RHS[row] = rhs_cell_data[sym_var_idx(c_A11, i, j)];
        ++row;
    }

    RHS[matrix_dim - 2] = rhs.K;
    RHS[matrix_dim - 1] = rhs.Pi;

    solve_linear_system(matrix_dim, (&LHS[0][0]), RHS);

    row = 0;
    FOR(i, j)
    {
        rhs_cell_data[sym_var_idx(c_A11, i, j)] = RHS[row];
        ++row;
    }
    rhs_cell_data[c_K] = RHS[N - 2];
    rhs_cell_data[c_Pi] = RHS[N - 1];

}

#endif /* FOURDERIVSCALARTENSOR_IMPL_HPP_ */
