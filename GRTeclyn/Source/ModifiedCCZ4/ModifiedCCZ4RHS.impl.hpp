/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(MODIFIEDCCZ4RHS_HPP_)
#error "This file should only be included through ModifiedCCZ4RHS.hpp"
#endif

#ifndef MODIFIEDCCZ4RHS_IMPL_HPP_
#define MODIFIEDCCZ4RHS_IMPL_HPP_
#include "DimensionDefinitions.hpp"

template <class theory_t, class deriv_t>
ModifiedCCZ4RHS<theory_t, deriv_t>::ModifiedCCZ4RHS(amrex::Real a_dx)
    : CCZ4RHS<deriv_t>(a_dx, 0.0 /*No cosmological constant*/)
{
    GRParmParse mod_gauge_pp("mod_gauge");
    mod_gauge_pp.get("mod_b", m_mod_b);
}

template <class theory_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
ModifiedCCZ4RHS<theory_t, deriv_t>::add_b_rhs(
    const int ix, const int iy, const int iz,
    const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<const amrex::Real> &state) const
{
    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);
    const typename theory_t::Vars vars(state_cell_data);
    
    // Construct derivatives
    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    auto d1_h = this->m_deriv.d1_sym_tensor(ix, iy, iz, state, c_h11);
    const auto chris = CCZ4Geometry::compute_christoffel(d1_h, h_UU);
    
    auto d1_A = this->m_deriv.d1_sym_tensor(ix, iy, iz, state, c_A11);
    auto d1_K     = this->m_deriv.d1_scalar(ix, iy, iz, state, c_K);
    auto d1_chi   = this->m_deriv.d1_scalar(ix, iy, iz, state, c_chi);
    auto d1_lapse = this->m_deriv.d1_scalar(ix, iy, iz, state, c_lapse);
    auto d1_Theta = this->m_deriv.d1_scalar(ix, iy, iz, state, c_Theta);
    auto d1_Gamma = this->m_deriv.d1_vector(ix, iy, iz, state, c_Gamma1);
    auto d2_chi   = this->m_deriv.d2_scalar(ix, iy, iz, state, c_chi);
    auto d2_h = this->m_deriv.d2_sym_tensor(ix, iy, iz, state, c_h11);

    // Compute ricci
    Tensor::Rank1 zero_Z{};
    auto ricci = CCZ4Geometry::compute_ricci_Z(
        vars, d1_chi, d1_Gamma, d1_h, d2_h, d2_chi, h_UU, chris, zero_Z);

    // Compute Z4
    const amrex::Real non_covariant_z4 = 1.0 - this->m_params.covariant_z4_coeff;
    const amrex::Real kappa1_times_lapse =
        this->m_params.covariant_z4_coeff * this->m_params.kappa1 +
        non_covariant_z4 * this->m_params.kappa1 * vars.lapse();

    Tensor::Rank1 Z_over_chi;

    // Select CCZ4 without introducing a branch into the GPU kernel.
    const amrex::Real ccz4_coeff = 1.0 - this->m_params.bssn_coeff;

    FOR (i)
    {
        Z_over_chi(i) =
            ccz4_coeff * 0.5 * (vars.Gamma(i) - chris.contracted(i));
    }

    // This is A_ij A^ij
    amrex::Real Aij_squared = CCZ4Geometry::compute_Aij_squared(vars, h_UU);

    // Compute Hamiltonian constraint
    amrex::Real Ham = ricci.scalar +
                  (GR_SPACEDIM - 1.0) * vars.K() * vars.K() / GR_SPACEDIM -
                  Aij_squared;
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
    // Compute momentum constraint
    Tensor::Rank1 Mom{};
    FOR (i)
    {
        Mom(i) = -(GR_SPACEDIM - 1.0) * d1_K(i) / GR_SPACEDIM;
    }
    FOR (i, j, k)
    {
        Mom(i) += h_UU(j, k) *
                  (covd_A(k, j, i) - GR_SPACEDIM * vars.A(i, j) *
                       d1_chi(k) / (2.0 * vars.chi()));
    }

    // Update evolution equations (pending to include BSSN option as well)
    amrex::Real factor_mod_b = m_mod_b / (1.0 + m_mod_b);
    //amrex::Real factor_mod_a = m_mod_a / (1. + m_mod_a);
    rhs_cell_data[c_K] += GR_SPACEDIM * factor_mod_b * 
	    (-0.5 / (GR_SPACEDIM - 1.) * vars.lapse() * Ham + 
	     kappa1_times_lapse * vars.Theta() *
                 (1.0 + 0.5 * this->m_params.kappa2));

    rhs_cell_data[c_Theta] += 0.5 * factor_mod_b * (-vars.lapse() * Ham +
         vars.Theta() * kappa1_times_lapse * 
	      ((GR_SPACEDIM - 3.0) / (2.0 + m_mod_b) +
	      (GR_SPACEDIM + 1.0) + this->m_params.kappa2 * (GR_SPACEDIM - 1.)));

    Tensor::Rank2 A_UU = CCZ4Geometry::compute_A_UU(vars, h_UU);
    FOR (i)
    {
	amrex::Real mod_gauge_term_Gamma = 2.0 * factor_mod_b * Z_over_chi(i) * 
		(vars.lapse() * vars.K() / GR_SPACEDIM  + 
		kappa1_times_lapse);
	FOR (j)
	{
	    mod_gauge_term_Gamma += 
                -factor_mod_b * 2.0 * h_UU(i, j) * vars.lapse() * 
		        (d1_Theta(j) + Mom(j));
	    FOR (k)
	    {
	       mod_gauge_term_Gamma += factor_mod_b * 2.0 * vars.lapse() * 
		       A_UU(i, j) * vars.h(j, k) * Z_over_chi(k);
	    }        
	}
        rhs_cell_data[c_Gamma1 + i] += mod_gauge_term_Gamma;
    }
}

// Function to add in EM Tensor matter terms to CCZ4 RHS
template <class theory_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
ModifiedCCZ4RHS<theory_t, deriv_t>::add_emtensor_rhs(
    const int ix, const int iy, const int iz,
    const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<const amrex::Real> &state) const
{
    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);

    const typename theory_t::Vars vars(state_cell_data);

    const auto h_UU = CCZ4Geometry::compute_inverse_metric(vars);

    const auto source = m_theory.compute_einstein_sources(ix, iy, iz, state,
                                                          this->m_deriv, h_UU);

    // Select the matter source terms without branching in the GPU kernel.
    const amrex::Real ccz4_coeff = 1.0 - this->m_params.bssn_coeff;

    const amrex::Real ccz4_K_matter_rhs =
        0.5 * vars.lapse() * (source.trS - 3.0 * source.rho / (1.0 + m_mod_b));
    const amrex::Real bssn_K_matter_rhs =
        0.5 * vars.lapse() * (source.trS + source.rho / (1.0 + m_mod_b));
    rhs_cell_data[c_K] += ccz4_coeff * ccz4_K_matter_rhs +
                          this->m_params.bssn_coeff * bssn_K_matter_rhs;

    const amrex::Real ccz4_Theta_matter_rhs = -vars.lapse() * source.rho / (1.0 + m_mod_b);
    rhs_cell_data[c_Theta] =
        ccz4_coeff * (rhs_cell_data[c_Theta] + ccz4_Theta_matter_rhs);

    // Update RHS for other variables

    FOR2_SYM(i, j)
    {

        rhs_cell_data[sym_var_idx(c_A11, i, j)] -=
            vars.chi() * vars.lapse() * source.S_TF(i, j);
    }

    FOR (i)
    {
        amrex::Real matter_term_Gamma = 0.0;
        FOR (j)
        {
            matter_term_Gamma -= 2.0 * vars.lapse() * h_UU(i, j) * source.j(j) / (1.0 + m_mod_b);
        }
        rhs_cell_data[c_Gamma1 + i] += matter_term_Gamma;
    }
}

template <class theory_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
ModifiedCCZ4RHS<theory_t, deriv_t>::add_theory_rhs(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<const amrex::Real> &state) const
{
    m_theory.add_theory_rhs(ix, iy, iz, rhs_state, state, this->m_deriv);
}

template <class theory_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
ModifiedCCZ4RHS<theory_t, deriv_t>::solve_lhs(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<const amrex::Real> &state) const
{
    m_theory.solve_lhs(ix, iy, iz, rhs_state, state, this->m_deriv);
}

template <class theory_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
ModifiedCCZ4RHS<theory_t, deriv_t>::apply_dissipation(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<const amrex::Real> &state) const
{
    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    this->m_deriv.add_dissipation(ix, iy, iz, rhs_cell_data, state,
                                  this->m_sigma, NUM_VARS);
}

template <class theory_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Rank2
ModifiedCCZ4RHS<theory_t, deriv_t>::get_full_kappa_Sij_TF(
    int ix, int iy, int iz, 
    const amrex::Array4<const amrex::Real> &state) const
{
    amrex::Real rhs_data[NUM_VARS]{};
    const amrex::Dim3 cell_begin{ix, iy, iz};
    const amrex::Dim3 cell_end{ix + 1, iy + 1, iz + 1};
    const amrex::Array4<amrex::Real> rhs_state(
        rhs_data, cell_begin, cell_end, NUM_VARS);

    this->compute_A_ij_and_Theta_and_Gamma(ix, iy, iz, rhs_state, state);
    add_b_rhs(ix, iy, iz, rhs_state, state);

    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    Tensor::Rank2 A_rhs_GR{};
    FOR2_SYM(i, j)
    {
        A_rhs_GR(i, j) = rhs_cell_data[sym_var_idx(c_A11, i, j)];
    }

    add_emtensor_rhs(ix, iy, iz, rhs_state, state);
    add_theory_rhs(ix, iy, iz, rhs_state, state);
    solve_lhs(ix, iy, iz, rhs_state, state);

    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);
    const typename theory_t::Vars vars(state_cell_data);

    Tensor::Rank2 out{};
    FOR2_SYM(i, j)
    {
        const amrex::Real component =
            (A_rhs_GR(i, j) - rhs_cell_data[sym_var_idx(c_A11, i, j)]) /
            vars.chi();
        out(i, j) = component;
        out(j, i) = component;
    }

    return out;
}

#endif /* MODIFIEDCCZ4RHS_IMPL_HPP_ */
