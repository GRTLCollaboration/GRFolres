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
}

template <class theory_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
ModifiedCCZ4RHS<theory_t, deriv_t>::add_a_and_b_rhs(
    const int ix, const int iy, const int iz,
    const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<const amrex::Real> &state) const
{
    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    const typename theory_t::Vars vars(state_cell_data);
    
    // Construct derivatives
    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    auto d1_h = this->m_deriv.d1_sym_tensor(ix, iy, iz, state, c_h11);
    const auto chris = CCZ4Geometry::compute_christoffel(d1_h, h_UU);
    
    auto d1_A = this->m_deriv.d1_sym_tensor(ix, iy, iz, state, c_A11);
    auto d1_chi   = this->m_deriv.d1_scalar(ix, iy, iz, state, c_chi);
    auto d1_Theta = this->m_deriv.d1_scalar(ix, iy, iz, state, c_Theta);
    auto d1_Gamma = this->m_deriv.d1_vector(ix, iy, iz, state, c_Gamma1);
    auto d2_chi   = this->m_deriv.d2_scalar(ix, iy, iz, state, c_chi);
    auto d2_h = this->m_deriv.d2_sym_tensor(ix, iy, iz, state, c_h11);

    auto ricci = CCZ4Geometry::compute_ricci(vars, d1_chi, d1_Gamma, d1_h,
                                                 d2_chi, d2_h, h_UU, chris);

    amrex::Real factor_b = m_params.modb / (1. + m_params.modb);

    // This is A_ij A^ij
    amrex::Real Aij_squared = CCZ4Geometry::compute_Aij_squared(vars, h_UU);

    amrex::Real Ham = ricci.scalar +
                  (GR_SPACEDIM - 1.) * vars.K() * vars.K() / GR_SPACEDIM -
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
    Tensor::Rank1 Mom{};
    FOR (i)
    {
        Mom(i) = -(GR_SPACEDIM - 1.) * d1_K(i) / GR_SPACEDIM;
    }

    rhs_cell_data[c_K] += -0.5 * GR_SPACEDIM / (GR_SPACEDIM - 1.) * factor_b * vars.lapse() * Ham;
    rhs_cell_data[c_Theta] += -0.5 * factor_b * vars.lapse() * Ham;

    FOR2_SYM(i, j) rhs_cell_data[c_Gamma1 + i] += -factor_b * 2. * h_UU(i, j) * vars.lapse() * (d1_Theta(j) + Mom(j)) ;

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
        0.5 * vars.lapse() * (source.trS - 3.0 * source.rho);
    const amrex::Real bssn_K_matter_rhs =
        0.5 * vars.lapse() * (source.trS + source.rho);
    rhs_cell_data[c_K] += ccz4_coeff * ccz4_K_matter_rhs +
                          this->m_params.bssn_coeff * bssn_K_matter_rhs;

    const amrex::Real ccz4_Theta_matter_rhs = -vars.lapse() * source.rho;
    rhs_cell_data[c_Theta] =
        ccz4_coeff * (rhs_cell_data[c_Theta] + ccz4_Theta_matter_rhs);

    // Update RHS for other variables
    Tensor::Rank2 S_TF = source.S;

    CCZ4Geometry::make_trace_free(S_TF, vars, h_UU);

    FOR2_SYM(i, j)
    {

        rhs_cell_data[sym_var_idx(c_A11, i, j)] -=
            vars.chi() * vars.lapse() * S_TF(i, j);
    }

    FOR (i)
    {
        amrex::Real matter_term_Gamma = 0.0;
        FOR (j)
        {
            matter_term_Gamma -= 2.0 * vars.lapse() * h_UU(i, j) * source.j(j);
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
ModifiedCCZ4RHS<theory_t, deriv_t>::apply_dissipation(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &rhs_state,
    const amrex::Array4<const amrex::Real> &state) const
{
    const amrex::CellData<amrex::Real> &rhs_cell_data =
        rhs_state.cellData(ix, iy, iz);
    this->m_deriv.add_dissipation(ix, iy, iz, rhs_cell_data, state,
                                  this->m_sigma, NUM_VARS);
}

#endif /* MODIFIEDCCZ4RHS_IMPL_HPP_ */
