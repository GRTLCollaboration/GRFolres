/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(MODIFIEDGRAVITYWEYL4_HPP_)
#error "This file should only be included through ModifiedGravityWeyl4.hpp"
#endif

#ifndef MODIFIEDGRAVITYWEYL4_IMPL_HPP_
#define MODIFIEDGRAVITYWEYL4_IMPL_HPP_

#include "DimensionDefinitions.hpp"

// AMReX includes
#include <AMReX_AmrLevel.H>

template <class theory_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void ModifiedGravityWeyl4<theory_t>::
operator()(int ix, int iy, int iz,
           const amrex::Array4<amrex::Real> &weyl_scalars,
           const amrex::Array4<amrex::Real const> &state) const
{
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);
    CCZ4Vars vars(state_cell_data);

    // we need d1 h and d2 chi, h for the geometry
    const Tensor::Sym12Rank3 d1_h =
        m_deriv.d1_sym_tensor(ix, iy, iz, state, c_h11);
    const Tensor::Rank2 h_UU = CCZ4Geometry::compute_inverse_metric(vars);
    const auto chris         = CCZ4Geometry::compute_christoffel(d1_h, h_UU);

    const Tensor::Sym12Rank2 d2_chi =
        m_deriv.d2_scalar(ix, iy, iz, state, c_chi);
    const Tensor::Sym12Sym34Rank4 d2_h =
        m_deriv.d2_sym_tensor(ix, iy, iz, state, c_h11);

    // Get the coordinates
    const Coordinates coords(amrex::IntVect{AMREX_D_DECL(ix, iy, iz)}, m_dx,
                             m_center);

    // Compute the spatial volume element
    const auto epsilon3_LUU = compute_epsilon3_LUU(vars, h_UU);

    // Compute the E and B fields (needs d1 chi, K, h, A)
    const Tensor::Rank1 d1_chi = m_deriv.d1_scalar(ix, iy, iz, state, c_chi);
    const Tensor::Rank2 d1_Gamma =
        m_deriv.d1_vector(ix, iy, iz, state, c_Gamma1);
    const Tensor::Rank1 d1_K = m_deriv.d1_scalar(ix, iy, iz, state, c_K);
    const Tensor::Sym12Rank3 d1_A =
        m_deriv.d1_sym_tensor(ix, iy, iz, state, c_A11);

    EBFields_t ebfields =
        compute_EB_fields(vars, d1_chi, d1_Gamma, d1_h, d1_K, d1_A, d2_chi,
                          d2_h, epsilon3_LUU, h_UU, chris);

    // Add in the theory terms to the E and B fields
    add_theory_EB(ebfields, ix, iy, iz, state, h_UU);

    // work out the Newman Penrose scalar
    weyl_scalar_t out = compute_Weyl4(ebfields, vars, h_UU, coords);

    // store the result
    weyl_scalars(ix, iy, iz, m_dcomp)     = out.Real;
    weyl_scalars(ix, iy, iz, m_dcomp + 1) = out.Im;
}

template <class theory_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
ModifiedGravityWeyl4<theory_t>::add_theory_EB(
    EBFields_t &ebfields, int ix, int iy, int iz,
    const amrex::Array4<const amrex::Real> &state,
    const Tensor::Rank2 &h_UU) const
{
    const Tensor::Rank2 kappa_times_Sij_TF =
        get_full_kappa_times_Sij_TF(ix, iy, iz, state, h_UU);

    // As we made the vacuum expression of Bij explicitly symmetric and Eij
    // explicitly trace-free, only Eij has theory terms
    FOR (i, j)
    {
        ebfields.E(i, j) += -0.5 * kappa_times_Sij_TF(i, j);
    }
}

template <class theory_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Rank2
ModifiedGravityWeyl4<theory_t>::get_full_kappa_times_Sij_TF(
    int ix, int iy, int iz, const amrex::Array4<const amrex::Real> &state,
    const Tensor::Rank2 &h_UU) const
{
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);
    const typename theory_t::Vars vars(state_cell_data);

    // Stress-energy sources, factor 8 pi G already applied
    const auto source =
        m_theory.compute_einstein_sources(ix, iy, iz, state, m_deriv, h_UU);

    Tensor::Rank2 S_TF = source.S;
    CCZ4Geometry::make_trace_free(S_TF, vars, h_UU);

    // Contribution of add_emtensor_rhs to the (conformal) bar-A_ij RHS is
    //   rhs.A_ij -= chi * lapse * S_TF_ij
    // and get_full_kappa_times_Sij_TF = (A_rhs_vacuum - A_rhs_full) / chi, so
    // the EM-tensor piece is lapse * S_TF_ij.
    Tensor::Rank2 out{};
    FOR (i, j)
    {
        out(i, j) = vars.lapse() * S_TF(i, j);
    }

    // TODO(4dST principal part): GRChombo's get_full_kappa_times_Sij_TF also
    // folds in the change to the bar-A_ij RHS made by
    // FourDerivScalarTensor::solve_lhs (the modified principal part). Add that
    // contribution here once solve_lhs is ported; .

    return out;
}

template <class theory_t>
void ModifiedGravityWeyl4<theory_t>::set_up(int a_state_index)
{
    const int num_ghosts = 2;

    auto &derive_lst     = amrex::AmrLevel::get_derive_lst();
    const auto &desc_lst = amrex::AmrLevel::get_desc_lst();

    // Register under the same name/vars as vacuum Weyl4 so the stock
    // WeylExtraction keeps working unchanged
    derive_lst.add(
        Weyl4::name, amrex::IndexType::TheCellType(),
        static_cast<int>(Weyl4::var_names.size()), Weyl4::var_names,
        ModifiedGravityWeyl4::compute_mf, [=](const amrex::Box &box)
        { return amrex::grow(box, num_ghosts); }, &amrex::cell_quartic_interp);

    derive_lst.addComponent(Weyl4::name, desc_lst, a_state_index, 0, NUM_VARS);
}

template <class theory_t>
void ModifiedGravityWeyl4<theory_t>::compute_mf(
    amrex::MultiFab &out_mf, int dcomp, int /*ncomp*/,
    const amrex::MultiFab &src_mf, const amrex::Geometry &geomdata,
    amrex::Real /*time*/, const int * /*bcrec*/, int /*level*/)
{
    const auto &out_arrays = out_mf.arrays();
    const auto &src_arrays = src_mf.const_arrays();

    ModifiedGravityWeyl4<theory_t> weyl4(geomdata.CellSize(0), dcomp);

    amrex::ParallelFor(
        out_mf,
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz) noexcept
        { weyl4(ix, iy, iz, out_arrays[box_no], src_arrays[box_no]); });
}

#endif /* MODIFIEDGRAVITYWEYL4_IMPL_HPP_ */
