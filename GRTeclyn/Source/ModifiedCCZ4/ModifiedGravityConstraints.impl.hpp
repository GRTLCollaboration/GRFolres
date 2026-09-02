/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(MODIFIEDGRAVITYCONSTRAINTS_HPP_)
#error                                                                         \
    "This file should only be included through ModifiedGravityConstraints.hpp"
#endif

#ifndef MODIFIEDGRAVITYCONSTRAINTS_IMPL_HPP_
#define MODIFIEDGRAVITYCONSTRAINTS_IMPL_HPP_

#include "DimensionDefinitions.hpp"

// AMReX includes
#include <AMReX_AmrLevel.H>

template <class theory_t>
ModifiedGravityConstraints<theory_t>::ModifiedGravityConstraints(
    amrex::Real dx, int a_c_Ham, const Interval &a_c_Moms,
    int a_c_Ham_abs_terms /* defaulted*/,
    const Interval &a_c_Moms_abs_terms /*defaulted*/)
    : Constraints(dx, a_c_Ham, a_c_Moms, a_c_Ham_abs_terms, a_c_Moms_abs_terms,
                  0.0 /*No cosmological constant*/)
{
}

template <class theory_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
ModifiedGravityConstraints<theory_t>::operator()(
    int ix, int iy, int iz, const amrex::Array4<amrex::Real> &constraints,
    const amrex::Array4<amrex::Real const> &state) const
{
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);
    typename theory_t::Vars vars(state_cell_data);

    // we need d1 chi, K, h, A and d2 chi, h
    const auto d1_chi   = m_deriv.d1_scalar(ix, iy, iz, state, c_chi);
    const auto d1_Gamma = m_deriv.d1_vector(ix, iy, iz, state, c_Gamma1);
    const auto d1_K     = m_deriv.d1_scalar(ix, iy, iz, state, c_K);
    const auto d1_A     = m_deriv.d1_sym_tensor(ix, iy, iz, state, c_A11);
    const auto d1_h     = m_deriv.d1_sym_tensor(ix, iy, iz, state, c_h11);

    const auto d2_chi = m_deriv.d2_scalar(ix, iy, iz, state, c_chi);
    const auto d2_h   = m_deriv.d2_sym_tensor(ix, iy, iz, state, c_h11);

    // Inverse metric and Christoffel symbol
    const auto h_UU  = CCZ4Geometry::compute_inverse_metric(vars);
    const auto chris = CCZ4Geometry::compute_christoffel(d1_h, h_UU);

    // Get the vacuum (geometric) terms for the constraints
    constraints_t out = constraint_equations(vars, d1_chi, d1_Gamma, d1_h, d1_K,
                                             d1_A, d2_chi, d2_h, h_UU, chris);

    // Stress-energy sources, with the factor 8 pi G already applied. For a
    // modified-gravity theory these also carry the extra (e.g. Gauss-Bonnet)
    // contributions to rho and j_i.
    const auto source =
        m_theory.compute_einstein_sources(ix, iy, iz, state, m_deriv, h_UU);

    // Hamiltonian constraint
    if (m_c_Ham >= 0 || m_c_Ham_abs_terms >= 0)
    {
        out.Ham           += -2.0 * source.rho;
        out.Ham_abs_terms += 2.0 * std::abs(source.rho);
    }

    // Momentum constraints
    if (m_c_Moms.size() > 0 || m_c_Moms_abs_terms.size() > 0)
    {
        FOR (i)
        {
            out.Mom(i)           += -source.j(i);
            out.Mom_abs_terms(i) += std::abs(source.j(i));
        }
    }

    // Write the constraints into the output FArrayBox
    store_vars(out, constraints.cellData(ix, iy, iz));
}

template <class theory_t>
void ModifiedGravityConstraints<theory_t>::set_up(int a_state_index,
                                                  bool a_calc_mom_norm)
{
    const int num_ghosts = 2; // no advection terms so only need 2 ghost cells

    auto &derive_lst     = amrex::AmrLevel::get_derive_lst();
    const auto &desc_lst = amrex::AmrLevel::get_desc_lst();

    const auto &comp_names = (a_calc_mom_norm) ? Constraints::var_names_norm
                                               : Constraints::var_names;

    // Add the constraints to the derive list (under the same name as the
    // vacuum Constraints so existing plot_vars / norm helpers keep working)
    derive_lst.add(
        Constraints::name, amrex::IndexType::TheCellType(),
        static_cast<int>(comp_names.size()), comp_names,
        ModifiedGravityConstraints::compute_mf, [=](const amrex::Box &box)
        { return amrex::grow(box, num_ghosts); }, &amrex::cell_quartic_interp);

    // We need all of the state variables to build the theory Vars class
    derive_lst.addComponent(Constraints::name, desc_lst, a_state_index, 0,
                            NUM_VARS);
}

template <class theory_t>
void ModifiedGravityConstraints<theory_t>::compute_mf(
    amrex::MultiFab &out_mf, int dcomp, int ncomp,
    const amrex::MultiFab &src_mf, const amrex::Geometry &geomdata,
    amrex::Real /*time*/, const int * /*bcrec*/, int /*level*/)
{
    const auto &out_arrays = out_mf.arrays();
    const auto &src_arrays = src_mf.const_arrays();

    const amrex::Real dx = geomdata.CellSize(0);
    const int iham       = dcomp; // Ham
    const Interval imom =
        Interval(dcomp + 1, dcomp + AMREX_SPACEDIM); // Mom1, Mom2, Mom3

    AMREX_ALWAYS_ASSERT(ncomp == (1 + AMREX_SPACEDIM));

    ModifiedGravityConstraints<theory_t> constraints(dx, iham, imom);

    amrex::ParallelFor(
        out_mf,
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz) noexcept
        { constraints(ix, iy, iz, out_arrays[box_no], src_arrays[box_no]); });
}

#endif /* MODIFIEDGRAVITYCONSTRAINTS_IMPL_HPP_ */
