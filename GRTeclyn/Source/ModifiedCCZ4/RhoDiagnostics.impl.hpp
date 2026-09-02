/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(RHODIAGNOSTICS_HPP_)
#error "This file should only be included through RhoDiagnostics.hpp"
#endif

#ifndef RHODIAGNOSTICS_IMPL_HPP_
#define RHODIAGNOSTICS_IMPL_HPP_

#include "DimensionDefinitions.hpp"

// AMReX includes
#include <AMReX_AmrLevel.H>

template <class theory_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void RhoDiagnostics<theory_t>::operator()(
    int ix, int iy, int iz,
    const amrex::Array4<amrex::Real> &rho_diagnostics,
    const amrex::Array4<amrex::Real const> &state) const
{
    const amrex::CellData<const amrex::Real> &state_cell_data =
        state.cellData(ix, iy, iz);
    const typename theory_t::Vars vars(state_cell_data);

    const auto h_UU = CCZ4Geometry::compute_inverse_metric(vars);

    // Compute all the rho contributions for this theory
    const AllRhos all_rhos =
        m_theory.compute_all_rhos(ix, iy, iz, state, m_deriv, h_UU);

    rho_diagnostics(ix, iy, iz, m_dcomp + 0) = all_rhos.phi;
    rho_diagnostics(ix, iy, iz, m_dcomp + 1) = all_rhos.g2;
    rho_diagnostics(ix, iy, iz, m_dcomp + 2) = all_rhos.g3;
    rho_diagnostics(ix, iy, iz, m_dcomp + 3) = all_rhos.GB;
}

template <class theory_t>
void RhoDiagnostics<theory_t>::set_up(int a_state_index)
{
    const int num_ghosts = 2; // 2nd derivatives need 2 ghost cells

    auto &derive_lst     = amrex::AmrLevel::get_derive_lst();
    const auto &desc_lst = amrex::AmrLevel::get_desc_lst();

    derive_lst.add(
        name, amrex::IndexType::TheCellType(),
        static_cast<int>(var_names.size()), var_names,
        RhoDiagnostics::compute_mf, [=](const amrex::Box &box)
        { return amrex::grow(box, num_ghosts); }, &amrex::cell_quartic_interp);

    // We need all of the state variables to build the theory Vars class
    derive_lst.addComponent(name, desc_lst, a_state_index, 0, NUM_VARS);
}

template <class theory_t>
void RhoDiagnostics<theory_t>::compute_mf(amrex::MultiFab &out_mf, int dcomp,
                                          int /*ncomp*/,
                                          const amrex::MultiFab &src_mf,
                                          const amrex::Geometry &geomdata,
                                          amrex::Real /*time*/,
                                          const int * /*bcrec*/, int /*level*/)
{
    const auto &out_arrays = out_mf.arrays();
    const auto &src_arrays = src_mf.const_arrays();

    RhoDiagnostics<theory_t> rho_diagnostics(geomdata.CellSize(0), dcomp);

    amrex::ParallelFor(
        out_mf,
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz) noexcept
        {
            rho_diagnostics(ix, iy, iz, out_arrays[box_no],
                            src_arrays[box_no]);
        });
}

#endif /* RHODIAGNOSTICS_IMPL_HPP_ */
