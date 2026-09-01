/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(FOURDERIVSCALARTENSOR_HPP_)
#error "This file should only be included through FourDerivScalarTensor.hpp"
#endif

#ifndef FOURDERIVSCALARTENSOR_IMPL_HPP_
#define FOURDERIVSCALARTENSOR_IMPL_HPP_

// Calculate the stress energy tensor elements
template <class potential_t, class deriv_t>
AMREX_GPU_DEVICE emtensor_t FourDerivScalarTensor<potential_t, deriv_t>::compute_emtensor(
    const int ix, const int iy, const int iz,
    const amrex::Array4<const amrex::Real> &state, const deriv_t &a_deriv,
    const Tensor::Rank2 &h_UU) const
{
    emtensor_t out;

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
    // S = T_ij
    FOR (i, j)
    {
        out.S(i, j) = -0.5 * vars.h(i, j) * Vt / vars.chi() +
                      d1_phi(i) * d1_phi(j) -
                      vars.h(i, j) * V_of_phi / vars.chi();
    }

    // rho = n^a n^b T_ab
    out.rho = vars.Pi() * vars.Pi() + 0.5 * Vt + V_of_phi;

    // trS = Tr_S_ij
    out.trS = vars.chi() * TensorAlgebra::compute_trace(out.S, h_UU);

    //    j_i (note lower index) = - n^a T_ai
    FOR (i)
    {
        out.j(i) = -d1_phi(i) * vars.Pi();
    }

    return out;
}

template <class potential_t, class deriv_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE einstein_sources_t
FourDerivScalarTensor<coupling_and_potential_t, deriv_t>::compute_einstein_sources(
    int ix, int iy, int iz, const amrex::Array4<const amrex::Real> &state,
    const deriv_t &a_deriv, const Tensor::Rank2 &h_UU) const
{
    const emtensor_t emtensor =
        compute_emtensor(ix, iy, iz, state, a_deriv, h_UU);
    const amrex::Real coupling = 8.0 * M_PI * m_G_Newton;

    einstein_sources_t out;
    out.rho = coupling * emtensor.rho;
    out.trS = coupling * emtensor.trS;
    FOR (i)
    {
        out.j(i) = coupling * emtensor.j(i);
    }
    FOR (i, j)
    {
        out.S(i, j) = coupling * emtensor.S(i, j);
    }
    return out;
}

// Adds in the RHS for the theory vars
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

#endif /* FOURDERIVSCALARTENSOR_IMPL_HPP_ */
