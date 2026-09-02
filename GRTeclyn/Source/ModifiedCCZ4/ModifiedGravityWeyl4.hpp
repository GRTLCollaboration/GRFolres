/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef MODIFIEDGRAVITYWEYL4_HPP_
#define MODIFIEDGRAVITYWEYL4_HPP_

#include "CCZ4Geometry.hpp"
#include "Weyl4.hpp"

//!  Calculates the Weyl4 scalar for spacetimes in a modified-gravity theory
/*!
   This class calculates the real and imaginary parts of the Weyl4 scalar. It
   inherits from the Weyl4 class (vacuum expression, decomposed into electric
   and magnetic parts) and adds in the theory contribution to the electric
   part. It is the modified-gravity analogue of GRTeclyn's Weyl4WithMatter.

   For a modified-gravity theory the trace-free stress that enters the electric
   field is not simply 8 pi G S_ij^TF: it also picks up the change to the
   \bar A_ij evolution coming from the modified principal part (the per-cell
   linear solve, FourDerivScalarTensor::solve_lhs). See
   get_full_kappa_times_Sij_TF() and GRChombo's
   ModifiedCCZ4RHS::get_full_kappa_times_Sij_TF().

*/
template <class theory_t> class ModifiedGravityWeyl4 : public Weyl4
{
  public:
    //! Constructor
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    ModifiedGravityWeyl4(amrex::Real a_dx, int a_dcomp)
        : Weyl4(a_dx, a_dcomp)
    {
    }

    //! The compute member which calculates the wave quantities at each point on
    //! the grid
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &weyl_scalars,
               const amrex::Array4<amrex::Real const> &state) const;
    // NOLINTEND(bugprone-easily-swappable-parameters)

    static void set_up(int a_state_index);

    // Has signature of DeriveFuncMF so that it can be stored in the derive_lst
    static void compute_mf(amrex::MultiFab &out_mf, int dcomp, int ncomp,
                           const amrex::MultiFab &src_mf,
                           const amrex::Geometry &geomdata,
                           amrex::Real /*time*/, const int * /*bcrec*/,
                           int /*level*/);

  protected:
    theory_t m_theory; //!< The theory object, e.g. 4dST

    //! Add the theory terms to the electric and magnetic parts
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_theory_EB(EBFields_t &eb_fields, int ix, int iy, int iz,
                  const amrex::Array4<const amrex::Real> &state,
                  const Tensor::Rank2 &h_UU) const;
    // NOLINTEND(bugprone-easily-swappable-parameters)

    //! Full \kappa S_{ij}^{TF} that enters the electric part, i.e. the change
    //! the theory makes to the (conformal) \bar A_ij RHS, divided by chi.
    //! Port of GRChombo ModifiedCCZ4RHS::get_full_kappa_times_Sij_TF.
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE Tensor::Rank2
    get_full_kappa_times_Sij_TF(
        int ix, int iy, int iz,
        const amrex::Array4<const amrex::Real> &state,
        const Tensor::Rank2 &h_UU) const;
    // NOLINTEND(bugprone-easily-swappable-parameters)
};

#include "ModifiedGravityWeyl4.impl.hpp"

#endif /* MODIFIEDGRAVITYWEYL4_HPP_ */
