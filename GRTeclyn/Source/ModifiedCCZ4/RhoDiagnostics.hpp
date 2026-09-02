/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef RHODIAGNOSTICS_HPP_
#define RHODIAGNOSTICS_HPP_

#include "CCZ4Geometry.hpp"
#include "FourthOrderDerivatives.hpp"
#include "ModifiedCCZ4RHS.hpp" // for the AllRhos struct
#include "Tensor.hpp"

// AMReX includes
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <array>
#include <string>

//! Calculates all the individual contributions to rho = T_ab n^a n^b for a
//! given modified-gravity theory, stored as an AMReX derived record so they
//! can be written to plotfiles on demand.
/*!
   The theory class theory_t must be default constructible (it reads its own
   runtime parameters) and must provide
       AllRhos compute_all_rhos(int ix, int iy, int iz,
                                const amrex::Array4<const amrex::Real> &state,
                                const deriv_t &a_deriv,
                                const Tensor::Rank2 &h_UU) const;
   For an example of a theory_t class see FourDerivScalarTensor.
*/
template <class theory_t> class RhoDiagnostics
{
  public:
    /// derive record name
    static inline const std::string name = "rho_diagnostics";

    /// Variable names
    static inline const amrex::Vector<std::string> var_names = {
        "rho_phi", "rho_g2", "rho_g3", "rho_GB"};

    //! Constructor
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    RhoDiagnostics(amrex::Real a_dx, int a_dcomp)
        : m_deriv(a_dx), m_dcomp(a_dcomp)
    {
    }

    //! Computes the rho diagnostics in a cell
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &rho_diagnostics,
               const amrex::Array4<amrex::Real const> &state) const;

    //! Adds the rho diagnostics to the derive list. Call in variableSetUp().
    static void set_up(int a_state_index);

    // Has signature of DeriveFuncMF so that it can be stored in the derive_lst
    static void compute_mf(amrex::MultiFab &out_mf, int dcomp, int ncomp,
                           const amrex::MultiFab &src_mf,
                           const amrex::Geometry &geomdata,
                           amrex::Real /*time*/, const int * /*bcrec*/,
                           int /*level*/);

  protected:
    theory_t m_theory;              //!< The theory object, e.g. 4dST
    FourthOrderDerivatives m_deriv; //!< for calculating derivs of vars
    int m_dcomp; //!< first component in the MultiFab to store into
};

#include "RhoDiagnostics.impl.hpp"

#endif /* RHODIAGNOSTICS_HPP_ */
