/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef MODIFIEDGRAVITYCONSTRAINTS_HPP_
#define MODIFIEDGRAVITYCONSTRAINTS_HPP_

#include "CCZ4Geometry.hpp"
#include "Constraints.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"

#include <array>

//!  Calculates the Hamiltonian and Momentum constraints with the theory fields
/*!
     The class calculates the Hamiltonian and Momentum constraints at each point
   in a box. It inherits from the Constraints class which calculates the
   constraints without the theory terms. It adds in the theory terms for a given
   theory class theory_t, which must provide the stress-energy sources with
   their gravitational coupling already applied (compute_einstein_sources).
   For an example of a theory_t class see FourDerivScalarTensor.
   \sa Constraints(), ConstraintsWithMatter(), FourDerivScalarTensor()

   This is the modified-gravity analogue of GRTeclyn's ConstraintsWithMatter:
   the Hamiltonian and Momentum constraints are geometric (gauge independent),
   so the only difference from vacuum is the stress-energy source, which for a
   modified-gravity theory carries the extra (e.g. Gauss-Bonnet) contributions.
*/
template <class theory_t>
class ModifiedGravityConstraints : public Constraints
{
  public:
    //! Constructor of class ModifiedGravityConstraints
    /*!
        Can specify the vars of the constraint vars instead of using the
        hardcoded ones.
    */
    ModifiedGravityConstraints(
        amrex::Real dx, int a_c_Ham, const Interval &a_c_Moms,
        int a_c_Ham_abs_terms              = -1,
        const Interval &a_c_Moms_abs_terms = Interval());

    //! The compute member which calculates the constraints at each point in the
    //! box
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &constraints,
               const amrex::Array4<amrex::Real const> &state) const;

    //! Adds the constraints to the derive list. Call in variableSetUp().
    static void set_up(int a_state_index, bool a_calc_mom_norm = false);

    // Has signature of DeriveFuncMF so that it can be stored in the derive_lst
    static void compute_mf(amrex::MultiFab &out_mf, int dcomp, int ncomp,
                           const amrex::MultiFab &src_mf,
                           const amrex::Geometry &geomdata,
                           amrex::Real /*time*/, const int * /*bcrec*/,
                           int /*level*/);

  protected:
    theory_t m_theory; //!< The theory object, e.g. 4dST
};

#include "ModifiedGravityConstraints.impl.hpp"

#endif /* MODIFIEDGRAVITYCONSTRAINTS_HPP_ */
