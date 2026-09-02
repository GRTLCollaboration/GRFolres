/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef MODIFIEDCCZ4RHS_HPP_
#define MODIFIEDCCZ4RHS_HPP_

#include "CCZ4Geometry.hpp"
#include "CCZ4RHS.hpp"
#include "CCZ4Vars.hpp"
#include "FourthOrderDerivatives.hpp"
#include "StateVariables.hpp" //This files needs NUM_VARS - total number of components
#include "TensorAlgebra.hpp"

//!  Calculates RHS using CCZ4 including matter terms, and matter variable
//!  evolution
/*!
     The class calculates the RHS evolution for all the variables. It inherits
   from the CCZ4RHS class, which it uses to do the non matter evolution of
   variables. It then adds in the additional matter terms to the CCZ4 evolution
   (those including the stress energy tensor), and calculates the evolution of
   the matter variables. It does not assume a specific form of matter but is
   templated over a matter class matter_t. Please see the class ScalarField as
   an example of a matter_t. \sa CCZ4RHS(), ScalarField()
*/

struct RhoAndSi
{
    Tensor::Rank1 Si; //!< S_i = T_ia_n^a
    amrex::Real rho;           //!< rho = T_ab n^a n^b
};

struct SijTFAndS
{
    Tensor::Rank2 Sij_TF; //!< S_ij_TF = (T_ab\gamma_i^a\gamma_j^b)^TF
    amrex::Real S;                 //!< S = \gamma^ijT_ab\gamma_i^a\gamma_j^b
};

struct ScalarVectorTensor
{
    amrex::Real scalar;
    Tensor::Rank1 vector;
    Tensor::Rank2 tensor;
};


template <class theory_t, class deriv_t = FourthOrderDerivatives>
class ModifiedCCZ4RHS : public CCZ4RHS<deriv_t>
{
  public:
    // Use this alias for the same template instantiation as this class
    using CCZ4 = CCZ4RHS<deriv_t>;

    using params_t = CCZ4_params_t;

    //! Constructor of class CCZ4RHSWithMatter
    /*!
       The evolution parameters are read from the inputs. The
       default-constructed matter object supplies stress-energy sources with
       their gravitational coupling already applied.
    */
    ModifiedCCZ4RHS(amrex::Real a_dx);

    //! Add the modified gauge terms to the CCZ4 RHS
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void add_a_and_b_rhs(
        const int ix, const int iy, const int iz,
	const amrex::Array4<amrex::Real>
            &rhs_state,
        const amrex::Array4<const amrex::Real>
            &state)
        const;	    
    
    //! Add stress-energy terms to the CCZ4 RHS, including the Gamma RHS.
    /** Call this before a gauge update that derives the B-field RHS from the
     * time derivative of the conformal connection functions.
     */
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void add_emtensor_rhs(
        const int ix, const int iy, const int iz,
        const amrex::Array4<amrex::Real>
            &rhs_state, //!< the RHS data for each variable at that point.
        const amrex::Array4<const amrex::Real>
            &state) //!< the current value of the variables at the point.
        const;

    //! Add the evolution equations for the matter variables.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_theory_rhs(const int ix, const int iy, const int iz,
                   const amrex::Array4<amrex::Real> &rhs_state,
                   const amrex::Array4<const amrex::Real> &state) const;

    //! Add dissipation to the CCZ4 and matter variables.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    apply_dissipation(int ix, int iy, int iz,
                      const amrex::Array4<amrex::Real> &rhs_state,
                      const amrex::Array4<const amrex::Real> &state) const;

  protected:
    // Class members
    theory_t m_theory; //!< The matter object, e.g. a scalar field.
};

#include "MODIFIEDCCZ4RHS.impl.hpp"

#endif /* MODIFIEDCCZ4RHS_HPP_ */
