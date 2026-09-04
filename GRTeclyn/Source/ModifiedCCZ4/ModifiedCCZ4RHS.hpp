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
#include "GRParmParse.hpp"
#include "StateVariables.hpp" //This files needs NUM_VARS - total number of components
#include "TensorAlgebra.hpp"

//!  Calculates RHS using CCZ4 (in the modified gauge) including theory terms, 
//!  and matter variable evolution
/*!
     The class calculates the RHS evolution for all the variables. It inherits
   from the CCZ4RHS class, which it uses to do the vacuum GR evolution of
   variables. It then adds in the additional terms from the modified CCZ4 gauge.
   Next, it adds the source terms of the metric sector corresponding to the RHS
   of Einstein's equations (those including via an effective stress energy tensor),
   which may include solving a linear system. Finally it calculates the evolution 
   of the non-metric variables, if any. It does not assume a specific theory but is
   templated over a theory class theory_t. Please see the class FourDerivScalarTensor 
   as an example of a theory_t. \sa CCZ4RHS(), FourDerivScalarTensor()
*/

struct RhoAndJ
{
    Tensor::Rank1 j; //!< S_i = T_ia_n^a
    amrex::Real rho;           //!< rho = T_ab n^a n^b
};

struct S_TFAndTrS
{
    Tensor::Rank2 S_TF; //!< S_ij_TF = (T_ab\gamma_i^a\gamma_j^b)^TF
    amrex::Real trS;                 //!< S = \gamma^ijT_ab\gamma_i^a\gamma_j^b
};

struct ScalarVectorTensor
{
    amrex::Real scalar;
    Tensor::Rank1 vector;
    Tensor::Rank2 tensor;
};

struct einstein_sources_TF
{
    amrex::Real rho;
    Tensor::Rank1 j;
    Tensor::Rank2 S_TF;
    amrex::Real trS;
};

// Individual contributions to rho = T_ab n^a n^b, stored as diagnostics.
// Note: additional components may be needed for some theories.
struct AllRhos
{
    amrex::Real phi; //!< minimally-coupled scalar contribution
    amrex::Real g2;  //!< k_2(phi) X^2 (k-essence) contribution
    amrex::Real g3;  //!< g_3(phi) box(phi) X contribution
    amrex::Real GB;  //!< Gauss-Bonnet contribution
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

    static void check_params()
    {
        // To be added checker for mod_a and mod_b params
    }

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

    //! Solve the LHS (if any)
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    solve_lhs(const int ix, const int iy, const int iz,
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
    amrex::Real m_mod_a;
    amrex::Real m_mod_b;
};

#include "ModifiedCCZ4RHS.impl.hpp"

#endif /* MODIFIEDCCZ4RHS_HPP_ */
