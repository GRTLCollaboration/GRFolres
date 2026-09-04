/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef FOURDERIVSCALARTENSOR_HPP_
#define FOURDERIVSCALARTENSOR_HPP_

#include "CCZ4Geometry.hpp"
#include "DefaultCouplingAndPotential.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRParmParse.hpp"
#include "LinearSolver.hpp"
#include "ModifiedCCZ4RHS.hpp"
#include "ScalarFieldVars.hpp"
#include "StateVariables.hpp" //This files needs NUM_VARS, total num of components
#include "TensorAlgebra.hpp"

//!  Calculates the matter type specific elements such as the EMTensor and
//   matter evolution
/*!
     This class is an example of a theory_t object which calculates the
     theory type specific elements for the RHS update and the evaluation
     of the constraints. This includes the source terms of the, and
     the theory evolution terms. In this case, a scalar field,
     the theory elements are phi and (minus) its conjugate momentum, Pi.
     It is templated over a coupling and potential function coupling_and_potential_t 
     which the user must specify in a class, although a default is provided which
     sets dfdphi, d2fdphi2, dVdphi, V_of_phi, g2 and dg2dphi to zero.
     It assumes minimal coupling of the field to gravity.
     Matter classes used by the diagnostic callbacks must be default
     constructible. Their default constructors should therefore read all
     runtime parameters needed to construct a fully configured matter object,
     as FourDerivScalarTensor and its coupling_and_potential do here.
     \sa ModifiedCCZ4(), ConstraintsMatter()
*/
template <class coupling_and_potential_t = DefaultCouplingAndPotential,
          class deriv_t     = FourthOrderDerivatives>
class FourDerivScalarTensor
{
  protected:
    coupling_and_potential_t m_coupling_and_potential;
    //! The local copy of the potential
    amrex::Real m_G_Newton{1.0};

  public:

    struct params_t
    {
        amrex::Real G_Newton{1.0};

        static void check_params()
        {
            GRParmParse four_deriv_scalar_tensor_pp("four_deriv_scalar_tensor");
            amrex::Real G_Newton{1.0};
            four_deriv_scalar_tensor_pp.queryAdd("G_Newton", G_Newton);
            if (G_Newton < 0.0)
            {
                four_deriv_scalar_tensor_pp.error("G_Newton", "must be >= 0.0");
            }
        }

        void fill_params()
        {
            GRParmParse four_deriv_scalar_tensor_pp("four_deriv_scalar_tensor");
            four_deriv_scalar_tensor_pp.query("G_Newton", G_Newton);
        }
    };

    //!  Constructor of class FourDerivScalarTensor, inputs are the theory parameters.
    FourDerivScalarTensor()
    {
        params_t params;
        params.fill_params();
        m_G_Newton = params.G_Newton;
    }

    AMREX_FORCE_INLINE explicit FourDerivScalarTensor(coupling_and_potential_t a_coupling_and_potential)
        : m_coupling_and_potential(a_coupling_and_potential)
    {
        params_t params;
        params.fill_params();
        m_G_Newton = params.G_Newton;
    }

    using Vars = ScalarFieldVars;

    //! The function which computes:
    //! M_{ij} = R_{ij} + KK_{ij} - K_{ik}K_j^{~k}
    //! N_i = D^jK_{ij}-D_iK (GR momentum constraint)
    //! M = \gamma^{ij}M_{ij} (GR Hamiltonian constraint)
    [[nodiscard]]
    AMREX_GPU_DEVICE ScalarVectorTensor compute_M_Ni_and_Mij(
        const int ix, const int iy, const int iz, //!< grid indicies
        const amrex::Array4<const amrex::Real>
            &state,             //!< the current value of state variables
        const deriv_t &a_deriv, //!< the object that calculates the derivative
        const Tensor::Rank2 &h_UU) //!< the inverse metric (raised indices)
        const;

    //! The function which computes the decomposition of Omega_{\mu\nu},
    //! \Omega_{\mu\nu} = \nabla_{\mu}\nabla_{\nu}f(\phi)
    [[nodiscard]]
    AMREX_GPU_DEVICE ScalarVectorTensor compute_Omega_munu(
        const int ix, const int iy, const int iz, //!< grid indicies
        const amrex::Array4<const amrex::Real>
            &state,             //!< the current value of state variables
        const deriv_t &a_deriv, //!< the object that calculates the derivative
        const Tensor::Rank2 &h_UU) //!< the inverse metric (raised indices)
        const;

    //! The function which calculates the rho and j components of the effective 
    //! EM Tensor, given the vars and derivatives, including the potential
    [[nodiscard]]
    AMREX_GPU_DEVICE RhoAndJ compute_rho_and_j(
        const int ix, const int iy, const int iz, //!< grid indicies
        const amrex::Array4<const amrex::Real>
            &state,             //!< the current value of state variables
        const deriv_t &a_deriv, //!< the object that calculates the derivative
        const Tensor::Rank2 &h_UU) //!< the inverse metric (raised indices)
        const;

    //! The function which calculates the S_TF and trS components of the effective
    //! EM Tensor, given the vars and derivatives, including the potential
    [[nodiscard]]
    AMREX_GPU_DEVICE S_TFAndTrS compute_S_TF_and_trS(
        const int ix, const int iy, const int iz, //!< grid indicies
        const amrex::Array4<const amrex::Real>
            &state,             //!< the current value of state variables
        const deriv_t &a_deriv, //!< the object that calculates the derivative
        const Tensor::Rank2 &h_UU) //!< the inverse metric (raised indices)
        const;

    //! Calculate the stress-energy sources including the factor 8 pi G.
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE einstein_sources_TF
    compute_einstein_sources(int ix, int iy, int iz,
                             const amrex::Array4<const amrex::Real> &state,
                             const deriv_t &a_deriv,
                             const Tensor::Rank2 &h_UU) const;

    //! The function which adds in the RHS for the theory field vars,
    //! including the potential
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void add_theory_rhs(
        const int ix, const int iy, const int iz, //!< grid indices
        const amrex::Array4<amrex::Real>
            &rhs_state, //!< the next value of state variables (rhs update)
        const amrex::Array4<const amrex::Real>
            &state, //!< the current value of state variables
        const deriv_t &a_deriv)
        const; //!< the object for calculating derivatives

    //! The function which computes the LHS matrix for some of the vars, which
    //! are A, K and Pi for this 4dST example
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void compute_lhs(
        const int ix, const int iy, const int iz, //!< grid indices
	const amrex::Array4<const amrex::Real> 
	    &state, //!< the current value of state variables
	const deriv_t &a_deriv, //!< the object for calculating derivates
	const int matrix_dim, //!< the dimension of the LHS matrix
	amrex::Real *LHS)
	const; //!< the LHS matrix itself

    //! The function which solves the linear system using the LHS matrix
    //! computed in compute_lhs and the RHS calculated before
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void solve_lhs(
	const int ix, const int iy, const int iz, //!< grid indices
	const amrex::Array4<amrex::Real> 
	    &rhs_state, //!< the next value of state variables (rhs update)
	const amrex::Array4<const amrex::Real> 
	    &state, //!< the current value of state variables
	const deriv_t &a_deriv)
        const; //!< the object for calculating derivatives

};

#include "FourDerivScalarTensor.impl.hpp"

#endif /* FOURDERIVSCALARTENSOR_HPP_ */
