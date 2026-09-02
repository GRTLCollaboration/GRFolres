/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef FOURDERIVSCALARTENSOR_HPP_
#define FOURDERIVSCALARTENSOR_HPP_

#include "CCZ4Geometry.hpp"
#include "DefaultPotential.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRParmParse.hpp"
#include "ScalarFieldVars.hpp"
#include "StateVariables.hpp" //This files needs NUM_VARS, total num of components
#include "TensorAlgebra.hpp"

//!  Calculates the matter type specific elements such as the EMTensor and
//   matter evolution
/*!
     This class is an example of a matter_t object which calculates the
     matter type specific elements for the RHS update and the evaluation
     of the constraints. This includes the Energy Momentum Tensor, and
     the matter evolution terms. In this case, a scalar field,
     the matter elements are phi and (minus) its conjugate momentum, Pi.
     It is templated over a potential function potential_t which the
     user must specify in a class, although a default is provided which
     sets dVdphi and V_of_phi to zero.
     It assumes minimal coupling of the field to gravity.
     Matter classes used by the diagnostic callbacks must be default
     constructible. Their default constructors should therefore read all
     runtime parameters needed to construct a fully configured matter object,
     as ScalarField and its potential do here.
     \sa MatterCCZ4(), ConstraintsMatter()
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

    //!  Constructor of class FourDerivScalarTensor, inputs are the matter parameters.
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

    //! The function which calculates the EM Tensor, given the vars and
    //! derivatives, including the potential
    [[nodiscard]]
    AMREX_GPU_DEVICE emtensor_t compute_emtensor(
        const int ix, const int iy, const int iz, //!< grid indicies
        const amrex::Array4<const amrex::Real>
            &state,             //!< the current value of state variables
        const deriv_t &a_deriv, //!< the object that calculates the derivative
        const Tensor::Rank2 &h_UU) //!< the inverse metric (raised indices)
        const;

    //! Calculate the stress-energy sources including the factor 8 pi G.
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE einstein_sources_t
    compute_einstein_sources(int ix, int iy, int iz,
                             const amrex::Array4<const amrex::Real> &state,
                             const deriv_t &a_deriv,
                             const Tensor::Rank2 &h_UU) const;

    // ! The function which adds in the RHS for the matter field vars,
    // ! including the potential
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void add_theory_rhs(
        const int ix, const int iy, const int iz, //!< grid indicies
        const amrex::Array4<amrex::Real>
            &rhs_state, //!< the next value of state variables (rhs update)
        const amrex::Array4<const amrex::Real>
            &state, //!< the current value of state variables
        const deriv_t &a_deriv)
        const; //!< the object for calculating derivatives

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void compute_lhs(
        const int ix, const int iy, const int iz,
	const amrex::Array4<const amrex::Real> &state,
	const deriv_t &a_deriv,
	const int matrix_dim, amrex::Real *LHS)
	const;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void solve_lhs(
	const int ix, const int iy, cons int iz,
	const amrex::Array4<amrex::Real> &rhs_state,
	const amrex::Array4<const amrex::Real> &state,
	const deriv_t &a_deriv)
        const;

};

#include "FourDerivScalarTensor.impl.hpp"

#endif /* FOURDERIVSCALARTENSOR_HPP_ */
