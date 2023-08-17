/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef TestField4dST_HPP_
#define TestField4dST_HPP_

#include "CCZ4Geometry.hpp"
#include "Coordinates.hpp"
#include "DefaultCouplingAndPotential.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "LinearSolver.hpp"
#include "ModifiedCCZ4RHS.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS, total num of components
#include "VarsTools.hpp"

//!  Calculates the theory type specific elements such as the EMTensor and
//   theory evolution
/*!
     This class is an example of a theory_t object which calculates the
     theory type specific elements for the RHS update and the evaluation
     of the constraints. This includes the Energy Momentum Tensor, and
     the theory evolution terms. In this case, a scalar field coupled to
     the Gauss-Bonnet invariant, the theory elements are phi and its
     conjugate momentum, Pi.
     It is templated over a coupling function coupling_t which the
     user must specify in a class, although a default is provided which
     sets g2, dfdphi and df2dphi2 to zero.
     It assumes a scalar field, coming from the action
     S=\int d^4x\left(R/16\pi G + X + V(\phi) + g_2(\phi) X^2
     + V(\phi) + f(\phi)/4{\mathcal R}_{GB}\right)
     where X = - 1/2\nabla_{\mu}\phi\nabla^{\mu}\phi
     \sa ModifiedCCZ4RHS(), ModifiedGravityConstraints()
*/

template <class coupling_and_potential_t = DefaultCouplingAndPotential>
class TestField4dST
{
  protected:
    //! The local copy of the coupling
    coupling_and_potential_t my_coupling_and_potential;


  public:
    //!  Constructor of class FourDerivScalarTensor, inputs are the theory
    //!  parameters.
    TestField4dST(
        const coupling_and_potential_t a_coupling_and_potential)
        : my_coupling_and_potential(a_coupling_and_potential)
    {
    }

    //! Structure containing the rhs variables for the theory fields
    template <class data_t> struct Vars
    {
        data_t phi; // the scalar field
        data_t Pi;  // conjugate momentum of the scalar field

        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_phi, phi);
            VarsTools::define_enum_mapping(mapping_function, c_Pi, Pi);
        }
    };

    //! Structure containing the rhs variables for the theory fields requiring
    //!  2nd derivs
    template <class data_t> struct Diff2Vars
    {
        data_t phi;

        /// Defines the mapping between members of Vars and Chombo grid
        ///  variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_phi, phi);
        }
    };

    //! The function which computes:
    //! M_{ij} = R_{ij} + KK_{ij} - K_{ik}K_j^{~k}
    //! N_i = D^jK_{ij}-D_iK (GR momentum constraint)
    //! M = \gamma^{ij}M_{ij} (GR Hamiltonian constraint)
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    ScalarVectorTensor<data_t> compute_M_Ni_and_Mij(
        const vars_t<data_t> &vars, //!< the value of the variables
        const vars_t<Tensor<1, data_t>>
            &d1, //!< the value of the 1st derivatives
        const diff2_vars_t<Tensor<2, data_t>> &d2)
        const; //!< the value of the 2nd derivatives



    //! The function which calculates rho and Si, given the vars and
    //! derivatives, for the given coupling function
    //! NOTE: this is computed in a separate function from the other
    //! elements, namely S and Sij, so that we can call directly
    //! this function to calculate the constraint without needing to
    //! specify the gauge variables
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    RhoAndSi<data_t> compute_rho_and_Si(
        const vars_t<data_t> &vars,          //!< the value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>>
            &d2, //!< the value of the 2nd derivs
        const Coordinates<data_t> &coords)
        const; //!< the value of the coordinates

    //! The function which calculates SijTF and S, given the vars and
    //! derivatives, for the given coupling function
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    SijTFAndS<data_t> compute_Sij_TF_and_S(
        const vars_t<data_t> &vars,          //!< the value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>>
            &d2,                     //!< the value of the 2nd derivs
        const vars_t<data_t> &advec, //!< the value of the advection terms
        const Coordinates<data_t> &coords)
        const; //!< the value of the coordinates

    //! The function which adds in the RHS for the theory field vars,
    //! given the coupling function
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t,
              template <typename> class rhs_vars_t>
    void add_theory_rhs(
        rhs_vars_t<data_t> &rhs,    //!< the value of the RHS for all vars
        const vars_t<data_t> &vars, //!< the value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< the value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>>
            &d2,                     //!< the value of the 2nd derivs
        const vars_t<data_t> &advec, //!< the value of the advection terms
        const Coordinates<data_t> &coords)
        const; //!< the value of the coordinates


    //! The function which solves the lhs of Pi  from rhs of Aij and K
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t,
              template <typename> class rhs_vars_t>
    void solve_lhs(
        rhs_vars_t<data_t> &rhs,             //!< value of the RHS for all vars
        const vars_t<data_t> &vars,          //!< value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>> &d2, //!< value of the 2nd derivs
        const vars_t<data_t> &advec, //!< the value of the advection terms
        const Coordinates<data_t> &coords)
        const; //!< the value of the coordinates
};

#include "TestField4dST.impl.hpp"

#endif /* TestField4dST_HPP_ */
