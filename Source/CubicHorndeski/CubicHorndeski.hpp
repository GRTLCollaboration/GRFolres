/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CUBICHORNDESKI_HPP_
#define CUBICHORNDESKI_HPP_

#include "CCZ4Geometry.hpp"
#include "Coordinates.hpp"
#include "DefaultCouplingAndPotential.hpp"
#include "FourthOrderDerivatives.hpp"
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
     the theory evolution terms. In this case, the quadratic and cubic
     Horneski terms, the theory elements are phi and its conjugate
     momentum, Pi.
     It is templated over a coupling and potential function
     coupling_and_potential_t which the user must specify in a class,
     although a default is provided which sets G_2(\phi, X) = g_2X^2
     and G_3(\phi, X) = g_3X.
     It assumes a scalar field, coming from the action
     S=\int d^4x\left(R/(16\pi G) + X + G_2(\phi, X)
     + G_3(\phi, X)\Box\phi\right),
     where X=-1/2\nabla_{\mu}\phi\nabla^{\mu}\phi and the potential
     of the scalar field is considered inside G_2(\phi) function
     \sa ModifiedCCZ4RHS(), ModifiedGravityConstraints()
*/

template <class coupling_and_potential_t = DefaultCouplingAndPotential>
class CubicHorndeski
{
  public:
    template <class data_t> struct expressions
    {
        data_t lie_deriv_Pi_no_lapse;

        data_t g2;
        data_t dg2_dX, dg3_dX;
        data_t dg2_dphi, dg3_dphi;
        data_t d2g2_dXX, d2g3_dXX;
        data_t d2g2_dXphi, d2g3_dXphi;
        data_t d2g3_dphiphi;

        data_t dcommon;

        data_t Pi2, X, Xplus;

        data_t tau;
        Tensor<1, data_t> tau_i;
        Tensor<2, data_t> tau_ij;

        data_t tau_i_dot_dphi;
        Tensor<1, data_t> tau_ij_dot_dphi;
        data_t tau_ij_dot_dphi2;

        data_t exprA, exprB;
        data_t ders1, ders2;
    };

    CubicHorndeski(const coupling_and_potential_t a_coupling_and_potential)
        : my_coupling_and_potential(a_coupling_and_potential)
    {
    }

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

    //! The function which calculates rho and Si, given the vars and
    //! derivatives, for the given coupling function
    //! NOTE: this is computed in a separate function from the other
    //! elements, namely S and Sij, so that we can call directly
    //! this function to calculate the constraint without needing to
    //! specify the gauge variables
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    RhoAndSi<data_t>
    compute_rho_and_Si(const vars_t<data_t> &vars,
                       const vars_t<Tensor<1, data_t>> &d1,
                       const diff2_vars_t<Tensor<2, data_t>> &d2,
                       const Coordinates<data_t> &coords) const;

    //! The function which calculates SijTF and S, given the vars and
    //! derivatives, for the given coupling function
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    SijTFAndS<data_t> compute_Sij_TF_and_S(
        const vars_t<data_t> &vars, const vars_t<Tensor<1, data_t>> &d1,
        const diff2_vars_t<Tensor<2, data_t>> &d2, const vars_t<data_t> &advec,
        const Coordinates<data_t> &coords) const;

    //! The function which adds in the RHS for the theory field vars,
    //! including the potential
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t,
              template <typename> class rhs_vars_t>
    void add_theory_rhs(
        rhs_vars_t<data_t> &total_rhs,       //!< value of the RHS for all vars
        const vars_t<data_t> &vars,          //!< value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>> &d2, //!< value of the 2nd derivs
        const vars_t<data_t> &advec,
        const Coordinates<data_t> &coords) //!< the value of the advection terms
        const;                             //!< the value of the coordinates

    //! The function which solves the linear system using the LHS matrix
    //! for those theories where necessary (not Cubic Horndeski)
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t,
              template <typename> class rhs_vars_t>
    void solve_lhs(
        rhs_vars_t<data_t> &total_rhs,       //!< value of the RHS for all vars
        const vars_t<data_t> &vars,          //!< value of the variables
        const vars_t<Tensor<1, data_t>> &d1, //!< value of the 1st derivs
        const diff2_vars_t<Tensor<2, data_t>> &d2, //!< value of the 2nd derivs
        const vars_t<data_t> &advec,
        const Coordinates<data_t> &coords) //!< the value of the advection terms
        const;                             //!< the value of the coordinates

    //! The function to compute NON-GAUGE dependent quantities
    //! (called in 'add_theory_rhs')
    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    void pre_compute_no_gauge(expressions<data_t> &E,
                              const vars_t<data_t> &vars,
                              const vars_t<Tensor<1, data_t>> &d1,
                              const diff2_vars_t<Tensor<2, data_t>> &d2,
                              const Tensor<2, data_t> &h_UU,
                              const Tensor<3, data_t> &chris_ULL) const;

  protected:
    coupling_and_potential_t my_coupling_and_potential;
    // const std::array<double, CH_SPACEDIM> &m_center;
};

#include "CubicHorndeski.impl.hpp"

#endif /* CUBICHORNDESKI_HPP_ */
