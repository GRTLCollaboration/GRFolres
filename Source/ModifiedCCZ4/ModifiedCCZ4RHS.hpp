/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MODIFIEDCCZ4RHS_HPP_
#define MODIFIEDCCZ4RHS_HPP_

#include "CCZ4Geometry.hpp"
#include "CCZ4RHS.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "ModifiedPunctureGauge.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//!  Calculates RHS using CCZ4 in the modified gauge including theory terms,
//!  and theory variable evolution
/*!
     The class calculates the RHS evolution for all the variables. It inherits
   from the CCZ4RHS class, which it uses to do the evolution of the CCZ4
   variables. It then adds in the terms depending on the a(x) and b(x) functions
   of the modified gauge (which are included through the class gauge_t).
   Next, it included the additional theory terms to the CCZ4 evolution
   (those including the stress energy tensor), and calculates the evolution of
   the theory variables. It does not assume a specific form of theory but is
   templated over a theory class theory_t. Please see the class
   FourDerivScalarTensor as an example of a theory_t. \sa ModifiedCCZ4RHS(),
   FourDerivScalarTensor()
*/

template <class data_t> struct RhoAndSi
{
    Tensor<1, data_t> Si; //!< S_i = T_ia_n^a
    data_t rho;           //!< rho = T_ab n^a n^b
};

template <class data_t> struct SijTFAndS
{
    Tensor<2, data_t> Sij_TF; //!< S_ij_TF = (T_ab\gamma_i^a\gamma_j^b)^TF
    data_t S;                 //!< S = \gamma^ijT_ab\gamma_i^a\gamma_j^b
};

template <class data_t> struct ScalarVectorTensor
{
    data_t scalar;
    Tensor<1, data_t> vector;
    Tensor<2, data_t> tensor;
};

// Note: It may be needed to add other components for some theories
template <class data_t> struct AllRhos
{
    data_t phi;
    data_t g2;
    data_t g3;
    data_t GB;
};

template <class data_t> struct WeakCouplingConditions
{
    data_t g2;
    data_t g3;
    data_t GB;
};

template <class theory_t, class gauge_t = ModifiedPunctureGauge,
          class deriv_t = FourthOrderDerivatives>
class ModifiedCCZ4RHS : public CCZ4RHS<gauge_t, deriv_t>
{
  public:
    // Use this alias for the same template instantiation as this class
    using CCZ4 = CCZ4RHS<gauge_t, deriv_t>;

    using modified_params_t = CCZ4_params_t<typename gauge_t::params_t>;

    template <class data_t>
    using TheoryVars = typename theory_t::template Vars<data_t>;

    template <class data_t>
    using TheoryDiff2Vars = typename theory_t::template Diff2Vars<data_t>;

    template <class data_t>
    using CCZ4Vars = typename CCZ4::template Vars<data_t>;

    template <class data_t>
    using CCZ4Diff2Vars = typename CCZ4::template Diff2Vars<data_t>;

    // Inherit the variable definitions from CCZ4RHS + theory_t
    template <class data_t>
    struct Vars : public CCZ4Vars<data_t>, public TheoryVars<data_t>
    {
        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            CCZ4Vars<data_t>::enum_mapping(mapping_function);
            TheoryVars<data_t>::enum_mapping(mapping_function);
        }
    };

    template <class data_t>
    struct Diff2Vars : public CCZ4Diff2Vars<data_t>,
                       public TheoryDiff2Vars<data_t>
    {
        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            CCZ4Diff2Vars<data_t>::enum_mapping(mapping_function);
            TheoryDiff2Vars<data_t>::enum_mapping(mapping_function);
        }
    };

    //!  Constructor of class ModifiedCCZ4RHS
    /*!
       Inputs are the grid spacing, plus the CCZ4 evolution parameters, the
       modified gauge functions and a theory object. It also takes the
       dissipation parameter sigma. It allows the user to set the value of
       Newton's constant, which is set to one by default.
    */
    ModifiedCCZ4RHS(theory_t a_theory, modified_params_t a_params,
                    gauge_t a_gauge, double a_dx, double a_sigma,
                    const std::array<double, CH_SPACEDIM> a_center,
                    double a_G_Newton = 1.0);

    //!  The compute member which calculates the RHS at each point in the box
    //!  \sa theory_rhs_equation()
    template <class data_t> void compute(Cell<data_t> current_cell) const;

    //! Same change as in CCZ4RHS.hpp
    // protected:

    //! The function which add in the modified gauge terms to the CCZ4 rhs
    template <class data_t>
    void add_a_and_b_rhs(
        Vars<data_t>
            &theory_rhs, //!< the RHS data for each variable at that point.
        const Vars<data_t> &vars, //!< the value of the variables at the point.
        const Vars<Tensor<1, data_t>>
            &d1, //!< the value of the first derivatives of the variables.
        const Diff2Vars<Tensor<2, data_t>>
            &d2, //!< the value of the second derivatives of the variables.
        const Vars<data_t> &advec, //!< the value of the advection terms.
        const Coordinates<data_t> &coords)
        const; //!< the value of the coordinates.

    //! The function which adds in the EM Tensor terms to the CCZ4 rhs \sa
    //! compute()
    template <class data_t>
    void add_emtensor_rhs(
        Vars<data_t>
            &theory_rhs, //!< the RHS data for each variable at that point.
        const Vars<data_t> &vars, //!< the value of the variables at the point.
        const Vars<Tensor<1, data_t>>
            &d1, //!< the value of the first derivatives of the variables.
        const Diff2Vars<Tensor<2, data_t>>
            &d2, //!< the value of the second derivatives of the variables.
        const Vars<data_t> &advec, //!< the value of the advection terms.
        const Coordinates<data_t> &coords)
        const; //!< the value of the coordinates.

    //! Function to get full \kappa S_{ij}^{TF} (including LHS if needed) so as
    //! to be called in ModifiedGravityWeyl4 class
    template <class data_t>
    Tensor<2, data_t> get_full_kappa_times_Sij_TF(
        const Vars<data_t>
            &theory_vars, //!< the value of the variables at the point.
        const Vars<Tensor<1, data_t>>
            &d1, //!< the value of the first derivatives of the variables.
        const Diff2Vars<Tensor<2, data_t>>
            &d2, //!< the value of the second derivatives of the variables.
        const Vars<data_t> &advec, //!< the value of the advection terms.
        const Coordinates<data_t> &coords)
        const; //!< the value of the coordinates.

    // Class members
    theory_t my_theory; //!< The theory object, e.g. 4dST.
    gauge_t my_gauge;   //!< The gauge object, which includes a(x) and b(x)
    const std::array<double, CH_SPACEDIM> m_center; //!< The center of the grid
    double m_G_Newton;
};

#include "ModifiedCCZ4RHS.impl.hpp"

#endif /* MODIFIEDCCZ4RHS_HPP_ */
