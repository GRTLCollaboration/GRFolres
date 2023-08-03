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
#include "DefaultModifiedGauge.hpp"
#include "FourthOrderDerivatives.hpp"
#include "MovingPunctureGauge.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//!  Calculates RHS using CCZ4 in the modified gauge including matter terms,
//!  and matter variable evolution
/*!
     The class calculates the RHS evolution for all the variables. It inherits
   from the CCZ4RHS class, which it uses to do the non matter evolution of
   variables. It then adds in the terms depending on the a(x) and b(x) functions
   of the modified gauge (which are included through the class mod_gauge_t).
   Next, it included the additional matter terms to the CCZ4 evolution
   (those including the stress energy tensor), and calculates the evolution of
   the matter variables. It does not assume a specific form of matter but is
   templated over a matter class matter_t. Please see the class EsGB as
   an example of a matter_t. \sa CCZ4RHS(), EsGB()
*/

//! NOTE: In contrast with MatterCCZ4RHS, the quantity 8\pi G has been to 1 and
//! the coupling for the matter terms are taken into account into the matter
//! class itself

template <class matter_t, class gauge_t = MovingPunctureGauge,
          class deriv_t = FourthOrderDerivatives,
          class modified_gauge_t = DefaultModifiedGauge>
class ModifiedCCZ4RHS : public CCZ4RHS<gauge_t, deriv_t>
{
  public:
    // Use this alias for the same template instantiation as this class
    using CCZ4 = CCZ4RHS<gauge_t, deriv_t>;

    using params_t = CCZ4_params_t<typename gauge_t::params_t>;

    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    template <class data_t>
    using MatterDiff2Vars = typename matter_t::template Diff2Vars<data_t>;

    template <class data_t>
    using CCZ4Vars = typename CCZ4::template Vars<data_t>;

    template <class data_t>
    using CCZ4Diff2Vars = typename CCZ4::template Diff2Vars<data_t>;

    // Inherit the variable definitions from CCZ4RHS + matter_t
    template <class data_t>
    struct Vars : public CCZ4Vars<data_t>, public MatterVars<data_t>
    {
        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            CCZ4Vars<data_t>::enum_mapping(mapping_function);
            MatterVars<data_t>::enum_mapping(mapping_function);
        }
    };

    template <class data_t>
    struct Diff2Vars : public CCZ4Diff2Vars<data_t>,
                       public MatterDiff2Vars<data_t>
    {
        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            CCZ4Diff2Vars<data_t>::enum_mapping(mapping_function);
            MatterDiff2Vars<data_t>::enum_mapping(mapping_function);
        }
    };

    //!  Constructor of class ModifiedCCZ4RHS
    /*!
       Inputs are the grid spacing, plus the CCZ4 evolution parameters, the
       modified gauge functions and a matter object. It also takes the
       dissipation parameter sigma, and allows the formulation to be toggled
       between CCZ4 and BSSN. The default is CCZ4.
    */
    ModifiedCCZ4RHS(matter_t a_matter, params_t a_params,
                     modified_gauge_t a_modified_gauge, double a_dx, double a_sigma,
		     const std::array<double, CH_SPACEDIM> a_center, 
		     int a_formulation = CCZ4RHS<>::USE_CCZ4, 
		     double a_G_Newton = 1.0);

    //!  The compute member which calculates the RHS at each point in the box
    //!  \sa matter_rhs_equation()
    template <class data_t> void compute(Cell<data_t> current_cell) const;

    //! Same change as in CCZ4RHS.hpp
    // protected:

    //! The function which add in the modified gauge terms to the CCZ4 rhs
    template <class data_t>
    void add_a_b_rhs(
        Vars<data_t>
            &matter_rhs, //!< the RHS data for each variable at that point.
        const Vars<data_t> &vars, //!< the value of the variables at the point.
        const Vars<Tensor<1, data_t>>
            &d1, //!< the value of the first derivatives of the variables.
        const Diff2Vars<Tensor<2, data_t>>
            &d2, //!< the value of the second derivatives of the variables.
        const Vars<data_t> &advec, //!< the value of the advection terms.
        const Coordinates<data_t> &coords)
        const; //!< the value of the coordinates for computing a(x), b(x).

    //! The function which adds in the EM Tensor terms to the CCZ4 rhs \sa
    //! compute()
    template <class data_t>
    void add_emtensor_rhs(
        Vars<data_t>
            &matter_rhs, //!< the RHS data for each variable at that point.
        const Vars<data_t> &vars, //!< the value of the variables at the point.
        const Vars<Tensor<1, data_t>>
            &d1, //!< the value of the first derivatives of the variables.
        const Diff2Vars<Tensor<2, data_t>>
            &d2, //!< the value of the second derivatives of the variables.
        const Vars<data_t> &advec, //!< the value of the advection terms.
        const Coordinates<data_t> &coords)
        const; //!< the value of the coordinates for computing a(x), b(x).

    // Class members
    matter_t my_matter;       //!< The matter object, e.g. EsGB.
    modified_gauge_t my_modified_gauge; //!< The modified gauge object, i.e a(x) and b(x)
    const std::array<double, CH_SPACEDIM> m_center; //!< The center of the grid
    double m_G_Newton;
};

#include "ModifiedCCZ4RHS.impl.hpp"

#endif /* MODIFIEDCCZ4RHS_HPP_ */
