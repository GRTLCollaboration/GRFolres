/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef RHSDIAGNOSTICS_HPP_
#define RHSDIAGNOSTICS_HPP_

#include "BSSNVars.hpp"
#include "CCZ4Geometry.hpp"
#include "CCZ4RHS.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DiagnosticVariables.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "ModifiedCCZ4RHS.hpp"
#include "ModifiedPunctureGauge.hpp"
#include "Tensor.hpp"
#include "simd.hpp"
#include <array>

//! Calculates all the diagnostics in a specific theory, which need
//! the rhs of the variables, such as the weak coupling condition

template <class theory_t, class gauge_t = ModifiedPunctureGauge,
          class deriv_t = FourthOrderDerivatives>
class RHSDiagnostics : public ModifiedCCZ4RHS<theory_t, gauge_t, deriv_t>
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

    template <class data_t>
    using Vars = typename ModifiedCCZ4RHS<theory_t, gauge_t,
                                          deriv_t>::template Vars<data_t>;

    template <class data_t>
    using Diff2Vars =
        typename ModifiedCCZ4RHS<theory_t, gauge_t,
                                 deriv_t>::template Diff2Vars<data_t>;

    /*    // Inherit the variable definitions from CCZ4 + theory_t
        template <class data_t>
        struct BSSNTheoryVars : public BSSNVars::VarsNoGauge<data_t>,
                                public TheoryVars<data_t>
        {
            /// Defines the mapping between members of Vars and Chombo grid
            /// variables (enum in User_Variables)
            template <typename mapping_function_t>
            void enum_mapping(mapping_function_t mapping_function)
            {
                BSSNVars::VarsNoGauge<data_t>::enum_mapping(mapping_function);
                TheoryVars<data_t>::enum_mapping(mapping_function);
            }
        };*/

    //! Constructor of class RHSDiagnostics
    RHSDiagnostics(theory_t a_theory, modified_params_t a_params,
                   gauge_t a_gauge, double a_dx, double a_sigma,
                   const std::array<double, CH_SPACEDIM> a_center,
                   double a_G_Newton = 1.0);

    //! The compute member which calculates the RHS at each point in the
    //! box
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    theory_t my_theory; //!< The theory object, e.g. 4dST
    gauge_t my_gauge;   //!< The gauge object, which includes a(x) and b(x)
    const std::array<double, CH_SPACEDIM> m_center; //!< The center of the grid
    // const FourthOrderDerivatives m_deriv;
    double m_G_Newton;
    // modified_params_t m_params;
};

#include "RHSDiagnostics.impl.hpp"

#endif /* RHSDIAGNOSTICS_HPP_ */