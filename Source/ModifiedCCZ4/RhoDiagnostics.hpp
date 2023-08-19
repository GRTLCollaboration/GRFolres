/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef RHODIAGNOSTICS_HPP_
#define RHODIAGNOSTICS_HPP_

#include "BSSNVars.hpp"
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "simd.hpp"
#include <array>

//! Calculates all the rho components in a specific theory, which
//! are stored as diagnostics

template <class theory_t> class RhoDiagnostics
{
  public:
    template <class data_t>
    using TheoryVars = typename theory_t::template Vars<data_t>;

    // Inherit the variable definitions from CCZ4 + theory_t
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
    };

    //! Constructor of class RhoDiagnostics
    RhoDiagnostics(const theory_t a_theory, double dx,
                   const std::array<double, CH_SPACEDIM> a_center,
                   int a_c_rho_phi, int a_c_rho_g2, int a_c_rho_g3,
                   int a_c_rho_GB);

    //! The compute member which calculates the constraints at each point in the
    //! box
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    theory_t my_theory; //!< The theory object, e.g. 4dST
    const std::array<double, CH_SPACEDIM> m_center; //!< The center of the grid
    const FourthOrderDerivatives m_deriv;
    const int m_c_rho_phi;
    const int m_c_rho_g2;
    const int m_c_rho_g3;
    const int m_c_rho_GB;
};

#include "RhoDiagnostics.impl.hpp"

#endif /* RHODIAGNOSTICS_HPP_ */
