/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MODIFIEDGRAVITYWEYL4_HPP_
#define MODIFIEDGRAVITYWEYL4_HPP_

#include "CCZ4RHS.hpp"
#include "Coordinates.hpp"
#include "Weyl4.hpp"

//!  Calculates the Weyl4 scalar for spacetimes in a specific theory
/*!
   This class calculates the Weyl4 scalar real and im parts. It inherits from
   the Weyl4 class and ModifiedCCZ4RHS classes and adds in the additional
   corresponding terms
*/
template <class theory_t, class gauge_t, class deriv_t>
class ModifiedGravityWeyl4 : public Weyl4,
                             public ModifiedCCZ4RHS<theory_t, gauge_t, deriv_t>
{
  public:
    template <class data_t>
    using Vars = typename ModifiedCCZ4RHS<theory_t, gauge_t,
                                          deriv_t>::template Vars<data_t>;
    template <class data_t>
    using Diff2Vars =
        typename ModifiedCCZ4RHS<theory_t, gauge_t,
                                 deriv_t>::template Diff2Vars<data_t>;
    using modified_params_t = CCZ4_params_t<typename gauge_t::params_t>;

    //! Constructor
    ModifiedGravityWeyl4(theory_t a_theory, modified_params_t a_params,
                         gauge_t a_gauge,
                         const std::array<double, CH_SPACEDIM> a_center,
                         const double a_dx, double a_sigma,
                         const int a_formulation = CCZ4RHS<>::USE_CCZ4)
        : Weyl4(a_center, a_dx, a_formulation),
          ModifiedCCZ4RHS<theory_t, gauge_t, deriv_t>(
              a_theory, a_params, a_gauge, a_dx, a_sigma, a_center),
          my_theory(a_theory)
    {
    }

    //! The compute member which calculates the wave quantities at each point on
    //! the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    theory_t my_theory; //!< The theory object, e.g. 4dST

    //! Add modified gravity terms to electric and magnetic parts
    template <class data_t>
    void add_theory_EB(
        EBFields_t<data_t>
            &eb_fields, //!< the value of the E and B Fields at the point.
        const Vars<data_t> &vars, //!< the value of the variables at the point.
        const Vars<Tensor<1, data_t>>
            &d1, //!< the value of the first derivatives of the variables.
        const Diff2Vars<Tensor<2, data_t>>
            &d2, //!< the value of the second derivatives of the variables.
        const Vars<data_t> &advec, //!< the value of the advection terms.
        const Coordinates<data_t> &coords)
        const; //!< the value of the coordinates.
};

#include "ModifiedGravityWeyl4.impl.hpp"

#endif /* MODIFIEDGRAVITYWEYL4_HPP_ */
