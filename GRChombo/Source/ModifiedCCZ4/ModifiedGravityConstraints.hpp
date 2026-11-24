/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MODIFIEDGRAVITYCONSTRAINTS_HPP_
#define MODIFIEDGRAVITYCONSTRAINTS_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "NewConstraints.hpp"
#include "Tensor.hpp"
#include "simd.hpp"
#include <array>

//!  Calculates the Hamiltonian and Momentum constraints with the theory fields
/*!
     The class calculates the Hamiltonian and Momentum constraints at each point
   in a box. It inherits from the Constraints class which calculates the
   constraints without the theory terms. It adds in the theory terms for a given
   theory class theory_t, which must provide it with the Energy Momentum Tensor.
   For an example of a theory_t class see 4dST. \sa Constraints(),
   FourDerivScalarTensor()
*/

template <class theory_t> class ModifiedGravityConstraints : public Constraints
{
  public:
    template <class data_t>
    using TheoryVars = typename theory_t::template Vars<data_t>;

    // Inherit the variable definitions from CCZ4 + theory_t
    template <class data_t>
    struct BSSNTheoryVars : public Constraints::MetricVars<data_t>,
                            public TheoryVars<data_t>
    {
        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            Constraints::MetricVars<data_t>::enum_mapping(mapping_function);
            TheoryVars<data_t>::enum_mapping(mapping_function);
        }
    };

    //! Constructor of class ModifiedGravityConstraints
    /*!
        Can specify the vars of the constraint vars instead of using the
        hardcoded ones.
    */
    ModifiedGravityConstraints(const theory_t a_theory, double dx,
                               const std::array<double, CH_SPACEDIM> a_center,
                               double G_Newton, int a_c_Ham,
                               const Interval &a_c_Moms,
                               int a_c_Ham_abs_terms = -1,
                               const Interval &a_c_Moms_abs_terms = Interval());

    //! The compute member which calculates the constraints at each point in the
    //! box
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    theory_t my_theory; //!< The theory object, e.g. 4dST
    const std::array<double, CH_SPACEDIM> m_center; //!< The center of the grid
    double m_G_Newton;
};

#include "ModifiedGravityConstraints.impl.hpp"

#endif /* MODIFIEDGRAVITYCONSTRAINTS_HPP_ */
