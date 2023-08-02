/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MODIFIEDGRAVITYCONSTRAINTS_HPP_
#define MODIFIEDGRAVITYCONSTRAINTS_HPP_

#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "NewConstraints.hpp"
#include "Tensor.hpp"
#include "simd.hpp"
#include <array>

//!  Calculates the Hamiltonian and Momentum constraints with matter fields
/*!
     The class calculates the Hamiltonian and Momentum constraints at each point
   in a box. It inherits from the Constraints class which calculates the
   constraints without the matter terms. It adds in the matter terms for a given
   matter class matter_t, which must provide it with the Energy Momentum Tensor.
   For an example of a matter_t class see EsGB. \sa Constraints(),
   EsGB()
*/

template <class matter_t> class ModifiedGravityConstraints : public Constraints
{
  public:
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    // Inherit the variable definitions from CCZ4 + matter_t
    template <class data_t>
    struct BSSNMatterVars : public Constraints::MetricVars<data_t>,
                            public MatterVars<data_t>
    {
        /// Defines the mapping between members of Vars and Chombo grid
        /// variables (enum in User_Variables)
        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            Constraints::MetricVars<data_t>::enum_mapping(mapping_function);
            MatterVars<data_t>::enum_mapping(mapping_function);
        }
    };

    //! Constructor of class ModifiedGravityConstraints
    /*!
        Can specify the vars of the constraint vars instead of using the
        hardcoded ones.
    */
    ModifiedGravityConstraints(const matter_t a_matter, double dx, double G_Newton, int a_c_Ham,
                      const Interval &a_c_Moms, int a_c_Ham_abs_terms = -1,
                      const Interval &a_c_Moms_abs_terms = Interval());

    //! The compute member which calculates the constraints at each point in the
    //! box
    template <class data_t> void compute(Cell<data_t> current_cell) const;

  protected:
    matter_t my_matter; //!< The matter object, e.g. EsGB
    double m_G_Newton;
};

#include "ModifiedGravityConstraints.impl.hpp"

#endif /* MODIFIEDGRAVITYCONSTRAINTS_HPP_ */
