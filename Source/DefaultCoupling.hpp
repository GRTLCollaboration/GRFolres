/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DEFAULTCOUPLING_HPP_
#define DEFAULTCOUPLING_HPP_

#include "Tensor.hpp"
#include "simd.hpp"

class DefaultCoupling
{
  public:
    //! The constructor
    DefaultCoupling() {}

    //! Set the EsGB coupling function here to zero
    template <class data_t, template <typename> class vars_t>
    void compute_coupling(data_t &dfGB, data_t &df2GB,
                          const vars_t<data_t> &vars) const
    {
        // The first derivative of the coupling function
        dfGB = 0.;

        // The second derivative of the coupling function
        df2GB = 0.;
    }
};

#endif /* DEFAULTCOUPLING_HPP_ */
