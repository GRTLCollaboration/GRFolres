/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MODGAUGE_HPP_
#define MODGAUGE_HPP_

#include "simd.hpp"

class ModGauge
{
  public:
    struct params_t
    {
        double a0; // constant value of a(x)
        double b0; // constant value of b(x)
        //!< This class has been prepared so that a and b can be
        //!< dependent of the coordinates
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    ModGauge(params_t a_params) : m_params(a_params) {}

    //! Set the modified gauge functions here
    template <class data_t, template <typename> class coords_t>
    void compute_mod_gauge(data_t &a_of_x, data_t &b_of_x,
                           const coords_t<data_t> &coords) const
    {
        // a(x) from \tilde{g}^{\mu\nu} = g^{\mu\nu} - a(x)n^{\mu}n^{\nu}
        a_of_x = m_params.a0;

        // b(x) from \hat{g}^{\mu\nu} = g^{\mu\nu} - b(x)n^{\mu}n^{\nu}
        b_of_x = m_params.b0;
    }
};

#endif /* MODGAUGE_HPP_ */
