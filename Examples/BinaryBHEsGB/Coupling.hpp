/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COUPLING_HPP_
#define COUPLING_HPP_

#include "simd.hpp"

class Coupling
{
  public:
    struct params_t
    {
        double lambda_GB; // Gauss-Bonnet coupling
        double cutoff_GB; // cutoff for switching off the Gauss-Bonnet terms
                          // inside the BH
        double factor_GB; // factor for the function smoothening the GB cutoff
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    Coupling(params_t a_params) : m_params(a_params) {}

    //! Set the EsGB coupling function fhere
    template <class data_t, template <typename> class vars_t>
    void compute_coupling(data_t &dfGB, data_t &df2GB,
                          const vars_t<data_t> &vars) const
    {
        // The first derivative of the coupling
        // It is set to 0 to the interior of the BH with a smooth function
        dfGB =
            m_params.lambda_GB /
            (1. + exp(-m_params.factor_GB * (vars.chi - m_params.cutoff_GB)));

        // The second derivative of the coupling
        df2GB = 0.;
    }
};

#endif /* COUPLING_HPP_ */
