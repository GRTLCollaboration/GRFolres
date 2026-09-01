/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COUPLINGANDPOTENTIAL_HPP_
#define COUPLINGANDPOTENTIAL_HPP_

#include "simd.hpp"

class CouplingAndPotential
{
  public:
    struct params_t
    {
        double dfdphi;
        double d2fdphi2;
        double g2;
        double dg2dphi;
        double V_of_phi;
        double dVdphi;
    };

  private:
    params_t m_params;

  public:
    //! The constructor
    CouplingAndPotential(params_t a_params) : m_params(a_params) {}

    //! Set the EsGB coupling function fhere
    template <class data_t, template <typename> class vars_t>
    void compute_coupling_and_potential(data_t &dfdphi, data_t &d2fdphi2,
                                        data_t &g2, data_t &dg2dphi,
                                        data_t &V_of_phi, data_t &dVdphi,
                                        const vars_t<data_t> &vars,
                                        const Coordinates<data_t> &coords) const
    {
        // The first derivative of the GB coupling function
        dfdphi = m_params.dfdphi;
        // The second derivative of the GB coupling function
        d2fdphi2 = m_params.d2fdphi2;
        // The coupling to the square of the kinetic term
        g2 = m_params.g2;
        // The first derivative of the g2 coupling
        dg2dphi = m_params.dg2dphi;
        // The potential of the scalar field
        V_of_phi = m_params.V_of_phi;
        // The first derivative of the potential
        dVdphi = m_params.dVdphi;
    }
};

#endif /* COUPLINGANDPOTENTIAL_HPP_ */
