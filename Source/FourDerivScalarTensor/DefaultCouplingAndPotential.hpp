/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DEFAULTCOUPLINGANDPOTENTIAL_HPP_
#define DEFAULTCOUPLINGANDPOTENTIAL_HPP_

#include "Coordinates.hpp"
#include "Tensor.hpp"
#include "simd.hpp"

class DefaultCouplingAndPotential {
public:
  //! The constructor
  DefaultCouplingAndPotential() {}

  //! Set the 4dST coupling function here to zero
  template <class data_t, template <typename> class vars_t>
  void compute_coupling_and_potential(data_t &dfdphi, data_t &d2fdphi2,
                                      data_t &g2, data_t &dg2dphi,
                                      data_t &V_of_phi, data_t &dVdphi,
                                      const vars_t<data_t> &vars,
                                      const Coordinates<data_t> &coords) const {
    // The first derivative of the GB coupling function
    dfdphi = 0.;
    // The second derivative of the GB coupling function
    d2fdphi2 = 0.;
    // The coupling to the square of the kinetic term
    g2 = 0.;
    // The first derivative of the g2 coupling
    dg2dphi = 0.;
    // The potential of the scalar field
    V_of_phi = 0.;
    // The first derivative of the potential
    dVdphi = 0.;
  }
};

#endif /* DEFAULTCOUPLINGANDPOTENTIAL_HPP_ */
