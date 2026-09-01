/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef DEFAULTCOUPLINGANDPOTENTIAL_HPP_
#define DEFAULTCOUPLINGANDPOTENTIAL_HPP_

#include "ScalarFieldVars.hpp"
#include "Tensor.hpp"
#include <AMReX_REAL.H>

class DefaultCouplingAndPotential
{
  public:
    //! The constructor
    DefaultCouplingAndPotential() = default;

    //! Set the potential function for the scalar field here to zero
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    compute_coupling_and_potential(amrex::Real &dfdphi, 
        amrex::Real &d2fdphi2, amrex::Real &V_of_phi, 
        amrex::Real &dVdphi, amrex::Real &g2, amrex::Real &dg2dphi,
        const ScalarFieldVars &vars)
    {
	dfdphi = 0.0;
	d2fdphi2 = 0.0;
	        
	// The potential value at phi
        V_of_phi = 0.0;

        // The potential gradient at phi
        dVdphi = 0.0;

	g2 = 0.0;
	dg2dphi = 0.0;
    }
};

#endif /* DEFAULTCOUPLINGANDPOTENTIAL_HPP_ */
