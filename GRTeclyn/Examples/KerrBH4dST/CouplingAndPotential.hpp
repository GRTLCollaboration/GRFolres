/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef COUPLINGANDPOTENTIAL_HPP_
#define COUPLINGANDPOTENTIAL_HPP_

#include "GRParmParse.hpp"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>

// Shift-symmetric or exponential quadratic Einstein-scalar-Gauss-Bonnet 
// coupling plus an optional quadratic potential:
// f(phi) = lambda*phi or f(phi) = lambda/(2 * beta)*(1-exp(-beta*phi^2)),
// with a smooth interior excision based on chi.
// g2(phi) = g2 (constant coupling to the square of the kinetic term X).
// V(phi)  = 1/2 (scalar_mass * phi)^2.
class CouplingAndPotential
{
  public:
    struct params_t
    {
        amrex::Real lambda{};      // Gauss-Bonnet coupling constant
        amrex::Real cutoff{0.07};     // chi cutoff for the interior excision
        amrex::Real factor{100.0};    // sharpness of the excision transition
        amrex::Real scalar_mass{}; // mass of the scalar field
	int coupling_type{}; // Type of coupling function to use
                // 0: Shift symmetric Gauss-Bonnet coupling function f(phi) =
                // lambda * phi
                // 1: Exponential quadratic Gauss-Bonnet coupling function
                // f(phi) = lambda / (2 * beta) * (1 - exp(-beta * phi^2))
        amrex::Real beta{
            100.0}; // parameter for the exponential quadratic Gauss-Bonnet 
		    // coupling function. Only used if coupling_type == 1
        amrex::Real g2{}; // coupling to the square of the kinetic term

        static void check_params()
        {
            GRParmParse fdst_pp("four_deriv_scalar_tensor");
            amrex::Real lambda{};
	    amrex::Real cutoff{0.07};
	    amrex::Real factor{100.0};
	    amrex::Real scalar_mass{};
	    int coupling_type{};
	    amrex::Real beta{100.0};
	    amrex::Real g2{};

            fdst_pp.queryAdd("lambda", lambda);
            fdst_pp.queryAdd("cutoff", cutoff);
            fdst_pp.queryAdd("factor", factor);
            fdst_pp.queryAdd("scalar_mass", scalar_mass);
	    fdst_pp.queryAdd("coupling_type", coupling_type);
	    fdst_pp.queryAdd("beta", beta);
	    fdst_pp.queryAdd("g2", g2);

            // there's not yet a GRParmParse for kerr_spin
	    /*if (cutoff >=
                0.8 * 0.2666 *
                    std::sqrt(1.0 - m_params.kerr_spin * m_params.spin))
            {
                fdst.warning(
                    "cutoff", "Gauss-Bonnet cutoff may be too large.");
            }*/ 
	    
	    if (scalar_mass < 0.0)
            {
                fdst_pp.error("scalar_mass", "must be >= 0.0");
            }

            GRParmParse geometry_pp("geometry");
            amrex::Real coarsest_dx{};
            geometry_pp.get("coarsest_dx", coarsest_dx);

            GRParmParse evolution_pp("evolution");
            amrex::Real dt_multiplier{};
            evolution_pp.get("dt_multiplier", dt_multiplier);
            if (scalar_mass >= 0.2 / coarsest_dx / dt_multiplier)
            {
                fdst_pp.warning(
                    "scalar_mass",
                    "oscillations of the scalar field may not be resolved on "
                    "the coarsest level");
            }

	    if (coupling_type != 0 && coupling_type != 1)
            {
                fdst_pp.error(
                    "coupling_type", "only 0 (shift symmetric) or 1 "
                                     "(exponential) currently supported");
            }
        }

        void fill_params()
        {
            GRParmParse fdst_pp("four_deriv_scalar_tensor");
            fdst_pp.get("lambda", lambda);
            fdst_pp.get("cutoff", cutoff);
            fdst_pp.get("factor", factor);
            fdst_pp.get("scalar_mass", scalar_mass);
	    fdst_pp.get("coupling_type", coupling_type);
	    fdst_pp.get("beta", beta);
	    fdst_pp.get("g2", g2);
        }
    };

    CouplingAndPotential() { m_params.fill_params(); }

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE explicit CouplingAndPotential(params_t a_params)
        : m_params(a_params)
    {
    }

    // Set the EsGB coupling function and the scalar potential.
    // vars must provide vars.chi() and vars.phi().
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    template <class vars_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void compute_coupling_and_potential(
        amrex::Real &dfdphi, amrex::Real &d2fdphi2, amrex::Real &g2,
        amrex::Real &dg2dphi, amrex::Real &V_of_phi, amrex::Real &dVdphi,
        const vars_t &vars) const
    {
        // excision setting the coupling to 0 in the interior of the BH with a
	// smooth function
        const amrex::Real cutoff_factor =
            1.0 +
            std::exp(-m_params.factor * (vars.chi() - m_params.cutoff));

        // Shift-symmetric or exponential quadratic coupling
	// The first derivative of the GB coupling function
	const amrex::phi_squared = vars.phi() * vars.phi();
        dfdphi   = m_params.lambda / cutoff_factor * 
		(1 - coupling_type + coupling_type * vars.phi() *
	         std::exp(-m_params.beta * phi_squared));
	// The second derivative of the GB coupling function
        d2fdphi2 = coupling_type * m_params.lambda / cutoff_factor *
		(1.0 - 2.0 * m_params.beta * phi_squared) *
		std::exp(-m_params.beta * phi_squared);

        // coupling to the square of the kinetic term
        g2 = m_params.g2;
	// The first derivative of the g2 coupling
        dg2dphi = 0.0;

        // quadratic potential
        const amrex::Real mass_times_phi = m_params.scalar_mass * vars.phi();
        V_of_phi = 0.5 * mass_times_phi * mass_times_phi;
	// The first derivative of the potential
        dVdphi   = m_params.scalar_mass * m_params.scalar_mass * vars.phi();
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

  private:
    params_t m_params{};
};

#endif /* COUPLINGANDPOTENTIAL_HPP_ */
