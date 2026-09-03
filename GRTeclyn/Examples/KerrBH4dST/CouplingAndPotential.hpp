/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef COUPLINGANDPOTENTIAL_HPP_
#define COUPLINGANDPOTENTIAL_HPP_

#include "GRParmParse.hpp"
#include "ScalarFieldVars.hpp"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>

class CouplingAndPotential
{
  public:
    struct params_t
    {
        amrex::Real lambda{0.0};   // Gauss-Bonnet coupling constant
        amrex::Real cutoff{0.07};  // Cut-off for switching off the Gauss-Bonnet
                                   // terms inside the BH
        amrex::Real factor{100.0}; // Factor inside the smoothing function for
                                   // the Gauss-Bonnet cut-off
        amrex::Real scalar_mass{1.0}; // Mass of the scalar field

        int coupling_type{
            1}; // Type of coupling function to use
                // 1: Shift symmetric Gauss-Bonnet coupling function f(phi) =
                // lambda * phi 2: Exponential Gauss-Bonnet coupling function
                // f(phi) = ( lambda / (2 * beta) ) * ( 1 - exp( - beta * phi^2
                // ) )

        amrex::Real beta{
            100.0}; // Parameter for the exponential Gauss-Bonnet coupling
                    // function. Only used if coupling_type == 2

        static void check_params()
        {
            GRParmParse scalar_field_pp("scalar_field");
            amrex::Real scalar_mass{1.0};
            scalar_field_pp.queryAdd("scalar_mass", scalar_mass);
            if (scalar_mass < 0.0)
            {
                scalar_field_pp.error("scalar_mass", "must be >= 0.0");
            }

            GRParmParse gauss_bonnet_pp("gauss_bonnet");
            int coupling_type{1};
            amrex::Real lambda{0.};
            amrex::Real cutoff{0.07};
            amrex::Real factor{100.0};
            amrex::Real beta{100.0};
            gauss_bonnet_pp.queryAdd("coupling_type", coupling_type);
            gauss_bonnet_pp.queryAdd("lambda", lambda);
            gauss_bonnet_pp.queryAdd("cutoff", cutoff);
            gauss_bonnet_pp.queryAdd("factor", factor);
            gauss_bonnet_pp.queryAdd("beta", beta);
            if (lambda < 0.0)
            {
                gauss_bonnet_pp.error("lambda", "must be >= 0.0");
            }
            if (cutoff < 0.0)
            {
                gauss_bonnet_pp.error("cutoff", "must be >= 0.0");
            }
            // TODO: ONCE KERR_SPIN IS ADDED, RE-ADD THIS
            // Impose that cutoff is at most 80% of chi_average at the horizon.
            // <chi>_H ~ 0.2666*sqrt(1 - j^2) - C.1 of
            // https://iopscience.iop.org/article/10.1088/1361-6382/ac6fa9
            // if (cutoff >= 0.8 * 0.2666 * std::sqrt(1.0 - m_params.kerr_spin *
            // m_params.kerr_spin))
            //{
            //    gauss_bonnet_pp.warning(
            //        "cutoff",
            //        "Gauss-Bonnet cutoff may be too large."
            //    );
            //}
            if (factor < 0.0)
            {
                gauss_bonnet_pp.error("factor", "must be >= 0.0");
            }
            if (coupling_type != 1 && coupling_type != 2)
            {
                gauss_bonnet_pp.error("coupling_type",
                                      "only 1 (shift symmetric) or 2 "
                                      "(exponential) currently supported");
            }
            if (beta <= 0.0)
            {
                gauss_bonnet_pp.error("beta", "must be > 0.0");
            }

            GRParmParse geometry_pp("geometry");
            amrex::Real coarsest_dx{};
            geometry_pp.get("coarsest_dx", coarsest_dx);

            GRParmParse evolution_pp("evolution");
            amrex::Real dt_multiplier{};
            evolution_pp.get("dt_multiplier", dt_multiplier);
            if (scalar_mass >= 0.2 / coarsest_dx / dt_multiplier)
            {
                scalar_field_pp.warning(
                    "scalar_mass",
                    "oscillations of the scalar field may not be resolved on "
                    "the coarsest level");
            }
        }

        void fill_params()
        {
            GRParmParse gauss_bonnet_pp("gauss_bonnet");
            gauss_bonnet_pp.get("coupling_type", coupling_type);
            gauss_bonnet_pp.get("lambda", lambda);
            gauss_bonnet_pp.get("cutoff", cutoff);
            gauss_bonnet_pp.get("factor", factor);
            gauss_bonnet_pp.get("beta", beta);

            GRParmParse scalar_field_pp("scalar_field");
            scalar_field_pp.get("scalar_mass", scalar_mass);
        }
    };

    CouplingAndPotential() { m_params.fill_params(); }

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE explicit CouplingAndPotential(params_t a_params)
        : m_params(a_params)
    {
    }

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute_coupling_and_potential(amrex::Real &dfdphi, amrex::Real &d2fdphi2,
                                   amrex::Real &V_of_phi, amrex::Real &dVdphi,
                                   amrex::Real &g2, amrex::Real &dg2dphi,
                                   const ScalarFieldVars &vars) const
    {

        const amrex::Real cutoff_factor =
            1. + std::exp(-m_params.factor * (vars.chi() - m_params.cutoff));

        if (m_params.coupling_type == 1)
        {
            // Shift symmetric Gauss-Bonnet coupling function f(phi) = lambda *
            // phi

            // Compute first and second derivative of coupling function.
            dfdphi = m_params.lambda / cutoff_factor;
            d2fdphi2 = 0.0;
        }

        else if (m_params.coupling_type == 2)
        {
            // Exponential Gauss-Bonnet coupling function f(phi) = ( lambda / (2
            // * beta) ) * ( 1 - exp( - beta * phi^2 ) )

            // Compute first and second derivative of coupling function.
            const amrex::Real phi_squared = vars.phi() * vars.phi();
            dfdphi = m_params.lambda * vars.phi() *
                     std::exp(-m_params.beta * phi_squared) / cutoff_factor;
            d2fdphi2 =
                m_params.lambda * std::exp(-m_params.beta * phi_squared) *
                (1.0 - 2.0 * m_params.beta * phi_squared) / cutoff_factor;
        }
        else
        {
            // We should never reach this block, but set to minimal coupling
            // values just in case.
            dfdphi = 0.0;
            d2fdphi2 = 0.0;
        }

        // Compute potential and its derivative.
        const amrex::Real mass_times_phi = m_params.scalar_mass * vars.phi();
        V_of_phi = 0.5 * mass_times_phi * mass_times_phi;
        dVdphi = m_params.scalar_mass * m_params.scalar_mass * vars.phi();

        // Compute coupling to the square of the kinetic term and its first
        // derivative.
        g2 = 0.0;
        dg2dphi = 0.0;
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

  private:
    params_t m_params{};
};

#endif /* COUPLINGANDPOTENTIAL_HPP_ */
