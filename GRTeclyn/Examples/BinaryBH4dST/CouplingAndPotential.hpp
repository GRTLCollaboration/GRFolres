/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef COUPLINGANDPOTENTIAL_HPP_
#define COUPLINGANDPOTENTIAL_HPP_

#include "Coordinates.hpp"
#include "GRParmParse.hpp"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>

//! Shift-symmetric (linear) Einstein-scalar-Gauss-Bonnet coupling plus an
//! optional quadratic potential, as used in the GRChombo BinaryBH4dST example.
//! f(phi) = lambda_GB * phi, with a smooth interior excision based on chi.
//! g2(phi) = g2 (constant coupling to the square of the kinetic term X).
//! V(phi)  = 1/2 (scalar_mass * phi)^2.
class CouplingAndPotential
{
  public:
    struct params_t
    {
        amrex::Real lambda_GB{0.0};   //!< Gauss-Bonnet coupling
        amrex::Real g2{0.0};          //!< coupling to the square of X
        amrex::Real cutoff_GB{0.15};  //!< chi cutoff for the interior excision
        amrex::Real factor_GB{100.0}; //!< sharpness of the excision transition
        amrex::Real scalar_mass{0.0}; //!< mass in the potential

        static void check_params()
        {
            GRParmParse fdst_pp("four_deriv_scalar_tensor");

            amrex::Real lambda_GB{0.0};
            fdst_pp.queryAdd("lambda_GB", lambda_GB);

            amrex::Real g2{0.0};
            fdst_pp.queryAdd("g2", g2);

            amrex::Real cutoff_GB{0.15};
            fdst_pp.queryAdd("cutoff_GB", cutoff_GB);

            amrex::Real factor_GB{100.0};
            fdst_pp.queryAdd("factor_GB", factor_GB);

            amrex::Real scalar_mass{0.0};
            fdst_pp.queryAdd("scalar_mass", scalar_mass);
            if (scalar_mass < 0.0)
            {
                fdst_pp.error("scalar_mass", "must be >= 0.0");
            }
        }

        void fill_params()
        {
            GRParmParse fdst_pp("four_deriv_scalar_tensor");
            fdst_pp.get("lambda_GB", lambda_GB);
            fdst_pp.get("g2", g2);
            fdst_pp.get("cutoff_GB", cutoff_GB);
            fdst_pp.get("factor_GB", factor_GB);
            fdst_pp.get("scalar_mass", scalar_mass);
        }
    };

    CouplingAndPotential() { m_params.fill_params(); }

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE explicit CouplingAndPotential(params_t a_params)
        : m_params(a_params)
    {
    }

    //! Set the EsGB coupling function and the scalar potential.
    //! vars must provide vars.chi() and vars.phi().
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    template <class vars_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void compute_coupling_and_potential(
        amrex::Real &dfdphi, amrex::Real &d2fdphi2, amrex::Real &g2,
        amrex::Real &dg2dphi, amrex::Real &V_of_phi, amrex::Real &dVdphi,
        const vars_t &vars, const Coordinates & /*coords*/) const
    {
        // excision: switch the coupling off smoothly in the BH interior
        const amrex::Real cutoff_factor =
            1.0 +
            std::exp(-m_params.factor_GB * (vars.chi() - m_params.cutoff_GB));

        // Shift-symmetric coupling f(phi) = lambda_GB * phi
        dfdphi   = m_params.lambda_GB / cutoff_factor;
        d2fdphi2 = 0.0;

        // coupling to the square of the kinetic term
        g2      = m_params.g2;
        dg2dphi = 0.0;

        // quadratic potential
        const amrex::Real mass_times_phi = m_params.scalar_mass * vars.phi();
        V_of_phi = 0.5 * mass_times_phi * mass_times_phi;
        dVdphi   = m_params.scalar_mass * m_params.scalar_mass * vars.phi();
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

  private:
    params_t m_params{};
};

#endif /* COUPLINGANDPOTENTIAL_HPP_ */
