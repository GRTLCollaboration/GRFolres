/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARFIELDINITIALDATA_HPP_
#define SCALARFIELDINITIALDATA_HPP_

#include "Coordinates.hpp"
#include "GRParmParse.hpp"
#include "StateVariables.hpp"

#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>

#include <array>
#include <cmath>

//! Class which sets the initial scalar field config
class ScalarFieldInitialData
{
  public:
    struct params_t
    {
        amrex::Real amplitude{0.0}; //!< Amplitude of the Gaussian
        std::array<amrex::Real, AMREX_SPACEDIM>
            center{};           //!< Centre of perturbation in the initial SF
        amrex::Real width{1.0}; //!< Width of the Gaussian
        amrex::Real r0{
            30.0}; //!< The radial coordinate of the position of the Gaussian

        static void check_params()
        {
            GRParmParse scalar_field_pp("scalar_field");
            amrex::Real scalar_amplitude{0.0};
            amrex::Real scalar_width{1.0};
            amrex::Real scalar_r0{30.0};
            scalar_field_pp.queryAdd("scalar_amplitude", scalar_amplitude);
            scalar_field_pp.queryAdd("scalar_width", scalar_width);
            scalar_field_pp.queryAdd("scalar_r0", scalar_r0);
            if (scalar_width <= 0.0)
            {
                scalar_field_pp.error("scalar_width", "must be > 0.0");
            }

            if (scalar_r0 < 0.0)
            {
                scalar_field_pp.error("scalar_r0", "must be >= 0.0");
            }
        }

        void fill_params()
        {
            GRParmParse geometry_pp("geometry");
            geometry_pp.get("center", center);

            GRParmParse scalar_field_pp("scalar_field");
            scalar_field_pp.get("scalar_amplitude", amplitude);
            scalar_field_pp.get("scalar_width", width);
            scalar_field_pp.get("scalar_r0", r0);
        }
    };

    AMREX_FORCE_INLINE explicit ScalarFieldInitialData(amrex::Real a_dx)
        : m_dx(a_dx)
    {
        m_params.fill_params();
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &state) const
    {
        const amrex::CellData<amrex::Real> &cell = state.cellData(ix, iy, iz);

        const Coordinates coords(amrex::IntVect(ix, iy, iz), m_dx,
                                 m_params.center);

        const amrex::Real rr = coords.get_radius();

        const amrex::Real arg = (rr - m_params.r0) / m_params.width;
        const amrex::Real arg2 = arg * arg;

        cell[c_phi] = m_params.amplitude * std::exp(-arg2);
        cell[c_Pi] = 0.0;
    }

  private:
    amrex::Real m_dx;
    params_t m_params{};
};

#endif /* SCALARFIELDINITIALDATA_HPP_ */