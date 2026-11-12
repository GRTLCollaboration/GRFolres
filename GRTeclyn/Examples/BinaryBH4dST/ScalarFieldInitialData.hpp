/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARFIELDINITIALDATA_HPP_
#define SCALARFIELDINITIALDATA_HPP_

#include "Coordinates.hpp"
#include "GRParmParse.hpp"
#include "StateVariables.hpp" // needs c_phi, c_Pi

#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>

#include <array>
#include <cmath>

// Sets a spherically-symmetric Gaussian bump for the scalar field phi with
// zero conjugate momentum Pi, matching the GRChombo BinaryBH4dST example.
class ScalarFieldInitialData
{
  public:
    struct params_t
    {
        amrex::Real amplitude{}; //!< amplitude of the Gaussian
        amrex::Real width{1.0};  //!< width of the Gaussian
        amrex::Real r0{30.0};    //!< radius of the shell of the Gaussian
        std::array<amrex::Real, AMREX_SPACEDIM> center{};

        static void check_params()
        {
            GRParmParse scalar_field_pp("scalar_field");
            amrex::Real scalar_amplitude{};
            scalar_field_pp.queryAdd("scalar_amplitude", scalar_amplitude);
            amrex::Real scalar_width{1.0};
            scalar_field_pp.queryAdd("scalar_width", scalar_width);
            if (scalar_width <= 0.0)
            {
                scalar_field_pp.error("scalar_width", "must be > 0.0");
            }
            amrex::Real scalar_r0{30.0};
            scalar_field_pp.queryAdd("scalar_r0", scalar_r0);

            GRParmParse geom_pp("geometry");
            std::array<amrex::Real, AMREX_SPACEDIM> scalar_center{};
            geom_pp.get("center", scalar_center);
            scalar_field_pp.queryAdd("scalar_center", scalar_center);
        }

        void fill_params()
        {
            GRParmParse scalar_field_pp("scalar_field");
            scalar_field_pp.get("scalar_amplitude", amplitude);
            scalar_field_pp.get("scalar_width", width);
            scalar_field_pp.get("scalar_r0", r0);

            GRParmParse geom_pp("geometry");
            geom_pp.get("center", center);
            scalar_field_pp.queryAdd("scalar_center", center);
        }
    };

    AMREX_FORCE_INLINE explicit ScalarFieldInitialData(amrex::Real a_dx)
        : m_dx(a_dx)
    {
        m_params.fill_params();
    }

    //! Function to compute the value of all the initial vars on the grid
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &state) const
    {
        const amrex::CellData<amrex::Real> &cell = state.cellData(ix, iy, iz);

        const Coordinates coords(amrex::IntVect(ix, iy, iz), m_dx,
                                 m_params.center);
        const amrex::Real rr = coords.get_radius();

        const amrex::Real arg = (rr - m_params.r0) / m_params.width;
        cell[c_phi] = m_params.amplitude * std::exp(-arg * arg);
        cell[c_Pi] = 0.0;
    }

  private:
    amrex::Real m_dx;
    params_t m_params{};
};

#endif /* SCALARFIELDINITIALDATA_HPP_ */
