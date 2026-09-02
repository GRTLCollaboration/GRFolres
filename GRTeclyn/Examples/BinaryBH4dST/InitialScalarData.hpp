/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INITIALSCALARDATA_HPP_
#define INITIALSCALARDATA_HPP_

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
class InitialScalarData
{
  public:
    struct params_t
    {
        amrex::Real amplitude{0.0}; //!< amplitude of the Gaussian
        amrex::Real width{1.0};     //!< width of the Gaussian
        amrex::Real r0{0.0};        //!< radius of the shell of the Gaussian
        std::array<amrex::Real, AMREX_SPACEDIM> center{};

        static void check_params()
        {
            GRParmParse pp("initial_scalar_data");
            amrex::Real amplitude{0.0};
            pp.queryAdd("amplitude", amplitude);
            amrex::Real width{1.0};
            pp.queryAdd("width", width);
            if (width <= 0.0)
            {
                pp.error("width", "must be > 0.0");
            }
            amrex::Real r0{0.0};
            pp.queryAdd("r0", r0);

            GRParmParse geom_pp("geometry");
            std::array<amrex::Real, AMREX_SPACEDIM> center{};
            geom_pp.get("center", center);
            pp.queryAdd("center", center);
        }

        void fill_params()
        {
            GRParmParse pp("initial_scalar_data");
            pp.get("amplitude", amplitude);
            pp.get("width", width);
            pp.get("r0", r0);

            GRParmParse geom_pp("geometry");
            geom_pp.get("center", center);
            pp.queryAdd("center", center);
        }
    };

    AMREX_FORCE_INLINE explicit InitialScalarData(amrex::Real a_dx) : m_dx(a_dx)
    {
        m_params.fill_params();
    }

    //! Function to compute the value of all the initial vars on the grid
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &state) const
    {
        const amrex::CellData<amrex::Real> &cell =
            state.cellData(ix, iy, iz);

        const Coordinates coords(amrex::IntVect(ix, iy, iz), m_dx,
                                 m_params.center);
        const amrex::Real rr = coords.get_radius();

        const amrex::Real arg = (rr - m_params.r0) / m_params.width;
        cell[c_phi] = m_params.amplitude * std::exp(-arg * arg);
        cell[c_Pi]  = 0.0;
    }

  private:
    amrex::Real m_dx;
    params_t m_params{};
};

#endif /* INITIALSCALARDATA_HPP_ */
