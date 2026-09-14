/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALAREXTRACTION_HPP_
#define SCALAREXTRACTION_HPP_

#include "SphericalExtraction.hpp"

/*!
   The class allows the user to extract data from the grid for the scalar
   components over spherical shells at specified radii. The values may then be
   written to an output file, or integrated across the surfaces.
*/
class ScalarExtraction : public SphericalExtraction<1>
{
  public:
    //! The constructor
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    ScalarExtraction(const spherical_extraction_params_t &a_params,
                     amrex::Real a_dt, amrex::Real a_time, bool a_first_step,
                     amrex::Real a_restart_time = 0.0)
        : SphericalExtraction<1>(a_params, a_dt, a_time, a_first_step,
                                 a_restart_time)
    {
        amrex::Vector<BCParity> parities = {BCParity::even};
        this->add_derived_vars({0}, parities, "phi");
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    amrex::Vector<BCParity> parities = {BCParity::even};

    //! The old constructor which assumes it is called in specific_post_timestep
    //! so the first time step is when m_time == m_dt
    ScalarExtraction(const spherical_extraction_params_t &a_params,
                     amrex::Real a_dt, amrex::Real a_time,
                     amrex::Real a_restart_time = 0.0)
        : ScalarExtraction(a_params, a_dt, a_time, (a_dt == a_time),
                           a_restart_time)
    {
    }

    //! Execute the query
    void execute_query(ParticleInterpolator<1> *a_interpolator)
    {
        // extract the values of the Weyl scalars on the spheres
        this->extract(a_interpolator);

        if (this->m_params.write_extraction)
        {
            this->write_extraction(this->m_params.extraction_file_prefix);
        }

        // now calculate and write the requested spherical harmonic modes
        std::vector<std::pair<std::vector<amrex::ParticleReal>,
                              std::vector<amrex::ParticleReal>>>
            mode_integrals(m_num_modes);

        // note that this is normalised by multiplying by radius
        // NOLINTBEGIN(bugprone-easily-swappable-parameters)
        auto normalised_scalar_complex =
            [](std::vector<amrex::ParticleReal> scalar_reim_parts,
               amrex::ParticleReal r, amrex::ParticleReal, amrex::ParticleReal)
        {
            // here the std::vector<amrex::ParticleReal> passed will just have
            // the real part of the scalar as its only component
            return std::pair(r * scalar_reim_parts[0], 0.0);
        };
        // NOLINTEND(bugprone-easily-swappable-parameters)

        // add the modes that will be integrated
        for (int imode = 0; imode < this->m_num_modes; ++imode)
        {
            const auto &mode = this->m_modes[imode];
            constexpr int spin_quantum_number = 0;
            this->add_mode_integrand(spin_quantum_number, mode.first,
                                     mode.second, normalised_scalar_complex,
                                     mode_integrals[imode]);
        }

        // do the integration over the surface
        this->integrate();

        // write the integrals
        for (int imode = 0; imode < this->m_num_modes; ++imode)
        {
            const auto &mode = this->m_modes[imode];
            std::string integrals_filename =
                this->m_params.integral_file_prefix +
                std::to_string(mode.first) + std::to_string(mode.second);
            std::vector<std::vector<amrex::ParticleReal>>
                integrals_for_writing = {
                    std::move(mode_integrals[imode].first),
                    std::move(mode_integrals[imode].second)};
            std::vector<std::string> labels = {"integral Re", "integral Im"};
            this->write_integrals(integrals_filename, integrals_for_writing,
                                  labels);
        }
    }
};

#endif /* SCALAREXTRACTION_HPP_ */
