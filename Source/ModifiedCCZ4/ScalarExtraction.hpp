/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SCALAREXTRACTION_HPP_
#define SCALAREXTRACTION_HPP_

#include "SphericalExtraction.hpp"
//! This class allows extraction of the values of a real scalar field
//! (normalised by radius) on spherical shells at specified radii, and
//! integration over those shells.
//! Optionally one can also extract and integrate the Noether flux too
class ScalarExtraction : public SphericalExtraction
{
  private:
    static constexpr int m_scalar = 0;

  public:
    //! The constructor
    ScalarExtraction(const spherical_extraction_params_t &a_params, double a_dt,
                     double a_time, bool a_first_step,
                     double a_restart_time = 0.0)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                              a_restart_time)
    {
        add_var(c_phi, VariableType::evolution);
    }

    //! The old constructor which assumes it is called in specificPostTimeStep
    //! so the first time step is when m_time == m_dt
    ScalarExtraction(const spherical_extraction_params_t &a_params, double a_dt,
                     double a_time, double a_restart_time = 0.0)
        : ScalarExtraction(a_params, a_dt, a_time, (a_dt == a_time),
                           a_restart_time)
    {
    }

    //! Execute the query
    void execute_query(AMRInterpolator<Lagrange<4>> *a_interpolator)
    {
        // extract the values of the scalar field on the spheres
        extract(a_interpolator);

        if (m_params.write_extraction)
            write_extraction(m_params.extraction_file_prefix);

        // now calculate and write the requested spherical harmonic modes
        std::vector<std::pair<std::vector<double>, std::vector<double>>>
            mode_integrals(m_num_modes);

        // note that this is normalised by multiplying by radius
        auto normalised_scalar =
            [](std::vector<double> scalar_reim_parts, double r, double, double)
        {
            // here the std::vector<double> passed will have the
            // values of scalar_Re = scalar and scalar_Im = 0.0 as its first two
            // components
            return std::make_pair(r * scalar_reim_parts[m_scalar], 0.0);
        };

        // add the modes that will be integrated
        for (int imode = 0; imode < m_num_modes; ++imode)
        {

            const auto &mode = m_modes[imode];
            // scalars have spin weight 0
            constexpr int es = 0;
            add_mode_integrand(es, mode.first, mode.second, normalised_scalar,
                               mode_integrals[imode]);
        }

        // do the integration
        integrate();

        // write the mode integrals
        for (int imode = 0; imode < m_num_modes; ++imode)
        {
            const auto &mode = m_modes[imode];
            std::string integrals_filename = m_params.integral_file_prefix +
                                             std::to_string(mode.first) +
                                             std::to_string(mode.second);
            std::vector<std::vector<double>> integrals_for_writing = {
                std::move(mode_integrals[imode].first),
                std::move(mode_integrals[imode].second)};
            std::vector<std::string> labels = {"integral Re", "integral Im"};
            write_integrals(integrals_filename, integrals_for_writing, labels);
        }
    }
};

#endif /* COMPLEXSCALAREXTRACTION_HPP_ */
