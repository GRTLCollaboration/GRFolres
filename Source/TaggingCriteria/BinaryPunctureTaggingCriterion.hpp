/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BINARYPUNCTURETAGGINGCRITERION_HPP_
#define BINARYPUNCTURETAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SphericalExtraction.hpp"
#include "Tensor.hpp"
#include "VarsTools.hpp"

//! This class tags cells based on three criteria - the
//! value of the second derivs, the extraction regions
//! and the puncture horizons (which must be covered to
//! a given level
template <class derivative_t = FourthOrderDerivatives>
class BinaryPunctureTaggingCriterion
{
  protected:
    const double m_dx;
    const int m_level;
    const std::array<int, 2> m_horizon_max_levels;
    const std::array<int, 2> m_puncture_max_levels;
    const bool m_activate_extraction;
    const bool m_track_punctures;
    // const std::array<double, 2> m_puncture_masses;
    const std::vector<double> m_puncture_masses;
    // const SphericalExtraction::params_t m_params;
    const spherical_extraction_params_t m_params;
    const derivative_t m_deriv;
    const std::vector<std::array<double, CH_SPACEDIM>> m_puncture_coords;
    const std::array<double, 2> m_buffers;
    const double m_puncture_min_separation; // the minimum coordinate separation
                                            // to tag around the punctures
                                            // separately
    double m_max_buffer;

  public:
    template <class data_t> struct Vars
    {
        data_t chi; //!< Conformal factor

        template <typename mapping_function_t>
        void enum_mapping(mapping_function_t mapping_function)
        {
            using namespace VarsTools; // define_enum_mapping is part of
                                       // VarsTools
            define_enum_mapping(mapping_function, c_chi, chi);
        }
    };

    // The constructor
    BinaryPunctureTaggingCriterion(
        const double dx, const int a_level,
        const std::array<int, 2> a_horizon_max_levels,
        const std::array<int, 2> a_puncture_max_levels,
        // const SphericalExtraction::params_t a_params,
        const spherical_extraction_params_t a_params,
        const std::vector<std::array<double, CH_SPACEDIM>> a_puncture_coords,
        const bool activate_extraction = false,
        const bool track_punctures = false,
        // const std::array<double, 2> a_puncture_masses = {1.0, 1.0},
	    const std::vector<double> a_puncture_masses = {1.0, 1.0},
        const std::array<double, 2> a_buffers = {0.5, 0.5},
        const double a_puncture_min_separation = 1e-3)
        : m_dx(dx), m_deriv(dx), m_params(a_params), m_level(a_level),
          m_horizon_max_levels(a_horizon_max_levels),
          m_puncture_max_levels(a_puncture_max_levels),
          m_track_punctures(track_punctures),
          m_activate_extraction(activate_extraction),
          m_puncture_masses(a_puncture_masses),
          m_puncture_coords(a_puncture_coords), m_buffers(a_buffers),
          m_puncture_min_separation(a_puncture_min_separation)
    {
        // check that the number of punctures is consistent
        CH_assert(m_puncture_coords.size() == 2);
        m_max_buffer = max(m_buffers[0], m_buffers[1]);
    };

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // first test the gradients for regions of high curvature
        const auto d2 = m_deriv.template diff2<Vars>(current_cell);
        data_t mod_d2_chi = 0;
        FOR2(idir, jdir)
        {
            mod_d2_chi += d2.chi[idir][jdir] * d2.chi[idir][jdir];
        }
        data_t criterion = m_dx * sqrt(mod_d2_chi);

        // if extracting weyl data at a given radius, enforce a given resolution
        // there
        if (m_activate_extraction)
        {
            for (int iradius = 0; iradius < m_params.num_extraction_radii;
                 ++iradius)
            {
                // regrid if within extraction level and not at required
                // refinement
                if (m_level < m_params.extraction_levels[iradius])
                {
                    const Coordinates<data_t> coords(current_cell, m_dx,
                                                     m_params.center);
                    const data_t r = coords.get_radius();
                    // add a 20% buffer to extraction zone so not too near to
                    // boundary
                    auto regrid = simd_compare_lt(
                        r, 1.2 * m_params.extraction_radii[iradius]);
                    criterion = simd_conditional(regrid, 100.0, criterion);
                }
            }
        }

        // next use the punctures to decide where to tag
        if (m_track_punctures)
        {
            double puncture_separation2 = 0.0;
            FOR1(idir)
            {
                double displacement =
                    m_puncture_coords[0][idir] - m_puncture_coords[1][idir];
                puncture_separation2 += displacement * displacement;
            }
            double puncture_separation = sqrt(puncture_separation2);

            const int merger_horizon_max_level =
                min(m_horizon_max_levels[0], m_horizon_max_levels[1]);

            const int merger_puncture_max_level =
                min(m_puncture_max_levels[0], m_puncture_max_levels[1]);

            if (puncture_separation > m_puncture_min_separation)
            {
                // punctures still far enough apart so tag around each one
                // separately
                for (int ipuncture = 0; ipuncture < 2; ++ipuncture)
                {
                    const Coordinates<data_t> coords(
                        current_cell, m_dx, m_puncture_coords[ipuncture]);
                    const data_t max_abs_xy =
                        simd_max(abs(coords.x), abs(coords.y));
                    const data_t max_abs_xyz =
                        simd_max(max_abs_xy, abs(coords.z));

                    if (m_level < m_horizon_max_levels[ipuncture])
                    {
                        // we want the 2nd and 3rd levels above
                        // puncture_max_level to be twice the size of the next
                        // finest level
                        const double factor =
                            pow(2.0, min(m_horizon_max_levels[ipuncture] -
                                             m_level - 1,
                                         2));

                        auto regrid = simd_compare_lt(
                            max_abs_xyz,
                            factor * (m_puncture_masses[ipuncture] +
                                      m_buffers[ipuncture]));

                        criterion = simd_conditional(regrid, 100.0, criterion);
                    }
                    else if (m_level < m_puncture_max_levels[ipuncture])
                    {
                        // make the first level below horizon_max_level
                        // 1/4 the length across and levels below that
                        // 1/2 as small each time
                        const double factor = pow(
                            2.0, m_horizon_max_levels[ipuncture] - m_level - 2);

                        auto regrid = simd_compare_lt(
                            max_abs_xyz, factor * m_puncture_masses[ipuncture]);

                        criterion = simd_conditional(regrid, 100.0, criterion);
                    }
                    else
                    {
                        // remove any finer levels for BHs with
                        // puncture_max_level < max_level
                        auto dont_regrid = simd_compare_lt(
                            max_abs_xyz, m_puncture_masses[ipuncture] +
                                             m_buffers[ipuncture]);
                        criterion =
                            simd_conditional(dont_regrid, 0.0, criterion);
                    }
                }
            }
            else
            {
                if (m_level >= merger_horizon_max_level &&
                    m_level >= merger_puncture_max_level)
                {
                    // drop any finer levels after merger
                    // tagging for merger BH handled below
                    criterion = 0.0;
                }
            }

            double sum_masses = m_puncture_masses[0] + m_puncture_masses[1];
            // if punctures are close enough together tag cells at the
            // center of mass for the merger BH
            if (puncture_separation < sum_masses + m_max_buffer)
            {
                std::array<double, CH_SPACEDIM> center_of_mass;
                FOR1(idir)
                {
                    center_of_mass[idir] =
                        (m_puncture_masses[0] * m_puncture_coords[0][idir] +
                         m_puncture_masses[1] * m_puncture_coords[1][idir]) /
                        sum_masses;
                }
                Coordinates<data_t> coords(current_cell, m_dx, center_of_mass);
                const data_t max_abs_xy =
                    simd_max(abs(coords.x), abs(coords.y));
                const data_t max_abs_xyz = simd_max(max_abs_xy, abs(coords.z));

                if (m_level < merger_horizon_max_level)
                {
                    const double factor = pow(
                        2.0, min(merger_horizon_max_level - m_level - 1, 2));
                    auto regrid2 = simd_compare_lt(
                        max_abs_xyz, factor * (sum_masses + m_max_buffer));
                    criterion = simd_conditional(regrid2, 100.0, criterion);
                }
                else if (m_level < merger_puncture_max_level)
                {
                    // make the first level below horizon_max_level
                    // 1/4 the length across and levels below that
                    // 1/2 as small each time
                    const double factor =
                        pow(2.0, merger_horizon_max_level - m_level - 2);

                    auto regrid2 =
                        simd_compare_lt(max_abs_xyz, factor * sum_masses);

                    criterion = simd_conditional(regrid2, 100.0, criterion);
                }
            }
        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* BINARYPUNCTURETAGGINGCRITERION_HPP_ */
