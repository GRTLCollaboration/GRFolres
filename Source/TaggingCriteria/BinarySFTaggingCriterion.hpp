/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef BINARYSFTAGGINGCRITERION_HPP_
#define BINARYSFTAGGINGCRITERION_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "CubicHorndeski.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "SimulationParametersBase.hpp"
#include "Tensor.hpp"

class BinarySFTaggingCriterion
{
  public:
    struct params
    {
        double threshold_chi;
        double threshold_phi;
        double amplitudeSF;

        bool activate_extraction;
        extraction_params_t extraction_params;

        bool do_min_chi;
        Vector<double> chi_contours;

        Vector<double> tag_force_times, tag_force_radii;
        Vector<int> tag_force_levels;
        Vector<std::array<double, GR_SPACEDIM>> centers;
    };

  protected:
    const double m_dx;
    const FourthOrderDerivatives m_deriv;
    const int m_level, m_max_level;
    const double m_time;
    const params m_params;

    template <class data_t>
    using MatterVars = typename CubicHorndeski<>::template Vars<data_t>;

    /// Vars object for chi
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

  public:
    BinarySFTaggingCriterion(const double dx, const int a_level,
                             const int a_max_level, double a_time,
                             const params &a_params)
        : m_dx(dx), m_deriv(dx), m_level(a_level), m_max_level(a_max_level),
          m_time(a_time), m_params(a_params){};

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        const auto d2 = m_deriv.template diff2<MatterVars>(current_cell);
        const auto d2chi = m_deriv.template diff2<Vars>(current_cell);

        data_t mod_d2_chi = 0.;
        data_t mod_d2_phi = 0.;
        data_t mod_d2_Pi = 0.;

        FOR2(idir, jdir)
        {
            mod_d2_chi += d2chi.chi[idir][jdir] * d2chi.chi[idir][jdir];
            mod_d2_phi += d2.phi[idir][jdir] * d2.phi[idir][jdir];
            mod_d2_Pi += d2.Pi[idir][jdir] * d2.Pi[idir][jdir];
        }

        data_t criterion_chi = m_dx / m_params.threshold_chi * sqrt(mod_d2_chi);
        data_t criterion_phi = m_dx / m_params.threshold_phi *
                               (sqrt(mod_d2_phi) + sqrt(mod_d2_Pi)) /
                               m_params.amplitudeSF;

        data_t criterion = simd_max(criterion_chi, criterion_phi);

        // if extracting weyl data at a given radius, enforce a given resolution
        // there
        if (m_params.activate_extraction)
        {
            for (int iradius = 0;
                 iradius < m_params.extraction_params.num_extraction_radii;
                 ++iradius)
            {
                // regrid if within extraction level and not at required
                // refinement
                if (m_level <
                    m_params.extraction_params.extraction_levels[iradius])
                {
                    const Coordinates<data_t> coords(
                        current_cell, m_dx,
                        m_params.extraction_params.extraction_center);
                    const data_t r = coords.get_radius();
                    // add a 20% buffer to extraction zone so not too near to
                    // boundary
                    auto regrid = simd_compare_lt(
                        r, 1.2 * m_params.extraction_params
                                     .extraction_radii[iradius]);
                    criterion = simd_conditional(regrid, 100.0, criterion);
                }
            }
        }

        // ensure that the horizons of the punctures are covered
        // by the max level - for this we need
        // only check the chi contours on the top 4 levels
        // which regrid (ie, max_level - 1 to max_level - 4)
        // (just the top level would be ok, but doing more ensures
        // the top levels are well spaced)
        if ((m_level > (m_max_level - 1 - (int)m_params.chi_contours.size())) &&
            m_params.do_min_chi)
        {
            // decide whether to tag if chi < 0.35 (ensures AH, around 0.2 to
            // 0.3, is well within)
            auto chi = current_cell.load_vars(c_chi);
            // static data_t min_chi[] = {0.35, 0.45, 0.55, 0.65};
            // to ensure a gap between levels and avoid regridding issues
            auto regrid = simd_compare_lt(
                chi, m_params.chi_contours[(m_max_level - 1) - m_level]);
            criterion = simd_conditional(regrid, 100.0, criterion);
        }

        if (m_params.tag_force_times.size() > 0)
        {
            CH_assert(m_params.tag_force_times.size() ==
                      m_params.tag_force_levels.size());
            int what_time_index = 0;
            const double eps = 1.e-7;
            while (what_time_index < m_params.tag_force_times.size() &&
                   m_time > m_params.tag_force_times[what_time_index] - eps)
                ++what_time_index;
            if (what_time_index < m_params.tag_force_times.size())
            {
                if (m_level >= m_params.tag_force_levels[what_time_index])
                    criterion = 0.;
                else
                {
                    for (int i = 0; i < m_params.centers.size(); ++i)
                    {
                        const Coordinates<data_t> coords(current_cell, m_dx,
                                                         m_params.centers[i]);
                        const data_t r = coords.get_radius();
                        auto inside = simd_compare_lt(
                            r, m_params.tag_force_radii[what_time_index]);
                        criterion = simd_conditional(inside, 100., criterion);
                    }
                }
            }
        }

        // Write back into the flattened Chombo box
        current_cell.store_vars(criterion, 0);
    }
};

#endif /* BINARYSFTAGGINGCRITERION_HPP_ */
