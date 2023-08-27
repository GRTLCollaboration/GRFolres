/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "ModifiedGravitySimulationParametersBase.hpp"

// Problem specific includes:
#include "BinarySFTaggingCriterion.hpp"
#include "BoostedBH.hpp"
#include "CCZ4.hpp"
#include "CouplingAndPotential.hpp"
#include "CubicHorndeski.hpp"
#include "InitialScalarData.hpp"

class SimulationParameters : public ModifiedGravitySimulationParametersBase<
                                 CubicHorndeski<CouplingAndPotential>>
{
  public:
    SimulationParameters(GRParmParse &pp)
        : ModifiedGravitySimulationParametersBase(pp)
    {
        read_params(pp);
        check_params();
    }

    /// Read parameters
    void read_params(GRParmParse &pp)
    {
        // Do we want puncture tracking and constraint norm calculation?
        pp.load("track_punctures", track_punctures, false);
        pp.load("puncture_tracking_level", puncture_tracking_level, max_level);
        // Two above not used (should we add them in Level.cpp?)
        pp.load("calculate_constraint_norms", calculate_constraint_norms,
                false);

        // Coupling and potential
        pp.load("scalar_mass", coupling_and_potential_params.scalar_mass);
        pp.load("g3", coupling_and_potential_params.g3);
        pp.load("g2", coupling_and_potential_params.g2);

        // Initial data
        pp.load("G_Newton", G_Newton, 1.0);

        // Initial scalar field data
        initial_params.center =
            center; // already read in SimulationParametersBase
        pp.load("scalar_amplitude", initial_params.amplitude, 0.);
        pp.load("scalar_width", initial_params.width, 1.0);
        pp.load("scalar_r0", initial_params.r0, 0.);

        // Initial BH data
        pp.load("massA", bh1_params.mass);
        pp.load("momentumA", bh1_params.momentum);
        pp.load("massB", bh2_params.mass);
        pp.load("momentumB", bh2_params.momentum);

        // Get the centers of the BHs either explicitly or as
        // an offset (not both, or they will be offset from center
        // provided)
        std::array<double, CH_SPACEDIM> centerA, centerB;
        std::array<double, CH_SPACEDIM> offsetA, offsetB;
        pp.load("centerA", centerA, center);
        pp.load("centerB", centerB, center);
        pp.load("offsetA", offsetA, {0.0, 0.0, 0.0});
        pp.load("offsetB", offsetB, {0.0, 0.0, 0.0});
        FOR(idir)
        {
            bh1_params.center[idir] = centerA[idir] + offsetA[idir];
            bh2_params.center[idir] = centerB[idir] + offsetB[idir];
        }

        // Other Tiago's parameters (to be revised)
        pp.load("make_SF_zero_at_restart", make_SF_zero_at_restart, false);
        pp.load("excise_time", excise_time, 1.e8); // very big
        pp.load("calculate_constraint_norm_interval",
                calculate_constraint_norm_interval, 1);

        // Tagging parameters
        pp.load("regrid_threshold_chi", tag_params.threshold_chi);
        pp.load("regrid_threshold_phi", tag_params.threshold_phi);
        // tag_params.amplitudeSF =
        //    std::max(initial_params_1.amplitude, initial_params_2.amplitude);

        if (activate_extraction)
        {
            // tag_params.activate_extraction = activate_extraction;
            // tag_params.extraction_params = extraction_params;
        }

        tag_params.do_min_chi = true;
        int num_chi_contours;
        pp.load("num_chi_contours", num_chi_contours);
        pp.getarr("chi_contours", tag_params.chi_contours, 0, num_chi_contours);

        int num_tag_force_times;
        pp.load("num_tag_force_times", num_tag_force_times);
        pp.getarr("tag_force_times", tag_params.tag_force_times, 0,
                  num_tag_force_times);
        pp.getarr("tag_force_radii", tag_params.tag_force_radii, 0,
                  num_tag_force_times);
        pp.getarr("tag_force_levels", tag_params.tag_force_levels, 0,
                  num_tag_force_times);

        tag_params.centers = {bh1_params.center, bh2_params.center};

#ifdef USE_AHFINDER
        pp.load("AH_1_initial_guess", AH_1_initial_guess,
                0.5 * bh1_params.mass);
        pp.load("AH_2_initial_guess", AH_2_initial_guess,
                0.5 * bh2_params.mass);
        pp.load("AH_set_origins_to_punctures", AH_set_origins_to_punctures,
                false);
#endif
    }

    void check_params()
    {
        warn_parameter("massA", bh1_params.mass, bh1_params.mass >= 0,
                       "should be >= 0");
        warn_parameter("massB", bh2_params.mass, bh2_params.mass >= 0,
                       "should be >= 0");
        warn_array_parameter(
            "momentumA", bh1_params.momentum,
            std::sqrt(ArrayTools::norm2(bh1_params.momentum)) <
                0.3 * bh1_params.mass,
            "approximation used for boosted BH only valid for small boosts");
        warn_array_parameter(
            "momentumB", bh2_params.momentum,
            std::sqrt(ArrayTools::norm2(bh2_params.momentum)) <
                0.3 * bh1_params.mass,
            "approximation used for boosted BH only valid for small boosts");
        FOR(idir)
        {
            std::string nameA = "centerA[" + std::to_string(idir) + "]";
            std::string nameB = "centerB[" + std::to_string(idir) + "]";
            double center_A_dir = bh1_params.center[idir];
            double center_B_dir = bh2_params.center[idir];
            warn_parameter(nameA, center_A_dir,
                           (center_A_dir >= 0.0) &&
                               (center_A_dir <= (ivN[idir] + 1) * coarsest_dx),
                           "should be within the computational domain");
            warn_parameter(nameB, center_B_dir,
                           (center_B_dir >= 0.0) &&
                               (center_B_dir <= (ivN[idir] + 1) * coarsest_dx),
                           "should be within the computational domain");
        }
        check_parameter("puncture_tracking_level", puncture_tracking_level,
                        (puncture_tracking_level >= 0) &&
                            (puncture_tracking_level <= max_level),
                        "must be between 0 and max_level (inclusive)");
    }

    bool track_punctures, calculate_constraint_norms;
    int puncture_tracking_level;

    // Tagging parameters
    BinarySFTaggingCriterion::params tag_params;

    // Collection of parameters necessary for initial conditions
    double G_Newton;
    InitialScalarData::params_t initial_params;
    BoostedBH::params_t bh2_params;
    BoostedBH::params_t bh1_params;
    CouplingAndPotential::params_t coupling_and_potential_params;

    // Other Tiago's parameters
    bool make_SF_zero_at_restart;
    double excise_time;
    int calculate_constraint_norm_interval;

#ifdef USE_AHFINDER
    double AH_1_initial_guess;
    double AH_2_initial_guess;
    bool AH_set_origins_to_punctures;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
