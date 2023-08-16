/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "BinarySFTaggingCriterion.hpp"
#include "BoostedBH.hpp"
#include "CCZ4.hpp"
#include "CouplingAndPotential.hpp"
#include "CubicHorndeski.hpp"
#include "InitialScalarData.hpp"
#include "ModifiedPunctureGauge.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
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

        // Initial scalar field data
        initial_params.center =
            center; // already read in SimulationParametersBase
        pp.load("G_Newton", G_Newton, 1.0);
        pp.load("scalar_amplitude", initial_params.amplitude, 0.);
        pp.load("scalar_width", initial_params.width, 1.0);
        pp.load("scalar_r0", initial_params.r0, 0.);

        // Coupling and potential
        pp.load("scalar_mass", coupling_and_potential_params.scalar_mass);
        pp.load("g3", coupling_and_potential_params.g3);
        pp.load("g2", coupling_and_potential_params.g2);

        // Modified gauge
        pp.load("a0", modified_ccz4_params.a0, 0.);
        pp.load("b0", modified_ccz4_params.b0, 0.);
        pp.load("lapse_advec_coeff", modified_ccz4_params.lapse_advec_coeff,
                0.);
        pp.load("lapse_power", modified_ccz4_params.lapse_power, 1.);
        pp.load("lapse_coeff", modified_ccz4_params.lapse_coeff, 2.);
        pp.load("shift_Gamma_coeff", modified_ccz4_params.shift_Gamma_coeff,
                0.75);
        pp.load("shift_advec_coeff", modified_ccz4_params.shift_advec_coeff,
                0.);
        pp.load("eta", modified_ccz4_params.eta, 1.);
        modified_ccz4_params.kappa1 = ccz4_base_params.kappa1;
        modified_ccz4_params.kappa2 = ccz4_base_params.kappa2;
        modified_ccz4_params.kappa3 = ccz4_base_params.kappa3;
        modified_ccz4_params.covariantZ4 = ccz4_base_params.covariantZ4;

        // Initial data
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

        pp.load("activate_extraction", activate_extraction, 0);

        pp.load("make_SF_zero_at_restart", make_SF_zero_at_restart, false);

        pp.load("excise_time", excise_time, 1.e8); // very big

        pp.load("calculate_constraint_norm_interval",
                calculate_constraint_norm_interval, 1);

        // TAGGING PARAMETERS
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
    }

    void check_params()
    {
        check_parameter("a(x)", modified_ccz4_params.a0,
                        modified_ccz4_params.a0 > -1, "should be >-1");
        warn_parameter("a(x)", modified_ccz4_params.a0,
                       modified_ccz4_params.a0 != 0, "should be !=0");
        warn_parameter("b(x)", modified_ccz4_params.b0,
                       (modified_ccz4_params.b0 > 0) &&
                           (modified_ccz4_params.b0 != modified_ccz4_params.a0),
                       "should be >0 and !=a(x)");
        check_parameter("kappa1", modified_ccz4_params.kappa1,
                        modified_ccz4_params.kappa1 > 0, "should be > 0");
        check_parameter("kappa2", modified_ccz4_params.kappa2,
                        modified_ccz4_params.kappa2 >
                            -2. / (2. + modified_ccz4_params.b0),
                        "should be > -2/(2+b(x))");
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

    BinarySFTaggingCriterion::params tag_params;
    int activate_extraction;

    double G_Newton;
    InitialScalarData::params_t initial_params;
    BoostedBH::params_t bh2_params;
    BoostedBH::params_t bh1_params;
    CouplingAndPotential::params_t coupling_and_potential_params;
    ModifiedCCZ4RHS<CubicHorndeski<CouplingAndPotential>, ModifiedPunctureGauge,
                    FourthOrderDerivatives>::modified_params_t
        modified_ccz4_params;

    bool make_SF_zero_at_restart;
    double excise_time;
    int calculate_constraint_norm_interval;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
