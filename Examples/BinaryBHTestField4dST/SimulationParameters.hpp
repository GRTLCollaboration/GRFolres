/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "ArrayTools.hpp"
#include "BoostedBH.hpp"
#include "CouplingAndPotential.hpp"
#include "InitialScalarData.hpp"
#include "ModifiedPunctureGauge.hpp"
#include "SphericalExtraction.hpp"
#include "TestField4dST.hpp"

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
        pp.load("lambda_GB", coupling_and_potential_params.lambda_GB, 0.);
        pp.load("g2", coupling_and_potential_params.g2, 0.);
        pp.load("cutoff_GB", coupling_and_potential_params.cutoff_GB, 0.15);
        pp.load("factor_GB", coupling_and_potential_params.factor_GB, 100.);
        pp.load("scalar_mass", coupling_and_potential_params.scalar_mass, 0.);

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

        std::string extraction_path;
        if (pp.contains("extraction_subpath"))
        {
            pp.load("extraction_subpath", extraction_path);
            if (!extraction_path.empty() && extraction_path.back() != '/')
                extraction_path += "/";
            if (output_path != "./" && !output_path.empty())
                extraction_path = output_path + extraction_path;
        }
        else
            extraction_path = data_path;

        scalar_extraction_params.data_path = data_path;
        scalar_extraction_params.extraction_path = extraction_path;
        pp.load("scalar_integral_file_prefix",
                scalar_extraction_params.integral_file_prefix,
                std::string("scalar_mode_"));
        // by default just extract on the same spheres and with the same
        // parameters as for WeylExtraction
        pp.load("scalar_num_extraction_radii",
                scalar_extraction_params.num_extraction_radii,
                extraction_params.num_extraction_radii);
        pp.load("scalar_extraction_levels",
                scalar_extraction_params.extraction_levels,
                scalar_extraction_params.num_extraction_radii,
                extraction_params.extraction_levels);
        pp.load("scalar_extraction_radii",
                scalar_extraction_params.extraction_radii,
                scalar_extraction_params.num_extraction_radii,
                extraction_params.extraction_radii);
        pp.load("scalar_num_points_phi",
                scalar_extraction_params.num_points_phi,
                extraction_params.num_points_phi);
        pp.load("scalar_num_points_theta",
                scalar_extraction_params.num_points_theta,
                extraction_params.num_points_theta);
        pp.load("scalar_extraction_center", scalar_extraction_params.center,
                extraction_params.center);
        pp.load("scalar_write_extraction",
                scalar_extraction_params.write_extraction,
                extraction_params.write_extraction);

        // if scalar extraction is activated then the modes must be
        // specified
        pp.load("scalar_num_modes", scalar_extraction_params.num_modes);
        std::vector<int> scalar_extraction_modes_vect(
            2 * scalar_extraction_params.num_modes);
        pp.load("scalar_modes", scalar_extraction_modes_vect,
                2 * scalar_extraction_params.num_modes);
        scalar_extraction_params.modes.resize(
            scalar_extraction_params.num_modes);
        for (int imode; imode < scalar_extraction_params.num_modes; ++imode)
        {
            scalar_extraction_params.modes[imode] = {
                scalar_extraction_modes_vect[2 * imode],
                scalar_extraction_modes_vect[2 * imode + 1]};
        }

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
    spherical_extraction_params_t scalar_extraction_params;

    // Collection of parameters necessary for initial conditions
    double G_Newton;
    InitialScalarData::params_t initial_params;
    CouplingAndPotential::params_t coupling_and_potential_params;
    BoostedBH::params_t bh2_params;
    BoostedBH::params_t bh1_params;
    ModifiedCCZ4RHS<TestField4dST<CouplingAndPotential>, ModifiedPunctureGauge,
                    FourthOrderDerivatives>::modified_params_t
        modified_ccz4_params;

#ifdef USE_AHFINDER
    double AH_1_initial_guess;
    double AH_2_initial_guess;
    bool AH_set_origins_to_punctures;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP */
