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
        readParams(pp);
    }

    /// Read parameters from the parameter file
    void readParams(GRParmParse &pp) { auto_read_params(pp); }

    void auto_read_params(GRParmParse &pp)
    {
        pp.load("G_Newton", G_Newton, 1.);

        // Initial Data
        pp.load("scalar_amplitude", initial_params.amplitude, 0.);
        pp.load("scalar_width", initial_params.width, 1.0);
        pp.load("scalar_r0", initial_params.r0, 0.);

        // Coupling and potential
        pp.load("scalar_mass", coupling_and_potential_params.scalar_mass);
        pp.load("g3", coupling_and_potential_params.g3);
        pp.load("g2", coupling_and_potential_params.g2);
        // gr_potential_params.scalar_mass = potential_params.scalar_mass;
        pout() << "Using g3=" << coupling_and_potential_params.g3
               << " && g2=" << coupling_and_potential_params.g2 << std::endl;

        FOR1(i)
        {
            bh1_params.center[i] += center[i];
            bh2_params.center[i] += center[i];
        }

        pp.load("activate_extraction", activate_extraction, 0);

        pp.load("make_SF_zero_at_restart", make_SF_zero_at_restart, false);

        pp.load("excise_time", excise_time, 1.e8); // very big

        pp.load("calculate_constraint_norms", calculate_constraint_norms);
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
    bool calculate_constraint_norms;
    int calculate_constraint_norm_interval;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
