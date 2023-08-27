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
#include "CouplingAndPotential.hpp"
#include "FourDerivScalarTensor.hpp"
#include "InitialScalarData.hpp"
#include "KerrBH.hpp"
#include "ModifiedPunctureGauge.hpp"

class SimulationParameters : public ModifiedGravitySimulationParametersBase<
                                 FourDerivScalarTensor<CouplingAndPotential>>
{
  public:
    SimulationParameters(GRParmParse &pp)
        : ModifiedGravitySimulationParametersBase(pp)
    {
        read_params(pp);
        check_params();
    }

    /// Read parameters from the parameter file
    void read_params(GRParmParse &pp)
    {
        // Coupling and potential
        pp.load("lambda_GB", coupling_and_potential_params.lambda_GB, 0.);
        pp.load("quadratic_factor",
                coupling_and_potential_params.quadratic_factor, 0.);
        pp.load("quartic_factor", coupling_and_potential_params.quartic_factor,
                0.);
        pp.load("cutoff_GB", coupling_and_potential_params.cutoff_GB, 0.15);
        pp.load("factor_GB", coupling_and_potential_params.factor_GB, 100.);
        pp.load("scalar_mass", coupling_and_potential_params.scalar_mass, 0.);

        // Initial data
        pp.load("G_Newton", G_Newton, 1.0);

        // Initial scalar field data
        pp.load("scalar_amplitude", initial_params.amplitude, 0.);
        pp.load("scalar_width", initial_params.width, 1.0);
        pp.load("scalar_r0", initial_params.r0, 0.);

        // Initial Kerr data
        pp.load("kerr_mass", kerr_params.mass);
        pp.load("kerr_spin", kerr_params.spin);
        pp.load("kerr_center", kerr_params.center, center);

#ifdef USE_AHFINDER
        pp.load("AH_initial_guess", AH_initial_guess, 0.5 * kerr_params.mass);
#endif
    }

    void check_params()
    {
        warn_parameter("kerr_mass", kerr_params.mass, kerr_params.mass >= 0.0,
                       "should be >= 0.0");
        check_parameter("kerr_spin", kerr_params.spin,
                        std::abs(kerr_params.spin) <= kerr_params.mass,
                        "must satisfy |a| <= M = " +
                            std::to_string(kerr_params.mass));
        FOR(idir)
        {
            std::string name = "kerr_center[" + std::to_string(idir) + "]";
            warn_parameter(
                name, kerr_params.center[idir],
                (kerr_params.center[idir] >= 0) &&
                    (kerr_params.center[idir] <= (ivN[idir] + 1) * coarsest_dx),
                "should be within the computational domain");
        }
    }

    double G_Newton;
    InitialScalarData::params_t initial_params;
    CouplingAndPotential::params_t coupling_and_potential_params;
    KerrBH::params_t kerr_params;

#ifdef USE_AHFINDER
    double AH_initial_guess;
#endif
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
