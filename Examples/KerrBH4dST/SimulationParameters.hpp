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
#include "CouplingAndPotential.hpp"
#include "FourDerivScalarTensor.hpp"
#include "InitialScalarData.hpp"
#include "KerrBH.hpp"
#include "ModifiedPunctureGauge.hpp"

class SimulationParameters : public SimulationParametersBase {
public:
  SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp) {
    read_params(pp);
    check_params();
  }

  /// Read parameters from the parameter file
  void read_params(GRParmParse &pp) {
    pp.load("G_Newton", G_Newton, 1.0);
    // Initial Data
    pp.load("scalar_amplitude", initial_params.amplitude, 0.);
    pp.load("scalar_width", initial_params.width, 1.0);
    // Coupling and potential
    pp.load("lambda_GB", coupling_and_potential_params.lambda_GB, 0.);
    pp.load("cutoff_GB", coupling_and_potential_params.cutoff_GB, 0.15);
    pp.load("factor_GB", coupling_and_potential_params.factor_GB, 100.);
    pp.load("scalar_mass", coupling_and_potential_params.scalar_mass, 0.);
    // Modified gauge
    pp.load("a0", modified_ccz4_params.a0, 0.);
    pp.load("b0", modified_ccz4_params.b0, 0.);
    pp.load("lapse_advec_coeff", modified_ccz4_params.lapse_advec_coeff, 0.);
    pp.load("lapse_power", modified_ccz4_params.lapse_power, 1.);
    pp.load("lapse_coeff", modified_ccz4_params.lapse_coeff, 2.);
    pp.load("shift_Gamma_coeff", modified_ccz4_params.shift_Gamma_coeff, 0.75);
    pp.load("shift_advec_coeff", modified_ccz4_params.shift_advec_coeff, 0.);
    pp.load("eta", modified_ccz4_params.eta, 1.);
    modified_ccz4_params.kappa1 = ccz4_base_params.kappa1;
    modified_ccz4_params.kappa2 = ccz4_base_params.kappa2;
    modified_ccz4_params.kappa3 = ccz4_base_params.kappa3;
    modified_ccz4_params.covariantZ4 = ccz4_base_params.covariantZ4;

    // Initial Kerr data
    pp.load("kerr_mass", kerr_params.mass);
    pp.load("kerr_spin", kerr_params.spin);
    pp.load("kerr_center", kerr_params.center, center);
  }

  void check_params() {
    warn_parameter("kerr_mass", kerr_params.mass, kerr_params.mass >= 0.0,
                   "should be >= 0.0");
    check_parameter("kerr_spin", kerr_params.spin,
                    std::abs(kerr_params.spin) <= kerr_params.mass,
                    "must satisfy |a| <= M = " +
                        std::to_string(kerr_params.mass));
    FOR(idir) {
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
  ModifiedCCZ4RHS<
      FourDerivScalarTensor<CouplingAndPotential>, ModifiedPunctureGauge,
      FourthOrderDerivatives>::modified_params_t modified_ccz4_params;
  KerrBH::params_t kerr_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
