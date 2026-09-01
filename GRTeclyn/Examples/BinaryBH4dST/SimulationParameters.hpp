/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "BaseParameterChecker.hpp"

// Problem specific includes:
#include "BoostedBHInitialData.hpp"
#include "CCZ4RHS.hpp"
#include "ExtractionTagger.hpp"
#include "FourthOrderDerivatives.hpp"
#include "PunctureTagger.hpp"
#include "PunctureTracker.hpp"
#include "SphericalExtractionParameters.hpp"
#ifdef USE_TWOPUNCTURES
#include "TwoPuncturesInitialData.hpp"
#endif

#include "CouplingAndPotential.hpp"
#include "FourDerivScalarTensor.hpp"
#include "InitialScalarData.hpp"
#include "ModifiedPunctureGauge.hpp"

class SimulationParameters
{
  public:
    // NOLINTNEXTLINE(readability-identifier-length)
    SimulationParameters() = delete;

    static void check_params()
    {
        BaseParameterChecker::check_params();

        CCZ4_params_t::check_params();
        ModifiedPunctureGauge<FourthOrderDerivatives>::params_t::check_params();
        ExtractionTagger::check_params();
        PunctureTagger<2>::check_params();
        puncture_tracker_params_t::check_params();

#ifndef USE_TWOPUNCTURES
        BoostedBHInitialData::params_t::check_params(1);
        BoostedBHInitialData::params_t::check_params(2);
#else
        TwoPuncturesInitialData::check_params();
#endif

        // 4dST theory parameters
        CouplingAndPotential::params_t::check_params();
        FourDerivScalarTensor<CouplingAndPotential>::params_t::check_params();
        InitialScalarData::params_t::check_params();

        spherical_extraction_params_t::check_params("weyl_extraction");
        spherical_extraction_params_t::check_params("scalar_extraction");
    }
};

#endif /* SIMULATIONPARAMETERS_HPP */
