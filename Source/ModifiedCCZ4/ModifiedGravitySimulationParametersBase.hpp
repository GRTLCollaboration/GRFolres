/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MODIFIEDGRAVITYSIMULATIONPARAMETERSBASE_HPP_
#define MODIFIEDGRAVITYSIMULATIONPARAMETERSBASE_HPP_

// General includes
#include "BoundaryConditions.hpp"
#include "FilesystemTools.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRParmParse.hpp"
#include "ModifiedCCZ4RHS.hpp"
#include "ModifiedPunctureGauge.hpp"
#include "SimulationParametersBase.hpp"
#include "SphericalExtraction.hpp"

//! Class to handle the simulations params that are always required in
//! simulations with a fixed background
template <class theory_t>
class ModifiedGravitySimulationParametersBase : public SimulationParametersBase
{
  public:
    ModifiedGravitySimulationParametersBase(GRParmParse &pp)
        : SimulationParametersBase(pp)
    {
        read_params(pp);
        check_params();
    }

  private:
    void read_params(GRParmParse &pp)
    {
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

        // Scalar extraction
        if (activate_extraction)
        {
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
        }
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
    }

  public:
    // Collection of parameters necessary for the scalar extraction
    spherical_extraction_params_t scalar_extraction_params;

    typename ModifiedCCZ4RHS<theory_t, ModifiedPunctureGauge,
                             FourthOrderDerivatives>::modified_params_t
        modified_ccz4_params;
};

#endif /* MODIFIEDGRAVITYSIMULATIONPARAMETERSBASE_HPP_ */
