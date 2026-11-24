/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BINARYBH4DSTLEVEL_HPP_
#define BINARYBH4DSTLEVEL_HPP_

#include "BinaryBH4dSTAmr.hpp"
#include "CouplingAndPotential.hpp"
#include "DefaultLevelBld.hpp"
#include "FourDerivScalarTensor.hpp"
#include "GRAmrLevel.hpp"
#include "PunctureTracker.hpp"

/// Evolution level for a binary black hole in 4-derivative scalar-tensor
/// (shift-symmetric Einstein-scalar-Gauss-Bonnet) gravity.
class BinaryBH4dSTLevel : public GRAmrLevel
{
  public:
    static void variableSetUp();

    // Inherit the constructors from GRAmrLevel
    using GRAmrLevel::GRAmrLevel;

    static constexpr int num_punctures = 2;

    // The 4dST theory, templated over the coupling-and-potential function.
    // (The derivative order is selected at runtime in specific_eval_rhs.)
    template <class deriv_t = FourthOrderDerivatives>
    using FourDerivScalarTensorWithCouplingAndPotential =
        FourDerivScalarTensor<CouplingAndPotential, deriv_t>;

    BinaryBH4dSTAmr *get_bh_amr_ptr();

    /// Get a reference to the PunctureTracker object stored by BHAmr
    PunctureTracker<num_punctures> &get_puncture_tracker();

    /// Things to do at every full timestep
    void specific_advance() override;

    /// Initial data calculation
    void initData() override;

    /// Calculation of the right hand side for the time stepping
    void specific_eval_rhs(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                           amrex::Real a_time) override;

    /// Things to do after dt*rhs has been added to the solution
    void specific_update_ode(amrex::MultiFab &a_soln) override;

    /// Things to do after each time step (extraction, tracking, diagnostics)
    void specific_post_timestep() override;

    /// Things to do before tagging cells for regridding
    void pre_tag_cells() final;

    /// Tag cells for regridding
    void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                   amrex::Real a_regrid_threshold) final;

    //! Things to do after a restart
    void specific_post_restart() override;

    //! Things to do after init
    void specific_post_init() override;

    //! Things to do after writing a plotfile
    void specific_post_plotfile(const std::string &a_dir,
                                std::ostream & /*a_os*/) override;

    //! Things to do after writing a checkpoint
    void specific_post_checkpoint(const std::string &a_dir,
                                  std::ostream & /*a_os*/) override;
};

#endif /* BINARYBH4DSTLEVEL_HPP_ */
