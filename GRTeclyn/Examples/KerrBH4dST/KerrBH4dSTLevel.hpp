/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef KERRBH4DSTLEVEL_HPP_
#define KERRBH4DSTLEVEL_HPP_

// Evolutionlevel for a Kerr black hole in 4dST.
#include "CouplingAndPotential.hpp"
#include "DefaultLevelBld.hpp"
#include "FourDerivScalarTensor.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRAmrLevel.hpp"
#include "KerrBH4dSTAmr.hpp"
#include "SixthOrderDerivatives.hpp"

class KerrBH4dSTLevel : public GRAmrLevel
{
  public:
    using GRAmrLevel::GRAmrLevel;

    template <class deriv_t = FourthOrderDerivatives>
    using KerrBH4dSTWithCouplingAndPotential =
        FourDerivScalarTensor<CouplingAndPotential, deriv_t>;

    KerrBH4dSTAmr *get_kerrbh4dst_amr_ptr();

    static void variableSetUp();

    void specific_advance() override;

    void initData() override;

    void specific_eval_rhs(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                           amrex::Real a_time) override;

    void specific_update_ode(amrex::MultiFab &a_soln) override;

    void specific_post_timestep() override;

    void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                   amrex::Real a_regrid_threshold) final;
};

#endif /* KERRBH4DSTLEVEL_HPP_ */
