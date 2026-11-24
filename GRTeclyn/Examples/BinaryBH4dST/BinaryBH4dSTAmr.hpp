/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BINARYBH4DSTAMR_HPP_
#define BINARYBH4DSTAMR_HPP_

#include "BHAmr.hpp"
#include "ParticleInterpolator.hpp"

//! AMR hierarchy for the BinaryBH4dST example.
/*!
 * Inherits the puncture tracker and the Weyl4 interpolator from BHAmr and adds
 * a separate interpolator used for the spherical extraction of the scalar
 * field.
 */
class BinaryBH4dSTAmr : public BHAmr<2>
{
  public:
    static constexpr int scalar_num_components = 1;
    ParticleInterpolator<scalar_num_components> m_scalar_interpolator;

    explicit BinaryBH4dSTAmr(amrex::LevelBld *a_level_bld)
        : BHAmr<2>(a_level_bld)
    {
    }

    void init(amrex::Real a_start_time, amrex::Real a_stop_time) override
    {
        BHAmr<2>::init(a_start_time, a_stop_time);

        m_scalar_interpolator.setup(this);
    }
};

#endif /* BINARYBH4DSTAMR_HPP_ */
