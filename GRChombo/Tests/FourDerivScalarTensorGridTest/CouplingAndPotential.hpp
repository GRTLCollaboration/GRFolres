/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COUPLINGANDPOTENTIAL_HPP_
#define COUPLINGANDPOTENTIAL_HPP_

#include "simd.hpp"

// Constant coupling / potential for FourDerivScalarTensorGridTest.
//
// The numbers below must stay identical to
//   GRTeclyn/Tests/FourDerivScalarTensorTest/TestCouplingAndPotential.hpp
// (note the different argument order of compute_coupling_and_potential between
// the two codes - only the argument order differs, not the values).
namespace TestCouplingConstants
{
static constexpr double dfdphi   = 0.21649827531440613;
static constexpr double d2fdphi2 = 0.12356664969697423;
static constexpr double g2       = 0.794128266835908;
static constexpr double dg2dphi  = 0.523505716367042;
static constexpr double V_of_phi = 0.7841125845230876;
static constexpr double dVdphi   = 0.18899643989058168;
} // namespace TestCouplingConstants

class CouplingAndPotential
{
  public:
    struct params_t
    {
    };

    CouplingAndPotential() = default;
    CouplingAndPotential(const params_t & /*a_params*/) {}

    //! GRChombo argument order: (dfdphi, d2fdphi2, g2, dg2dphi, V_of_phi, dVdphi)
    template <class data_t, template <typename> class vars_t>
    void compute_coupling_and_potential(data_t &dfdphi, data_t &d2fdphi2,
                                        data_t &g2, data_t &dg2dphi,
                                        data_t &V_of_phi, data_t &dVdphi,
                                        const vars_t<data_t> & /*vars*/,
                                        const Coordinates<data_t> & /*coords*/)
        const
    {
        dfdphi   = TestCouplingConstants::dfdphi;
        d2fdphi2 = TestCouplingConstants::d2fdphi2;
        g2       = TestCouplingConstants::g2;
        dg2dphi  = TestCouplingConstants::dg2dphi;
        V_of_phi = TestCouplingConstants::V_of_phi;
        dVdphi   = TestCouplingConstants::dVdphi;
    }
};

#endif /* COUPLINGANDPOTENTIAL_HPP_ */
