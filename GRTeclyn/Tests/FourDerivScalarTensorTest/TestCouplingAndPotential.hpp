/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TESTCOUPLINGANDPOTENTIAL_HPP_
#define TESTCOUPLINGANDPOTENTIAL_HPP_

#include "ScalarFieldVars.hpp"
#include <AMReX_REAL.H>

// Constant coupling / potential used only by FourDerivScalarTensorTest.
//
// Every quantity is a fixed number so that the modified-CCZ4 + 4dST RHS is a
// closed expression of the (polynomial) grid data, and the GRTeclyn result can
// be compared against the GRChombo one produced with the *same* constants (see
// GRChombo/Tests/FourDerivScalarTensorGridTest/TestCouplingAndPotential.hpp).
// Keep the values below in sync with the GRChombo copy.
namespace TestCouplingConstants
{
static constexpr amrex::Real dfdphi   = 0.21649827531440613;
static constexpr amrex::Real d2fdphi2 = 0.12356664969697423;
static constexpr amrex::Real g2       = 0.794128266835908;
static constexpr amrex::Real dg2dphi  = 0.523505716367042;
static constexpr amrex::Real V_of_phi = 0.7841125845230876;
static constexpr amrex::Real dVdphi   = 0.18899643989058168;
} // namespace TestCouplingConstants

class TestCouplingAndPotential
{
  public:
    TestCouplingAndPotential() = default;

    // NB: GRTeclyn's argument order is
    //   (dfdphi, d2fdphi2, V_of_phi, dVdphi, g2, dg2dphi)
    // which differs from GRChombo's - the GRChombo copy of this class orders
    // the arguments to match GRChombo. The assigned values are identical.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    compute_coupling_and_potential(amrex::Real &dfdphi, amrex::Real &d2fdphi2,
                                   amrex::Real &V_of_phi, amrex::Real &dVdphi,
                                   amrex::Real &g2, amrex::Real &dg2dphi,
                                   const ScalarFieldVars & /*vars*/)
    {
        dfdphi   = TestCouplingConstants::dfdphi;
        d2fdphi2 = TestCouplingConstants::d2fdphi2;
        V_of_phi = TestCouplingConstants::V_of_phi;
        dVdphi   = TestCouplingConstants::dVdphi;
        g2       = TestCouplingConstants::g2;
        dg2dphi  = TestCouplingConstants::dg2dphi;
    }
};

#endif /* TESTCOUPLINGANDPOTENTIAL_HPP_ */
