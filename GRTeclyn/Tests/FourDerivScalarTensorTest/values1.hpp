/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Reference RHS values for the FourDerivScalarTensorTest probe cell.
//
// This file is #included inside run_four_deriv_scalar_tensor_test(), where
//   std::array<amrex::Real, NUM_VARS> known;   // zero-initialised
//   bool have_reference;
// are in scope. Fill one "known[c_<var>] = <number>;" line per component and
// set have_reference = true.
//
// HOW TO REGENERATE
// -----------------
// 1. Build GRChombo/Tests/FourDerivScalarTensorGridTest (needs a GRChombo +
//    GRFolres checkout with CHOMBO_HOME / GRCHOMBO_SOURCE set), e.g.
//       cd GRFolres/GRChombo/Tests/FourDerivScalarTensorGridTest && make all
// 2. Run the executable; it prints a "known[c_...] = ...;" block followed by
//    "have_reference = true;" for the same grid, parameters and constant
//    coupling used by this test.
// 3. Replace everything below this line with that output.
//
// Until then the test still builds and runs: it reports the computed RHS and
// only checks that it contains no NaNs.

have_reference = false;
// --- paste GRChombo output below and delete the line above ------------------
