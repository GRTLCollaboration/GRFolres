/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TESTCASES_HPP_
#define TESTCASES_HPP_

// doctest header
#include "doctest.h"

// AMReX includes
#include <AMReX.H>

// Test cases
#include "FourDerivScalarTensorTest.hpp"

TEST_CASE("FourDerivScalarTensor") { run_four_deriv_scalar_tensor_test(); }

#endif /* TESTCASES_HPP_ */
