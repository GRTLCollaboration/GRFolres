/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef LINEARSOLVER_HPP_
#define LINEARSOLVER_HPP_

class LinearSolver {
public:
  LinearSolver() {}
  template <class data_t>
  static void solve_linear_system(const int N, data_t *LHS, data_t *RHS);
};

//#include "LinearSolver.impl.hpp"

#endif /* LINEARSOLVER_HPP_ */
