/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

/*The functions implemented in this file calculate the solution of the linear
system of equations A * x=b through the following steps:
1. decompose A into the format A = L * U (where L is lower triangular and U is
upper triangular). This is the so-called LU decomposition
2. solve A * y = b for y using forward substitution - y is stored in b
3. solve A * x = y for x using back substitution - in the calculation we use b
in place of y (allowed due to the fact that we store y in b in the previous
step); then the solution x is stored in b. Notice that due to this "storing"
strategy (which saves a few floating point operations), the variables x and y
need not be defined in the code below.*/

#ifndef LINEARSOLVER_HPP_
#define LINEARSOLVER_HPP_

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void LU_decompose(const int &N, amrex::Real *A)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            *(A + j * N + i) /= *(A + i * N + i);
            for (int k = i + 1; k < N; k++)
            {
                *(A + j * N + k) -= *(A + j * N + i) * (*(A + i * N + k));
            }
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void forward_substitution(const int &N, const amrex::Real *L, amrex::Real *b)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < i; j++)
        {
            *(b + i) -= *(L + i * N + j) * (*(b + j));
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void back_substitution(const int &N, const amrex::Real *U, amrex::Real *b)
{
    for (int i = N - 1; i >= 0; i--)
    {
        for (int j = i + 1; j < N; j++)
        {
            *(b + i) -= *(U + i * N + j) * (*(b + j));
        }
        *(b + i) /= *(U + i * N + i);
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void solve_linear_system(const int &N, amrex::Real *A, amrex::Real *b)
{
    LU_decompose(N, A);
    forward_substitution(N, A, b);
    back_substitution(N, A, b);
}

#endif // LINEARSOLVER_HPP_
