// __has_include requires C++17
// This tests whether or not MKL library was used in compilation
// (for example when using Intel compilers)

#if __has_include("mkl_lapacke.h")
#define USE_MKL
#endif

#ifdef USE_MKL
#include "mkl_lapacke.h"
#else
#include <lapacke.h>
#endif

#include "EsGB.hpp"

template <>
void EsGB_lapack::solve_system_lapack(const int N, double *LHS, double *RHS)
{
    int nrhs = 1;
    int lda = N;
    int ldb = N;
    int ipiv[N];
    int info;

#ifdef USE_MKL
    LAPACKE_dgesv(LAPACK_COL_MAJOR, N, nrhs, &LHS[0], lda, ipiv, &RHS[0], ldb);
#else
    dgesv_(&N, &nrhs, &LHS[0], &lda, ipiv, &RHS[0], &ldb, &info);
#endif

    if (info < 0)
    {
        pout() << "Illegal value in argument " << info << std::endl;
        MayDay::Error("Illegal value for lapack function");
    }
    else if (info > 0)
    {
        pout() << "Solution could not be computed" << std::endl;
        MayDay::Error("Solution could not be computed");
    }
}

template <>
void EsGB_lapack::solve_system_lapack(const int N, simd<double> *LHS,
                                      simd<double> *RHS)
{
    int simd_len = simd_traits<double>::simd_len;

    double LHS_arr[simd_len][N][N];
    double RHS_arr[simd_len][N];

    for (int n1 = 0; n1 < N; ++n1)
    {
        for (int n2 = 0; n2 < N; ++n2)
        {
            double in_arr[simd_len];
            simd<double>::store(in_arr, LHS[n1 * N + n2]);
            for (int s = 0; s < simd_len; ++s)
                LHS_arr[s][n1][n2] = in_arr[s];
        }

        double in_arr[simd_len];
        simd<double>::store(in_arr, RHS[n1]);
        for (int s = 0; s < simd_len; ++s)
            RHS_arr[s][n1] = in_arr[s];
    }

    for (int s = 0; s < simd_len; ++s)
        solve_system_lapack(N, (double *)(&LHS_arr[s][0][0]),
                            (double *)(RHS_arr[s]));

    for (int n1 = 0; n1 < N; ++n1)
    {
        double out_arr[simd_len];
        for (int s = 0; s < simd_len; ++s)
            out_arr[s] = RHS_arr[s][n1];
        RHS[n1] = simd<double>::load(out_arr);
    }
}
