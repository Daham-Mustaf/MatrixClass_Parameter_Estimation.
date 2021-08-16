#ifndef f77blas_h_
#define f77blas_h_
#include "openblas_config.h"
#define f77blas_dgemm dgemm_
#define f77blas_dpotrf dpotrf_
#define f77blas_dtrsm dtrsm_
extern "C" {
void f77blas_dpotrf(const char &UL, const int &N, double *A, const int &lda,
                    int &info);
void f77blas_dgemm(const char &transA, const char &transB, const int &M,
                   const int &N, const int &K, const double &alpha, const double *A,
                   const int &lda, const double *B, const int &ldb,
                   const double &beta, double *C, const int &ldc);
void f77blas_dtrsm(const char &LR, const char &UL, const char &transA,
                   const char &UN, const int &M, const int &N,
                   const double &alpha, double *A, const int &lda, double *B,
                   const int &ldb);
}
#endif