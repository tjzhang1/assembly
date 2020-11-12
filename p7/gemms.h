#ifndef GEMMS_H
   #define GEMMS_H

void dgemm_a1b0
   (int ta, int tb, int M, int N, int K, double alpha, const double *A,
   int lda, const double *B, int ldb, double beta, double *C, int ldc);
void dgemm_a1b1
   (int ta, int tb, int M, int N, int K, double alpha, const double *A,
   int lda, const double *B, int ldb, double beta, double *C, int ldc);
void dgemm_a1bX
   (int ta, int tb, int M, int N, int K, double alpha, const double *A,
   int lda, const double *B, int ldb, double beta, double *C, int ldc);
void dgemm_aXb0
   (int ta, int tb, int M, int N, int K, double alpha, const double *A,
   int lda, const double *B, int ldb, double beta, double *C, int ldc);
void dgemm_aXb1
   (int ta, int tb, int M, int N, int K, double alpha, const double *A,
   int lda, const double *B, int ldb, double beta, double *C, int ldc);
void dgemm_aXbX
   (int ta, int tb, int M, int N, int K, double alpha, const double *A,
   int lda, const double *B, int ldb, double beta, double *C, int ldc);

void dgemm
   (int ta, int tb, int M, int N, int K, double alpha, const double *A,
   int lda, const double *B, int ldb, double beta, double *C, int ldc);

typedef void (*DGEMM_PTR)
   (int ta, int tb, int M, int N, int K, double alpha, const double *A,
   int lda, const double *B, int ldb, double beta, double *C, int ldc);

#endif
