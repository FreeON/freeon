/* Interface to library.
 */

#include "matrix.h"

int
mal_allocate (const int M, const int N, matrix_t **A);

double
mal_get (const int i, const int j, const matrix_t *A);

void
mal_set (const int i, const int j, const double Aij, matrix_t *A);

void
mal_zero (matrix_t *A);

void
mal_dgemm (const char *transA, const char *transB, const int M, const int N,
    const int K, const double alpha, const matrix_t A, const int lda,
    const matrix_t B, const int ldb, const double beta, matrix_t C,
    const int ldc);
