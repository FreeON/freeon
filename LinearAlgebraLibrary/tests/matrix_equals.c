#include <lal.h>

#include <stdlib.h>

int
main ()
{
  int M = 10;
  int N = 10;

  int i, j;

  matrix_t *A = NULL;
  matrix_t *B = NULL;

  lal_allocate(M, N, &A);
  lal_allocate(M, N, &B);

  lal_rand(A);
  lal_zero(B);

  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j)
    {
      lal_set(i, j, lal_get(i, j, A), B);
    }
  }

  return lal_equals(A, B);
}
