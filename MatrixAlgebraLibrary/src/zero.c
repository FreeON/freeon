#include "mal.h"

void
mal_zero (matrix_t *A)
{
  int i, j;

  for (i = 0; i < A->M; ++i) {
    for (j = 0; j < A->N; ++j)
    {
      mal_set(i, j, 0, A);
    }
  }
}
