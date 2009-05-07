#include "lal.h"

#include <stdlib.h>

void
lal_rand (matrix_t *A)
{
  int i, j;

  for (i = 0; i < A->M; ++i) {
    for (j = 0; j < A->N; ++j)
    {
      lal_set(i, j, rand()/(double)(RAND_MAX), A);
    }
  }
}
