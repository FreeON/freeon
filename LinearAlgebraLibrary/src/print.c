#include "lal.h"

#include <stdio.h>

void
lal_print (matrix_t *A)
{
  int i, j;

  for (i = 0; i < A->M; ++i) {
    for (j = 0; j < A->N; ++j)
    {
      printf("%f", lal_get(i, j, A));
      if (j < A->N-1) { printf(" "); }
    }
    printf("\n");
  }
}
