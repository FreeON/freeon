#include "lal.h"

#include <assert.h>

void
lal_set (const int i, const int j, const double Aij, matrix_t *A)
{
  assert(i >= 0 && i < A->M);
  assert(j >= 0 && j < A->N);

  A->data[j+i*A->N] = Aij;
}
