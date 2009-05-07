#include "mal.h"

#include <assert.h>

double
mal_get (const int i, const int j, const matrix_t *A)
{
  assert(i >= 0 && i < A->M);
  assert(j >= 0 && j < A->N);

  return A->data[j+i*A->N];
}
