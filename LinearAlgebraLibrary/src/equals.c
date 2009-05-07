#include "lal.h"

/* Compare 2 matrices and return:
 *
 * 0  - if the 2 are identical.
 * -1 - if they are different.
 */

int
lal_equals (matrix_t *A, matrix_t *B)
{
  int i, j;

  if (A->M != B->M) { return -1; }
  if (A->N != B->N) { return -1; }

  for (i = 0; i < A->M; ++i) {
    for (j = 0; j < A->N; ++j)
    {
      if (lal_get(i, j, A) != lal_get(i, j, B)) { return -1; }
    }
  }

  return 0;
}
