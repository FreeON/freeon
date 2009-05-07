#include "lal.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

int
lal_allocate (const int M, const int N, matrix_t **A)
{
  assert(M > 0);
  assert(N > 0);

  /* Allocate space for matrix object structure. */
  *A = (matrix_t*) (malloc(sizeof(matrix_t)));

  /* Set dimensions. */
  (*A)->M = M;
  (*A)->N = N;

  /* Allocate space for matrix data. */
  (*A)->data = (double*) (malloc(sizeof(double)*M*N));

  if ((*A)->data == NULL)
  {
    fprintf(stderr, "[lal_allocate]", "error allocating matrix\n");
    return -1;
  }

  else { return 0; }
}
