#include <lal.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int
main ()
{
  int M, N;
  matrix_t *A = NULL;

  M = (int) (100+rand()/((double) (RAND_MAX))*100);
  N = (int) (100+rand()/((double) (RAND_MAX))*100);

  return lal_allocate(M, N, &A);
}
