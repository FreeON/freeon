#include "mal.h"

#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  matrix_t *A = NULL;

  mal_allocate(10, 10, &A);

  mal_zero(A);

  if (mal_get(0, 0, A) == 0.0) { return 0; }
  else { return 1; }
}
