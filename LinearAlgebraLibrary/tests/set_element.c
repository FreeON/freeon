#include <lal.h>

#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  matrix_t *A = NULL;

  lal_allocate(10, 10, &A);

  lal_zero(A);

  lal_set(1, 2, 5, A);

  if (lal_get(1, 2, A) == 5.0) { return 0; }
  else { return 1; }
}
