#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

void
temporarydirectory (char *path, int *max_length)
{
  int result = 0;

#if defined HAVE_MKOSTEMP
  result = mkostemp(path, 0);
#elif defined HAVE_MKSTEMP
  result = mkstemp(path);
#endif
  if (result < 0)
  {
    printf("[temporarydirectory] error creating temporary directory: %s", strerror(errno));
    exit(1);
  }
}
