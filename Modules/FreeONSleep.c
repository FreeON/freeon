#include "config.h"

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

void
freeonsleep (float sleeptime)
{
  struct timespec nanoseconds;

  /* Split the sleeptime into seconds and nanoseconds. */
  nanoseconds.tv_sec = (time_t) floor(sleeptime);
  nanoseconds.tv_nsec = (long) floor((sleeptime-nanoseconds.tv_sec)*1e9);

  if (nanosleep(&nanoseconds, NULL) != 0)
  {
    printf("[FreeONSleep] could not sleep, dying out of sleep deprevation: %s\n", strerror(errno));
    exit(1);
  }
}

void
F77_FUNC(freeonsleep_single, FREEONSLEEP_SINGLE) (float *sleeptime)
{
  freeonsleep(*sleeptime);
}

void
F77_FUNC(freeonsleep_integer, FREEONSLEEP_INTEGER) (int *sleeptime_int)
{
  float sleeptime = (float) *sleeptime_int;
  freeonsleep(sleeptime);
}
