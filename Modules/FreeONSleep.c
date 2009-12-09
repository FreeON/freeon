#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void
freeonsleep_single (float *sleeptime)
{
  int seconds = (int) floor(*sleeptime);
  int microseconds = (int) floor((*sleeptime-seconds)*1e6);

  if (sleep(seconds) != 0)
  {
    printf("[FreeONSleep] could not sleep, dying out of sleep deprevation.\n");
    exit(1);
  }
  if (usleep(microseconds) != 0)
  {
    printf("[FreeONSleep] could not microsleep, dying out of sleep deprevation.\n");
    exit(1);
  }
}

void
freeonsleep_single_ (float *sleeptime)
{
  freeonsleep_single(sleeptime);
}

void
freeonsleep_single__ (float *sleeptime)
{
  freeonsleep_single(sleeptime);
}

void
freeonsleep_integer (int *sleeptime_int)
{
  float sleeptime = (float) *sleeptime_int;
  freeonsleep_single(&sleeptime);
}

void
freeonsleep_integer_ (int *sleeptime_int)
{
  freeonsleep_integer(sleeptime_int);
}

void
freeonsleep_integer__ (int *sleeptime_int)
{
  freeonsleep_integer(sleeptime_int);
}
