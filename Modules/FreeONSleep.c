#include "config.h"

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

void
freeonsleep_single (float *sleeptime)
{
  struct timespec nanoseconds;

  /* Split the sleeptime into seconds and nanoseconds. */
  nanoseconds.tv_sec = (time_t) floor(*sleeptime);
  nanoseconds.tv_nsec = (long) floor((*sleeptime-nanoseconds.tv_sec)*1e9);

  if (nanosleep(&nanoseconds, NULL) != 0)
  {
    printf("[FreeONSleep] could not sleep, dying out of sleep deprevation: %s\n", strerror(errno));
    exit(1);
  }

//#ifdef HAVE_USLEEP
//  useconds_t microseconds = (useconds_t) floor((*sleeptime-seconds)*1e6);
//#endif

//if (sleep(seconds) != 0)
//{
//  printf("[FreeONSleep] could not sleep, dying out of sleep deprevation.\n");
//  exit(1);
//}

//#ifdef HAVE_USLEEP
//  if (usleep(microseconds) != 0)
//  {
//    printf("[FreeONSleep] could not microsleep, dying out of sleep deprevation.\n");
//    exit(1);
//  }
//#else
//  if (seconds == 0)
//  {
//    /* Sleep a little bit more. */
//    if (sleep(1) != 0)
//    {
//      printf("[FreeONSleep] could not sleep, dying out of sleep deprevation.\n");
//      exit(1);
//    }
//  }
//#endif
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
