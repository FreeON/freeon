#include <unistd.h>

void
freeonsleep (int *time)
{
  sleep(*time);
}

void
freeonsleep_ (int *time)
{
  freeonsleep(time);
}

void
freeonsleep__ (int *time)
{
  freeonsleep(time);
}
