/*
     This code is part of the MondoSCF suite of programs for linear scaling
     electronic structure theory and ab initio molecular dynamics.

     Copyright (2004). The Regents of the University of California. This
     material was produced under U.S. Government contract W-7405-ENG-36
     for Los Alamos National Laboratory, which is operated by the University
     of California for the U.S. Department of Energy. The U.S. Government has
     rights to use, reproduce, and distribute this software.  NEITHER THE
     GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
     OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.

     This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by the
     Free Software Foundation; either version 2 of the License, or (at your
     option) any later version. Accordingly, this program is distributed in
     the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE. See the GNU General Public License at www.gnu.org for details.

     While you may do as you like with this software, the GNU license requires
     that you clearly mark derivative software.  In addition, you are encouraged
     to return derivative works to the MondoSCF group for review, and possible
     disemination in future releases.
*/

/* This code provides some wrappers to get resource usage.
 *
 * Author: Nicolas Bock
 */

#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

/* getmemoryusage:
 *
 * This function returns the virtual memory size of the current process in
 * kBytes.
 */
void
getmemoryusage_wrapper (int * size)
{
  pid_t our_pid;
  FILE * file = NULL;
  char path [2048];
  char token [256];
  unsigned long temp;
  int i;

  /* Get PID. */
  our_pid = getpid();

  /* Open stat in the proc filesystem. */
  sprintf(path, "/proc/%i/stat", (int) our_pid);
  file = fopen(path, "r");

  if (file != NULL)
  {
    for (i = 0; i < 22; ++i)
    {
      if (fscanf(file, "%s", & token[0]) <= 0)
      {
        printf("[getmemoryusage] error reading\n");
        exit(1);
      }
    }
    if (fscanf(file, "%lu", & temp) <= 0)
    {
      printf("[getmemoryusage] error reading\n");
      exit(1);
    }
    *size = (int) round(temp/1024.);

    /* Close file again. */
    fclose(file);
  }

  else { *size = -1; }
}

void
getmemoryusage_wrapper_ (int * size)
{
  getmemoryusage_wrapper(size);
}

void
getmemoryusage_wrapper__ (int * size)
{
  getmemoryusage_wrapper(size);
}

double
getcputime ()
{
  struct rusage current_time;

  if (getrusage(RUSAGE_SELF, &current_time) < 0)
  {
    printf("[getCPUTime] error\n");
    exit(1);
  }

  return current_time.ru_utime.tv_sec+current_time.ru_utime.tv_usec/(double) 1e6;
}

double
getcputime_ ()
{
  return getcputime();
}

double
getcputime__ ()
{
  return getcputime();
}
