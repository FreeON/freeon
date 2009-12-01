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

#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

/* This code creates a unique scratch directory.
 *
 * Author: Nicolas Bock
 */

void
temporarydirectory (char *path, int *length)
{
  int i;
  int result = 0;
  char* temp = NULL;

  temp = (char*) malloc(sizeof(char)*((*length)+1));
  for (i = 0; i < *length; ++i) { temp[i] = path[i]; }
  temp[*length] = '\0';

#if defined HAVE_MKOSTEMP
  result = mkostemp(temp, 0);
#elif defined HAVE_MKSTEMP
  result = mkstemp(temp);
#else
  result = -1;
  printf("[temporarydirectory] I have no working mkstemp() or mkostemp()\n");
  free(temp);
  temp = (char*) malloc(sizeof(char)*(strlen(HAVE_MONDO_SCRATCH)+1));
  if(strcpy(temp, HAVE_MONDO_SCRATCH) == NULL)
  {
    printf("[temporarydirectory] error copying HAVE_MONDO_SCRATCH\n");
    exit(1);
  }
#endif
  if (result < 0)
  {
    printf("[temporarydirectory] error creating temporary file (%s): %s\n", temp, strerror(errno));
  }

  /* Create directory. */
  if (unlink(temp) != 0)
  {
    printf("[temporarydirectory] error unlinking temporary file (%s): %s\n", temp, strerror(errno));
    exit(1);
  }

  if (mkdir(temp, S_IRUSR | S_IWUSR | S_IXUSR) != 0)
  {
    printf("[temporarydirectory] error creating temporary directory (%s): %s\n", temp, strerror(errno));
    exit(1);
  }

  /* Copy result back. */
  for (i = 0; i < *length; ++i) { path[i] = temp[i]; }

  /* Free memory. */
  free(temp);
}

void
temporarydirectory_ (char *path, int *max_length)
{
  temporarydirectory(path, max_length);
}

void
temporarydirectory__ (char *path, int *max_length)
{
  temporarydirectory(path, max_length);
}
