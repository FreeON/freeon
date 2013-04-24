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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

/* The number of bytes read in one chunk when copying a file. */
#define CHUNKSIZE 1024

/* This function copies a file. It does not check whether the target file exists
 * already and overwrites it in case it does.
 *
 * Author: Nicolas Bock <nbock@lanl.gov>
 */

void
F77_FUNC(filecopywrapper, FILECOPYWRAPPER) (int *lenA, char *fileA_arg, int *lenB, char *fileB_arg)
{
  char *buffer = NULL;
  char *fileA  = NULL;
  char *fileB  = NULL;

  size_t bytes_read;
  size_t bytes_written;

  FILE *fdA;
  FILE *fdB;

  int i;

  /* Fix strings. */
  fileA = (char*) malloc(sizeof(char)*(*lenA+1));
  fileB = (char*) malloc(sizeof(char)*(*lenB+1));

  for (i = 0; i < *lenA; ++i)
  {
    fileA[i] = fileA_arg[i];
  }

  for (i = 0; i < *lenB; ++i)
  {
    fileB[i] = fileB_arg[i];
  }

  fileA[*lenA] = '\0';
  fileB[*lenB] = '\0';

  /* Read file and copy. */
  fdA = fopen(fileA, "r");
  if (fdA == NULL)
  {
    printf("[filecopy] error opening file %s\n", fileA);
    printf("[filecopy] %s\n", strerror(errno));
    free(fileA);
    free(fileB);
    return;
  }

  fdB = fopen(fileB, "w");
  if (fdB == NULL)
  {
    printf("[filecopy] error opening file %s\n", fileB);
    printf("[filecopy] %s\n", strerror(errno));
    free(fileA);
    free(fileB);
    return;
  }

  /* Allocate buffer memory. */
  buffer = (char*) malloc(CHUNKSIZE);
  while (1)
  {
    bytes_read = fread((void*) buffer, 1, CHUNKSIZE, fdA);
    if (bytes_read < CHUNKSIZE && ferror(fdA))
    {
      printf("[filecopy] error reading from %s: %s\n", fileA, strerror(errno));
      break;
    }

    if (bytes_read > 0)
    {
      bytes_written = fwrite((void*) buffer, 1, bytes_read, fdB);
      if (bytes_written != bytes_read)
      {
        printf("[filecopy] I read %u bytes, but wrote %u bytes\n", (unsigned int) bytes_read, (unsigned int) bytes_written);
        if (ferror(fdB))
        {
          printf("[filecopy]   %s\n", strerror(errno));
        }
        break;
      }
    }

    if (ferror(fdB) || feof(fdA)) { break; }
  }

  /* Release memory. */
  free(buffer);

  /* Free memory for filenames. */
  free(fileA);
  free(fileB);

  /* Close files. */
  fclose(fdA);
  fclose(fdB);
}
