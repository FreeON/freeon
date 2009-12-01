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
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "config.h"

double
random_wrapper ()
{
  static int first_time = 1;
  static unsigned int seed = 0;

  if (first_time == 1)
  {
#ifdef HAVE_RANDOM_SEED
    seed = time(NULL);
#else
    seed = 0;

    /* Print out a warning. */
    printf("[Random] seed random number generator with 0. This makes the random number\n");
    printf("         generator produce a reproducible sequence or random numbers.\n");
#endif
    srand(seed);
    first_time = 0;
  }

  return rand()/(double) RAND_MAX;
}

double
random_ ()
{
  return random_wrapper();
}

double
random__ ()
{
  return random_wrapper();
}
