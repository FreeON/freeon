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

/* This code provides some wrappers to get resource limits.
 *
 * Author: Nicolas Bock
 */

#include <sys/time.h>
#include <sys/resource.h>

#include <stdio.h>

void
getstacksizelimit (int * current_limit, int * max_limit)
{
  struct rlimit limits;

  getrlimit (RLIMIT_STACK, & limits);
  if (limits.rlim_cur == RLIM_INFINITY)
  {
    *current_limit = -1;
  }

  else
  {
    *current_limit = limits.rlim_cur;
  }

  if (limits.rlim_cur == RLIM_INFINITY)
  {
    *max_limit = -1;
  }

  else
  {
    *max_limit = limits.rlim_max;
  }
}

void
getstacksizelimit_ (int * current_limit, int * max_limit)
{
  getstacksizelimit (current_limit, max_limit);
}

void
getstacksizelimit__ (int * current_limit, int * max_limit)
{
  getstacksizelimit (current_limit, max_limit);
}
