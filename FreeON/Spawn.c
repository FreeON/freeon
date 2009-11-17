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

/*   FORKS A CHILD PROCESS    */
/*   Author: Matt Challacombe */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>

int
spawn_ (int *nc, int *maxlen, int *ichr)
{
  int i, j, k, status;
  char** argv;
  pid_t pid, wpid;

  int ZERO_ERROR=0;
  int EXIT_ERROR=-120384;
  int FORK_ERROR=-320498;
  int DUMP_ERROR=-580234;
  int SGNL_ERROR=-674034;
  int MISC_ERROR=-843503;
  int ierr=MISC_ERROR;

  /* Allocate argument list. */
  argv = (char**) malloc(sizeof(char*)*(*nc+1));

  for(i = 0, k = 0; i < *nc; i++)
  {
    argv[i] = (char*) malloc(sizeof(char*)*(*maxlen+1));
    for(j = 0; j < *maxlen; j++, k++)
    {
      argv[i][j] = (char) ichr[k];
      if (argv[i][j] == ' ') { argv[i][j] = '\0'; }
    }
    argv[i][*maxlen] = '\0';
  }

  /* Terminate argument list. */
  argv[*nc] = NULL;

  /* Fork process. */
  pid = fork();

  if(pid == 0)
  {
    /* This is the child process. */
    status = execvp(argv[0], argv);
    printf("error with execvp: %s\n", strerror(errno));
    exit(status);
  }

  else if(pid < 0)
  {
    /* Error on fork(). */
    printf("unable to fork() \"%s\".\n",argv[0]);
    for(i = 0; i < *nc; i++)
    {
      printf("argv[%i] = \"%s\"\n", i, argv[i]);
    }
    ierr = FORK_ERROR;
  }

  /* This is the parent process. */
  printf("waiting for child to finish\n");
  wpid = waitpid(pid, &status, 0);
  printf("child has returned\n");

  if (WEXITSTATUS(status) == 0)
  {
    ierr = 0;
  }

  else
  {
    if(WSTOPSIG(status)!=0) ierr = SGNL_ERROR;
    else if(WCOREDUMP(status)!=0) ierr = DUMP_ERROR;
    else ierr = EXIT_ERROR;
  }

  /* Free memory. */
  for (i = 0; i < *nc+1; ++i) { free(argv[i]); }
  free(argv);

  printf("returning from spawn\n");
  return ierr;
}

int
spawn (int *nc, int *maxlen, int *ichr)
{
  return spawn_(nc, maxlen, ichr);
}
