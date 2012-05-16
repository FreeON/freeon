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

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <sys/select.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>

char *
signal_name (int signal_number)
{
  char *name;

  switch(signal_number)
  {
    case 9:
      name = strdup("SIGKILL: Kill signal");
      break;

    case 11:
      name = strdup("SIGSEGV: Invalid memory reference");
      break;

    default:
      name = strdup("unknown signal");
      break;
  }

  return name;
}

#define MAX_BUFFER_LENGTH 2000
int
F77_FUNC(spawn, SPAWN) (int *nc, int *maxlen, int *ichr)
{
  int i, j, k, status;
  char** argv;
  pid_t pid, wpid;

  int EXIT_ERROR=-120384;
  int FORK_ERROR=-320498;
  int DUMP_ERROR=-580234;
  int SGNL_ERROR=-674034;
  int MISC_ERROR=-843503;
  int ierr=MISC_ERROR;

#ifdef DUMP_PIPE
  /* Create a pipe to standard output and error so we can capture those in the
   * main process. Otherwise we won't get any useful failure message from the
   * front-end.
   */
  char pipe_buffer[MAX_BUFFER_LENGTH];
  int child_pipe[2];

  if ((status = pipe(child_pipe)) != 0)
  {
    printf("error creating pipe for standard output\n");
    exit(status);
  }

  int pipe_read  = child_pipe[0];
  int pipe_write = child_pipe[1];
#endif

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
#ifdef DUMP_PIPE
    /* This is the child process. */
    int stdout_file = fileno(stdout);
    int stderr_file = fileno(stderr);

    /* Replace stdout and stderr with pipe. */
    if ((status = dup2(pipe_write, stdout_file)) < 0)
    {
      printf("error attaching pipe\n");
      exit(status);
    }

    if ((status = dup2(pipe_write, stderr_file)) < 0)
    {
      printf("error attaching pipe\n");
      exit(status);
    }
#endif

    /* Replace this process with backend. */
    status = execvp(argv[0], argv);

    /* If we are here, then execvp could not start. */
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

#ifdef DUMP_PIPE
  /* Dump standard output and error from back-end process. */
  fd_set read_set;
  FD_ZERO(&read_set);
  FD_SET(pipe_read, &read_set);
  struct timeval timeout;
  while (1)
  {
    timeout.tv_sec = 1;
    timeout.tv_usec = 0;
    status = select(pipe_read+1, &read_set, NULL, NULL, &timeout);

    if (status == 0)
    {
      /* We timed out. */
      printf("pipe timed out\n");
      break;
    }

    else if (status > 0)
    {
      /* There is something to read from the pipe. */
      if (read(pipe_read, pipe_buffer, MAX_BUFFER_LENGTH) > 0)
      {
        printf("%s", pipe_buffer);
      }
      else { break; }
    }

    else {
      /* An error occurred. */
      printf("error on select()\n");
      exit(status);
    }
  }
#endif

  /* This is the parent process. */
  if ((wpid = waitpid(pid, &status, 0)) != pid)
  {
    printf("[Spawn] error calling wpid() (which returned %i)\n", wpid);
    exit(1);
  }

  if (WIFSIGNALED(status))
  {
    printf("[Spawn] child (%s) terminated by signal %i (%s)\n", argv[0], WTERMSIG(status), signal_name(WTERMSIG(status)));
    ierr = SGNL_ERROR;
#ifdef WCOREDUMP
    if (WCOREDUMP(status))
    {
      ierr = DUMP_ERROR;
    }
#endif
  }

  else if (WIFSTOPPED(status))
  {
    printf("[Spawn] child (%s) stopped by signal %i (%s)\n", argv[0], WSTOPSIG(status), signal_name(WSTOPSIG(status)));
    ierr = SGNL_ERROR;
  }

  else if (WIFEXITED(status))
  {
    ierr = EXIT_ERROR;

    if (WEXITSTATUS(status) == 0)
    {
      ierr = 0;
    }

    else
    {
      printf("[Spawn] child (%s) exited due to general error (return code %i)\n", argv[0], WEXITSTATUS(status));
    }
  }

  else
  {
    printf("[Spawn] FIXME\n");
    exit(1);
  }

  /* Free memory. */
  for (i = 0; i < *nc+1; ++i) { free(argv[i]); }
  free(argv);

  return ierr;
}
