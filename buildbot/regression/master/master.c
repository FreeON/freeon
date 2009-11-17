#include "config.h"

#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>

int
invoke (char *program, char *path, char **argv)
{
  pid_t pid;
  int status;

  printf("[master:invoke] starting %s in %s\n", program, path);

  pid = fork();
  if (pid == 0)
  {
    /* This is the child process. */
    if ((status = execvp(program, argv)) < 0)
    {
      printf("[master:invoke] error with execvp: %s\n", strerror(errno));
      return 1;
    }
  }
  printf("[master:invoke] this is the parent process\n");

  /* This is the parent process. */
  if (wait(&status) != pid)
  {
    printf("[master:invoke] error waiting for child\n");
    return 1;
  }

  if (WIFEXITED(status))
  {
    if (WEXITSTATUS(status) == 0)
    {
      printf("[master:invoke] normal exit of child\n");
    }

    else
    {
      printf("[master:invoke] abnormal exit of child: %i\n", WEXITSTATUS(status));
      return WEXITSTATUS(status);
    }
  }

  else
  {
    printf("[master:invoke] error\n");
    return 1;
  }

  return 0;
}

int
main (int argc, char **argv)
{
  char *exe_path = EXEPATH;

  if (invoke("slave_1", exe_path, argv) != 0) { return 1; }
  if (invoke("slave_2", exe_path, argv) != 0) { return 1; }
  if (invoke("slave_3", exe_path, argv) != 0) { return 1; }

  return 0;
}
