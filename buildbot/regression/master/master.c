#include "config.h"

#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>

#define MAX_LENGTH 10000

int
invoke (char *program, char *path, char **argv)
{
  pid_t pid;
  int status;
  char program_path[MAX_LENGTH];

  printf("[master:invoke] starting %s in %s\n", program, path);

  /* Construct absolute path to program. */
  snprintf(program_path, MAX_LENGTH, "%s/%s", path, program);

  pid = fork();
  if (pid == 0)
  {
    /* This is the child process.
     *
     * First try to use the hardcoded path, i.e. the install path of the
     * binaries. If that fails, fall back to using simply the binary name,
     * relying on the PATH environment variable and the shell.
     */
    printf("[master:invoke] invoking %s\n", program_path);
    if ((status = execvp(program_path, argv)) < 0)
    {
      printf("[master:invoke] error with execvp: %s\n", strerror(errno));

      printf("[master:invoke] invoking %s\n", program);
      if ((status = execvp(program, argv)) < 0)
      {
        printf("[master:invoke] error with execvp: %s\n", strerror(errno));
        return 1;
      }
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
