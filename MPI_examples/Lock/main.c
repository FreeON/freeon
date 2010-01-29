#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

int
main (int argc, char **argv)
{
  int rank, size;
  int *error_codes;
  int i;
  char *buffer;
  int *flag;
  int all_done;
  MPI_Comm intercomm;
  MPI_Request *request;
  MPI_Status *status;
  struct timespec nanoseconds;
  char **child_argv;

  MPI_Init(&argc, &argv);

  /* Get some information. */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0)
  {
    error_codes = (int*) malloc(sizeof(int)*size);
    buffer = (char*) malloc(sizeof(char)*size);
    request = (MPI_Request*) malloc(sizeof(MPI_Request)*size);
    status = (MPI_Status*) malloc(sizeof(MPI_Status)*size);
    flag = (int*) malloc(sizeof(int)*size);

    /* Construct arguments for children. */
    child_argv = (char**) malloc(sizeof(char*)*(argc+1));
    for (i = 1; i < argc; ++i)
    {
      child_argv[i-1] = (char*) malloc(sizeof(char)*(strlen(argv[i])+1));
      strcpy(child_argv[i-1], argv[i]);
    }
    child_argv[argc-1] = (char*) malloc(sizeof(char)*3);
    snprintf(child_argv[argc-1], 3, "%i", size);
    child_argv[argc] = NULL;

    while (1)
    {
      printf("[master] spawning %i processes\n", size);
      MPI_Comm_spawn("./other", child_argv, size, MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomm, error_codes);

      /* Wait for children to finish. */
      for (i = 0; i < size; ++i)
      {
        MPI_Irecv(buffer, 1, MPI_CHAR, MPI_ANY_SOURCE, 1, intercomm, &request[i]);
        flag[i] = 0;
      }

      all_done = 0;
      while (! all_done)
      {
        all_done = 1;
        for (i = 0; i < size; ++i)
        {
          if (! flag[i])
          {
            MPI_Test(&request[i], &flag[i], &status[i]);
            if (! flag[i])
            {
              all_done = 0;
            }
          }
        }

        /* Sleep a little. */
        nanoseconds.tv_sec = 0;
        nanoseconds.tv_nsec = 0.5e9;
        nanosleep(&nanoseconds, NULL);
      }

      printf("[master] all children are done\n");
    }
  }

  printf("[master (rank: %i)] waiting at barrier\n", rank);
  MPI_Barrier(MPI_COMM_WORLD);
  printf("[master (rank: %i)] done\n", rank);

  MPI_Finalize();
}
