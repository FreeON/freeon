#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

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

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0)
  {
    error_codes = (int*) malloc(sizeof(int)*size);
    buffer = (char*) malloc(sizeof(char)*size);
    request = (MPI_Request*) malloc(sizeof(MPI_Request)*size);
    status = (MPI_Status*) malloc(sizeof(MPI_Status)*size);
    flag = (int*) malloc(sizeof(int)*size);

    while (1)
    {
      printf("[master] spawning %i processes\n", size);
      MPI_Comm_spawn("./other", argv, size, MPI_INFO_NULL, 0, MPI_COMM_SELF, &intercomm, error_codes);

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
    }
  }

  printf("[master (%i)] waiting at barrier\n", rank);
  MPI_Barrier(MPI_COMM_WORLD);
  printf("[master (%i)] done\n", rank);

  MPI_Finalize();
}
