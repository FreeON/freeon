#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>

int
main (int argc, char **argv)
{
  int rank;
  int size;
  int *error_codes;
  int spawn_counter = 0;
  char *slave_argv[] = { "arg1", "arg2", 0 };
  char *buffer;
  MPI_Request *request;
  MPI_Status *status;
  int *flag;
  int buffer_index;
  int number_children;
  int i;
  MPI_Comm spawn;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0)
  {
    printf("[master] running on %i processors\n", size);

    while (1)
    {
      printf("[master] (%i) forking processes\n", spawn_counter++);
      error_codes = (int*) malloc(sizeof(int)*size);
      MPI_Comm_spawn("./slave", slave_argv, size, MPI_INFO_NULL, 0, MPI_COMM_SELF, &spawn, error_codes);
      MPI_Comm_remote_size(spawn, &number_children);
      printf("[master] waiting for messages from %i children\n", number_children);

      /* Prepare receive buffers. */
      buffer = (char*) malloc(sizeof(char)*number_children);
      request = (MPI_Request*) malloc(sizeof(MPI_Request)*number_children);
      flag = (int*) malloc(sizeof(int)*number_children);
      status = (MPI_Status*) malloc(sizeof(MPI_Status)*number_children);
      for (buffer_index = 0; buffer_index < number_children; ++buffer_index)
      {
        MPI_Irecv(buffer+buffer_index, 1, MPI_CHAR, MPI_ANY_SOURCE, 0, spawn, request+buffer_index);
        *(flag+buffer_index) = 0;
      }

      buffer_index = 0;
      while (buffer_index < number_children)
      {
        for (i = 0; i < number_children; ++i)
        {
          if (! *(flag+i))
          {
            MPI_Test(request+i, flag+i, status+i);
            if (*(flag+i))
            {
              printf("[master] child %i finished\n", i+1);
              buffer_index++;
            }
          }
        }

        sleep(2);
      }

      printf("[master] everybody finished, freeing memory\n");
      free(status);
      free(flag);
      free(request);
      free(buffer);
      free(error_codes);
    }
  }

  MPI_Finalize();
}
