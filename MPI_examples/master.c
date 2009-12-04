#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int
main (int argc, char **argv)
{
  int rank;
  int size;
  int *error_codes;
  int spawn_counter = 0;
  char *slave_argv[] = { "arg1", "arg2", 0 };
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
      printf("[master] waiting at barrier\n");
      MPI_Barrier(spawn);
      free(error_codes);
    }
  }

  MPI_Finalize();
}
