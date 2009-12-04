#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>

int
main (int argc, char **argv)
{
  int rank;
  int size;
  int i, j;
  double temp;
  MPI_Comm spawn;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  printf("[slave %i] working\n", rank);
  for (i = 0; i < 10000; ++i) {
    for (j = 0; j < 500000; ++j)
    {
      temp = rand();
    }
  }

  printf("[slave %i] waiting at barrier\n", rank);
  MPI_Comm_get_parent(&spawn);
  MPI_Barrier(spawn);

  MPI_Finalize();
}
