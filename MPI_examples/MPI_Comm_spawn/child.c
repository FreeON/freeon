#include <mpi.h>
#include <stdio.h>
#include <unistd.h>

int
main (int argc, char **argv)
{
  int rank;
  char buffer[1];
  MPI_Comm parent;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_get_parent(&parent);

  printf("[slave (%i)] starting up, sleeping...\n", rank);
  sleep(5);
  printf("[slave (%i)] done sleeping, signalling master\n", rank);
  MPI_Send(buffer, 1, MPI_CHAR, 0, 1, parent);

  MPI_Finalize();
}
