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
  int workload;
  char buffer;
  MPI_Comm spawn;
  FILE *random_dev;
  int random_seed;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /* Create some random workload. */
  random_dev = fopen("/dev/urandom", "r");
  fgets((char*) &random_seed, sizeof(int), random_dev);
  srand(random_seed);
  fclose(random_dev);

  workload = 20000+(rand()/(double) RAND_MAX-0.5)*3000;
  printf("[slave %i] working with workload %i\n", rank, workload);

  for (i = 0; i < workload; ++i) {
    for (j = 0; j < 100000; ++j)
    {
      temp = rand();
    }
  }

  printf("[slave %i] sending message to master\n", rank);
  MPI_Comm_get_parent(&spawn);
  MPI_Send(&buffer, 1, MPI_CHAR, 0, 0, spawn);

  MPI_Finalize();
}
