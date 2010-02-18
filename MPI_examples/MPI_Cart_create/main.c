#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>

int
main (int argc, char **argv)
{
  int rank, cart_rank, size;
  int *dims, *periods, *coords, *remain_dims;
  MPI_Comm cart_comm;
  MPI_Comm mondo_comm;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  dims = (int*) malloc(sizeof(int)*2);
  periods = (int*) malloc(sizeof(int)*2);
  coords = (int*) malloc(sizeof(int)*2);
  remain_dims = (int*) malloc(sizeof(int)*2);

  dims[0] = 1;
  dims[1] = size;

  periods[0] = 0;
  periods[1] = 0;

  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart_comm);

  if (cart_comm != MPI_COMM_NULL)
  {
    MPI_Cart_coords(cart_comm, rank, 2, coords);
  }

  else
  {
    printf("dead node\n");
    exit(1);
  }

  remain_dims[0] = 1;
  remain_dims[1] = 0;

  MPI_Cart_sub(cart_comm, remain_dims, &mondo_comm);

  MPI_Comm_rank(mondo_comm, &cart_rank);

  printf("[rank %i] coords = [ %i, %i ], cartesian rank = %i\n", rank, coords[0], coords[1], cart_rank);

  MPI_Finalize();
}
