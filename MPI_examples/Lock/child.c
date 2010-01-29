/* Exclusive lock with MPI.
 *
 * This is adapted from:
 *
 * ``Implementing MPI-IO Atomic Mode and Shared File Pointers Using MPI
 * One-Sided Communication'', Robert Latham, Robert Ross and Rajeev Thakur,
 * International Journal of High Performance Computing Applications 21, 132
 * (2007).
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define LOCKTAG 1

int
main (int argc, char **argv)
{
  int rank, cartrank;
  int size;
  int i;
  char buffer[1];
  int dimensions[2];
  int periods[2];
  int remaining_dimensions[2];
  MPI_Comm parent;
  MPI_Comm cartcomm;
  MPI_Comm subcomm;
  MPI_Info win_info;
  MPI_Win win;
  int *waitflag = NULL;
  int *array_of_blocklengths = NULL;
  int *array_of_displacements = NULL;
  MPI_Datatype waitflagcopy_t;
  int *waitflagcopy = NULL;
  int val;
  int nextrank;

  MPI_Init(&argc, &argv);

  /* Get some information. */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_get_parent(&parent);

  /* Create a cartesian communicator. */
  printf("[child (rank: %i)] called with \"", rank);
  for (i = 0; i < argc; ++i)
  {
    printf("%s", argv[i]);
    if (i < argc-1) { printf(" "); }
  }
  printf("\"\n");
  fflush(stdout);

  dimensions[0] = 1;
  dimensions[1] = strtol(argv[argc-1], NULL, 10);

  periods[0] = 0;
  periods[1] = 0;

  remaining_dimensions[0] = 1;
  remaining_dimensions[1] = 0;

  MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, periods, 1, &cartcomm);

  /* Partition cartesian communicator into subgroups. */
  MPI_Cart_sub(cartcomm, remaining_dimensions, &subcomm);

  /* Get new rank within cartesian subgroup. */
  MPI_Comm_rank(subcomm, &cartrank);

  /* Print something. */
  printf("[child (rank: %i, cart: %i)] starting up...\n", rank, cartrank);
  fflush(stdout);

  /* Get window object. */
  MPI_Info_create(&win_info);
  if (rank == 0)
  {
    waitflag = (int*) malloc(sizeof(int)*size);
    for (i = 0; i < size; ++i)
    {
      waitflag[i] = 0;
    }
    MPI_Win_create((void*) waitflag, size, sizeof(int), win_info, MPI_COMM_WORLD, &win);
  }

  else
  {
    MPI_Win_create(NULL, 0, 1, win_info, MPI_COMM_WORLD, &win);
  }

  /* Create local copy of waitflag array. This array contains size-1 flags, it
   * excludes the flag for rank. */
  array_of_blocklengths = (int*) malloc(sizeof(int)*(size-1));
  array_of_displacements = (int*) malloc(sizeof(int)*(size-1));

  for (i = 0; i < size; ++i)
  {
    array_of_blocklengths[i] = 1;
    if (i < rank)
    {
      array_of_displacements[i] = i;
    }

    else if (i > rank)
    {
      array_of_displacements[i-1] = i;
    }
  }

  /* Create new datatype to get all waitflags except for ours. */
  MPI_Type_indexed(size-1, array_of_blocklengths, array_of_displacements, MPI_INT, &waitflagcopy_t);
  MPI_Type_commit(&waitflagcopy_t);

  /* Allocate memory for local copy of waitflags. */
  waitflagcopy = (int*) malloc(sizeof(int)*(size-1));

  /* Get lock on window. */
  MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);

  /* Get waitflags for everyone but us. */
  MPI_Get((void*) waitflagcopy, size-1, MPI_INT, 0, 0, 1, waitflagcopy_t, win);

  /* Put our waitflag. */
  val = 1;
  MPI_Put(&val, 1, MPI_INT, 0, rank, 1, MPI_INT, win);

  /* Unlock window. */
  MPI_Win_unlock(0, win);

  /* Check whether some other rank already holds the lock. */
  for (i = 0; i < size-1 && waitflagcopy[i] == 0; ++i) {}
  if (i < size-1)
  {
    /* Wait for notification from whoever holds the lock. */
    printf("[child (rank: %i, cart: %i)] waiting to acquire lock...\n", rank, cartrank);
    fflush(stdout);
    MPI_Recv(NULL, 0, MPI_BYTE, MPI_ANY_SOURCE, LOCKTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  printf("[child (rank: %i, cart: %i)] acquired lock...\n", rank, cartrank);
  fflush(stdout);

  /* And sleep. */
  printf("[child (rank: %i, cart: %i)] sleeping...\n", rank, cartrank);
  fflush(stdout);
  sleep(10);

  /* Release lock. */
  MPI_Win_lock(MPI_LOCK_EXCLUSIVE, 0, 0, win);

  /* Get waitflags. */
  MPI_Get((void*) waitflagcopy, size-1, MPI_INT, 0, 0, 1, waitflagcopy_t, win);

  /* Put our waitflag. */
  val = 0;
  MPI_Put(&val, 1, MPI_INT, 0, rank, 1, MPI_INT, win);

  /* Unlock window. */
  MPI_Win_unlock(0, win);

  /* In case someone else is waiting on the lock, notify them that we are
   * done. */
  for (i = 0; i < size-1 && waitflagcopy[i] == 0; ++i) {}
  if (i < size-1)
  {
    for (nextrank = rank; nextrank < size-1 && waitflagcopy[nextrank] == 0; ++nextrank) {}

    /* nextrank is off by one, because we dind't get our own waitflag. */
    if (nextrank < size-1) { nextrank++; }
    else
    {
      for (nextrank = 0; nextrank < rank && waitflagcopy[nextrank] == 0; nextrank++) {}
    }

    /* Notify the other rank. */
    printf("[child (rank: %i, cart: %i)] telling %i to take lock\n", rank, cartrank, nextrank);
    fflush(stdout);
    MPI_Send(NULL, 0, MPI_BYTE, nextrank, LOCKTAG, MPI_COMM_WORLD);
  }

  /* We are done. */
  printf("[child (rank: %i, cart: %i)] done sleeping, signalling master\n", rank, cartrank);
  fflush(stdout);
  MPI_Send(buffer, 1, MPI_CHAR, 0, 1, parent);

  /* Close window. */
  MPI_Win_free(&win);

  MPI_Finalize();
}
