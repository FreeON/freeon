#include "MondoConfig.h"

PROGRAM locktest

#if defined(PARALLEL_CLONES)
  USE MondoMPI
  USE DerivedTypes
  USE MondoLogger
  USE Parse
  USE Utilities

  IMPLICIT NONE

  TYPE(FreeONLock) :: lock
  INTEGER          :: rank
  INTEGER          :: IErr
  INTEGER          :: pid

  ! Initialize MPI.
  CALL MPI_INIT(IErr)

  ! Get some information on world communicator.
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, IErr)
  pid = GetPIDWrapper()
  CALL MondoLog(DEBUG_NONE, "lock "//TRIM(IntToChar(pid)), "starting up...", "rank "//TRIM(IntToChar(rank)))

  ! Allocate lock.
  lock%Alloc = ALLOCATED_FALSE
  CALL AllocateLock(lock, MPI_COMM_WORLD)

  CALL MondoLog(DEBUG_NONE, "lock "//TRIM(IntToChar(pid)), "trying to get lock", "rank "//TRIM(IntToChar(rank)))
  !CALL AcquireLock(lock, FreeONLockExclusive)
  CALL AcquireLock(lock, FreeONLockShared)
  CALL MondoLog(DEBUG_NONE, "lock "//TRIM(IntToChar(pid)), "sleeping", "rank "//TRIM(IntToChar(rank)))
  CALL FreeONSleep(10)
  CALL MondoLog(DEBUG_NONE, "lock "//TRIM(IntToChar(pid)), "releasing lock", "rank "//TRIM(IntToChar(rank)))
  CALL ReleaseLock(lock)

  ! Wait for other ranks.
  CALL MondoLog(DEBUG_NONE, "lock "//TRIM(IntToChar(pid)), "waiting for other ranks to finish", "rank "//TRIM(IntToChar(rank)))
  CALL MPI_Barrier(MPI_COMM_WORLD, IErr)
  CALL MondoLog(DEBUG_NONE, "lock "//TRIM(IntToChar(pid)), "everybody finished", "rank "//TRIM(IntToChar(rank)))

  ! Free lock.
  CALL MondoLog(DEBUG_NONE, "lock "//TRIM(IntToChar(pid)), "freeing lock", "rank "//TRIM(IntToChar(rank)))
  CALL MondoLog(DEBUG_NONE, "lock "//TRIM(IntToChar(pid)), "communicator = "//TRIM(IntToChar(lock%communicator)), "rank "//TRIM(IntToChar(rank)))
  CALL FreeLock(lock)

  CALL MPI_FINALIZE(IErr)

#endif

END PROGRAM locktest
