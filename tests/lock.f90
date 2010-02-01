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

  CALL InitMPI()

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, IErr)
  CALL MondoLog(DEBUG_NONE, "lock", "starting up...", "rank "//TRIM(IntToChar(MyID)))

  lock%Alloc = ALLOCATED_FALSE
  CALL AllocateLock(lock, MPI_COMM_WORLD)

  CALL MondoLog(DEBUG_NONE, "lock", "trying to get lock", "rank "//TRIM(IntToChar(MyID)))
  CALL AcquireLock(lock, FreeONLockExclusive)
  CALL MondoLog(DEBUG_NONE, "lock", "sleeping", "rank "//TRIM(IntToChar(MyID)))
  CALL FreeONSleep(5)
  CALL MondoLog(DEBUG_NONE, "lock", "relasing lock", "rank "//TRIM(IntToChar(MyID)))
  CALL ReleaseLock(lock)

  CALL FreeLock(lock)

  CALL FiniMPI()

#endif

END PROGRAM locktest
