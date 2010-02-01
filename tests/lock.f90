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
  CALL InitMPI()

  ! Get some information on worl communicator.
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, IErr)
  pid = GetPIDWrapper()
  CALL MondoLog(DEBUG_NONE, "lock "//TRIM(IntToChar(pid)), "starting up...", "rank "//TRIM(IntToChar(MyID)))

  ! Allocate lock.
  lock%Alloc = ALLOCATED_FALSE
  CALL AllocateLock(lock, MPI_COMM_WORLD)

  CALL MondoLog(DEBUG_NONE, "lock "//TRIM(IntToChar(pid)), "trying to get lock", "rank "//TRIM(IntToChar(MyID)))
  CALL AcquireLock(lock, FreeONLockExclusive)
  CALL MondoLog(DEBUG_NONE, "lock "//TRIM(IntToChar(pid)), "sleeping", "rank "//TRIM(IntToChar(MyID)))
  CALL FreeONSleep(5)
  CALL MondoLog(DEBUG_NONE, "lock "//TRIM(IntToChar(pid)), "relasing lock", "rank "//TRIM(IntToChar(MyID)))
  CALL ReleaseLock(lock)

  ! Wait for other ranks.
  CALL MondoLog(DEBUG_NONE, "lock "//TRIM(IntToChar(pid)), "waiting for other ranks to finish", "rank "//TRIM(IntToChar(MyID)))
  CALL MPI_Barrier(MPI_COMM_WORLD, IErr)
  CALL MondoLog(DEBUG_NONE, "lock "//TRIM(IntToChar(pid)), "everybody finished", "rank "//TRIM(IntToChar(MyID)))

  ! Free lock.
  CALL MondoLog(DEBUG_NONE, "lock "//TRIM(IntToChar(pid)), "freeing lock", "rank "//TRIM(IntToChar(MyID)))
  CALL FreeLock(lock)

  CALL FiniMPI()

#endif

END PROGRAM locktest
