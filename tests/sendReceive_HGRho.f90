#include "MondoConfig.h"

PROGRAM SendReceive

#if defined(PARALLEL_CLONES)
  USE DerivedTypes
  USE MondoMPI
  USE PrettyPrint

  INTEGER     :: rank, size
  INTEGER     :: IErr
  TYPE(HGRho) :: A, B

  CALL MPI_Init(IErr)

  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, IErr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, size, IErr)

  ! Allocate object.
  CALL Initialize(A)
  CALL Initialize(B)

  ! Allocate array.
  CALL New(A, (/ 3, 3, 3, 1 /))

  IF(rank == 0) THEN
    ! Receive objects from other ranks.
    CALL PPrint(A, "A", Unit_O = 6)
    DO rank = 1, size-1
      CALL Recv(B, rank, 1, comm_O = MPI_COMM_WORLD)
      CALL PPrint(B, "B", Unit_O = 6)
    ENDDO
  ELSE
    CALL Send(A, 0, 1, comm_O = MPI_COMM_WORLD)
  ENDIF

  CALL MPI_Finalize(IErr)
#endif

END PROGRAM SendReceive
