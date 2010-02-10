#include "MondoConfig.h"

PROGRAM SendReceive

#if defined(PARALLEL_CLONES)
  USE DerivedTypes
  USE MondoMPI

  INTEGER               :: rank, size
  INTEGER               :: IErr
  INTEGER               :: i, j
  INTEGER, DIMENSION(2) :: NLower, NUpper
  TYPE(DBL_RNK2)        :: A, B

  CALL MPI_Init(IErr)

  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, IErr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, size, IErr)

  ! Allocate object.
  CALL Initialize(A)
  CALL Initialize(B)

  ! Set dimensions.
  NLower(1) = 1
  NUpper(1) = 3
  NLower(2) = 1
  NUpper(2) = 6

  ! Allocate array.
  CALL New(A, NUpper, M_O = NLower)

  ! Fill object with random data.
  DO i = NLower(1), NUpper(1)
    DO j = NLower(2), NUpper(2)
      A%D(i, j) = (i-1)+(j-1)*(NUpper(2)-NLower(2)+1)
    ENDDO
  ENDDO

  IF(rank == 0) THEN
    ! Receive objects from other ranks.
    DO rank = 1, size-1
      CALL Recv(B, rank, 1, comm_O = MPI_COMM_WORLD)
      DO i = NLower(1), NUpper(1)
        DO j = NLower(2), NUpper(2)
          IF(ABS(A%D(i, j)-B%D(i, j)) > 1.0D-12) THEN
            WRITE(*,*) "test failed"
            STOP
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ELSE
    CALL Send(A, 0, 1, comm_O = MPI_COMM_WORLD)
  ENDIF

  CALL MPI_Finalize(IErr)
#endif

END PROGRAM SendReceive
