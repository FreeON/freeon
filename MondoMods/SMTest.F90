!    COMPUTES THE MIN AND MAX (PREFETCHED) MFLOP RATE OF THE SMALL MATRIX 
!    DGEMM LIBRARIES CREATED BY PHiPAC FOR VAROUS BLOCK SIZES AND THE CURRENT PLATFORM.
!    Author: Matt Challacombe
!-------------------------------------------------------------------------
PROGRAM SMTest
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE Parse
  USE MemMan
  USE Order
  IMPLICIT NONE
  TYPE(TIME)                     :: Perf
  INTEGER                        :: I,J,K,N,N2,NBloks
  INTEGER, DIMENSION(3)          :: BlokSize=(/3,9,15/)
  INTEGER, PARAMETER             :: NNon0=10000*15**2
  REAL(DOUBLE),DIMENSION(NNon0)  :: A,B,C
  CHARACTER(LEN=6),PARAMETER     :: Prog='SMTest'
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: String
!----------------------------------------------------------------------------------------------------------- 
!
!
  A=One; B=Two;C=Zero
  DO I=1,3
     N=BlokSize(I); N2=N**2; NBloks=NNon0/N**2-1
!
     CALL Elapsed_TIME(Perf,'Init')
     J=1
     DO K=1,NBloks 
        CALL DGEMM_NN(N,N,N,One,A(J),B(J),C(J))
        J=J+N2
     ENDDO 
     CALL Elapsed_TIME(Perf,'Accum')
     Perf%FLOP=DBLE(NBloks*N2*(2*N-1))
     String='('//TRIM(IntToChar(N))//'x'//TRIM(IntToChar(N))//')'
     CALL PPrint(Perf,String,Unit_O=6)
!
     CALL Elapsed_TIME(Perf,'Init')
     J=1
     DO K=1,NBloks 
        CALL DGEMM_NN(N,N,N,One,A(1),B(1),C(1))
        J=J+N2
     ENDDO 
     CALL Elapsed_TIME(Perf,'Accum')
     Perf%FLOP=DBLE(NBloks*N2*(2*N-1))
     String='('//TRIM(IntToChar(N))//'x'//TRIM(IntToChar(N))//')[prefetch]'
     CALL PPrint(Perf,String,Unit_O=6)
!
  ENDDO

END PROGRAM
