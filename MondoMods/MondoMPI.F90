MODULE MondoMPI
   USE DerivedTypes  
   USE GlobalScalars
   USE GlobalObjects
   USE ProcessControl
   USE MPIInclude
   USE MemMan
   USE Parse
   IMPLICIT NONE   
#ifdef PARALLEL
!    NOTE DIFFERENCES BETWEEN DERIVED TYPES, EG. INT_VECT, WHICH INCLUDE
!    AUXILIARY INFO SUCH AS ARRAY BOUNDS, AND THE ARRAYS THEMSELVES GIVEN
!    BY THE FULLER NAME, EG. INT_VECTOR, WHICH SPECIFIES THE ARRAY ONLY.
!
   INTEGER,SAVE :: MONDO_COMM=MPI_COMM_WORLD
   INTERFACE BCast          ! Wrappers for MPI_BCAST
      MODULE PROCEDURE BCast_DBL_SCLR, BCast_DBL_VECT,  BCast_DBL_RNK2, &
                       BCast_DBL_RNK3, BCast_DBL_RNK4,  BCast_DBL_RNK6, & 
                       BCast_INT_SCLR,                                  &
                       BCast_INT_VECT, BCast_INT_RNK2,  BCast_INT_RNK3, &
                       BCast_INT_RNK4, BCast_CHR_SCLR,  BCast_LOG_SCLR, &
                       BCast_DEBG                        
   END INTERFACE
   INTERFACE AllReduce   ! Wrappers for MPI_ALLREDUCE
      MODULE PROCEDURE AllReduce_DBL_SCLR
   END INTERFACE
   INTERFACE Reduce       ! Wrappers for MPI_REDUCE
      MODULE PROCEDURE Reduce_DBL_SCLR, Reduce_INT_SCLR
   END INTERFACE
   INTERFACE Gather      ! Wrappers for MPI_GATHER and MPI GATHERV
      MODULE PROCEDURE Gather_INT_SCLR, Gather_INT_VECT, &
                       Gather_DBL_SCLR, Gather_DBL_VECT
   END INTERFACE
   INTERFACE Send       ! Wrappers for MPI_SEND (Blocking)
      MODULE PROCEDURE Send_INT_SCLR, Send_INT_VECT, & 
                       Send_DBL_SCLR, Send_DBL_VECT, &
                       Send_INT_VECTOR 
   END INTERFACE
   INTERFACE Recv       ! Wrappers for MPI_RECV (Blocking)
      MODULE PROCEDURE Recv_INT_SCLR, Recv_INT_VECT, &
                       Recv_DBL_SCLR, Recv_DBL_VECT, &
                       Recv_INT_VECTOR 
   END INTERFACE
   INTERFACE ISend       ! Wrappers for MPI_ISEND (Non-Blocking)
      MODULE PROCEDURE ISend_INT_SCLR, ISend_INT_VECT, & 
                       ISend_DBL_SCLR, ISend_DBL_VECT
   END INTERFACE
   INTERFACE IRecv       ! Wrappers for MPI_IRECV (Non-Blocking)
      MODULE PROCEDURE IRecv_INT_SCLR, IRecv_INT_VECT, &
                       IRecv_DBL_SCLR, IRecv_DBL_VECT
   END INTERFACE

   INTERFACE PSpew
      MODULE PROCEDURE PSpew_INT_VECT, PSpew_DBL_VECT
   END INTERFACE
   CONTAINS
!===============================================================
!
!     WRAPPERS FOR MPI INITIALIZATION AND FINALIZATION
!
!===============================================================
!
      SUBROUTINE InitMPI()
         INTEGER :: IErr
         CHARACTER(LEN=7),PARAMETER :: Sub='InitMPI'
!------------------------------------------------
         CALL MPI_INIT(IErr)
         CALL ErrChk(IErr,Sub)
!        Load global MPI variables
         MyID=MRank()
         NPrc=MSize()
         InParallel=.TRUE.
         CALL CreateComm()
      END SUBROUTINE InitMPI

      SUBROUTINE FiniMPI()
         INTEGER                    :: IErr
         CHARACTER(LEN=7),PARAMETER :: Sub='FiniMPI'
         CALL MPI_FINALIZE(IErr)
         CALL ErrChk(IErr,Sub)
      END SUBROUTINE FiniMPI

      FUNCTION MRank()
         INTEGER :: IErr,MRank
         CHARACTER(LEN=5),PARAMETER :: Sub='MRank'
         CALL MPI_COMM_RANK(MONDO_COMM,MRank,IErr)
         CALL ErrChk(IErr,Sub)
      END FUNCTION MRANK

      FUNCTION MSize()
         INTEGER :: IErr,MSize
         CHARACTER(LEN=5),PARAMETER :: Sub='MSize'
         CALL MPI_COMM_SIZE(MONDO_COMM,MSize,IErr)
         CALL ErrChk(IErr,Sub)
      END FUNCTION MSize
!
!===============================================================
!
!     BCAST WRAPPERS
!
!===============================================================
!--------------------------------------------------------BCAST CHARACTERS
      SUBROUTINE BCast_CHR_SCLR(A)
         CHARACTER(LEN=*)                   :: A
         CHARACTER(LEN=1),DIMENSION(LEN(A)) :: LocalChars
         INTEGER                            :: IErr,I,L
         L=LEN(A)
         DO I=1,L
            LocalChars(I)=A(I:I)
         ENDDO
         CALL MPI_BCAST(LocalChars,L,MPI_CHARACTER,ROOT,MONDO_COMM,IErr)
         IF(IErr/=MPI_SUCCESS) &
            CALL HaltMPI(' MPI_BCAST failed in BCast_STR',IErr)
         DO I=1,L
            A(I:I)=LocalChars(I)
         ENDDO
      END SUBROUTINE BCast_CHR_SCLR

      SUBROUTINE BCast_CHR_VECT(A)
         TYPE(CHR_VECT), INTENT(INOUT)  :: A
         INTEGER                        :: I,N
         N=SIZE(A%C)
         DO I=1,N
            CALL BCast_CHR_SCLR(A%C(I))
         ENDDO
      END SUBROUTINE BCast_CHR_VECT
!--------------------------------------------------------BCAST CHARACTERS
      SUBROUTINE BCast_LOG_SCLR(A)
         LOGICAL,INTENT(INOUT)      :: A
         INTEGER                    :: IErr,I
         I=0; IF(A)I=1
         CALL MPI_BCAST(I,1,MPI_INTEGER,ROOT,MONDO_COMM,IErr)
         IF(IErr/=MPI_SUCCESS) &
            CALL HaltMPI(' MPI_BCAST failed in BCast_STR',IErr)
         A=.FALSE.; IF(I==1)A=.TRUE.
      END SUBROUTINE BCast_LOG_SCLR
!--------------------------------------------------------BCAST DOUBLES
      SUBROUTINE BCast_DBL_SCLR(D,Id_0)
         REAL(DOUBLE), INTENT(INOUT) :: D
         INTEGER, OPTIONAL           :: Id_0
         INTEGER                     :: IErr,Id
         ID=ROOT
         IF(PRESENT(Id_0))Id=Id_0
         CALL MPI_BCAST(D,1,MPI_DOUBLE_PRECISION,Id,MONDO_COMM,IErr)
         IF(IErr/=MPI_SUCCESS) &
            CALL HaltMPI(' MPI_BCAST failed in BCast_DBL',IErr)
      END SUBROUTINE BCast_DBL_SCLR

      SUBROUTINE BCast_DBL_VECT(A,N_O)
         TYPE(DBL_VECT),INTENT(INOUT) :: A
         INTEGER,OPTIONAL,INTENT(IN) :: N_O

         INTEGER :: IErr,L
         L=SIZE(A%D); IF(PRESENT(N_O)) L=N_O
         CALL MPI_BCAST(A%D,L,MPI_DOUBLE_PRECISION,ROOT,MONDO_COMM,IErr)
         IF(IErr/=MPI_SUCCESS) &
            CALL HaltMPI(' MPI_BCAST failed in BCast_DBL_VECT',IErr)
      END SUBROUTINE BCast_DBL_VECT

      SUBROUTINE BCast_DBL_RNK2(A)
         TYPE(DBL_RNK2), INTENT(INOUT) :: A
         INTEGER :: L,M
         L=SIZE(A%D,1); M=SIZE(A%D,2)
         CALL CRNDAV(L*M,A%D)
      END SUBROUTINE BCast_DBL_RNK2

      SUBROUTINE BCast_DBL_RNK3(A)
         TYPE(DBL_RNK3), INTENT(INOUT) :: A
         INTEGER :: L,M,N
         L=SIZE(A%D,1); M=SIZE(A%D,2); N=SIZE(A%D,3)
         CALL CRNDAV(L*M*N,A%D)
      END SUBROUTINE BCast_DBL_RNK3

      SUBROUTINE BCast_DBL_RNK4(A)
         TYPE(DBL_RNK4), INTENT(INOUT) :: A
         INTEGER :: L,M,N,K
         L=SIZE(A%D,1); M=SIZE(A%D,2)
         N=SIZE(A%D,3); K=SIZE(A%D,4)
         CALL CRNDAV(L*M*N*K,A%D)
      END SUBROUTINE BCast_DBL_RNK4

      SUBROUTINE BCast_DBL_RNK6(A)
         TYPE(DBL_RNK6), INTENT(INOUT) :: A
         INTEGER :: L
         L=SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3)* &
           SIZE(A%D,4)*SIZE(A%D,5)*SIZE(A%D,6)
         CALL CRNDAV(L,A%D)
      END SUBROUTINE BCast_DBL_RNK6
!--------------------------------------------------------BCAST INTEGERS
      SUBROUTINE BCast_INT_SCLR(A,Id_O)
         INTEGER, INTENT(INOUT) :: A
         INTEGER, OPTIONAL      :: Id_O
         INTEGER                :: IErr,Id
         Id=ROOT; IF(PRESENT(Id_O))Id=Id_O   
         CALL MPI_BCAST(A,1,MPI_INTEGER,Id,MONDO_COMM,IErr)
         IF(IErr/=MPI_SUCCESS) &
            CALL HaltMPI(' MPI_BCAST failed in BCast_INT',IErr)
      END SUBROUTINE BCast_INT_SCLR

      SUBROUTINE BCast_INT_VECT(A,N_O)
         TYPE(INT_VECT),    INTENT(INOUT) :: A
         INTEGER, OPTIONAL, INTENT(IN)    :: N_O
         INTEGER                          :: IErr,N
         N=SIZE(A%I); IF(PRESENT(N_O))N=N_O
         CALL MPI_BCAST(A%I,N,MPI_INTEGER,ROOT,MONDO_COMM,IErr)
         IF(IErr/=MPI_SUCCESS) &
            CALL HaltMPI(' MPI_BCAST failed in BCast_INT_VECT',IErr)
      END SUBROUTINE BCast_INT_VECT

      SUBROUTINE BCast_INT_RNK2(A)
         TYPE(INT_RNK2), INTENT(INOUT) :: A
         INTEGER :: L,M
         L=SIZE(A%I,1) ; M=SIZE(A%I,2)
         CALL CRNIAV(L*M,A%I)
      END SUBROUTINE BCast_INT_RNK2

      SUBROUTINE BCast_INT_RNK3(A)
         TYPE(INT_RNK3), INTENT(INOUT) :: A
         INTEGER :: L,M,N
         L=SIZE(A%I,1) ; M=SIZE(A%I,2); N=SIZE(A%I,3)
         CALL CRNIAV(L*M*N,A%I)
      END SUBROUTINE BCast_INT_RNK3

      SUBROUTINE BCast_INT_RNK4(A)
         TYPE(INT_RNK4),INTENT(INOUT) :: A
         INTEGER :: I,J,L,M
         I=SIZE(A%I,1); J=SIZE(A%I,2)
         L=SIZE(A%I,3) ; M=SIZE(A%I,4)
         CALL CRNIAV(I*J*L*M,A%I)
      END SUBROUTINE BCast_INT_RNK4
!
!=================================================================
!
!     RANK N -> RANK 1 BCAST
!  
!=================================================================
!
      SUBROUTINE CRNIAV(L,A)               ! BCast RANK N Integer As VECTOR
         INTEGER :: L,IErr
         INTEGER, DIMENSION(L), INTENT(INOUT)  :: A
         CALL MPI_BCAST(A,L,MPI_INTEGER,ROOT,MONDO_COMM,IErr)
         IF(IErr/=MPI_SUCCESS) &
            CALL HaltMPI(' MPI_BCAST failed for integer vector ',IErr)
      END SUBROUTINE CRNIAV

      SUBROUTINE CRNDAV(L,A)               ! BCast RANK N Double As VECTOR
         INTEGER :: L,IErr
         REAL(DOUBLE), DIMENSION(L), INTENT(INOUT)  :: A
         CALL MPI_BCAST(A,L,MPI_DOUBLE_PRECISION,ROOT,MONDO_COMM,IErr)
         IF(IErr/=MPI_SUCCESS) &
            CALL HaltMPI(' MPI_BCAST failed for double vector ',IErr)
      END SUBROUTINE CRNDAV
!-------------------------------------------------- 
       SUBROUTINE BCast_DEBG(A)
         TYPE(DEBG) :: A
         CALL BCast(A%Key)
         CALL BCast(A%Chk)
         CALL BCast(A%Mat)
         CALL BCast(A%Set)
         CALL BCast(A%Int)
         CALL BCast(A%Rho)
         CALL BCast(A%Fmt)
       END SUBROUTINE BCast_DEBG
!===============================================================
!
!     REDUCE WRAPPERS
!
!===============================================================
!
      FUNCTION AllReduce_DBL_SCLR(D,Op_O)
         REAL(DOUBLE), INTENT(INOUT)  :: D
         INTEGER,OPTIONAL,INTENT(IN)  :: Op_O
         REAL(DOUBLE)                 :: AllReduce_DBL_SCLR
         INTEGER                      :: Op,IErr
         CHARACTER(LEN=18), PARAMETER :: Sub='AllReduce_DBL_SCLR'
         IF(PRESENT(Op_O))THEN
            Op=Op_O
         ELSE
            Op=MPI_SUM
         ENDIF
         CALL MPI_ALLREDUCE(D,AllReduce_DBL_SCLR,1, &
              MPI_DOUBLE_PRECISION,Op,MONDO_COMM,IErr)
         CALL ErrChk(IErr,Sub)
      END FUNCTION AllReduce_DBL_SCLR
!-------------------------------------------------- DOUBLES
!
!
      FUNCTION Reduce_DBL_SCLR(D,Op_O)
         REAL(DOUBLE), INTENT(INOUT)  :: D
         INTEGER,OPTIONAL,INTENT(IN)  :: Op_O
         REAL(DOUBLE)                 :: Reduce_DBL_SCLR
         INTEGER                      :: Op,IErr
         CHARACTER(LEN=15), PARAMETER :: Sub='Reduce_DBL_SCLR'
         IF(PRESENT(Op_O))THEN
            Op=Op_O
         ELSE
            Op=MPI_SUM
         ENDIF
         CALL MPI_REDUCE(D,Reduce_DBL_SCLR,1,MPI_DOUBLE_PRECISION, &
                         Op,ROOT,MONDO_COMM,IErr)
         CALL ErrChk(IErr,Sub)
      END FUNCTION Reduce_DBL_SCLR
!-------------------------------------------------- INTEGERS
!
!
      FUNCTION Reduce_INT_SCLR(I)
         INTEGER                      :: I,Reduce_INT_SCLR
         INTEGER                      :: IErr
         CHARACTER(LEN=15), PARAMETER :: Sub='Reduce_INT_SCLR'
         CALL MPI_REDUCE(I,Reduce_INT_SCLR,1,MPI_INTEGER, &
                         MPI_SUM,ROOT,MONDO_COMM,IErr)
         CALL ErrChk(IErr,Sub)
      END FUNCTION Reduce_INT_SCLR
!
!======================================================================
!     GATHERS (Note, be carefull to avoid 
!     reference to uninitialised POINTERs
!======================================================================
!
      SUBROUTINE Gather_INT_SCLR(A,B)
         INTEGER,       INTENT(IN)    :: A
         TYPE(INT_VECT),INTENT(INOUT) :: B
         INTEGER                      :: IErr
         CHARACTER(LEN=15), PARAMETER :: Sub='Gather_INT_SCLR'
         IF(MyId==ROOT)THEN 
            CALL MPI_GATHER(  A,1,MPI_INTEGER, &
                            B%I,1,MPI_INTEGER, &
                            ROOT,MONDO_COMM,IErr)
         ELSE
            CALL MPI_GATHER(  A,1,MPI_INTEGER, &        ! Non-root nodes may not have
                              1,1,MPI_INTEGER, &        ! B allocated, need to protect
                            ROOT,MONDO_COMM,IErr)   ! against refrencing a pointer
         ENDIF            
         CALL ErrChk(IErr,Sub)
      END SUBROUTINE Gather_INT_SCLR

      SUBROUTINE Gather_INT_VECT(Send,Recv,NSnd,NRcv,Disp)
         TYPE(INT_VECT),INTENT(IN)    :: Send,NRcv,Disp
         INTEGER,       INTENT(IN)    :: NSnd        
         TYPE(INT_VECT),INTENT(INOUT) :: Recv
         INTEGER                      :: IErr
         CHARACTER(LEN=15), PARAMETER :: Sub='Gather_INT_VECT'
         IF(MyId==ROOT)THEN
            CALL MPI_GATHERV(Send%I,NSnd,         MPI_INTEGER, &
                             Recv%I,NRcv%I,Disp%I,MPI_INTEGER, &
                             ROOT,MONDO_COMM,IErr)
         ELSE
            CALL MPI_GATHERV(Send%I,NSnd,         MPI_INTEGER, &  
                                  1,   1,       1,MPI_INTEGER, &
                             ROOT,MONDO_COMM,IErr)
         ENDIF
         CALL ErrChk(IErr,Sub)
      END SUBROUTINE Gather_INT_VECT

      SUBROUTINE Gather_DBL_SCLR(A,B)
         REAL(DOUBLE),  INTENT(IN)    :: A
         TYPE(DBL_VECT),INTENT(INOUT) :: B
         INTEGER                      :: IErr
         CHARACTER(LEN=15),PARAMETER  :: Sub='Gather_DBL_SCLR'
         IF(MyId==ROOT)THEN 
            CALL MPI_GATHER(  A,1,MPI_DOUBLE_PRECISION, &
                            B%D,1,MPI_DOUBLE_PRECISION, &
                            ROOT,MONDO_COMM,IErr)
         ELSE
            CALL MPI_GATHER(    A,1,MPI_DOUBLE_PRECISION, &
                            1.0D0,1,MPI_DOUBLE_PRECISION, &
                            ROOT,MONDO_COMM,IErr)
         ENDIF
         CALL ErrChk(IErr,Sub)
      END SUBROUTINE Gather_DBL_SCLR

      SUBROUTINE Gather_DBL_VECT(Send,Recv,NSnd,NRcv,Disp)
         TYPE(DBL_VECT),INTENT(IN)    :: Send
         TYPE(INT_VECT),INTENT(IN)    :: NRcv,Disp
         INTEGER,       INTENT(IN)    :: NSnd        
         TYPE(DBl_VECT),INTENT(INOUT) :: Recv
         INTEGER                      :: IErr
         CHARACTER(LEN=15),PARAMETER  :: Sub='Gather_DBL_VECT'
         IF(MyId==ROOT)THEN 
            CALL MPI_GATHERV(Send%D,NSnd,         MPI_DOUBLE_PRECISION, &
                             Recv%D,NRcv%I,Disp%I,MPI_DOUBLE_PRECISION, &
                             ROOT,MONDO_COMM,IErr)
         ELSE
            CALL MPI_GATHERV(Send%D,NSnd,         MPI_DOUBLE_PRECISION, &
                              1.0D0,   1,       1,MPI_DOUBLE_PRECISION, &
                             ROOT,MONDO_COMM,IErr)
         ENDIF
         CALL ErrChk(IErr,Sub)
      END SUBROUTINE Gather_DBL_VECT
!======================================================================
!
!     BLOCKING SENDS AND RECIEVES
!
!---------------------------------------------------------------------
!
!
      SUBROUTINE Send_DBL_SCLR(Snd,To,Tag)
         REAL(DOUBLE),INTENT(IN)      :: Snd
         INTEGER,     INTENT(IN)      :: To,Tag
         INTEGER                      :: IErr
         CHARACTER(LEN=13), PARAMETER :: Sub='Send_DBL_SCLR'
         CALL MPI_SEND(Snd,1,MPI_DOUBLE_PRECISION,To,Tag, &
                       MONDO_COMM,IErr)
         CALL ErrChk(IErr,Sub)
      END SUBROUTINE Send_DBL_SCLR
!---------------------------------------------------------------------
!
!
      SUBROUTINE Send_DBL_VECT(Snd,N,To,Tag,M_O)
         TYPE(DBL_VECT),  INTENT(IN) :: Snd
         INTEGER,         INTENT(IN) :: N,To,Tag
         INTEGER                     :: L,M,IErr
         INTEGER,OPTIONAL,INTENT(IN) :: M_O
         CHARACTER(LEN=13),PARAMETER :: Sub='Send_DBL_SCLR'
         M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
!         WRITE(*,11)M,N,L
!        11 format('m = ',i4,' n = ',i4,' l = ',i4)
         CALL MPI_SEND(Snd%D(M:N),L,MPI_DOUBLE_PRECISION,To,Tag, &
                       MONDO_COMM,IErr)
         CALL ErrChk(IErr,Sub)
      END SUBROUTINE Send_DBL_VECT
!---------------------------------------------------------------------
!
!
      SUBROUTINE Send_INT_SCLR(Snd,To,Tag)
         INTEGER, INTENT(IN)          :: Snd,To,Tag
         INTEGER                      :: IErr
         CHARACTER(LEN=13), PARAMETER :: Sub='Send_INT_SCLR'
         CALL MPI_SEND(Snd,1,MPI_INTEGER,To,Tag, &
                       MONDO_COMM,IErr)
         CALL ErrChk(IErr,Sub)
      END SUBROUTINE Send_INT_SCLR
!---------------------------------------------------------------------
!
!
      SUBROUTINE Send_INT_VECT(Snd,N,To,Tag,M_O)
         TYPE(INT_VECT),  INTENT(IN) :: Snd
         INTEGER,         INTENT(IN) :: N,To,Tag
         INTEGER                     :: L,M,IErr
         INTEGER,OPTIONAL,INTENT(IN) :: M_O
         CHARACTER(LEN=13),PARAMETER :: Sub='Send_INT_VECT'
         M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
!         WRITE(*,*)' send send send send send send '
!         WRITE(*,*)' M = ',M,' N = ',N,' L = ',L
!         WRITE(*,*)' To = ',To,' Tag = ',Tag
!         WRITE(*,*)' SIZE SND = ',SIZE(Snd%I)
!         WRITE(*,*)' Snd = ',Snd%I(M:N)
         CALL MPI_SEND(Snd%I(M:N),L,MPI_INTEGER,To,Tag, &
                       MONDO_COMM,IErr)
         CALL ErrChk(IErr,Sub)
      END SUBROUTINE Send_INT_VECT
!---------------------------------------------------------------------
!
!
      SUBROUTINE Send_INT_VECTOR(Snd,N,To,Tag,M_O)
         INTEGER,DIMENSION(:),INTENT(IN) :: Snd
         INTEGER,             INTENT(IN) :: N,To,Tag
         INTEGER                         :: L,M,IErr
         INTEGER,OPTIONAL,INTENT(IN)     :: M_O
         CHARACTER(LEN=15),PARAMETER     :: Sub='Send_INT_VECTOR'
         M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
         CALL MPI_SEND(Snd(M:N),L,MPI_INTEGER,To,Tag, &
                       MONDO_COMM,IErr)
         CALL ErrChk(IErr,Sub)
      END SUBROUTINE Send_INT_VECTOR
!---------------------------------------------------------------------
!
!
      SUBROUTINE Recv_DBL_SCLR(Rec,From,Tag)
         REAL(DOUBLE),INTENT(INOUT)          :: Rec
         INTEGER,     INTENT(IN)             :: From,Tag
         INTEGER                             :: IErr
         INTEGER, DIMENSION(MPI_STATUS_SIZE) :: Status
         CHARACTER(LEN=13),PARAMETER         :: Sub='Recv_DBL_SCLR'
         CALL MPI_RECV(Rec,1,MPI_DOUBLE_PRECISION,From,Tag, &
                       MONDO_COMM,Status,IErr)
         CALL ErrChk(IErr,Sub)
      END SUBROUTINE Recv_DBL_SCLR
!---------------------------------------------------------------------
!
!
      SUBROUTINE Recv_DBL_VECT(Rec,N,From,Tag,M_O)
         TYPE(DBL_VECT),INTENT(INOUT)        :: Rec
         INTEGER,       INTENT(IN)           :: N,From,Tag
         INTEGER                             :: L,M,IErr
         INTEGER, DIMENSION(MPI_STATUS_SIZE) :: Status
         INTEGER, OPTIONAL, INTENT(IN)       :: M_O
         CHARACTER(LEN=13),PARAMETER         :: Sub='Recv_DBL_VECT'
         M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
         CALL MPI_RECV(Rec%D(M:N),L,MPI_DOUBLE_PRECISION,From,Tag, &
                       MONDO_COMM,Status,IErr)
         CALL ErrChk(IErr,Sub)
      END SUBROUTINE Recv_DBL_VECT

!---------------------------------------------------------------------
!
!
      SUBROUTINE Recv_INT_SCLR(Rec,From,Tag)
         INTEGER, INTENT(INOUT)              :: Rec
         INTEGER, INTENT(IN)                 :: From,Tag
         INTEGER                             :: IErr
         INTEGER, DIMENSION(MPI_STATUS_SIZE) :: Status
         CHARACTER(LEN=13),PARAMETER         :: Sub='Recv_INT_SCLR'
         CALL MPI_RECV(Rec,1,MPI_INTEGER,From,Tag, &
                       MONDO_COMM,Status,IErr)
         CALL ErrChk(IErr,Sub)
      END SUBROUTINE Recv_INT_SCLR
!---------------------------------------------------------------------
!
!
      SUBROUTINE Recv_INT_VECT(Rec,N,From,Tag,M_O)
         TYPE(INT_VECT),INTENT(INOUT)        :: Rec
         INTEGER,       INTENT(IN)           :: N,From,Tag
         INTEGER                             :: L,M,IErr
         INTEGER, DIMENSION(MPI_STATUS_SIZE) :: Status
         INTEGER, OPTIONAL, INTENT(IN)       :: M_O
         CHARACTER(LEN=13),PARAMETER         :: Sub='Recv_INT_VECT'
         M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
!         WRITE(*,*)' recv recv recv recv recv recv recv '
!         WRITE(*,*)' M = ',M,' N = ',N,' L = ',L
!         WRITE(*,*)' From = ',From,' Tag = ',Tag
!         WRITE(*,*)' SIZE RCV = ',SIZE(Rec%I)
         CALL MPI_RECV(Rec%I(M:N),L,MPI_INTEGER,From,Tag, &
                       MONDO_COMM,Status,IErr)
!         WRITE(*,*)' SIZE REC = ',SIZE(Rec%I)
!         WRITE(*,*)' Rec = ',Rec%I(M:N)
!         WRITE(*,*)' IErr = ',IErr
         CALL ErrChk(IErr,Sub)
      END SUBROUTINE Recv_INT_VECT
!---------------------------------------------------------------------
!
!
      SUBROUTINE Recv_INT_VECTOR(Rec,N,From,Tag,M_O)
         INTEGER,DIMENSION(:),INTENT(INOUT)  :: Rec
         INTEGER,             INTENT(IN)     :: N,From,Tag
         INTEGER                             :: L,M,IErr
         INTEGER, DIMENSION(MPI_STATUS_SIZE) :: Status
         INTEGER, OPTIONAL, INTENT(IN)       :: M_O
         CHARACTER(LEN=15),PARAMETER         :: Sub='Recv_INT_VECTOR'
         M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
         CALL MPI_RECV(Rec(M:N),L,MPI_INTEGER,From,Tag, &
                       MONDO_COMM,Status,IErr)
         CALL ErrChk(IErr,Sub)
      END SUBROUTINE Recv_INT_VECTOR
!======================================================================
!
!     WRAPPERS FOR NON-BLOCKING (IMEDIATE) SENDS AND RECIEVES
!
!--------------------------------------------------------------------
!
!
      FUNCTION ISend_DBL_SCLR(Snd,To,Tag)
         REAL(DOUBLE),INTENT(INOUT)     :: Snd
         INTEGER,     INTENT(IN)     :: To,Tag
         INTEGER                     :: ISend_DBL_SCLR,IErr
         CHARACTER(LEN=14),PARAMETER :: Sub='ISend_DBL_SCLR'
         CALL MPI_ISEND(Snd,1,MPI_DOUBLE_PRECISION,To,Tag, &
                        MONDO_COMM,ISend_DBL_SCLR,IErr)
         CALL ErrChk(IErr,Sub)
      END FUNCTION ISend_DBL_SCLR
!--------------------------------------------------------------------
!
!
      FUNCTION ISend_DBL_VECT(Snd,N,To,Tag,M_O)
         TYPE(DBL_VECT),INTENT(INOUT)    :: Snd
         INTEGER,       INTENT(IN)    :: N,To,Tag
         INTEGER                      :: ISend_DBL_VECT,L,M,IErr
         INTEGER,OPTIONAL,INTENT(IN)  :: M_O
         CHARACTER(LEN=14),PARAMETER  :: Sub='ISend_DBL_VECT'
         M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
         CALL MPI_ISEND(Snd%D(M:N),L,MPI_DOUBLE_PRECISION,To,Tag, &
                        MONDO_COMM,ISend_DBL_VECT,IErr)
         CALL ErrChk(IErr,Sub)
      END FUNCTION ISend_DBL_VECT
!--------------------------------------------------------------------
!
!
      FUNCTION ISend_INT_SCLR(Snd,To,Tag)
         INTEGER,INTENT(INOUT)          :: Snd
         INTEGER,INTENT(IN)          :: To,Tag
         INTEGER                     :: ISend_INT_SCLR,IErr
         CHARACTER(LEN=14),PARAMETER :: Sub='ISend_INT_SCLR'
!         WRITE(*,*)' isend isend isend isend isend isend isend isend isend isend '
!         WRITE(*,*)'MyId = ',MyId,' To = ',To,' Tag = ',Tag
         CALL MPI_ISEND(Snd,1,MPI_INTEGER,To,Tag, &
                        MONDO_COMM,ISend_INT_SCLR,IErr)
         CALL ErrChk(IErr,Sub)
      END FUNCTION ISend_INT_SCLR
!--------------------------------------------------------------------
!
!
      FUNCTION ISend_INT_VECT(Snd,N,To,Tag,M_O)
         TYPE(INT_VECT),INTENT(INOUT) :: Snd
         INTEGER,       INTENT(IN)    :: N,To,Tag
         INTEGER                      :: ISend_INT_VECT,L,M,IErr
         INTEGER,OPTIONAL,INTENT(IN)  :: M_O
         CHARACTER(LEN=14),PARAMETER  :: Sub='ISend_INT_VECT'
         M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
!if(myid==root)then
!         WRITE(*,*)' isend isend isend isend isend isend isend isend isend isend '
!         WRITE(*,*)' M = ',M,' N = ',N,' L = ',L
!         WRITE(*,*)'MyId = ',MyId,' To = ',To,' Tag = ',Tag
!         WRITE(*,*)' SIZE SND = ',SIZE(Snd%I)
!         WRITE(*,*)' Snd = ',Snd%I(M:N)
!endif
         CALL MPI_ISEND(Snd%I(M:N),L,MPI_INTEGER,To,Tag, &
                        MONDO_COMM,ISend_INT_VECT,IErr)
         CALL ErrChk(IErr,Sub)
      END FUNCTION ISend_INT_VECT
!-----------------------------------------------------------------------
!
!
      FUNCTION IRecv_DBL_SCLR(Rec,From,Tag)
         REAL(DOUBLE),INTENT(INOUT) :: Rec
         INTEGER,INTENT(IN)         :: From,Tag
         INTEGER                    :: IRecv_DBL_SCLR,IErr
         CHARACTER(LEN=14)          :: Sub='IRecv_DBL_SCLR'
         CALL MPI_IRECV(Rec,1,MPI_DOUBLE_PRECISION,From,Tag, &
                        MONDO_COMM,IRecv_DBL_SCLR,IErr)
         CALL ErrChk(IErr,Sub)
      END FUNCTION IRecv_DBL_SCLR
!-----------------------------------------------------------------------
!
!
      FUNCTION IRecv_DBL_VECT(Rec,N,From,Tag,M_O)
         TYPE(DBL_VECT),INTENT(INOUT)  :: Rec
         INTEGER,INTENT(IN)            :: From,N,Tag
         INTEGER                       :: IRecv_DBL_VECT,L,M,IErr
         INTEGER,INTENT(IN),OPTIONAL   :: M_O
         CHARACTER(LEN=14)             :: Sub='IRecv_DBL_VECT'
         M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
         CALL MPI_IRECV(Rec%D(M:N),L,MPI_DOUBLE_PRECISION,From,Tag, &
                        MONDO_COMM,IRecv_DBL_VECT,IErr)
         CALL ErrChk(IErr,Sub)
      END FUNCTION IRecv_DBL_VECT
!-----------------------------------------------------------------------
!
!
      FUNCTION IRecv_INT_SCLR(Rec,From,Tag)
         INTEGER,INTENT(INOUT) :: Rec
         INTEGER,INTENT(IN)    :: From,Tag
         INTEGER               :: IRecv_INT_SCLR,IErr
         CHARACTER(LEN=14)     :: Sub='IRecv_INT_SCLR'
!         WRITE(*,*)' irecv irecv irecv irecv irecv irecv irecv irecv irecv irecv irecv '
!         WRITE(*,*)' MyId = ',MyId,' From = ',From,' Tag = ',Tag
         CALL MPI_IRECV(Rec,1,MPI_INTEGER,From,Tag, &
                        MONDO_COMM,IRecv_INT_SCLR,IErr)
         CALL ErrChk(IErr,Sub)
      END FUNCTION IRecv_INT_SCLR
!-----------------------------------------------------------------------
!
!
      FUNCTION IRecv_INT_VECT(Rec,N,From,Tag,M_O)
         TYPE(INT_VECT),INTENT(INOUT)  :: Rec
         INTEGER,INTENT(IN)            :: From,N,Tag
         INTEGER                       :: IRecv_INT_VECT,L,M,IErr
         INTEGER,INTENT(IN),OPTIONAL   :: M_O
         CHARACTER(LEN=14)             :: Sub='IRecv_INT_VECT'
         M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
!if(myid==root)then
!         WRITE(*,*)' irecv irecv irecv irecv irecv irecv irecv irecv irecv irecv irecv '
!         WRITE(*,*)' M = ',M,' N = ',N,' L = ',L
!         WRITE(*,*)' MyId = ',MyId,' From = ',From,' Tag = ',Tag
!         WRITE(*,*)' SIZE RCV = ',SIZE(Rec%I)
!endif
         CALL MPI_IRECV(Rec%I(M:N),L,MPI_INTEGER,From,Tag, &
                        MONDO_COMM,IRecv_INT_VECT,IErr)
         CALL ErrChk(IErr,Sub)
      END FUNCTION IRecv_INT_VECT




!===============================================================
!
!     WRAPER FOR CART_CREATE ROUTINE
      
      SUBROUTINE CreateComm()
         INTEGER                       :: IErr
         CHARACTER(LEN=10)             :: Sub='CreateComm'
         CALL MPI_CART_CREATE(MPI_COMM_WORLD,1,NPrc,.FALSE.,.TRUE.,MONDO_COMM,IErr)
         CALL ErrChk(IErr,Sub)
      END SUBROUTINE CreateComm
!===============================================================
!
!     WRAPERS FOR WAIT ROUTINES
!
      SUBROUTINE WaitAll(Reqs,Mssg)
         TYPE(INT_VECT),INTENT(INOUT)    :: Reqs     
         TYPE(INT_RNK2)                  :: Stats
         INTEGER                         :: I,J,N,IErr
         CHARACTER(LEN=*), OPTIONAL      :: Mssg
         CHARACTER(LEN=8)                :: Sub='MWaitAll'
!------------------------------------------------------------
         IF(PRESENT(Mssg))THEN
            WRITE(*,*)' WAITALL:'
            CALL PSpew(Mssg,Reqs)
            WRITE(*,*)' ...aiting waiting waiting waiti... '
         ENDIF
!        Check to see if all requests are NULL
         N=SIZE(Reqs%I)
         DO I=1,N
            IF(Reqs%I(I)/=MPI_REQUEST_NULL)  &
               GOTO 10 ! Old habits die hard ...
         ENDDO
         RETURN
      10 CONTINUE 
!        Wait ...
         CALL New(Stats,(/MPI_STATUS_SIZE,N/))
         CALL MPI_WAITALL(N,Reqs%I(1:N),Stats%I,IErr)
         CALL ErrChk(IErr,Sub)
         CALL Delete(Stats)
      END SUBROUTINE WaitAll
!---------------------------------------------------------------
!
!
      SUBROUTINE WaitSome(Reqs,ToDo)
         TYPE(INT_VECT),             INTENT(INOUT)  :: Reqs,ToDo            
         TYPE(INT_VECT)                             :: Indxs
         TYPE(INT_RNK2)                             :: Stats
         INTEGER                                    :: I,J,L,N,IErr
         CHARACTER(LEN=8),PARAMETER                 :: Sub='WaitSome'
!-----------------------------------------------------------------------------
!        Allocation check
         IF(AllocQ(ToDo%Alloc))CALL Delete(ToDo)
!        Check to see if all requests are NULL
         N=SIZE(Reqs%I)
         DO I=1,N
            IF(Reqs%I(I)/=MPI_REQUEST_NULL)  &
               GOTO 10 ! Old habits die hard ...
         ENDDO
         CALL New(ToDo,1)
         ToDo%I=FAIL
         RETURN
      10 CONTINUE 
!        Wait ...
         CALL New(Indxs,N)
         CALL New(Stats,(/MPI_STATUS_SIZE,N/))
         CALL MPI_WAITSOME(N,Reqs%I,L,Indxs%I,Stats%I,IErr)         
         CALL ErrChk(IErr,Sub)
         CALL New(ToDo,L)
         ToDo%I(1:L)=Indxs%I(1:L)
         CALL Delete(Indxs)
         CALL Delete(Stats)
      END SUBROUTINE WaitSome
!===============================================================
!
!     WRAPER FOR MPI_TYPE_FREE
!
      SUBROUTINE FreeType(Type)
         INTEGER, INTENT(INOUT) :: Type
         INTEGER                :: IErr
         CHARACTER(LEN=8)       :: Sub='FreeType'
         IF(Type/=MPI_DATATYPE_NULL)THEN
            CALL MPI_TYPE_FREE(Type,IErr)
            CALL ErrChk(IErr,Sub)
         ENDIF
      END SUBROUTINE FreeType
!===============================================================
!
!     NODE ALLINGMENT
!
      SUBROUTINE AlignNodes(String_O)
         INTEGER                        :: IErr
         CHARACTER(LEN=*),OPTIONAL      :: String_O
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: String 
         CHARACTER(LEN=10),PARAMETER :: Sub='AlignNodes'
         CALL MPI_BARRIER(MONDO_COMM,IErr)
         CALL ErrChk(IErr,Sub)
         IF(PRESENT(String_O))THEN
            String="Node#"//TRIM(IntToChar(MyId))//' :: '//TRIM(String_O)
            WRITE(*,*)TRIM(String)
         ENDIF
      END SUBROUTINE AlignNodes
!===============================================================
!
!     Error checking
!
      SUBROUTINE ErrChk(IErr,Caller,Seqnce_O)
         INTEGER                              :: IErr
         CHARACTER(LEN=*),INTENT(IN)          :: Caller
         CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Seqnce_O
         IF(IErr/=MPI_SUCCESS)THEN
            IF(PRESENT(Seqnce_O))THEN
               CALL HaltMPI('>>>MondoMPI error:  '  & 
                         //TRIM(Seqnce_O)//'$'    &
                         //TRIM(Caller)//'.',IErr)

            ELSE
               CALL HaltMPI('>>>MondoMPI error:  '  & 
                         //TRIM(Caller)//'.',IErr)

            ENDIF
         ENDIF
      END SUBROUTINE


      SUBROUTINE InTurn(Mssg)
         CHARACTER(LEN=*) :: Mssg
         INTEGER          :: I,IErr
         CHARACTER(LEN=6),PARAMETER :: Sub='InTurn'
         CALL MPI_BARRIER(MONDO_COMM,IErr)
         CALL ErrChk(IErr,Sub)
         DO I=0,NPrc-1
            CALL MPI_BARRIER(MONDO_COMM,IErr)
            CALL ErrChk(IErr,Sub)
            IF(MyId==I)WRITE(*,*)' MyId = ',MyId,'$',TRIM(Mssg)
            CALL MPI_BARRIER(MONDO_COMM,IErr)
            CALL ErrChk(IErr,Sub)
         ENDDO  
      END SUBROUTINE InTurn


      SUBROUTINE PSpew_INT_VECT(Name,A)
         CHARACTER(LEN=*) :: Name
         TYPE(INT_VECT)   :: A
         INTEGER          :: I
         DO I=0,NPrc-1
            CALL AlignNodes()         
            IF(MyId==I)WRITE(*,55)MyId,TRIM(Name),A%I
         ENDDO
      55 FORMAT(' MyId = ',I2,', ',A,' = ',1000(I20,', '))
      ENDSUBROUTINE PSpew_INT_VECT

      SUBROUTINE PSpew_DBL_VECT(Name,A)
         CHARACTER(LEN=*) :: Name
         TYPE(DBL_VECT)   :: A
         INTEGER          :: I
         DO I=0,NPrc-1
            CALL AlignNodes()         
            IF(MyId==I)WRITE(*,55)MyId,TRIM(Name),A%D
         ENDDO
      55 FORMAT(' MyId = ',I2,', ',A,' = ',1000(D9.3,', '))
      ENDSUBROUTINE PSpew_DBL_VECT
#endif
END MODULE
