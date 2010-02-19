! vim: tw=0:sw=3
!------------------------------------------------------------------------------
!    This code is part of the FreeON suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the FreeON group for review, and possible
!    dissemination in future releases.
!------------------------------------------------------------------------------

#include "MondoConfig.h"

MODULE MondoMPI
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalObjects
   USE ProcessControl
   USE MemMan
   USE Parse
   USE MondoLogger
   USE Utilities

#if defined(PARALLEL) || defined(PARALLEL_CLONES)
   USE MPI
#endif

   IMPLICIT NONE

#if defined(PARALLEL) || defined(PARALLEL_CLONES)
   ! NOTE DIFFERENCES BETWEEN DERIVED TYPES, EG. INT_VECT, WHICH INCLUDE
   ! AUXILIARY INFO SUCH AS ARRAY BOUNDS, AND THE ARRAYS THEMSELVES GIVEN BY THE
   ! FULLER NAME, EG. INT_VECTOR, WHICH SPECIFIES THE ARRAY ONLY.
   INTEGER,SAVE :: MONDO_COMM

   INTERFACE BCast ! Wrappers for MPI_BCAST
      MODULE PROCEDURE BCast_DBL_SCLR, BCast_DBL_VECT,  BCast_DBL_RNK2, &
                       BCast_DBL_RNK3, BCast_DBL_RNK4,  BCast_DBL_RNK6, &
                       BCast_INT_SCLR,                                  &
                       BCast_INT_VECT, BCast_INT_RNK2,  BCast_INT_RNK3, &
                       BCast_INT_RNK4, BCast_CHR_SCLR,  BCast_LOG_SCLR, &
                       BCast_DEBG, Bcast_LOG_VECT
   END INTERFACE

   INTERFACE AllReduce ! Wrappers for MPI_ALLREDUCE
      MODULE PROCEDURE AllReduce_DBL_SCLR,AllReduce_INT_SCLR
   END INTERFACE

   INTERFACE Reduce ! Wrappers for MPI_REDUCE
      MODULE PROCEDURE Reduce_DBL_SCLR, Reduce_INT_SCLR
   END INTERFACE

   INTERFACE Gather ! Wrappers for MPI_GATHER and MPI GATHERV
      MODULE PROCEDURE Gather_INT_SCLR, Gather_INT_VECT, &
                       Gather_DBL_SCLR, Gather_DBL_VECT
   END INTERFACE

   INTERFACE Send ! Wrappers for MPI_SEND (Blocking)
      MODULE PROCEDURE Send_INT_SCLR, Send_LOG_SCLR,   &
                       Send_INT_VECT, Send_DBL_SCLR,   &
                       Send_DBL_VECT, Send_INT_VECTOR, &
                       Send_DBL_RNK2, Send_HGRho,      &
                       Send_CMPoles, Send_CRDS,        &
                       Send_CHR10_VECT, Send_PBCInfo,  &
                       Send_CellSet
   END INTERFACE

   INTERFACE Recv ! Wrappers for MPI_RECV (Blocking)
      MODULE PROCEDURE Recv_INT_SCLR, Recv_LOG_SCLR,   &
                       Recv_INT_VECT, Recv_DBL_SCLR,   &
                       Recv_DBL_VECT, Recv_INT_VECTOR, &
                       Recv_DBL_RNK2, Recv_HGRho,      &
                       Recv_CMPoles, Recv_CRDS,        &
                       Recv_CHR10_VECT, Recv_PBCInfo,  &
                       Recv_CellSet
   END INTERFACE

   INTERFACE ISend ! Wrappers for MPI_ISEND (Non-Blocking)
      MODULE PROCEDURE ISend_INT_SCLR, ISend_INT_VECT, &
                       ISend_DBL_SCLR, ISend_DBL_VECT
   END INTERFACE

   INTERFACE IRecv ! Wrappers for MPI_IRECV (Non-Blocking)
      MODULE PROCEDURE IRecv_INT_SCLR, IRecv_INT_VECT, &
                       IRecv_DBL_SCLR, IRecv_DBL_VECT
   END INTERFACE

   INTERFACE PSpew
      MODULE PROCEDURE PSpew_INT_VECT, PSpew_DBL_VECT
   END INTERFACE

CONTAINS

   !===============================================================================
   !     WRAPPERS FOR MPI INITIALIZATION AND FINALIZATION
   !===============================================================================
   SUBROUTINE InitMPI()
      INTEGER :: IErr
      CHARACTER(LEN=*),PARAMETER :: Sub='InitMPI'

      CALL MPI_INIT(IErr)
      CALL ErrChk(IErr,Sub)
      MONDO_COMM=MPI_COMM_WORLD

      ! Load global MPI variables
      MyID=MRank()
      NPrc=MSize()
      InParallel=.TRUE.
   END SUBROUTINE InitMPI

   SUBROUTINE FiniMPI()
      INTEGER                    :: IErr
      CHARACTER(LEN=*),PARAMETER :: Sub='FiniMPI'
#ifdef MPI2
      INTEGER                    :: PARENT
      !INTEGER                    :: INTRA_PARENT
      CHARACTER                  :: message_buffer

      ! Get the parent communicator
      CALL MPI_COMM_GET_PARENT(PARENT, IErr)

      ! Merge the parent and current local communicators
      !CALL MPI_INTERCOMM_MERGE(PARENT, .TRUE., INTRA_PARENT, IErr)

      ! Sychronize parent and child. We do that by sending a message to the
      ! parent.
      CALL MondoLog(DEBUG_MAXIMUM, "FiniMPI", "sending message to rank 0", "rank "//TRIM(IntToChar(MyID)))
      CALL MPI_SEND(message_buffer, 1, MPI_CHARACTER, 0, BARRIER_TAG, PARENT, IErr)

      ! Free the communicator.
      CALL MPI_COMM_FREE(PARENT, IErr)
      !CALL MPI_COMM_FREE(INTRA_PARENT, IErr)

      ! Wait before finalizing...
      CALL MPI_Barrier(MPI_COMM_WORLD, IErr)
#endif
      CALL MPI_FINALIZE(IErr)
   END SUBROUTINE FiniMPI

   FUNCTION MRank(communicator_O)
      INTEGER                     :: IErr,MRank
      INTEGER, OPTIONAL           :: communicator_O
      CHARACTER(LEN=*), PARAMETER :: Sub = 'MRank'

      IF(PRESENT(communicator_O)) THEN
         CALL MPI_COMM_RANK(communicator_O, MRank, IErr)
      ELSE
         CALL MPI_COMM_RANK(MONDO_COMM, MRank, IErr)
      ENDIF
      CALL ErrChk(IErr, Sub)
   END FUNCTION MRank

   FUNCTION MSize(communicator_O)
      INTEGER                     :: IErr,MSize
      INTEGER, OPTIONAL           :: communicator_O
      CHARACTER(LEN=*), PARAMETER :: Sub = 'MSize'

      IF(PRESENT(communicator_O)) THEN
         CALL MPI_COMM_SIZE(communicator_O, MSize, IErr)
      ELSE
         CALL MPI_COMM_SIZE(MONDO_COMM, MSize, IErr)
      ENDIF
      CALL ErrChk(IErr, Sub)
   END FUNCTION MSize

   !===============================================================================
   !     BCAST WRAPPERS
   !===============================================================================
   SUBROUTINE Bcast_CHR_SCLR(A)
      CHARACTER(LEN=*)                   :: A
      CHARACTER(LEN=1),DIMENSION(LEN(A)) :: LocalChars
      INTEGER                            :: IErr,I,L

      L=LEN(A)
      DO I=1,L
         LocalChars(I)=A(I:I)
      ENDDO
      CALL MPI_BCAST(LocalChars,L,MPI_CHARACTER,ROOT,MONDO_COMM,IErr)
      IF(IErr/=MPI_SUCCESS) THEN
         CALL HaltMPI('MPI_BCAST failed in Bcast_STR',IErr)
      ENDIF
      DO I=1,L
         A(I:I)=LocalChars(I)
      ENDDO
   END SUBROUTINE Bcast_CHR_SCLR

   SUBROUTINE Bcast_CHR_VECT(A)
      TYPE(CHR_VECT), INTENT(INOUT)  :: A
      INTEGER                        :: I,N

      N=SIZE(A%C)
      DO I=1,N
         CALL Bcast_CHR_SCLR(A%C(I))
      ENDDO
   END SUBROUTINE Bcast_CHR_VECT

   SUBROUTINE Bcast_LOG_SCLR(A)
      LOGICAL,INTENT(INOUT)      :: A
      INTEGER                    :: IErr,I

      I=0; IF(A)I=1
      CALL MPI_BCAST(I,1,MPI_INTEGER,ROOT,MONDO_COMM,IErr)
      IF(IErr/=MPI_SUCCESS) THEN
         CALL HaltMPI('MPI_BCAST failed in Bcast_STR',IErr)
      ENDIF
      A=.FALSE.; IF(I==1)A=.TRUE.
   END SUBROUTINE Bcast_LOG_SCLR

   SUBROUTINE Bcast_LOG_VECT(A)
      TYPE(LOG_VECT),INTENT(INOUT)            :: A
      TYPE(INT_VECT)                          :: IA
      INTEGER                                 :: IErr,I,N

      N=SIZE(A%L)
      CALL New(IA,N)
      IA%I=0
      DO I=1,N
         IF(A%L(I)) IA%I(I)=1
      ENDDO
      CALL MPI_BCAST(IA%I,N,MPI_INTEGER,ROOT,MONDO_COMM,IErr)
      IF(IErr/=MPI_SUCCESS) THEN
         CALL HaltMPI('MPI_BCAST failed in Bcast_LOG_VECT',IErr)
      ENDIF
      A%L=.FALSE.
      DO I=1,N
         IF(IA%I(I)==1) A%L(I)=.TRUE.
      ENDDO
      CALL Delete(IA)
   END SUBROUTINE Bcast_LOG_VECT

   SUBROUTINE Bcast_DBL_SCLR(D,Id_0)
      REAL(DOUBLE), INTENT(INOUT) :: D
      INTEGER, OPTIONAL           :: Id_0
      INTEGER                     :: IErr,Id

      ID=ROOT
      IF(PRESENT(Id_0))Id=Id_0
      CALL MPI_BCAST(D,1,MPI_DOUBLE_PRECISION,Id,MONDO_COMM,IErr)
      IF(IErr/=MPI_SUCCESS) THEN
         CALL HaltMPI('MPI_BCAST failed in Bcast_DBL',IErr)
      ENDIF
   END SUBROUTINE Bcast_DBL_SCLR

   SUBROUTINE Bcast_DBL_VECT(A,N_O)
      TYPE(DBL_VECT),INTENT(INOUT) :: A
      INTEGER,OPTIONAL,INTENT(IN)  :: N_O
      INTEGER                      :: IErr,L

      L=SIZE(A%D); IF(PRESENT(N_O)) L=N_O
      CALL MPI_BCAST(A%D,L,MPI_DOUBLE_PRECISION,ROOT,MONDO_COMM,IErr)
      IF(IErr/=MPI_SUCCESS) THEN
         CALL HaltMPI('MPI_BCAST failed in Bcast_DBL_VECT',IErr)
      ENDIF
   END SUBROUTINE Bcast_DBL_VECT

   SUBROUTINE Bcast_DBL_RNK2(A)
      TYPE(DBL_RNK2), INTENT(INOUT) :: A
      INTEGER :: L,M

      L=SIZE(A%D,1); M=SIZE(A%D,2)
      CALL CRNDAV(L*M,A%D)
   END SUBROUTINE Bcast_DBL_RNK2

   SUBROUTINE Bcast_DBL_RNK3(A)
      TYPE(DBL_RNK3), INTENT(INOUT) :: A
      INTEGER :: L,M,N

      L=SIZE(A%D,1); M=SIZE(A%D,2); N=SIZE(A%D,3)
      CALL CRNDAV(L*M*N,A%D)
   END SUBROUTINE Bcast_DBL_RNK3

   SUBROUTINE Bcast_DBL_RNK4(A)
      TYPE(DBL_RNK4), INTENT(INOUT) :: A
      INTEGER :: L,M,N,K

      L=SIZE(A%D,1); M=SIZE(A%D,2)
      N=SIZE(A%D,3); K=SIZE(A%D,4)
      CALL CRNDAV(L*M*N*K,A%D)
   END SUBROUTINE Bcast_DBL_RNK4

   SUBROUTINE Bcast_DBL_RNK6(A)
      TYPE(DBL_RNK6), INTENT(INOUT) :: A
      INTEGER :: L

      L=SIZE(A%D,1)*SIZE(A%D,2)*SIZE(A%D,3)*SIZE(A%D,4)*SIZE(A%D,5)*SIZE(A%D,6)
      CALL CRNDAV(L,A%D)
   END SUBROUTINE Bcast_DBL_RNK6

   SUBROUTINE Bcast_INT_SCLR(A,Id_O)
      INTEGER, INTENT(INOUT) :: A
         INTEGER, OPTIONAL      :: Id_O
         INTEGER                :: IErr,Id

         Id=ROOT; IF(PRESENT(Id_O))Id=Id_O
         CALL MPI_BCAST(A,1,MPI_INTEGER,Id,MONDO_COMM,IErr)
         IF(IErr/=MPI_SUCCESS) THEN
            CALL HaltMPI('MPI_BCAST failed in Bcast_INT',IErr)
         ENDIF
   END SUBROUTINE Bcast_INT_SCLR

   SUBROUTINE Bcast_INT_VECT(A,N_O)
      TYPE(INT_VECT),    INTENT(INOUT) :: A
      INTEGER, OPTIONAL, INTENT(IN)    :: N_O
      INTEGER                          :: IErr,N

      N=SIZE(A%I); IF(PRESENT(N_O))N=N_O
      CALL MPI_BCAST(A%I,N,MPI_INTEGER,ROOT,MONDO_COMM,IErr)
      IF(IErr/=MPI_SUCCESS) THEN
         CALL HaltMPI('MPI_BCAST failed in Bcast_INT_VECT',IErr)
      ENDIF
   END SUBROUTINE Bcast_INT_VECT

   SUBROUTINE Bcast_INT_RNK2(A)
      TYPE(INT_RNK2), INTENT(INOUT) :: A
      INTEGER :: L,M

      L=SIZE(A%I,1) ; M=SIZE(A%I,2)
      CALL CRNIAV(L*M,A%I)
   END SUBROUTINE Bcast_INT_RNK2

   SUBROUTINE Bcast_INT_RNK3(A)
      TYPE(INT_RNK3), INTENT(INOUT) :: A
      INTEGER :: L,M,N

      L=SIZE(A%I,1) ; M=SIZE(A%I,2); N=SIZE(A%I,3)
      CALL CRNIAV(L*M*N,A%I)
   END SUBROUTINE Bcast_INT_RNK3

   SUBROUTINE Bcast_INT_RNK4(A)
      TYPE(INT_RNK4),INTENT(INOUT) :: A
      INTEGER :: I,J,L,M

      I=SIZE(A%I,1); J=SIZE(A%I,2)
      L=SIZE(A%I,3) ; M=SIZE(A%I,4)
      CALL CRNIAV(I*J*L*M,A%I)
   END SUBROUTINE Bcast_INT_RNK4

   !===============================================================================
   !     RANK N -> RANK 1 BCAST
   !===============================================================================
   SUBROUTINE CRNIAV(L,A) ! Bcast RANK N Integer As VECTOR
      INTEGER :: L,IErr
      INTEGER, DIMENSION(L), INTENT(INOUT)  :: A

      CALL MPI_BCAST(A,L,MPI_INTEGER,ROOT,MONDO_COMM,IErr)
      IF(IErr/=MPI_SUCCESS) THEN
         CALL HaltMPI('MPI_BCAST failed for integer vector',IErr)
      ENDIF
   END SUBROUTINE CRNIAV

   SUBROUTINE CRNDAV(L,A) ! Bcast RANK N Double As VECTOR
      INTEGER :: L,IErr
      REAL(DOUBLE), DIMENSION(L), INTENT(INOUT)  :: A

      CALL MPI_BCAST(A,L,MPI_DOUBLE_PRECISION,ROOT,MONDO_COMM,IErr)
      IF(IErr/=MPI_SUCCESS) THEN
         CALL HaltMPI('MPI_BCAST failed for double vector',IErr)
      ENDIF
   END SUBROUTINE CRNDAV

   SUBROUTINE Bcast_DEBG(A)
      TYPE(DEBG) :: A

      CALL Bcast(A%Key)
      CALL Bcast(A%Chk)
      CALL Bcast(A%Mat)
      CALL Bcast(A%Set)
      CALL Bcast(A%Int)
      CALL Bcast(A%Rho)
      CALL Bcast(A%Fmt)
   END SUBROUTINE Bcast_DEBG

   !===============================================================================
   !     REDUCE WRAPPERS
   !===============================================================================
   FUNCTION AllReduce_INT_SCLR(D,Op_O)
      INTEGER, INTENT(INOUT)       :: D
      INTEGER,OPTIONAL,INTENT(IN)  :: Op_O
      INTEGER                      :: AllReduce_INT_SCLR
      INTEGER                      :: Op,IErr
      CHARACTER(LEN=*), PARAMETER  :: Sub='AllReduce_INT_SCLR'

      IF(PRESENT(Op_O))THEN
         Op=Op_O
      ELSE
         Op=MPI_SUM
      ENDIF
      CALL MPI_ALLREDUCE(D,AllReduce_INT_SCLR,1,MPI_INTEGER,Op,MONDO_COMM,IErr)
      CALL ErrChk(IErr,Sub)
   END FUNCTION AllReduce_INT_SCLR

   FUNCTION AllReduce_DBL_SCLR(D,Op_O)
      REAL(DOUBLE), INTENT(INOUT)  :: D
      INTEGER,OPTIONAL,INTENT(IN)  :: Op_O
      REAL(DOUBLE)                 :: AllReduce_DBL_SCLR
      INTEGER                      :: Op,IErr
      CHARACTER(LEN=*), PARAMETER  :: Sub='AllReduce_DBL_SCLR'

      IF(PRESENT(Op_O))THEN
         Op=Op_O
      ELSE
         Op=MPI_SUM
      ENDIF
      CALL MPI_ALLREDUCE(D,AllReduce_DBL_SCLR,1,MPI_DOUBLE_PRECISION,Op,MONDO_COMM,IErr)
      CALL ErrChk(IErr,Sub)
   END FUNCTION AllReduce_DBL_SCLR

   FUNCTION Reduce_DBL_SCLR(D,Op_O)
      REAL(DOUBLE), INTENT(INOUT)  :: D
      INTEGER,OPTIONAL,INTENT(IN)  :: Op_O
      REAL(DOUBLE)                 :: Reduce_DBL_SCLR
      INTEGER                      :: Op,IErr
      CHARACTER(LEN=*), PARAMETER  :: Sub='Reduce_DBL_SCLR'

      IF(PRESENT(Op_O))THEN
         Op=Op_O
      ELSE
         Op=MPI_SUM
      ENDIF
      CALL MPI_REDUCE(D,Reduce_DBL_SCLR,1,MPI_DOUBLE_PRECISION,Op,ROOT,MONDO_COMM,IErr)
      CALL ErrChk(IErr,Sub)
   END FUNCTION Reduce_DBL_SCLR

   FUNCTION Reduce_INT_SCLR(I)
      INTEGER                      :: I,Reduce_INT_SCLR
      INTEGER                      :: IErr
      CHARACTER(LEN=*), PARAMETER  :: Sub='Reduce_INT_SCLR'

      CALL MPI_REDUCE(I,Reduce_INT_SCLR,1,MPI_INTEGER,MPI_SUM,ROOT,MONDO_COMM,IErr)
      CALL ErrChk(IErr,Sub)
   END FUNCTION Reduce_INT_SCLR

   !===============================================================================
   !     GATHERS (Note, be carefull to avoid
   !     reference to uninitialised POINTERs
   !===============================================================================
   SUBROUTINE Gather_INT_SCLR(A,B)
      INTEGER,       INTENT(IN)    :: A
      TYPE(INT_VECT),INTENT(INOUT) :: B
      INTEGER                      :: IErr
      CHARACTER(LEN=*), PARAMETER  :: Sub='Gather_INT_SCLR'

      IF(MyId==ROOT)THEN
         CALL MPI_GATHER(A,1,MPI_INTEGER,B%I,1,MPI_INTEGER,ROOT,MONDO_COMM,IErr)
      ELSE
         ! Non-root nodes may not have B allocated, need to protect against
         ! refrencing a pointer
         CALL MPI_GATHER(A,1,MPI_INTEGER,1,1,MPI_INTEGER,ROOT,MONDO_COMM,IErr)
      ENDIF
      CALL ErrChk(IErr,Sub)
   END SUBROUTINE Gather_INT_SCLR

   SUBROUTINE Gather_INT_VECT(Send,Recv,NSnd,NRcv,Disp)
      TYPE(INT_VECT),INTENT(IN)    :: Send,NRcv,Disp
      INTEGER,       INTENT(IN)    :: NSnd
      TYPE(INT_VECT),INTENT(INOUT) :: Recv
      INTEGER                      :: IErr
      CHARACTER(LEN=*), PARAMETER  :: Sub='Gather_INT_VECT'

      IF(MyId==ROOT)THEN
         CALL MPI_GATHERV(Send%I,NSnd,MPI_INTEGER, Recv%I,NRcv%I,Disp%I,MPI_INTEGER,ROOT,MONDO_COMM,IErr)
      ELSE
         CALL MPI_GATHERV(Send%I,NSnd,MPI_INTEGER,1,1,1,MPI_INTEGER,ROOT,MONDO_COMM,IErr)
      ENDIF
      CALL ErrChk(IErr,Sub)
   END SUBROUTINE Gather_INT_VECT

   SUBROUTINE Gather_DBL_SCLR(A,B)
      REAL(DOUBLE),  INTENT(IN)    :: A
      TYPE(DBL_VECT),INTENT(INOUT) :: B
      INTEGER                      :: IErr
      CHARACTER(LEN=*),PARAMETER   :: Sub='Gather_DBL_SCLR'

      IF(MyId==ROOT)THEN
         CALL MPI_GATHER(A,1,MPI_DOUBLE_PRECISION,B%D,1,MPI_DOUBLE_PRECISION,ROOT,MONDO_COMM,IErr)
      ELSE
         CALL MPI_GATHER(A,1,MPI_DOUBLE_PRECISION,1.0D0,1,MPI_DOUBLE_PRECISION,ROOT,MONDO_COMM,IErr)
      ENDIF
      CALL ErrChk(IErr,Sub)
   END SUBROUTINE Gather_DBL_SCLR

   SUBROUTINE Gather_DBL_VECT(Send,Recv,NSnd,NRcv,Disp)
      TYPE(DBL_VECT),INTENT(IN)    :: Send
      TYPE(INT_VECT),INTENT(IN)    :: NRcv,Disp
      INTEGER,       INTENT(IN)    :: NSnd
      TYPE(DBl_VECT),INTENT(INOUT) :: Recv
      INTEGER                      :: IErr
      CHARACTER(LEN=*),PARAMETER   :: Sub='Gather_DBL_VECT'

      IF(MyId==ROOT)THEN
         CALL MPI_GATHERV(Send%D,NSnd,MPI_DOUBLE_PRECISION,Recv%D,NRcv%I,Disp%I,MPI_DOUBLE_PRECISION,ROOT,MONDO_COMM,IErr)
      ELSE
         CALL MPI_GATHERV(Send%D,NSnd,MPI_DOUBLE_PRECISION,1.0D0,1,1,MPI_DOUBLE_PRECISION,ROOT,MONDO_COMM,IErr)
      ENDIF
      CALL ErrChk(IErr,Sub)
   END SUBROUTINE Gather_DBL_VECT

   !===============================================================================
   !     BLOCKING SENDS AND RECIEVES
   !===============================================================================
   SUBROUTINE Send_DBL_SCLR(Snd,To,Tag, comm_O)
      REAL(DOUBLE),INTENT(IN)       :: Snd
      INTEGER,     INTENT(IN)       :: To,Tag
      INTEGER, OPTIONAL, INTENT(IN) :: comm_O
      INTEGER                       :: IErr
      INTEGER                       :: communicator
      CHARACTER(LEN=*), PARAMETER   :: Sub='Send_DBL_SCLR'

      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      CALL MPI_SEND(Snd,1,MPI_DOUBLE_PRECISION,To,Tag,MPI_STATUS_IGNORE,IErr)
      CALL ErrChk(IErr,Sub)
   END SUBROUTINE Send_DBL_SCLR

   SUBROUTINE Send_DBL_VECT(Snd,N,To,Tag,M_O, comm_O)
      TYPE(DBL_VECT),  INTENT(IN)   :: Snd
      INTEGER,         INTENT(IN)   :: N,To,Tag
      INTEGER, OPTIONAL, INTENT(IN) :: comm_O
      INTEGER                       :: L,M,IErr
      INTEGER,OPTIONAL,INTENT(IN)   :: M_O
      INTEGER                       :: communicator
      CHARACTER(LEN=*),PARAMETER    :: Sub='Send_DBL_SCLR'

      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
      CALL MPI_SEND(Snd%D(M:N),L,MPI_DOUBLE_PRECISION,To,Tag,communicator,IErr)
      CALL ErrChk(IErr,Sub)
   END SUBROUTINE Send_DBL_VECT

   SUBROUTINE Send_INT_SCLR(Snd, To, Tag, comm_O)
      INTEGER, INTENT(IN)           :: Snd,To,Tag
      INTEGER, OPTIONAL, INTENT(IN) :: comm_O
      INTEGER                       :: IErr
      INTEGER                       :: communicator
      CHARACTER(LEN=*), PARAMETER   :: Sub='Send_INT_SCLR'

      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      CALL MPI_SEND(Snd, 1, MPI_INTEGER, To, Tag, communicator, IErr)
      CALL ErrChk(IErr, Sub)
   END SUBROUTINE Send_INT_SCLR

   SUBROUTINE Send_LOG_SCLR(var, rank, tag, comm_O)
      LOGICAL, INTENT(IN)           :: var
      INTEGER, INTENT(IN)           :: rank, tag
      INTEGER, OPTIONAL, INTENT(IN) :: comm_O
      INTEGER                       :: IErr
      INTEGER                       :: communicator
      CHARACTER(LEN=*), PARAMETER   :: Sub='Send_LOG_SCLR'

      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      CALL MPI_SEND(var, 1, MPI_LOGICAL, rank, tag, communicator, IErr)
      CALL ErrChk(IErr, Sub)
   END SUBROUTINE Send_LOG_SCLR

   SUBROUTINE Send_INT_VECT(Snd,N,To,Tag,M_O, comm_O)
      TYPE(INT_VECT),  INTENT(IN)     :: Snd
      INTEGER,         INTENT(IN)     :: N,To,Tag
      INTEGER, OPTIONAL, INTENT(IN)   :: comm_O
      INTEGER                         :: L,M,IErr
      INTEGER,OPTIONAL,INTENT(IN)     :: M_O
      INTEGER                         :: communicator
      CHARACTER(LEN=*),PARAMETER      :: Sub='Send_INT_VECT'

      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
      CALL MPI_SEND(Snd%I(M:N),L,MPI_INTEGER,To,Tag,communicator,IErr)
      CALL ErrChk(IErr,Sub)
   END SUBROUTINE Send_INT_VECT

   SUBROUTINE Send_INT_VECTOR(Snd,N,To,Tag,M_O, comm_O)
      INTEGER,DIMENSION(:),INTENT(IN) :: Snd
      INTEGER,             INTENT(IN) :: N,To,Tag
      INTEGER, OPTIONAL, INTENT(IN)   :: comm_O
      INTEGER                         :: L,M,IErr
      INTEGER,OPTIONAL,INTENT(IN)     :: M_O
      INTEGER                         :: communicator
      CHARACTER(LEN=*),PARAMETER      :: Sub='Send_INT_VECTOR'

      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
      CALL MPI_SEND(Snd(M:N),L,MPI_INTEGER,To,Tag,communicator,IErr)
      CALL ErrChk(IErr,Sub)
   END SUBROUTINE Send_INT_VECTOR

   SUBROUTINE Send_DBL_RNK2(var, rank, tag, comm_O)
      TYPE(DBL_RNK2), INTENT(IN)    :: var
      INTEGER, INTENT(IN)           :: rank, tag
      INTEGER, OPTIONAL, INTENT(IN) :: comm_O
      INTEGER, DIMENSION(2)         :: NLower, NUpper
      INTEGER                       :: communicator
      INTEGER                       :: IErr
      CHARACTER(LEN=*),PARAMETER    :: Sub = 'Send_DBL_RNK2'

      IF(.NOT. AllocQ(var%Alloc)) THEN
         CALL Halt("["//TRIM(Sub)//"] variable not allocated")
      ENDIF
      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      NLower(1) = LBOUND(var%D, 1)
      NUpper(1) = UBOUND(var%D, 1)
      NLower(2) = LBOUND(var%D, 2)
      NUpper(2) = UBOUND(var%D, 2)
      CALL Send(NLower(1), rank, tag, comm_O = communicator)
      CALL Send(NUpper(1), rank, tag, comm_O = communicator)
      CALL Send(NLower(2), rank, tag, comm_O = communicator)
      CALL Send(NUpper(2), rank, tag, comm_O = communicator)
      CALL MPI_Send(var%D, (NUpper(1)-NLower(1)+1)*(NUpper(2)-NLower(2)+1), MPI_DOUBLE_PRECISION, rank, PUT_TAG, communicator, IErr)
      CALL ErrChk(IErr,Sub)
   END SUBROUTINE Send_DBL_RNK2

   SUBROUTINE Send_HGRho(var, rank, tag, comm_O)
      TYPE(HGRho), INTENT(IN)       :: var
      INTEGER, INTENT(IN)           :: rank, tag
      INTEGER, OPTIONAL, INTENT(IN) :: comm_O
      INTEGER                       :: communicator
      CHARACTER(LEN=*),PARAMETER    :: Sub = 'Send_HGRho'

      IF(.NOT. AllocQ(var%Alloc)) THEN
         CALL Halt("["//TRIM(Sub)//"] variable not allocated")
      ENDIF
      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      CALL Send(var%NSDen, rank, tag, comm_O = communicator)
      CALL Send(var%NExpt, rank, tag, comm_O = communicator)
      CALL Send(var%NDist, rank, tag, comm_O = communicator)
      CALL Send(var%NCoef, rank, tag, comm_O = communicator)
      CALL Send(var%NQ, UBOUND(var%NQ%I, 1), rank, tag, M_O = LBOUND(var%NQ%I, 1), comm_O = communicator)
      CALL Send(var%Lndx, UBOUND(var%Lndx%I, 1), rank, tag, M_O = LBOUND(var%Lndx%I, 1), comm_O = communicator)
      CALL Send(var%OffQ, UBOUND(var%OffQ%I, 1), rank, tag, M_O = LBOUND(var%OffQ%I, 1), comm_O = communicator)
      CALL Send(var%OffR, UBOUND(var%OffR%I, 1), rank, tag, M_O = LBOUND(var%OffR%I, 1), comm_O = communicator)
      CALL Send(var%Expt, UBOUND(var%Expt%D, 1), rank, tag, M_O = LBOUND(var%Expt%D, 1), comm_O = communicator)
      CALL Send(var%Qx, UBOUND(var%Qx%D, 1), rank, tag, M_O = LBOUND(var%Qx%D, 1), comm_O = communicator)
      CALL Send(var%Qy, UBOUND(var%Qy%D, 1), rank, tag, M_O = LBOUND(var%Qy%D, 1), comm_O = communicator)
      CALL Send(var%Qz, UBOUND(var%Qz%D, 1), rank, tag, M_O = LBOUND(var%Qz%D, 1), comm_O = communicator)
      CALL Send(var%Co, UBOUND(var%Co%D, 1), rank, tag, M_O = LBOUND(var%Co%D, 1), comm_O = communicator)
   END SUBROUTINE Send_HGRho

   SUBROUTINE Send_CMPoles(var, rank, tag, comm_O)
      TYPE(CMPoles), INTENT(IN)     :: var
      INTEGER, INTENT(IN)           :: rank, tag
      INTEGER, OPTIONAL, INTENT(IN) :: comm_O
      INTEGER                       :: communicator
      CHARACTER(LEN=*),PARAMETER    :: Sub = 'Send_CMPoles'

      IF(.NOT. AllocQ(var%Alloc)) THEN
         CALL Halt("["//TRIM(Sub)//"] variable not allocated")
      ENDIF
      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      CALL Send(var%MPole, rank, tag, comm_O = communicator)
      CALL Send(var%DPole, 3, rank, tag, comm_O = communicator)
      CALL Send(var%QPole, 6, rank, tag, comm_O = communicator)
   END SUBROUTINE Send_CMPoles

   SUBROUTINE Send_CRDS(var, rank, tag, comm_O)
      TYPE(CRDS), INTENT(IN)        :: var
      INTEGER, INTENT(IN)           :: rank, tag
      INTEGER, OPTIONAL, INTENT(IN) :: comm_O
      INTEGER                       :: communicator
      CHARACTER(LEN=*),PARAMETER    :: Sub = 'Send_CRDS'

      IF(.NOT. AllocQ(var%Alloc)) THEN
         CALL Halt("["//TRIM(Sub)//"] variable not allocated")
      ENDIF
      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      CALL Send(var%NAtms, rank, tag, comm_O = communicator)

      CALL Send(var%InAU, rank, tag, comm_O = communicator)
      CALL Send(var%Ordrd, rank, tag, comm_O = communicator)
      CALL Send(var%Confg, rank, tag, comm_O = communicator)
      CALL Send(var%NElec, rank, tag, comm_O = communicator)
      CALL Send(var%Multp, rank, tag, comm_O = communicator)
      CALL Send(var%TotCh, rank, tag, comm_O = communicator)
      CALL Send(var%NAlph, rank, tag, comm_O = communicator)
      CALL Send(var%NBeta, rank, tag, comm_O = communicator)
      CALL Send(var%ETotal, rank, tag, comm_O = communicator)
      CALL Send(var%ETotalPerSCF, UBOUND(var%ETotalPerSCF%D, 1), rank, tag, M_O = LBOUND(var%ETotalPerSCF%D, 1), comm_O = communicator)
      CALL Send(var%GradRMS, rank, tag, comm_O = communicator)
      CALL Send(var%GradMax, rank, tag, comm_O = communicator)
      CALL Send(var%Unstable, rank, tag, comm_O = communicator)
      CALL Send(var%BndBox, rank, tag, comm_O = communicator)
      CALL Send(var%PBC, rank, tag, comm_O = communicator)
      CALL Send(var%OvCells, rank, tag, comm_O = communicator)
      CALL Send(var%InCells, rank, tag, comm_O = communicator)
      CALL Send(var%NKind, rank, tag, comm_O = communicator)
      CALL Send(var%AtNum, UBOUND(var%AtNum%D, 1), rank, tag, M_O = LBOUND(var%AtNum%D, 1), comm_O = communicator)
      CALL Send(var%AtTyp, UBOUND(var%AtTyp%I, 1), rank, tag, M_O = LBOUND(var%AtTyp%I, 1), comm_O = communicator)
      CALL Send(var%AtNam, rank, tag, comm_O = communicator)
      CALL Send(var%AtMss, UBOUND(var%AtMss%D, 1), rank, tag, M_O = LBOUND(var%AtMss%D, 1), comm_O = communicator)
      CALL Send(var%CConstrain, UBOUND(var%CConstrain%I, 1), rank, tag, M_O = LBOUND(var%CConstrain%I, 1), comm_O = communicator)
      CALL Send(var%DoFreq, UBOUND(var%DoFreq%I, 1), rank, tag, M_O = LBOUND(var%DoFreq%I, 1), comm_O = communicator)
      CALL Send(var%Carts, rank, tag, comm_O = communicator)
      CALL Send(var%BoxCarts, rank, tag, comm_O = communicator)
      CALL Send(var%Velocity, rank, tag, comm_O = communicator)
      CALL Send(var%Gradients, rank, tag, comm_O = communicator)
      CALL Send(var%Fext, rank, tag, comm_O = communicator)
      CALL Send(var%Displ, rank, tag, comm_O = communicator)
      CALL Send(var%PBCDispl, rank, tag, comm_O = communicator)
      CALL Send(var%LatticeOnly, rank, tag, comm_O = communicator)
      CALL Send(var%AltCount, rank, tag, comm_O = communicator)
   END SUBROUTINE Send_CRDS

   SUBROUTINE Send_CHR10_VECT(var, rank, tag, comm_O)
      TYPE(CHR10_VECT), INTENT(IN)  :: var
      INTEGER, INTENT(IN)           :: rank, tag
      INTEGER, OPTIONAL, INTENT(IN) :: comm_O
      INTEGER                       :: communicator, i, IErr
      CHARACTER(LEN=*),PARAMETER    :: Sub = 'Send_CHR10_VECT'

      IF(.NOT. AllocQ(var%Alloc)) THEN
         CALL Halt("["//TRIM(Sub)//"] variable not allocated")
      ENDIF
      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      CALL Send(LBOUND(var%C, 1), rank, tag, comm_O = communicator)
      CALL Send(UBOUND(var%C, 1), rank, tag, comm_O = communicator)
      DO i = LBOUND(var%C, 1), UBOUND(var%C, 1)
         CALL MPI_Send(var%C(i), 10, MPI_CHARACTER, rank, tag, communicator, IErr)
         CALL ErrChk(IErr,Sub)
      ENDDO
   END SUBROUTINE Send_CHR10_VECT

   SUBROUTINE Send_PBCInfo(var, rank, tag, comm_O)
      TYPE(PBCInfo), INTENT(IN)     :: var
      INTEGER, INTENT(IN)           :: rank, tag
      INTEGER, OPTIONAL, INTENT(IN) :: comm_O
      INTEGER                       :: communicator, i, IErr
      CHARACTER(LEN=*),PARAMETER    :: Sub = 'Send_PBCInfo'

      IF(.NOT. AllocQ(var%Alloc)) THEN
         CALL Halt("["//TRIM(Sub)//"] variable not allocated")
      ENDIF
      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      CALL Send(var%Dimen, rank, tag, comm_O = communicator)
      CALL Send(var%PFFMaxEll, rank, tag, comm_O = communicator)
      CALL Send(var%PFFWelSep, rank, tag, comm_O = communicator)
      CALL Send(var%InVecForm, rank, tag, comm_O = communicator)
      CALL Send(var%InAtomCrd, rank, tag, comm_O = communicator)
      CALL Send(var%Translate, rank, tag, comm_O = communicator)
      CALL Send(var%CellVolume, rank, tag, comm_O = communicator)
      CALL Send(var%Epsilon, rank, tag, comm_O = communicator)
      CALL Send(var%DipoleFAC, rank, tag, comm_O = communicator)
      CALL Send(var%AutoW, 3, rank, tag, comm_O = communicator)
      CALL Send(var%SuperCell, 3, rank, tag, comm_O = communicator)
      CALL Send(var%CellCenter, 3, rank, tag, comm_O = communicator)
      CALL Send(var%TransVec, 3, rank, tag, comm_O = communicator)
      CALL Send(var%BoxShape, rank, tag, comm_O = communicator)
      CALL Send(var%InvBoxSh, rank, tag, comm_O = communicator)
      CALL Send(var%LatFrc, rank, tag, comm_O = communicator)
   END SUBROUTINE Send_PBCInfo

   SUBROUTINE Send_CellSet(var, rank, tag, comm_O)
      TYPE(CellSet), INTENT(IN)     :: var
      INTEGER, INTENT(IN)           :: rank, tag
      INTEGER, OPTIONAL, INTENT(IN) :: comm_O
      INTEGER                       :: communicator, i, IErr
      CHARACTER(LEN=*),PARAMETER    :: Sub = 'Send_CellSet'

      IF(.NOT. AllocQ(var%Alloc)) THEN
         CALL Halt("["//TRIM(Sub)//"] variable not allocated")
      ENDIF
      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      CALL Send(var%NCells, rank, tag, comm_O = communicator)
      CALL Send(var%Radius, rank, tag, comm_O = communicator)
      CALL Send(var%CellCarts, rank, tag, comm_O = communicator)
   END SUBROUTINE Send_CellSet

   SUBROUTINE Recv_DBL_SCLR(Rec,From,Tag, comm_O)
      REAL(DOUBLE),INTENT(INOUT)          :: Rec
      INTEGER,     INTENT(IN)             :: From,Tag
      INTEGER                             :: IErr
      INTEGER, OPTIONAL, INTENT(IN)       :: comm_O
      INTEGER                             :: communicator
      CHARACTER(LEN=*),PARAMETER          :: Sub='Recv_DBL_SCLR'

      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      CALL MPI_RECV(Rec,1,MPI_DOUBLE_PRECISION,From,Tag,communicator,MPI_STATUS_IGNORE,IErr)
      CALL ErrChk(IErr,Sub)
   END SUBROUTINE Recv_DBL_SCLR

   SUBROUTINE Recv_DBL_VECT(Rec,N,From,Tag,M_O, comm_O)
      TYPE(DBL_VECT),INTENT(INOUT)        :: Rec
      INTEGER,       INTENT(IN)           :: N,From,Tag
      INTEGER, OPTIONAL, INTENT(IN)       :: comm_O
      INTEGER                             :: L,M,IErr
      INTEGER, OPTIONAL, INTENT(IN)       :: M_O
      INTEGER                             :: communicator
      CHARACTER(LEN=*),PARAMETER          :: Sub='Recv_DBL_VECT'

      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
      CALL MPI_RECV(Rec%D(M:N),L,MPI_DOUBLE_PRECISION,From,Tag,communicator,MPI_STATUS_IGNORE,IErr)
      CALL ErrChk(IErr,Sub)
   END SUBROUTINE Recv_DBL_VECT

   SUBROUTINE Recv_INT_SCLR(Rec, From, Tag, comm_O)
      INTEGER, INTENT(INOUT)              :: Rec
      INTEGER, INTENT(IN)                 :: From,Tag
      INTEGER, OPTIONAL, INTENT(IN)       :: comm_O
      INTEGER                             :: IErr
      INTEGER                             :: communicator
      CHARACTER(LEN=*),PARAMETER          :: Sub='Recv_INT_SCLR'

      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      CALL MPI_RECV(Rec, 1, MPI_INTEGER, From, Tag, communicator, MPI_STATUS_IGNORE, IErr)
      CALL ErrChk(IErr, Sub)
   END SUBROUTINE Recv_INT_SCLR

   SUBROUTINE Recv_LOG_SCLR(var, rank, tag, comm_O)
      LOGICAL, INTENT(INOUT)              :: var
      INTEGER, INTENT(IN)                 :: rank, tag
      INTEGER, OPTIONAL, INTENT(IN)       :: comm_O
      INTEGER                             :: IErr
      INTEGER                             :: communicator
      CHARACTER(LEN=*),PARAMETER          :: Sub='Recv_LOG_SCLR'

      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      CALL MPI_RECV(var, 1, MPI_LOGICAL, rank, tag, communicator, MPI_STATUS_IGNORE, IErr)
      CALL ErrChk(IErr, Sub)
   END SUBROUTINE Recv_LOG_SCLR

   SUBROUTINE Recv_INT_VECT(Rec,N,From,Tag,M_O, comm_O)
      TYPE(INT_VECT),INTENT(INOUT)        :: Rec
      INTEGER,       INTENT(IN)           :: N,From,Tag
      INTEGER, OPTIONAL, INTENT(IN)       :: comm_O
      INTEGER                             :: L,M,IErr
      INTEGER, OPTIONAL, INTENT(IN)       :: M_O
      INTEGER                             :: communicator
      CHARACTER(LEN=*),PARAMETER          :: Sub='Recv_INT_VECT'

      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
      CALL MPI_RECV(Rec%I(M:N),L,MPI_INTEGER,From,Tag,communicator,MPI_STATUS_IGNORE,IErr)
      CALL ErrChk(IErr,Sub)
   END SUBROUTINE Recv_INT_VECT

   SUBROUTINE Recv_INT_VECTOR(Rec,N,From,Tag,M_O)
      INTEGER,DIMENSION(:),INTENT(INOUT)  :: Rec
      INTEGER,             INTENT(IN)     :: N,From,Tag
      INTEGER                             :: L,M,IErr
      INTEGER, OPTIONAL, INTENT(IN)       :: M_O
      CHARACTER(LEN=*),PARAMETER          :: Sub='Recv_INT_VECTOR'

      M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
      CALL MPI_RECV(Rec(M:N),L,MPI_INTEGER,From,Tag,MONDO_COMM,MPI_STATUS_IGNORE,IErr)
      CALL ErrChk(IErr,Sub)
   END SUBROUTINE Recv_INT_VECTOR

   SUBROUTINE Recv_DBL_RNK2(var, rank, tag, comm_O)
      TYPE(DBL_RNK2), INTENT(INOUT)       :: var
      INTEGER, INTENT(IN)                 :: rank, tag
      INTEGER, OPTIONAL, INTENT(IN)       :: comm_O
      INTEGER, DIMENSION(2)               :: NLower, NUpper
      INTEGER                             :: communicator
      INTEGER                             :: IErr
      CHARACTER(LEN=*),PARAMETER          :: Sub = 'Recv_DBL_RNK2'

      IF(AllocQ(var%Alloc)) THEN
         CALL Delete(var)
      ENDIF
      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      CALL Recv(NLower(1), rank, tag, comm_O = communicator)
      CALL Recv(NUpper(1), rank, tag, comm_O = communicator)
      CALL Recv(NLower(2), rank, tag, comm_O = communicator)
      CALL Recv(NUpper(2), rank, tag, comm_O = communicator)
      CALL New(var, NUpper, M_O = NLower)
      CALL MPI_Recv(var%D, (NUpper(1)-NLower(1)+1)*(NUpper(2)-NLower(2)+1), MPI_DOUBLE_PRECISION, rank, PUT_TAG, communicator, MPI_STATUS_IGNORE, IErr)
      CALL ErrChk(IErr,Sub)
   END SUBROUTINE Recv_DBL_RNK2

   SUBROUTINE Recv_HGRho(var, rank, tag, comm_O)
      TYPE(HGRho), INTENT(INOUT)    :: var
      INTEGER, INTENT(IN)           :: rank, tag
      INTEGER, OPTIONAL, INTENT(IN) :: comm_O
      INTEGER                       :: communicator
      INTEGER                       :: NSDen, NExpt, NDist, NCoef
      CHARACTER(LEN=*),PARAMETER    :: Sub = 'Send_HGRho'

      IF(AllocQ(var%Alloc)) THEN
         CALL Delete(var)
      ENDIF
      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      CALL Recv(NSDen, rank, tag, comm_O = communicator)
      CALL Recv(NExpt, rank, tag, comm_O = communicator)
      CALL Recv(NDist, rank, tag, comm_O = communicator)
      CALL Recv(NCoef, rank, tag, comm_O = communicator)

      ! Allocate density.
      CALL NEW(var, N_O = (/ NExpt, NDist, NCoef, NSDen /))

      CALL Recv(var%NQ, UBOUND(var%NQ%I, 1), rank, tag, M_O = LBOUND(var%NQ%I, 1), comm_O = communicator)
      CALL Recv(var%Lndx, UBOUND(var%Lndx%I, 1), rank, tag, M_O = LBOUND(var%Lndx%I, 1), comm_O = communicator)
      CALL Recv(var%OffQ, UBOUND(var%OffQ%I, 1), rank, tag, M_O = LBOUND(var%OffQ%I, 1), comm_O = communicator)
      CALL Recv(var%OffR, UBOUND(var%OffR%I, 1), rank, tag, M_O = LBOUND(var%OffR%I, 1), comm_O = communicator)
      CALL Recv(var%Expt, UBOUND(var%Expt%D, 1), rank, tag, M_O = LBOUND(var%Expt%D, 1), comm_O = communicator)
      CALL Recv(var%Qx, UBOUND(var%Qx%D, 1), rank, tag, M_O = LBOUND(var%Qx%D, 1), comm_O = communicator)
      CALL Recv(var%Qy, UBOUND(var%Qy%D, 1), rank, tag, M_O = LBOUND(var%Qy%D, 1), comm_O = communicator)
      CALL Recv(var%Qz, UBOUND(var%Qz%D, 1), rank, tag, M_O = LBOUND(var%Qz%D, 1), comm_O = communicator)
      CALL Recv(var%Co, UBOUND(var%Co%D, 1), rank, tag, M_O = LBOUND(var%Co%D, 1), comm_O = communicator)
   END SUBROUTINE Recv_HGRho

   SUBROUTINE Recv_CMPoles(var, rank, tag, comm_O)
      TYPE(CMPoles), INTENT(INOUT)  :: var
      INTEGER, INTENT(IN)           :: rank, tag
      INTEGER, OPTIONAL, INTENT(IN) :: comm_O
      INTEGER                       :: communicator
      CHARACTER(LEN=*),PARAMETER    :: Sub = 'Send_CMPoles'

      IF(AllocQ(var%Alloc)) THEN
         CALL Delete(var)
      ENDIF
      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      CALL New(var)
      CALL Recv(var%MPole, rank, tag, comm_O = communicator)
      CALL Recv(var%DPole, 3, rank, tag, comm_O = communicator)
      CALL Recv(var%QPole, 6, rank, tag, comm_O = communicator)
   END SUBROUTINE Recv_CMPoles

   SUBROUTINE Recv_CRDS(var, rank, tag, comm_O)
      TYPE(CRDS), INTENT(INOUT)     :: var
      INTEGER, INTENT(IN)           :: rank, tag
      INTEGER, OPTIONAL, INTENT(IN) :: comm_O
      INTEGER                       :: communicator
      CHARACTER(LEN=*),PARAMETER    :: Sub = 'Recv_CRDS'

      IF(AllocQ(var%Alloc)) THEN
         CALL Delete(var)
      ENDIF
      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      CALL Recv(var%NAtms, rank, tag, comm_O = communicator)
      CALL New(var)

      CALL Recv(var%InAU, rank, tag, comm_O = communicator)
      CALL Recv(var%Ordrd, rank, tag, comm_O = communicator)
      CALL Recv(var%Confg, rank, tag, comm_O = communicator)
      CALL Recv(var%NElec, rank, tag, comm_O = communicator)
      CALL Recv(var%Multp, rank, tag, comm_O = communicator)
      CALL Recv(var%TotCh, rank, tag, comm_O = communicator)
      CALL Recv(var%NAlph, rank, tag, comm_O = communicator)
      CALL Recv(var%NBeta, rank, tag, comm_O = communicator)
      CALL Recv(var%ETotal, rank, tag, comm_O = communicator)
      CALL Recv(var%ETotalPerSCF, UBOUND(var%ETotalPerSCF%D, 1), rank, tag, M_O = LBOUND(var%ETotalPerSCF%D, 1), comm_O = communicator)
      CALL Recv(var%GradRMS, rank, tag, comm_O = communicator)
      CALL Recv(var%GradMax, rank, tag, comm_O = communicator)
      CALL Recv(var%Unstable, rank, tag, comm_O = communicator)
      CALL Recv(var%BndBox, rank, tag, comm_O = communicator)
      CALL Recv(var%PBC, rank, tag, comm_O = communicator)
      CALL Recv(var%OvCells, rank, tag, comm_O = communicator)
      CALL Recv(var%InCells, rank, tag, comm_O = communicator)
      CALL Recv(var%NKind, rank, tag, comm_O = communicator)
      CALL Recv(var%AtNum, UBOUND(var%AtNum%D, 1), rank, tag, M_O = LBOUND(var%AtNum%D, 1), comm_O = communicator)
      CALL Recv(var%AtTyp, UBOUND(var%AtTyp%I, 1), rank, tag, M_O = LBOUND(var%AtTyp%I, 1), comm_O = communicator)
      CALL Recv(var%AtNam, rank, tag, comm_O = communicator)
      CALL Recv(var%AtMss, UBOUND(var%AtMss%D, 1), rank, tag, M_O = LBOUND(var%AtMss%D, 1), comm_O = communicator)
      CALL Recv(var%CConstrain, UBOUND(var%CConstrain%I, 1), rank, tag, M_O = LBOUND(var%CConstrain%I, 1), comm_O = communicator)
      CALL Recv(var%DoFreq, UBOUND(var%DoFreq%I, 1), rank, tag, M_O = LBOUND(var%DoFreq%I, 1), comm_O = communicator)
      CALL Recv(var%Carts, rank, tag, comm_O = communicator)
      CALL Recv(var%BoxCarts, rank, tag, comm_O = communicator)
      CALL Recv(var%Velocity, rank, tag, comm_O = communicator)
      CALL Recv(var%Gradients, rank, tag, comm_O = communicator)
      CALL Recv(var%Fext, rank, tag, comm_O = communicator)
      CALL Recv(var%Displ, rank, tag, comm_O = communicator)
      CALL Recv(var%PBCDispl, rank, tag, comm_O = communicator)
      CALL Recv(var%LatticeOnly, rank, tag, comm_O = communicator)
      CALL Recv(var%AltCount, rank, tag, comm_O = communicator)
   END SUBROUTINE Recv_CRDS

   SUBROUTINE Recv_CHR10_VECT(var, rank, tag, comm_O)
      TYPE(CHR10_VECT), INTENT(INOUT) :: var
      INTEGER, INTENT(IN)             :: rank, tag
      INTEGER, OPTIONAL, INTENT(IN)   :: comm_O
      INTEGER                         :: communicator
      INTEGER                         :: M, N, i, IErr
      CHARACTER(LEN=*),PARAMETER      :: Sub = 'Send_CHR10_VECT'

      IF(AllocQ(var%Alloc)) THEN
         CALL Delete(var)
      ENDIF
      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      CALL Recv(M, rank, tag, comm_O = communicator)
      CALL Recv(N, rank, tag, comm_O = communicator)
      CALL New(var, N, M_O = M)
      DO i = M, N
         CALL MPI_Recv(var%C(i), 10, MPI_CHARACTER, rank, tag, communicator, MPI_STATUS_IGNORE, IErr)
         CALL ErrChk(IErr,Sub)
      ENDDO
   END SUBROUTINE Recv_CHR10_VECT

   SUBROUTINE Recv_PBCInfo(var, rank, tag, comm_O)
      TYPE(PBCInfo), INTENT(INOUT)    :: var
      INTEGER, INTENT(IN)             :: rank, tag
      INTEGER, OPTIONAL, INTENT(IN)   :: comm_O
      INTEGER                         :: communicator
      INTEGER                         :: M, N, i, IErr
      CHARACTER(LEN=*),PARAMETER      :: Sub = 'Send_PBCInfo'

      IF(AllocQ(var%Alloc)) THEN
         CALL Delete(var)
      ENDIF
      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      CALL New(var)
      CALL Recv(var%Dimen, rank, tag, comm_O = communicator)
      CALL Recv(var%PFFMaxEll, rank, tag, comm_O = communicator)
      CALL Recv(var%PFFWelSep, rank, tag, comm_O = communicator)
      CALL Recv(var%InVecForm, rank, tag, comm_O = communicator)
      CALL Recv(var%InAtomCrd, rank, tag, comm_O = communicator)
      CALL Recv(var%Translate, rank, tag, comm_O = communicator)
      CALL Recv(var%CellVolume, rank, tag, comm_O = communicator)
      CALL Recv(var%Epsilon, rank, tag, comm_O = communicator)
      CALL Recv(var%DipoleFAC, rank, tag, comm_O = communicator)
      CALL Recv(var%AutoW, 3, rank, tag, comm_O = communicator)
      CALL Recv(var%SuperCell, 3, rank, tag, comm_O = communicator)
      CALL Recv(var%CellCenter, 3, rank, tag, comm_O = communicator)
      CALL Recv(var%TransVec, 3, rank, tag, comm_O = communicator)
      CALL Recv(var%BoxShape, rank, tag, comm_O = communicator)
      CALL Recv(var%InvBoxSh, rank, tag, comm_O = communicator)
      CALL Recv(var%LatFrc, rank, tag, comm_O = communicator)
   END SUBROUTINE Recv_PBCInfo

   SUBROUTINE Recv_CellSet(var, rank, tag, comm_O)
      TYPE(CellSet), INTENT(INOUT)    :: var
      INTEGER, INTENT(IN)             :: rank, tag
      INTEGER, OPTIONAL, INTENT(IN)   :: comm_O
      INTEGER                         :: communicator
      INTEGER                         :: M, N, i, IErr
      CHARACTER(LEN=*),PARAMETER      :: Sub = 'Send_CellSet'

      IF(AllocQ(var%Alloc)) THEN
         CALL Delete(var)
      ENDIF
      IF(PRESENT(comm_O)) THEN
        communicator = comm_O
      ELSE
        communicator = MONDO_COMM
      ENDIF
      CALL Recv(var%NCells, rank, tag, comm_O = communicator)
      CALL New(var, var%NCells)
      CALL Recv(var%Radius, rank, tag, comm_O = communicator)
      CALL Recv(var%CellCarts, rank, tag, comm_O = communicator)
   END SUBROUTINE Recv_CellSet

   !===============================================================================
   !     WRAPPERS FOR NON-BLOCKING (IMEDIATE) SENDS AND RECIEVES
   !===============================================================================
   FUNCTION ISend_DBL_SCLR(Snd,To,Tag)
      REAL(DOUBLE),INTENT(INOUT)  :: Snd
      INTEGER,     INTENT(IN)     :: To,Tag
      INTEGER                     :: ISend_DBL_SCLR,IErr
      CHARACTER(LEN=*),PARAMETER  :: Sub='ISend_DBL_SCLR'

      CALL MPI_ISEND(Snd,1,MPI_DOUBLE_PRECISION,To,Tag,MONDO_COMM,ISend_DBL_SCLR,IErr)
      CALL ErrChk(IErr,Sub)
   END FUNCTION ISend_DBL_SCLR

   FUNCTION ISend_DBL_VECT(Snd,N,To,Tag,M_O)
      TYPE(DBL_VECT),INTENT(INOUT) :: Snd
      INTEGER,       INTENT(IN)    :: N,To,Tag
      INTEGER                      :: ISend_DBL_VECT,L,M,IErr
      INTEGER,OPTIONAL,INTENT(IN)  :: M_O
      CHARACTER(LEN=*),PARAMETER   :: Sub='ISend_DBL_VECT'

      M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
      CALL MPI_ISEND(Snd%D(M:N),L,MPI_DOUBLE_PRECISION,To,Tag,MONDO_COMM,ISend_DBL_VECT,IErr)
      CALL ErrChk(IErr,Sub)
   END FUNCTION ISend_DBL_VECT

   FUNCTION ISend_INT_SCLR(Snd,To,Tag)
      INTEGER,INTENT(INOUT)       :: Snd
      INTEGER,INTENT(IN)          :: To,Tag
      INTEGER                     :: ISend_INT_SCLR,IErr
      CHARACTER(LEN=*),PARAMETER  :: Sub='ISend_INT_SCLR'

      !WRITE(*,*)' isend isend isend isend isend isend isend isend isend isend '
      !WRITE(*,*)'MyId = ',MyId,' To = ',To,' Tag = ',Tag
      CALL MPI_ISEND(Snd,1,MPI_INTEGER,To,Tag,MONDO_COMM,ISend_INT_SCLR,IErr)
      CALL ErrChk(IErr,Sub)
   END FUNCTION ISend_INT_SCLR

   FUNCTION ISend_INT_VECT(Snd,N,To,Tag,M_O)
      TYPE(INT_VECT),INTENT(INOUT) :: Snd
      INTEGER,       INTENT(IN)    :: N,To,Tag
      INTEGER                      :: ISend_INT_VECT,L,M,IErr
      INTEGER,OPTIONAL,INTENT(IN)  :: M_O
      CHARACTER(LEN=*),PARAMETER   :: Sub='ISend_INT_VECT'

      M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
      !if(myid==root)then
      !   WRITE(*,*)' isend isend isend isend isend isend isend isend isend isend '
      !   WRITE(*,*)' M = ',M,' N = ',N,' L = ',L
      !   WRITE(*,*)'MyId = ',MyId,' To = ',To,' Tag = ',Tag
      !   WRITE(*,*)' SIZE SND = ',SIZE(Snd%I)
      !   WRITE(*,*)' Snd = ',Snd%I(M:N)
      !endif
      CALL MPI_ISEND(Snd%I(M:N),L,MPI_INTEGER,To,Tag,MONDO_COMM,ISend_INT_VECT,IErr)
      CALL ErrChk(IErr,Sub)
   END FUNCTION ISend_INT_VECT

   FUNCTION IRecv_DBL_SCLR(Rec,From,Tag)
      REAL(DOUBLE),INTENT(INOUT) :: Rec
      INTEGER,INTENT(IN)         :: From,Tag
      INTEGER                    :: IRecv_DBL_SCLR,IErr
      CHARACTER(LEN=*),PARAMETER :: Sub='IRecv_DBL_SCLR'

      CALL MPI_IRECV(Rec,1,MPI_DOUBLE_PRECISION,From,Tag,MONDO_COMM,IRecv_DBL_SCLR,IErr)
      CALL ErrChk(IErr,Sub)
   END FUNCTION IRecv_DBL_SCLR

   FUNCTION IRecv_DBL_VECT(Rec,N,From,Tag,M_O)
      TYPE(DBL_VECT),INTENT(INOUT)  :: Rec
      INTEGER,INTENT(IN)            :: From,N,Tag
      INTEGER                       :: IRecv_DBL_VECT,L,M,IErr
      INTEGER,INTENT(IN),OPTIONAL   :: M_O
      CHARACTER(LEN=*),PARAMETER    :: Sub='IRecv_DBL_VECT'

      M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
      CALL MPI_IRECV(Rec%D(M:N),L,MPI_DOUBLE_PRECISION,From,Tag,MONDO_COMM,IRecv_DBL_VECT,IErr)
      CALL ErrChk(IErr,Sub)
   END FUNCTION IRecv_DBL_VECT

   FUNCTION IRecv_INT_SCLR(Rec,From,Tag)
      INTEGER,INTENT(INOUT)       :: Rec
      INTEGER,INTENT(IN)          :: From,Tag
      INTEGER                     :: IRecv_INT_SCLR,IErr
      CHARACTER(LEN=*),PARAMETER  :: Sub='IRecv_INT_SCLR'

      !WRITE(*,*)' irecv irecv irecv irecv irecv irecv irecv irecv irecv irecv irecv '
      !WRITE(*,*)' MyId = ',MyId,' From = ',From,' Tag = ',Tag
      CALL MPI_IRECV(Rec,1,MPI_INTEGER,From,Tag,MONDO_COMM,IRecv_INT_SCLR,IErr)
      CALL ErrChk(IErr,Sub)
   END FUNCTION IRecv_INT_SCLR

   FUNCTION IRecv_INT_VECT(Rec,N,From,Tag,M_O)
      TYPE(INT_VECT),INTENT(INOUT)  :: Rec
      INTEGER,INTENT(IN)            :: From,N,Tag
      INTEGER                       :: IRecv_INT_VECT,L,M,IErr
      INTEGER,INTENT(IN),OPTIONAL   :: M_O
      CHARACTER(LEN=*),PARAMETER    :: Sub='IRecv_INT_VECT'

      M=1; IF(PRESENT(M_O))M=M_O; L=N-M+1
      !if(myid==root)then
      !   WRITE(*,*)' irecv irecv irecv irecv irecv irecv irecv irecv irecv irecv irecv '
      !   WRITE(*,*)' M = ',M,' N = ',N,' L = ',L
      !   WRITE(*,*)' MyId = ',MyId,' From = ',From,' Tag = ',Tag
      !   WRITE(*,*)' SIZE RCV = ',SIZE(Rec%I)
      !endif
      CALL MPI_IRECV(Rec%I(M:N),L,MPI_INTEGER,From,Tag,MONDO_COMM,IRecv_INT_VECT,IErr)
      CALL ErrChk(IErr,Sub)
   END FUNCTION IRecv_INT_VECT

   !===============================================================================
   !     WRAPERS FOR WAIT ROUTINES
   !===============================================================================
   SUBROUTINE WaitAll(Reqs,Mssg)
      TYPE(INT_VECT),INTENT(INOUT)    :: Reqs
      TYPE(INT_RNK2)                  :: Stats
      INTEGER                         :: I,N,IErr
      CHARACTER(LEN=*), OPTIONAL      :: Mssg
      CHARACTER(LEN=*),PARAMETER      :: Sub='MWaitAll'

      IF(PRESENT(Mssg))THEN
         WRITE(*,*)'WAITALL:'
         CALL PSpew(Mssg,Reqs)
         WRITE(*,*)'...aiting waiting waiting waiti...'
      ENDIF
      ! Check to see if all requests are NULL
      N=SIZE(Reqs%I)
      DO I=1,N
         IF(Reqs%I(I)/=MPI_REQUEST_NULL)  &
         GOTO 10 ! Old habits die hard ...
      ENDDO
      RETURN
      10 CONTINUE
      ! Wait ...
      CALL New(Stats,(/MPI_STATUS_SIZE,N/))
      CALL MPI_WAITALL(N,Reqs%I(1:N),Stats%I,IErr)
      CALL ErrChk(IErr,Sub)
      CALL Delete(Stats)
   END SUBROUTINE WaitAll

   SUBROUTINE WaitSome(Reqs,ToDo)
      TYPE(INT_VECT),             INTENT(INOUT)  :: Reqs,ToDo
      TYPE(INT_VECT)                             :: Indxs
      TYPE(INT_RNK2)                             :: Stats
      INTEGER                                    :: I,L,N,IErr
      CHARACTER(LEN=*),PARAMETER                 :: Sub='WaitSome'

      ! Allocation check
      IF(AllocQ(ToDo%Alloc))CALL Delete(ToDo)
      ! Check to see if all requests are NULL
      N=SIZE(Reqs%I)
      DO I=1,N
         IF(Reqs%I(I)/=MPI_REQUEST_NULL)  &
         GOTO 10 ! Old habits die hard ...
      ENDDO
      CALL New(ToDo,1)
      ToDo%I=FAIL
      RETURN
      10 CONTINUE
      ! Wait ...
      CALL New(Indxs,N)
      CALL New(Stats,(/MPI_STATUS_SIZE,N/))
      CALL MPI_WAITSOME(N,Reqs%I,L,Indxs%I,Stats%I,IErr)
      CALL ErrChk(IErr,Sub)
      CALL New(ToDo,L)
      ToDo%I(1:L)=Indxs%I(1:L)
      CALL Delete(Indxs)
      CALL Delete(Stats)
   END SUBROUTINE WaitSome

   !===============================================================================
   !     WRAPER FOR MPI_TYPE_FREE
   !===============================================================================
   SUBROUTINE FreeType(Type)
      INTEGER, INTENT(INOUT)      :: Type
      INTEGER                     :: IErr
      CHARACTER(LEN=*),PARAMETER  :: Sub='FreeType'

      IF(Type/=MPI_DATATYPE_NULL)THEN
         CALL MPI_TYPE_FREE(Type,IErr)
         CALL ErrChk(IErr,Sub)
      ENDIF
   END SUBROUTINE FreeType

   !===============================================================================
   !     NODE ALLINGMENT FOR MONDO_COMM
   !===============================================================================
   SUBROUTINE AlignNodes(String_O)
      INTEGER                        :: IErr
      CHARACTER(LEN=*),OPTIONAL      :: String_O
      CHARACTER(LEN=*),PARAMETER     :: Sub='AlignNodes'

      CALL MondoLog(DEBUG_NONE, "AlignNodes", "waiting at barrier, MONDO_COMM = "//TRIM(IntToChar(MONDO_COMM)))
      CALL MPI_BARRIER(MONDO_COMM,IErr)
      CALL ErrChk(IErr,Sub)
      IF(PRESENT(String_O))THEN
         CALL MondoLog(DEBUG_NONE, "FreeON", String_O, "Node # "//TRIM(IntToChar(MyId)))
      ENDIF
   END SUBROUTINE AlignNodes

   !===============================================================================
   !     NODE ALLINGMENT FOR MPI_COMM_WORLD
   !===============================================================================
   SUBROUTINE AlignClones(String_O)
      INTEGER                        :: IErr
      CHARACTER(LEN=*),OPTIONAL      :: String_O
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: String
      CHARACTER(LEN=*),PARAMETER     :: Sub='AlignNodes'

      CALL MPI_BARRIER(MPI_COMM_WORLD,IErr)
      CALL ErrChk(IErr,Sub)
      IF(PRESENT(String_O))THEN
         String="Node#"//TRIM(IntToChar(MyId))//' :: '//TRIM(String_O)
         WRITE(*,*)TRIM(String)
      ENDIF
   END SUBROUTINE AlignClones

   !===============================================================================
   !     Error checking
   !===============================================================================
   SUBROUTINE ErrChk(IErr,Caller,Seqnce_O)
      INTEGER                              :: IErr
      CHARACTER(LEN=*),INTENT(IN)          :: Caller
      CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Seqnce_O

      IF(IErr/=MPI_SUCCESS)THEN
         IF(PRESENT(Seqnce_O))THEN
            CALL HaltMPI('>>>MondoMPI error: '//TRIM(Seqnce_O)//'$'//TRIM(Caller)//'.',IErr)
         ELSE
            CALL HaltMPI('>>>MondoMPI error: '//TRIM(Caller)//'.',IErr)
         ENDIF
      ENDIF
   END SUBROUTINE

   SUBROUTINE InTurn(Mssg)
      CHARACTER(LEN=*)            :: Mssg
      INTEGER                     :: I,IErr
      CHARACTER(LEN=*),PARAMETER  :: Sub='InTurn'

      CALL MPI_BARRIER(MONDO_COMM,IErr)
      CALL ErrChk(IErr,Sub)
      DO I=0,NPrc-1
         CALL MPI_BARRIER(MONDO_COMM,IErr)
         CALL ErrChk(IErr,Sub)
         IF(MyId==I)WRITE(*,*)'MyId = ',MyId,'$',TRIM(Mssg)
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
      55 FORMAT('MyId = ',I2,', ',A,' = ',1000(I20,', '))
   ENDSUBROUTINE PSpew_INT_VECT

   SUBROUTINE PSpew_DBL_VECT(Name,A)
      CHARACTER(LEN=*) :: Name
      TYPE(DBL_VECT)   :: A
      INTEGER          :: I

      DO I=0,NPrc-1
         CALL AlignNodes()
         IF(MyId==I)WRITE(*,55)MyId,TRIM(Name),A%D
      ENDDO
      55 FORMAT('MyId = ',I2,', ',A,' = ',1000(D9.3,', '))
   ENDSUBROUTINE PSpew_DBL_VECT

   SUBROUTINE AlignFrontends()
      CHARACTER                             :: shutdownBuffer
      INTEGER, DIMENSION(:), ALLOCATABLE    :: shutdownRequest, IErr
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: shutdownStatus
      INTEGER                               :: I
      LOGICAL                               :: allDone, done

      ! Allocate memory for shutdown negotiations.
      ALLOCATE(shutdownRequest(NPrc))
      ALLOCATE(shutdownStatus(NPrc, MPI_STATUS_SIZE))
      ALLOCATE(IErr(NPrc))

      IF(MyID == 0) THEN
         ! Set up non-blocking send to all other ranks.
         DO I = 1, NPrc
            CALL MondoLog(DEBUG_MAXIMUM, "FreeON", "sending shutdown message to frontend rank "//TRIM(IntToChar(I)))
            CALL MPI_Isend(shutdownBuffer, 1, MPI_CHARACTER, I, FRONTEND_TAG, MPI_COMM_WORLD, shutdownRequest(I), IErr(I))
         ENDDO

         allDone = .FALSE.
         DO WHILE(.NOT. allDone)
            ! Check for the other ranks.
            allDone = .TRUE.
            DO I = 1, NPrc
               CALL MPI_Test(shutdownRequest(I), done, shutdownStatus(I,:), IErr(I))
               IF(.NOT. done) THEN
                  allDone = .FALSE.
                  EXIT
               ENDIF
            ENDDO
         ENDDO
         CALL MondoLog(DEBUG_MAXIMUM, "FreeON", "all frontends are done")
      ELSE
         ! Wait a bit so the log message ends up in the out file.
         CALL FreeONSleep(2)

         ! The log messages currently only end up in the out file if the
         ! logLevel is DEBUG_NONE because only the rank 0 frontend calls
         ! ParsePrintFlags(). We would have to communicate the PrintFlags to
         ! rank > 0 frontends to get MondoLog() to print messages here.

         ! Set up non-blocking receive from rank 0 front-end.
         CALL MondoLog(DEBUG_MAXIMUM, "FreeON", "waiting for shutdown message from Frontend rank 0", "Frontend rank "//TRIM(IntToChar(MyID)))
         CALL MPI_Irecv(shutdownBuffer, 1, MPI_CHARACTER, 0, FRONTEND_TAG, MPI_COMM_WORLD, shutdownRequest(1), IErr(1))

         done = .FALSE.
         DO WHILE(.NOT. done)
            CALL FreeONSleep(2)
            CALL MPI_Test(shutdownRequest(1), done, shutdownStatus(1,:), IErr(1))
         ENDDO
         CALL MondoLog(DEBUG_MAXIMUM, "FreeON", "received shutdown message from rank 0", "Frontend rank "//TRIM(IntToChar(MyID)))
      ENDIF

      ! Free memory.
      DEALLOCATE(shutdownRequest)
      DEALLOCATE(shutdownStatus)
      DEALLOCATE(IErr)

      ! Shut down MPI.
      CALL MPI_FINALIZE(IErr(1))
   END SUBROUTINE AlignFrontends

   ! Initialize a lock object.
   SUBROUTINE AllocateLock (lock, communicator)
      TYPE(FreeONLock), INTENT(INOUT)     :: lock
      INTEGER, INTENT(IN)                 :: communicator
      INTEGER                             :: IErr
      INTEGER                             :: i
      INTEGER, DIMENSION(:), ALLOCATABLE  :: arrayOfBlocklengths, arrayOfDisplacements

      IF(AllocQ(lock%Alloc)) THEN
         CALL MondoHalt(MPIS_ERROR, "lock is already allocated")
      ENDIF

      ! Get information on size of current communicator.
      lock%communicator = communicator
      CALL MPI_COMM_SIZE(lock%communicator, lock%lockSize, IErr)
      CALL MPI_COMM_RANK(lock%communicator, lock%lockRank, IErr)

      CALL MondoLog(DEBUG_NONE, "AllocateLock", "allocating lock", "rank "//TRIM(IntToChar(lock%lockRank)))

      ! Create window object.
      CALL MPI_INFO_CREATE(lock%windowInfo, IErr)
      IF(lock%lockRank == ROOT) THEN
         ALLOCATE(lock%waitflag(lock%lockSize))
         DO i = 1, lock%lockSize
            lock%waitflag(i) = 0
         ENDDO
         CALL MondoLog(DEBUG_NONE, "AllocateLock", "creating window", "rank "//TRIM(IntToChar(lock%lockRank)))
         CALL MondoLog(DEBUG_NONE, "AllocateLock", "size = "//TRIM(IntToChar(lock%lockSize)), "rank "//TRIM(IntToChar(lock%lockRank)))
         CALL MPI_WIN_CREATE(lock%waitflag, lock%lockSize, 1, lock%windowInfo, lock%communicator, lock%window, IErr)
      ELSE
         ! Only rank 0 stores the actual flags. If we are rank > 0 we create a
         ! window of 0 size.
         CALL MondoLog(DEBUG_NONE, "AllocateLock", "creating empty window", "rank "//TRIM(IntToChar(lock%lockRank)))
         CALL MPI_WIN_CREATE(lock%waitflag, 0, 1, lock%windowInfo, lock%communicator, lock%window, IErr)
      ENDIF

      ! Create local copy of waitflag array. This array contains size-1 flags,
      ! it excludes the flag for rank.
      ALLOCATE(arrayOfBlocklengths(lock%lockSize-1))
      ALLOCATE(arrayOfDisplacements(lock%lockSize-1))
      DO i = 1, lock%lockSize
         IF(i-1 < lock%lockRank) THEN
            arrayOfBlocklengths(i) = 1
            arrayOfDisplacements(i) = i-1
         ELSE IF(i-1 > lock%lockRank) THEN
            arrayOfBlocklengths(i-1) = 1
            arrayOfDisplacements(i-1) = i-1
         ENDIF
      ENDDO

      ! Create new datatype to represent all waitflags except for ours.
      CALL MondoLog(DEBUG_NONE, "AllocateLock", "creating new datatype", "rank "//TRIM(IntToChar(lock%lockRank)))
      CALL MPI_TYPE_INDEXED(lock%lockSize-1, arrayOfBlocklengths, arrayOfDisplacements, MPI_BYTE, lock%waitflagCopyType, IErr)
      CALL MPI_TYPE_COMMIT(lock%waitflagCopyType, IErr)

      ! Allocate local copy of waitflags.
      ALLOCATE(lock%waitflagCopy(lock%lockSize-1))

      ! Free some memory.
      DEALLOCATE(arrayOfBlocklengths)
      DEALLOCATE(arrayOfDisplacements)

      ! We are done allocating.
      lock%Alloc = ALLOCATED_TRUE
      lock%lockAcquired = .FALSE.

      END SUBROUTINE AllocateLock

      SUBROUTINE FreeLock (lock)
      TYPE(FreeONLock), INTENT(INOUT) :: lock
      INTEGER                         :: IErr

      IF(.NOT. AllocQ(lock%Alloc)) THEN
         CALL MondoHalt(MPIS_ERROR, "lock not allocated")
      ENDIF

      ! Free window.
      CALL MPI_WIN_FREE(lock%window, IErr)

      ! Free the waitflag datatype.
      CALL MPI_TYPE_FREE(lock%waitflagCopyType, IErr)

      ! Deallocate memory.
      IF(lock%lockRank == ROOT) THEN
         DEALLOCATE(lock%waitflag)
      ENDIF
      DEALLOCATE(lock%waitflagCopy)

      ! We are done freeing.
      lock%Alloc = ALLOCATED_FALSE

   END SUBROUTINE FreeLock

   SUBROUTINE AcquireLock (lock, lockType)
      TYPE(FreeONLock), INTENT(INOUT) :: lock
      INTEGER, INTENT(IN)             :: lockType
      INTEGER(KIND = INT1)            :: value
      INTEGER                         :: i
      INTEGER                         :: rank
      INTEGER                         :: IErr
      INTEGER                         :: emptyBuffer
      CHARACTER(LEN=200)              :: outputFormat
      CHARACTER(LEN=200)              :: outputBuffer

      ! Sanity check.
      IF(lock%lockAcquired) THEN
         CALL MondoHalt(MPIS_ERROR, "this lock is already acquired")
      ENDIF

      IF(lockType == FreeONLockExclusive) THEN
         CALL MondoLog(DEBUG_NONE, "AcquireLock", "acquiring exclusive lock", "rank "//TRIM(IntToChar(lock%lockRank)))
      ELSE
         CALL MondoLog(DEBUG_NONE, "AcquireLock", "acquiring shared lock", "rank "//TRIM(IntToChar(lock%lockRank)))
      ENDIF

      ! Store the lock type.
      lock%lockType = lockType

      ! Get exclusive lock on window.
      CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE, ROOT, 0, lock%window, IErr)

      ! Get waitflags.
      CALL MPI_GET(lock%waitflagCopy, lock%lockSize-1, MPI_BYTE, ROOT, 0, 1, lock%waitflagCopyType, lock%window, IErr)

      ! Put our waitflag, indicating that we would like to acquire this lock.
      value = 1
      CALL MPI_PUT(value, 1, MPI_BYTE, ROOT, lock%lockRank, 1, MPI_BYTE, lock%window, IErr)

      ! Unlock window.
      CALL MPI_WIN_UNLOCK(ROOT, lock%window, IErr)

      ! Debugging.
      WRITE(outputFormat, "(A,I3,A)") "(", lock%lockSize-1, "I2)"
      WRITE(outputBuffer, outputFormat) lock%waitflagCopy
      CALL MondoLog(DEBUG_NONE, "AcquireLock", "waitflagCopy = ["//TRIM(outputBuffer)//" ]", "rank "//TRIM(IntToChar(lock%lockRank)))

      IF(lock%lockType == FreeONLockExclusive) THEN
         ! Check whether another rank is already holding this lock.
         rank = -1
         DO i = 1, lock%lockSize-1
            IF(lock%waitflagCopy(i) == 1) THEN
               ! Remember that ranks start counting with 0.
               rank = i-1
               CALL MondoLog(DEBUG_NONE, "AcquireLock", "rank "//TRIM(IntToChar(rank))//" is holding the lock already", "rank "//TRIM(IntToChar(lock%lockRank)))
               EXIT
            ENDIF
         ENDDO

         IF(rank >= 0) THEN
            ! Somebody else is holding this lock. We have to wait.
            CALL MondoLog(DEBUG_NONE, "AcquireLock", "another rank is holding the lock, I am going to wait", "rank "//TRIM(IntToChar(lock%lockRank)))
            CALL MPI_RECV(emptyBuffer, 0, MPI_BYTE, MPI_ANY_SOURCE, LOCK_TAG, lock%communicator, MPI_STATUS_IGNORE, IErr)
            CALL MondoLog(DEBUG_NONE, "AcquireLock", "received notification, acquiring lock", "rank "//TRIM(IntToChar(lock%lockRank)))
         ELSE
            CALL MondoLog(DEBUG_NONE, "AcquireLock", "no rank is holding the lock, done", "rank "//TRIM(IntToChar(lock%lockRank)))
         ENDIF
      ENDIF

      ! We successfully acquired this lock.
      CALL MondoLog(DEBUG_NONE, "AcquireLock", "lock acquired", "rank "//TRIM(IntToChar(lock%lockRank)))
      lock%lockAcquired = .TRUE.

   END SUBROUTINE AcquireLock

   SUBROUTINE ReleaseLock (lock)
      TYPE(FreeONLock), INTENT(INOUT) :: lock
      INTEGER(KIND = INT1)            :: value
      INTEGER                         :: i
      INTEGER                         :: nextrank
      INTEGER                         :: IErr
      INTEGER                         :: emptyBuffer
      CHARACTER(LEN=200)              :: outputFormat
      CHARACTER(LEN=200)              :: outputBuffer

      ! Sanity check.
      IF(.NOT. lock%lockAcquired) THEN
         CALL MondoHalt(MPIS_ERROR, "this lock has not been acquired")
      ENDIF

      IF(lock%lockType == FreeONLockExclusive) THEN
         CALL MondoLog(DEBUG_NONE, "ReleaseLock", "releasing exclusive lock", "rank "//TRIM(IntToChar(lock%lockRank)))
      ELSE
         CALL MondoLog(DEBUG_NONE, "ReleaseLock", "releasing shared lock", "rank "//TRIM(IntToChar(lock%lockRank)))
      ENDIF

      ! Get exclusive lock on window.
      CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE, ROOT, 0, lock%window, IErr)

      ! Get waitflags.
      CALL MPI_GET(lock%waitflagCopy, lock%lockSize-1, MPI_BYTE, ROOT, 0, 1, lock%waitflagCopyType, lock%window, IErr)

      ! Put our waitflag, indicating that we would like to acquire this lock.
      value = 0
      CALL MPI_PUT(value, 1, MPI_BYTE, ROOT, lock%lockRank, 1, MPI_BYTE, lock%window, IErr)

      ! Unlock window.
      CALL MPI_WIN_UNLOCK(ROOT, lock%window, IErr)

      ! Debugging.
      WRITE(outputFormat, "(A,I3,A)") "(", lock%lockSize-1, "I2)"
      WRITE(outputBuffer, outputFormat) lock%waitflagCopy
      CALL MondoLog(DEBUG_NONE, "ReleaseLock", "waitflagCopy = ["//TRIM(outputBuffer)//" ]", "rank "//TRIM(IntToChar(lock%lockRank)))

      ! We need to notify anybody only in the case of holding an exclusive
      ! lock.
      IF(lock%lockType == FreeONLockExclusive) THEN
         ! In case someone else is waiting on the lock, notify them that we are
         ! done.
         nextrank = -1
         DO i = 1, lock%lockSize-1
            IF(lock%waitflagCopy(i) == 1) THEN
               ! Remember that ranks start counting with 0.
               nextrank = i-1
               EXIT
            ENDIF
         ENDDO

         IF(nextrank >= 0) THEN
            ! Somebody else is holding this lock. We should notify them that we
            ! are releasing the lock.
            nextrank = -1
            DO i = lock%lockRank+1, lock%lockSize-1
               IF(lock%waitflagCopy(i) == 1) THEN
                  ! Remember that ranks start counting with 0.
                  nextrank = (i-1)+1
                  EXIT
               ENDIF
            ENDDO

            IF(nextrank < 0) THEN
               DO i = 1, lock%lockRank
                  IF(lock%waitflagCopy(i) == 1) THEN
                     ! Remember that ranks start counting with 0.
                     nextrank = i-1
                     EXIT
                  ENDIF
               ENDDO
            ENDIF

            CALL MondoLog(DEBUG_NONE, "ReleaseLock", "another rank is waiting for the lock, I am going to notify rank "//TRIM(IntToChar(nextrank)), "rank "//TRIM(IntToChar(lock%lockRank)))
            CALL MPI_SEND(emptyBuffer, 0, MPI_BYTE, nextrank, LOCK_TAG, lock%communicator, IErr)
         ELSE
            CALL MondoLog(DEBUG_NONE, "ReleaseLock", "no other rank is waiting for the lock, done", "rank "//TRIM(IntToChar(lock%lockRank)))
         ENDIF
      ENDIF

      ! We successfully released this lock.
      CALL MondoLog(DEBUG_NONE, "ReleaseLock", "lock released", "rank "//TRIM(IntToChar(lock%lockRank)))
      lock%lockAcquired = .FALSE.

   END SUBROUTINE ReleaseLock

#endif

END MODULE
