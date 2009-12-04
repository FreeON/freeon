!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
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
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------

#include "MondoConfig.h"

MODULE FastMatrices
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalObjects
   USE MemMan
   USE Order
   USE LinAlg
#ifdef PARALLEL
   USE MondoMPI
#endif
   USE Thresholding
   IMPLICIT NONE
!======================================================================
!   LINKED ROW LIST WITH COLOUMN INDEXED SEARCH TREE
!   A FAST O(N LG N) SPARSE MATRIX DATA STRUCTURE FOR ALL PROCCEDURES
!======================================================================
   TYPE FASTMAT
      INTEGER                  :: Alloc   !-- Allocation key
      INTEGER                  :: Nodes   !-- Number of nodes in this SRST
      INTEGER                  :: Row     !-- Row index of this link
      INTEGER                  :: NSMat   !-- Nbr spin matrices
      TYPE(SRST),POINTER   :: RowRoot !-- Row link to a sparse row search tree
      TYPE(FASTMAT), POINTER   :: Next    !-- Next row in linked list
   END TYPE FASTMAT
!======================================================================
!   SPARSE ROW SEARCH TREE: A FAST (LG N) SPARSE VECTOR DS
!======================================================================
    TYPE SRST
       INTEGER                  :: Alloc  !-- Allocation key
       INTEGER                  :: Row    !-- Row number
       INTEGER                  :: Tier   !-- Tree depth
       INTEGER                  :: Number !-- This nodes number in the tree
       INTEGER                  :: L,R    !-- Left/right coloumn interval bounds
       real(double) :: Part
       TYPE(SRST), POINTER      :: Left
       TYPE(SRST), POINTER      :: Right
       TYPE(SRST), POINTER      :: Next
       REAL(DOUBLE), POINTER, DIMENSION(:,:)  :: MTrix  !-- Matrix block
    END TYPE SRST
    INTEGER :: SRSTCount
    TYPE(SRST),POINTER :: GlobalP
    TYPE(SRST),POINTER :: GlobalP1
!======================================================================
#ifdef PARALLEL
  CONTAINS
!======================================================================

!======================================================================
! COMPUTE BCSR MATRIX DIMENSIONS CORESPONDING TO A FAST MATRIX
!======================================================================
  FUNCTION MatDimensions_1(A,RowLimits,ColLimits) RESULT(Dim_res)
    TYPE(FASTMAT),POINTER :: A,C
    TYPE(SRST),POINTER :: P
    INTEGER,DIMENSION(2) :: RowLimits,ColLimits
    INTEGER,DIMENSION(3) :: Dim_Res
    INTEGER :: RowNum,BlkNum,MtxEleNum,M,N,Row,Col
    !
    ! Set some pointers.
    NULLIFY(C,P)
    !
    RowNum = 0
    BlkNum = 0
    MtxEleNum = 0
    C => A%Next
    DO
      IF(.NOT. ASSOCIATED(C)) EXIT
      Row = C%Row
      IF(Row >= RowLimits(1) .AND. Row <= RowLimits(2)) THEN
        RowNum = RowNum + 1
        M = BSiz%I(Row)
        !! calculate the number of blks and number of matrix elements
        P => C%RowRoot
        DO
          IF(.NOT. ASSOCIATED(P)) EXIT
          Col = P%L
          IF(Col == P%R .AND. Col >= ColLimits(1) .AND. &
             Col <= ColLimits(2) ) THEN
            BlkNum = BlkNum + 1
            N = BSiz%I(Col)
            MtxEleNum = MtxEleNum + M*N
          ENDIF
          P => P%Next
        ENDDO
      ENDIF
      C => C%Next
    ENDDO
    Dim_Res = (/RowNum+1,BlkNum,MtxEleNum/)
  END FUNCTION MatDimensions_1

!======================================================================
  !
  SUBROUTINE Set_BCSR_EQ_DFASTMAT(C,A)
    TYPE(FASTMAT),POINTER :: A
    TYPE(BCSR) :: B,C
    INTEGER :: &
      PrevColSize,LocalNBlks,PrevBlkSize,NewPt,NAtms,IErr,GBNBlks,GBNNon0,I
    TYPE(INT_VECT) :: CA,CB,CN,DispN,DispB,DispA

    CALL GetLocalBCSR(B,A,(/Beg%I(MyID),End%I(MyID)/))

    NAtms = End%I(MyID)-Beg%I(MyID)+1
    IF(MyID == 0) THEN
      CALL New(CA,NPrc-1,M_O=0)
      CALL New(CB,NPrc-1,M_O=0)
      CALL New(CN,NPrc-1,M_O=0)
      CALL New(DispN,NPrc-1,M_O=0)
      CALL New(DispB,NPrc-1,M_O=0)
      CALL New(DispA,NPrc-1,M_O=0)
    ELSE
      ! Allocate something to allow full debug...
      CALL New(CA,0,0)
      CALL New(CB,0,0)
      CALL New(CN,0,0)
    ENDIF

    CALL MPI_Gather(NAtms,1,MPI_INTEGER,CA%I(0),1,MPI_INTEGER,0,MONDO_COMM,IErr)
    CALL MPI_Gather(B%NBlks,1,MPI_INTEGER,CB%I(0),1,MPI_INTEGER,0,MONDO_COMM,IErr)
    CALL MPI_Gather(B%NNon0,1,MPI_INTEGER,CN%I(0),1,MPI_INTEGER,0,MONDO_COMM,IErr)

    IF(MyID == 0) THEN
      GBNBlks = 0
      GBNNon0 = 0
      DO I = 0, NPrc-1
        GBNBlks = GBNBlks + CB%I(I)
        GBNNon0 = GBNNon0 + CN%I(I)
      ENDDO
      ! Messy output is unsightly and unessesary.  Protect with if defs in future...
      !      CALL OpenASCII(OutFile,Out)
      !      WRITE(Out,*) 'Set_BCSR_EQ_DFASTMAT: GBNBlks = ', GBNBlks,', GBNNon0 = ', GBNNon0
      !      CLOSE(Out,STATUS='KEEP')
      CALL New(C,(/NAtoms,GBNBlks,GBNNon0/))
      DispN%I(0) = 0
      DispB%I(0) = 0
      DispA%I(0) = 0
      DO I = 1, NPrc-1
        DispN%I(I) = DispN%I(I-1) + CN%I(I-1)
        DispB%I(I) = DispB%I(I-1) + CB%I(I-1)
        DispA%I(I) = DispA%I(I-1) + CA%I(I-1)
      ENDDO
    ENDIF

    IF(MyID == 0) THEN
      CALL MPI_GatherV(B%MTrix%D,B%NNon0,MPI_DOUBLE_PRECISION,C%MTrix%D,CN%I(0),DispN%I(0),MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
    ELSE
      CALL MPI_GatherV(B%MTrix%D,B%NNon0,MPI_DOUBLE_PRECISION,1.0D0,1,1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
    ENDIF

    IF(MyID == 0) THEN
      CALL MPI_GatherV(B%ColPt%I(1),B%NBlks,MPI_INTEGER,C%ColPt%I(1),CB%I(0),DispB%I(0),MPI_INTEGER,0,MONDO_COMM,IErr)
    ELSE
      CALL MPI_GatherV(B%ColPt%I(1),B%NBlks,MPI_INTEGER,1,1,1,MPI_INTEGER,0,MONDO_COMM,IErr)
    ENDIF


    ! take differences
    LocalNBlks = B%NBlks
    DO I = 1, LocalNBlks-1
      B%BlkPt%I(I) = B%BlkPt%I(I+1)-B%BlkPt%I(I)
    ENDDO
    B%BlkPt%I(LocalNBlks) = B%NNon0-B%BlkPt%I(LocalNBlks)+1

    IF(MyID == 0) THEN
      CALL MPI_GatherV(B%BlkPt%I(1),LocalNBlks,MPI_INTEGER,C%BlkPt%I(1),CB%I(0),DispB%I(0),MPI_INTEGER,0,MONDO_COMM,IErr)
      PrevBlkSize = C%BlkPt%I(1)
      C%BlkPt%I(1) = 1
      DO I = 2, GBNBlks
        ! new pointer = previous cumulative pointer + prev size
        NewPt = C%BlkPt%I(I-1) + PrevBlkSize
        PrevBlkSize = C%BlkPt%I(I)
        C%BlkPt%I(I) = NewPt
      ENDDO
    ELSE
      CALL MPI_GatherV(B%BlkPt%I(1),B%NBlks,MPI_INTEGER,1,1,1,MPI_INTEGER,0,MONDO_COMM,IErr)
    ENDIF

    DO I = 1, NAtms
      B%RowPt%I(I) = B%RowPt%I(I+1)-B%RowPt%I(I)
    ENDDO
    IF(MyID == 0) THEN
      CALL MPI_GatherV(B%RowPt%I(1),NAtms,MPI_INTEGER,C%RowPt%I(1),CA%I(0),DispA%I(0),MPI_INTEGER,0,MONDO_COMM,IErr)
      PrevColSize = C%RowPt%I(1)
      C%RowPt%I(1) = 1
      DO I = 2, NAtoms+1
        NewPt = C%RowPt%I(I-1) + PrevColSize
        PrevColSize = C%RowPt%I(I)
        C%RowPt%I(I) = NewPt
      ENDDO
    ELSE
      CALL MPI_GatherV(B%RowPt%I(1),NAtms,MPI_INTEGER,1,1,1,MPI_INTEGER,0,MONDO_COMM,IErr)
    ENDIF

    CALL Delete(CA)
    CALL Delete(CB)
    CALL Delete(CN)

    IF(MyID == 0) THEN
      CALL Delete(DispA)
      CALL Delete(DispB)
      CALL Delete(DispN)
    ENDIF
    CALL Delete(B)

  END SUBROUTINE Set_BCSR_EQ_DFASTMAT

!======================================================================
! GetLocalBCSR fails if there is an empty row in the FastMat.
!======================================================================
!   CONVERT A FAST MATRIX INTO A BCSR MATRIX
!======================================================================
  SUBROUTINE GetLocalBCSR(B,A,RowLimits_O)
    TYPE(BCSR) :: B
    TYPE(FASTMAT),POINTER :: A,R
    TYPE(SRST),POINTER :: C
    INTEGER,OPTIONAL,DIMENSION(2) :: RowLimits_O
    INTEGER,DIMENSION(2) :: RowLimits
    INTEGER,DIMENSION(3) :: MatDims
    INTEGER :: MN,P,RowOffSt,MtxBegInd,NAtms,M,N,Row,Col,RowNum
    !
    ! Set some pointers.
    NULLIFY(C,R)
    !
    IF(PRESENT(RowLimits_O)) THEN
      RowLimits = RowLimits_O
    ELSE
      RowLimits = (/1,NAtoms/)
    ENDIF
    MatDims = MatDimensions_1(A,RowLimits,(/1,NAtoms/))
    NAtms = End%I(MyID) - Beg%I(MyID) + 1
    NAtms = NAtms + 1
    IF(NAtms /= MatDims(1)) THEN
      WRITE(*,*) 'NAtms = ', NAtms
      WRITE(*,*) 'MatDims(1) = ', MatDims(1)
      WRITE(*,*) 'WARNING: NAtms is not equal to MatDims(1)'
      STOP 'ERR: Basic assumption of GetLocalBCSR is not valid!'
    ENDIF
    CALL New(B,N_O=(/MatDims(1)-1,MatDims(2),MatDims(3)/),OnAll_O=.TRUE.)
    RowNum = 0
    R => A%Next
    MtxBegInd = 1
    P = 1
    B%RowPt%I(1) = 1
    RowOffSt = Beg%I(MyID)-1
    DO
      IF(.NOT. ASSOCIATED(R)) EXIT
      Row = R%Row
      M = BSiz%I(Row)
      RowNum = RowNum+1
      IF(Row-RowOffSt /= RowNum) THEN
        WRITE(*,*) 'Row = ', Row, ', RowOffSt = ', RowOffSt
        WRITE(*,*) 'Row-RowOffSt = ', Row-RowOffSt
        WRITE(*,*) 'RowNum = ', RowNum
        STOP 'ERR: Row not equal to RowNum in GetLocalBCSR!'
      ENDIF
      C => R%RowRoot
      DO
        IF(.NOT. ASSOCIATED(C)) EXIT
        Col = C%L
        IF(Col == C%R) THEN
          N = BSiz%I(Col)
          MN = M*N
          B%MTrix%D(MtxBegInd:MtxBegInd+MN-1) = PACK(C%MTrix,.TRUE.)
          B%BlkPt%I(P) = MtxBegInd
          B%ColPt%I(P) = Col
          P = P + 1
          MtxBegInd = MtxBegInd + MN
          B%RowPt%I(RowNum+1) = P
        ENDIF
        C => C%Next
      ENDDO
      R => R%Next
    ENDDO
    B%NAtms = RowNum
    B%NBlks = P-1
    B%NNon0 = MtxBegInd-1
  END SUBROUTINE GetLocalBCSR
!=================================================================
  SUBROUTINE Delete_FastMat1(A,RowLimits_O)
    TYPE(FASTMAT),POINTER         :: A,R,NextR,P
    INTEGER,OPTIONAL,DIMENSION(2) :: RowLimits_O
    INTEGER,         DIMENSION(2) :: RowLimits
    INTEGER :: Row
    !
    ! Set some pointers.
    NULLIFY(R,NextR,P)
    !
    IF(PRESENT(RowLimits_O))THEN
      RowLimits=RowLimits_O
    ELSE
      RowLimits=(/1,1000000000/)
    ENDIF
    R => A
    NULLIFY(P)
    DO
      IF(.NOT. ASSOCIATED(R)) EXIT
      Row = R%Row
      !IF(MyID == 0) THEN
        !! WRITE(*,*) 'MyID = ',MyID, ' Testing row ', Row
      !ENDIF
      IF(Row >= RowLimits(1) .AND. Row <= RowLimits(2)) THEN
         !IF(MyID == 0) THEN
          !! WRITE(*,*) 'MyID = ',MyID, ' deleting row ', Row
         !ENDIF

        NextR => R%Next
        !! CALL Print_SRST_1(R%RowRoot)
        IF(ASSOCIATED(R%RowRoot)) CALL Delete_SRST_1(R%RowRoot)   !vw add that: if(associated(R%RowRoot))
        !!CALL Delete_SRST_1(R%RowRoot)   !vw add that: if(associated(R%RowRoot))
        DEALLOCATE(R)
        P%Next => NextR
        R => NextR
      ELSE
        P => R
        R => R%Next
      ENDIF
    ENDDO
  END SUBROUTINE Delete_FastMat1
!=================================================================
  RECURSIVE SUBROUTINE Print_SRST_1(A)
    TYPE(SRST),POINTER :: A
    IF(.NOT. ASSOCIATED(A%Left) .AND. .NOT. ASSOCIATED(A%Right)) THEN
      !IF(MyID == 0) THEN
        WRITE(*,*) 'Col = ', A%L
      !ENDIF

      IF(ASSOCIATED(A%Left) .OR. ASSOCIATED(A%Right)) THEN
        STOP 'ERR: Left and right must be NULL! '
      ENDIF
    ELSE
      IF(ASSOCIATED(A%Left)) THEN
        CALL Print_SRST_1(A%Left)
      ENDIF
      IF(ASSOCIATED(A%Right)) THEN
        CALL Print_SRST_1(A%Right)
      ENDIF
    ENDIF
  END SUBROUTINE Print_SRST_1

!=================================================================
  RECURSIVE SUBROUTINE Delete_SRST_1(A)
    TYPE(SRST),POINTER :: A

    IF(.NOT. ASSOCIATED(A)) THEN
      WRITE(*,*) 'ERR: A is null in Delete_SRST_1!'
      STOP
    ENDIF
    IF(.NOT. ASSOCIATED(A%Left) .AND. .NOT. ASSOCIATED(A%Right)) THEN
      IF(ASSOCIATED(A%Left) .OR. ASSOCIATED(A%Right)) THEN
        STOP 'ERR: Left and right must be NULL! '
      ENDIF
      !IF(MyID == 0) THEN
        !! WRITE(*,*) 'Delete_SRST_1 , Col = ', A%L, ' Col 2 = ', A%R
      !ENDIF
      IF(ASSOCIATED(A%MTrix)) THEN
        DEALLOCATE(A%MTrix)
      ENDIF
      DEALLOCATE(A)
      NULLIFY(A)
    ELSE
      IF(ASSOCIATED(A%Left)) THEN
        CALL Delete_SRST_1(A%Left)
      ENDIF
      IF(ASSOCIATED(A%Right)) THEN
        CALL Delete_SRST_1(A%Right)
      ENDIF
      DEALLOCATE(A)
      NULLIFY(A)
    ENDIF
  END SUBROUTINE Delete_SRST_1

!======================================================================
  SUBROUTINE FlattenAllRows(S)
    TYPE(FastMat),POINTER :: S,P
    !
    ! Set some pointers.
    NULLIFY(P)
    !
    P => S%Next
    DO
      IF(.NOT. ASSOCIATED(P)) EXIT
        IF(ASSOCIATED(P%RowRoot)) THEN
          NULLIFY(GlobalP)
          CALL Flatten(P%RowRoot)
          !! Terminates the tail
          NULLIFY(GlobalP%Next)
        ENDIF
      P => P%Next
    ENDDO
  END SUBROUTINE FlattenAllRows
!======================================================================
  RECURSIVE SUBROUTINE Flatten(A)
    TYPE(SRST),POINTER :: A

    IF(.NOT. ASSOCIATED(GlobalP)) THEN
      GlobalP => A
    ELSE
      GlobalP%Next => A
      GlobalP => A
    ENDIF

    IF(A%L == A%R) THEN
      !! do nothing
      IF(ASSOCIATED(A%Left) .OR. ASSOCIATED(A%Right)) THEN
        STOP 'ERR in Flatten: either left or right is not null!'
      ELSE
        !! IF(MyID == 0) THEN
        !!  WRITE(*,*) 'Flatten: assertion okay!'
        !! ENDIF
      ENDIF
    ELSE
      IF(ASSOCIATED(A%Left)) THEN
        CALL Flatten(A%Left)
      ENDIF
      IF(ASSOCIATED(A%Right)) THEN
        CALL Flatten(A%Right)
      ENDIF
    ENDIF
  END SUBROUTINE Flatten
!======================================================================
!======================================================================
!======================================================================
  SUBROUTINE Redistribute_FastMat(A)
    TYPE(FastMat),POINTER :: A
    TYPE(INT_RNK2)        :: LocalDims,RemoteDims
    TYPE(INT_VECT)        :: SndBeg,SndAbsRow,SndRow,SndCol,SndBlk,SndMtx,&
                             SndEnd,RcvBeg,SF,Dest,RcvAbsRow,RcvRow,&
                             RcvCol,RcvBlk,RcvMtx,RcvEnd,SendToQ,RecvFrQ
    REAL(DOUBLE),ALLOCATABLE,DIMENSION(:)  :: SendBuffer,RecvBuffer
    REAL(DOUBLE),ALLOCATABLE,DIMENSION(:) :: DimTmArr
    INTEGER,ALLOCATABLE,DIMENSION(:) :: IntArr
    INTEGER,DIMENSION(2)  :: MyRow,AllCol
    INTEGER               :: B,E,Num,AbsRow,Row,Col,Blk,Mtx,N,I,J,K,Tag,IErr, &
                             RE,To,From,NRecvs,NSends,SendTo,RecvFr
    CHARACTER(LEN=20)      :: Sub='Redistribute_FastMat'
    REAL(DOUBLE)           :: StartTm,EndTm,TotTm
    REAL(DOUBLE)           :: DimBegTm,DimEndTm,DimTotTm
    REAL(DOUBLE)           :: AllToAllBegTm,AllToAllEndTm,AllToAllTotTm
    INTEGER,DIMENSION(MPI_STATUS_SIZE) :: Status
    INTEGER :: TotDblSent,DblSentMax,NumDblSent,LDblSendMaxSize,LDblRecvMaxSize,ActDblRecvAmt
    REAL(DOUBLE),EXTERNAL :: MondoTimer

    StartTm = MondoTimer()

    ! Messy output is unsightly and unessesary.  Protect with if defs in future...
    !    IF(MyID == 0) THEN
    !     CALL OpenASCII(OutFile,Out)
    !    WRITE(Out,*) 'MyID=0, FastMat_redistribute is entered...'
    !   CLOSE(Out,STATUS='KEEP')
    !ENDIF

    CALL FlattenAllRows(A)

    ALLOCATE(DimTmArr(0:NPrc-1))
    ALLOCATE(IntArr(0:NPrc-1))

    CALL New(LocalDims,(/3,NPrc-1/),(/1,0/))
    CALL New(RemoteDims,(/3,NPrc-1/),(/1,0/))

    AllCol=(/1,NAtoms/)
    DO N=0,NPrc-1
       LocalDims%I(:,N)=MatDimensions_1(A,(/Beg%I(N),End%I(N)/),AllCol)
    ENDDO
    LocalDims%I(:,MyId)=(/0,0,0/)

    CALL AlignNodes()
    AllToAllBegTm = MondoTimer()

    CALL MPI_ALLTOALL( LocalDims%I(1,0),3,MPI_INTEGER,RemoteDims%I(1,0),3,MPI_INTEGER,MONDO_COMM,IErr)
    AllToAllEndTm = MondoTimer()
    AllToAllTotTm = AllToAllEndTm - AllToAllBegTm
    CALL MPI_Allgather(AllToAllTotTm,1,MPI_DOUBLE_PRECISION,DimTmArr(0),1,MPI_DOUBLE_PRECISION,MONDO_COMM,IErr)
    CALL ErrChk(IErr,Sub)

    CALL New(SendToQ,NPrc-1,0)
    CALL New(RecvFrQ,NPrc-1,0)
    NRecvs=0
    RecvFrQ%I(:) = 0
    DO N=0,NPrc-1
       IF(N /= MyId .AND. RemoteDims%I(2,N) /= 0)THEN
          NRecvs = NRecvs+1
          RecvFrQ%I(N) = 1
       ENDIF
    ENDDO
    NSends=0
    SendToQ%I(:) = 0
    DO N=0,NPrc-1
       IF(N /= MyId .AND. LocalDims%I(2,N) /= 0)THEN
          NSends = NSends+1
          SendToQ%I(N) = 1
       ENDIF
    ENDDO
    ! Check for no work
    IF(NRecvs==0.AND.NSends==0)THEN
       CALL Delete(LocalDims)
       CALL Delete(RemoteDims)
       CALL AlignNodes()
       ! Messy output is unsightly and unessesary.  Protect with if defs in future...
       !       IF(MyID == ROOT) THEN
       !         CALL OpenASCII(OutFile,Out)
       !         WRITE(Out,*) 'No sending and receiving needed in FastMat_redistribute!'
       !         CLOSE(Out,STATUS='KEEP')
       !       ENDIF
       RETURN
    ENDIF

    CALL New(SndBeg,NPrc-1,0)
    CALL New(SndRow,NPrc-1,0)
    CALL New(SndAbsRow,NPrc-1,0)
    CALL New(SndCol,NPrc-1,0)
    CALL New(SndBlk,NPrc-1,0)
    CALL New(SndMtx,NPrc-1,0)
    CALL New(SndEnd,NPrc-1,0)

    LDblSendMaxSize = -1
    DO N = 0, NPrc-1
      SndBeg%I(N) = 1
      SndAbsRow%I(N) = SndBeg%I(N)+3
      SndRow%I(N) = SndAbsRow%I(N)+LocalDims%I(1,N)+1
      SndCol%I(N)=SndRow%I(N)+LocalDims%I(1,N)+1
      SndBlk%I(N)=SndCol%I(N)+LocalDims%I(2,N)+1
      SndMtx%I(N)=SndBlk%I(N)+LocalDims%I(2,N)+1
      SndEnd%I(N)=SndMtx%I(N)+LocalDims%I(3,N)+1
      LDblSendMaxSize=Max(LDblSendMaxSize,SndEnd%I(N))
    ENDDO
    ALLOCATE(SendBuffer(1:LDblSendMaxSize),STAT=MemStatus)
    CALL IncMem(MemStatus,0,LDblSendMaxSize)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    CALL New(RcvBeg,NPrc-1,0)
    CALL New(RcvRow,NPrc-1,0)
    CALL New(RcvAbsRow,NPrc-1,0)
    CALL New(RcvCol,NPrc-1,0)
    CALL New(RcvBlk,NPrc-1,0)
    CALL New(RcvMtx,NPrc-1,0)
    CALL New(RcvEnd,NPrc-1,0)

    LDblRecvMaxSize = -1
    DO N=0,NPrc-1
       RcvBeg%I(N)=1
       RcvAbsRow%I(N) = RcvBeg%I(N)+3
       RcvRow%I(N) = RcvAbsRow%I(N)+RemoteDims%I(1,N)+1
       RcvCol%I(N)=RcvRow%I(N)+RemoteDims%I(1,N)+1
       RcvBlk%I(N)=RcvCol%I(N)+RemoteDims%I(2,N)+1
       RcvMtx%I(N)=RcvBlk%I(N)+RemoteDims%I(2,N)+1
       RcvEnd%I(N)=RcvMtx%I(N)+RemoteDims%I(3,N)+1
       LDblRecvMaxSize=Max(LDblRecvMaxSize,RcvEnd%I(N))
    ENDDO
    ! Allocate contigous memory for non-blocking recieves
    RE=LDblRecvMaxSize
    ALLOCATE(RecvBuffer(1:RE),STAT=MemStatus)
    CALL IncMem(MemStatus,0,RE)

    ! sending and receiving
    CALL New(SF,NPrc-1,0)
    CALL New(Dest,NPrc-1,0)
    NumDblSent = 0
    DO I = 1, NPrc-1
      CALL AlignNodes()
      DO J = 0, NPrc-1
        SendTo = MODULO(J+I,NPrc)
        Dest%I(J) = SendTo
        SF%I(J) = 0
      ENDDO
      DO J = 0, NPrc-1
        IF(SF%I(J) == 0 .AND. SF%I(Dest%I(J)) == 0) THEN
          SF%I(J) = 1
          SF%I(Dest%I(J)) = 2
        ENDIF
      ENDDO
      SendTo = MODULO(MyID+I,NPrc)
      RecvFr = MODULO(MyID-I,NPrc)
      IF(SF%I(MyID) == 1) THEN
        IF(SendToQ%I(SendTo) == 1) THEN
          To = SendTo
          B = SndBeg%I(To)
          E = SndEnd%I(To)
          AbsRow = SndAbsRow%I(To)
          Row = SndRow%I(To)
          Col = SndCol%I(To)
          Blk = SndBlk%I(To)
          Mtx = SndMtx%I(To)
          Num = E-B+1
          Tag = MyID
          CALL PackFastMat(A,(/Beg%I(To),End%I(To)/), &
               SendBuffer(B),SendBuffer(B+1),SendBuffer(B+2),&
               SendBuffer(AbsRow:Row-1),SendBuffer(Row:Col-1),&
               SendBuffer(Col:Blk-1),SendBuffer(Blk:MTx-1),SendBuffer(Mtx:E))
          CALL Delete_FastMat1(A,(/Beg%I(To),End%I(To)/))
          NumDblSent = NumDblSent + Num
          CALL MPI_Send(SendBuffer(B),Num,MPI_DOUBLE_PRECISION,To,Tag,MONDO_COMM,IErr)
          CALL ErrChk(IErr,Sub)
        ENDIF

        IF(RecvFrQ%I(RecvFr) == 1) THEN
          From = RecvFr
          B = RcvBeg%I(From)
          E = RcvEnd%I(From)
          Num = E-B+1
          Tag = RecvFr
          CALL MPI_Recv(RecvBuffer(B),Num,MPI_DOUBLE_PRECISION,From,Tag,MONDO_COMM,Status,IErr)
          CALL MPI_Get_Count(Status,MPI_DOUBLE_PRECISION,ActDblRecvAmt,IErr)
          IF(Num /= ActDblRecvAmt) THEN
            WRITE(*,*) 'Receive later : MyID = ',MyID, ' Num = ',Num, ',ActDblRecvAmt = ',ActDblRecvAmt
            WRITE(*,*) 'ERR: Num is not the same as ActDblRecvAmt'
            STOP
          ELSE
            !! WRITE(*,*) 'Receive later : MyID = ',MyID, ', recv from ', From, ' num = ',num
          ENDIF
          B = RcvBeg%I(From)
          AbsRow = RcvAbsRow%I(From)
          Row = RcvRow%I(From)
          Col = RcvCol%I(From)
          Blk = RcvBlk%I(From)
          Mtx = RcvMtx%I(From)
          E = RcvEnd%I(From)
          CALL UnPackNSumFastMat(A,OffSt%I(MyId),RecvBuffer(B),&
               RecvBuffer(B+1),RecvBuffer(B+2),&
               RecvBuffer(AbsRow:Row-1),RecvBuffer(Row:Col-1),&
               RecvBuffer(Col:Blk-1),RecvBuffer(Blk:Mtx-1),RecvBuffer(Mtx:E))
        ENDIF
      ELSE
        IF(RecvFrQ%I(RecvFr) == 1) THEN
          From = RecvFr
          B = RcvBeg%I(From)
          E = RcvEnd%I(From)
          Num = E-B+1
          Tag = RecvFr
          CALL MPI_Recv(RecvBuffer(B),Num,MPI_DOUBLE_PRECISION,From,Tag,MONDO_COMM,Status,IErr)
          CALL MPI_Get_Count(Status,MPI_DOUBLE_PRECISION,ActDblRecvAmt,IErr)
          IF(Num /= ActDblRecvAmt) THEN
            WRITE(*,*) 'Receive first: MyID = ',MyID, ' Num = ',Num, ',ActDblRecvAmt = ',ActDblRecvAmt
            WRITE(*,*) 'ERR: Num is not the same as ActDblRecvAmt'
            STOP
          ELSE
            !! WRITE(*,*) 'Receive first : MyID = ',MyID, ', recv from ', From, ' num = ',num
          ENDIF
          B = RcvBeg%I(From)
          AbsRow = RcvAbsRow%I(From)
          Row = RcvRow%I(From)
          Col = RcvCol%I(From)
          Blk = RcvBlk%I(From)
          Mtx = RcvMtx%I(From)
          E = RcvEnd%I(From)
          CALL UnPackNSumFastMat(A,OffSt%I(MyId),RecvBuffer(B),&
               RecvBuffer(B+1),RecvBuffer(B+2),&
               RecvBuffer(AbsRow:Row-1),RecvBuffer(Row:Col-1),&
               RecvBuffer(Col:Blk-1),RecvBuffer(Blk:Mtx-1),RecvBuffer(Mtx:E))
        ENDIF

        IF(SendToQ%I(SendTo) == 1) THEN
          To = SendTo
          B = SndBeg%I(To)
          E = SndEnd%I(To)
          AbsRow = SndAbsRow%I(To)
          Row = SndRow%I(To)
          Col = SndCol%I(To)
          Blk = SndBlk%I(To)
          Mtx = SndMtx%I(To)
          Num = E-B+1
          Tag = MyID
          CALL PackFastMat(A,(/Beg%I(To),End%I(To)/), &
               SendBuffer(B),SendBuffer(B+1),SendBuffer(B+2),  &
               SendBuffer(AbsRow:Row-1),SendBuffer(Row:Col-1),&
               SendBuffer(Col:Blk-1),SendBuffer(Blk:MTx-1),SendBuffer(Mtx:E))
          CALL Delete_FastMat1(A,(/Beg%I(To),End%I(To)/))
          NumDblSent = NumDblSent + Num
          CALL MPI_Send(SendBuffer(B),Num,MPI_DOUBLE_PRECISION,To,Tag,MONDO_COMM,IErr)
          CALL ErrChk(IErr,Sub)
        ENDIF
      ENDIF

    ENDDO

    CALL MPI_Gather(NumDblSent,1,MPI_INTEGER,IntArr(0),1,MPI_INTEGER,0,MONDO_COMM,IErr)
    IF(MyID == 0) THEN
      DblSentMax = -1
      TotDblSent = 0
      DO I = 0, NPrc-1
        DblSentMax = Max(IntArr(I),DblSentMax)
        TotDblSent = TotDblSent + IntArr(I)
      ENDDO
      ! Messy output is unsightly and unessesary.  Protect with if defs in future...
      !      CALL OpenASCII(OutFile,Out)
      !      WRITE(Out,*) 'Redistribute : DblSentMax = ',DblSentMax, ', TotDblSent = ',TotDblSent
      !      CLOSE(Out,STATUS='KEEP')
    ENDIF

    DEALLOCATE(SendBuffer,STAT=MemStatus)
    DEALLOCATE(RecvBuffer,STAT=MemStatus)

    CALL DecMem(MemStatus,LDblSendMaxSize,0)
    CALL DecMem(MemStatus,LDblRecvMaxSize,0)

    CALL Delete(LocalDims)
    CALL Delete(RemoteDims)
    CALL Delete(SndBeg)
    CALL Delete(SndAbsRow)
    CALL Delete(SndRow)
    CALL Delete(SndCol)
    CALL Delete(SndBlk)
    CALL Delete(SndMtx)
    CALL Delete(SndEnd)
    CALL Delete(RcvBeg)
    CALL Delete(RcvAbsRow)
    CALL Delete(RcvRow)
    CALL Delete(RcvCol)
    CALL Delete(RcvBlk)
    CALL Delete(RcvMtx)
    CALL Delete(RcvEnd)

    CALL FlattenAllRows(A)

    CALL AlignNodes()
    EndTm = MondoTimer()
    TotTm = EndTm - StartTm
    ! Messy output is unsightly and unessesary.  Protect with if defs in future...
    !    IF(MyID == ROOT) THEN
    !       CALL OpenASCII(OutFile,Out)
    !       WRITE(Out,*) 'Total time to Redistribute_FastMat is ', TotTm
    !       CLOSE(Out,STATUS='KEEP')
    !    ENDIF
  END SUBROUTINE Redistribute_FastMat

!======================================================================
  SUBROUTINE PackFastMat(A,RowLimits,NAtms,NBlks,Non0s,  &
                         AbsRowPt,RowPt,ColPt,BlkPt,MTrix)
    TYPE(FastMat),POINTER     :: A,R
    TYPE(SRST),POINTER        :: C
    REAL(DOUBLE),DIMENSION(:) :: AbsRowPt,RowPt,ColPt,BlkPt,MTrix
    REAL(DOUBLE)              :: NAtms,NBlks,Non0s
    INTEGER,DIMENSION(2)      :: RowLimits
    INTEGER :: INAtms,Row,M,N,MN,Col,Non0BlkP,MtxBegInd
    !
    ! Set some pointers.
    NULLIFY(R,C)
    !
    R => A%Next
    NAtms=0.0D0
    INAtms = 0
    MtxBegInd = 1 !! for block 1
    Non0BlkP = 1 !! atom-atom block
    RowPt(1) = 1
    DO
      IF(.NOT. ASSOCIATED(R)) EXIT
      Row = R%Row
      IF(Row >= RowLimits(1) .AND. Row <= RowLimits(2)) THEN
        NAtms=NAtms+1.0D0
        INAtms = INAtms + 1
        AbsRowPt(INAtms) = DBLE(Row)
        M = BSiz%I(Row)
        C => R%RowRoot
        DO
          IF(.NOT. ASSOCIATED(C)) EXIT
          Col = C%L
          IF(Col == C%R) THEN
            N = BSiz%I(Col)
            MN = M*N
            MTrix(MtxBegInd:MtxBegInd+MN-1)=PACK(C%MTrix,.TRUE.)
            BlkPt(Non0BlkP) = DBLE(MtxBegInd)
            ColPt(Non0BlkP) = DBLE(Col)
            Non0BlkP = Non0BlkP + 1
            MtxBegInd = MtxBegInd + MN
          ENDIF
          RowPt(INAtms+1) = DBLE(Non0BlkP)
          C => C%Next
        ENDDO
      ENDIF
      R => R%Next
    ENDDO
    NBlks=Non0BlkP-1
    Non0s=MtxBegInd-1
  END SUBROUTINE PackFastMat

!======================================================================
  SUBROUTINE UnPackNSumFastMat(A,OffSt,NAtms,NBlks,Non0s,  &
                               AbsRowPt,RowPt,ColPt,BlkPt,MTrix,Perf_O)
    TYPE(FastMat),POINTER     :: A,P
    TYPE(SRST),POINTER  :: U
    REAL(DOUBLE)              :: NAtms,NBlks,Non0s,Op
    REAL(DOUBLE),DIMENSION(:) :: AbsRowPt,RowPt,ColPt,BlkPt,MTrix
    INTEGER                   :: AbsRow,I,J,M,N,MN,OffSt,MtxBegInd,Col, &
                                 IAtms,BegBlk,EndBlk
    TYPE(TIME),   OPTIONAL :: Perf_O
    !
    ! Set some pointers.
    NULLIFY(P,U)
    !
    Op = ZERO
    IAtms = NAtms
    DO I = 1, IAtms
      AbsRow = AbsRowPt(I)
      M=BSiz%I(AbsRow)
      P=>FindFASTMATRow_1(A,AbsRow,SoftFind_O=.FALSE.)
      BegBlk=RowPt(I)
      EndBlk=RowPt(I+1)-1
      DO J = BegBlk, EndBlk
        Col=ColPt(J)
        MtxBegInd=BlkPt(J)
        N=BSiz%I(Col)
        MN = M*N
        U=>InsertSRSTNode(P%RowRoot,Col)
        IF(ASSOCIATED(U%MTrix))THEN
          ! Resum the block
          U%MTrix=U%MTrix+RESHAPE(MTrix(MtxBegInd:MtxBegInd+MN-1),(/M,N/))
          Op=Op+MN
        ELSE
          ! Initialize the block
          ALLOCATE(U%MTrix(M,N),STAT=MemStatus)
          CALL IncMem(MemStatus,0,MN,'AddFASTMATBlok')
          U%MTrix=RESHAPE(MTrix(MtxBegInd:MtxBegInd+MN-1),(/M,N/))
        ENDIF
      ENDDO
    ENDDO
    IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+Op
  END SUBROUTINE UnPackNSumFastMat
!=================================================================
!   ADD A BLOCK TO THE FAST MATRIX DATA STRUCTURE
!=================================================================
    SUBROUTINE AddFASTMATBlok(A,Row,Col,M,N,B)
      TYPE(FASTMAT),POINTER       :: A,C,D
      REAL(DOUBLE),DIMENSION(:,:) :: B
      TYPE(SRST), POINTER         :: P,Q
      INTEGER                     :: Row,Col,I,J,M,N
!-----------------------------------------------------------------
      !
      ! Set some pointer.
      NULLIFY(C,D,P,Q)
      !
      ! Find the current row in the fast matrix
      C=>FindFastMatRow_1(A,Row)
      ! Init global counter
      SRSTCount=C%Nodes

      ! Find/add node(s) in this rows search tree
      P=>InsertSRSTNode(C%RowRoot,Col)
      IF(ASSOCIATED(P%MTrix))THEN
         ! Resum the block
         P%MTrix=P%MTrix+B
      ELSE
         !M=BSiz%I(Row)
         !N=BSiz%I(Col)*NSMat !<<<SPIN
         ALLOCATE(P%MTrix(M,N),STAT=MemStatus)
         CALL IncMem(MemStatus,0,N*M,'AddFASTMATBlok')
         P%MTrix(1:M,1:N)=B(1:M,1:N)
      ENDIF
      ! Update counter
      C%Nodes=SRSTCount
    END SUBROUTINE AddFASTMATBlok
!=================================================================
  FUNCTION FindFastMatRow_1(A,Row,SoftFind_O) RESULT(C)
    TYPE(FASTMAT),POINTER       :: A,P,C
    LOGICAL, OPTIONAL           :: SoftFind_O
    INTEGER                     :: Row
    LOGICAL                     :: SoftFind
    !
    ! Set some pointers.
    NULLIFY(P,C)
    !
    IF(PRESENT(SoftFind_O))THEN
      SoftFind=SoftFind_O
    ELSE
      SoftFind=.FALSE.
    ENDIF

    C => A
    NULLIFY(P)
    DO
      IF(.NOT. ASSOCIATED(C)) THEN
        IF(.NOT. ASSOCIATED(P)) THEN
          WRITE(*,*) 'ERR: P is null in FindFastMatRow_1 (Logic A) !'
          STOP
        ENDIF
        IF(SoftFind) THEN
          !! do nothing, C is null anyway
        ELSE
          CALL New_FASTMAT(P%Next,Row)
          C => P%Next
        ENDIF
        RETURN
      ELSE
        IF(Row > C%Row) THEN
          P => C
          C => C%Next
        ELSEIF(Row == C%Row) THEN
          RETURN
        ELSE
          IF(SoftFind) THEN
            NULLIFY(C)
            RETURN
          ELSE
            IF(.NOT. ASSOCIATED(P)) THEN
              WRITE(*,*) 'ERR: P is null in FindFastMatRow_1 (Logic B)!'
              STOP
            ELSE
              CALL New_FASTMAT(P%Next,Row)
              P%Next%Next => C
              C => P%Next
              RETURN
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  END FUNCTION FindFastMatRow_1


    SUBROUTINE Set_FASTMAT_EQ_BCSR(B,A)
!---------------------------------------------------------------------
      TYPE(FASTMAT),POINTER :: B,C
      TYPE(SRST),POINTER    :: S
      TYPE(BCSR)            :: A
      INTEGER               :: I,J,JP,P,N,M
!---------------------------------------------------------------------
      !
      ! Set some pointers.
      NULLIFY(C,S)
      !
      !write(*,*) 'A%NAtms',A%NAtms,'A%NNon0',A%NNon0
      !write(*,*) 'A%RowPt%I',A%RowPt%I
      !write(*,*) 'A%ColPt%I',A%ColPt%I
      ! Check for prior allocation
      IF(ASSOCIATED(B))THEN
         CALL Delete_FASTMAT1(B)
         ! Begin with a new header Node
      ENDIF
      CALL New_FASTMAT(B,0,(/0,0/),NSMat_O=A%NSMat)
      !
      DO I=1,A%NAtms
         M=BSiz%I(I)
         IF(A%RowPt%I(I+1)-A%RowPt%I(I)>0)THEN
            ! Set current row link
            !write(*,*) 'I',I
            C=>FindFastMatRow_1(B,I)
            DO JP=A%RowPt%I(I),A%RowPt%I(I+1)-1
               J=A%ColPt%I(JP)
               P=A%BlkPt%I(JP)
               N=BSiz%I(J)*A%NSMat
               ! Add bloks to this sparse row search tree
               CALL AddFASTMATBlok(C,I,J,M,N,RESHAPE(A%MTrix%D(P:P+M*N-1),(/M,N/)))
               !CALL AddFASTMATBlok(C,I,J,VectToBlock(M,N,A%MTrix%D(P:P+M*N-1)))
            ENDDO
         ENDIF
      ENDDO
      !
      CALL FlattenAllRows(B)
      !
    END SUBROUTINE Set_FASTMAT_EQ_BCSR

!======================================================================
   SUBROUTINE Set_BCSR_EQ_FASTMAT(B,A)
      TYPE(FASTMAT),POINTER :: A,R
      TYPE(SRST),POINTER    :: U
      TYPE(BCSR)            :: B
      INTEGER               :: I,J,JP,P,N,M,Q,OI,OJ,IC,JC
!---------------------------------------------------------------------
      !
      ! Set some pointers.
      NULLIFY(R,U)
      !
      ! Check for prior allocation
      !write(*,*) 'In Set_BCSR_EQ_FASTMAT -1'
      !write(*,*) ASSOCIATED(A),ASSOCIATED(A%Next)
      IF(.NOT.ASSOCIATED(A%Next)) CALL Halt(' A not associated in Set_BCSR_EQ_FASTMAT')
      !
      !write(*,*) 'In Set_BCSR_EQ_FASTMAT 0'
      CALL FlattenAllRows(A)
      !
      IF(AllocQ(B%Alloc)) CALL Delete(B)
      CALL New(B)
      !
      P=1
      Q=1
      OI=0
      B%RowPt%I(1)=1
      !
      R => A%Next
      DO
         IF(.NOT.ASSOCIATED(R)) EXIT
         I = R%Row
         M = BSiz%I(I)
         OJ=0
         U => R%RowRoot
         DO
            IF(.NOT.ASSOCIATED(U)) EXIT
            IF(U%L.EQ.U%R) THEN
               IF(ASSOCIATED(U%MTrix)) THEN
                  J = U%L
                  N = BSiz%I(J)
                  DO JC=1,N
                     DO IC=1,M
                        B%MTrix%D(Q+IC-1+(JC-1)*M) = U%MTrix(IC,JC)
                     ENDDO
                  ENDDO
                  B%BlkPt%I(P)=Q
                  B%ColPt%I(P)=J
                  Q=Q+M*N
                  P=P+1
                  B%RowPt%I(I+1)=P
                  OJ=OJ+N
               ENDIF
               OI=OI+M
            ENDIF
            U => U%Next
            !write(*,*) 'In Set_BCSR_EQ_FASTMAT 7'
         ENDDO
         R => R%Next
         !write(*,*) 'In Set_BCSR_EQ_FASTMAT 8'
      ENDDO
      B%NAtms=NAtoms
      B%NBlks=P-1
      B%NNon0=Q-1
      !write(*,*) 'In Set_BCSR_EQ_FASTMAT 9'
      !
    END SUBROUTINE Set_BCSR_EQ_FASTMAT
!======================================================================
!
!======================================================================
   SUBROUTINE Set_DFASTMAT_EQ_DBCSR(B,A)
      TYPE(FASTMAT),POINTER :: B,C
      TYPE(SRST),POINTER    :: S
      TYPE(DBCSR)           :: A
      INTEGER               :: I,IRow,J,JP,P,N,M
integer::iii
!---------------------------------------------------------------------
      !
      ! Set some pointers.
      NULLIFY(C,S)
      !
      ! Check for prior allocation
      IF(ASSOCIATED(B))THEN
         CALL Delete_FASTMAT1(B)
         ! Begin with a new header Node
         CALL New_FASTMAT(B,0,(/0,0/),NSMat_O=A%NSMat)
      ENDIF
      !
      DO I = Beg%I(MyID),End%I(MyID)
         IRow = I-Beg%I(MyID)+1
         M = BSiz%I(I)
         IF(A%RowPt%I(IRow+1)-A%RowPt%I(IRow)>0) THEN
            ! Set current row link
            C => FindFastMatRow_1(B,I)
            DO JP = A%RowPt%I(IRow),A%RowPt%I(IRow+1)-1
               J = A%ColPt%I(JP)
               P = A%BlkPt%I(JP)
               N = BSiz%I(J)*A%NSMat
               ! Add bloks to this sparse row search tree
               CALL AddFASTMATBlok(C,I,J,M,N,RESHAPE(A%MTrix%D(P:P+M*N-1),(/M,N/)))
            ENDDO
         ENDIF
      ENDDO
      !
      CALL FlattenAllRows(B)
      !
    END SUBROUTINE Set_DFASTMAT_EQ_DBCSR
    !
    !
!======================================================================
    SUBROUTINE Set_DBCSR_EQ_DFASTMAT(B,A)
      TYPE(FASTMAT), POINTER :: A,R
      TYPE(SRST   ), POINTER :: U
      TYPE(DBCSR  )          :: B
      INTEGER                :: I,J,JP,P,N,M,Q,OI,OJ,IC,JC
!---------------------------------------------------------------------
      !
      ! Set some pointers.
      NULLIFY(R,U)
      !
      ! Check for prior allocation
      !write(*,*) 'In Set_BCSR_EQ_FASTMAT -1'
      IF(.NOT.ASSOCIATED(A%Next)) CALL Halt(' A not associated in Set_DBCSR_EQ_DFASTMAT')
      !
      !write(*,*) 'In Set_BCSR_EQ_FASTMAT 0'
      CALL FlattenAllRows(A)
      !
      IF(AllocQ(B%Alloc)) CALL Delete(B)
      CALL New(B)
      !
      P=1
      Q=1
      OI=0
      B%RowPt%I(1)=1
      !
      R => A%Next
      DO
         IF(.NOT.ASSOCIATED(R)) EXIT
         I = R%Row
         M = BSiz%I(I)
         OJ=0
         U => R%RowRoot
         DO
            IF(.NOT.ASSOCIATED(U)) EXIT
            IF(U%L.EQ.U%R) THEN
               IF(ASSOCIATED(U%MTrix)) THEN
                  J = U%L
                  N = BSiz%I(J)
                  DO JC=1,N                                                       ! I can create a BlockToBlock3 with that
                     DO IC=1,M                                                    !---->
                        B%MTrix%D(Q+IC-1+(JC-1)*M) = U%MTrix(IC,JC)               !
                     ENDDO                                                        !
                  ENDDO                                                           !<----
                  B%BlkPt%I(P)=Q
                  B%ColPt%I(P)=J
                  Q=Q+M*N
                  P=P+1
                  B%RowPt%I(I+1)=P
                  OJ=OJ+N
               ENDIF
               OI=OI+M
            ENDIF
            U => U%Next
            !write(*,*) 'In Set_DBCSR_EQ_DFASTMAT 7'
         ENDDO
         R => R%Next
         !write(*,*) 'In Set_DBCSR_EQ_DFASTMAT 8'
      ENDDO
      B%NAtms=NAtoms
      B%NBlks=P-1
      B%NNon0=Q-1
      !write(*,*) 'In Set_DBCSR_EQ_DFASTMAT 9'
      !
    END SUBROUTINE Set_DBCSR_EQ_DFASTMAT
!=================================================================
!   DEALLOCATE A SPARSE ROW SEARCH TREE
!=================================================================
    RECURSIVE SUBROUTINE Delete_SRST(A)
      TYPE(SRST),POINTER :: A
      INTEGER            :: MN
!-----------------------------------------------------------------
      IF(ASSOCIATED(A%Left))              &
         CALL Delete_SRST(A%Left)
      IF(ASSOCIATED(A%Right))             &
         CALL Delete_SRST(A%Right)
      IF(.NOT.ASSOCIATED(A%Left).AND.     &
         .NOT.ASSOCIATED(A%Right))THEN
         IF(ASSOCIATED(A%MTrix))THEN
            MN=SIZE(A%MTrix,1)*SIZE(A%MTrix,2)
            DEALLOCATE(A%MTrix,STAT=MemStatus)
            ! Deallocated M*N doubles
            CALL DecMem(MemStatus,0,MN)
         ENDIF
         DEALLOCATE(A,STAT=MemStatus)
         CALL DecMem(MemStatus,6,0)
         SRSTCount=SRSTCount-1
      ENDIF
    END SUBROUTINE Delete_SRST

!=================================================================
!     Allocate a new Search Tree Sparse Row Block
!=================================================================
    SUBROUTINE New_SRST(A,L,R,T)
      TYPE(SRST),POINTER :: A
      INTEGER            :: L,R,T
!-----------------------------------------------------------------
      ALLOCATE(A,STAT=MemStatus)
      ! 4 integers + 3 pointers
      CALL IncMem(MemStatus,7,0,'New_SRST')
      A%L=L
      A%R=R
      A%Tier=T
      SRSTCount=SRSTCount+1
      A%Number=SRSTCount
      NULLIFY(A%Left) !=>NULL()
      NULLIFY(A%Right)!=>NULL()
      NULLIFY(A%MTrix)!=>NULL()
    END SUBROUTINE New_SRST
!======================================================================
!     Allocate a new Search Tree Block Sparse Matrix
!======================================================================
    SUBROUTINE New_FASTMAT(A,Row,Cols_O,NSMat_O)
      TYPE(FASTMAT),POINTER         :: A
      INTEGER                       :: Row
      INTEGER,OPTIONAL,DIMENSION(2) :: Cols_O
      INTEGER,OPTIONAL              :: NSMat_O
      INTEGER,DIMENSION(2)          :: Cols
      INTEGER                       :: I
!----------------------------------------------------------------------
      IF(PRESENT(Cols_O))THEN
         Cols=Cols_O
      ELSE
         Cols=(/1,NAtoms/)
      ENDIF
      ALLOCATE(A,STAT=MemStatus)
      ! 3 integers + 2 pointers
      CALL IncMem(MemStatus,5,0,'New_FASTMAT')
      ! Initialize skip pointer flag
      A%Alloc=ALLOCATED_SKIP_POINTERS_OFF
      A%Row=Row
      A%Nodes=1
      A%NSMat=1
      IF(PRESENT(NSMat_O))A%NSMat=NSMat_O
      SRSTCount=0
      NULLIFY(A%Next)
      CALL New_SRST(A%RowRoot,Cols(1),Cols(2),0)
    END SUBROUTINE New_FASTMAT

!=================================================================
!   FIND OR ADD A COLUMN BLOCK INTO A SPARSE ROW SEARCH TREE
!=================================================================
    RECURSIVE FUNCTION InsertSRSTNode(A,Col) RESULT(C)
      TYPE(SRST), POINTER :: A,C
      INTEGER             :: Col,Split
!---------------------------------------------------------
      !
      ! Set some pointer.
      NULLIFY(C)
      !
      IF(A%Tier==0.AND.(Col<A%L.OR.Col>A%R)) THEN
         write(*,*) 'Col',Col
         write(*,*) 'A%Tier==0',A%Tier==0
         write(*,*) 'A%L',A%L
         write(*,*) 'A%R',A%R
         CALL Halt(' Logic error in InsertSRSTNode ')
      ENDIF
!     Halt recursion if we've found a node with our column
      IF(Col==A%L.AND.Col==A%R)THEN
         C=>A
         RETURN
      ENDIF
!     Halve interval
      Split=IntervalSplit(A%L,A%R)
      IF(Col>=A%L.AND.Col<=Split)THEN
         ! Go left
         IF(.NOT.ASSOCIATED(A%Left)) &
              CALL New_SRST(A%Left,A%L,Split,A%Tier+1)
         C=>InsertSRSTNode(A%Left,Col)
      ELSE
         ! Go right
         IF(.NOT.ASSOCIATED(A%Right)) &
              CALL New_SRST(A%Right,Split+1,A%R,A%Tier+1)
         C=>InsertSRSTNode(A%Right,Col)
      ENDIF
    END FUNCTION InsertSRSTNode
!=================================================================
!     SPLIT AN INTERVAL INTO
!=================================================================
      FUNCTION IntervalSplit(L,R) RESULT(Split)
         INTEGER L,R,N,Split
         N=(R-L+1)/2
         IF(N==1)N=0
         IF(L+N>R)THEN
            Split=R
         ELSEIF(R-N<L)THEN
            Split=L
         ELSE
            Split=L+N
         ENDIF
      END FUNCTION IntervalSplit

!======================================================================
!======================================================================
    SUBROUTINE SetFASTMATPart(A)
      IMPLICIT NONE
      TYPE(FASTMAT), POINTER :: A
      TYPE(FASTMAT), POINTER :: P
      TYPE(SRST   ), POINTER :: U
      !
      !
      ! Set some pointers.
      NULLIFY(P,U)
      !
      IF(.NOT.ASSOCIATED(A)) STOP 'A is null in SetFASTMATPart, STOP.' !CALL Halt(' A is null in SetFASTMATPart, STOP. ')
      P => A%Next
      DO
         IF(.NOT.ASSOCIATED(P)) EXIT
         U => P%RowRoot
         DO
            IF(.NOT.ASSOCIATED(U)) EXIT
            IF(U%L == U%R) THEN
               U%Part = Zero
            ENDIF
            U => U%Next
         ENDDO
         P => P%Next
      ENDDO
      !
    END SUBROUTINE SetFASTMATPart

  SUBROUTINE Reduce_FASTMAT(A,AFM)
    IMPLICIT NONE
    TYPE(BCSR)              :: A
    TYPE(FASTMAT), POINTER  :: AFM
    TYPE(BCSR)              :: B,C
    INTEGER                 :: LMax,NCom,L,i,To,From,iErr,it
    REAL(DOUBLE), PARAMETER :: INVLOG2=1.44269504088896D0
    !
    ! Depth of the tree.
    LMax=CEILING(LOG(DBLE(NPrc))*INVLOG2)
    !
    CALL New_BCSR(A,OnAll_O=.TRUE.,NSMat_O=AFM%NSMat)
    CALL New_BCSR(B,OnAll_O=.TRUE.,NSMat_O=AFM%NSMat)
    CALL New_BCSR(C,OnAll_O=.TRUE.,NSMat_O=AFM%NSMat)
    CALL Set_LBCSR_EQ_DFASTMAT(A,AFM)
    !
    ! Run over the shells in the tree.
    DO L=LMax,1,-1
       ! Compute the number max of send/recive.
       NCom=2**(L-1)
       !IF(MyID==0)WRITE(*,*) 'L=',L,' NCom=',NCom
       !
       To=0
       DO i=1,NCom
          From=To+2**(LMax-l);
          !IF(MyID==0)WRITE(*,*) 'To =',To,'; From =',From,'; i =',i,MyID
          !
          ! Take into account If NPrc is not a power of 2.
          IF(From.GT.NPrc-1) THEN
             !IF(MyID==0)WRITE(*,*) 'This one is grather than NPrc, we exit.'
             EXIT
          ENDIF
          !
          IF(MyID.EQ.To) THEN
             CALL Recv_LBCSR(B,From)
             IF(    A%NSMat.EQ.1.AND.B%NSMat.EQ.1)THEN
                CALL Add_LBCSR(A,B,C,1,1,1,1)
             ELSEIF(A%NSMat.EQ.2.AND.B%NSMat.EQ.2)THEN
                CALL Add_LBCSR(A,B,C,1,1,1,2)
                CALL Add_LBCSR(A,B,C,2,2,2,2)
             ELSEIF(A%NSMat.EQ.4.AND.B%NSMat.EQ.4)THEN
                CALL Add_LBCSR(A,B,C,1,1,1,4)
                CALL Add_LBCSR(A,B,C,2,2,2,4)
                CALL Add_LBCSR(A,B,C,3,3,3,4)
                CALL Add_LBCSR(A,B,C,4,4,4,4)
             ELSE
                write(*,*) 'A%NSMat=',A%NSMat,' B%NSMat=',B%NSMat,' C%NSMat=',C%NSMat,MyID
                CALL Halt('Add_BCSR: Error with NSMat!')
             ENDIF
             CALL Set_LBCSR_EQ_LBCSR(A,C)
          ENDIF
          IF(MyID.EQ.From) CALL Send_LBCSR(A,To)
          To=To+INT(2**(LMax-L+1))
       ENDDO
    ENDDO
    !
    CALL Delete_BCSR(B,OnAll_O=.TRUE.)
    CALL Delete_BCSR(C,OnAll_O=.TRUE.)
  END SUBROUTINE Reduce_FASTMAT

  SUBROUTINE Set_LBCSR_EQ_LBCSR(B,A)
    TYPE(BCSR), INTENT(INOUT) :: A
    TYPE(BCSR), INTENT(INOUT) :: B
    IF(AllocQ(B%Alloc) .AND. &
         (B%NAtms<A%NAtms.OR.B%NBlks<A%NBlks.OR.B%NNon0<A%NNon0) )THEN
       CALL Delete(B,OnAll_O=.TRUE.)
       CALL New(B,OnAll_O=.TRUE.,NSMat_O=A%NSMat)
       B%NAtms=A%NAtms; B%NBlks=A%NBlks; B%NNon0=A%NNon0
    ELSE
       IF(.NOT.AllocQ(B%Alloc))CALL New(B,OnAll_O=.TRUE.,NSMat_O=A%NSMat)
       B%NAtms=A%NAtms; B%NBlks=A%NBlks; B%NNon0=A%NNon0
    ENDIF
    CALL INT_VECT_EQ_INT_VECT(A%NAtms+1,B%RowPt%I(1),A%RowPt%I(1))
    CALL INT_VECT_EQ_INT_VECT(A%NBlks  ,B%ColPt%I(1),A%ColPt%I(1))
    CALL INT_VECT_EQ_INT_VECT(A%NBlks  ,B%BlkPt%I(1),A%BlkPt%I(1))
    CALL DBL_VECT_EQ_DBL_VECT(A%NNon0  ,B%MTrix%D(1),A%MTrix%D(1))
    !B%RowPt%I(1:A%NAtms+1)=A%RowPt%I(1:A%NAtms+1)
    !B%ColPt%I(1:A%NBlks)  =A%ColPt%I(1:A%NBlks)
    !B%BlkPt%I(1:A%NBlks)  =A%BlkPt%I(1:A%NBlks)
    !B%MTrix%D(1:A%NNon0)  =A%MTrix%D(1:A%NNon0)
  END SUBROUTINE Set_LBCSR_EQ_LBCSR

  SUBROUTINE Add_LBCSR(A,B,C,ASMat,BSMat,CSMat,NSMat)
    IMPLICIT NONE
    TYPE(BCSR), INTENT(INOUT) :: A,B
    TYPE(BCSR), INTENT(INOUT) :: C
    INTEGER                   :: ASMat,BSMat,CSMat,NSMat,Status,CLen
    REAL(DOUBLE)              :: FlOp
    IF(.NOT.AllocQ(C%Alloc)) CALL New(C,OnAll_O=.TRUE.)
    CALL New(Flag,NAtoms)
    CALL SetEq(Flag,0)
    Flop=Zero
    Status=Add_GENERIC(ASMat,BSMat,CSMat,NSMat,                 &
         &             SIZE(C%ColPt%I,1),SIZE(C%MTrix%D,1),     &
         &             A%NAtms,0,                               &
         &             C%NAtms,C%NBlks,C%NNon0,                 &
         &             A%RowPt%I,A%ColPt%I,A%BlkPt%I,A%MTrix%D, &
         &             B%RowPt%I,B%ColPt%I,B%BlkPt%I,B%MTrix%D, &
         &             C%RowPt%I,C%ColPt%I,C%BlkPt%I,C%MTrix%D, &
         &             BSiz%I,Flag%I,Flop)
    CALL Delete(Flag)
    IF(Status==FAIL) THEN
       WRITE(*,*) "Status = ",Status
       WRITE(*,*) A%NAtms,A%NBlks,A%NNon0,SIZE(A%MTrix%D)
       WRITE(*,*) B%NAtms,B%NBlks,B%NNon0,SIZE(B%MTrix%D)
       WRITE(*,*) C%NAtms,C%NBlks,C%NNon0,SIZE(C%MTrix%D)
       CALL Halt('Dimensions in Add_LBCSR')
    ENDIF
  END SUBROUTINE Add_LBCSR

  SUBROUTINE Send_LBCSR(A,To)
    IMPLICIT NONE
    TYPE(BCSR)     :: A
    TYPE(INT_VECT) :: V
    INTEGER        :: To,IErr,N,DUM(3)
    !
    DUM(1)=A%NAtms
    DUM(2)=A%NBlks
    DUM(3)=A%NNon0
    write(*,'(A,3I10,A,I3,A,I3)') 'We send',DUM(1),DUM(2),DUM(3),'  MyID',MyID,' To  ',To
    CALL MPI_SEND(DUM(1),3,MPI_INTEGER,To,10,MONDO_COMM,IErr)
    !
    N=A%NAtms+1+2*A%NBlks
    CALL New(V,N)
    CALL INT_VECT_EQ_INT_VECT(A%NAtms+1,V%I(1)                ,A%RowPt%I(1))
    CALL INT_VECT_EQ_INT_VECT(A%NBlks  ,V%I(A%NAtms+2)        ,A%ColPt%I(1))
    CALL INT_VECT_EQ_INT_VECT(A%NBlks  ,V%I(A%NAtms+2+A%NBlks),A%BlkPt%I(1))
    CALL MPI_SEND(V%I(1),N,MPI_INTEGER,To,11,MONDO_COMM,IErr)
    CALL Delete(V)
    !
    !CALL MPI_SEND(A%RowPt%I(1),A%NAtms+1,MPI_INTEGER,To,11,MONDO_COMM,IErr)
    !CALL MPI_SEND(A%ColPt%I(1),A%NBlks  ,MPI_INTEGER,To,12,MONDO_COMM,IErr)
    !CALL MPI_SEND(A%BlkPt%I(1),A%NBlks  ,MPI_INTEGER,To,13,MONDO_COMM,IErr)
    !
    CALL MPI_SEND(A%MTrix%D(1),A%NNon0  ,MPI_DOUBLE_PRECISION,To,14,MONDO_COMM,IErr)
  END SUBROUTINE Send_LBCSR

  SUBROUTINE Recv_LBCSR(A,From)
    IMPLICIT NONE
    TYPE(BCSR)     :: A
    TYPE(INT_VECT) :: V
    INTEGER        :: From,IErr,N,DUM(3),Status(MPI_STATUS_SIZE)
    ! May use that MPI_STATUS_IGNORE or MPI_STATUSES_IGNORE
    !
    CALL MPI_RECV(DUM(1),3,MPI_INTEGER,From,10,MONDO_COMM,Status,IErr)
    !CALL MPI_RECV(DUM(1),3,MPI_INTEGER,From,10,MONDO_COMM,MPI_STATUS_IGNORE,IErr)
    write(*,'(A,3I10,A,I3,A,I3)') 'We recv',DUM(1),DUM(2),DUM(3),'  MyID',MyID,' From',From
    A%NAtms=DUM(1)
    A%NBlks=DUM(2)
    A%NNon0=DUM(3)
    !
    N=A%NAtms+1+2*A%NBlks
    CALL New(V,N)
    CALL MPI_RECV(V%I(1),N,MPI_INTEGER,From,11,MONDO_COMM,Status,IErr)
    !CALL MPI_RECV(V%I(1),N,MPI_INTEGER,From,11,MONDO_COMM,MPI_STATUS_IGNORE,IErr)
    CALL INT_VECT_EQ_INT_VECT(A%NAtms+1,A%RowPt%I(1),V%I(1)                )
    CALL INT_VECT_EQ_INT_VECT(A%NBlks  ,A%ColPt%I(1),V%I(A%NAtms+2)        )
    CALL INT_VECT_EQ_INT_VECT(A%NBlks  ,A%BlkPt%I(1),V%I(A%NAtms+2+A%NBlks))
    CALL Delete(V)
    !
    !CALL MPI_RECV(A%RowPt%I(1),A%NAtms+1,MPI_INTEGER,From,11, &
    !     &        MONDO_COMM,Status,IErr)
    !CALL MPI_RECV(A%ColPt%I(1),A%NBlks  ,MPI_INTEGER,From,12, &
    !     &        MONDO_COMM,Status,IErr)
    !CALL MPI_RECV(A%BlkPt%I(1),A%NBlks  ,MPI_INTEGER,From,13, &
    !     &        MONDO_COMM,Status,IErr)
    !
    CALL MPI_RECV(A%MTrix%D(1),A%NNon0  ,MPI_DOUBLE_PRECISION, &
         &        From,14,MONDO_COMM,Status,IErr)
    !CALL MPI_RECV(A%MTrix%D(1),A%NNon0  ,MPI_DOUBLE_PRECISION, &
    !     &        From,14,MONDO_COMM,MPI_STATUS_IGNORE,IErr)
  END SUBROUTINE Recv_LBCSR

  SUBROUTINE Set_LBCSR_EQ_DFASTMAT(B,A)
    TYPE(FASTMAT),POINTER :: A,R
    TYPE(SRST),POINTER    :: U
    TYPE(BCSR)            :: B
    INTEGER               :: I,J,JP,P,N,M,MN,Q,OI,OJ,IC,JC,At,OldR
    NULLIFY(R,U)
    ! Check for prior allocation
    !write(*,*) 'In Set_BCSR_EQ_FASTMAT -1'
    IF(.NOT.ASSOCIATED(A%Next)) CALL Halt(' A not associated in Set_BCSR_EQ_FASTMAT')
    !write(*,*) 'In Set_BCSR_EQ_FASTMAT 0',MyID
    CALL FlattenAllRows(A)
    !write(*,*) 'In Set_BCSR_EQ_FASTMAT 0.1'
    !IF(AllocQ(B%Alloc)) CALL Delete(B,OnAll_O=.TRUE.)
    !CALL New(B,OnAll_O=.TRUE.)
    !write(*,*) 'In Set_BCSR_EQ_FASTMAT 0.2'
    CALL INT_VECT_EQ_INT_SCLR(NAtoms+1,B%RowPt%I(1),-100000)
    P=1
    Q=1
    OI=0
    B%RowPt%I(1)=1
    R => A%Next
    !write(*,*) 'In Set_BCSR_EQ_FASTMAT 0.3'
    DO
       IF(.NOT.ASSOCIATED(R)) EXIT
       I = R%Row
       M = BSiz%I(I)
       OJ=0
       U => R%RowRoot
       !write(*,*) 'In Set_BCSR_EQ_FASTMAT 0.4'
       DO
          IF(.NOT.ASSOCIATED(U)) EXIT
          IF(U%L.EQ.U%R) THEN
             IF(ASSOCIATED(U%MTrix)) THEN
                J = U%L
                N = BSiz%I(J)
                MN=M*N*B%NSMat !<<< SPIN
                !DO JC=1,N
                !   DO IC=1,M
                !      B%MTrix%D(Q+IC-1+(JC-1)*M) = U%MTrix(IC,JC)
                !   ENDDO
                !ENDDO
                CALL DBL_VECT_EQ_DBL_VECT(MN,B%MTrix%D(Q),U%MTrix(1,1))
                B%BlkPt%I(P)=Q
                B%ColPt%I(P)=J
                Q=Q+MN
                P=P+1
                B%RowPt%I(I+1)=P
                OJ=OJ+N
             ENDIF
             OI=OI+M
          ENDIF
          U => U%Next
          !write(*,*) 'In Set_BCSR_EQ_FASTMAT 7'
       ENDDO
       R => R%Next
       !write(*,*) 'In Set_BCSR_EQ_FASTMAT 8'
    ENDDO
    B%NAtms=NAtoms
    B%NBlks=P-1
    B%NNon0=Q-1
    !write(*,*) 'In Set_BCSR_EQ_FASTMAT 9',MyID
    OldR=B%RowPt%I(1)
    DO At=2,NAtoms+1
       IF(B%RowPt%I(At).EQ.-100000) B%RowPt%I(At)=OldR
       OldR=B%RowPt%I(At)
    ENDDO
  END SUBROUTINE Set_LBCSR_EQ_DFASTMAT
#endif
   END MODULE FASTMATRICES


#ifdef BROCKEN_CODE
!======================================================================
!!$  SUBROUTINE FlattenAllRows1(S)
!!$    TYPE(FastMat),POINTER :: S,P
!!$    !
!!$    ! Set some pointers.
!!$    NULLIFY(P)
!!$    !
!!$    P => S%Next
!!$    DO
!!$      IF(.NOT. ASSOCIATED(P)) EXIT
!!$        IF(ASSOCIATED(P%RowRoot)) THEN
!!$          NULLIFY(GlobalP1)
!!$          CALL Flatten1(P%RowRoot)
!!$          !! Terminates the tail
!!$          NULLIFY(GlobalP1%Next)
!!$        ENDIF
!!$      P => P%Next
!!$    ENDDO
!!$  END SUBROUTINE FlattenAllRows1
!!$  RECURSIVE SUBROUTINE Flatten1(A)
!!$    TYPE(SRST),POINTER :: A
!!$    IF(.NOT. ASSOCIATED(GlobalP1)) THEN
!!$      GlobalP1 => A
!!$    ELSE
!!$      GlobalP1%Next => A
!!$      GlobalP1 => A
!!$    ENDIF
!!$    IF(A%L == A%R) THEN
!!$      !! do nothing
!!$      IF(ASSOCIATED(A%Left) .OR. ASSOCIATED(A%Right)) THEN
!!$        STOP 'ERR in Flatten: either left or right is not null!'
!!$      ELSE
!!$        !! IF(MyID == 0) THEN
!!$        !!  WRITE(*,*) 'Flatten: assertion okay!'
!!$        !! ENDIF
!!$      ENDIF
!!$    ELSE
!!$      IF(ASSOCIATED(A%Left)) THEN
!!$        CALL Flatten1(A%Left)
!!$      ENDIF
!!$      IF(ASSOCIATED(A%Right)) THEN
!!$        CALL Flatten1(A%Right)
!!$      ENDIF
!!$    ENDIF
!!$  END SUBROUTINE Flatten1
!!$!======================================================================
!!$!======================================================================
!!$    SUBROUTINE Delete_FASTMAT_Node(A)
!!$!H---------------------------------------------------------------------------------
!!$!H SUBROUTINE Delete_FASTMAT_Node(A)
!!$!H
!!$!H---------------------------------------------------------------------------------
!!$      IMPLICIT NONE
!!$      TYPE(FASTMAT), POINTER :: A
!!$      IF(ASSOCIATED(A%RowRoot)) CALL Delete_SRST_1(A%RowRoot)
!!$      DEALLOCATE(A)
!!$    END SUBROUTINE Delete_FASTMAT_Node
!!$!======================================================================
!!$  SUBROUTINE PrintAllLinearRows(S)
!!$    TYPE(FastMat),POINTER :: S,P
!!$    TYPE(SRST),POINTER :: C
!!$    INTEGER :: Row
!!$    !
!!$    ! Set some pointers.
!!$    NULLIFY(P,C)
!!$    !
!!$    P => S%Next
!!$    DO
!!$      IF(.NOT. ASSOCIATED(P)) EXIT
!!$      !!  P is a valid row head
!!$      Row = P%Row
!!$      C => P%RowRoot
!!$      DO
!!$        IF(.NOT. ASSOCIATED(C)) EXIT
!!$        IF(C%L == C%R) THEN
!!$          !IF(MyID == 0) THEN
!!$             !!WRITE(*,*) 'Row = ',Row, ' Col = ', C%L
!!$          !ENDIF
!!$          IF(ASSOCIATED(C%Left) .OR. ASSOCIATED(C%Right)) THEN
!!$            STOP 'ERR in PrintAllLinearRows.. Left or Right is not null!'
!!$          ELSE
!!$            !! WRITE(*,*) 'ASSERTION in PrintAllLinearRows is okay!'
!!$          ENDIF
!!$        ENDIF
!!$        C =>  C%Next
!!$      ENDDO
!!$      P => P%Next
!!$    ENDDO
!!$  END SUBROUTINE PrintAllLinearRows
!!$!======================================================================
!!$  RECURSIVE SUBROUTINE PrintLeaf(A)
!!$    TYPE(SRST),POINTER :: A
!!$    IF(A%L == A%R) THEN
!!$      !! do nothing!
!!$    ELSE
!!$      IF(ASSOCIATED(A%Left)) THEN
!!$        CALL PrintLeaf(A%Left)
!!$      ENDIF
!!$      IF(ASSOCIATED(A%Right)) THEN
!!$        CALL PrintLeaf(A%Right)
!!$      ENDIF
!!$    ENDIF
!!$  END SUBROUTINE PrintLeaf
!!$!======================================================================
!!$  SUBROUTINE Multiply_FASTMAT_SCALAR(A,Alpha,Perf_O)
!!$!H---------------------------------------------------------------------------------
!!$!H SUBROUTINE Multiply_FASTMAT_SCALAR(A,Alpha,Perf_O)
!!$!H  This routine multilply a FASTMAT by a scalar.
!!$!H
!!$!H---------------------------------------------------------------------------------
!!$    TYPE(FASTMAT), POINTER    :: A,R
!!$    TYPE(SRST   ), POINTER    :: C
!!$    REAL(DOUBLE ), INTENT(IN) :: Alpha
!!$    TYPE(TIME   ), OPTIONAL   :: Perf_O
!!$    INTEGER                   :: Row,Col,M,N
!!$    REAL(DOUBLE )             :: Op
!!$    !
!!$    ! Set some pointers.
!!$    NULLIFY(R,C)
!!$    !
!!$    ! Build up the Skip list.
!!$    CALL FlattenAllRows(A)
!!$    !
!!$    Op=Zero                                                             !vw I should put the FlattenAllRows ?
!!$    IF(.NOT.ASSOCIATED(A)) &
!!$         & CALL Halt('ERR: A is null in Multiply_FASTMAT_SCALAR!')
!!$    R => A%Next
!!$    DO
!!$       IF(.NOT. ASSOCIATED(R)) EXIT
!!$       Row = R%Row
!!$       M = BSiz%I(Row)
!!$       C => R%RowRoot
!!$       DO
!!$          IF(.NOT. ASSOCIATED(C)) EXIT
!!$          Col = C%L
!!$          IF(Col == C%R) THEN
!!$             N = BSiz%I(Col)
!!$             Op = Op + M*N
!!$             C%MTrix = C%MTrix*Alpha
!!$          ENDIF
!!$          C => C%Next
!!$       ENDDO
!!$       R => R%Next
!!$    ENDDO
!!$    IF(PRESENT(Perf_O)) Perf_O%FLOP = Perf_O%FLOP+Op
!!$  END SUBROUTINE Multiply_FASTMAT_SCALAR
!!$  !
!!$  REAL(DOUBLE) FUNCTION Max_FASTMAT(A,Perf_O)
!!$!H---------------------------------------------------------------------------------
!!$!H FUNCTION Max_FASTMAT(A,Alpha,Perf_O)
!!$!H  This routine find the biggest element of a FASTMAT.
!!$!H
!!$!H---------------------------------------------------------------------------------
!!$    TYPE(FASTMAT), POINTER    :: A,R
!!$    TYPE(SRST   ), POINTER    :: C
!!$    TYPE(TIME   ), OPTIONAL   :: Perf_O
!!$    INTEGER                   :: Row,Col,M,N
!!$    REAL(DOUBLE )             :: Op
!!$    REAL(DOUBLE), EXTERNAL    :: DDOT
!!$    !
!!$    ! Set some pointers.
!!$    NULLIFY(R,C)
!!$    !
!!$    ! Build up the Skip list.
!!$    CALL FlattenAllRows(A)
!!$    !
!!$    Max_FASTMAT=Zero
!!$    Op=Zero                                                             !vw I should put the FlattenAllRows ?
!!$    IF(.NOT.ASSOCIATED(A)) &
!!$         & CALL Halt('ERR: A is null in Max_FASTMAT!')
!!$    R => A%Next
!!$    DO
!!$       IF(.NOT. ASSOCIATED(R)) EXIT
!!$       Row = R%Row
!!$       M = BSiz%I(Row)
!!$       C => R%RowRoot
!!$       DO
!!$          IF(.NOT. ASSOCIATED(C)) EXIT
!!$          Col = C%L
!!$          IF(Col == C%R) THEN
!!$             N = BSiz%I(Col)
!!$             Op = Op + M*N
!!$             Max_FASTMAT = MAX(Max_FASTMAT,SQRT(DDOT(M*N,C%MTrix(1,1),1,C%MTrix(1,1),1)))
!!$          ENDIF
!!$          C => C%Next
!!$       ENDDO
!!$       R => R%Next
!!$    ENDDO
!!$    IF(PRESENT(Perf_O)) Perf_O%FLOP = Perf_O%FLOP+Op
!!$  END FUNCTION Max_FASTMAT
!!$
!!$  REAL(DOUBLE) FUNCTION FNorm_FASTMAT(A,Perf_O)
!!$!H---------------------------------------------------------------------------------
!!$!H FUNCTION FNorm_FASTMAT(A,Alpha,Perf_O)
!!$!H  This routine computes the Frobenus norm of a FASTMAT.
!!$!H
!!$!H---------------------------------------------------------------------------------
!!$    TYPE(FASTMAT), POINTER    :: A,R
!!$    TYPE(SRST   ), POINTER    :: C
!!$    TYPE(TIME   ), OPTIONAL   :: Perf_O
!!$    INTEGER                   :: Row,Col,M,N
!!$    REAL(DOUBLE )             :: Op
!!$    REAL(DOUBLE ), EXTERNAL   :: DDOT
!!$    !
!!$    ! Set some pointers.
!!$    NULLIFY(R,C)
!!$    !
!!$    ! Build up the Skip list.
!!$    CALL FlattenAllRows(A)
!!$    !
!!$    FNorm_FASTMAT=Zero
!!$    Op=Zero                                                             !vw I should put the FlattenAllRows ?
!!$    IF(.NOT.ASSOCIATED(A)) &
!!$         & CALL Halt('ERR: A is null in FNorm_FASTMAT!')
!!$    R => A%Next
!!$    DO
!!$       IF(.NOT. ASSOCIATED(R)) EXIT
!!$       Row = R%Row
!!$       M = BSiz%I(Row)
!!$       C => R%RowRoot
!!$       DO
!!$          IF(.NOT. ASSOCIATED(C)) EXIT
!!$          Col = C%L
!!$          IF(Col == C%R) THEN
!!$             N = BSiz%I(Col)
!!$             Op = Op + 2*M*N
!!$             FNorm_FASTMAT = FNorm_FASTMAT+DDOT(M*N,C%MTrix(1,1),1,C%MTrix(1,1),1)
!!$          ENDIF
!!$          C => C%Next
!!$       ENDDO
!!$       R => R%Next
!!$    ENDDO
!!$    FNorm_FASTMAT=SQRT(FNorm_FASTMAT)
!!$    IF(PRESENT(Perf_O)) Perf_O%FLOP = Perf_O%FLOP+Op
!!$  END FUNCTION FNorm_FASTMAT
!!$!======================================================================
!!$!   ADD TWO DIFFERENT FAST MATRICES TOGETHER TO YEILD A THIRD:
!!$!   NEED THIS ROUTINE TO MAINTAIN OO EQUIVALENCE WITH BCSR CALLS
!!$!======================================================================
!!$    SUBROUTINE Add_FASTMAT(A,B,C,Perf_O)
!!$      TYPE(FASTMAT),POINTER  :: A,B,C
!!$      TYPE(TIME),   OPTIONAL :: Perf_O
!!$      IF(.NOT.ASSOCIATED(A))  &
!!$           CALL Halt(' A not associated in Add_FASTMAT ')
!!$      IF(.NOT.ASSOCIATED(B))  &
!!$           CALL Halt(' B not associated in Add_FASTMAT ')
!!$      IF(.NOT.ASSOCIATED(C)) &
!!$           CALL New_FASTMAT(C,0,(/0,0/))
!!$      CALL AddIn_FASTMAT(A,C,Perf_O=Perf_O)
!!$      CALL AddIn_FASTMAT(B,C,Perf_O=Perf_O)
!!$    END SUBROUTINE Add_FASTMAT
!!$!======================================================================
!!$!======================================================================
!!$    SUBROUTINE AddIn_FASTMAT(A,B,Alpha_O,Perf_O)
!!$!H---------------------------------------------------------------------------------
!!$!H SUBROUTINE AddIn_FASTMAT(A,B,Alpha_O,Perf_O)
!!$!H  This routine add two fast matrices together: A=Alpha*A+B.
!!$!H
!!$!H Comments:
!!$!H  - The B fast matrix IS FLATTENED in this routine!
!!$!H  - vw removed some bugs (See after).
!!$!H---------------------------------------------------------------------------------
!!$      TYPE(FASTMAT),POINTER  :: A,B
!!$      TYPE(FASTMAT),POINTER  :: P,Q
!!$      TYPE(SRST),   POINTER  :: U,V
!!$      REAL(DOUBLE), OPTIONAL :: Alpha_O
!!$      REAL(DOUBLE)           :: Alpha,Op
!!$      TYPE(TIME),   OPTIONAL :: Perf_O
!!$      INTEGER                :: Row,Col,M,N,MN
!!$!---------------------------------------------------------------------
!!$    !
!!$    ! Set some pointers.
!!$    NULLIFY(P,Q,U,V)
!!$    !
!!$      IF(PRESENT(Alpha_O))THEN
!!$         Alpha=Alpha_O
!!$      ELSE
!!$         Alpha=One
!!$      ENDIF
!!$      Op=Zero
!!$!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$      ! Skip pointers of A should be off since its getting added to
!!$      !vw CALL SkipsOffQ(A%Alloc,'AddIn_Fastmat :: A ')
!!$      ! B is walked upon, turn on skip pointers
!!$      !vw CALL SkipsOffQ(B%Alloc,'AddIn_FASTMAT : B')
!!$      ! +++++ Build up the Skip list +++++
!!$      CALL FlattenAllRows(B)
!!$!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$      Q=>B%Next
!!$      ! Go over rows of B
!!$      DO WHILE(ASSOCIATED(Q))
!!$         ! V is a Sparse Row ST belonging to B
!!$         V=>Q%RowRoot
!!$         Row=Q%Row
!!$         M=BSiz%I(Row)
!!$         ! Add/find P, a Sparse Row ST corresponding to Row of A
!!$         P=>FindFastMatRow_1(A,Row)               !vw there was a bug here B -> A
!!$         ! Initialize column nodes counter for Row of A
!!$         SRSTCount=P%Nodes
!!$         ! Go over columns of B
!!$         DO
!!$            IF(.NOT.ASSOCIATED(V)) EXIT
!!$            ! Check for leaf nodes
!!$            IF(V%L==V%R)THEN
!!$               Col=V%L
!!$               N=BSiz%I(Col)
!!$               ! Find/add Col node in this Rows search tree
!!$               U=>InsertSRSTNode(P%RowRoot,Col)
!!$               IF(ASSOCIATED(U%MTrix))THEN
!!$                  ! Resum the block
!!$                  U%MTrix=Alpha*U%MTrix+V%MTrix
!!$                  Op=Op+M*N
!!$               ELSE
!!$                  ! Initialize the block
!!$                  ALLOCATE(U%MTrix(M,N),STAT=MemStatus)
!!$                  CALL IncMem(MemStatus,0,N*M,'AddFASTMATBlok')
!!$                  U%MTrix=V%MTrix
!!$               ENDIF
!!$            ENDIF
!!$            !IF(ASSOCIATED(V%Left))THEN            !vw problems come from here
!!$            !    V=>V%Left                         !vw problems come from here
!!$            !ELSEIF(ASSOCIATED(V%Right))THEN       !vw problems come from here
!!$            !    V=>V%Right                        !vw problems come from here
!!$            !ELSE                                  !vw problems come from here
!!$            !    EXIT                              !vw problems come from here
!!$            !ENDIF                                 !vw problems come from here
!!$            V=>V%Next
!!$         ENDDO
!!$         ! Update counter
!!$         P%Nodes=SRSTCount
!!$         Q=>Q%Next
!!$      ENDDO
!!$!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$      ! +++++ Build up the Skip list +++++
!!$      CALL FlattenAllRows(A)
!!$!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$      !vw CALL SkipsOffFASTMAT(B)                 !vw We do not use this any more.
!!$      IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+Two*Op
!!$    END SUBROUTINE AddIn_FASTMAT
!!$!======================================================================
!!$!======================================================================
!!$    RECURSIVE FUNCTION FilterOut_SRST(A,Tol) RESULT(Op)
!!$      TYPE(SRST),POINTER  :: A
!!$      REAL(DOUBLE)        :: Tol,Op,FNorm
!!$      INTEGER             :: MN
!!$!---------------------------------------------------------------------
!!$      IF(ASSOCIATED(A%Left )) Op = Op+FilterOut_SRST(A%Left,Tol)
!!$      IF(ASSOCIATED(A%Right)) Op = Op+FilterOut_SRST(A%Right,Tol)
!!$      ! Found a leaf node?
!!$      IF(.NOT.ASSOCIATED(A%Left).AND..NOT.ASSOCIATED(A%Right))THEN
!!$         ! Check for leaf nodes with blocks
!!$         IF(ASSOCIATED(A%MTrix))THEN
!!$            ! Compute Frobinious norm
!!$            FNorm = SQRT(DOT_PRODUCT(PACK(A%MTrix,.TRUE.),PACK(A%MTrix,.TRUE.)))
!!$            MN = SIZE(A%MTrix,1)*SIZE(A%MTrix,2)
!!$            Op = MN+6 ! Approximately 6 flops for the SQRT
!!$            IF(FNorm<Tol)THEN
!!$               ! Block is under threshold, delete the parent SRST node
!!$               DEALLOCATE(A%MTrix,STAT=MemStatus)
!!$               ! Deallocated M*N doubles
!!$               CALL DecMem(MemStatus,0,MN)
!!$               DEALLOCATE(A,STAT=MemStatus)
!!$               CALL DecMem(MemStatus,6,0)
!!$               ! Decrement SRST counter
!!$               SRSTCount = SRSTCount-1
!!$            ENDIF
!!$         ELSE
!!$            ! This must be a dangling node, and we are now backtracking,
!!$            ! deleting links leading to the matrix containg node.
!!$            DEALLOCATE(A,STAT=MemStatus)
!!$            CALL DecMem(MemStatus,6,0)
!!$            SRSTCount = SRSTCount-1
!!$         ENDIF
!!$      ENDIF
!!$    END FUNCTION FilterOut_SRST
!!$!======================================================================
!!$    SUBROUTINE PChkSum_FASTMAT(A,Name,Unit_O,Proc_O,ChkInPar_O)
!!$!H---------------------------------------------------------------------------------
!!$!H SUBROUTINE  PChkSum_FASTMAT(A,Name,Unit_O,Proc_O,ChkInPar_O)
!!$!H
!!$!H---------------------------------------------------------------------------------
!!$      TYPE(FASTMAT)    , POINTER              :: A
!!$      CHARACTER(LEN=*) , INTENT(IN)           :: Name
!!$      CHARACTER(LEN=*) , INTENT(IN), OPTIONAL :: Proc_O
!!$      INTEGER          , INTENT(IN), OPTIONAL :: Unit_O
!!$      LOGICAL          , INTENT(IN), OPTIONAL :: ChkInPar_O
!!$      !-------------------------------------------------------------------
!!$      TYPE(FASTMAT)    , POINTER              :: R
!!$      TYPE(SRST   )    , POINTER              :: U
!!$      INTEGER                                 :: I,PU,J,M,N
!!$      REAL(DOUBLE)                            :: Chk
!!$      LOGICAL                                 :: InPara
!!$      CHARACTER(LEN=DEFAULT_CHR_LEN)       :: ChkStr
!!$      !-------------------------------------------------------------------
!!$      REAL(DOUBLE), EXTERNAL               :: DDOT
!!$      !-------------------------------------------------------------------
!!$      !
!!$      ! IF(PrintFlags%Key/=DEBUG_MAXIMUM.AND. &
!!$      !  PrintFlags%Chk/=DEBUG_CHKSUMS)RETURN
!!$      !
!!$      NULLIFY(R,U)
!!$      Chk=Zero
!!$      !
!!$      IF(PRESENT(ChkInPar_O)) THEN
!!$         InPara=ChkInPar_O
!!$      ELSE
!!$         InPara=.FALSE.
!!$      ENDIF
!!$      ! Flatten A.
!!$      CALL FlattenAllRows(A)
!!$      !
!!$      R => A%Next
!!$      DO
!!$         IF(.NOT.ASSOCIATED(R)) EXIT
!!$         I = R%Row
!!$         M = BSiz%I(I)
!!$         U => R%RowRoot
!!$         DO
!!$            IF(.NOT.ASSOCIATED(U)) EXIT
!!$            IF(U%L.EQ.U%R) THEN
!!$               IF(ASSOCIATED(U%MTrix)) THEN
!!$                  J = U%L
!!$                  N = BSiz%I(J)
!!$                  Chk=Chk+DDOT(M*N,U%MTrix(1,1),1,U%MTrix(1,1),1)
!!$               ENDIF
!!$            ENDIF
!!$            U => U%Next
!!$         ENDDO
!!$         R => R%Next
!!$      ENDDO
!!$#ifdef PARALLEL
!!$      IF(InPara) Chk=Reduce(Chk)
!!$      IF(MyID==ROOT)THEN
!!$#endif
!!$         Chk=SQRT(Chk)
!!$         ChkStr=CheckSumString(Chk,Name,Proc_O)
!!$         ! Write check string
!!$         ! PU=OpenPU(Unit_O=Unit_O)
!!$         ! WRITE(PU,'(1x,A)')TRIM(ChkStr)
!!$         WRITE(*,'(1x,A)')TRIM(ChkStr)
!!$         ! CALL ClosePU(PU)
!!$#ifdef PARALLEL
!!$      ENDIF
!!$#endif
!!$      !
!!$    END SUBROUTINE PChkSum_FASTMAT
!!$!=================================================================
!!$!
!!$!=================================================================
!!$    SUBROUTINE Delete_FASTMAT(A,RowLimits_O)
!!$      TYPE(FASTMAT),POINTER         :: A,B,C
!!$      INTEGER,OPTIONAL,DIMENSION(2) :: RowLimits_O
!!$      INTEGER,         DIMENSION(2) :: RowLimits
!!$!-----------------------------------------------------------------
!!$      STOP 'ERR: This is the old Delete_FASTMAT. Do not use this!'
!!$    END SUBROUTINE Delete_FASTMAT
!!$!======================================================================
!!$!  ADD A SCALAR TO A FAST MATRIX
!!$!======================================================================
!!$    SUBROUTINE Add_FASTMAT_SCLR(A,Alpha,Perf_O)
!!$      TYPE(FASTMAT),POINTER  :: A,P
!!$      TYPE(SRST),   POINTER  :: U
!!$      REAL(DOUBLE)           :: Alpha,Op
!!$      TYPE(TIME),   OPTIONAL :: Perf_O
!!$      INTEGER                :: M
!!$!---------------------------------------------------------------------
!!$      !
!!$      ! Set some pointers.
!!$      NULLIFY(P,U)
!!$      !
!!$      IF(.NOT.ASSOCIATED(A))  &
!!$           CALL Halt(' A not associated in Add_FASTMAT_SCLR ')
!!$      !
!!$      ! Build up the Skip list.
!!$      CALL FlattenAllRows(A)
!!$      !
!!$      ! CALL SkipsOffQ(A%Alloc,'Add_Fastmat_SCLR :: A ')
!!$      ! CALL SkipsOnFASTMAT(A)
!!$      Op=Zero
!!$!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$      P=>A%Next
!!$      DO
!!$         IF(.NOT.ASSOCIATED(P)) EXIT
!!$         U=>P%RowRoot
!!$         M=BSiz%I(P%Row)
!!$         DO
!!$            IF(.NOT.ASSOCIATED(U)) EXIT
!!$            IF(U%L==U%R)THEN
!!$               U%MTrix=U%MTrix+Alpha
!!$               Op=Op+M*BSiz%I(U%L)
!!$            ENDIF
!!$            U=>U%Next
!!$            !IF(ASSOCIATED(U%Left))THEN
!!$            !   U=>U%Left
!!$            !ELSEIF(ASSOCIATED(U%Right))THEN
!!$            !   U=>U%Right
!!$            !ELSE
!!$            !   EXIT
!!$            !ENDIF
!!$         ENDDO
!!$         P=>P%Next
!!$      ENDDO
!!$!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$      IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+Op
!!$   END SUBROUTINE Add_FASTMAT_SCLR
!!$!======================================================================
!!$!  MULTIPLY TWO FAST MATRICES TOGETHER: C=Alpha*A*B+Beta*C
!!$!======================================================================
!!$    SUBROUTINE Multiply_FASTMAT(A,B,C,Alpha_O,Beta_O,Perf_O)
!!$      TYPE(FASTMAT),POINTER  :: A,B,C
!!$      TYPE(FASTMAT),POINTER  :: P,Q,R
!!$      TYPE(SRST),   POINTER  :: U,V,W
!!$      REAL(DOUBLE), OPTIONAL :: Alpha_O,Beta_O
!!$      REAL(DOUBLE)           :: Alpha,Beta,Op
!!$      TYPE(TIME),   OPTIONAL :: Perf_O
!!$      INTEGER                :: Row,Col,K,MA,MB,NB,MN
!!$!----------------------------------------------------------------------
!!$      !
!!$      ! Set some pointers.
!!$      NULLIFY(P,Q,R,U,V,W)
!!$      !
!!$
!!$      IF(PRESENT(Alpha_O))THEN
!!$         Alpha=Alpha_O
!!$      ELSE
!!$         Alpha=One
!!$      ENDIF
!!$      IF(PRESENT(Beta_O))THEN
!!$         Beta=Beta_O
!!$      ELSE
!!$         Beta=One
!!$      ENDIF
!!$      Op=Zero
!!$!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$      IF(.NOT.ASSOCIATED(A))  &
!!$           CALL Halt(' A not associated in Multiply_FASTMAT ')
!!$      IF(.NOT.ASSOCIATED(B))  &
!!$           CALL Halt(' B not associated in Multiply_FASTMAT ')
!!$      IF(.NOT.ASSOCIATED(C)) &
!!$           CALL New_FASTMAT(C,0,(/0,0/))
!!$      !
!!$      ! Build up the Skip lists.
!!$      CALL FlattenAllRows(A)
!!$      CALL FlattenAllRows1(B)
!!$      !
!!$      ! Skip pointers of A and B are set on for column traversal
!!$      !CALL SkipsOffQ(A%Alloc,'Multiply_FASTMAT : A')
!!$      !CALL SkipsOffQ(B%Alloc,'Multiply_FASTMAT : B')
!!$      !CALL SkipsOnFASTMAT(A)
!!$      !CALL SkipsOnFASTMAT(B)
!!$      !Skip pointers of C should be off as its getting resummed
!!$      !CALL SkipsOffQ(C%Alloc,'Multiply_FASTMAT : C')
!!$!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$      P=>A%Next
!!$      ! Go over rows of A
!!$      DO
!!$         IF(.NOT.ASSOCIATED(P)) EXIT
!!$         ! U is the Sparse Row ST containing A
!!$         U=>P%RowRoot
!!$         Row=P%Row
!!$         MA=BSiz%I(Row)
!!$         ! R is the Sparse Row ST containing
!!$         ! cooresponding coloumns of C
!!$         R=>FindFastMatRow_1(C,Row)
!!$         ! Initialize column nodes counter for Row of C
!!$         SRSTCount=R%Nodes
!!$         ! Go over columns of A
!!$         DO
!!$            IF(.NOT.ASSOCIATED(U)) EXIT
!!$            ! Check for leaf nodes
!!$            IF(U%L==U%R)THEN
!!$               ! K is the kontraction index: cols of A, rows of B
!!$               K=U%L
!!$               ! Perform a hard find (no add) of row K in B,
!!$               ! proceed only if it exists
!!$               Q=>FindFastMatRow_1(B,K,SoftFind_O=.TRUE.)
!!$               IF(ASSOCIATED(Q))THEN
!!$                  MB=BSiz%I(K)
!!$                  MN=MA*MB
!!$                  ! V is the Sparse Row ST containing the columns of B
!!$                  V=>Q%RowRoot
!!$                  DO
!!$                     IF(.NOT.ASSOCIATED(V)) EXIT
!!$                     ! Check for leaf node
!!$                     IF(V%L==V%R)THEN
!!$                        Col=V%L
!!$                        NB=BSiz%I(Col)
!!$                        ! Fast lookup/add of W=>C(A%Row,B%Col)
!!$                        W=>InsertSRSTNode(R%RowRoot,Col)
!!$                        IF(ASSOCIATED(W%MTrix))THEN
!!$                           ! Resum the block
!!$                           CALL DGEMM_NNC(MA,MB,NB,Alpha,Beta,U%MTrix(1,1),V%MTrix(1,1),W%MTrix(1,1))
!!$                        ELSE
!!$                           ! Initialize this block
!!$                           ALLOCATE(W%MTrix(MA,NB),STAT=MemStatus)
!!$                           CALL IncMem(MemStatus,0,MA*NB,'AddFASTMATBlok')
!!$                           CALL DGEMM_NNC(MA,MB,NB,Alpha,Zero,U%MTrix(1,1),V%MTrix(1,1),W%MTrix(1,1))
!!$                        ENDIF
!!$                        Op=Op+DBLE(MN*NB)
!!$                     ENDIF
!!$                     ! Tracing skip pointers over V, the coloumn index of B
!!$                     V=>V%Next
!!$                     !IF(ASSOCIATED(V%Left))THEN
!!$                     !   V=>V%Left
!!$                     !ELSEIF(ASSOCIATED(V%Right))THEN
!!$                     !   V=>V%Right
!!$                     !ELSE
!!$                     !   EXIT
!!$                     !ENDIF
!!$                  ENDDO
!!$               ENDIF
!!$            ENDIF
!!$            ! Tracing skip pointers over U, the kontraction index K
!!$            U=>U%Next
!!$            !IF(ASSOCIATED(U%Left))THEN          !vw this does not work anymore!
!!$            !   U=>U%Left                        !vw this does not work anymore!
!!$            !ELSEIF(ASSOCIATED(U%Right))THEN     !vw this does not work anymore!
!!$            !   U=>U%Right                       !vw this does not work anymore!
!!$            !ELSE                                !vw this does not work anymore!
!!$            !   EXIT                             !vw this does not work anymore!
!!$            !ENDIF                               !vw this does not work anymore!
!!$         ENDDO ! Cols of A                      !vw this does not work anymore!
!!$         ! Update C column nodes counter
!!$         R%Nodes=SRSTCount
!!$         ! Keep following As row links
!!$         P=>P%Next
!!$      ENDDO ! Rows of A
!!$!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$      !CALL SkipsOffFASTMAT(A)
!!$      !CALL SkipsOffFASTMAT(B)
!!$      IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+Two*Op
!!$    END SUBROUTINE Multiply_FASTMAT
!!$
!!$    FUNCTION Trace_FASTMAT(A,Perf_O) RESULT(Trace)
!!$!H---------------------------------------------------------------------------------
!!$!H SUBROUTINE Trace_FASTMAT(A,Perf_O)
!!$!H  This routine computes the trace of a FASTMAT.
!!$!H
!!$!H---------------------------------------------------------------------------------
!!$      TYPE(FASTMAT), POINTER  :: A,P
!!$      TYPE(TIME   ), OPTIONAL :: Perf_O
!!$      TYPE(SRST   ), POINTER  :: U
!!$      INTEGER                 :: Row,Col,M,I
!!$      REAL(DOUBLE)            :: Trace,Op
!!$      !
!!$      ! Set some pointers.
!!$      NULLIFY(P,U)
!!$      !
!!$      IF(.NOT.ASSOCIATED(A))  &
!!$           & CALL Halt(' A not associated in Trace_FASTMAT ')
!!$      !
!!$      CALL FlattenAllRows(A)
!!$      !
!!$      Trace = Zero
!!$      Op = Zero
!!$      !
!!$      P=>A%Next
!!$      DO
!!$         IF(.NOT.ASSOCIATED(P)) EXIT
!!$         Row = P%Row
!!$         M = BSiz%I(Row)
!!$         U => P%RowRoot
!!$         DO
!!$            IF(.NOT.ASSOCIATED(U)) EXIT
!!$            IF(U%L==U%R) THEN
!!$               Col = U%L
!!$               IF(Row==Col) THEN
!!$                  IF(ASSOCIATED(U%MTrix)) THEN
!!$                     !
!!$                     ! Get the trace of the sub-block.
!!$                     DO I=1,M
!!$                        Trace=Trace+U%MTrix(I,I)
!!$                     ENDDO
!!$                     Op = Op+M
!!$                  ENDIF
!!$               ENDIF
!!$            ENDIF
!!$            U => U%Next
!!$         ENDDO
!!$         P=>P%Next
!!$      ENDDO
!!$      !
!!$    END FUNCTION Trace_FASTMAT
!!$    !
!!$    FUNCTION TraceMM_FASTMAT(A,B,Perf_O) RESULT(Trace)
!!$!H---------------------------------------------------------------------------------
!!$!H SUBROUTINE TraceMM_FASTMAT(A,B,Perf_O)
!!$!H  This routine computes the trace of a product of two FASTMAT.
!!$!H
!!$!H---------------------------------------------------------------------------------
!!$      TYPE(FASTMAT), POINTER  :: A,B,P,Q
!!$      TYPE(TIME   ), OPTIONAL :: Perf_O
!!$      TYPE(SRST   ), POINTER  :: U,V
!!$      INTEGER                 :: Row,Col,M,N,Split
!!$      REAL(DOUBLE)            :: Trace,Op
!!$
!!$      REAL(DOUBLE), EXTERNAL :: BlkTrace_2
!!$      !
!!$      ! Set some pointers.
!!$      NULLIFY(P,Q,U,V)
!!$      !
!!$      !
!!$      IF(.NOT.ASSOCIATED(A))  &
!!$           & CALL Halt(' A not associated in TraceMM_FASTMAT ')
!!$      IF(.NOT.ASSOCIATED(B))  &
!!$           & CALL Halt(' B not associated in TraceMM_FASTMAT ')
!!$      !
!!$      CALL FlattenAllRows(A)
!!$      !
!!$      Trace = Zero
!!$      Op = Zero
!!$      !
!!$      P=>A%Next
!!$      DO !i
!!$         IF(.NOT.ASSOCIATED(P)) EXIT
!!$         Row = P%Row
!!$         M = BSiz%I(Row)
!!$         U => P%RowRoot
!!$         DO !j
!!$            IF(.NOT.ASSOCIATED(U)) EXIT
!!$            IF(U%L==U%R) THEN
!!$               Col = U%L
!!$               N = BSiz%I(Col)
!!$               IF(ASSOCIATED(U%MTrix)) THEN
!!$                  Q=>FindFastMatRow_1(B,Col,SoftFind_O=.TRUE.)!j
!!$                  V=>Q%RowRoot
!!$                  Tree:DO !i
!!$                     IF(Row==V%L.AND.Row==V%R) then
!!$                        IF(ASSOCIATED(V%MTrix)) THEN
!!$                           ! Get the trace of the sub-block.
!!$                           Trace=Trace+BlkTrace_2(M,N,U%MTrix(1,1),V%MTrix(1,1))
!!$                           Op = Op+M
!!$                        ENDIF
!!$                        EXIT Tree
!!$                     ENDIF
!!$                     Split=IntervalSplit(V%L,V%R)
!!$                     IF(Row>=V%L.AND.Row<=Split.AND.ASSOCIATED(V%Left))THEN
!!$                        V=>V%Left
!!$                     ELSEIF(Row>=Split.AND.Row<=V%R.AND.ASSOCIATED(V%Right))THEN
!!$                        V=>V%Right
!!$                     ELSE
!!$                        EXIT Tree
!!$                     ENDIF
!!$                  ENDDO Tree
!!$               ENDIF
!!$            ENDIF
!!$            U => U%Next
!!$         ENDDO
!!$         P=>P%Next
!!$      ENDDO
!!$      !
!!$    END FUNCTION TraceMM_FASTMAT
!!$!======================================================================
!!$!  IN PLACE FILTERATION OF SMALL BLOCKS FROM A FAST MATRIX
!!$!======================================================================
!!$    SUBROUTINE FilterOut_FASTMAT(A,Tol_O,Perf_O)
!!$!H---------------------------------------------------------------------------------
!!$!H SUBROUTINE FilterOut_FASTMAT(A,Tol_O,Perf_O)
!!$!H
!!$!H---------------------------------------------------------------------------------
!!$      TYPE(FASTMAT), POINTER  :: A,P,Q
!!$      REAL(DOUBLE ), OPTIONAL :: Tol_O
!!$      TYPE(TIME   ), OPTIONAL :: Perf_O
!!$      REAL(DOUBLE )           :: Tol,Op
!!$      TYPE(SRST   ), POINTER  :: U,V
!!$      INTEGER                 :: Row,Col,M,N,MN
!!$!----------------------------------------------------------------------
!!$      !
!!$      ! Set some pointers.
!!$      NULLIFY(P,Q,U,V)
!!$      !
!!$      IF(PRESENT(Tol_O))THEN
!!$         Tol=Tol_O
!!$      ELSE
!!$         Tol=Thresholds%Trix
!!$      ENDIF
!!$      Op=Zero
!!$!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$      ! Cols of A are walked on as an ordinary binary tree,  ! We do not need this any more
!!$      ! skip pointers should be off                          ! We do not need this any more
!!$      !vw CALL SkipsOffQ(A%Alloc,'FilterOut_FASTMAT : A')    ! We do not need this any more
!!$!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$      ! Q is the parent link of P, to allow excission of a row
!!$      Q=>A
!!$      P=>A%Next
!!$      ! Go over rows of A
!!$      DO
!!$         IF(.NOT.ASSOCIATED(P)) EXIT
!!$         ! Initialize column nodes counter for Row of A
!!$         SRSTCount = P%Nodes
!!$         ! Filter out columns for this row, incrementing FLOP count
!!$         Op = Op+FilterOut_SRST(P%RowRoot,Tol)
!!$         ! Reset SRST node counter
!!$         P%Nodes=SRSTCount
!!$         ! Check to see if the SRST is now empty for this row
!!$         IF(P%Nodes==0)THEN
!!$            ! Skip over P in the LL
!!$            Q%Next=>P%Next
!!$            ! Delete P
!!$            !CALL Delete_FASTMAT(P) !vw changed ""->1
!!$            CALL Delete_FASTMAT_Node(P)      ! vw add that
!!$            !write(*,*) 'go out delete',myid
!!$            P=>Q                             ! vw add that
!!$            P%Next=>Q%Next                   ! vw add that
!!$         ENDIF
!!$         Q=>P
!!$         P=>P%Next
!!$      ENDDO
!!$!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$      IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+Op
!!$    END SUBROUTINE FilterOut_FASTMAT
!!$!======================================================================
!!$!   SET SKIP POINTERS FOR A FAST MATRIX
!!$!======================================================================
!!$    RECURSIVE SUBROUTINE SkipsOnFASTMAT(A)
!!$      TYPE(FASTMAT),POINTER :: A,C
!!$!----------------------------------------------------------------------
!!$      ! Check to see if skip pointers are already on. This is not
!!$      ! nessesarily an error, as if already on its probably because
!!$      ! the same pointer was passed twice to a procedure as in B=A*A.
!!$      STOP 'Err: SkipsOnFASTMAT should not be called at all!'
!!$    END SUBROUTINE SkipsOnFASTMAT
!!$!======================================================================
!!$!   SET SKIP POINTERS FOR THE SPARSE ROW SEARCH TREE (SRST),
!!$!   ALLOWING IN-PROCEDURE RECURSION OVER LEAF NODES
!!$!======================================================================
!!$  RECURSIVE SUBROUTINE SetSRSTSkipPtrs(A,B)
!!$    TYPE(SRST),POINTER          :: A
!!$    TYPE(SRST),POINTER,OPTIONAL :: B
!!$    LOGICAL                     :: AssocB,AssocAL,AssocAR, &
!!$                                   AssocBL,AssocBR
!!$!----------------------------------------------------------------------
!!$    STOP 'ERR: SetSRSTSkipPtrs should not be called!'
!!$  END SUBROUTINE SetSRSTSkipPtrs
!!$!======================================================================
!!$!   QUERY TO SEE IF SKIP POINTERS ARE ON
!!$!======================================================================
!!$    SUBROUTINE SkipsOnQ(Key,Proc)
!!$       INTEGER          :: Key
!!$       CHARACTER(LEN=*) :: Proc
!!$       IF(Key==ALLOCATED_SKIP_POINTERS_ON)RETURN
!!$       IF(Key==ALLOCATED_SKIP_POINTERS_OFF)THEN
!!$          CALL Halt(' Skip pointers incorrectly set in '//TRIM(Proc))
!!$       ELSE
!!$          CALL Halt(' Uninitialized allocation key in '//TRIM(Proc))
!!$       ENDIF
!!$     END SUBROUTINE SkipsOnQ
!!$!======================================================================
!!$!   UNSET SKIP POINTERS FOR A FAST MATRIX
!!$!======================================================================
!!$    RECURSIVE SUBROUTINE SkipsOffFASTMAT(A)
!!$      TYPE(FASTMAT),POINTER :: A,C
!!$!----------------------------------------------------------------------
!!$      ! Check to see if skip pointers are already off. This is not
!!$      ! nessesarily an error, as if already on its probably because
!!$      ! the same pointer was passed twice to a procedure as in B=A*A.
!!$
!!$      STOP 'ERR: SkipsOffFASTMAT is no longer needed!'
!!$    END SUBROUTINE SkipsOffFASTMAT
!!$!======================================================================
!!$!   REMOVE SKIP POINTERS, YEILDING A PLAIN 1-D BINARY SEARCH TREE
!!$!======================================================================
!!$    SUBROUTINE UnSetSRSTSkipPtrs(A)
!!$      TYPE(SRST),POINTER :: A,B,C
!!$      !
!!$      ! Set some pointers.
!!$      NULLIFY(B,C)
!!$      !
!!$        C=>A
!!$        DO ! Recur, always following the left link
!!$           IF(ASSOCIATED(C%Left))THEN
!!$              C=>C%Left
!!$           ELSEIF(ASSOCIATED(C%Right))THEN
!!$              ! Check right node for skip pointer
!!$              IF(C%Right%Tier<=C%Tier)THEN
!!$                 B=>C%Right
!!$                 ! Remove pointer
!!$                 NULLIFY(C%Right)!=>NULL()
!!$                 C=>B
!!$              ELSE
!!$                 C=>C%Right
!!$              ENDIF
!!$           ELSE
!!$              ! All done
!!$              EXIT
!!$           ENDIF
!!$        ENDDO
!!$    END SUBROUTINE UnSetSRSTSkipPtrs
!!$!======================================================================
!!$!   QUERY TO SEE IF SKIP POINTERS ARE OFF
!!$!======================================================================
!!$    SUBROUTINE SkipsOffQ(Key,Proc)
!!$       INTEGER          :: Key
!!$       CHARACTER(LEN=*) :: Proc
!!$       IF(Key==ALLOCATED_SKIP_POINTERS_OFF)RETURN
!!$       IF(Key==ALLOCATED_SKIP_POINTERS_ON)THEN
!!$          CALL Halt(' Skip pointers incorrectly set in '//TRIM(Proc))
!!$       ELSE
!!$          CALL Halt(' Uninitialized allocation key in '//TRIM(Proc))
!!$       ENDIF
!!$     END SUBROUTINE SkipsOffQ
!!$
!!$!======================================================================
!!$!======================================================================
!!$     SUBROUTINE Symmetrized_FASMAT(A,MatFormat,Perf_O)
!!$!H---------------------------------------------------------------------------------
!!$!H SUBROUTINE Symmetrized_FASMAT(A,MatFormat,Perf_O)
!!$!H
!!$!H
!!$!H---------------------------------------------------------------------------------
!!$       IMPLICIT NONE
!!$       TYPE(FASTMAT)   , POINTER    :: A,B,C
!!$       TYPE(TIME   )   , OPTIONAL   :: Perf_O
!!$       INTEGER                      :: Col,Row,M
!!$       REAL(DOUBLE)                 :: Op
!!$       CHARACTER(LEN=1), INTENT(IN) :: MatFormat !="L" (Lower) or "U" (Upper)
!!$!                                                ! block triangular matrix.
!!$!----------------------------------------------------------------------
!!$       !
!!$       ! Set some pointers.
!!$       NULLIFY(B,C)
!!$       !
!!$       IF(.NOT.ASSOCIATED(A))  &
!!$            CALL Halt(' A not associated in Symmetrized_FASTMAT ')
!!$!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$       ! +++++ Aii=0.5*Aii and Ai>j=0 or Ai<j=0 +++++
!!$       CALL CleanDiag_FASTMAT(A,MatFormat,Alpha_O=0.5d0,Perf_O=Perf_O)
!!$       ! +++++ B=A^t +++++
!!$       B => Transpose_FASTMAT(A,Perf_O=Perf_O)
!!$       ! +++++ A=A+B +++++
!!$       CALL AddIn_FASTMAT(A,B,Perf_O=Perf_O)
!!$       ! +++++ Delete B +++++
!!$       CALL Delete_FASTMAT1(B)
!!$!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$     END SUBROUTINE Symmetrized_FASMAT
!!$!======================================================================
!!$!  TRANSPOSE (LOCALY) A FAST MATT-RIX
!!$!======================================================================
!!$     FUNCTION Transpose_FASTMAT(A,Perf_O) RESULT(B)
!!$!H---------------------------------------------------------------------------------
!!$!H FUNCTION Transpose_FASTMAT(A,Perf_O) RESULT(B)
!!$!H  This routine transpose (localy) a fast matrix.
!!$!H
!!$!H Comments:
!!$!H
!!$!H---------------------------------------------------------------------------------
!!$       IMPLICIT NONE
!!$       TYPE(FASTMAT), POINTER  :: A,P,B
!!$       TYPE(SRST   ), POINTER  :: U
!!$       TYPE(TIME   ), OPTIONAL :: Perf_O
!!$       INTEGER                 :: Col,Row,M
!!$       REAL(DOUBLE)            :: Op
!!$!----------------------------------------------------------------------
!!$    !
!!$    ! Set some pointers.
!!$    NULLIFY(P,B,U)
!!$    !
!!$       IF(.NOT.ASSOCIATED(A)) CALL Halt(' A not associated in Transpose_FASTMAT ')
!!$       CALL New_FASTMAT(B,0,(/0,0/))
!!$!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$       CALL FlattenAllRows(A)
!!$
!!$       Op = Zero
!!$       P => A%Next
!!$       DO
!!$          IF(.NOT.ASSOCIATED(P)) EXIT
!!$          Row = P%Row
!!$          M = BSiz%I(Row)
!!$          U => P%RowRoot
!!$          DO
!!$             IF(.NOT.ASSOCIATED(U)) EXIT
!!$             IF(U%L == U%R) THEN
!!$                IF(ASSOCIATED(U%MTrix)) THEN
!!$                   Col = U%L
!!$                   ! +++++ Copy the transposed block in the new FASTMAT +++++
!!$                   CALL AddFASTMATBlok(B,Col,Row,TRANSPOSE(U%MTrix(:,:)))
!!$                   Op = Op+M*BSiz%I(U%L)
!!$                ENDIF
!!$             ENDIF
!!$             U => U%Next
!!$          ENDDO
!!$          P => P%Next
!!$       ENDDO
!!$
!!$       CALL FlattenAllRows(A)
!!$!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$       IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+Op
!!$     END FUNCTION Transpose_FASTMAT
!!$!======================================================================
!!$     SUBROUTINE CleanDiag_FASTMAT(A,MatFormat,Alpha_O,Perf_O)
!!$!H---------------------------------------------------------------------------------
!!$!H SUBROUTINE CleanDiag_FASTMAT(A,MatFormat,Perf_O)
!!$!H  This routine clean the diagonal part of a fast matrix.
!!$!H  e.g.: Aii=Alpha*Aii and Ai>j=0 or Ai<j=0
!!$!H
!!$!H  MatFormat = 'L', The matrix is in LOWER triangular part -> Clean UPPER part.
!!$!H  MatFormat = 'U', The matrix is in UPPER triangular part -> Clean LOWER part.
!!$!H
!!$!H Comments:
!!$!H
!!$!H---------------------------------------------------------------------------------
!!$       IMPLICIT NONE
!!$       TYPE(FASTMAT)   , POINTER    :: A,P
!!$       TYPE(SRST   )   , POINTER    :: U
!!$       TYPE(TIME   )   , OPTIONAL   :: Perf_O
!!$       REAL(DOUBLE )   , OPTIONAL   :: Alpha_O
!!$       CHARACTER(LEN=1), INTENT(IN) :: MatFormat
!!$       INTEGER                      :: Row,M,I,J,Split
!!$       REAL(DOUBLE)                 :: Op,Alpha
!!$       LOGICAL                      :: IsLower
!!$!----------------------------------------------------------------------
!!$       !
!!$       ! Set some pointers.
!!$       NULLIFY(P,U)
!!$       !
!!$       IF(MatFormat.NE.'L'.AND.MatFormat.NE.'U')  &
!!$            CALL Halt(' Cannot recognize the MatFormat in CleanDiag_FASTMAT ')
!!$!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$       IF(PRESENT(Alpha_O)) THEN
!!$          Alpha=Alpha_O
!!$       ELSE
!!$          Alpha=1.0d0
!!$       ENDIF
!!$       IF(MatFormat.EQ.'L') IsLower=.TRUE.
!!$       IF(MatFormat.EQ.'U') IsLower=.FALSE.
!!$       Op=Zero
!!$       P=>A%Next
!!$       DO WHILE(ASSOCIATED(P))
!!$          Row=P%Row
!!$          U=>P%RowRoot
!!$          M=BSiz%I(Row)
!!$          Tree: DO
!!$             IF(Row==U%L.AND.Row==U%R) then
!!$                IF(ASSOCIATED(U%MTrix)) THEN
!!$                   IF(IsLower) THEN
!!$                      ! +++++ Clean Upper Part +++++
!!$                      U%MTrix(1,1)=U%MTrix(1,1)*Alpha
!!$                      DO J=2,M
!!$                         U%MTrix(J,J)=U%MTrix(J,J)*Alpha
!!$                         DO I=1,J-1
!!$                            U%MTrix(I,J)=Zero
!!$                         ENDDO
!!$                      ENDDO
!!$                   ELSE
!!$                      ! +++++ Clean Lower Part +++++
!!$                      DO J=1,M-1
!!$                         U%MTrix(J,J)=U%MTrix(J,J)*Alpha
!!$                         DO I=J+1,M
!!$                            U%MTrix(I,J)=Zero
!!$                         ENDDO
!!$                      ENDDO
!!$                      U%MTrix(M,M)=U%MTrix(M,M)*Alpha
!!$                   ENDIF
!!$                   Op=Op+DBLE(M*(M-1))
!!$                ENDIF
!!$                EXIT Tree
!!$             ENDIF
!!$             Split=IntervalSplit(U%L,U%R)
!!$             IF(Row>=U%L.AND.Row<=Split.AND.ASSOCIATED(U%Left))THEN
!!$                U=>U%Left
!!$             ELSEIF(Row>=Split.AND.Row<=U%R.AND.ASSOCIATED(U%Right))THEN
!!$                U=>U%Right
!!$             ELSE
!!$                EXIT Tree
!!$             ENDIF
!!$          ENDDO Tree
!!$          P=>P%Next
!!$       ENDDO
!!$!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$       IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+Op
!!$     END SUBROUTINE CleanDiag_FASTMAT
!!$
!!$!======================================================================
!!$!   FAST MATT-RIX
!!$!======================================================================
!!$     SUBROUTINE DepthTree_FASTMAT(A)
!!$       IMPLICIT NONE
!!$       TYPE(FASTMAT), POINTER :: A,P
!!$       TYPE(SRST   ), POINTER :: U
!!$       INTEGER                :: Col,Row,DepthMin,DepthMax
!!$       INTEGER                :: Tier
!!$!----------------------------------------------------------------------
!!$       !
!!$       ! Set some pointers.
!!$       NULLIFY(P,U)
!!$       !
!!$       P=>A%Next
!!$       DO WHILE(ASSOCIATED(P))
!!$          Row=P%Row
!!$          U=>P%RowRoot
!!$          DepthMin=0
!!$          DepthMax=0
!!$          call DepthTree(U,DepthMin,DepthMax)
!!$          write(*,*) "Row=",Row," DepthMin",DepthMin," DepthMax",DepthMax
!!$          P=>P%Next
!!$       ENDDO
!!$!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$     END SUBROUTINE DepthTree_FASTMAT
!!$
!!$     RECURSIVE SUBROUTINE DepthTree(U,DepthMin,DepthMax)
!!$       TYPE(SRST   ), POINTER :: U
!!$       INTEGER                :: DepthMin,DepthMax
!!$       IF(ASSOCIATED(U%Left )) call depthTree(U%Left ,DepthMin,DepthMax)
!!$       IF(ASSOCIATED(U%Right)) call depthTree(U%Right,DepthMin,DepthMax)
!!$       IF(.NOT.ASSOCIATED(U%Left).AND..NOT.ASSOCIATED(U%Right)) THEN
!!$          DepthMin=MIN(DepthMin,U%Tier)
!!$          DepthMax=MAX(DepthMax,U%Tier)
!!$       ENDIF
!!$     END SUBROUTINE DepthTree
#endif
