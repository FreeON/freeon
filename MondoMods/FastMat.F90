MODULE FastMatrices
#ifdef PARALLEL_DEVELOPMENT
   USE DerivedTypes
   USE GlobalScalars   
   USE GlobalObjects
   USE MemMan
   USE Order
   USE MondoMPI
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
       TYPE(SRST), POINTER      :: Left 
       TYPE(SRST), POINTER      :: Right 
       TYPE(SRST), POINTER      :: Next
       REAL(DOUBLE), POINTER, DIMENSION(:,:)  :: MTrix  !-- Matrix block
    END TYPE SRST
    INTEGER :: SRSTCount
    TYPE(SRST),POINTER :: GlobalP
!======================================================================
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
    ENDIF

    CALL MPI_Gather(NAtms,1,MPI_INTEGER,CA%I(0),1,MPI_INTEGER,0,&
           MPI_COMM_WORLD,IErr)
    CALL MPI_Gather(B%NBlks,1,MPI_INTEGER,CB%I(0),1,MPI_INTEGER,0,&
           MPI_COMM_WORLD,IErr)
    CALL MPI_Gather(B%NNon0,1,MPI_INTEGER,CN%I(0),1,MPI_INTEGER,0,&
           MPI_COMM_WORLD,IErr)
    IF(MyID == 0) THEN
      GBNBlks = 0
      GBNNon0 = 0
      DO I = 0, NPrc-1
        GBNBlks = GBNBlks + CB%I(I)
        GBNNon0 = GBNNon0 + CN%I(I)
      ENDDO
      WRITE(*,*) 'GBNBlks = ', GBNBlks
      WRITE(*,*) 'GBNNon0 = ', GBNNon0
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
      CALL MPI_GatherV(B%MTrix%D,B%NNon0,MPI_DOUBLE_PRECISION,&
             C%MTrix%D,CN%I(0),DispN%I(0),MPI_DOUBLE_PRECISION,&
             0,MPI_COMM_WORLD,IErr)
    ELSE
      CALL MPI_GatherV(B%MTrix%D,B%NNon0,MPI_DOUBLE_PRECISION,&
             1.0D0,1,1,MPI_DOUBLE_PRECISION,&
             0,MPI_COMM_WORLD,IErr)
    ENDIF

    IF(MyID == 0) THEN
      CALL MPI_GatherV(B%ColPt%I(1),B%NBlks,MPI_INTEGER,&
             C%ColPt%I(1),CB%I(0),DispB%I(0),MPI_INTEGER,&
             0,MPI_COMM_WORLD,IErr)
    ELSE
      CALL MPI_GatherV(B%ColPt%I(1),B%NBlks,MPI_INTEGER,&
             1,1,1,MPI_INTEGER,&
             0,MPI_COMM_WORLD,IErr)
    ENDIF
 
    ! take differences
    LocalNBlks = B%NBlks
    DO I = 1, LocalNBlks-1
      B%BlkPt%I(I) = B%BlkPt%I(I+1)-B%BlkPt%I(I)
    ENDDO
    B%BlkPt%I(LocalNBlks) = B%NNon0-B%BlkPt%I(LocalNBlks)+1

    IF(MyID == 0) THEN
      CALL MPI_GatherV(B%BlkPt%I(1),LocalNBlks,MPI_INTEGER,&
             C%BlkPt%I(1),CB%I(0),DispB%I(0),MPI_INTEGER,&
             0,MPI_COMM_WORLD,IErr)
      PrevBlkSize = C%BlkPt%I(1)
      C%BlkPt%I(1) = 1
      DO I = 2, GBNBlks
        ! new pointer = previous cumulative pointer + prev size
        NewPt = C%BlkPt%I(I-1) + PrevBlkSize
        PrevBlkSize = C%BlkPt%I(I)
        C%BlkPt%I(I) = NewPt
      ENDDO
    ELSE
      CALL MPI_GatherV(B%BlkPt%I(1),B%NBlks,MPI_INTEGER,&
             1,1,1,MPI_INTEGER,&
             0,MPI_COMM_WORLD,IErr)
    ENDIF

    DO I = 1, NAtms
      B%RowPt%I(I) = B%RowPt%I(I+1)-B%RowPt%I(I)
    ENDDO
    IF(MyID == 0) THEN
      CALL MPI_GatherV(B%RowPt%I(1),NAtms,MPI_INTEGER,&
             C%RowPt%I(1),CA%I(0),DispA%I(0),MPI_INTEGER,&
             0,MPI_COMM_WORLD,IErr)
      PrevColSize = C%RowPt%I(1)
      C%RowPt%I(1) = 1
      DO I = 2, NAtoms+1
        NewPt = C%RowPt%I(I-1) + PrevColSize
        PrevColSize = C%RowPt%I(I)
        C%RowPt%I(I) = NewPt
      ENDDO
    ELSE
      CALL MPI_GatherV(B%RowPt%I(1),NAtms,MPI_INTEGER,&
             1,1,1,MPI_INTEGER,&
             0,MPI_COMM_WORLD,IErr)
    ENDIF
    
    IF(MyID == 0) THEN
      CALL Delete(CA)
      CALL Delete(CB)
      CALL Delete(CN)
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
      IF(MyID == 0) THEN
        !! WRITE(*,*) 'MyID = ',MyID, ' Testing row ', Row
      ENDIF
      IF(Row >= RowLimits(1) .AND. Row <= RowLimits(2)) THEN
        IF(MyID == 0) THEN
          !! WRITE(*,*) 'MyID = ',MyID, ' deleting row ', Row
        ENDIF

        NextR => R%Next
        !! CALL Print_SRST_1(R%RowRoot)
        CALL Delete_SRST_1(R%RowRoot)
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
      IF(MyID == 0) THEN
        WRITE(*,*) 'Col = ', A%L
      ENDIF
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
      IF(MyID == 0) THEN
        !! WRITE(*,*) 'Delete_SRST_1 , Col = ', A%L, ' Col 2 = ', A%R
      ENDIF
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
  SUBROUTINE PrintAllLinearRows(S)
    TYPE(FastMat),POINTER :: S,P
    TYPE(SRST),POINTER :: C
    INTEGER :: Row
    P => S%Next
    DO 
      IF(.NOT. ASSOCIATED(P)) EXIT
      !!  P is a valid row head
      Row = P%Row
      C => P%RowRoot
      DO 
        IF(.NOT. ASSOCIATED(C)) EXIT
        IF(C%L == C%R) THEN
          IF(MyID == 0) THEN
            !! WRITE(*,*) 'Row = ',Row, ' Col = ', C%L
          ENDIF
          IF(ASSOCIATED(C%Left) .OR. ASSOCIATED(C%Right)) THEN
            STOP 'ERR in PrintAllLinearRows.. Left or Right is not null!'
          ELSE
            !! WRITE(*,*) 'ASSERTION in PrintAllLinearRows is okay!'
          ENDIF
        ENDIF
        C =>  C%Next
      ENDDO
      P => P%Next
    ENDDO
  END SUBROUTINE PrintAllLinearRows

!======================================================================
  SUBROUTINE FlattenAllRows(S)
    TYPE(FastMat),POINTER :: S,P
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
  RECURSIVE SUBROUTINE PrintLeaf(A)
    TYPE(SRST),POINTER :: A
    IF(A%L == A%R) THEN
      !! do nothing!
    ELSE
      IF(ASSOCIATED(A%Left)) THEN
        CALL PrintLeaf(A%Left)
      ENDIF
      IF(ASSOCIATED(A%Right)) THEN
        CALL PrintLeaf(A%Right)
      ENDIF
    ENDIF
  END SUBROUTINE PrintLeaf
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

    StartTm = MPI_Wtime()
    IF(MyID == 0) THEN 
      WRITE(*,*) 'MyID=', MyID, ', FastMat_redistribute is entered...'
    ENDIF

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
    AllToAllBegTm = MPI_Wtime()
    
    CALL MPI_ALLTOALL( LocalDims%I(1,0),3,MPI_INTEGER, &
         RemoteDims%I(1,0),3,MPI_INTEGER,MONDO_COMM,IErr)
    AllToAllEndTm = MPI_Wtime()
    AllToAllTotTm = AllToAllEndTm - AllToAllBegTm
    CALL MPI_Allgather(AllToAllTotTm,1,MPI_DOUBLE_PRECISION,DimTmArr(0),1,&
      MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,IErr)
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
          CALL MPI_Send(SendBuffer(B),Num,MPI_DOUBLE_PRECISION,To,Tag,&
                 MONDO_COMM,IErr)
          CALL ErrChk(IErr,Sub)
        ENDIF
  
        IF(RecvFrQ%I(RecvFr) == 1) THEN
          From = RecvFr
          B = RcvBeg%I(From)
          E = RcvEnd%I(From)
          Num = E-B+1
          Tag = RecvFr
          CALL MPI_Recv(RecvBuffer(B),Num,MPI_DOUBLE_PRECISION,From,&
                 Tag,MONDO_COMM,Status,IErr)
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
          CALL MPI_Recv(RecvBuffer(B),Num,MPI_DOUBLE_PRECISION,From,&
                  Tag,MONDO_COMM,Status,IErr)
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
          CALL MPI_Send(SendBuffer(B),Num,MPI_DOUBLE_PRECISION,To,&
                 Tag,MONDO_COMM,IErr)
          CALL ErrChk(IErr,Sub)
        ENDIF
      ENDIF
        
    ENDDO

    CALL MPI_Gather(NumDblSent,1,MPI_INTEGER,IntArr(0),1,MPI_INTEGER,0,&
           MPI_COMM_WORLD,IErr)
    IF(MyID == 0) THEN
      DblSentMax = -1
      TotDblSent = 0
      DO I = 0, NPrc-1
        DblSentMax = Max(IntArr(I),DblSentMax)
        TotDblSent = TotDblSent + IntArr(I)
      ENDDO
      WRITE(*,*) 'Redistribute : DblSentMax = ',DblSentMax, ', TotDblSent = ',TotDblSent
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
    EndTm = MPI_Wtime()
    TotTm = EndTm - StartTm
    IF(MyID == ROOT) THEN
       WRITE(*,*) 'Total time to Redistribute_FastMat is ', TotTm
    ENDIF
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
        AbsRowPt(INAtms) = Row
        M = BSiz%I(Row)
        C => R%RowRoot
        DO
          IF(.NOT. ASSOCIATED(C)) EXIT
          Col = C%L
          IF(Col == C%R) THEN
            N = BSiz%I(Col)
            MN = M*N
            MTrix(MtxBegInd:MtxBegInd+MN-1)=PACK(C%MTrix,.TRUE.)
            BlkPt(Non0BlkP) = MtxBegInd
            ColPt(Non0BlkP) = Col
            Non0BlkP = Non0BlkP + 1
            MtxBegInd = MtxBegInd + MN 
          ENDIF
          RowPt(INAtms+1) = Non0BlkP 
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

!======================================================================
  SUBROUTINE Multiply_FASTMAT_SCALAR(A,Alpha,Perf_O)
    TYPE(FASTMAT),POINTER  :: A,R
    TYPE(SRST),POINTER :: C
    REAL(DOUBLE) :: Alpha,Op
    TYPE(TIME),OPTIONAL :: Perf_O
    INTEGER :: Row,Col,M,N
    Op=Zero
    IF(.NOT. ASSOCIATED(A)) THEN
      WRITE(*,*) 'ERR: A is null in Multiply_FASTMAT_SCALAR!'
      STOP
    ENDIF
    R => A%Next
    DO 
      IF(.NOT. ASSOCIATED(R)) EXIT
      Row = R%Row
      M = BSiz%I(Row)
      C => R%RowRoot
      DO 
        IF(.NOT. ASSOCIATED(C)) EXIT
        Col = C%L
        IF(Col == C%R) THEN
          N = BSiz%I(Col)
          Op = Op + M*N
          C%MTrix = C%MTrix*Alpha
        ENDIF
        C => C%Next 
      ENDDO
      R => R%Next
    ENDDO
    IF(PRESENT(Perf_O)) Perf_O%FLOP = Perf_O%FLOP+Op
  END SUBROUTINE Multiply_FASTMAT_SCALAR

!=================================================================
!   ADD A BLOCK TO THE FAST MATRIX DATA STRUCTURE
!=================================================================    
    SUBROUTINE AddFASTMATBlok(A,Row,Col,B)
      TYPE(FASTMAT),POINTER       :: A,C,D
      REAL(DOUBLE),DIMENSION(:,:) :: B
      TYPE(SRST), POINTER         :: P,Q
      INTEGER                     :: Row,Col,I,J,M,N
!-----------------------------------------------------------------
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
         M=BSiz%I(Row)
         N=BSiz%I(Col)
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

!======================================================================
!   ADD TWO DIFFERENT FAST MATRICES TOGETHER TO YEILD A THIRD:
!   NEED THIS ROUTINE TO MAINTAIN OO EQUIVALENCE WITH BCSR CALLS
!======================================================================
    SUBROUTINE Add_FASTMAT(A,B,C,Perf_O)
      TYPE(FASTMAT),POINTER  :: A,B,C
      TYPE(TIME),   OPTIONAL :: Perf_O         
      IF(.NOT.ASSOCIATED(A))  &
           CALL Halt(' A not associated in Add_FASTMAT ')
      IF(.NOT.ASSOCIATED(B))  &
           CALL Halt(' B not associated in Add_FASTMAT ')
      IF(.NOT.ASSOCIATED(C)) &
           CALL New_FASTMAT(C,0,(/0,0/))
      CALL AddIn_FASTMAT(A,C,Perf_O=Perf_O)
      CALL AddIn_FASTMAT(B,C,Perf_O=Perf_O)
    END SUBROUTINE Add_FASTMAT
!======================================================================
!   ADD TWO FAST MATRICES TOGETHER: A=Alpha*A+B
!======================================================================
    SUBROUTINE AddIn_FASTMAT(A,B,Alpha_O,Perf_O)
      TYPE(FASTMAT),POINTER  :: A,B
      TYPE(FASTMAT),POINTER  :: P,Q
      TYPE(SRST),   POINTER  :: U,V
      REAL(DOUBLE), OPTIONAL :: Alpha_O
      REAL(DOUBLE)           :: Alpha,Op
      TYPE(TIME),   OPTIONAL :: Perf_O         
      INTEGER                :: Row,Col,M,N,MN
!---------------------------------------------------------------------
      IF(PRESENT(Alpha_O))THEN
         Alpha=Alpha_O
      ELSE
         Alpha=One
      ENDIF 
      Op=Zero
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Skip pointers of A should be off since 
      ! its getting added to
      CALL SkipsOffQ(A%Alloc,'AddIn_Fastmat :: A ')
      ! B is walked upon, turn on skip pointers
      CALL SkipsOffQ(B%Alloc,'AddIn_FASTMAT : B')
      CALL SkipsOnFASTMAT(B)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      Q=>B%Next
      ! Go over rows of B
      DO WHILE(ASSOCIATED(Q))
         ! V is a Sparse Row ST belonging to B
         V=>Q%RowRoot
         Row=Q%Row
         M=BSiz%I(Row)
         ! Add/find P, a Sparse Row ST
         ! corresponding to Row of A
         P=>FindFastMatRow_1(B,Row)
         ! Initialize column nodes counter for Row of A
         SRSTCount=P%Nodes
         ! Go over columns of B
         DO
            ! Check for leaf nodes
            IF(V%L==V%R)THEN
               Col=V%L
               N=BSiz%I(Col)
               ! Find/add Col node in this Rows search tree
               U=>InsertSRSTNode(P%RowRoot,Col)
               IF(ASSOCIATED(U%MTrix))THEN
                  ! Resum the block
                  U%MTrix=Alpha*U%MTrix+V%MTrix
                  Op=Op+M*N
               ELSE
                  ! Initialize the block
                  ALLOCATE(U%MTrix(M,N),STAT=MemStatus)
                  CALL IncMem(MemStatus,0,N*M,'AddFASTMATBlok')
                  U%MTrix=V%MTrix
               ENDIF
            ENDIF
            IF(ASSOCIATED(V%Left))THEN
               V=>V%Left
            ELSEIF(ASSOCIATED(V%Right))THEN
               V=>V%Right
            ELSE
               EXIT
            ENDIF
         ENDDO
         ! Update counter
         P%Nodes=SRSTCount
         Q=>Q%Next
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      CALL SkipsOffFASTMAT(B)
      IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+Two*Op
    END SUBROUTINE AddIn_FASTMAT
!======================================================================
!  ADD A SCALAR TO A FAST MATRIX
!======================================================================
    SUBROUTINE Add_FASTMAT_SCLR(A,B,Perf_O)
      TYPE(FASTMAT),POINTER  :: A,P
      TYPE(SRST),   POINTER  :: U
      REAL(DOUBLE)           :: B,Op
      TYPE(TIME),   OPTIONAL :: Perf_O         
      INTEGER                :: M
!---------------------------------------------------------------------
      IF(.NOT.ASSOCIATED(A))  &
           CALL Halt(' A not associated in Add_FASTMAT_SCLR ')
      CALL SkipsOffQ(A%Alloc,'Add_Fastmat_SCLR :: A ')
      CALL SkipsOnFASTMAT(A)
      Op=Zero
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      P=>A%Next
      DO WHILE(ASSOCIATED(P))
         U=>P%RowRoot
         M=BSiz%I(P%Row)
         DO
            IF(U%L==U%R)THEN
               U%MTrix=U%MTrix+B
               Op=Op+M*BSiz%I(U%L)
            ENDIF
            IF(ASSOCIATED(U%Left))THEN
               U=>U%Left
            ELSEIF(ASSOCIATED(U%Right))THEN
               U=>U%Right
            ELSE
               EXIT
            ENDIF
         ENDDO
         P=>P%Next
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+Op      
   END SUBROUTINE Add_FASTMAT_SCLR
!======================================================================
!  MULTIPLY TWO FAST MATRICES TOGETHER: C=Alpha*A*B+Beta*C
!======================================================================
    SUBROUTINE Multiply_FASTMAT(A,B,C,Alpha_O,Beta_O,Perf_O)
      TYPE(FASTMAT),POINTER  :: A,B,C
      TYPE(FASTMAT),POINTER  :: P,Q,R
      TYPE(SRST),   POINTER  :: U,V,W
      REAL(DOUBLE), OPTIONAL :: Alpha_O,Beta_O
      REAL(DOUBLE)           :: Alpha,Beta,Op
      TYPE(TIME),   OPTIONAL :: Perf_O         
      INTEGER                :: Row,Col,K,MA,MB,NB,MN
!----------------------------------------------------------------------
      IF(PRESENT(Alpha_O))THEN
         Alpha=Alpha_O
      ELSE
         Alpha=One
      ENDIF
      IF(PRESENT(Beta_O))THEN
         Beta=Beta_O
      ELSE
         Beta=One
      ENDIF
      Op=Zero
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      IF(.NOT.ASSOCIATED(A))  &
         CALL Halt(' A not associated in Multiply_FASTMAT ')
      IF(.NOT.ASSOCIATED(B))  &
         CALL Halt(' B not associated in Multiply_FASTMAT ')
      IF(.NOT.ASSOCIATED(C)) &
         CALL New_FASTMAT(C,0,(/0,0/))
      ! Skip pointers of A and B are set on for column traversal
      CALL SkipsOffQ(A%Alloc,'Multiply_FASTMAT : A')
      CALL SkipsOffQ(B%Alloc,'Multiply_FASTMAT : B')
      CALL SkipsOnFASTMAT(A)
      CALL SkipsOnFASTMAT(B)
      ! Skip pointers of C should be off as its getting resummed
      CALL SkipsOffQ(C%Alloc,'Multiply_FASTMAT : C')
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      P=>A%Next
      ! Go over rows of A
      DO WHILE(ASSOCIATED(P))
         ! U is the Sparse Row ST containing A
         U=>P%RowRoot
         Row=P%Row
         MA=BSiz%I(Row)
         ! R is the Sparse Row ST containing 
         ! cooresponding coloumns of C
         R=>FindFastMatRow_1(C,Row)
         ! Initialize column nodes counter for Row of C
         SRSTCount=R%Nodes
         ! Go over columns of A
         DO
            ! Check for leaf nodes
            IF(U%L==U%R)THEN
               ! K is the kontraction index: cols of A, rows of B
               K=U%L
               ! Perform a hard find (no add) of row K in B, 
               ! proceed only if it exists
               Q=>FindFastMatRow_1(B,K,SoftFind_O=.TRUE.)         
               IF(ASSOCIATED(Q))THEN
                  MB=BSiz%I(K)
                  MN=MA*MB
                  ! V is the Sparse Row ST containing the columns of B
                  V=>Q%RowRoot
                  DO WHILE(ASSOCIATED(V))
                     ! Check for leaf node
                     IF(V%L==V%R)THEN
                        Col=V%L
                        NB=BSiz%I(Col)
                        ! Fast lookup/add of W=>C(A%Row,B%Col) 
                        W=>InsertSRSTNode(R%RowRoot,Col)
                        IF(ASSOCIATED(W%MTrix))THEN
                           ! Resum the block 
                           CALL DGEMM_NNC(MA,MB,NB,Alpha,Beta,U%MTrix(1,1),V%MTrix(1,1),W%MTrix(1,1))
                        ELSE
                           ! Initialize this block
                           ALLOCATE(W%MTrix(MA,NB),STAT=MemStatus)
                           CALL IncMem(MemStatus,0,MA*NB,'AddFASTMATBlok')
                           CALL DGEMM_NNC(MA,MB,NB,Alpha,Zero,U%MTrix(1,1),V%MTrix(1,1),W%MTrix(1,1))
                        ENDIF
                        Op=Op+DBLE(MN*NB) 
                     ENDIF
                     ! Tracing skip pointers over V, the coloumn index of B
                     IF(ASSOCIATED(V%Left))THEN
                        V=>V%Left
                     ELSEIF(ASSOCIATED(V%Right))THEN
                        V=>V%Right
                     ELSE
                        EXIT
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF
            ! Tracing skip pointers over U, the kontraction index K
            IF(ASSOCIATED(U%Left))THEN
               U=>U%Left
            ELSEIF(ASSOCIATED(U%Right))THEN
               U=>U%Right
            ELSE
               EXIT
            ENDIF
         ENDDO ! Cols of A
         ! Update C column nodes counter
         R%Nodes=SRSTCount
         ! Keep following As row links
         P=>P%Next
      ENDDO ! Rows of A
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      CALL SkipsOffFASTMAT(A)
      CALL SkipsOffFASTMAT(B)
      IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+Two*Op
    END SUBROUTINE Multiply_FASTMAT

!======================================================================
!  IN PLACE FILTERATION OF SMALL BLOCKS FROM A FAST MATRIX
!======================================================================
   SUBROUTINE FilterOut_FASTMAT(A,Tol_O,Perf_O)
      TYPE(FASTMAT),POINTER  :: A,P,Q
      REAL(DOUBLE),OPTIONAL  :: Tol_O         
      TYPE(TIME)  ,OPTIONAL  :: Perf_O         
      REAL(DOUBLE)           :: Tol,Op 
      TYPE(SRST),   POINTER  :: U,V
      INTEGER                :: Row,Col,M,N,MN
!----------------------------------------------------------------------
      IF(PRESENT(Tol_O))THEN
         Tol=Tol_O
      ELSE
         Tol=Thresholds%Trix
      ENDIF
      Op=Zero
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Cols of A are walked on as an ordinary binary tree, 
      ! skip pointers should be off
      CALL SkipsOffQ(A%Alloc,'FilterOut_FASTMAT : A')
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Q is the parent link of P, to allow excission of a row
      Q=>A
      P=>A%Next
      ! Go over rows of A
      DO WHILE(ASSOCIATED(P))
         ! Initialize column nodes counter for Row of A
         SRSTCount=P%Nodes
         ! Filter out columns for this row, incrementing FLOP count
         Op=Op+FilterOut_SRST(P%RowRoot,Tol)
         ! Reset SRST node counter 
         P%Nodes=SRSTCount
         ! Check to see if the SRST is now empty for this row
         IF(P%Nodes==0)THEN
            ! Skip over P in the LL 
            Q%Next=>P%Next
            ! Delete P
            CALL Delete_FASTMAT(P)
         ENDIF
         Q=>P
         P=>P%Next
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+Op
    END SUBROUTINE FilterOut_FASTMAT
!======================================================================
!
!======================================================================
    RECURSIVE FUNCTION FilterOut_SRST(A,Tol) RESULT(Op)
      TYPE(SRST),POINTER  :: A
      REAL(DOUBLE)        :: Tol,Op,FNorm
      INTEGER             :: MN
!---------------------------------------------------------------------
      IF(ASSOCIATED(A%Left))                 &
         Op=Op+FilterOut_SRST(A%Left,Tol)
      IF(ASSOCIATED(A%Right))                &      
         Op=Op+FilterOut_SRST(A%Right,Tol)
      ! Found a leaf node?
      IF(.NOT.ASSOCIATED(A%Left).AND.        &
         .NOT.ASSOCIATED(A%Right))THEN
         ! Check for leaf nodes with blocks
         IF(ASSOCIATED(A%MTrix))THEN
            ! Compute Frobinious norm
            FNorm=SQRT(DOT_PRODUCT(PACK(A%MTrix,.TRUE.),PACK(A%MTrix,.TRUE.)))
            MN=SIZE(A%MTrix,1)*SIZE(A%MTrix,2)
            Op=MN+6 ! Approximately 6 flops for the SQRT
            IF(FNorm<Tol)THEN
               ! Block is under threshold, delete the parent SRST node
               DEALLOCATE(A%MTrix,STAT=MemStatus)
               ! Deallocated M*N doubles
               CALL DecMem(MemStatus,0,MN)
               DEALLOCATE(A,STAT=MemStatus)
               CALL DecMem(MemStatus,6,0)
               ! Decrement SRST counter
               SRSTCount=SRSTCount-1
            ENDIF
         ELSE
            ! This must be a dangling node, and we are now backtracking, 
            ! deleting links leading to the matrix containg node.
            DEALLOCATE(A,STAT=MemStatus)
            CALL DecMem(MemStatus,6,0)
            SRSTCount=SRSTCount-1
         ENDIF
      ENDIF
    END FUNCTION FilterOut_SRST

!======================================================================
   SUBROUTINE Set_FASTMAT_EQ_BCSR(B,A)
      TYPE(FASTMAT),POINTER :: B,C
      TYPE(SRST),POINTER    :: S
      TYPE(BCSR)            :: A
      INTEGER               :: I,J,JP,P,N,M
!---------------------------------------------------------------------
      ! Check for prior allocation
      IF(ASSOCIATED(B%Next))THEN
         CALL Delete_FASTMAT(B)
         ! Begin with a new header Node     
         CALL New_FASTMAT(B,0,(/0,0/))
      ENDIF
      !
      DO I=1,A%NAtms
         M=BSiz%I(I)
         IF(A%RowPt%I(I+1)-A%RowPt%I(I)>1)THEN
            ! Set current row link 
            C=>FindFastMatRow_1(B,I)         
            DO JP=A%RowPt%I(I),A%RowPt%I(I+1)-1
               J=A%ColPt%I(JP)
               P=A%BlkPt%I(JP)
               N=BSiz%I(J)
               ! Add bloks to this sparse row search tree
               CALL AddFASTMATBlok(C,I,J,RESHAPE(A%MTrix%D(P:P+M*N-1),(/M,N/)))
!              CALL AddFASTMATBlok(C,I,J,VectToBlock(M,N,A%MTrix%D(P:P+M*N-1)))
            ENDDO
         ENDIF
      ENDDO
    END SUBROUTINE Set_FASTMAT_EQ_BCSR

!=================================================================
!
!=================================================================    
    SUBROUTINE Delete_FASTMAT(A,RowLimits_O)
      TYPE(FASTMAT),POINTER         :: A,B,C
      INTEGER,OPTIONAL,DIMENSION(2) :: RowLimits_O
      INTEGER,         DIMENSION(2) :: RowLimits
!-----------------------------------------------------------------

      STOP 'ERR: This is the old Delete_FASTMAT. Do not use this!'
    END SUBROUTINE Delete_FASTMAT

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
    SUBROUTINE New_FASTMAT(A,Row,Cols_O)
      TYPE(FASTMAT),POINTER         :: A
      INTEGER                       :: Row
      INTEGER,OPTIONAL,DIMENSION(2) :: Cols_O                       
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
        IF(A%Tier==0.AND.(Col<A%L.OR.Col>A%R)) &
           CALL Halt(' Logic error in InsertSRSTNode ')
!       Halt recursion if we've found a node with our column 
        IF(Col==A%L.AND.Col==A%R)THEN
           C=>A
           RETURN
        ENDIF
!       Halve interval
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
!   SET SKIP POINTERS FOR A FAST MATRIX
!======================================================================
    RECURSIVE SUBROUTINE SkipsOnFASTMAT(A)
      TYPE(FASTMAT),POINTER :: A,C
!----------------------------------------------------------------------
      ! Check to see if skip pointers are already on. This is not 
      ! nessesarily an error, as if already on its probably because
      ! the same pointer was passed twice to a procedure as in B=A*A.
      STOP 'Err: SkipsOnFASTMAT should not be called at all!'
    END SUBROUTINE SkipsOnFASTMAT
!======================================================================
!   SET SKIP POINTERS FOR THE SPARSE ROW SEARCH TREE (SRST), 
!   ALLOWING IN-PROCEDURE RECURSION OVER LEAF NODES
!======================================================================
  RECURSIVE SUBROUTINE SetSRSTSkipPtrs(A,B)
    TYPE(SRST),POINTER          :: A
    TYPE(SRST),POINTER,OPTIONAL :: B
    LOGICAL                     :: AssocB,AssocAL,AssocAR, &
                                   AssocBL,AssocBR
!----------------------------------------------------------------------
    STOP 'ERR: SetSRSTSkipPtrs should not be called!'
  END SUBROUTINE SetSRSTSkipPtrs
!======================================================================
!   QUERY TO SEE IF SKIP POINTERS ARE ON
!======================================================================
    SUBROUTINE SkipsOnQ(Key,Proc)
       INTEGER          :: Key
       CHARACTER(LEN=*) :: Proc
       IF(Key==ALLOCATED_SKIP_POINTERS_ON)RETURN
       IF(Key==ALLOCATED_SKIP_POINTERS_OFF)THEN
          CALL Halt(' Skip pointers incorrectly set in '//TRIM(Proc))
       ELSE
          CALL Halt(' Uninitialized allocation key in '//TRIM(Proc))
       ENDIF
     END SUBROUTINE SkipsOnQ
!======================================================================
!   UNSET SKIP POINTERS FOR A FAST MATRIX
!======================================================================
    RECURSIVE SUBROUTINE SkipsOffFASTMAT(A)
      TYPE(FASTMAT),POINTER :: A,C
!----------------------------------------------------------------------
      ! Check to see if skip pointers are already off. This is not 
      ! nessesarily an error, as if already on its probably because
      ! the same pointer was passed twice to a procedure as in B=A*A.

      STOP 'ERR: SkipsOffFASTMAT is no longer needed!'
    END SUBROUTINE SkipsOffFASTMAT
!======================================================================
!   REMOVE SKIP POINTERS, YEILDING A PLAIN 1-D BINARY SEARCH TREE
!======================================================================
    SUBROUTINE UnSetSRSTSkipPtrs(A)
      TYPE(SRST),POINTER :: A,B,C
        C=>A        
        DO ! Recur, always following the left link 
           IF(ASSOCIATED(C%Left))THEN
              C=>C%Left
           ELSEIF(ASSOCIATED(C%Right))THEN
              ! Check right node for skip pointer
              IF(C%Right%Tier<=C%Tier)THEN
                 B=>C%Right
                 ! Remove pointer
                 NULLIFY(C%Right)!=>NULL()
                 C=>B
              ELSE
                 C=>C%Right
              ENDIF
           ELSE
              ! All done
              EXIT
           ENDIF
        ENDDO
    END SUBROUTINE UnSetSRSTSkipPtrs
!======================================================================
!   QUERY TO SEE IF SKIP POINTERS ARE OFF
!======================================================================
    SUBROUTINE SkipsOffQ(Key,Proc)
       INTEGER          :: Key
       CHARACTER(LEN=*) :: Proc
       IF(Key==ALLOCATED_SKIP_POINTERS_OFF)RETURN
       IF(Key==ALLOCATED_SKIP_POINTERS_ON)THEN
          CALL Halt(' Skip pointers incorrectly set in '//TRIM(Proc))
       ELSE
          CALL Halt(' Uninitialized allocation key in '//TRIM(Proc))
       ENDIF
     END SUBROUTINE SkipsOffQ
#endif
  END MODULE FASTMATRICES
