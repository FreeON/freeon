MODULE FastMatrices
#ifdef PARALLEL_DEVELOPMENT
   USE DerivedTypes
   USE GlobalScalars   
   USE GlobalObjects
   USE MemMan
   USE Order
   USE MondoMPI
   USE Thresholding
!======================================================================
!   LINKED ROW LIST WITH COLOUMN INDEXED SEARCH TREE
!   A FAST O(N LG N) SPARSE MATRIX DATA STRUCTURE FOR ALL PROCCEDURES
!======================================================================
   TYPE FASTMAT
      INTEGER                  :: Alloc   !-- Allocation key
      INTEGER                  :: Nodes   !-- Number of nodes in this SRST
      INTEGER                  :: Row     !-- Row index of this link
      TYPE(SRST),    POINTER   :: RowRoot !-- Row link to a sparse row search tree  
      TYPE(FASTMAT), POINTER   :: Next    !-- Next row in linked list
   END TYPE FASTMAT
!======================================================================
!   SPARSE ROW SEARCH TREE: A FAST (LG N) SPARSE VECTOR DS
!======================================================================
    TYPE SRST
       INTEGER                                :: Alloc  !-- Allocation key
       INTEGER                                :: Row    !-- Row number
       INTEGER                                :: Tier   !-- Tree depth
       INTEGER                                :: Number !-- This nodes number in the tree
       INTEGER                                :: L,R    !-- Left and right coloumn interval bounds 
       TYPE(SRST), POINTER                    :: Left 
       TYPE(SRST), POINTER                    :: Right
       REAL(DOUBLE), POINTER, DIMENSION(:,:)  :: MTrix  !-- Matrix block
    END TYPE SRST

! the following line is from MemMan.F90
    INTEGER :: SRSTCount
  CONTAINS
!======================================================================
!
!======================================================================
  SUBROUTINE Redistribute_FastMat(A)
    TYPE(FastMat),    POINTER :: A
    TYPE(INT_RNK2)             :: LocalDims,RemoteDims
    TYPE(INT_VECT)             :: SndBeg,SndRow,SndCol,SndBlk,SndMtx,SndEnd, &
         RcvBeg,RcvRow,RcvCol,RcvBlk,RcvMtx,RcvEnd, &
         RecvReqst,SendReqst,ToDo,                  &
         SendOrder,RecvOrder,SendSched,RecvSched
    REAL(DOUBLE),ALLOCATABLE, &
         DIMENSION(:)  :: SendBuffer,RecvBuffer
    INTEGER,DIMENSION(2)       :: MyRow,AllCol
    INTEGER                    :: B,E,Num,Row,Col,Blk,Mtx,N,I,K,Tag,IErr,    &
         RE,To,From,NRecvs,NSends
    CHARACTER(LEN=20)          :: Sub='Redistribute_FastMat'
    REAL(DOUBLE)               :: StartTm,EndTm,TotTm
    StartTm = MPI_Wtime()
!----------------------------------------------------------------------
    CALL SkipsOffQ(A%Alloc,'Redistribute_FASTMAT : A')
    CALL SkipsOnFASTMAT(A)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    CALL New(LocalDims,(/3,NPrc-1/),(/1,0/))
    CALL New(RemoteDims,(/3,NPrc-1/),(/1,0/))
    AllCol=(/1,NAtoms/)
    DO N=0,NPrc-1 
       LocalDims%I(:,N)=MatDimensions(A,(/Beg%I(N),End%I(N)/),AllCol)
    ENDDO
    LocalDims%I(:,MyId)=(/0,0,0/)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    CALL MPI_ALLTOALL( LocalDims%I(1,0),3,MPI_INTEGER, &
         RemoteDims%I(1,0),3,MPI_INTEGER, &
         MONDO_COMM,IErr)
    CALL ErrChk(IErr,Sub)            
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Randomize communications
    CALL New(RecvOrder,NPrc)
    CALL New(SendOrder,NPrc)
    CALL New(RecvSched,NPrc)
    CALL New(SendSched,NPrc)
    NRecvs=0
    DO N=0,NPrc-1
       IF(N/=MyId.AND.RemoteDims%I(2,N)/=0)THEN
          NRecvs=NRecvs+1
          RecvSched%I(NRecvs)=N
          RecvOrder%I(NRecvs)=RANDOM((/1,100000/))
       ENDIF
    ENDDO
    NSends=0
    DO N=0,NPrc-1
       IF(N/=MyId.AND.RemoteDims%I(2,N)/=0)THEN
          NSends=NSends+1
          SendSched%I(NSends)=N
          SendOrder%I(NSends)=RANDOM((/1,100000/))
       ENDIF
    ENDDO
    ! Check for no work 
    IF(NRecvs==0.AND.NSends==0)THEN
       CALL Delete(LocalDims)
       CALL Delete(RemoteDims)
       CALL Delete(RecvOrder)
       CALL Delete(SendOrder)
       CALL Delete(RecvSched)
       CALL Delete(SendSched)
       CALL SkipsOffFASTMAT(A)
       CALL AlignNodes()
       RETURN
    ENDIF
!
    CALL Sort(RecvOrder,RecvSched,NRecvs)
    CALL Sort(SendOrder,SendSched,NSends)
    CALL Delete(RecvOrder)
    CALL Delete(SendOrder)
    CALL New(ToDo,NRecvs)
    CALL New(RecvReqst,NRecvs)
    CALL New(SendReqst,NSends)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    CALL New(SndBeg,NPrc-1,0)
    CALL New(SndRow,NPrc-1,0)
    CALL New(SndCol,NPrc-1,0)
    CALL New(SndBlk,NPrc-1,0)
    CALL New(SndMtx,NPrc-1,0)
    CALL New(SndEnd,NPrc-1,0)
    SndBeg%I(0)=1
    DO N=0,NPrc-1
       IF(N>0)  &
            SndBeg%I(N)=SndEnd%I(N-1)+1
       SndRow%I(N)=SndBeg%I(N)+3
       SndCol%I(N)=SndRow%I(N)+LocalDims%I(1,N)+1
       SndBlk%I(N)=SndCol%I(N)+LocalDims%I(2,N)+1
       SndMtx%I(N)=SndBlk%I(N)+LocalDims%I(2,N)+1
       SndEnd%I(N)=SndMtx%I(N)+LocalDims%I(3,N)+1
    ENDDO
    ! Allocate contigous memory for non-blocking sends
    ALLOCATE(SendBuffer(1:SndEnd%I(NPrc-1)),STAT=MemStatus)
    CALL IncMem(MemStatus,0,SndEnd%I(NPrc-1))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    CALL New(RcvBeg,NPrc-1,0)
    CALL New(RcvRow,NPrc-1,0)
    CALL New(RcvCol,NPrc-1,0)
    CALL New(RcvBlk,NPrc-1,0)
    CALL New(RcvMtx,NPrc-1,0)
    CALL New(RcvEnd,NPrc-1,0)
    RcvBeg%I(0)=1
    DO N=0,NPrc-1
       IF(N>0)  &
       RcvBeg%I(N)=RcvEnd%I(N-1)+1
       RcvRow%I(N)=RcvBeg%I(N)+3
       RcvCol%I(N)=RcvRow%I(N)+RemoteDims%I(1,N)+1
       RcvBlk%I(N)=RcvCol%I(N)+RemoteDims%I(2,N)+1
       RcvMtx%I(N)=RcvBlk%I(N)+RemoteDims%I(2,N)+1
       RcvEnd%I(N)=RcvMtx%I(N)+RemoteDims%I(3,N)+1
    ENDDO
    ! Allocate contigous memory for non-blocking recieves
    RE=RcvEnd%I(NPrc-1)
    ALLOCATE(RecvBuffer(1:RE),STAT=MemStatus)
    CALL IncMem(MemStatus,0,RcvEnd%I(NPrc-1))
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Post recieves 
    DO N=1,NRecvs
       From=RecvSched%I(N)
       B=RcvBeg%I(From)
       E=RcvEnd%I(From)
       Num=E-B+1
       Tag=From*MaxProc+MyId
       CALL MPI_IRECV(RecvBuffer(B),Num,MPI_DOUBLE_PRECISION,  &
            From,Tag,MONDO_COMM,RecvReqst%I(N),IErr)
       CALL ErrChk(IErr,Sub)            
    ENDDO
    ! Wait for all recieves to post before executing sends.
    ! This way blocking sends should be fast
    CALL AlignNodes()
    ! Post sends
    DO N=1,NSends
       To=SendSched%I(N)
       B=SndBeg%I(To)
       E=SndEnd%I(To)
       Row=SndRow%I(To)
       Col=SndCol%I(To)
       Blk=SndBlk%I(To)
       Mtx=SndMtx%I(To)
       Num=E-B+1
       Tag=MyId*MaxProc+To
       CALL PackFastMat(A,(/Beg%I(To),End%I(To)/),          &
            SendBuffer(B),SendBuffer(B+1),SendBuffer(B+2),  &
            SendBuffer(Row:Col-1),SendBuffer(Col:Blk-1),    &
            SendBuffer(Blk:MTx-1),SendBuffer(Mtx:E))
       CALL MPI_SEND(SendBuffer(B),Num,MPI_DOUBLE_PRECISION,To,Tag,MONDO_COMM,IErr)
!       CALL MPI_ISEND(SendBuffer(B),Num,MPI_DOUBLE_PRECISION,To, &
!            Tag,MONDO_COMM,SendReqst%I(N),IErr)        
       CALL ErrChk(IErr,Sub)            
    ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Delete rows from A that have been sent to another processor
    CALL SkipsOffFASTMAT(A)
    CALL Delete_FASTMAT(A,(/Beg%I(MyId),End%I(MyId)/))
    CALL SkipsOnFASTMAT(A)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    DO 
       ! Wait for some buffers to fill
       CALL WaitSome(RecvReqst,ToDo)
       ! If done, exit
       IF(ToDo%I(1)==FAIL)EXIT
       ! Go over filled buffers
       DO N=1,SIZE(ToDo%I)
          From=RecvSched%I(ToDo%I(N))
          B=RcvBeg%I(From)
          Row=RcvRow%I(From)
          Col=RcvCol%I(From)
          Blk=RcvBlk%I(From)
          Mtx=RcvMtx%I(From)
          E=RcvEnd%I(From)
          CALL UnPackNSumFastMat(A,OffSt%I(MyId),RecvBuffer(B),                &
               RecvBuffer(B+1),RecvBuffer(B+2),              &
               RecvBuffer(Row:Col-1),RecvBuffer(Col:Blk-1),  &
               RecvBuffer(Blk:Mtx-1),RecvBuffer(Mtx:E))
       ENDDO
    ENDDO
    CALL AlignNodes()      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Delete contigous memory
    DEALLOCATE(SendBuffer,STAT=MemStatus)
    CALL DecMem(MemStatus,SndEnd%I(NPrc-1),0)
    DEALLOCATE(RecvBuffer,STAT=MemStatus) 
    ! Delete the rest
    CALL DecMem(MemStatus,RcvEnd%I(NPrc-1),0)
    CALL Delete(LocalDims)
    CALL Delete(RemoteDims)
    CALL Delete(RecvReqst)
    CALL Delete(SendReqst)  
    CALL Delete(ToDo)
    CALL Delete(SndBeg)
    CALL Delete(SndRow)
    CALL Delete(SndCol)
    CALL Delete(SndBlk)
    CALL Delete(SndMtx)
    CALL Delete(SndEnd)
    CALL Delete(RcvBeg)
    CALL Delete(RcvRow)
    CALL Delete(RcvCol)
    CALL Delete(RcvBlk)
    CALL Delete(RcvMtx)
    CALL Delete(RcvEnd)
    CALL Delete(RecvSched)
    CALL Delete(SendSched)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    CALL SkipsOffFASTMAT(A)
    CALL AlignNodes()
    EndTm = MPI_Wtime()
    TotTm = EndTm - StartTm
    IF(MyID == ROOT) THEN
       WRITE(*,*) 'Total time to Redistribute_FastMat is ', TotTm
    ENDIF
  END SUBROUTINE Redistribute_FastMat
!======================================================================
!
!======================================================================
    SUBROUTINE PackFastMat(A,RowLimits,NAtms,NBlks,Non0s,  &
                           RowPt,ColPt,BlkPt,MTrix)
      TYPE(FastMat),POINTER     :: A,C
      TYPE(SRST),   POINTER     :: S
      INTEGER,DIMENSION(2)      :: RowLimits
      REAL(DOUBLE)              :: NAtms,NBlks,Non0s
      REAL(DOUBLE),DIMENSION(:) :: RowPt,ColPt,BlkPt,MTrix
      INTEGER                   :: OffSt,I,Row,Col,Blk,M,N,MN,P
!----------------------------------------------------------------------
      CALL SkipsOnQ(A%Alloc,'PackFastmat')
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      C=>A%Next
      P=1
      Blk=1
      NAtms=0
      RowPt(1)=1 
      ! Go over Rows of the LL
      DO WHILE(ASSOCIATED(C))
         Row=C%Row
         IF(Row>=RowLimits(1).AND.Row<=RowLimits(2))THEN
            ! Offset by the first row to compact the BCSR
            I=Row-RowLimits(1)+1
            NAtms=NAtms+1
            S=>C%RowRoot
            M=BSiz%I(Row) 
            DO
               IF(S%L==S%R)THEN
                  Col=S%L
                  N=BSiz%I(Col)
                  MN=M*N
                  MTrix(Blk:Blk+MN-1)=PACK(S%MTrix,.TRUE.)
                  BlkPt(P)=Blk
                  ColPt(P)=Col             
                  P=P+1
                  Blk=Blk+MN
                  RowPt(I+1)=P
               ENDIF
               IF(ASSOCIATED(S%Left))THEN
                  S=>S%Left
               ELSEIF(ASSOCIATED(S%Right))THEN
                  S=>S%Right
               ELSE
                  EXIT
               ENDIF
            ENDDO
         ENDIF
         C=>C%Next
      ENDDO
      NBlks=P-1
      Non0s=Blk-1
!!$      WRITE(*,*)' NAts = ',NAtms,' NBlks = ',NBlks,' NNon0s = ',Non0s
    END SUBROUTINE PackFastMat
!======================================================================
!
!======================================================================
    SUBROUTINE UnPackNSumFastMat(A,OffSt,NAtms,NBlks,Non0s,  &
                                 RowPt,ColPt,BlkPt,MTrix,Perf_O)
      TYPE(FastMat),POINTER     :: A,P
      TYPE(SRST)   ,POINTER     :: U
      REAL(DOUBLE)              :: NAtms,NBlks,Non0s,Op
      REAL(DOUBLE),DIMENSION(:) :: RowPt,ColPt,BlkPt,MTrix
      INTEGER                   :: I,J,M,N,MN,OffSt,Blk,Row,Col, &
                                   IAtms,BegRow,EndRow
      TYPE(TIME),   OPTIONAL :: Perf_O      
!----------------------------------------------------------------------
      Op=Zero
      IAtms=NAtms
      ! Go over the number of rows (I) and the absolute row position (Row)
      DO I=1,IAtms; 
         Row=I+OffSt
         M=BSiz%I(Row) 
         ! Add/find P, a Sparse Row ST
         ! corresponding to Row of A
         P=>FindFASTMATRow(A,Row,HardFind_O=.TRUE.)
         ! Go over columns
         BegRow=RowPt(I)
         EndRow=RowPt(I+1)-1
         DO J=BegRow,EndRow
            ! Initialize column nodes counter for Row
            SRSTCount=P%Nodes
            Col=ColPt(J)
            Blk=BlkPt(J)
            N=BSiz%I(Col)
            ! Find/add Col node in this Rows search tree
            U=>InsertSRSTNode(P%RowRoot,Col)
            IF(ASSOCIATED(U%MTrix))THEN
               ! Resum the block
               U%MTrix=U%MTrix+RESHAPE(MTrix(Blk:Blk+M*N-1),(/M,N/))
               Op=Op+M*N
            ELSE
               ! Initialize the block
               ALLOCATE(U%MTrix(M,N),STAT=MemStatus)
               CALL IncMem(MemStatus,0,N*M,'AddFASTMATBlok')
               U%MTrix=RESHAPE(MTrix(Blk:Blk+M*N-1),(/M,N/))
            ENDIF
         ENDDO
      ENDDO
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+Op  
    END SUBROUTINE UnPackNSumFastMat
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
         P=>FindFASTMATRow(B,Row)
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
         R=>FindFASTMATRow(C,Row)
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
               Q=>FindFASTMATRow(B,K,HardFind_O=.TRUE.)         
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
!  MULTIPLY A FAST MATRIX BY A SCALAR: A=Alpha*A
!======================================================================
    SUBROUTINE Multiply_FASTMAT_SCALAR(A,Alpha,Perf_O)
      TYPE(FASTMAT),POINTER  :: A
      TYPE(FASTMAT),POINTER  :: P
      TYPE(SRST),   POINTER  :: U
      REAL(DOUBLE)           :: Alpha,Op
      TYPE(TIME),   OPTIONAL :: Perf_O         
      INTEGER                :: Row,MA
!----------------------------------------------------------------------
      IF(.NOT.ASSOCIATED(A))  &
         CALL Halt(' A not associated in Multiply_FASTMAT ')
      ! Skip pointers of A are set on for column traversal
      CALL SkipsOffQ(A%Alloc,'Multiply_FASTMAT_SCALAR : A')
      CALL SkipsOnFASTMAT(A)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      Op=Zero
      P=>A%Next
      ! Go over rows of A
      DO WHILE(ASSOCIATED(P))
         U=>P%RowRoot  ! U is the Sparse Row ST containing A
         Row=P%Row
         MA=BSiz%I(Row)
         ! Go over columns of A
         DO
            ! Check for leaf nodes
            IF(U%L==U%R)THEN
               Op=Op+MA*BSiz%I(U%L)
               U%MTrix=U%MTrix*Alpha
            ENDIF
            ! Tracing skip pointers over U, the Col index
            IF(ASSOCIATED(U%Left))THEN
               U=>U%Left
            ELSEIF(ASSOCIATED(U%Right))THEN
               U=>U%Right
            ELSE
               EXIT
            ENDIF
         ENDDO ! Cols of A
         ! Keep following row links in A
         P=>P%Next
      ENDDO ! Rows of A
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      CALL SkipsOffFASTMAT(A)
      IF(PRESENT(Perf_O))Perf_O%FLOP=Perf_O%FLOP+Op
    END SUBROUTINE Multiply_FASTMAT_SCALAR
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
!   CONVERT A DISTRIBUTED FAST MATRIX INTO A ROOTED BCSR MATRIX
!======================================================================
    SUBROUTINE Set_BCSR_EQ_DFASTMAT(C,A)
      TYPE(FASTMAT),POINTER :: A
      TYPE(BCSR)            :: B,C
      TYPE(INT_VECT)        :: MA,MB,MN,NA,NB,NN 
      INTEGER               :: I,K,NAtms
!-----------------------------------------------------------------------
      ! Obtain local, intermediate BCSR matrices from the kosher rows of A
      CALL Set_BCSR_EQ_FASTMAT(B,A,OffSt%I(MyId),(/Beg%I(MyId),End%I(MyId)/))
      ! Global limits for rooted BCSR matrix
      C%NAtms=NAtoms
      C%NBlks=Reduce(B%NBlks)
      C%NNon0=Reduce(B%NNon0)
      ! Allocate rooted BCSR matrix
      IF(MyId==ROOT) &
           CALL New(C,(/C%NAtms,C%NBlks,C%NNon0/))
      ! Allocate displacement indecies for gather
      IF(MyId==ROOT)THEN
         CALL New(MA,NPrc,M_O=0)
         CALL New(MB,NPrc,M_O=0)
         CALL New(MN,NPrc,M_O=0)
         CALL New(NA,NPrc,M_O=0)
         CALL New(NB,NPrc,M_O=0)
         CALL New(NN,NPrc,M_O=0)
      ENDIF
      ! Number of atoms for this node
      NAtms=End%I(MyId)-Beg%I(MyId)+1
      IF(MyID==NPrc-1)NAtms=NAtms+1
      ! Gather the indecies to root
      CALL Gather(NAtms  ,MA)
      CALL Gather(B%NBlks,MB)
      CALL Gather(B%NNon0,MN)
      ! Calculate displacement indeces (offsets)
      IF(MyID==ROOT)THEN
         NA%I(0)=  NAtms
         NB%I(0)=B%NBlks
         NN%I(0)=B%NNon0
         DO I=1,NPrc-1
            NA%I(I)=MA%I(I)+NA%I(I-1)
            NB%I(I)=MB%I(I)+NB%I(I-1)
            NN%I(I)=MN%I(I)+NN%I(I-1)
         ENDDO
         DO I=NPrc,1,-1
            NA%I(I)=NA%I(I-1)
            NB%I(I)=NB%I(I-1)
            NN%I(I)=NN%I(I-1)
         ENDDO
         NA%I(0)=0
         NB%I(0)=0
         NN%I(0)=0
         ! Copy local portions of the matirx
         CALL SetEq(C%RowPt,B%RowPt,  NAtms)
         CALL SetEq(C%ColPt,B%ColPt,B%NBlks)
         CALL SetEq(C%BlkPt,B%BlkPt,B%NBlks)
         CALL SetEq(C%MTrix,B%MTrix,B%NNon0)           
      ENDIF
      ! Gather the rest of the matrix to ROOT
      CALL Gather(B%RowPt,C%RowPt,  NAtms,MA,NA)
      CALL Gather(B%ColPt,C%ColPt,B%NBlks,MB,NB)
      CALL Gather(B%BlkPt,C%BlkPt,B%NBlks,MB,NB)
      CALL Gather(B%MTrix,C%MTrix,B%NNon0,MN,NN)
      ! Add in offsets to achieve correct indexing
      IF(MyId==ROOT)THEN
         DO I=1,NPrc-1
            DO K=NA%I(I)+1,NA%I(I+1)
               C%RowPt%I(K)=C%RowPt%I(K)+NB%I(I)
            ENDDO
            DO K=NB%I(I)+1,NB%I(I+1)
               C%BlkPt%I(K)=C%BlkPt%I(K)+NN%I(I)
            ENDDO
         ENDDO
      ENDIF
      ! Clean up 
      IF(MyId==ROOT)THEN
         CALL Delete(MA)
         CALL Delete(MB)
         CALL Delete(MN)
         CALL Delete(NA)
         CALL Delete(NB)
         CALL Delete(NN)
         CALL Delete(B)
      ENDIF
!!$IF(MyId==0)THEN
!!$WRITE(*,*)' C%RowPt = ',C%RowPt%I
!!$WRITE(*,*)' C%ColPt = ',C%ColPt%I
!!$WRITE(*,*)' C%BlkPt = ',C%BlkPt%I
!!$WRITE(*,*)' C%MTrix = ',C%MTrix%D
!!$ENDIF
    END SUBROUTINE Set_BCSR_EQ_DFASTMAT
!======================================================================
!   CONVERT A FAST MATRIX INTO A BCSR MATRIX
!======================================================================
    SUBROUTINE Set_BCSR_EQ_FASTMAT(B,A,OffSt,RowLimits_O)
      TYPE(FASTMAT),POINTER         :: A,C
      TYPE(SRST),POINTER            :: S
      TYPE(BCSR)                    :: B
      INTEGER,OPTIONAL,DIMENSION(2) :: RowLimits_O
      INTEGER,DIMENSION(2)          :: RowLimits
      INTEGER,DIMENSION(3)          :: MatDims
      INTEGER                       :: OffSt,RowOff
      INTEGER                       :: P,Q,I,J,M,N,MN
!---------------------------------------------------------------------
      IF(PRESENT(RowLimits_O))THEN
         RowLimits=RowLimits_O
      ELSE
         RowLimits=(/1,NAtoms/)
      ENDIF
      CALL SkipsOffQ(A%Alloc,'Set_BCSR_EQ_FASTMAT')
      CALL SkipsOnFASTMAT(A)
      MatDims=MatDimensions(A,RowLimits,(/1,NAtoms/))
      IF(MatDims(1)==0)THEN
         B%NAtms=Zero
         RETURN
      END IF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(AllocQ(B%Alloc))THEN 
         CALL Delete(B)
      ENDIF
      CALL New(B,N_O=MatDims,OnAll_O=.TRUE.)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      C=>A%Next
      P=1
      Q=1
      B%RowPt%I(1)=1 
      DO WHILE(ASSOCIATED(C))
         RowOff=C%Row-OffSt
         IF(C%Row>=RowLimits(1).AND.C%Row<=RowLimits(2))THEN
            S=>C%RowRoot
            M=BSiz%I(C%Row) 
            DO
               IF(S%L==S%R)THEN
                  J=S%L
                  N=BSiz%I(J)
                  MN=M*N
                  B%MTrix%D(Q:Q+MN-1)=PACK(S%MTrix,.TRUE.)
!                  CALL BlockToBlock(M,N,0,0,S%MTrix,B%MTrix%D(Q:Q+MN-1))
                  B%BlkPt%I(P)=Q
                  B%ColPt%I(P)=J               
                  P=P+1
                  Q=Q+MN
                  B%RowPt%I(RowOff+1)=P
               ENDIF
               IF(ASSOCIATED(S%Left))THEN
                  S=>S%Left
               ELSEIF(ASSOCIATED(S%Right))THEN
                  S=>S%Right
               ELSE
                  EXIT
               ENDIF
            ENDDO
         ENDIF
         C=>C%Next
      ENDDO
      B%NAtms=NAtoms
      B%NBlks=P-1
      B%NNon0=Q-1
!      WRITE(*,*)MyId,' B%NBlks = ',B%NBlks
!      WRITE(*,*)' B%RowPt = ',B%RowPt%I
!      WRITE(*,*)' B%ColPt = ',B%ColPt%I
!      WRITE(*,*)' B%BlkPt = ',B%BlkPt%I
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CALL SkipsOffFASTMAT(A)
    END SUBROUTINE Set_BCSR_EQ_FASTMAT
!======================================================================
!   COMPUTE BCSR MATRIX DIMENSIONS CORESPONDING TO A FAST MATRIX
!======================================================================
    FUNCTION MatDimensions(A,RowLimits,ColLimits) RESULT(MatDim)
      TYPE(FASTMAT),POINTER :: A,C
      INTEGER               :: Atms,Blks,Non0s             
      INTEGER,DIMENSION(2)  :: RowDim,RowLimits,ColLimits
      INTEGER,DIMENSION(3)  :: MatDim
!-----------------------------------------------------------------------
      CALL SkipsOnQ(A%Alloc,'MatDimensions')
      Atms=0
      RowDim=0
      C=>A%Next
      DO WHILE(ASSOCIATED(C))
         IF(C%Row>=RowLimits(1).AND.C%Row<=RowLimits(2))THEN
            Atms=Atms+1
            RowDim=RowDim+RowDimensions(C%RowRoot,C%Row,ColLimits)
         ENDIF
         C=>C%Next
      ENDDO
      ! Atms+1 since CSR row pts addressing goes like RowPt(I+1)-1
      MatDim=(/Atms+1,RowDim(1),RowDim(2)/)
    END FUNCTION MatDimensions
!======================================================================
!   COMPUTE SPARSE ROW DIMENSIONS FROM A SPARSE ROW SEARCH TREE (SRST)
!======================================================================
    FUNCTION RowDimensions(A,Row,ColLimits) RESULT(RowDim)
      TYPE(SRST),POINTER   :: A,C
      INTEGER              :: Row,M,N
      INTEGER,DIMENSION(2) :: RowDim,ColLimits
!---------------------------------------------------------------------
      M=BSiz%I(Row)
      C=>A
      RowDim=0
      DO 
         IF(ASSOCIATED(C%Left))THEN
            C=>C%Left
         ELSEIF(ASSOCIATED(C%Right))THEN
            IF(C%L==C%R.AND.C%L>=ColLimits(1).AND.C%L<=ColLimits(2))THEN
               RowDim(1)=RowDim(1)+1
               RowDim(2)=RowDim(2)+M*BSiz%I(C%L)
            ENDIF
            C=>C%Right
         ELSE
            IF(C%L==C%R.AND.C%L>=ColLimits(1).AND.C%L<=ColLimits(2))THEN
               RowDim(1)=RowDim(1)+1
               RowDim(2)=RowDim(2)+M*BSiz%I(C%L)
            ENDIF
            EXIT
         ENDIF
      ENDDO
    END FUNCTION RowDimensions
!======================================================================
!
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
            C=>FindFASTMATRow(B,I)         
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
      CALL SkipsOffQ(A%Alloc,'Delete_FASTMAT')
      IF(PRESENT(RowLimits_O))THEN
         RowLimits=RowLimits_O
      ELSE
         RowLimits=(/1000000000,-1/)
      ENDIF
      B=>A
      DO 
         IF(ASSOCIATED(A%Next%Next))THEN
            ! We are somewhere in the midle of the list
            IF(A%Next%Row<RowLimits(1).OR.  &
               A%Next%Row>RowLimits(2))THEN
               C=>A%Next%Next
               CALL Delete_SRST(A%Next%RowRoot)
               DEALLOCATE(A%Next)
               A%Next=>C
            ELSE
               A=>A%Next
            ENDIF
         ELSEIF(ASSOCIATED(A%Next))THEN
            ! We are one before the end of the list 
            IF(A%Next%Row<RowLimits(1).OR.  &
               A%Next%Row>RowLimits(2))THEN
               CALL Delete_SRST(A%Next%RowRoot)
               DEALLOCATE(A%Next)
               NULLIFY(A%Next)
            ENDIF
            EXIT
         ENDIF
      ENDDO
!!$      A=>B
!!$      DO WHILE(ASSOCIATED(A))
!!$         WRITE(*,*)MyId,'ASSOCIATEDA = ',ASSOCIATED(A),' A%Row = ',A%Row
!!$         IF(ASSOCIATED(A%Next))THEN
!!$            A=>A%Next
!!$         ELSE
!!$            EXIT
!!$         ENDIF
!!$      ENDDO
      A=>B
    END SUBROUTINE Delete_FASTMAT
!=================================================================
!
!=================================================================    
    SUBROUTINE Delete_FASTMAT_LINK(A)
      TYPE(FASTMAT),POINTER :: A,B
      B=>A%Next
      CALL Delete_SRST(A%RowRoot)
      DEALLOCATE(A)
      A=>B
    END SUBROUTINE DELETE_FASTMAT_LINK
!=================================================================
!
!=================================================================    
    RECURSIVE SUBROUTINE Delete_FASTMAT2(A)
      TYPE(FASTMAT),POINTER :: A
!-----------------------------------------------------------------
      CALL SkipsOffQ(A%Alloc,'Delete_FASTMAT')
      IF(ASSOCIATED(A%Next))THEN
         CALL Delete_FASTMAT2(A)
      ELSE
         CALL Delete_SRST(A%RowRoot)
         DEALLOCATE(A)
      ENDIF
    END SUBROUTINE Delete_FASTMAT2
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
!   ADD A BLOCK TO THE FAST MATRIX DATA STRUCTURE
!=================================================================    
    SUBROUTINE AddFASTMATBlok(A,Row,Col,B)
      TYPE(FASTMAT),POINTER       :: A,C,D
      REAL(DOUBLE),DIMENSION(:,:) :: B
      TYPE(SRST), POINTER         :: P,Q
      INTEGER                     :: Row,Col,I,J,M,N
!-----------------------------------------------------------------
      ! Find the current row in the fast matrix
      C=>FindFASTMATRow(A,Row) 
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
!   FIND OR CREATE A ROW IN A FAST MATRIX
!=================================================================    
    FUNCTION FindFASTMATRow(A,Row,HardFind_O) RESULT(C)
      TYPE(FASTMAT),POINTER       :: A,C,D
      LOGICAL, OPTIONAL           :: HardFind_O
      INTEGER                     :: Row
      LOGICAL                     :: HardFind
!-----------------------------------------------------------------
      IF(PRESENT(HardFind_O))THEN
         HardFind=HardFind_O
      ELSE
         HardFind=.FALSE.
      ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      C=>A
      DO ! Search through linked list for row
         IF(C%Row==Row)THEN
            EXIT
         ELSEIF(ASSOCIATED(C%Next))THEN
            IF(Row>C%Row.AND.Row<C%Next%Row)THEN
               IF(HardFind)THEN
                  ! Link does not exist, return null 
                  ! under HardFind directive
                  NULLIFY(C)!=>NULL()
                  RETURN
               ENDIF
               ! Insert link to preserve assention
               D=>C%Next
               NULLIFY(C%Next)
               CALL New_FASTMAT(C%Next,Row)
               C%Next=>D              
            ELSE
               C=>C%Next
            ENDIF
         ELSE
            IF(HardFind)THEN
               ! Link does not exist, return null 
               ! under HardFind directive
               NULLIFY(C)!=>NULL()
               RETURN
            ELSE
               ! If row doesnt show, add link to end of list
               CALL New_FASTMAT(C%Next,Row)
               C=>C%Next
            ENDIF
            EXIT
         ENDIF
      ENDDO
    END FUNCTION FindFASTMATRow
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
      NULLIFY(A%Next)!=>NULL()
      CALL New_SRST(A%RowRoot,Cols(1),Cols(2),0)
    END SUBROUTINE New_FASTMAT
!======================================================================
!   SET SKIP POINTERS FOR A FAST MATRIX
!======================================================================
    RECURSIVE SUBROUTINE SkipsOnFASTMAT(A)
      TYPE(FASTMAT),POINTER :: A,C
!----------------------------------------------------------------------
      ! Check to see if skip pointers are already on. This is not 
      ! nessesarily an error, as if already on its probably because
      ! the same pointer was passed twice to a procedure as in B=A*A.
      IF(A%Alloc==ALLOCATED_SKIP_POINTERS_ON)RETURN
      ! Start past head link
      C=>A%Next
      ! Go over each link 
      DO WHILE(ASSOCIATED(C))
         IF(ASSOCIATED(C%RowRoot))THEN 
            ! Set skip pointers for this sparse row search tree
            IF(ASSOCIATED(C%RowRoot%Left).AND.ASSOCIATED(C%RowRoot%Right))THEN
               CALL SetSRSTSkipPtrs(C%RowRoot%Left,C%RowRoot%Right)
            ELSEIF(ASSOCIATED(C%RowRoot%Left))THEN
               CALL SetSRSTSkipPtrs(C%RowRoot%Left)
            ELSEIF(ASSOCIATED(C%RowRoot%Right))THEN
               CALL SetSRSTSkipPtrs(C%RowRoot%Right)
            ENDIF
         ENDIF
         C=>C%Next
      ENDDO
      A%Alloc=ALLOCATED_SKIP_POINTERS_ON
    END SUBROUTINE SkipsOnFASTMAT
!======================================================================
!   SET SKIP POINTERS FOR THE SPARSE ROW SEARCH TREE (SRST), 
!   ALLOWING IN-PROCEDURE RECURSION OVER LEAF NODES
!======================================================================
  RECURSIVE SUBROUTINE SetSRSTSkipPtrs(A,B)
    TYPE(SRST),POINTER          :: A
    TYPE(SRST),POINTER,OPTIONAL :: B
    LOGICAL                     :: AssocAL,AssocAR, &
                                   AssocBL,AssocBR
!----------------------------------------------------------------------
    IF(PRESENT(B))THEN
       AssocB=ASSOCIATED(B)
       IF(A%L==A%R)THEN ! Done!
          A%Right=>B    ! Set the skip pointer
       ELSE
          AssocAL=ASSOCIATED(A%Left)
          AssocAR=ASSOCIATED(A%Right)
          ! Follow the left link (A) 
          IF(AssocAL.AND.AssocAR)THEN
             CALL SetSRSTSkipPtrs(A%Left,A%Right)
             CALL SetSRSTSkipPtrs(A%Right,B)
          ELSEIF(AssocAL)THEN
             CALL SetSRSTSkipPtrs(A%Left,B)
          ELSEIF(AssocAR)THEN
             CALL SetSRSTSkipPtrs(A%Right,B)
          ENDIF
       ENDIF
       ! Check for dead end
       IF(B%L==B%R)RETURN 
       AssocBL=ASSOCIATED(B%Left)
       AssocBR=ASSOCIATED(B%Right)
       ! Follow the right link (B)
       IF(AssocBL.AND.AssocBR)THEN
          CALL SetSRSTSkipPtrs(B%Left,B%Right)
       ELSEIF(AssocBL)THEN
          CALL SetSRSTSkipPtrs(B%Left)
       ELSEIF(AssocBR)THEN
          CALL SetSRSTSkipPtrs(B%Right)
       ENDIF
    ELSE ! There was no right link passed in ...
       ! Check for dead end
       IF(A%L==A%R)RETURN 
       AssocAL=ASSOCIATED(A%Left)
       AssocAR=ASSOCIATED(A%Right)
       ! Follow the left link (A)
       IF(AssocAL.AND.AssocAR)THEN
          CALL SetSRSTSkipPtrs(A%Left,A%Right)
       ELSEIF(AssocAL)THEN
          CALL SetSRSTSkipPtrs(A%Left)
       ELSEIF(AssocAR)THEN
          CALL SetSRSTSkipPtrs(A%Right)
       ENDIF
    ENDIF
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
      IF(A%Alloc==ALLOCATED_SKIP_POINTERS_OFF)RETURN
      ! Start past head link
      C=>A%Next
      ! Go over row links 
      DO WHILE(ASSOCIATED(C))
         IF(ASSOCIATED(C%RowRoot))THEN 
            ! Unset skip pointers for the sparse row search tree
            CALL UnSetSRSTSkipPtrs(C%RowRoot)
         ENDIF
         C=>C%Next
      ENDDO
      A%Alloc=ALLOCATED_SKIP_POINTERS_OFF
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
!
!======================================================================
      SUBROUTINE Print_FASTMAT(A)
         TYPE(FASTMAT),POINTER :: A,P
!----------------------------------------------------------------------
         P=>A%Next
         DO WHILE(ASSOCIATED(P))
            WRITE(*,*)'===================== Row = ',P%Row,'====================='
            CALL Print_SRST(P%RowRoot)
            P=>P%Next
         ENDDO
       END SUBROUTINE Print_FASTMAT
!======================================================================
!
!======================================================================
    RECURSIVE SUBROUTINE Print_SRST(A)
      TYPE(SRST),POINTER :: A
!----------------------------------------------------------------------
      IF(ASSOCIATED(A%Left))CALL Print_SRST(A%Left)
      IF(ASSOCIATED(A%Right))CALL Print_SRST(A%Right)      
      IF(ASSOCIATED(A%Left).AND.ASSOCIATED(A%Right))THEN
         WRITE(*,73)A%Tier,A%Number,A%L,A%R,A%Left%Number,A%Right%Number
      ELSEIF(ASSOCIATED(A%Left))THEN
         WRITE(*,74)A%Tier,A%Number,A%L,A%R,A%Left%Number
      ELSEIF(ASSOCIATED(A%Right))THEN
         WRITE(*,77)A%Tier,A%Number,A%L,A%R,A%Right%Number
      ELSE
         WRITE(*,75)A%Tier,A%Number,A%L,A%R
      ENDIF
73    FORMAT('xT=',I4,', Num=',I4,', [',I4,',',I4,'], L#=',I4,'R#=',I4)           
77    FORMAT('xT=',I4,', Num=',I4,', [',I4,',',I4,'], R#=',I4)           
74    FORMAT('xT=',I4,', Num=',I4,', [',I4,',',I4,'], L#=',I4)           
75    FORMAT('xT=',I4,', Num=',I4,', [',I4,',',I4,']')           
    END SUBROUTINE Print_SRST


    RECURSIVE SUBROUTINE Print_SkipPtrsSRST(A)
      TYPE(SRST),POINTER :: A,C
        C=>A        
        DO WHILE(ASSOCIATED(C%Left).OR.ASSOCIATED(C%Right))
           IF(ASSOCIATED(C%Left))THEN
              IF(ASSOCIATED(C%Left).AND.ASSOCIATED(C%Right))THEN
                 WRITE(*,74)C%Tier,C%Number,C%L,C%R,C%Left%Number,C%Right%Number
              ELSE
                 WRITE(*,79)C%Tier,C%Number,C%L,C%R,C%Left%Number
              ENDIF
              C=>C%Left
           ELSEIF(ASSOCIATED(C%Right))THEN
              IF(ASSOCIATED(C%Left).AND.ASSOCIATED(C%Right))THEN
                 WRITE(*,74)C%Tier,C%Number,C%L,C%R,C%Left%Number,C%Right%Number
              ELSE
                 WRITE(*,77)C%Tier,C%Number,C%L,C%R,C%Right%Number
              ENDIF
              C=>C%Right
           ENDIF
           IF(C%L==C%R)THEN
              IF(ASSOCIATED(C%Right))THEN
                 WRITE(*,73)C%Tier,C%Number,C%L,C%R,C%Right%Number
              ELSE
                 WRITE(*,72)C%Tier,C%Number,C%L,C%R
              ENDIF
           ENDIF              
72         FORMAT('T=',I2,', Num=',I3,', [',I3,',',I3,']')           
73         FORMAT('T=',I2,', Num=',I3,', [',I3,',',I3,'] R#=',I4)           
77         FORMAT('T=',I2,', Num=',I3,', [',I3,',',I3,'] R#=',I4)           
79         FORMAT('T=',I2,', Num=',I3,', [',I3,',',I3,'] L#=',I4)           
74         FORMAT('T=',I2,', Num=',I3,', [',I3,',',I3,'] L#=',I4,' R#=',I4)           
        ENDDO
        RETURN
    END SUBROUTINE Print_SkipPtrsSRST

#endif
  END MODULE FASTMATRICES
