MODULE ONXGet
!H=================================================================================
!H MODULE ONXGet
!H This MODULE contains:
!H PUBLIC:
!H  o SUB GetAdrB
!H  o SUB Get_Essential_RowCol
!H  o SUB 
!H PRIVATE:
!H
!H Comments:
!H
!H---------------------------------------------------------------------------------
  !
#ifndef PARALLEL
#undef ONX2_PARALLEL
#endif
  !
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE FastMatrices
  !
  USE GlobalCharacters
  USE InOut
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
!
  IMPLICIT NONE
!  PRIVATE
  !
!--------------------------------------------------------------------------------- 
! PUBLIC DECLARATIONS
!--------------------------------------------------------------------------------- 
  PUBLIC :: GetAdrB
  PUBLIC :: GetOffArr
#ifdef ONX2_PARALLEL
  PUBLIC :: Get_Essential_RowCol
  PUBLIC :: Set_DFASTMAT_EQ_DBCSR2
  PUBLIC :: Reduce_FASTMAT
#endif
  !
CONTAINS
  !
  SUBROUTINE GetAdrB(I,J,Ind,A,E_O)
!----------------------------------------------------------------------
!H Finds the CSR address corresponding to the dense matrix index pair.
!H Assumes that the collumn index array is in order (binary search)
!H Inputs:
!H   I     = row index of dense matrix
!H   J     = collumn index of dense matrix
!H   IType = 0 halt if element not found
!H         = 1 return index 0 if not found
!H   A     = A sparse matrix
!H   E_O   = 0 call Halt if I cant find the address
!H         = 1 return Ind=0 if I cant find the address
!H Outputs:
!H   Ind   = sparse matrix index corresponding to (I,J)
!----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN)          :: I,J
    INTEGER,INTENT(OUT)         :: Ind
    TYPE(BCSR),INTENT(IN)       :: A
    INTEGER,OPTIONAL,INTENT(IN) :: E_O
    INTEGER                     :: iMid,iBgn,iEnd,Err
    iBgn=A%RowPt%I(I)
    iEnd=A%RowPt%I(I+1)-1
    Ind=0
    Err=0
    IF (PRESENT(E_O)) Err=E_O
    SELECT CASE (Err)
    CASE (0)
1000   CONTINUE
       iMid=(iEnd+iBgn)/2
       IF (A%ColPt%I(iMid)==J) GOTO 2000
       IF (iBgn>=iEnd) THEN
          WRITE(*,*) "I=",I," J=",J
          CALL Halt(' Matrix element not found in ONX:GetAdrB')
       ENDIF
       IF (A%ColPt%I(iMid)>J) THEN
          iEnd=iMid-1
       ELSE
          iBgn=iMid+1
       ENDIF
       GOTO 1000
2000   CONTINUE
       Ind=iMid
    CASE (1)
3000   CONTINUE
       iMid=(iEnd+iBgn)/2
       IF (A%ColPt%I(iMid)==J) GOTO 4000
       IF (iBgn>=iEnd) RETURN
       IF (A%ColPt%I(iMid)>J) THEN
          iEnd=iMid-1
       ELSE
          iBgn=iMid+1
       ENDIF
       GOTO 3000
4000   CONTINUE
       Ind=iMid
    CASE DEFAULT
       CALL Halt(' Illegal switch in ONX:GetAdrB')
    END SELECT
  END SUBROUTINE GetAdrB
  !  
  !
  SUBROUTINE GetOffArr(OffArr,BS)
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(INT_RNK2) :: OffArr
    TYPE(BSET)     :: BS
    !-------------------------------------------------------------------
    INTEGER        :: K,I,Off
    !-------------------------------------------------------------------
    OffArr%I=BIG_INT
    DO K=1,BS%NKind
       Off=1
       DO I=1,BS%NCFnc%I(K)
          OffArr%I(I,K)=Off
          Off=Off+BS%LStop%I(I,K)-BS%LStrt%I(I,K)+1
       ENDDO
    ENDDO
  END SUBROUTINE GetOffArr
  !
  !
#ifdef ONX2_PARALLEL
  SUBROUTINE Get_Essential_RowCol(A,RowPt,NRow,RowMin,RowMax,ColPt,NCol,ColMin,ColMax)
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(FASTMAT ), POINTER       :: A
    TYPE(INT_VECT), INTENT(INOUT) :: RowPt,ColPt
    INTEGER       , INTENT(  OUT) :: NRow,NCol,RowMin,RowMax,ColMin,ColMax
    !-------------------------------------------------------------------
    TYPE(FASTMAT ), POINTER       :: P
    TYPE(SRST    ), POINTER       :: U
    TYPE(INT_VECT)                :: TmpRowPt,TmpColPt
    INTEGER                       :: AtA,Row,Col
    !-------------------------------------------------------------------
    ! Create temporary array.
    CALL New(TmpRowPt,NAtoms)
    CALL SetEq(TmpRowPt,0)
    CALL New(TmpColPt,NAtoms)
    CALL SetEq(TmpColPt,0)
    !
    IF(.NOT.ASSOCIATED(A)) CALL Halt(' A is null in PDrv_Collect, STOP. ')
    !
    ! Get the essential rows and columns.
    P => A%Next
    DO
       IF(.NOT.ASSOCIATED(P)) EXIT
       U => P%RowRoot
       Row = P%Row
       TmpRowPt%I(Row)=1
       DO
          IF(.NOT.ASSOCIATED(U)) EXIT
          IF(U%L == U%R) THEN
             Col = U%R
             TmpColPt%I(Col)=1
          ENDIF
          U => U%Next
       ENDDO
       P => P%Next
    ENDDO
    !
    ! Set the essentail number of rows and columns.
    NRow = SUM(TmpRowPt%I(:))
    NCol = SUM(TmpColPt%I(:))
    !
    ! Create the essential arrays.
    CALL New(RowPt,NRow)
    CALL SetEq(RowPt,0)
    CALL New(ColPt,NCol)
    CALL SetEq(ColPt,0)
    !
    ! Fill them.
    Row = 1
    Col = 1
    RowMin=0
    RowMax=0
    ColMin=0
    ColMax=0
    DO AtA = 1,NAtoms
       IF(TmpRowPt%I(AtA).EQ.1) THEN
          RowPt%I(Row) = AtA
          IF(RowMin.EQ.0) RowMin=AtA
          RowMax=AtA
          Row = Row+1
       ENDIF
       IF(TmpColPt%I(AtA).EQ.1) THEN
          ColPt%I(Col) = AtA
          IF(ColMin.EQ.0) ColMin=AtA
          ColMax=AtA
          Col = Col+1
       ENDIF
    ENDDO
    !
    ! Delete temorary arrays.
    CALL Delete(TmpRowPt)
    CALL Delete(TmpColPt)
    !
  END SUBROUTINE Get_Essential_RowCol
  !
  SUBROUTINE Set_DFASTMAT_EQ_DBCSR2(B,A,RowPt,RowNbr)
    TYPE(FASTMAT),POINTER :: B,C
    TYPE(SRST),POINTER    :: S
    TYPE(DBCSR)           :: A
    TYPE(INT_VECT)        :: RowPt
    INTEGER, INTENT(IN)   :: RowNbr
    INTEGER               :: I,IRow,J,JP,P,N,M
    !---------------------------------------------------------------------
    !
    ! Set some pointers.
    NULLIFY(C,S)
    !
    ! Check for prior allocation
    IF(ASSOCIATED(B))THEN
       CALL Delete_FASTMAT1(B)
       ! Begin with a new header Node
       CALL New_FASTMAT(B,0,(/0,0/))
    ENDIF
    !
    DO I = 1,RowNbr
       IRow=RowPt%I(I)
       M = BSiz%I(IRow)
       IF(A%RowPt%I(I+1)-A%RowPt%I(I)>1) THEN
          ! Set current row link
          C => FindFastMatRow_1(B,IRow)
          DO JP = A%RowPt%I(I),A%RowPt%I(I+1)-1
             J = A%ColPt%I(JP)
             P = A%BlkPt%I(JP)
             N = BSiz%I(J)
             ! Add bloks to this sparse row search tree
             CALL AddFASTMATBlok(C,IRow,J,RESHAPE(A%MTrix%D(P:P+M*N-1),(/M,N/)))
          ENDDO
       ENDIF
    ENDDO
    !
    CALL FlattenAllRows(B)
    !
  END SUBROUTINE Set_DFASTMAT_EQ_DBCSR2
#endif
  !
!#ifdef ONX2_PARALLEL
#ifdef 1
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
    CALL New_BCSR(A,OnAll_O=.TRUE.)
    CALL New_BCSR(B,OnAll_O=.TRUE.)
    CALL New_BCSR(C,OnAll_O=.TRUE.)
    CALL Set_LBCSR_EQ_DFASTMAT(A,AFM)
    !
    ! Run over the shells in the tree.
    DO L=LMax,1,-1
       ! Compute the number max of send/recive.
       NCom=2**(L-1)
       IF(MyID==0)WRITE(*,*) 'L=',L,' NCom=',NCom
       !
       To=0
       DO i=1,NCom
          From=To+2**(LMax-l);
          IF(MyID==0)WRITE(*,*) 'To =',To,'; From =',From,'; i =',i,MyID
          !
          ! Take into account If NPrc is not a power of 2.
          IF(From.GT.NPrc-1) THEN
             IF(MyID==0)WRITE(*,*) 'This one is grather than NPrc, we exit.'
             EXIT
          ENDIF
          !
          IF(MyID.EQ.To) THEN
             CALL Recv_LBCSR(B,From)
             CALL Add_LBCSR(A,B,C)
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
       CALL New(B,OnAll_O=.TRUE.)
       B%NAtms=A%NAtms; B%NBlks=A%NBlks; B%NNon0=A%NNon0
    ELSE
       IF(.NOT.AllocQ(B%Alloc))CALL New(B,OnAll_O=.TRUE.)
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
  
  SUBROUTINE Add_LBCSR(A,B,C)
    IMPLICIT NONE
    TYPE(BCSR), INTENT(INOUT) :: A,B
    TYPE(BCSR), INTENT(INOUT) :: C
    INTEGER                   :: Status
    REAL(DOUBLE)              :: FlOp
    IF(.NOT.AllocQ(C%Alloc)) CALL New(C,OnAll_O=.TRUE.)
    CALL New(Flag,NAtoms)
    CALL SetEq(Flag,0)
    Flop=Zero
    Status=Add_GENERIC(SIZE(C%ColPt%I),SIZE(C%MTrix%D),         &
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
    INTEGER               :: I,J,JP,P,N,M,Q,OI,OJ,IC,JC,At,OldR
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
                DO JC=1,N
                   DO IC=1,M
                      B%MTrix%D(Q+IC-1+(JC-1)*M) = U%MTrix(IC,JC)
                   ENDDO
                ENDDO
                !CALL DBL_VECT_EQ_DBL_VECT(M*N,B%MTrix%D(Q),U%MTrix(1,1))
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
    !write(*,*) 'In Set_BCSR_EQ_FASTMAT 9',MyID
    OldR=B%RowPt%I(1)
    DO At=2,NAtoms+1
       IF(B%RowPt%I(At).EQ.-100000) B%RowPt%I(At)=OldR
       OldR=B%RowPt%I(At)
    ENDDO
  END SUBROUTINE Set_LBCSR_EQ_DFASTMAT
#endif
END MODULE ONXGet
