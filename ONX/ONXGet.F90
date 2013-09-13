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
!  PUBLIC :: Reduce_FASTMAT
  PUBLIC :: GetDab
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
       CALL New_FASTMAT(B,0,(/0,0/),NSMat_O=A%NSMat)
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
             N = BSiz%I(J)*A%NSMat
             ! Add bloks to this sparse row search tree
             CALL AddFASTMATBlok(C,IRow,J,M,N,RESHAPE(A%MTrix%D(P:P+M*N-1),(/M,N/)))
          ENDDO
       ENDIF
    ENDDO
    !
    CALL FlattenAllRows(B)
    !
  END SUBROUTINE Set_DFASTMAT_EQ_DBCSR2
#endif
  !
#ifdef BLABLABLA
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
#ifdef ONX2_PARALLEL
  SUBROUTINE GetDab(DFMab,APt,ANbr,BPt,BNbr,Args)
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(FASTMAT), POINTER    :: DFMab
    TYPE(INT_VECT)            :: APt,BPt
    INTEGER      , INTENT(IN) :: ANbr,BNBr
    TYPE(ARGMT)               :: Args
    !-------------------------------------------------------------------
    TYPE(BCSR )     :: A
    TYPE(DBCSR)     :: B
    TYPE(INT_VECT)  :: ANBrArr,BNBrArr,APtArr,BPtArr,ADisp,BDisp
    INTEGER         :: iPrc,ANbrTot,BNbrTot,NSMat
    INTEGER         :: I,J,JG,M,MN,MN1,P,NAtms,NBlks,NNon0,ICol,IRow
    LOGICAL         :: ReAllocate
    !-------------------------------------------------------------------
    INTEGER, EXTERNAL :: IBinSrch
    !-------------------------------------------------------------------
    INTEGER :: TotBlk
    !
    !Send ANbr,BNbr to ROOT
    CALL New(ANBrArr,NPrc)
    CALL New(BNBrArr,NPrc)
    CALL Gather(ANBr,ANBrArr)
    CALL Gather(BNBr,BNBrArr)
    !
    IF(MyID==0) THEN
       !write(*,*) 'ANBrArr',ANBrArr%I
       !write(*,*) 'BNBrArr',BNBrArr%I
    ENDIF
    !
    ! Compute some indicies.
    IF(MyID.EQ.ROOT) THEN
       CALL New(ADisp,NPrc+1)
       CALL New(BDisp,NPrc+1)
       ANbrTot=ANBrArr%I(1)
       BNbrTot=BNBrArr%I(1)
       ADisp%I(1)=0
       BDisp%I(1)=0
       DO iPrc=2,NPrc
          ADisp%I(iPrc)=ANBrArr%I(iPrc-1)+ADisp%I(iPrc-1)
          BDisp%I(iPrc)=BNBrArr%I(iPrc-1)+BDisp%I(iPrc-1)
          ANbrTot=ANbrTot+ANBrArr%I(iPrc)
          BNbrTot=BNbrTot+BNBrArr%I(iPrc)
       ENDDO
       ADisp%I(NPrc+1)=ANbrTot
       BDisp%I(NPrc+1)=BNbrTot
       !write(*,*) 'ANbrTot',ANbrTot
       !write(*,*) 'ADisp',ADisp%I
       !write(*,*) 'BNbrTot',BNbrTot
       !write(*,*) 'BDisp',BDisp%I
       CALL New(APtArr,ANbrTot)
       CALL New(BPtArr,BNbrTot)
    ENDIF
    !
    !Alloc the arrays A,B on ROOT to get the atom lists
    CALL Gather(APt,APtArr,ANBr,ANBrArr,ADisp)
    CALL Gather(BPt,BPtArr,BNBr,BNBrArr,BDisp)
    !
    !IF(MyId==ROOT) write(*,*) 'APtArr',APtArr%I
    !IF(MyId==ROOT) write(*,*) 'BPtArr',BPtArr%I
    !
    ! Get the B matrix from disc.
    CALL Get_BCSR(A,TrixFile('D',Args,0))
    !
    !CALL PChkSum(A,'Density matrix on root',Unit_O=6)
    NSMat=A%NSMat
    CALL BCast(NSMat)
    !
    CALL New(B,NSMat_O=NSMat)
    !
!------------------------------------------------
!        Distribute to each processor
!
    DO iPrc=NPrc-1,0,-1
       IF(MyId==ROOT)THEN
          B%NAtms=0
          B%NBlks=1
          B%NNon0=1
          B%RowPt%I(1)=1
          !
          !vw beg and end must be set somewhere else.
          IRow=1
          DO I=1,NAtoms !Row min, Row max.
             ! Is it the right row?
             IF(IBinSrch(APtArr%I(ADisp%I(iPrc+1)+1),I,ADisp%I(iPrc+2)-ADisp%I(iPrc+1)).LT.1) CYCLE
             !if(iPrc==1) write(*,*) 'I',I
             ICol=1
             M=BSiz%I(I)
             DO J=A%RowPt%I(I),A%RowPt%I(I+1)-1
                JG=A%ColPt%I(J)
                ! Is it the right col?
                IF(IBinSrch(BPtArr%I(BDisp%I(iPrc+1)+1),JG,BDisp%I(iPrc+2)-BDisp%I(iPrc+1)).LT.1) CYCLE
                !if(iPrc==1) write(*,*) 'JG',JG
                B%ColPt%I(B%NBlks)=JG
                B%BlkPt%I(B%NBlks)=B%NNon0
                B%NBlks=B%NBlks+1
                MN=M*BSiz%I(JG)*NSMat
                MN1=MN-1
                P=A%BlkPt%I(J)
                B%MTrix%D(B%NNon0:B%NNon0+MN1)=A%MTrix%D(P:P+MN1)
                B%NNon0=B%NNon0+MN
             ENDDO
             B%NAtms=B%NAtms+1
             B%RowPt%I(B%NAtms+1)=B%NBlks
          ENDDO
          B%NBlks=B%NBlks-1
          B%NNon0=B%NNon0-1
! MINUS
          IF(iPrc/=ROOT)THEN
             CALL Send(B%NAtms,iPrc,1)
             CALL Send(B%NBlks,iPrc,2)
             CALL Send(B%NNon0,iPrc,3)
             CALL Send(B%RowPt,B%NAtms+1,iPrc,4)
             CALL Send(B%ColPt,B%NBlks,iPrc,5)
             CALL Send(B%BlkPt,B%NBlks,iPrc,6)
             CALL Send(B%MTrix,B%NNon0,iPrc,7)
          ENDIF
       ELSEIF(MyId==iPrc)THEN
          CALL Recv(B%NAtms,ROOT,1)
          CALL Recv(B%NBlks,ROOT,2)
          CALL Recv(B%NNon0,ROOT,3)
          ReAllocate=(SIZE(B%RowPt%I)<B%NAtms+1).OR. &
               (SIZE(B%ColPt%I)<B%NBlks)  .OR. &
               (SIZE(B%BlkPt%I)<B%NBlks)  .OR. &
               (SIZE(B%MTrix%D)<B%NNon0)
          IF(ReAllocate)THEN
             NAtms=B%NAtms; NBlks=B%NBlks; NNon0=B%NNon0
             CALL Delete(B)
             CALL New(B,(/NAtms,NBlks,NNon0/))
          ENDIF
          CALL Recv(B%RowPt,B%NAtms+1,ROOT,4)
          CALL Recv(B%ColPt,B%NBlks,ROOT,5)
          CALL Recv(B%BlkPt,B%NBlks,ROOT,6)
          CALL Recv(B%MTrix,B%NNon0,ROOT,7)
       ENDIF
    ENDDO
    !------------------------------------------------
    CALL BCast(A%NBlks)
    IF(MyId==ROOT)THEN
       CALL SetEq(B%GRwPt,A%RowPt,NAtoms+1)
       CALL SetEq(B%GClPt,A%ColPt,A%NBlks)
    ENDIF
    CALL BCast(B%GRwPt,NAtoms+1)
    CALL BCast(B%GClPt,A%NBlks)
    B%Node=MyId
    B%GUpDate=STATUS_TRUE
    !------------------------------------------------
    !
    !IF(MyID==1)write(*,*) 'B MyID=',MyID,B%MTrix%D
    !IF(MyID==1)write(*,*) 'B MyID=',MyID,B%RowPt%I
    !IF(MyID==1)write(*,*) 'B MyID=',MyID,B%ColPt%I
    !
    IF(MyID.EQ.ROOT) TotBlk=A%NBlks
    CALL BCast(TotBlk)
    !
    WRITE(*,100) DBLE(B%NBlks)/DBLE(TotBlk)*100d0,DBLE(B%NBlks)/(DBLE(NAtoms)**2)*100d0,MyID
100 FORMAT(' Remaining block, Rel=',F8.3,'% Abs=',F8.3,'% MyID=',I4)
    !
    ! Copy the DBCSR to a FastMat.
    CALL New_FASTMAT(DFMab,0,(/0,0/),NSMat_O=B%NSMat)
    CALL Set_DFASTMAT_EQ_DBCSR2(DFMab,B,APt,ANBr)
    !CALL PChkSum_FASTMAT2(DFMab,'Density matrix ab',Unit_O=6)
    !
    ! Delete arrays.
    IF(MyID.EQ.ROOT) THEN
       CALL Delete(APtArr)
       CALL Delete(BPtArr)
       CALL Delete(A)
    ENDIF
    !
    CALL Delete(B)
    CALL Delete(ANBrArr)
    CALL Delete(BNBrArr)
    !
  END SUBROUTINE GetDab
  !
  !
!!$  SUBROUTINE Set_DFASTMAT_EQ_DBCSR2(B,A,RowPt,RowNbr)
!!$    TYPE(FASTMAT),POINTER :: B,C
!!$    TYPE(SRST),POINTER    :: S
!!$    TYPE(DBCSR)           :: A
!!$    TYPE(INT_VECT)        :: RowPt
!!$    INTEGER, INTENT(IN)   :: RowNbr
!!$    INTEGER               :: I,IRow,J,JP,P,N,M
!!$    !---------------------------------------------------------------------
!!$    !
!!$    ! Set some pointers.
!!$    NULLIFY(C,S)
!!$    !
!!$    ! Check for prior allocation
!!$    IF(ASSOCIATED(B))THEN
!!$       CALL Delete_FASTMAT1(B)
!!$       ! Begin with a new header Node
!!$       CALL New_FASTMAT(B,0,(/0,0/))
!!$    ENDIF
!!$    !
!!$    DO I = 1,RowNbr
!!$       IRow=RowPt%I(I)
!!$       M = BSiz%I(IRow)
!!$       IF(A%RowPt%I(I+1)-A%RowPt%I(I)>1) THEN
!!$          ! Set current row link
!!$          C => FindFastMatRow_1(B,IRow)
!!$          DO JP = A%RowPt%I(I),A%RowPt%I(I+1)-1
!!$             J = A%ColPt%I(JP)
!!$             P = A%BlkPt%I(JP)
!!$             N = BSiz%I(J)
!!$             ! Add bloks to this sparse row search tree
!!$             CALL AddFASTMATBlok(C,IRow,J,RESHAPE(A%MTrix%D(P:P+M*N-1),(/M,N/)))
!!$          ENDDO
!!$       ENDIF
!!$    ENDDO
!!$    !
!!$    CALL FlattenAllRows(B)
!!$    !
!!$  END SUBROUTINE Set_DFASTMAT_EQ_DBCSR2
  !
#endif
END MODULE ONXGet
