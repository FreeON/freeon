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
  IMPLICIT NONE
  PRIVATE
  !
!--------------------------------------------------------------------------------- 
! PUBLIC DECLARATIONS
!--------------------------------------------------------------------------------- 
  PUBLIC :: GetAdrB
  PUBLIC :: GetOffArr
#ifdef ONX2_PARALLEL
  PUBLIC :: Get_Essential_RowCol
  PUBLIC :: Set_DFASTMAT_EQ_DBCSR2
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
END MODULE ONXGet
