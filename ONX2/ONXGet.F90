MODULE ONXGet
!H=================================================================================
!H MODULE ONXGet
!H This MODULE contains:
!H PUBLIC:
!H  o SUB GetIntCode
!H  o SUB GetIntSpace
!H  o SUB GetAdrB
!H  o SUB GetSubBlk
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
  PUBLIC :: GetIntCode
  PUBLIC :: GetIntSpace
  PUBLIC :: GetAdrB
  PUBLIC :: GetSubBlk
  PUBLIC :: GetOffArr
#ifdef ONX2_PARALLEL
  PUBLIC :: Get_Essential_RowCol
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
  SUBROUTINE GetIntCode(LTot,TBra,TKet,IntCodeV,IntCodeC,Explicit)
    IMPLICIT NONE
    INTEGER           :: LTot,TBra,TKet,IntCodeV,IntCodeC
    INTEGER           :: TB,TK
    LOGICAL           :: Explicit

    IF(LTot.GT.8) THEN
       Explicit=.FALSE.
       RETURN
    ELSE
       Explicit=.TRUE.
    ENDIF
    TB=TBra
    TK=TKet

    IF(TBra.eq.0301) TB=0201
    IF(TBra.eq.0303) TB=0202
    IF(TBra.eq.0302) TB=0202
    IF(TBra.eq.0603) TB=0602
    IF(TKet.eq.0301) TK=0201
    IF(TKet.eq.0303) TK=0202
    IF(TKet.eq.0302) TK=0202
    IF(TKet.eq.0603) TK=0602

    IntCodeC = TBra*10000+TKet
    IntCodeV = TB*10000+TK

!
! L=1 exceptions
!
!    IF (IntCodeC.eq.03010101) Explicit=.FALSE.
!    IF (IntCodeC.eq.01010301) Explicit=.FALSE.
!
! L=2 exceptions
!
!    IF (IntCodeC.eq.02010301) Explicit=.FALSE.
!    IF (IntCodeC.eq.03010201) Explicit=.FALSE.
!    IF (IntCodeC.eq.03010301) Explicit=.FALSE.
!    IF (IntCodeC.eq.06010101) Explicit=.FALSE.
!    IF (IntCodeC.eq.01010601) Explicit=.FALSE.
!
! L=3 exceptions
!
!    IF (IntCodeC.eq.02020201) Explicit=.FALSE.
!    IF (IntCodeC.eq.02010202) Explicit=.FALSE.  
!    IF (IntCodeC.eq.06020101) Explicit=.FALSE.
!    IF (IntCodeC.eq.01010602) Explicit=.FALSE.
!    IF (IntCodeC.eq.06010201) Explicit=.FALSE.
!    IF (IntCodeC.eq.02010601) Explicit=.FALSE.
!    IF (IntCodeC.eq.03020301) Explicit=.FALSE.  ! 02020201
!    IF (IntCodeC.eq.03030301) Explicit=.FALSE.  ! 02020201
!    IF (IntCodeC.eq.03020201) Explicit=.FALSE.  ! 02020201
!    IF (IntCodeC.eq.03030201) Explicit=.FALSE.  ! 02020201
!    IF (IntCodeC.eq.02020301) Explicit=.FALSE.  ! 02020201
!    IF (IntCodeC.eq.03010302) Explicit=.FALSE.  ! 02010202
!    IF (IntCodeC.eq.02010303) Explicit=.FALSE.  ! 02010202
!    IF (IntCodeC.eq.03010303) Explicit=.FALSE.  ! 02010202
!    IF (IntCodeC.eq.03010202) Explicit=.FALSE.  ! 02010202
!    IF (IntCodeC.eq.02010302) Explicit=.FALSE.  ! 02010202
!    IF (IntCodeC.eq.06030101) Explicit=.FALSE.  ! 06020101
!    IF (IntCodeC.eq.01010603) Explicit=.FALSE.  ! 01010602
!    IF (IntCodeC.eq.06010301) Explicit=.FALSE.  ! 06010201
!    IF (IntCodeC.eq.03010601) Explicit=.FALSE.  ! 02010601 
!
! L=4 exceptions
!
!    IF (IntCodeC.eq.02020202) Explicit=.FALSE.
!    IF (IntCodeC.eq.06060101) Explicit=.FALSE.
!    IF (IntCodeC.eq.01010606) Explicit=.FALSE.
!    IF (IntCodeC.eq.06010601) Explicit=.FALSE.
!    IF (IntCodeC.eq.06020201) Explicit=.FALSE.
!    IF (IntCodeC.eq.02010602) Explicit=.FALSE.
!    IF (IntCodeC.eq.06010202) Explicit=.FALSE.
!    IF (IntCodeC.eq.02020601) Explicit=.FALSE.
!    IF (IntCodeC.eq.03020202) Explicit=.FALSE.  ! 02020202
!    IF (IntCodeC.eq.02020302) Explicit=.FALSE.  ! 02020202
!    IF (IntCodeC.eq.03030202) Explicit=.FALSE.  ! 02020202
!    IF (IntCodeC.eq.02020303) Explicit=.FALSE.  ! 02020202
!    IF (IntCodeC.eq.03020302) Explicit=.FALSE.  ! 02020202
!    IF (IntCodeC.eq.03030302) Explicit=.FALSE.  ! 02020202
!    IF (IntCodeC.eq.03020303) Explicit=.FALSE.  ! 02020202
!    IF (IntCodeC.eq.03030303) Explicit=.FALSE.  ! 02020202
!    IF (IntCodeC.eq.03010602) Explicit=.FALSE.  ! 02010602
!    IF (IntCodeC.eq.02010603) Explicit=.FALSE.  ! 02010602
!    IF (IntCodeC.eq.03010603) Explicit=.FALSE.  ! 02010602
!    IF (IntCodeC.eq.06020301) Explicit=.FALSE.  ! 06020201
!    IF (IntCodeC.eq.06030201) Explicit=.FALSE.  ! 06020201
!    IF (IntCodeC.eq.06030301) Explicit=.FALSE.  ! 06020201
!    IF (IntCodeC.eq.06010302) Explicit=.FALSE.  ! 06010202
!    IF (IntCodeC.eq.06010303) Explicit=.FALSE.  ! 06010202
!    IF (IntCodeC.eq.03020601) Explicit=.FALSE.  ! 02020601
!    IF (IntCodeC.eq.03030601) Explicit=.FALSE.  ! 02020601
  END SUBROUTINE GetIntCode
  !
  !
  SUBROUTINE GetIntSpace(TBra,TKet,LBra,LKet,IS)
    IMPLICIT NONE
    INTEGER           :: TBra,TKet
    INTEGER           :: LBra,LKet
    TYPE(ISpc)        :: IS
    INTEGER           :: IndV,IndAC,IndBD
    INTEGER           :: ITypeA,ITypeB,ITypeC,ITypeD
    INTEGER,PARAMETER :: Space1(0:4) = (/1,4,17,24,36/)
    INTEGER,PARAMETER :: Space2(36)  = (/ 1,   4,    3,    0,    0,    6,  &
                                          4,  17,   12,    0,    0,   22,  &
                                          3,  12,    9,    0,    0,   16,  &
                                          0,   0,    0,    0,    0,    0,  &
                                          0,   0,    0,    0,    0,    0,  &
                                          6,  22,   16,    0,    0,   31   /)
    INTEGER,PARAMETER :: SpaceV(25)  = (/ 1,   4,   17,   24,   48, &
                                          4,  16,   68,   97,  197, &
                                         17,  68,  289,  408,  816, &
                                         24,  88,  408,  526,  983, &
                                         48, 140,  816,  863, 1708  /)
    INTEGER,PARAMETER :: SpaceF(28)  = (/ 1,                  &
                                          4, 3,               &
                                         10, 9, 6,            &
                                         20,19,16,10,         &
                                         35,34,31,25,15,      &
                                         56,55,52,46,36,21,   &
                                         84,83,80,74,64,49,28 /)

    IF(LBra>4.OR.LKet>4) THEN
      CALL Halt('Illegal LBra or LKet in GetIntSpace')
    ENDIF

    ITypeA = MOD(TBra,100)
    ITypeC = (TBra-ITypeA)/100
    ITypeB = MOD(TKet,100)
    ITypeD = (TKet-ITypeB)/100

    IndV = (LBra+1)+LKet*5
    IndAC = ITypeA+(ITypeC-1)*6
    IndBD = ITypeB+(ITypeD-1)*6

    IS%NB1  = Space1(LBra)
    IS%NK1  = Space1(LKet)
    IS%NB2  = Space2(IndAC)
    IS%NK2  = Space2(IndBD)  !vw seem to be not used.
    IS%L1   = SpaceF(ITypeC)
    IS%L2   = SpaceF(ITypeA)
    IS%L3   = SpaceF(ITypeD)
    IS%L4   = SpaceF(ITypeB)
    IS%NVRR = SpaceV(IndV)
    !
  END SUBROUTINE GetIntSpace
  !
  !
  SUBROUTINE GetSubBlk(NAR,NAC,NBR,NBC,RS,CS,A,B)
    IMPLICIT REAL*8 (a-h,o-z)
    IMPLICIT INTEGER (i-n)
    REAL*8 A(NAR,NAC)
    REAL*8 B(NBR,NBC)
    INTEGER NAR,NAC
    INTEGER NBR,NBC,RS,CS
    INTEGER I,J,RS1,CS1
    RS1=RS-1
    CS1=CS-1
    DO J=1,NBC
       DO I=1,NBR
          B(I,J)=A(I+RS1,J+CS1)
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE GetSubBlk
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
#endif
  !
END MODULE ONXGet
