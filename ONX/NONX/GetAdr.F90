SUBROUTINE GetAdrB(I,J,Ind,A,E_O)
!----------------------------------------------------------------------
! Finds the CSR address corresponding to the dense matrix index pair.
! Assumes that the collumn index array is in order (binary search)
!
! Inputs:
!   I     = row index of dense matrix
!   J     = collumn index of dense matrix
!   IType = 0 halt if element not found
!         = 1 return index 0 if not found
!   A     = A sparse matrix
!   E_O   = 0 call Halt if I cant find the address
!         = 1 return Ind=0 if I cant find the address
!
! Outputs:
!   Ind   = sparse matrix index corresponding to (I,J)
!
!----------------------------------------------------------------------
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
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
  1000 CONTINUE
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
  2000 CONTINUE
  Ind=iMid

  CASE (1)
  3000 CONTINUE
    iMid=(iEnd+iBgn)/2
    IF (A%ColPt%I(iMid)==J) GOTO 4000
    IF (iBgn>=iEnd) RETURN
    IF (A%ColPt%I(iMid)>J) THEN
      iEnd=iMid-1
    ELSE
      iBgn=iMid+1
    ENDIF
    GOTO 3000
  4000 CONTINUE
  Ind=iMid

  CASE DEFAULT
    CALL Halt(' Illegal switch in ONX:GetAdrB')
  END SELECT
END SUBROUTINE GetAdrB

