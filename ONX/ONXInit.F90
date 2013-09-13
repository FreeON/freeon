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
MODULE ONXInit
!H=================================================================================
!H MODULE ONXInit
!H This MODULE contains:
!H  PUBLIC:
!H  o MOD PROC InitK
!H  o SUB InitSubInd
!H  o SUB InitBfnInd
!H  PRIVATE:
!H  o SUB InitK_BCSR
!H
!H  Comments:
!H
!H---------------------------------------------------------------------------------
  !
#ifndef PARALLEL
#undef ONX2_PARALLEL
#endif
  !
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE PrettyPrint
  USE ONXParameters
  !
  IMPLICIT NONE
  PRIVATE
  !
!---------------------------------------------------------------------------------
! PUBLIC DECLARATIONS
!---------------------------------------------------------------------------------
  PUBLIC  :: InitK
  PUBLIC  :: InitSubInd
  PUBLIC  :: InitBfnInd
  !
!---------------------------------------------------------------------------------
! PRIVATE DECLARATIONS
!---------------------------------------------------------------------------------
  PRIVATE :: InitK_BCSR
  !
!---------------------------------------------------------------------------------
! MOD PROC DECLARATIONS
!---------------------------------------------------------------------------------
  INTERFACE InitK
     MODULE PROCEDURE InitK_BCSR
  END INTERFACE
  !
CONTAINS
  !
  SUBROUTINE InitK_BCSR(BS,GM,K)
!H---------------------------------------------------------------------------------
!H SUBROUTINE InitK_BCSR(BS,GM,K)
!H
!H---------------------------------------------------------------------------------
    TYPE(BSET), INTENT(IN   ) :: BS
    TYPE(CRDS), INTENT(IN   ) :: GM
    TYPE(BCSR), INTENT(INOUT) :: K
    INTEGER                   :: AtA,KA,NBFA
    INTEGER                   :: AtB,KB,NBFB
    REAL(DOUBLE)              :: Ax,Ay,Az,AB2
    REAL(DOUBLE)              :: Bx,By,Bz
    INTEGER                   :: J,iPnt
    !
    ! Set the number of Atoms in K.
    K%NAtms = NAtoms
    !
    j=1
    iPnt=1
    DO AtA=1,NAtoms
       KA=GM%AtTyp%I(AtA)
       NBFA=BS%BfKnd%I(KA)
       Ax=GM%Carts%D(1,AtA)
       Ay=GM%Carts%D(2,AtA)
       Az=GM%Carts%D(3,AtA)
       K%RowPt%I(AtA)=j
       DO AtB=1,NAtoms
          KB=GM%AtTyp%I(AtB)
          NBFB=BS%BfKnd%I(KB)
          Bx=GM%Carts%D(1,AtB)
          By=GM%Carts%D(2,AtB)
          Bz=GM%Carts%D(3,AtB)
          AB2= (Ax-Bx)*(Ax-Bx) + &
               (Ay-By)*(Ay-By) + &
               (Az-Bz)*(Az-Bz)
          IF (SQRT(AB2).LE.ONXRange) THEN
             K%ColPt%I(j)=AtB
             K%BlkPt%I(j)=iPnt
             j=j+1
             iPnt=iPnt+K%NSMat*NBFA*NBFB
          END IF
       END DO ! ci
    END DO ! ri
    K%RowPt%I(NRows+1)=j
  END SUBROUTINE InitK_BCSR
  !
  !
  SUBROUTINE InitSubInd(BS,GM,SubInd)
!H---------------------------------------------------------------------------------
!H SUBROUTINE InitSubInd(BS,GM,SubInd)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(BSET    ), INTENT(IN   ) :: BS
    TYPE(CRDS    ), INTENT(IN   ) :: GM
    TYPE(INT_RNK2), INTENT(INOUT) :: SubInd
    INTEGER                       :: IndexA1,IndexA2
    INTEGER                       :: AtA,KA,NBFA,CFA
    INTEGER                       :: StartLA,StopLA,StrideA
    !
    IndexA1=0
    DO AtA=1,NAtoms
       KA=GM%AtTyp%I(AtA)
       NBFA=BS%BfKnd%I(KA)
       IndexA2=0
       DO CFA=1,BS%NCFnc%I(KA)
          IndexA1=IndexA1+1
          StartLA=BS%LStrt%I(CFA,KA)
          StopLA=BS%LStop%I(CFA,KA)
          StrideA=StopLA-StartLA+1
          SubInd%I(1,IndexA1)=AtA
          SubInd%I(2,IndexA1)=NBFA
          SubInd%I(3,IndexA1)=IndexA2+1
          IndexA2=IndexA2+StrideA
       ENDDO
    ENDDO
    !
  END SUBROUTINE InitSubInd
  !
  !
#ifdef ONX2_PARALLEL
  SUBROUTINE InitBfnInd(DB,BS,GM,RCPtr,RCNbr,LBfnInd_O)
#else
  SUBROUTINE InitBfnInd(DB,BS,GM,BfnInd)
#endif
!H---------------------------------------------------------------------------------
!H SUBROUTINE InitBfnInd(DB,BS,GM,BfnInd)
!H
!H Comments:
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
#ifdef ONX2_PARALLEL
    TYPE(INT_VECT), INTENT(INOUT), OPTIONAL :: LBfnInd_O
    TYPE(INT_VECT), INTENT(INOUT) :: RCPtr
    INTEGER       , INTENT(IN   ) :: RCNbr
    INTEGER                       :: iRC
#else
    TYPE(INT_VECT), INTENT(INOUT) :: BfnInd                                            !vw new
#endif
    TYPE(BSET    ), INTENT(IN   ) :: BS
    TYPE(CRDS    ), INTENT(IN   ) :: GM
    TYPE(DBuf    ), INTENT(INOUT) :: DB
    INTEGER                       :: AtA,ShellA,KA,CFA
    !
    ShellA=0
#ifdef ONX2_PARALLEL
    IF(PRESENT(LBfnInd_O)) THEN
       DO iRC = 1,RCNbr
          AtA = RCPtr%I(iRC)
          LBfnInd_O%I(AtA) = ShellA
          KA=GM%AtTyp%I(AtA)
          DO CFA=1,BS%NCFnc%I(KA)
             ShellA=ShellA+1
          ENDDO
       ENDDO
    ELSE
       DO iRC = 1,RCNbr
          AtA = RCPtr%I(iRC)
          KA=GM%AtTyp%I(AtA)
          DO CFA=1,BS%NCFnc%I(KA)
             ShellA=ShellA+1
          ENDDO
       ENDDO
    ENDIF
#else
    DO AtA=1,NAtoms
       BfnInd%I(AtA) = ShellA
       KA=GM%AtTyp%I(AtA)
       DO CFA=1,BS%NCFnc%I(KA)
          ShellA=ShellA+1
       ENDDO
    ENDDO
#endif
    DB%NShells=ShellA
  END SUBROUTINE InitBfnInd
  !
END MODULE ONXInit
