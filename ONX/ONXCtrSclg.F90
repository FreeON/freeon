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
MODULE ONXCtrSclg
!H=================================================================================
!H MODULE ONXCtrSclg
!H This MODULE contains:
!H  PUBLIC:
!H  o MOD PRO TrnMatBlk
!H  PRIVATE:
!H  o SUB TrnMatBlk_BCSR
!H  o SUB TrnMatBlk_DBCSR
!H  o SUB TrnMatBlk_FASTMAT
!H  o SUB TrnSubBlk
!H  o FUN FacNrm
!H
!H Comments:
!H  - TrnMatBlk_DBCSR must be keeped, needed for Density matrix.
!H---------------------------------------------------------------------------------
  !
#ifndef PARALLEL
#undef ONX2_PARALLEL
#endif
  !
  USE DerivedTypes
  USE GlobalCharacters
  USE GlobalScalars
  USE GlobalObjects
  USE ProcessControl
  USE MemMan
  USE PrettyPrint
  USE ONXParameters
#ifdef ONX2_PARALLEL
  USE FastMatrices
#endif
  IMPLICIT NONE
  PRIVATE
!---------------------------------------------------------------------------------
! PUBLIC DECLARATIONS
!---------------------------------------------------------------------------------
  PUBLIC  :: TrnMatBlk
!---------------------------------------------------------------------------------
! PRIVATE DECLARATIONS
!---------------------------------------------------------------------------------
  PRIVATE :: TrnMatBlk_BCSR
#ifdef ONX2_PARALLEL
  PRIVATE :: TrnMatBlk_DBCSR
  PRIVATE :: TrnMatBlk_FASTMAT
#endif
  PRIVATE :: TrnSubBlk
  PRIVATE :: FacNrm
!---------------------------------------------------------------------------------
! MODULE PROCEDURE DECLARATIONS
!---------------------------------------------------------------------------------
  INTERFACE TrnMatBlk
     MODULE PROCEDURE TrnMatBlk_BCSR
#ifdef ONX2_PARALLEL
     MODULE PROCEDURE TrnMatBlk_DBCSR
     MODULE PROCEDURE TrnMatBlk_FASTMAT
#endif
  END INTERFACE
!
CONTAINS
!
#ifdef ONX2_PARALLEL
  SUBROUTINE TrnMatBlk_DBCSR(BS,GM,A)
!H---------------------------------------------------------------------------------
!H SUBROUTINE TrnMatBlk_DBCSR(BS,GM,A)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(BSET),INTENT(IN)     :: BS
    TYPE(CRDS),INTENT(IN)     :: GM
    TYPE(DBCSR),INTENT(INOUT) :: A
    INTEGER                   :: ri,ci,iPtr,N2
    INTEGER                   :: AtA,AtB,KA,KB,NBFA,NBFB
    DO AtA=Beg%I(MyId),End%I(MyId)
       ri=AtA-Beg%I(MyId)+1
       KA=GM%AtTyp%I(AtA)
       NBFA=BS%BfKnd%I(KA)
       DO ci=A%RowPt%I(ri),A%RowPt%I(ri+1)-1
          AtB=A%ColPt%I(ci)
          KB=GM%AtTyp%I(AtB)
          NBFB=BS%BfKnd%I(KB)
          N2=NBFA*NBFB
          iPtr=A%BlkPt%I(ci)
          CALL TrnSubBlk(BS,KA,KB,NBFA,NBFB,A%MTrix%D(iPtr:iPtr+N2-1))
       END DO
    END DO
  END SUBROUTINE TrnMatBlk_DBCSR

  SUBROUTINE TrnMatBlk_FASTMAT(BS,GM,A)
!H---------------------------------------------------------------------------------
!H SUBROUTINE TrnMatBlk_FASTMAT(BS,GM,A)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(BSET)   , INTENT(IN) :: BS
    TYPE(CRDS)   , INTENT(IN) :: GM
    TYPE(FASTMAT), POINTER    :: A,P
    TYPE(SRST   ), POINTER    :: U
    INTEGER                   :: Row,Col,N2
    INTEGER                   :: AtA,AtB,KA,KB,NBFA,NBFB,iSMat
!
    IF(.NOT.ASSOCIATED(A)) CALL Halt(' A not associated in TrnMatBlk_FASTMAT ')
    CALL FlattenAllRows(A)
!
    P => A%Next
    DO
       IF(.NOT. ASSOCIATED(P)) EXIT
       Row  = P%Row
       KA   = GM%AtTyp%I(Row)
       NBFA = BS%BfKnd%I(KA)
       U => P%RowRoot
       DO
          IF(.NOT. ASSOCIATED(U)) EXIT
          IF(U%L == U%R) THEN
             Col  = U%L
             KB   = GM%AtTyp%I(Col)
             NBFB = BS%BfKnd%I(KB)
             N2=NBFA*NBFB
             !if(myid==1)write(*,*) 'Row',Row,'Col',Col
             DO iSMat=1,A%NSMat
                CALL TrnSubBlk(BS,KA,KB,NBFA,NBFB,U%MTrix(1:NBFA, (iSMat-1)*NBFB+1:iSMat*NBFB  ))
             ENDDO
             !CALL TrnSubBlk(BS,KA,KB,NBFA,NBFB,U%MTrix)
          ENDIF
          U => U%Next
       ENDDO
       P => P%Next
    ENDDO
!
    CALL FlattenAllRows(A)
  END SUBROUTINE TrnMatBlk_FASTMAT
#endif
!
  SUBROUTINE TrnMatBlk_BCSR(BS,GM,A)
!H---------------------------------------------------------------------------------
!H SUBROUTINE TrnMatBlk_BCSR(BS,GM,A)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(BSET),INTENT(IN)    :: BS
    TYPE(CRDS),INTENT(IN)    :: GM
    TYPE(BCSR),INTENT(INOUT) :: A
    INTEGER                  :: ci,iPtr,N2,iSMat
    INTEGER                  :: AtA,AtB,KA,KB,NBFA,NBFB
!
    DO AtA=1,NAtoms
       KA=GM%AtTyp%I(AtA)
       NBFA=BS%BfKnd%I(KA)
       DO ci=A%RowPt%I(AtA),A%RowPt%I(AtA+1)-1
          AtB=A%ColPt%I(ci)
          KB=GM%AtTyp%I(AtB)
          NBFB=BS%BfKnd%I(KB)
          N2=NBFA*NBFB
          iPtr=A%BlkPt%I(ci)
          DO iSMat=1,A%NSMat
             CALL TrnSubBlk(BS,KA,KB,NBFA,NBFB,A%MTrix%D(iPtr:iPtr+N2-1))
             iPtr=iPtr+N2
          ENDDO
       END DO
    END DO
  END SUBROUTINE TrnMatBlk_BCSR
!
  SUBROUTINE TrnSubBlk(BS,KA,KB,NBFA,NBFB,A)
!H---------------------------------------------------------------------------------
!H SUBROUTINE TrnSubBlk(BS,KA,KB,NBFA,NBFB,A)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(BSET),INTENT(IN)              :: BS
    INTEGER,INTENT(IN)                 :: KA,KB,NBFA,NBFB
    REAL(DOUBLE),DIMENSION(NBFA*NBFB),INTENT(INOUT) :: A
    INTEGER                            :: Index
    INTEGER                            :: CFA,StartLA,StopLA,LMNA,LA,MA,NA
    INTEGER                            :: CFB,StartLB,StopLB,LMNB,LB,MB,NB
    REAL(DOUBLE)                       :: FacA,FacB
    Index=0
    DO CFB=1,BS%NCFnc%I(KB)
       StartLB=BS%LStrt%I(CFB,KB)
       StopLB=BS%LStop%I(CFB,KB)
       DO LMNB=StartLB,StopLB
          LB=BS%LxDex%I(LMNB)
          MB=BS%LyDex%I(LMNB)
          NB=BS%LzDex%I(LMNB)
          FacB=FacNrm(LB,MB,NB)
          DO CFA=1,BS%NCFnc%I(KA)
             StartLA=BS%LStrt%I(CFA,KA)
             StopLA=BS%LStop%I(CFA,KA)
             DO LMNA=StartLA,StopLA
                LA=BS%LxDex%I(LMNA)
                MA=BS%LyDex%I(LMNA)
                NA=BS%LzDex%I(LMNA)
                FacA=FacNrm(LA,MA,NA)
                Index=Index+1
                A(Index)=A(Index)*FacA*FacB
             END DO ! LMNA
          END DO ! CFA
       END DO ! LMNB
    END DO ! CFB
  END SUBROUTINE TrnSubBlk
!
  FUNCTION FacNrm(L,M,N)
!H---------------------------------------------------------------------------------
!H FUNCTION FacNrm(L,M,N)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(IN)     ::  L,M,N
    REAL(DOUBLE),PARAMETER :: Fct(0:8) = (/1D0,1D0,3D0,15D0,105D0,945D0, &
                                           10395D0,135135D0,2027025D0/)
    REAL(DOUBLE)           :: x,y,FacNrm
    X=Fct(L+M+N)
    Y=Fct(L)*Fct(M)*Fct(N)
    FacNrm=SQRT(X/Y)
  END FUNCTION FacNrm
!
END MODULE ONXCtrSclg
