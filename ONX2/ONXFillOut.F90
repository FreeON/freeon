MODULE ONXFillOut
!H=================================================================================
!H MODULE ONXFillOut
!H This MODULE contains:
!H  PUBLIC:
!H  o MOD PROC FillOut
!H  PRIVATE:
!H  o SUB FillOutBCSR
!H  o SUB FillOutFASTMAT
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
  USE ONXParameters
  USE ONXGet
#ifdef ONX2_PARALLEL
  USE FastMatrices
  USE MondoMPI
#endif
  IMPLICIT NONE
  PRIVATE
  !
!--------------------------------------------------------------------------------- 
! PUBLIC DECLARATIONS
!--------------------------------------------------------------------------------- 
  PUBLIC  :: FillOutBCSR
!!$#ifdef ONX2_PARALLEL
!!$  PUBLIC  :: FillOutFASTMAT
!!$#endif
  !
!--------------------------------------------------------------------------------- 
! PRIVATE DECLARATIONS
!--------------------------------------------------------------------------------- 
  !
CONTAINS
  !
  SUBROUTINE FillOutBCSR(BS,GM,A)
!H---------------------------------------------------------------------------------
!H SUBROUTINE FillOutBCSR(BS,GM,A)
!H 
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(BSET),INTENT(IN)    :: BS
    TYPE(CRDS),INTENT(IN)    :: GM
    TYPE(BCSR),INTENT(INOUT) :: A
    INTEGER                  :: AtA,KA,NBFA
    INTEGER                  :: AtB,KB,NBFB
    INTEGER                  :: iPnt1,iPnt2,iPnt1T,iPnt2T,ci,Ind,iSmat
!
!call print_bcsr(A,'FillOut: A1',unit_o=6)
#ifdef ONX2_PARALLEL
    IF (MyID==ROOT) THEN
#endif
    DO AtA=1,NAtoms
       KA=GM%AtTyp%I(AtA)
       NBFA=BS%BfKnd%I(KA)
       DO ci=A%RowPt%I(AtA),A%RowPt%I(AtA+1)-1
          AtB=A%ColPt%I(ci)
          KB=GM%AtTyp%I(AtB)
          NBFB=BS%BfKnd%I(KB)
          IF (AtA.GE.AtB) THEN
             CALL GetAdrB(AtB,AtA,Ind,A,0)
             iPnt1=A%BlkPt%I(ci)
             iPnt2=A%BlkPt%I(Ind)
             IF(iPnt1.EQ.iPnt2) CYCLE
             IF(A%NSMat.EQ.1)THEN
                !IF (iPnt1.EQ.iPnt2) THEN
                !   CALL XPose1C(NBFA,NBFB,A%MTrix%D(iPnt1),A%MTrix%D(ipnt2))
                !ELSE
                   CALL XPose2C(NBFA,NBFB,A%MTrix%D(iPnt1),A%MTrix%D(iPnt2))
                !END IF
             ELSEIF(A%NSMat.EQ.2)THEN
                !DO iSMat=1,2
                   !IF (iPnt1.EQ.iPnt2) THEN
                   !   CALL XPose1C(NBFA,NBFB,A%MTrix%D(iPnt1),A%MTrix%D(ipnt2))
                   !ELSE
                   CALL XPose2C(NBFA,NBFB,A%MTrix%D(iPnt1),A%MTrix%D(iPnt2))
                   !END IF
                   iPnt1=iPnt1+NBFA*NBFB
                   iPnt2=iPnt2+NBFA*NBFB
                   CALL XPose2C(NBFA,NBFB,A%MTrix%D(iPnt1),A%MTrix%D(iPnt2))
                !ENDDO
             ELSEIF(A%NSMat.EQ.4)THEN
                iPnt1T=iPnt1
                iPnt2T=iPnt2
                !A_aa
                CALL XPose2C(NBFA,NBFB,A%MTrix%D(iPnt1T),A%MTrix%D(iPnt2T))
                !A_ab
                iPnt1T=iPnt1+  NBFA*NBFB
                iPnt2T=iPnt2+2*NBFA*NBFB
                CALL XPose2C(NBFA,NBFB,A%MTrix%D(iPnt1T),A%MTrix%D(iPnt2T))
                !A_ba
                iPnt1T=iPnt1+2*NBFA*NBFB
                iPnt2T=iPnt2+  NBFA*NBFB
                CALL XPose2C(NBFA,NBFB,A%MTrix%D(iPnt1T),A%MTrix%D(iPnt2T))
                !A_bb
                iPnt1T=iPnt1+3*NBFA*NBFB
                iPnt2T=iPnt2+3*NBFA*NBFB
                CALL XPose2C(NBFA,NBFB,A%MTrix%D(iPnt1T),A%MTrix%D(iPnt2T))
             ELSE
                CALL Halt('ONXFillOut: wrong value for NSMat!')
             ENDIF
          END IF
       END DO
    END DO
#ifdef ONX2_PARALLEL
    END IF
#endif
!call print_bcsr(A,'FillOut: A2',unit_o=6)
  END SUBROUTINE FillOutBCSR
  !
  !
!!$#ifdef ONX2_PARALLEL
!!$  SUBROUTINE FillOutFASTMAT(BS,GM,KFastMat)
!!$!H---------------------------------------------------------------------------------
!!$!H SUBROUTINE FillOutDBCSR(BS,GM,KFastMat)
!!$!H 
!!$!H---------------------------------------------------------------------------------
!!$    IMPLICIT NONE
!!$    TYPE(BSET),INTENT(IN)     :: BS
!!$    TYPE(CRDS),INTENT(IN)     :: GM
!!$    TYPE(FastMat),POINTER          :: KFastMat
!!$    !
!!$    CALL Symmetrized_FASMAT(KFastMat,'L')
!!$    !
!!$  END SUBROUTINE FillOutFASTMAT
!!$#endif
  !
END MODULE ONXFillOut
