SUBROUTINE Fillout_BCSR(BS,GM,A)
  USE DerivedTypes
  USE GlobalScalars
  IMPLICIT NONE
  TYPE(BSET),INTENT(IN)    :: BS
  TYPE(CRDS),INTENT(IN)    :: GM
  TYPE(BCSR),INTENT(INOUT) :: A
  INTEGER                  :: AtA,KA,NBFA
  INTEGER                  :: AtB,KB,NBFB
  INTEGER                  :: iPnt1,iPnt2,ci,Ind

#ifdef PARALLEL
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
          IF (iPnt1.EQ.iPnt2) THEN
            CALL XPose1C(NBFA,NBFB,A%MTrix%D(iPnt1),A%MTrix%D(ipnt2))
          ELSE
            CALL XPose2C(NBFA,NBFB,A%MTrix%D(iPnt1),A%MTrix%D(ipnt2))
          END IF
        END IF
      END DO
    END DO
#ifdef PARALLEL
  END IF
#endif
END SUBROUTINE Fillout_BCSR

#ifdef PARALLEL
SUBROUTINE Fillout_DBCSR(BS,GM,A,NameBuf)
  USE DerivedTypes
  USE GlobalScalars
  USE ONXParameters
  IMPLICIT NONE
  TYPE(BSET),INTENT(IN)     :: BS
  TYPE(CRDS),INTENT(IN)     :: GM
  TYPE(DBCSR),INTENT(INOUT) :: A
  TYPE(INT_VECT),INTENT(IN) :: NameBuf
  INTEGER                   :: AtA,KA,NBFA
  INTEGER                   :: AtB,KB,NBFB
  INTEGER                   :: iPnt1,iPnt2,ci,Ind,ri
  DO ri=1,NRows
    AtA=NameBuf%I(ri)
    KA=GM%AtTyp%I(AtA)
    NBFA=BS%BfKnd%I(KA)
    DO ci=A%RowPt%I(AtA),A%RowPt%I(AtA+1)-1
      AtB=A%ColPt%I(ci)
      KB=GM%AtTyp%I(AtB)
      NBFB=BS%BfKnd%I(KB)
      IF (AtA.GE.AtB) THEN
        CALL GetAdrB(AtB,AtA,Ind,0,A%ColPt%I,A%RowPt%I)
        iPnt1=A%BlkPt%I(ci)
        iPnt2=A%BlkPt%I(Ind)
        IF (iPnt1.EQ.iPnt2) THEN
          CALL XPose1C(NBFA,NBFB,A%MTrix%D(iPnt1),A%MTrix%D(ipnt2))
        ELSE
          CALL XPose2C(NBFA,NBFB,A%MTrix%D(iPnt1),A%MTrix%D(ipnt2))
        END IF
      END IF
    END DO
  END DO

END SUBROUTINE Fillout_DBCSR
#endif
