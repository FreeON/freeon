  SUBROUTINE RangeOfDensity(D,NameBuf,BfnInd,DB,BS,GM)
    USE DerivedTypes
    USE GlobalScalars
    USE Macros
    USE ONXParameters
    IMPLICIT NONE
#ifdef PARALLEL
    TYPE(DBCSR),INTENT(IN)        :: D
#else
    TYPE(BCSR),INTENT(IN)         :: D
#endif
    TYPE(INT_VECT),INTENT(INOUT)  :: NameBuf
    TYPE(INT_RNK2),INTENT(INOUT)  :: BfnInd
    TYPE(BSET),INTENT(IN)         :: BS
    TYPE(CRDS),INTENT(IN)         :: GM
    TYPE(DBuf)                    :: DB
    INTEGER                       :: AtA,ShellA,KA,CFA
    INTEGER                       :: AtB
    INTEGER                       :: ri,ci
    NameBuf%I=0
#ifdef PARALLEL
  NameBuf%I=0
  DO AtA=Beg%I(MyID),End%I(MyID)
    ri=AtA-Beg%I(MyID)+1
    NameBuf%I(AtA)=1
    DO ci=D%RowPt%I(ri),D%RowPt%I(ri+1)-1
      AtB=D%ColPt%I(ci)
      NameBuf%I(AtB)=1
    END DO
  END DO
#else
  NameBuf%I=0
  DO AtA=1,NAtoms
    NameBuf%I(AtA)=1
    DO ci=D%RowPt%I(AtA),D%RowPt%I(AtA+1)-1
      AtB=D%ColPt%I(ci)
      NameBuf%I(AtB)=1
    END DO
  END DO
  ShellA=0
  DO AtA=1,NAtoms
    KA=GM%AtTyp%I(AtA)
    DO CFA=1,BS%NCFnc%I(KA)
      ShellA=ShellA+1
      BfnInd%I(AtA,CFA)=ShellA
    END DO
  END DO 
  DB%NShells=ShellA
#endif
  END SUBROUTINE RangeOfDensity

