  SUBROUTINE RangeOfDensity(D,NameBuf)
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
    INTEGER                       :: AtA
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
#endif
  END SUBROUTINE RangeOfDensity

