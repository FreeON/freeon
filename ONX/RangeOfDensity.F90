#ifdef PARALLEL
SUBROUTINE RangeOfDensity(D,NameBuf)
  USE DerivedTypes
  USE GlobalScalars
  USE Macros
  USE ONXParameters
  IMPLICIT NONE
  TYPE(DBCSR),INTENT(IN)        :: D
  TYPE(INT_VECT),INTENT(INOUT)  :: NameBuf
  INTEGER                       :: AtA
  INTEGER                       :: AtB
  INTEGER                       :: ri,ci
  NameBuf%I=0
  DO AtA=Beg%I(MyID),End%I(MyID)
    ri=AtA-Beg%I(MyID)+1
    NameBuf%I(AtA)=1
    DO ci=D%RowPt%I(ri),D%RowPt%I(ri+1)-1
      AtB=D%ColPt%I(ci)
      NameBuf%I(AtB)=1
    END DO 
  END DO
END SUBROUTINE RangeOfDensity
#else
SUBROUTINE RangeOfDensity(D,NameBuf)
  USE DerivedTypes
  USE GlobalScalars
  USE Macros
  USE ONXParameters
  IMPLICIT NONE
  TYPE(BCSR),INTENT(IN)         :: D
  TYPE(INT_VECT),INTENT(INOUT)  :: NameBuf
  INTEGER                       :: AtA
  INTEGER                       :: AtB
  INTEGER                       :: ri,ci
  NameBuf%I=0
  DO AtA=1,NAtoms
    NameBuf%I(AtA)=1
    DO ci=D%RowPt%I(AtA),D%RowPt%I(AtA+1)-1
      AtB=D%ColPt%I(ci)
      NameBuf%I(AtB)=1
    END DO
  END DO
END SUBROUTINE RangeOfDensity
#endif
