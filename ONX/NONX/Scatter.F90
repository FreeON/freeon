SUBROUTINE Scatter(N,NA,NB,IndexA,SB,SubInd,DB,KB,K)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE ONXParameters
  USE ONXMemory
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: N,NA,NB,IndexA
  TYPE(DSL),INTENT(IN)       :: SB
  TYPE(INT_RNK2),INTENT(IN)  :: SubInd
  TYPE(DBuf),INTENT(IN)      :: DB
  REAL(DOUBLE),INTENT(IN)    :: KB(N,NA,NB)
  TYPE(BCSR),INTENT(INOUT)   :: K
  INTEGER                    :: I,I0,Ind,IndexB,Ioff
  INTEGER                    :: AtA,NBFA,RS
  INTEGER                    :: AtB,NBFB,CS

  DO I=1,N
    I0     = SB%SLDis%I(I)-3
    IndexB = INT(ABS(DB%DisBuf%D(I0)))
    AtA    = SubInd%I(1,IndexA)
    NBFA   = SubInd%I(2,IndexA)
    RS     = SubInd%I(3,IndexA)
    AtB    = SubInd%I(1,IndexB)
    NBFB   = SubInd%I(2,IndexB)
    CS     = SubInd%I(3,IndexB)
    CALL GetAdrB(AtA,AtB,Ind,K,0)
    Ioff=K%BlkPt%I(Ind)
    IF (Ind > 0) CALL PutSubBlk(I,N,NBFA,NBFB,NA,NB,RS,CS, &
                                K%MTrix%D(Ioff),KB)
  END DO

END SUBROUTINE Scatter
