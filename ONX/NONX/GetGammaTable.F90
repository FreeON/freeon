SUBROUTINE GetGammaTable(L,IB)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE Thresholding
  USE ONXParameters
  USE ONXMemory
  USE Stats
  TYPE(IBuf),INTENT(INOUT)   :: IB
  INTEGER,INTENT(IN)         :: L

  IF (L/=IB%Lval) THEN
    CALL GammaTable(L,IB%Mesh,IB%FAsymp,IB%GT%D(0,0))
    IB%Lval=L
  ENDIF

END SUBROUTINE GetGammaTable



