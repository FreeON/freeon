MODULE GetTables
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE Thresholding
  USE ONXParameters
  USE ONXMemory
  USE Stats

  CONTAINS

  SUBROUTINE GetGammaTable(L,IB)
    TYPE(IBuf),INTENT(INOUT)   :: IB
    INTEGER,INTENT(IN)         :: L
    IF (L/=IB%Lval) THEN
      CALL GammaTable(L,IB%Mesh,IB%GT%D(0,0))
      IB%Lval=L
    ENDIF
  END SUBROUTINE GetGammaTable

  SUBROUTINE GetExpTable(IB)
    TYPE(IBuf),INTENT(INOUT)   :: IB
    CALL ExpTable(IB%Mesh,IB%ET%D(0,0))
  END SUBROUTINE GetExpTable

END MODULE GetTables

