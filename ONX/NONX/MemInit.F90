SUBROUTINE MemInit(DB,IB,BSc,BSp)
  USE InOut
  USE ONXMemory
  IMPLICIT NONE
  TYPE(DBuf),INTENT(INOUT) :: DB
  TYPE(IBuf),INTENT(INOUT) :: IB
  TYPE(BSet),INTENT(IN)    :: BSc,BSp
  IF (ErrorCode==eInit) THEN
    DB%MAXD    = 100
    DB%MAXC    = 10
    DB%MAXT    = 10
    DB%MAXP    = 8
    DB%NShells = 20   !  fix this
    DB%NTypes  = 20   !  fix this 
    DB%NCnts   = 20   !  fix this
    IB%MAXI    = 1000
    IB%NPrim   = MAX(BSc%NPrim,BSp%NPrim)
    IB%LngCC   = 60000
    CALL VRRLng(IB%LngVRR,IB%LngLoc)
    CALL New(DB)
    CALL New(IB)
    CALL CCDriver(IB%CDrv%I(1),IB%LngDrv)
    CALL VRRDriver(IB%VLOC%I(1),IB%SLOC%I(1,1),IB%LngVRR,IB%LngLoc)
  ELSEIF (ErrorCode==eMAXI) THEN
    CALL Delete(IB)
    IB%MAXI = IB%MAXI * 5
    CALL New(IB)
    CALL CCDriver(IB%CDrv%I(1),IB%LngDrv)
    CALL VRRDriver(IB%VLOC%I(1),IB%SLOC%I(1,1),IB%LngVRR,IB%LngLoc)
  ELSEIF (ErrorCode==eMAXD) THEN
    CALL Delete(DB)
    DB%MAXD = DB%MAXD * 2
    CALL New(DB)
  ELSEIF (ErrorCode==eMAXC) THEN
    CALL Delete(DB)
    DB%MAXC = DB%MAXC + 1
    CALL New(DB)
  ELSEIF (ErrorCode==eMAXT) THEN
    CALL Delete(DB)
    DB%MAXT = DB%MAXT + 1
    CALL New(DB)
  ELSE
    write(*,*) "Error code = ",ErrorCode
    CALL Halt(' Unknown error code in ONX')
  END IF

END SUBROUTINE MemInit
