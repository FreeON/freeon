SUBROUTINE MemInit(DB,IB,MB,Drv,BSc,BSp)
  USE InOut
  USE ONXMemory
  IMPLICIT NONE
  TYPE(DBuf),INTENT(INOUT) :: DB
  TYPE(IBuf),INTENT(INOUT) :: IB
  TYPE(DML),INTENT(INOUT)  :: MB
  TYPE(IDrv),INTENT(INOUT) :: Drv
  TYPE(BSet),INTENT(IN)    :: BSc,BSp
  write(*,*) "In MemInit: ErrorCode=",ErrorCode
  IF (ErrorCode==eInit) THEN
    DB%MAXDis  = 1000
    DB%MAXPrm  = 3000
    DB%MAXD    = 100
    DB%MAXC    = 10
    DB%MAXT    = 10
    DB%MAXP    = 8
    DB%NShells = 20   !  fix this
    DB%NTypes  = 20   !  fix this 
    DB%NCnts   = 20   !  fix this
    DB%NPrim   = MAX(BSc%NPrim,BSp%NPrim)
    DB%MInfo   = 0    
    MB%MAXML   = 10000
    IB%MAXI    = 30000
    IB%NPrim   = MAX(BSc%NPrim,BSp%NPrim)
    IB%Lval    = -1
    Drv%LngCC  = 60000
    CALL GammaHeader(IB%Mesh,IB%Switch,IB%Grid)
    CALL VRRLng(Drv%LngVRR,Drv%LngLoc)
    CALL New(DB)
    CALL New(IB)
    CALL New(MB)
    CALL New(Drv)
    CALL CCDriver(Drv%CDrv%I(1),Drv%LngDrv)
    CALL VRRDriver(Drv%VLOC%I(1),Drv%SLOC%I(1,1),Drv%LngVRR,Drv%LngLoc)
  ELSEIF (ErrorCode==eMAXI) THEN
    CALL Delete(IB)
    IB%MAXI = IB%MAXI * 5
    CALL New(IB)
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
  ELSEIF (ErrorCode==eMAXDis) THEN
    CALL Delete(DB)
    DB%MAXDis = DB%MAXDis * 5
    CALL New(DB)
  ELSEIF (ErrorCode==eMAXPrm) THEN
    CALL Delete(DB)
    DB%MAXPrm = DB%MAXPrm * 5
    CALL New(DB)
  ELSEIF (ErrorCode==eMAXML) THEN
    CALL Delete(MB)
    MB%MAXML = MB%MAXML * 2
    CALL New(MB)
  ELSE
    write(*,*) "Error code = ",ErrorCode
    CALL Halt(' Sorry, fatal error in ONX')
  END IF
  write(*,*) "Leaving MemInit"
END SUBROUTINE MemInit
