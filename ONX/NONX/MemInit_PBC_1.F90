SUBROUTINE MemInit(DB,IB,SB,Drv,BSc,BSp)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE ONXParameters
  USE ONXMemory
  IMPLICIT NONE
  TYPE(DBuf),INTENT(INOUT) :: DB
  TYPE(IBuf),INTENT(INOUT) :: IB
  TYPE(DSL),INTENT(INOUT)  :: SB
  TYPE(IDrv),INTENT(INOUT) :: Drv
  TYPE(BSet),INTENT(IN)    :: BSc,BSp
  INTEGER                  :: LR

!  write(*,*) "In MemInit: ",ErrorCode

  IF (ErrorCode==eInit) THEN
    DB%MAXDis  = 1000
    DB%MAXPrm  = 3000
    DB%MAXD    = 100
    DB%MAXT    = 10
    DB%MAXK    = 10
    IF (Gradient) THEN
      DB%MAXP    = 10
    ELSE
      DB%MAXP    = 8
    ENDIF
    DB%MAXC    = 11
    DB%NPrim   = MAX(BSc%NPrim,BSp%NPrim)
    DB%MInfo   = 0    
    SB%MAXSL   = 10000
    IB%MAXI    = 1000000
    IB%MaxInts = 1000
    IB%NPrim   = MAX(BSc%NPrim,BSp%NPrim)
    IB%Lval    = -1
    IB%MAXL    = 12
    Drv%LngCC  = 60000
    CALL GammaHeader(IB%Mesh,IB%Switch,IB%Grid)
    CALL VRRLng(Drv%LngVRR,Drv%LngLoc)
    IF(.NOT.AllocQ(DB%Alloc)) &
    CALL New(DB)
    IF(.NOT.AllocQ(IB%Alloc)) &
    CALL New(IB)
    IF(.NOT.AllocQ(SB%Alloc)) &
    CALL New(SB)
    IF(.NOT.AllocQ(Drv%Alloc)) &
    CALL New(Drv)
    CALL CCDriver(Drv%CDrv%I(1),Drv%LngDrv)
    CALL VRRDriver(Drv%VLOC%I(1),Drv%SLOC%I(1,1),Drv%LngVRR,Drv%LngLoc)
    CALL GammaAsymptotics(Lr,IB%GammaA%D)
  ELSEIF (ErrorCode==eMAXI) THEN
    CALL Delete(IB)
    IB%MAXI = IB%MAXI * 5
    CALL GammaHeader(IB%Mesh,IB%Switch,IB%Grid)
    CALL New(IB)
    IB%Lval=-1
    CALL GammaAsymptotics(Lr,IB%GammaA%D)
  ELSEIF (ErrorCode==eMAXD) THEN
    CALL Delete(DB)
    DB%MAXD = DB%MAXD * 2
    CALL New(DB)
  ELSEIF (ErrorCode==eMAXK) THEN
    CALL Delete(DB)
    DB%MAXK = DB%MAXK + 1
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
  ELSEIF (ErrorCode==eMAXSL) THEN
    CALL Delete(SB)
    SB%MAXSL = SB%MAXSL * 2
    CALL New(SB)
  ELSE
    write(*,*) "Error code = ",ErrorCode
    CALL Halt(' Sorry, fatal error in ONX')
  END IF
END SUBROUTINE MemInit
