SUBROUTINE MemInit(DB,IB,SB,Drv,BSc,BSp)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE ONXParameters
  USE ONXMemory
  USE DOrder
  IMPLICIT NONE
  TYPE(DBuf),INTENT(INOUT) :: DB
  TYPE(IBuf),INTENT(INOUT) :: IB
  TYPE(DSL),INTENT(INOUT)  :: SB
  TYPE(IDrv),INTENT(INOUT) :: Drv
  TYPE(BSet),INTENT(IN)    :: BSc,BSp
  INTEGER                  :: LR

!  write(*,*) "In MemInit: ",ErrorCode

  IF (ErrorCode==eInit) THEN
    DB%MAXDis  = 10000
    DB%MAXPrm  = 30000
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
    CALL VRRLng(Drv%LngVRR,Drv%LngLoc)
    IF(.NOT.AllocQ(DB%Alloc)) &
    CALL New(DB)
    IF(.NOT.AllocQ(IB%Alloc)) &
    CALL New(IB)
    IF(.NOT.AllocQ(SB%Alloc)) &
    CALL New(SB)
    IF(.NOT.AllocQ(Drv%Alloc)) &
    CALL New(Drv)
    IF(.NOT.AllocQ(BufN%Alloc)) &
    CALL New(BufN,(/DB%MAXT,DB%NPrim*DB%NPrim/))
    IF(.NOT.AllocQ(BufT%Alloc)) &
    CALL New(BufT,(/DB%MAXD,DB%MAXT,DB%MAXK/))
    IF(.NOT.AllocQ(SchT%Alloc)) &
    CALL New(SchT,(/DB%MAXD,DB%MAXT,DB%MAXK/))
    CALL CCDriver(Drv%CDrv%I(1),Drv%LngDrv)
    CALL VRRDriver(Drv%VLOC%I(1),Drv%SLOC%I(1,1),Drv%LngVRR,Drv%LngLoc)
  ELSEIF (ErrorCode==eMAXI) THEN
    CALL Delete(IB)
    IB%MAXI = IB%MAXI * 5
    CALL New(IB)
    IB%Lval=-1
  ELSEIF (ErrorCode==eMAXD) THEN
    CALL Delete(DB)
    CALL Delete(BufT)
    CALL Delete(SchT)
    DB%MAXD = DB%MAXD * 1.3D0
    CALL New(DB)
    CALL New(BufT,(/DB%MAXD,DB%MAXT,DB%MAXK/))
    CALL New(SchT,(/DB%MAXD,DB%MAXT,DB%MAXK/))
  ELSEIF (ErrorCode==eMAXK) THEN
    CALL Delete(DB)
    CALL Delete(BufT)
    CALL Delete(SchT)
    DB%MAXK = DB%MAXK + 1
    CALL New(DB)
    CALL New(BufT,(/DB%MAXD,DB%MAXT,DB%MAXK/))
    CALL New(SchT,(/DB%MAXD,DB%MAXT,DB%MAXK/))
  ELSEIF (ErrorCode==eMAXT) THEN
    CALL Delete(DB)
    CALL Delete(BufN)
    CALL Delete(BufT)
    CALL Delete(SchT)
    DB%MAXT = DB%MAXT + 1
    CALL New(DB)
    CALL New(BufN,(/DB%MAXT,DB%NPrim*DB%NPrim/))
    CALL New(BufT,(/DB%MAXD,DB%MAXT,DB%MAXK/))
    CALL New(SchT,(/DB%MAXD,DB%MAXT,DB%MAXK/))
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
