MODULE ONXMemInit
!H=================================================================================
!H MODULE MemoryInit
!H This MODULE contains:
!H  o SUB MemInit
!H  o SUB 
!H
!H Comments:
!H  - Should delete(BufT,BufP,BufN)
!H---------------------------------------------------------------------------------
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE ONXParameters
  USE ONXMemory
  USE ONXDOrder  !just need that for the buffer BufN, BufT, SchT.
  !
  IMPLICIT NONE
  PRIVATE
  !
!--------------------------------------------------------------------------------- 
! PUBLIC DECLARATIONS
!--------------------------------------------------------------------------------- 
  PUBLIC :: MemInit
  !
CONTAINS
  !
  SUBROUTINE MemInit(DB,IB,SB,Drv,BSc,BSp)
!H--------------------------------------------------------------------------------- 
!H SUBROUTINE MemInit(DB,IB,SB,Drv,BSc,BSp)
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(DBuf),INTENT(INOUT) :: DB
    TYPE(IBuf),INTENT(INOUT) :: IB
    TYPE(DSL ),INTENT(INOUT) :: SB
    TYPE(IDrv),INTENT(INOUT) :: Drv
    TYPE(BSet),INTENT(IN   ) :: BSc,BSp
    INTEGER                  :: LR
!
    IF (ErrorCode==eInit) THEN
       DB%MAXDis  = 10000
       DB%MAXPrm  = 30000
       DB%MAXD    = 100
       DB%MAXT    = 10
       DB%MAXK    = 100
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
       IF(.NOT.AllocQ(DB%Alloc))  CALL New(DB)
       IF(.NOT.AllocQ(IB%Alloc))  CALL New(IB)
       IF(.NOT.AllocQ(SB%Alloc))  CALL New(SB)
       IF(.NOT.AllocQ(Drv%Alloc)) CALL New(Drv)

       IF(.NOT.AllocQ(BufN%Alloc)) CALL New(BufN,(/DB%MAXT,DB%NPrim*DB%NPrim/))    !per
       IF(.NOT.AllocQ(BufT%Alloc)) CALL New(BufT,(/DB%MAXD,DB%MAXT,DB%MAXK/))      !per
       IF(.NOT.AllocQ(SchT%Alloc)) CALL New(SchT,(/DB%MAXD,DB%MAXT,DB%MAXK/))      !per

       CALL CCDriver(Drv%CDrv%I(1),Drv%LngDrv)
       CALL VRRDriver(Drv%VLOC%I(1),Drv%SLOC%I(1,1),Drv%LngVRR,Drv%LngLoc)
    ELSEIF (ErrorCode==eMAXI) THEN
       CALL Delete(IB)
       IB%MAXI = IB%MAXI * 5
       CALL New(IB)
       IB%Lval=-1
    ELSEIF (ErrorCode==eMAXD) THEN
       CALL Delete(DB)

       CALL Delete(BufT) !per
       CALL Delete(SchT) !per

!old       DB%MAXD = DB%MAXD * 2
       DB%MAXD = DB%MAXD * 2 !1.3D0 !per

       CALL New(DB)

       CALL New(BufT,(/DB%MAXD,DB%MAXT,DB%MAXK/))  !per
       CALL New(SchT,(/DB%MAXD,DB%MAXT,DB%MAXK/))  !per

    ELSEIF (ErrorCode==eMAXK) THEN
       CALL Delete(DB)

       CALL Delete(BufT) !per
       CALL Delete(SchT) !per

       DB%MAXK = DB%MAXK + 1
       CALL New(DB)

       CALL New(BufT,(/DB%MAXD,DB%MAXT,DB%MAXK/))   !per
       CALL New(SchT,(/DB%MAXD,DB%MAXT,DB%MAXK/))   !per

    ELSEIF (ErrorCode==eMAXT) THEN
       CALL Delete(DB)

       CALL Delete(BufN)   !per
       CALL Delete(BufT)   !per
       CALL Delete(SchT)   !per

       DB%MAXT = DB%MAXT + 1
       CALL New(DB)

       CALL New(BufN,(/DB%MAXT,DB%NPrim*DB%NPrim/))  !per
       CALL New(BufT,(/DB%MAXD,DB%MAXT,DB%MAXK/))    !per
       CALL New(SchT,(/DB%MAXD,DB%MAXT,DB%MAXK/))    !per

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
       WRITE(*,*) "Error code = ",ErrorCode
       CALL Halt(' Sorry, fatal error in ONX')
    END IF
!
  END SUBROUTINE MemInit
!
END MODULE ONXMemInit
