!    MAIN PROGRAM
!    Authors:  Matt Challacombe and C.J. Tymczak
!------------------------------------------------------------------------------
PROGRAM MondoSCF
  USE DerivedTypes
  USE GlobalScalars
  USE InOut
  USE Macros
  USE Parse  
  USE SCFLocals
  USE ParseInput
  USE PrettyPrint
  USE SetSCFs
  USE ChkSCFs
  USE DrvSCFs
  USE DrvFrcs
  USE ProcessControl    
  USE InOut
  IMPLICIT NONE
  TYPE(SCFControls)     :: Ctrl
  INTEGER               :: ICyc,ISet,IGeo,Bak,IChk,K
  INTEGER, DIMENSION(3) :: Begin
  INTEGER, DIMENSION(3) :: PrevState,CrntState
  CHARACTER(LEN=DCL)    :: RmAllPrv
  LOGICAL               :: DoForce
!------------------------------------------------------------
! Intialize 
  CALL Init(PerfMon)
  CALL Init(MemStats)
! Parse input
  CALL ParseInp(Ctrl)
! Set up the SCF
  CALL SetSCF(Ctrl)
! Set up the Restart
  IF(Ctrl%Rest)THEN
     Begin=(/0,Ctrl%Current(2),Ctrl%Current(3)/)
  ELSE
     Begin=(/0,1,1/)
  ENDIF     
  Ctrl%Previous=Begin
!
! Start with first geometry and converge 
! the SCF for each set 
!
  DO IGeo=Begin(3),Ctrl%NGeom
     DO ISet=Begin(2),Ctrl%NSet
        DO ICyc=Begin(1),MaxSCFs
           Ctrl%Current=(/ICyc,ISet,IGeo/)
           IF(ICyc==0) &
           CALL OneEMats(Ctrl)
           CALL SCFCycle(Ctrl,Begin)
           IF(ConvergedQ(Ctrl))GOTO 999
        ENDDO
        CALL MondoHalt(DRIV_ERROR,'Failed to converge SCF in ' &
                                  //TRIM(IntToChar(MaxSCFs))   &
                                  //' SCF iterations. ')
999     CONTINUE
        Begin(1) = 0
!       Summarize SCF stats
        CALL SCFSummry(Ctrl) 
     ENDDO
     Begin(2) = 1
!    Do Action
     SELECT CASE (Ctrl%ForceAction)
     CASE('MolecularDynamics')
        CALL Forces(Ctrl)
        CALL NextMDGeometry(Ctrl)
     CASE('GeometryOptimization')   
!        CALL Forces(Ctrl)        
!        CALL NextGOGeometry(Ctrl,IChk)
     CASE('Force')
        CALL Forces(Ctrl)
     END SELECT
  ENDDO
!
END PROGRAM MondoSCF
!
!!$  Ctrl%NGeom=50
!!$  IGeo=Begin(2)
!!$  Bak=1
!!$  IChk=1
!!$  CrntState=Ctrl%Current
!!$  PrevState=Ctrl%Current
!!$!  WRITE(*,*)'1 Previous',TRIM(StatsToChar(Ctrl%Previous))
!!$  CALL OpenHDF(InfFile)
!!$  CALL Get(GM,IntToChar(IGeo-1))
!!$  CALL Put(GM,IntToChar(IGeo))
!!$  CALL Put(One,'StepSize')
!!$  CALL CloseHDF()
!!$! Loop over geometries with the last set
!!$  DO K=1,50
!!$     Ctrl%Previous=PrevState
!!$     Ctrl%Current=CrntState        
!!$     CALL NewStep(Ctrl,IChk)
!!$!      STOP
!!$!------------------------------------------------
!!$!     CALL OpenHDF(InfFile)
!!$!     CALL Get(GM,IntToChar(IGeo-1))
!!$!     CALL PPrint(GM,Unit_O=6,PrintGeom_O='XYZ')
!!$!     CALL Get(GM,IntToChar(IGeo))
!!$!     CALL PPrint(GM,Unit_O=6,PrintGeom_O='XYZ')
!!$!     CALL CloseHDF()
!!$!------------------------------------------------
!!$     DO ICyc=0,20
!!$        Ctrl%Current=(/ICyc,Ctrl%NSet,IGeo/)
!!$        IF(ICyc==0) &
!!$        CALL OneEMats(Ctrl)
!!$        CALL SCFCycle(Ctrl)
!!$        IF(ConvergedQ(Ctrl))GOTO 9999
!!$     ENDDO
!!$     IF(BaK==4) &
!!$     CALL MondoHalt(DRIV_ERROR,'Failed to converge SCF in 10 SCF iterations. ')
!!$9999 CONTINUE
!!$     Ctrl%Previous=PrevState
!!$     IChk=ChkStep(Ctrl,GM,Bak)
!!$     IF(IChk==1)THEN
!!$!       Summarize SCF stats
!!$        CALL SCFSummry(Ctrl)
!!$!       Increment geometry counter
!!$        CALL OpenHDF(InfFile)
!!$        CALL Put(One,'StepSize')
!!$        CALL CloseHDF()
!!$        RmAllPrv=TRIM(Ctrl%Name)//'_Geom#'//TRIM(IntToChar(IGeo-1))//'*' 
!!$        CALL SYSTEM('/bin/rm '//TRIM(RmAllPrv))
!!$!        CALL Invoke('ls',(/RmAllPrv/),AbsPath_O=.TRUE.)
!!$        IGeo=IGeo+1
!!$        PrevState=Ctrl%Current 
!!$        CrntState=Ctrl%Current 
!!$     ELSEIF(IChk==-1)THEN
!!$        EXIT
!!$     ENDIF
!!$  ENDDO
!!$  CALL SCFSummry(Ctrl)
!!$END PROGRAM
!     CALL PPrint(GM,GeoFile,Geo,'XYZ')
!    Summarize SCF stats
