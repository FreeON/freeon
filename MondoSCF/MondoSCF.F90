!------------------------------------------------------------------------------
!--  This code is part of the MondoSCF suite of programs for linear scaling 
!    electronic structure theory and ab initio molecular dynamics.
!
!--  Copyright (c) 2001, the Regents of the University of California.  
!    This SOFTWARE has been authored by an employee or employees of the 
!    University of California, operator of the Los Alamos National Laboratory 
!    under Contract No. W-7405-ENG-36 with the U.S. Department of Energy.  
!    The U.S. Government has rights to use, reproduce, and distribute this 
!    SOFTWARE.  The public may copy, distribute, prepare derivative works 
!    and publicly display this SOFTWARE without charge, provided that this 
!    Notice and any statement of authorship are reproduced on all copies.  
!    Neither the Government nor the University makes any warranty, express 
!    or implied, or assumes any liability or responsibility for the use of 
!    this SOFTWARE.  If SOFTWARE is modified to produce derivative works, 
!    such modified SOFTWARE should be clearly marked, so as not to confuse 
!    it with the version available from LANL.  The return of derivative works
!    to the primary author for integration and general release is encouraged. 
!    The first publication realized with the use of MondoSCF shall be
!    considered a joint work.  Publication of the results will appear
!    under the joint authorship of the researchers nominated by their
!    respective institutions. In future publications of work performed
!    with MondoSCF, the use of the software shall be properly acknowledged,
!    e.g. in the form "These calculations have been performed using MondoSCF, 
!    a suite of programs for linear scaling electronic structure theory and
!    ab initio molecular dynamics", and given appropriate citation.  
!------------------------------------------------------------------------------
!    Author:  Matt Challacombe
!    MAIN DRIVER
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
!       Summarize SCF stats
        CALL SCFSummry(Ctrl) 
     ENDDO
!    Do Action
     SELECT CASE (Ctrl%ForceAction)
     CASE('Molecular-Dynamics')
        CALL Forces(Ctrl)
!        CALL NextMDGeometry(Ctrl)
     CASE('Geometry-Optimization') 
!        IChk=CheckStep(Ctrl)      
        CALL Forces(Ctrl)        
!        CALL NextGOGeometry(Ctrl,IChk)
     CASE('Forces')
        WRITE(*,*) 'Doing Force'
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
