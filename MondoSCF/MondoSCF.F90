 !
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!
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
  USE ProcessControl    
  IMPLICIT NONE
  TYPE(SCFControls)     :: Ctrl
  INTEGER               :: ISet,ICyc,IGeo
  INTEGER, DIMENSION(2) :: Begin
!------------------------------------------------------------
! Intialize 
  CALL Init(PerfMon)
  CALL Init(MemStats)
! Parse input
  CALL ParseInp(Ctrl)
! Set up the SCF
  CALL SetSCF(Ctrl)
!
  IF(Ctrl%Rest)THEN
     Begin=Ctrl%Current%I(2:3)
  ELSE
     Begin=(/1,1/)
  ENDIF
  IF(Begin(2)==1)THEN  
!    Start with first geometry and converge 
!    the SCF for each set
     Ctrl%Previous%I=(/0,Begin(1),1/)
     DO ISet=Begin(1),Ctrl%NSet
        DO ICyc=0,MaxSCFs
           Ctrl%Current%I=(/ICyc,ISet,1/)
           IF(ICyc==0) &
           CALL OneEMats(Ctrl)
           CALL SCFCycle(Ctrl)
           IF(ConvergedQ(Ctrl))EXIT
        ENDDO
!       Summarize SCF stats
        CALL SCFSummry(Ctrl)
     ENDDO
     Begin=(/Ctrl%NSet,2/)
  ENDIF 
! Loop over geometries with the last set
  DO IGeo=Begin(2),Ctrl%NGeom
     DO ICyc=0,MaxSCFs
        Ctrl%Current%I=(/ICyc,Ctrl%NSet,IGeo/)
        IF(ICyc==0) &
        CALL OneEMats(Ctrl)
        CALL SCFCycle(Ctrl)
        IF(ConvergedQ(Ctrl))EXIT
     ENDDO
!    Summarize SCF stats
     CALL SCFSummry(Ctrl)
  ENDDO
END PROGRAM
