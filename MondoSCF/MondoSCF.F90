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

#if defined(PARALLEL) && defined(MPI2)
  CALL InitMPI()
  InParallel=.FALSE.
  IF(MyID==0)THEN
#endif
! Initialize
  CALL Init(PerfMon)
  CALL Init(MemStats)
! Parse input
  CALL ParseInp(Ctrl)
! Set up the SCF
  CALL SetSCF(Ctrl)
  Ctrl%Current=(/0,1,1/)
  Ctrl%Previous=(/0,1,1/)
  CALL SetGlobalCtrlIndecies(Ctrl)           
! Decide about forces
  SELECT CASE(Ctrl%Grad)
  CASE(GRAD_ONE_FORCE)
     DO IGeo=1,Ctrl%NGeom
        Ctrl%Current(3)=IGeo
        DO ISet=1,Ctrl%NSet
           Ctrl%Current(2)=ISet
           CALL OneSCF(Ctrl)
        ENDDO
        CALL Forces(Ctrl)
        CALL NForce(Ctrl)
     ENDDO
  CASE(GRAD_MD)
     CALL MondoHalt(USUP_ERROR,' Look for MD in version 1.0b2. ')
  CASE(GRAD_QNEW_OPT)
     DO ISet=1,Ctrl%NSet
!       Optimize geometry for each basis set
        Ctrl%Current=(/0,ISet,CGeo/)
        CALL GeOp(Ctrl)
     ENDDO
  CASE(GRAD_QNEW_ONE_OPT)
!    Loop first over basis sets.
     DO ISet=1,Ctrl%NSet-1
        Ctrl%Current(2)=ISet
        CALL OneSCF(Ctrl)
     ENDDO
!    Optimize geometry only in last basis set
     Ctrl%Current=(/0,Ctrl%NSet,CGeo/)
     CALL GeOp(Ctrl)
  CASE(GRAD_TS_SEARCH)
     CALL MondoHalt(USUP_ERROR,' Look for transition state optimizer in version 1.0b2.')
  CASE(GRAD_NO_GRAD)
!    Loop first over basis sets 
     DO ISet=1,Ctrl%NSet
        Ctrl%Current(2)=ISet
        CALL OneSCF(Ctrl)
     ENDDO
     IF(Ctrl%NGeom>1)THEN
!       Go over additional geometries at last basis set
        DO IGeo=2,Ctrl%NGeom
           Ctrl%Current=(/0,Ctrl%NSet,IGeo/)
           CALL OneSCF(Ctrl)
        ENDDO
     ENDIF
  END SELECT
#if defined(PARALLEL) && defined(MPI2)
  ENDIF
  CALL FiniMPI()
#endif
  CALL TimeStamp('Succesfull MondoSCF run',.FALSE.)   
END PROGRAM MondoSCF
