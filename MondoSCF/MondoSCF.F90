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
#endif
! Initialize
  CALL Init(PerfMon)
  CALL Init(MemStats)
! Parse input
  CALL ParseInp(Ctrl)
! Set up the SCF
  CALL SetSCF(Ctrl)
! Set status with option for restart
  IF(Ctrl%Rest)THEN
     Begin=(/0,Ctrl%Current(2),Ctrl%Current(3)/)
  ELSE
     Begin=(/0,1,1/)
  ENDIF    
  Ctrl%Current=Begin
  Ctrl%Previous=Begin
  CALL SetGlobalCtrlIndecies(Ctrl)           
! Decide about forces
  SELECT CASE(Ctrl%Grad)
  CASE(GRAD_ONE_FORCE)
     DO IGeo=Begin(3),Ctrl%NGeom
        Ctrl%Current(3)=IGeo
        DO ISet=Begin(2),Ctrl%NSet
           Ctrl%Current(2)=ISet
           CALL OneSCF(Ctrl)
           Ctrl%Previous=Ctrl%Current
        ENDDO
        CALL Forces(Ctrl)
     ENDDO
  CASE(GRAD_VERLET_MD)
     CALL MondoHalt(USUP_ERROR,' Look for MD in version 1.0b2. ')
  CASE(GRAD_QNEW_OPT)
     DO ISet=Begin(2),Ctrl%NSet
        Ctrl%Current=(/0,ISet,CGeo/)
        CALL GeOp(Ctrl)
     ENDDO
  CASE(GRAD_TS_SEARCH)
     CALL MondoHalt(USUP_ERROR,' Look for transition state optimizer in version 1.0b2.')
  CASE(GRAD_NO_GRAD)
!    Single points.  Loop first over basis sets...
     DO ISet=Begin(2),Ctrl%NSet
        Ctrl%Current(2)=ISet
        CALL OneSCF(Ctrl)
        Ctrl%Previous=Ctrl%Current
     ENDDO
!    Then over geometries at last basis set
     DO IGeo=Begin(3),Ctrl%NGeom
        Ctrl%Current(3)=IGeo
        CALL OneSCF(Ctrl)
        Ctrl%Previous=Ctrl%Current
     ENDDO
  END SELECT
#if defined(PARALLEL) && defined(MPI2)
  CALL FiniMPI()
#endif
END PROGRAM MondoSCF
