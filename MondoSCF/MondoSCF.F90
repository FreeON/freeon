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
#ifdef MMech
  USE Mechanics
#endif
!
  IMPLICIT NONE
  TYPE(SCFControls)     :: Ctrl
  INTEGER               :: ICyc,ISet,IGeo,Bak,IChk,K
  INTEGER, DIMENSION(3) :: Begin
  INTEGER, DIMENSION(3) :: PrevState,CrntState
  CHARACTER(LEN=DCL)    :: RmAllPrv
  LOGICAL               :: DoForce
!------------------------------------------------------------

#if defined(PARALLEL) && defined(MPI2)
WRITE(*,*)' INITING '
  CALL InitMPI()
WRITE(*,*)' INITed '
  InParallel=.FALSE.
  IF(MyID==0)THEN
#endif
! Initialize
  CALL Init(PerfMon)
  CALL Init(MemStats)
! Parse input
  CALL ParseInp(Ctrl)
#ifdef MMech
  MechFlag=Ctrl%Mechanics
#endif
!
! Set up the SCF or MM
  Ctrl%Current=(/0,1,1/)
  Ctrl%Previous=(/0,1,1/)
#ifdef MMech
  IF(HasQM())THEN
#endif
     CALL SetSCF(Ctrl)
#ifdef MMech
  ENDIF
#endif
  CALL SetGlobalCtrlIndecies(Ctrl)           
! Decide about forces
  SELECT CASE(Ctrl%Grad)
  CASE(GRAD_ONE_FORCE)
!    Loop first over basis sets 
     DO ISet=1,Ctrl%NSet
        Ctrl%Current=(/0,ISet,1/)
        CALL OneSCF(Ctrl)
     ENDDO
!    Calculate the Force (last basis set)
     CALL Forces(Ctrl)
!    Go over additional geometries at last basis set
     IF(Ctrl%NGeom>1)THEN
        DO IGeo=2,Ctrl%NGeom
           Ctrl%Current=(/0,Ctrl%NSet,IGeo/)
           CALL OneSCF(Ctrl)
           CALL Forces(Ctrl)
        ENDDO
     ENDIF
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
!
  CASE(GRAD_NO_GRAD) !!!!! Energy only calculation
!
#ifdef MMech
   If(HasQM()) Then
#endif
     DO ISet=1,Ctrl%NSet
        Ctrl%Current=(/0,ISet,1/)
#ifdef MMech
! Ooopss... KN forgot to add an ifdef here.  Also, should get hidden in OneSCF probably
        IF(ISet==1.AND.HasMM()) CALL MM_ENERG(Ctrl)
#endif
        CALL OneSCF(Ctrl)
     ENDDO
     IF(Ctrl%NGeom>1)THEN
!       Go over additional geometries at last basis set
        DO IGeo=2,Ctrl%NGeom
           Ctrl%Current=(/0,Ctrl%NSet,IGeo/)
#ifdef MMech
! Ooopss... KN forgot to add an ifdef here.  Also, should get hidden in OneSCF probably
           IF(HasMM()) CALL MM_ENERG(Ctrl) !!! temporary; overwrites energy terms calcd at prev geoms
#endif
           CALL OneSCF(Ctrl)
        ENDDO
     ENDIF
#ifdef MMech
   EndIf
#endif
!
!    Loop first over basis sets 
#ifdef MMech
   If(MMOnly()) Then
     Ctrl%Current(2)=1
     CALL MM_ENERG(Ctrl)
   ENDIF
#endif
  END SELECT
#if defined(PARALLEL) && defined(MPI2)
  ENDIF
  CALL FiniMPI()
#endif
  CALL TimeStamp('Succesfull MondoSCF run',.FALSE.)   

END PROGRAM MondoSCF
