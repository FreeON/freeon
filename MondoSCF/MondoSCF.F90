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
  USE GeomOpt
  USE ProcessControl    
  USE InOut
  USE MD
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
!
! initialize Ctrl%Current before parsing!
!
  Ctrl%Current=(/0,1,1/)
  Ctrl%Previous=(/0,1,1/)
  CALL SetGlobalCtrlIndecies(Ctrl)
!
! Parse input
  CALL ParseInp(Ctrl)  
  CALL OpenHDF(InfFile)
#ifdef MMech
   Call InitMMech()
#endif
!
#ifdef MMech
  IF(HasQM())THEN
#endif
     CALL SetSCF(Ctrl)
#ifdef MMech
  ENDIF
#endif
! Decide about forces
!
  SELECT CASE(Ctrl%Grad)
  CASE(GRAD_ONE_FORCE) 
     CALL CALC_GRAD_ONE_FORCE(Ctrl) 
  CASE(GRAD_MD)
     CALL CALC_MD(Ctrl)
  CASE(GRAD_QNEW_OPT)
     CALL CALC_GRAD_QNEW_OPT(Ctrl)
  CASE(GRAD_STPDESC_OPT,GRAD_DIAGHESS_OPT)
     CALL OptimizeNew(Ctrl)
  CASE(GRAD_QNEW_ONE_OPT)
     CALL CALC_GRAD_QNEW_ONE_OPT(Ctrl)
  CASE(GRAD_TS_SEARCH)
     CALL CALC_TS_SEARCH(Ctrl)
  CASE(GRAD_NO_GRAD) !!!!! Energy only calculation
     CALL CALC_GRAD_NO_GRAD(Ctrl)
  END SELECT
!
#ifdef MMech
  IF(HasQM().AND.Ctrl%PopAnalysis(1:8)==OPT_MULLIKEN) THEN
#else
  IF(Ctrl%PopAnalysis(1:8)==OPT_MULLIKEN) THEN
#endif
    CALL Mulliken_Analysis(Ctrl)
  ENDIF
!
  CALL CloseHDF()
!
#if defined(PARALLEL) && defined(MPI2)
  ENDIF
  CALL FiniMPI()
#endif
!--------------------------------------------------------
  CALL TimeStamp('Successful MondoSCF run',.FALSE.)   
!--------------------------------------------------------
  END PROGRAM MondoSCF
