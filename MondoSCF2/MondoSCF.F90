!==============================================================
PROGRAM MondoSCF
  USE SCFs
  USE Macros
  USE Response
  USE PunchHDF
  USE Optimizer
  USE ParseInput
  USE ZippyQuote
  USE PrintParsed
  USE MDynamics
  USE MonteCarlo
  IMPLICIT NONE
  TYPE(Controls) :: C
  !------------------------------------------------------------!
  CALL Init(PerfMon)
  CALL Init(MemStats)
#if defined(PARALLEL) && defined(MPI2)
  CALL InitMPI()
  InParallel=.FALSE.
  IF(MyID==0)THEN
#endif
  CALL ParseTheInput(C)  
! Initialize controls
  CALL InitGlobal(C)
! Print startup 
  CALL PrintsStartUp(C%Nams)
! Much ado about gradients
  SELECT CASE(C%Opts%Grad)
  CASE(GRAD_NO_GRAD)
     CALL SinglePoints(C)
     CALL CPSCF(C)
  CASE(GRAD_GO_DOWNHILL)
     CALL Descender(C)
  CASE(GRAD_TS_SEARCH_NEB)  
!    Place holder for whatever
     CALL Descender(C)
  CASE(GRAD_DO_DYNAMICS)
     CALL MD(C)
  CASE(GRAD_DO_HYBRIDMC)
     CALL HybridMC(C)
  CASE(GRAD_DO_SCAN)
     CALL Halt('SCAN Not Implimented')
!     CALL ScanGeom(C)
  CASE(GRAD_ONE_FORCE)
     CALL SinglePoints(C)
     CALL Force(C%Sets%NBSets,1,C%Nams,C%Opts,C%Stat,C%Geos,C%Sets,C%MPIs)
  CASE(GRAD_DO_NHESSIAN)
     CALL NHessian(C)
  END SELECT
  ! Something surreal to celebrate this run
  CALL ZippySez(C)
  !--------------------------------------------------------
  CALL TimeStamp('Successful MondoSCF run',.FALSE.)   
  !--------------------------------------------------------
#if defined(PARALLEL) && defined(MPI2)
  ENDIF
  CALL AlignNodes("MONDOSCF")
  CALL FiniMPI()
#endif
END PROGRAM MondoSCF

