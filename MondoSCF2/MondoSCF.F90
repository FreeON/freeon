!==============================================================
! 
!==============================================================
PROGRAM MondoSCF
  USE SCFs
  USE Macros
  USE PunchHDF
  USE Optimizer
  USE ParseInput
  USE ZippyQuote
  USE PrintParsed
  IMPLICIT NONE
  TYPE(Controls) :: C
  !------------------------------------------------------------!
  CALL Init(PerfMon)
  CALL Init(MemStats)
#if defined(PARALLEL) && defined(MPI2)
  CALL InitMPI()
#endif
  CALL ParseTheInput(C)  
  ! Initialize the HDF archival file 
  CALL InitArchive(C%Nams)
  ! Print startup 
  CALL PrintsStartUp(C%Nams)
  ! Much ado about gradients
  SELECT CASE(C%Opts%Grad)
  CASE(GRAD_NO_GRAD)
     CALL SinglePoints(C)
  CASE(GRAD_GO_DOWNHILL)
      CALL Descender(C)
  CASE(GRAD_TS_SEARCH_NEB)
     ! Place holder for whatever
      CALL Descender(C)
   CASE(GRAD_DO_DYNAMICS)
      ! Needs work...
  END SELECT
#if defined(PARALLEL) && defined(MPI2)
  CALL FiniMPI()
#endif
  ! Something surreal to celibrate this run
  CALL ZippySez()
  !--------------------------------------------------------
  CALL TimeStamp('Successful MondoSCF run',.FALSE.)   
  !--------------------------------------------------------
END PROGRAM MondoSCF

