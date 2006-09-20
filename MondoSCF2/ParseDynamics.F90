! Authors : Hugh Nymeyer, CJ Tymczak and Matt Challacombe
MODULE ParseDynamics
  USE Parse
  USE InOut
  USE PrettyPrint
  USE DynamicsKeys
  USE ControlStructures
  USE OptionKeys
  IMPLICIT NONE
CONTAINS
!========================================================================================
! PARSE FOR MOLECULAR DYNAMICS DIRECTIVES AND VALUES
!========================================================================================
  SUBROUTINE LoadDynamics(N,O,G,D)
    TYPE(FileNames)    :: N
    TYPE(Options)      :: O
    TYPE(Geometries)   :: G
    TYPE(Dynamics)     :: D
    CHARACTER(LEN=DCL) :: Line
!---------------------------------------------------------------------------------------!
    CALL OpenASCII(N%IFile,Inp)       
!   Initialize
    D%DoingMD        = .TRUE.
    D%Initial_Temp   = .FALSE.
    D%Temp_Scaling   = .FALSE.
    D%Const_Temp     = .FALSE.
    D%Const_Press    = .FALSE.
    D%Parallel_Rep   = .FALSE.
    D%MDGeuss        = MD_DGeuss
!   Parse MD
    IF(O%Grad==GRAD_DO_DYNAMICS) THEN
!      Parse the Density Matrix Projection Algoithm
       IF(    OptKeyQ(Inp,MD_PM_OPTION,MD_DMVerlet0)) THEN
          D%MDGeuss=MD_DMVerlet0
       ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMVerlet1)) THEN
          D%MDGeuss=MD_DMVerlet1
       ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_FMVerlet0)) THEN
          D%MDGeuss=MD_FMVerlet0
       ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_FMVerlet1)) THEN
          D%MDGeuss=MD_FMVerlet1
       ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP0)    ) THEN
          D%MDGeuss=MD_DMP0
       ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP1)    ) THEN
          D%MDGeuss=MD_DMP1
       ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP2)    ) THEN
          D%MDGeuss=MD_DMP2
       ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP3)    ) THEN
          D%MDGeuss=MD_DMP3
       ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP4)    ) THEN
          D%MDGeuss=MD_DMP4
       ELSE
          WRITE(*,*) ' In input file, no MD DM Projection algorithm Defined'
          CALL Warn(' In input file, no MD DM Projection algorithm Defined')
       ENDIF
!      Parse MD Options: First MD Algorithmn
       IF(OptKeyQ(Inp,MD_AL_OPTION,MD_AL_VERLET))THEN
          D%MDAlgorithm = VERLET_MD_AL
       ELSEIF(OptKeyQ(Inp,MD_AL_OPTION,MD_AL_GEAR)) THEN
          D%MDAlgorithm = GEAR_MD_AL
       ELSE
          D%MDAlgorithm = VERLET_MD_AL
       ENDIF
!      MD Time Step
       IF(.NOT. OptDblQ(Inp,MD_TIME_STEP,D%DTime)) THEN
          CALL MondoHalt(PRSE_ERROR,MD_TIME_STEP//' not found in input.')
       ENDIF
!      MD MaxSteps
       IF(.NOT. OptIntQ(Inp,MD_MAX_STEP,D%MDMaxSteps)) THEN
          CALL MondoHalt(PRSE_ERROR,MD_MAX_STEP//' not found in input.')
       ENDIF
!      Parse for Number of SCF cycles
       IF(.NOT. OptIntQ(Inp,MD_NUM_SCF,D%MDNumSCF)) THEN
          CALL Warn('Number of SCF cycles for MD is not set: DEFAULT is eight')
          D%MDNumSCF=8
       ENDIF
!      Parse for Initial Temp, if any
       IF(OptDblQ(Inp,MD_INIT_TEMP,D%TempInit))THEN
          D%Initial_Temp   = .TRUE.
       ELSE
          D%Initial_Temp   = .FALSE.
          D%TempInit       = Zero
       ENDIF
!      Parse for Temperature Scaling
       IF(OptDblQ(Inp,MD_TEMP_SCALE,D%TargetTemp)) THEN
          D%Temp_Scaling = .TRUE.
          IF(.NOT. OptIntQ(Inp,MD_TSCALE_INT,D%RescaleInt)) THEN
             D%RescaleInt   = 100
          ENDIF
       ELSE
          D%Temp_Scaling = .FALSE.
          D%TargetTemp   = Zero
          D%RescaleInt   = 1          
       ENDIF
!   Parse Hybrid MC
    ELSEIF(O%Grad==GRAD_DO_HYBRIDMC) THEN
!      Parse the Density Matrix Projection Algoithm
       IF(    OptKeyQ(Inp,MD_PM_OPTION,MD_DMVerlet0)) THEN
          D%MDGeuss=MD_DMVerlet0
       ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMVerlet1)) THEN
          D%MDGeuss=MD_DMVerlet1
       ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_FMVerlet0)) THEN
          D%MDGeuss=MD_FMVerlet0
       ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_FMVerlet1)) THEN
          D%MDGeuss=MD_FMVerlet1
       ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP0)     ) THEN
          D%MDGeuss=MD_DMP0
       ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP1)     ) THEN
          D%MDGeuss=MD_DMP1
       ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP2)     ) THEN
          D%MDGeuss=MD_DMP2
       ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP3)     ) THEN
          D%MDGeuss=MD_DMP3
       ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP4)     ) THEN
          D%MDGeuss=MD_DMP4
       ELSE
          WRITE(*,*) ' In input file, no MD DM Projection algorithm Defined'
          CALL Warn('In input file, no MD DM Projection algorithm Defined')
       ENDIF
!      MD Algorithmn
       IF(OptKeyQ(Inp,MD_AL_OPTION,MD_AL_VERLET))THEN
          D%MDAlgorithm = VERLET_MD_AL
       ELSEIF(OptKeyQ(Inp,MD_AL_OPTION,MD_AL_GEAR)) THEN
          D%MDAlgorithm = GEAR_MD_AL
       ELSE
          D%MDAlgorithm = VERLET_MD_AL
       ENDIF
!      MD Time Step
       IF(.NOT. OptDblQ(Inp,MD_TIME_STEP,D%DTime)) THEN
          CALL MondoHalt(PRSE_ERROR,MD_TIME_STEP//' not found in input.')
       ENDIF
!      MD MaxSteps
       IF(.NOT. OptIntQ(Inp,MD_MAX_STEP,D%MDMaxSteps)) THEN
          CALL MondoHalt(PRSE_ERROR,MD_MAX_STEP//' not found in input.')
       ENDIF
!      Parse Hybrid MC Options: First  MC MaxSteps
       IF(.NOT. OptIntQ(Inp,MC_MAX_STEP,D%MCMaxSteps)) THEN
          CALL MondoHalt(PRSE_ERROR,MC_MAX_STEP//' not found in input.')
       ENDIF
!      MC Temp
       IF(.NOT. OptDblQ(Inp,MC_TEMP,D%MCTemp))THEN
          CALL MondoHalt(PRSE_ERROR,MC_TEMP//' not found in input.')
       ENDIF
    ENDIF
!
  END SUBROUTINE LoadDynamics
!---------------------------------------------------------------------------------------!
END MODULE ParseDynamics
