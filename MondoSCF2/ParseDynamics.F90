! Authors : Hugh Nymeyer and Matt Challacombe
MODULE ParseDynamics
  USE Parse
  USE InOut
  USE PrettyPrint
  USE DynamicsKeys
  USE ControlStructures
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
    D%Velcty_Scaling = .FALSE.
    D%Const_Temp     = .FALSE.
    D%Const_Press    = .FALSE.
    D%Parallel_Rep   = .FALSE.
!   Parse MD Options: First MD Algorithmn
    IF(OptKeyQ(Inp,MD_AL_OPTION,MD_AL_VERLET))THEN
       D%MDAlgorithm = VERLET_MD_AL
    ELSEIF(OptKeyQ(Inp,MD_AL_OPTION,MD_AL_GEAR)) THEN
       D%MDAlgorithm = GEAR_MD_AL
    ELSE
       D%MDAlgorithm = VERLET_MD_AL
    ENDIF
!   MD Time Step
    IF(.NOT. OptDblQ(Inp,MD_TIME_STEP,D%DTime)) THEN
       CALL MondoHalt(PRSE_ERROR,MD_TIME_STEP//' not found in input.')
    ENDIF
!   MD MaxSteps
    IF(.NOT. OptIntQ(Inp,MD_MAX_STEP,D%MDMaxSteps)) THEN
       CALL MondoHalt(PRSE_ERROR,MD_MAX_STEP//' not found in input.')
    ENDIF
!
    CLOSE(UNIT=Inp,STATUS='KEEP')
!
  END SUBROUTINE LoadDynamics
!---------------------------------------------------------------------------------------!
END MODULE ParseDynamics
