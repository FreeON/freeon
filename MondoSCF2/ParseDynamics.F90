! Authors : Hugh Nymeyer and Matt Challacombe
MODULE ParseDynamics
  USE Parse
  USE DynamicsKeys
  USE ControlStructures
  IMPLICIT NONE
CONTAINS
  !========================================================================================
  ! PARSE FOR MOLECULAR DYNAMICS DIRECTIVES AND VALUES
  !========================================================================================
  SUBROUTINE LoadDynamics(N,O,G,D)
    TYPE(FileNames)  :: N
    TYPE(Options)    :: O
    TYPE(Geometries) :: G
    TYPE(Dynamics)   :: D
    !---------------------------------------------------------------------------------------!
    CALL OpenASCII(N%IFile,Inp)
    D%MAX_STEPS=O%NSteps
    D%AtomWrap=.TRUE.
    IF(OptKeyQ(Inp,MOLDYNE,MD_VV))THEN
       D%MDAlgorithm=MD_SERIAL_VERLET
    ELSEIF(OptKeyQ(Inp,MOLDYNE,MD_PTEMPER))THEN
       D%MDAlgorithm=MD_PARALLEL_REP
       IF(.NOT.OptIntQ(Inp,MD_REPLICAS,G%Klones))THEN
          CALL MondoHalt(PRSE_ERROR,'Did not find input number of parallel replicas to run ')
       ENDIF
    ENDIF
    IF(OptKeyQ(Inp,MOLDYNE,MD_CLOBBER))THEN
       D%CLOBBER=.TRUE.
    ELSE
       D%CLOBBER=.FALSE.
    ENDIF
    ! freq. to write coords
    IF(.NOT.OptIntQ(Inp,MD_CRDFREQ,D%CRDfreq))D%CRDfreq=1
    ! freq to write velocity file
    IF(.NOT.OptIntQ(Inp,MD_VELFREQ,D%VELfreq))D%VELfreq=1
    ! freq to write energy file
    IF(.NOT.OptIntQ(Inp,MD_ENEFREQ,D%ENEfreq))D%ENEfreq=1
    ! freq to write restart
    IF(.NOT.OptIntQ(Inp,MD_RESFREQ,D%RESfreq))D%RESfreq=1
    IF(.NOT.OptDblQ(Inp,MD_TIME_STEP,D%DT))D%DT= 0.5D0
    IF(.NOT.OptDblQ(Inp,MD_TEMP0,D%TEMP0))D%TEMP0=Zero
    IF(.NOT.OptDblQ(Inp,MD_TEMP,D%TEMP))D%TEMP=300D0
    IF(.NOT.OptDblQ(Inp,MD_PRES,D%PRES))D%PRES=1.0D0
    IF(.NOT.OptDblQ(Inp,MD_TTAU,D%TTAU))D%TTAU=1.0D1
    IF(.NOT.OptDblQ(Inp,MD_PTAU,D%PTAU))D%PTAU=1.0D1
    IF(OptIntQ(Inp,MD_REM_TRANS,D%TRANSfreq))THEN
       D%REM_TRANS=.TRUE.
    ELSE
       D%TRANSfreq=0
       D%REM_TRANS=.FALSE.
    ENDIF
    IF(OptIntQ(Inp,MD_REM_ROTAT,D%ROTATfreq))THEN
       D%REM_ROTAT=.TRUE.
    ELSE
       D%ROTATfreq=0
       D%REM_ROTAT=.FALSE.
    ENDIF
    IF(.NOT.OptCharQ(Inp,MD_RESTRT_IN,D%RESTRT_IN))D%RESTRT_IN='Restart.in'
    IF(.NOT.OptCharQ(Inp,MD_RESTRT_OUT,D%RESTRT_OUT))D%RESTRT_OUT='Restart.out'
    IF(.NOT.OptCharQ(Inp,MD_CRD_OUT,D%CRD_NAME))D%CRD_NAME='crd'
    IF(.NOT.OptCharQ(Inp,MD_VEL_OUT,D%VEL_NAME))D%VEL_NAME='vel'
    IF(.NOT.OptCharQ(Inp,MD_ENE_OUT,D%ENE_NAME))D%ENE_NAME='ene'
    CLOSE(UNIT=Inp,STATUS='KEEP')
    ! Convert to internal units
    D%DT   = D%DT   * 1.0D-15 * SecondsToInternalTime
    D%TTAU = D%TTAU * 1.0D-15 * SecondsToInternalTime
    D%PTAU = D%PTAU * 1.0D-15 * SecondsToInternalTime
  END SUBROUTINE LoadDynamics
END MODULE ParseDynamics
