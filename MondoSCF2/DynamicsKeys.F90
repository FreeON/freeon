MODULE DynamicsKeys
  IMPLICIT NONE
! MD Algorithm
  CHARACTER(LEN=*),  PARAMETER :: MD_AL_OPTION       ='MDMethod' 
! Velocity Velet
  CHARACTER(LEN=*),  PARAMETER :: MD_AL_VERLET       ='Verlet'
  INTEGER,           PARAMETER :: VERLET_MD_AL       =34235423
! Gear PC
  CHARACTER(LEN=*),  PARAMETER :: MD_AL_GEAR         ='Gear'
  INTEGER,           PARAMETER :: GEAR_MD_AL         =34235424
! MD Inputs
  CHARACTER(LEN=*),  PARAMETER :: MD_TIME_STEP       ='DeltaTime'
  CHARACTER(LEN=*),  PARAMETER :: MD_MAX_STEP        ='MaxMDStep' 
!
END MODULE DynamicsKeys
