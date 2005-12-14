MODULE DynamicsKeys
  IMPLICIT NONE
! DM Projection Method
  CHARACTER(LEN=*),  PARAMETER :: MD_PM_OPTION       ='MDProjection' 
! DM Verlet
  CHARACTER(LEN=*),  PARAMETER :: MD_DMVerlet0       ='DMVerlet0'
  CHARACTER(LEN=*),  PARAMETER :: MD_DMVerlet1       ='DMVerlet1' 
! FM Verlet  
  CHARACTER(LEN=*),  PARAMETER :: MD_FMVerlet0       ='FMVerlet0'
  CHARACTER(LEN=*),  PARAMETER :: MD_FMVerlet1       ='FMVerlet1' 
! Stencil Denisty Matrix Projectors
  CHARACTER(LEN=*),  PARAMETER :: MD_DMP0            ='DMProj0'
  CHARACTER(LEN=*),  PARAMETER :: MD_DMP1            ='DMProj1'
  CHARACTER(LEN=*),  PARAMETER :: MD_DMP2            ='DMProj2'
  CHARACTER(LEN=*),  PARAMETER :: MD_DMP3            ='DMProj3'
  CHARACTER(LEN=*),  PARAMETER :: MD_DMP4            ='DMProj4'
! Diagonal Geuss
  CHARACTER(LEN=*),  PARAMETER :: MD_DGeuss          ='DMDGeuss'
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
! Initial Temp
  CHARACTER(LEN=*),  PARAMETER :: MD_INIT_TEMP       ='InitTemp'
! Temperature Rescaling
  CHARACTER(LEN=*),  PARAMETER :: MD_TEMP_SCALE      ='TempScaling'
  CHARACTER(LEN=*),  PARAMETER :: MD_TSCALE_INT      ='TempIntScaling'
! MC Inputs
  CHARACTER(LEN=*),  PARAMETER :: MC_MAX_STEP        ='MaxMCStep'
  CHARACTER(LEN=*),  PARAMETER :: MC_TEMP            ='MCTemperature'
!
END MODULE DynamicsKeys 
