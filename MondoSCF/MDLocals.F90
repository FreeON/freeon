!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This contains type definitions for holding some state 
! information for MD that I didn't want to go into Ctrl 
! or GM.
! Author: Hugh Nymeyer
!
! To do: Add Nose-Hoover information, etc.
!
MODULE MDLocals
!
  USE GlobalScalars
  USE GlobalCharacters
  IMPLICIT NONE
!
!---------------------------------------------------------
!
  TYPE MDState
     CHARACTER(LEN=DEFAULT_CHR_LEN)  :: TITLE
     INTEGER                         :: STEP
     REAL(DOUBLE)                    :: TIME
     REAL(DOUBLE)                    :: E_POT
     REAL(DOUBLE)                    :: E_KIN
     REAL(DOUBLE)                    :: E_TOT
     REAL(DOUBLE)                    :: E_TOT_START
     REAL(DOUBLE)                    :: DE_TOT
     REAL(DOUBLE),DIMENSION(3)       :: VIRIAL
     REAL(DOUBLE),DIMENSION(3)       :: PRESSURE
     REAL(DOUBLE)                    :: TEMPERATURE
     REAL(DOUBLE),DIMENSION(3)       :: CM
     REAL(DOUBLE),DIMENSION(3)       :: VCM
     REAL(DOUBLE)                    :: CM_KIN
     REAL(DOUBLE),DIMENSION(3)       :: OMEGA
     REAL(DOUBLE)                    :: OMEGA_KIN
  END TYPE MDState
!
!---------------------------------------------------------
!
END MODULE MDLocals
