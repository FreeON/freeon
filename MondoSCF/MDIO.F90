!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module handles i/o for special MD files: CRD,   !
! VEL, ENE, and RESTART files
! Author: Hugh Nymeyer                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
MODULE MDIO
!
  USE GlobalCharacters        ! DEFAULT_CHR_LEN
  USE GlobalScalars           ! Unit numbers
  USE SCFLocals               ! SCFControls, MDControls
  USE MDLocals                ! MDState
!
  IMPLICIT NONE
  CONTAINS
!
!--------------------------------------------------------
!
  SUBROUTINE OpenCRD(FILENAME,TITLE,CLOBBER)
!
    CHARACTER(LEN=DEFAULT_CHR_LEN)       :: FILENAME
    CHARACTER(LEN=DEFAULT_CHR_LEN)       :: TITLE
    LOGICAL                              :: CLOBBER
    INTEGER                              :: OPEN_ERROR
!
    IF (CLOBBER) THEN
       OPEN(UNIT   = CRD_UNIT,    &
            FILE   = FILENAME,    &
            STATUS = "REPLACE",   &
            FORM   = "FORMATTED", &
            IOSTAT = OPEN_ERROR)
    ELSE
       OPEN(UNIT   = CRD_UNIT,    &
            FILE   = FILENAME,    &
            STATUS = "NEW",       &
            FORM   = "FORMATTED", &
            IOSTAT = OPEN_ERROR)
    ENDIF
!
    IF (OPEN_ERROR .NE. 0) THEN
       CALL HALT('! OpenCRD failed')
    END IF
!
    WRITE(UNIT=CRD_UNIT,FMT='(A80)')TITLE
!
  END SUBROUTINE OpenCRD
!
!--------------------------------------------------------
!
  SUBROUTINE OpenVEL(FILENAME,TITLE,CLOBBER)
!
    CHARACTER(LEN=DEFAULT_CHR_LEN)       :: FILENAME
    CHARACTER(LEN=DEFAULT_CHR_LEN)       :: TITLE
    LOGICAL                              :: CLOBBER
    INTEGER                              :: OPEN_ERROR
!
    IF (CLOBBER) THEN
       OPEN(UNIT   = VEL_UNIT,    &
            FILE   = FILENAME,    &
            STATUS = "REPLACE",   &
            FORM   = "FORMATTED", &
            IOSTAT = OPEN_ERROR)
    ELSE
       OPEN(UNIT   = VEL_UNIT,    &
            FILE   = FILENAME,    &
            STATUS = "NEW",       &
            FORM   = "FORMATTED", &
            IOSTAT = OPEN_ERROR)
    ENDIF
!
    IF (OPEN_ERROR .NE. 0) THEN
       CALL HALT('! OpenVEL failed')
    END IF
!
    WRITE(UNIT=VEL_UNIT,FMT='(A80)')TITLE
!    
!
  END SUBROUTINE OpenVEL

!
!--------------------------------------------------------
!
  SUBROUTINE OpenENE(FILENAME,TITLE,CLOBBER)
!
    CHARACTER(LEN=DEFAULT_CHR_LEN)       :: FILENAME
    CHARACTER(LEN=DEFAULT_CHR_LEN)       :: TITLE
    LOGICAL                              :: CLOBBER
    INTEGER                              :: OPEN_ERROR
!
    IF (CLOBBER) THEN
       OPEN(UNIT   = ENE_UNIT,    &
            FILE   = FILENAME,    &
            STATUS = "REPLACE",   &
            FORM   = "FORMATTED", &
            IOSTAT = OPEN_ERROR)
    ELSE
       OPEN(UNIT   = ENE_UNIT,    &
            FILE   = FILENAME,    &
            STATUS = "NEW",       &
            FORM   = "FORMATTED", &
            IOSTAT = OPEN_ERROR)
    ENDIF
!
    IF (OPEN_ERROR .NE. 0) THEN
       CALL HALT('! OpenENE failed')
    END IF
!
    WRITE(UNIT=ENE_UNIT,FMT='(A80)')TITLE
    WRITE(UNIT=ENE_UNIT,FMT='(A80)') &
	"--------------------------------------------------------------------------------"
    WRITE(UNIT=ENE_UNIT,FMT='(" 0 ",6(A14,X))') &
         "STEP_NUMBER","TIME","E_POT","E_KIN","E_TOT","DE_TOT"
    WRITE(UNIT=ENE_UNIT,FMT='(" 1 ",4(A14,X))') &
         "BOX_X","BOX_Y","BOX_Z","VOLUME"
    WRITE(UNIT=ENE_UNIT,FMT='(" 2 ",5(A14,X))') &
         "VIRIAL_X","VIRIAL_Y","VIRIAL_Z","SET_PRESSURE","PRESSURE"
    WRITE(UNIT=ENE_UNIT,FMT='(" 3 ",2(A14,X))') &
         "SET_TEMP","TEMPERATURE"
    WRITE(UNIT=ENE_UNIT,FMT='(" 4 ",4(A14,X))') &
         "VCM_X","VCM_Y","VCM_Z","CM_KIN"
    WRITE(UNIT=ENE_UNIT,FMT='(" 5 ",4(A14,X))') &
         "OMEGA_X","OMEGA_Y","OMEGA_Z","OMEGA_KIN"
    WRITE(UNIT=ENE_UNIT,FMT='(A80)') &
	"--------------------------------------------------------------------------------"
!
  END SUBROUTINE OpenENE
!
!--------------------------------------------------------
!
  SUBROUTINE CloseCRD
    CLOSE(CRD_UNIT)
  END SUBROUTINE CloseCRD
!
!--------------------------------------------------------
!
  SUBROUTINE CloseVEL
    CLOSE(VEL_UNIT)
  END SUBROUTINE CloseVEL
!
!--------------------------------------------------------
!
  SUBROUTINE CloseENE
    CLOSE(ENE_UNIT)
  END SUBROUTINE CloseENE
!
!--------------------------------------------------------
!
  SUBROUTINE WriteCRD(GM,CTRL)
!
    TYPE(CRDS)          :: GM
    TYPE(SCFControls)   :: CTRL
    INTEGER             :: i
!
!  For now, no box info is written to crd file
!  Need to convert to external units (Angstroms)
!
#ifdef PERIODIC
    IF (CTRL%MDC%AtomWrap) THEN
       WRITE(CRD_UNIT,FMT='(10F8.3)') &
            (GM%Carts%D(1,i) * BohrsToAngstroms, &
             GM%Carts%D(2,i) * BohrsToAngstroms, &
             GM%Carts%D(3,i) * BohrsToAngstroms, &
             i=1,GM%NAtms)
    ELSE
       WRITE(CRD_UNIT,FMT='(10F8.3)') &
            (GM%Carts%D(1,i) * BohrsToAngstroms, &
             GM%Carts%D(2,i) * BohrsToAngstroms, &
             GM%Carts%D(3,i) * BohrsToAngstroms, &
             i=1,GM%NAtms)
    END IF
#else
    WRITE(CRD_UNIT,FMT='(10F8.3)') &
         (GM%Carts%D(1,i) * BohrsToAngstroms,&
          GM%Carts%D(2,i) * BohrsToAngstroms,&
          GM%Carts%D(3,i) * BohrsToAngstroms,&
          i=1,GM%NAtms)
#endif
!
  END SUBROUTINE WriteCRD
!
!--------------------------------------------------------
!
  SUBROUTINE WriteVEL(GM)
!
    TYPE(CRDS)          :: GM
    INTEGER             :: i
!
!  For now, no box info is written to vel file
!
    WRITE(VEL_UNIT,FMT='(10F8.3)') &
         (GM%Vects%D(1,i), &
          GM%Vects%D(2,i), &
          GM%Vects%D(3,i), &
          i=1,GM%NAtms)
!
  END SUBROUTINE WriteVEL
!
!--------------------------------------------------------
!
  SUBROUTINE WriteENE(GM,CTRL,MDS)
!
    TYPE(CRDS)          :: GM
    TYPE(SCFControls)   :: CTRL
    TYPE(MDState)       :: MDS
    INTEGER             :: i
!
! Need to convert to "outside" units before printing
!
    WRITE(UNIT=ENE_UNIT,FMT='(" 0 ",I14,1X,5(E14.6,2X))') &
         MDS%STEP                               , &
	 MDS%TIME * InternalTimeToSeconds       , &
	 MDS%E_POT                              , &
         MDS%E_KIN                              , &
         MDS%E_TOT                              , &
         MDS%DE_TOT
!
#ifdef PERIODIC
    WRITE(UNIT=ENE_UNIT,FMT='(" 1 ",4(E14.6,2X))') &
         (GM%PBC%TransVec(i),i=1,3),GM%PBC%CellVolume
#else
    WRITE(UNIT=ENE_UNIT,FMT='(" 1 ",4(E14.6,2X))') &
         0.0D0,0.0D0,0.0.0D0,0.0D0
#endif
!
    WRITE(UNIT=ENE_UNIT,FMT='(" 2 ",5(E14.6,2X))') &
         (MDS%VIRIAL(i),i=1,3),CTRL%MDC%PRES,MDS%PRESSURE
!
    WRITE(UNIT=ENE_UNIT,FMT='(" 3 ",2(E14.6,2X))') &
         CTRL%MDC%TEMP,MDS%TEMPERATURE
!
    WRITE(UNIT=ENE_UNIT,FMT='(" 4 ",4(E14.6,2X))') &
         (MDS%VCM(i),i=1,3),MDS%CM_KIN
!
    WRITE(UNIT=ENE_UNIT,FMT='(" 5 ",4(E14.6,2X))') &
         (MDS%OMEGA(i),i=1,3),MDS%OMEGA_KIN
!
    WRITE(UNIT=ENE_UNIT,FMT='(A80)') &
	"--------------------------------------------------------------------------------"
!
  END SUBROUTINE WriteENE
!
!--------------------------------------------------------
!
  SUBROUTINE OpenRESIN(FILENAME,TITLE)
!
    CHARACTER(LEN=DEFAULT_CHR_LEN)       :: FILENAME
    CHARACTER(LEN=DEFAULT_CHR_LEN)       :: TITLE
    INTEGER                              :: OPEN_ERROR
!
!   Do Nothing for now
!
  END SUBROUTINE OpenRESIN
!
!--------------------------------------------------------
!
  SUBROUTINE OpenRESOUT(FILENAME,TITLE,CLOBBER)
!
    CHARACTER(LEN=DEFAULT_CHR_LEN)       :: FILENAME
    CHARACTER(LEN=DEFAULT_CHR_LEN)       :: TITLE
    LOGICAL                              :: CLOBBER
    INTEGER                              :: OPEN_ERROR
!
!   Do Nothing for now
!
  END SUBROUTINE OpenRESOUT
!
!--------------------------------------------------------
!
  SUBROUTINE ReadRESIN(Ctrl,MDS,GM)
!
    TYPE(SCFControls)  :: Ctrl
    TYPE(MDState)      :: MDS
    TYPE(CRDS)         :: GM
!
!   Do Nothing for now
!    
  END SUBROUTINE ReadRESIN
!
!--------------------------------------------------------
!
  SUBROUTINE WriteRESOUT(Ctrl,MDS,GM)
!
    TYPE(SCFControls)  :: Ctrl
    TYPE(MDState)      :: MDS
    TYPE(CRDS)         :: GM
!
!   Do Nothing for now
!    
  END SUBROUTINE WriteRESOUT
!
!--------------------------------------------------------
!
  SUBROUTINE CloseRESIN
  END SUBROUTINE CloseRESIN
!
!
!--------------------------------------------------------
!
  SUBROUTINE CloseRESOUT
  END SUBROUTINE CloseRESOUT
!
END MODULE MDIO
