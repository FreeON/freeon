MODULE Response
!H=================================================================================
!H MODULE Response
!H
!H  OPTIONS:
!H  DEBUGING: 
!H  INFO    : 
!H
!H Comment:
!H
!H Ref:
!H
!H---------------------------------------------------------------------------------
  !
  USE Parse
  USE InOut
  USE LinAlg 
  USE GlobalObjects
  USE SCFKeys
  USE Overlay
  USE SCFKeys
  USE PunchHDF
  USE Numerics
  USE OptionKeys
  USE Functionals
  USE ControlStructures
  USE NEB
  USE SCFs
  USE SetXYZ 
  USE MLP !ONLY: MLPDriver
  !
  IMPLICIT NONE 
  PRIVATE
  !
!--------------------------------------------------------------------------------- 
! PUBLIC DECLARATIONS
!--------------------------------------------------------------------------------- 
  PUBLIC  :: CPSCF
  !
!--------------------------------------------------------------------------------- 
! PRIVATE DECLARATIONS
!--------------------------------------------------------------------------------- 
  PRIVATE :: Save_LastSCFCycleNbr
  !
CONTAINS
  !
  SUBROUTINE CPSCF(C)
!H---------------------------------------------------------------------------------
!H SUBROUTINE CPSCF(C)
!H  Main driver to solve the CPSCF equations.
!H
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(Controls)                            :: C
    !-------------------------------------------------------------------
    TYPE(DBL_RNK2)                            :: ETot,DMax,DIIS
    INTEGER                       , PARAMETER :: MaxCPSCFs=64 !32
    INTEGER                                   :: iCPSCF,iXYZ
    CHARACTER(LEN=*), DIMENSION(3), PARAMETER :: Cart=(/'X','Y','Z'/)
    !-------------------------------------------------------------------
    !
    ! Check if we need to compute a Response.
    IF(.NOT.C%POpt%Resp%StcAlpha) RETURN
    !
    ! Save last SCF cycle number.
    CALL Save_LastSCFCycleNbr(C)
    !
    CALL New(C%Stat%Action,3)
    !
    ! Compute dipole moments.
    C%Stat%Action%C(1)='DipoleBuild'
    C%Stat%Action%C(2)='Dipole'
    C%Stat%Action%C(3)='All'
    CALL Invoke('MakeM',C%Nams,C%Stat,C%MPIs)
    !
    ! Compute Multipole Moments.
    !CALL MLPDriver()
    !    
    ! Let's compute the response.
    DO iXYZ=1,3
       !
       ! Do we need to compute the polarizability along this axis?
       IF(.NOT.C%POpt%Resp%AlphaAxis(iXYZ)) CYCLE
       !
       ! Allocate space for convergence statistics
       CALL New(ETot,(/MaxCPSCFs,C%Geos%Clones/),(/0,1/))
       CALL New(DMax,(/MaxCPSCFs,C%Geos%Clones/),(/0,1/))
       CALL New(DIIS,(/MaxCPSCFs,C%Geos%Clones/),(/0,1/))
       !
       C%Opts%Guess=GUESS_EQ_DIPOLE
       C%Stat%Action%C(2)='Dipole'
       C%Stat%Action%C(3)=Cart(iXYZ)
       !
       DO iCPSCF=0,MaxCPSCFs
          !
          ! Do an SCF cycle.
          IF(SCFCycle(iCPSCF,C%Stat%Current%I(2),C%Stat%Current%I(3), &
               &      C%Nams,C%Stat,C%Opts,C%Geos,C%MPIs,ETot,DMax,   &
               &      DIIS,CPSCF_O=.TRUE.))THEN
             EXIT
          ENDIF
       ENDDO
       !
       ! Free memory.
       CALL Delete(ETot)
       CALL Delete(DMax)
       CALL Delete(DIIS)
       IF(iCPSCF.GT.MaxCPSCFs+1) THEN
          CALL MondoHalt(DRIV_ERROR,'Failed to converge CPSCF in ' &
               &         //TRIM(IntToChar(MaxCPSCFs))//' CPSCF iterations.')
       ENDIF
    ENDDO
    !
    !
    CALL Delete(C%Stat%Action)
    !
  END SUBROUTINE CPSCF
  !
  !
  SUBROUTINE Save_LastSCFCycleNbr(C)
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(Controls) :: C
    !-------------------------------------------------------------------
    INTEGER        :: iClone
    !-------------------------------------------------------------------
    !
    HDFFileID=OpenHDF(C%Nams%HFile)
    HDF_CurrentID=HDFFileID
    !DO iClone=1,C%Geos%Clones
       !HDF_CurrentID=InitHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
       HDF_CurrentID=OpenHDFGroup(HDFFileID,'Clone #'//TRIM(IntToChar(1)))
       CALL Put(C%Stat%Current%I(1),'lastscfcycle')
       !CALL Get(iCPSCF,'lastscfcycle')
       !write(*,*) 'LastSCFCycle=',iCPSCF
       CALL CloseHDFGroup(HDF_CurrentID)
    !ENDDO
    CALL CloseHDF(HDFFileID)
    !
  END SUBROUTINE Save_LastSCFCycleNbr
  !
  !
END MODULE Response




