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
  !USE MLP !ONLY: MLPDriver
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
  PRIVATE :: Save_LastCPSCFCycleNbr
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
    TYPE(Controls)                             :: C
    !-------------------------------------------------------------------
    TYPE(DBL_RNK2)                             :: ETot,DMax,DIIS
    INTEGER                        , PARAMETER :: MaxCPSCFs=64 !32
    INTEGER                                    :: iCPSCF,iXYZ,jXYZ,kXYZ,ijXYZ
    ! 1 2 3
    ! X Y Z
    CHARACTER(LEN=*), DIMENSION( 3), PARAMETER :: Cart =(/'X','Y','Z'/)
    !  1  2  3  4  5  6
    ! XX XY XZ YY YZ ZZ
    CHARACTER(LEN=*), DIMENSION( 6), PARAMETER :: Cart2=(/'XX','YY','ZZ', &
         &                                                'XY','YZ','XZ'/)
    !  1   2   3   4   5   6   7   8   9   10
    ! XXX,XXY,XXZ,XYY,XYZ,XZZ,YYY,YYZ,YZZ,ZZZ
    CHARACTER(LEN=*), DIMENSION(10), PARAMETER :: Cart3=(/'XXX','XXY','XXZ', &
         &                                                'XYY','XYZ','XZZ', &
         &                                                'YYY','YYZ','YZZ', &
         &                                                'ZZZ'/)
    integer :: iii
    !-------------------------------------------------------------------
    !
    IF(.NOT. C%POpt%Resp%StcAlpha) RETURN
    !-------------------------------------------------------------------
    ! Save last SCF cycle number.

    CALL Save_LastCPSCFCycleNbr(C,'lastscfcycle')
    !
    !-------------------------------------------------------------------
    ! COMPUTE DIPOLE COMPUTE DIPOLE COMPUTE DIPOLE COMPUTE DIPOLE COMPU
    !-------------------------------------------------------------------
    !
    CALL New(C%Stat%Action,3)
    !
    ! Compute dipole moments.
    C%Stat%Action%C(1)='DipoleBuild'
    C%Stat%Action%C(2)='Dipole'
    C%Stat%Action%C(3)='All'
    CALL Invoke('MakeM',C%Nams,C%Stat,C%MPIs)
    !
    CALL Delete(C%Stat%Action)
    !
    ! Compute Multipole Moments.
    !CALL MLPDriver()
    !
    !-------------------------------------------------------------------
    ! LINEAR RESPONSE LINEAR RESPONSE LINEAR RESPONSE LINEAR RESPONSE L
    !-------------------------------------------------------------------
    !
    CALL New(C%Stat%Action,4)
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
       C%Stat%Action%C(2)=Cart(iXYZ)
       C%Stat%Action%C(3)='Dipole'
       C%Stat%Action%C(4)=Cart(iXYZ)
       !
       DO iCPSCF=0,MaxCPSCFs
          !
          ! Do a SCF cycle.
          IF(SCFCycle(iCPSCF,C%Stat%Current%I(2),C%Stat%Current%I(3), &
               &      C%Nams,C%Stat,C%Opts,C%Geos,C%Dyns,C%MPIs,ETot,DMax,   &
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
          CALL MondoHalt(DRIV_ERROR,'Failed to converge linear CPSCF in ' &
               &         //TRIM(IntToChar(MaxCPSCFs))//' CPSCF iterations.')
       ENDIF
       !
       ! Save last SCF cycle number.
       CALL Save_LastCPSCFCycleNbr(C,'lastcpscfcycle'//Cart(iXYZ))
       !
    ENDDO
    !
    CALL Delete(C%Stat%Action)
    !
    !-------------------------------------------------------------------
    ! QUADRATIC RESPONSE QUADRATIC RESPONSE QUADRATIC RESPONSE QUADRATI
    !-------------------------------------------------------------------
    !
    ! Check if we need to compute the Quadratic Response.
    IF(.NOT.C%POpt%Resp%StcBeta) RETURN
    !
    CALL New(C%Stat%Action,4)
    !
    ! Let's compute the quadratic response.
    ijXYZ=0
    DO iXYZ=1,3
       !
       DO jXYZ=iXYZ,3
          !
          ijXYZ=ijXYZ+1
          !
          ! Do we need to compute the hyperpolarizability along this axis?
          IF(.NOT.C%POpt%Resp%BetaAxis(ijXYZ)) CYCLE
          !
          ! Allocate space for convergence statistics
          CALL New(ETot,(/MaxCPSCFs,C%Geos%Clones/),(/0,1/))
          CALL New(DMax,(/MaxCPSCFs,C%Geos%Clones/),(/0,1/))
          CALL New(DIIS,(/MaxCPSCFs,C%Geos%Clones/),(/0,1/))
          !
          C%Opts%Guess=GUESS_EQ_DIPOLE!GUESS_EQ_NOGUESS
          C%Stat%Action%C(2)=Cart(iXYZ)//Cart(jXYZ)
          C%Stat%Action%C(3)='Dipole'
          C%Stat%Action%C(4)=Cart(jXYZ)
          !
          DO iCPSCF=0,MaxCPSCFs
             !
             ! Do a SCF cycle.
             IF(SCFCycle(iCPSCF,C%Stat%Current%I(2),C%Stat%Current%I(3), &
                  &      C%Nams,C%Stat,C%Opts,C%Geos,C%Dyns,C%MPIs,ETot,DMax,   &
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
             CALL MondoHalt(DRIV_ERROR,'Failed to converge quadratic CPSCF in ' &
                  &         //TRIM(IntToChar(MaxCPSCFs))//' CPSCF iterations.')
          ENDIF
          !
          ! Save last CPSCF cycle number.
          CALL Save_LastCPSCFCycleNbr(C,'lastcpscfcycle'//Cart(iXYZ)//Cart(jXYZ))
          !
       ENDDO
    ENDDO
    !
    !
    CALL Delete(C%Stat%Action)
    !
    !-------------------------------------------------------------------
    ! CUBIC RESPONSE CUBIC RESPONSE CUBIC RESPONSE CUBIC RESPONSE CUBIC 
    !-------------------------------------------------------------------
    !
    ! Check if we need to compute the Quadratic Response.
    IF(.NOT.C%POpt%Resp%StcGamma) RETURN
    !
    CALL New(C%Stat%Action,4)
    !
    ! Let's compute the quadratic response.
    ijXYZ=0
    DO iXYZ=1,3
       !
       DO jXYZ=iXYZ,3
          !
          DO kXYZ=jXYZ,3
             !
             ijXYZ=ijXYZ+1
             !
             ! Do we need to compute the hyperpolarizability along this axis?
             IF(.NOT.C%POpt%Resp%GammaAxis(ijXYZ)) CYCLE
             !
             ! Allocate space for convergence statistics
             CALL New(ETot,(/MaxCPSCFs,C%Geos%Clones/),(/0,1/))
             CALL New(DMax,(/MaxCPSCFs,C%Geos%Clones/),(/0,1/))
             CALL New(DIIS,(/MaxCPSCFs,C%Geos%Clones/),(/0,1/))
             !
             C%Opts%Guess=GUESS_EQ_DIPOLE!GUESS_EQ_NOGUESS
             C%Stat%Action%C(2)=Cart(iXYZ)//Cart(jXYZ)//Cart(kXYZ)
             C%Stat%Action%C(3)='Dipole'
             C%Stat%Action%C(4)=Cart(kXYZ)
             !
             DO iCPSCF=0,MaxCPSCFs
                !
                ! Do a SCF cycle.
                IF(SCFCycle(iCPSCF,C%Stat%Current%I(2),C%Stat%Current%I(3), &
                     &      C%Nams,C%Stat,C%Opts,C%Geos,C%Dyns,C%MPIs,ETot,DMax,   &
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
                CALL MondoHalt(DRIV_ERROR,'Failed to converge cubic CPSCF in ' &
                     &         //TRIM(IntToChar(MaxCPSCFs))//' CPSCF iterations.')
             ENDIF
             !
             ! Save last CPSCF cycle number.
             CALL Save_LastCPSCFCycleNbr(C,'lastcpscfcycle'//Cart(iXYZ)//Cart(jXYZ)//Cart(kXYZ))
             !
          ENDDO
       ENDDO
    ENDDO
    !
    !
    CALL Delete(C%Stat%Action)
    !
  END SUBROUTINE CPSCF
  !
  !
  SUBROUTINE Save_LastCPSCFCycleNbr(C,Name)
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(Controls)  , INTENT(IN) :: C
    CHARACTER(LEN=*), INTENT(IN) :: Name
    !-------------------------------------------------------------------
    INTEGER                      :: iClone!,iCPSCF
    !-------------------------------------------------------------------
    !
    HDFFileID=OpenHDF(C%Nams%HFile)
    HDF_CurrentID=HDFFileID
    !DO iClone=1,C%Geos%Clones
       !HDF_CurrentID=InitHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
       HDF_CurrentID=OpenHDFGroup(HDFFileID,'Clone #'//TRIM(IntToChar(1)))
       CALL Put(C%Stat%Current%I(1),TRIM(Name))
       CALL CloseHDFGroup(HDF_CurrentID)
    !ENDDO
    CALL CloseHDF(HDFFileID)
    !
  END SUBROUTINE Save_LastCPSCFCycleNbr
  !
  !
END MODULE Response




