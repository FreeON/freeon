MODULE ParseProperties
!H=================================================================================
!H MODULE ParseProperties
!H This MODULE contains:
!H  PUBLIC:
!H  o SUB 
!H
!H  PRIVATE:
!H  o SUB 
!H
!H  OPTIONS:
!H  DEBUGING: 
!H  INFO    : 
!H
!H  Comments:
!H
!H---------------------------------------------------------------------------------
  USE InOut
  USE ControlStructures
#ifdef NAG
  USE F90_UNIX
#endif
  !
  IMPLICIT NONE
  PRIVATE
!---------------------------------------------------------------------------------
! PUBLIC DECLARATIONS
!---------------------------------------------------------------------------------
  PUBLIC :: LoadPropertyOptions
  !
!---------------------------------------------------------------------------------
! PRIVATE DECLARATIONS
!---------------------------------------------------------------------------------
  !PRIVATE :: 
  !
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: BEGIN_PROPERTIES  = '<BeginProperties>'
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: END_PROPERTIES    = '<EndProperties>'
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: PROP_RESPONSE     = 'Response'
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: PROP_STATIC_ALPHA = 'StcAlpha'
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: PROP_STATIC_BETA  = 'StcBeta'
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: PROP_STATIC_GAMMA = 'StcGamma'
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: PROP_X_DIRECTION  = 'X'
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: PROP_Y_DIRECTION  = 'Y'
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: PROP_Z_DIRECTION  = 'Z'
  !
CONTAINS
  !
  !
  SUBROUTINE LoadPropertyOptions(N,Prp)
!H---------------------------------------------------------------------------------
!H SUBROUTINE LoadPropertyOptions(N,Prp)
!H  
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(FileNames)    :: N
    TYPE(PropOpts )    :: Prp
    !-------------------------------------------------------------------

    !-------------------------------------------------------------------
    !
    CALL OpenASCII(N%IFile,Inp)
    !
    ! Initialize Properties
    CALL PPrp_Init(Prp)
    !
    !IF(.NOT.FindKey(BEGIN_PROPERTIES,Inp)) RETURN
    !
    ! Parse Response.
    CALL PPrp_Response(Prp%Resp)
    !
    ! Parse ...
    !
    !
    CLOSE(UNIT=Inp,STATUS='KEEP')
    !
  END SUBROUTINE LoadPropertyOptions
  !
  !
  SUBROUTINE PPrp_Response(R)
!H---------------------------------------------------------------------------------
!H SUBROUTINE PPrp_Response(R)
!H  
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(RespOpts) :: R
    !-------------------------------------------------------------------
    !
    ! Look for Static Polarizability.
    IF(OptKeyQ(Inp,PROP_RESPONSE,PROP_STATIC_ALPHA)) THEN
       R%StcAlpha = .TRUE.
       IF(    OptKeyQ(Inp,PROP_RESPONSE,PROP_X_DIRECTION).OR. &
            & OptKeyQ(Inp,PROP_RESPONSE,PROP_Y_DIRECTION).OR. &
            & OptKeyQ(Inp,PROP_RESPONSE,PROP_Z_DIRECTION)) THEN
          R%AlphaAxis(:) = .FALSE.
          IF(OptKeyQ(Inp,PROP_RESPONSE,PROP_X_DIRECTION)) R%AlphaAxis(1) = .TRUE.
          IF(OptKeyQ(Inp,PROP_RESPONSE,PROP_Y_DIRECTION)) R%AlphaAxis(2) = .TRUE.
          IF(OptKeyQ(Inp,PROP_RESPONSE,PROP_Z_DIRECTION)) R%AlphaAxis(3) = .TRUE.
       ENDIF
    ENDIF
    !
    ! Look for First Static HyperPolarizability.
    IF(OptKeyQ(Inp,PROP_RESPONSE,PROP_STATIC_BETA )) THEN 
       R%StcAlpha = .TRUE.
       R%StcBeta  = .TRUE.
       CALL MondoHalt(PRSE_ERROR,'Well tryed, but Quadratic Response has not been implemented yet.')
    ENDIF
    !
    ! Look for Second Static HyperPolarizability.
    IF(OptKeyQ(Inp,PROP_RESPONSE,PROP_STATIC_GAMMA)) THEN
       R%StcAlpha = .TRUE.
       R%StcBeta  = .TRUE.
       R%StcGamma = .TRUE.
       CALL MondoHalt(PRSE_ERROR,'Well tryed, but Cubic Response has not been implemented yet.')
    ENDIF
    !
    !write(*,*) 'PPrp_Response: R%StcAlpha=',R%StcAlpha
    !write(*,*) 'PPrp_Response: R%StcBeta =',R%StcBeta 
    !write(*,*) 'PPrp_Response: R%StcGamma=',R%StcGamma
    !
  END SUBROUTINE PPrp_Response
  !
  !
  SUBROUTINE PPrp_Init(Prp)
!H---------------------------------------------------------------------------------
!H SUBROUTINE PPrp_Init(Prp)
!H  Initialize the Property Options.
!H---------------------------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------------
    TYPE(PropOpts) :: Prp
    !-------------------------------------------------------------------
    ! Initialize the static polarizabilities.
    Prp%Resp%StcAlpha = .FALSE.
    Prp%Resp%StcBeta  = .FALSE.
    Prp%Resp%StcGamma = .FALSE.
    ! Initialize the polarizability axis.
    Prp%Resp%AlphaAxis(:) = .TRUE.
    ! 
    ! ...
    !
  END SUBROUTINE PPrp_Init
  !
  !

END MODULE ParseProperties
