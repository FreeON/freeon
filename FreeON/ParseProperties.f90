!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
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
  USE MondoLogger
#ifdef NAG
  USE F90_UNIX
#endif

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
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: BEGIN_PROPERTIES  = '<BeginProperties>'
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: END_PROPERTIES    = '<EndProperties>'
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: PROP_RESPONSE     = 'Response'
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: PROP_TD_SCF       = 'TD-SCF'
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: PROP_STATIC_ALPHA = 'StcAlpha'
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: PROP_STATIC_BETA  = 'StcBeta'
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: PROP_STATIC_GAMMA = 'StcGamma'
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
    ! Close the input file.
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
    INTEGER :: ijXYZ,ikXYZ,jkXYZ,ijkXYZ,iXYZ,jXYZ,kXYZ
    LOGICAL :: IsThere,Tmp
    LOGICAL, DIMENSION( 3) :: ArrTmpA
    LOGICAL, DIMENSION( 6) :: ArrTmpB
    LOGICAL, DIMENSION(10) :: ArrTmpG
    CHARACTER(LEN=*), DIMENSION(3), PARAMETER :: Cart=(/'X','Y','Z'/)
    !-------------------------------------------------------------------
    ! Look for TD-SCF
    IF(OptKeyQ(Inp,PROP_RESPONSE,PROP_TD_SCF)) THEN
       R%TD_SCF = .TRUE.
       R%StcAlpha = .FALSE.
    ENDIF
    !
    ! Look for Static Polarizability.
    IF(OptKeyQ(Inp,PROP_RESPONSE,PROP_STATIC_ALPHA)) THEN
       R%StcAlpha = .TRUE.
       !
       ArrTmpA=.FALSE.
       IsThere=.FALSE.
       DO iXYZ=1,3
          IF(OptKeyQ(Inp,PROP_RESPONSE,Cart(iXYZ))) THEN
             ArrTmpA(iXYZ)=.TRUE.
             IsThere=.TRUE.
          ENDIF
       ENDDO
       !
       IF(IsThere) THEN
          R%AlphaAxis = ArrTmpA
          CALL MondoLog(DEBUG_MAXIMUM, "Response", "calculating polarizability along the following axes: [ "// &
            TRIM(LogicalToChar(R%AlphaAxis(1)))//" "// &
            TRIM(LogicalToChar(R%AlphaAxis(2)))//" "// &
            TRIM(LogicalToChar(R%AlphaAxis(3)))//" ]")
       ENDIF
    ENDIF
    !
    ! Look for First Static HyperPolarizability.
    IF(OptKeyQ(Inp,PROP_RESPONSE,PROP_STATIC_BETA)) THEN
       R%StcAlpha = .TRUE.
       R%StcBeta  = .TRUE.
       !
       ArrTmpA=.FALSE.
       ArrTmpB=.FALSE.
       IsThere=.FALSE.
       ijXYZ=0
       DO iXYZ=1,3
          DO jXYZ=iXYZ,3
             ijXYZ=ijXYZ+1
             IF(OptKeyQ(Inp,PROP_RESPONSE,Cart(iXYZ)//Cart(jXYZ))) THEN
                ArrTmpB(ijXYZ)=.TRUE.
                ArrTmpA(iXYZ)=.TRUE.
                ArrTmpA(jXYZ)=.TRUE.
                IsThere=.TRUE.
             ENDIF
          ENDDO
       ENDDO
       !
       IF(IsThere) THEN
          R%AlphaAxis = ArrTmpA
          R%BetaAxis  = ArrTmpB
       ENDIF
    ENDIF
    !
    ! Look for Second Static HyperPolarizability.
    IF(OptKeyQ(Inp,PROP_RESPONSE,PROP_STATIC_GAMMA)) THEN
       R%StcAlpha = .TRUE.
       R%StcBeta  = .TRUE.
       R%StcGamma = .TRUE.
       !
       ArrTmpA=.FALSE.
       ArrTmpB=.FALSE.
       ArrTmpG=.FALSE.
       IsThere=.FALSE.
       ijkXYZ =0
       ijXYZ  =0
       DO iXYZ=1,3
          jkXYZ=(-iXYZ**2+9*iXYZ-8)/2
          DO jXYZ=iXYZ,3
             ijXYZ=ijXYZ+1
             DO kXYZ=jXYZ,3
                ijkXYZ=ijkXYZ+1
                jkXYZ=jkXYZ+1
                ikXYZ=(-iXYZ**2+7*iXYZ-6)/2+kXYZ
                !  1  2  3  4  5  6
                ! XX XY XZ YY YZ ZZ
                !  1   2   3   4   5   6   7   8   9   10
                ! XXX,XXY,XXZ,XYY,XYZ,XZZ,YYY,YYZ,YZZ,ZZZ
                IF(OptKeyQ(Inp,PROP_RESPONSE,Cart(iXYZ)//Cart(jXYZ)//Cart(kXYZ))) THEN
                   ArrTmpG(ijkXYZ)=.TRUE.
                   ArrTmpB(ijXYZ)=.TRUE.
                   ArrTmpB(jkXYZ)=.TRUE.
                   ArrTmpB(ikXYZ)=.TRUE.
                   ArrTmpA(iXYZ)=.TRUE.
                   ArrTmpA(jXYZ)=.TRUE.
                   ArrTmpA(kXYZ)=.TRUE.
                   IsThere=.TRUE.
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       !
       !write(*,*) 'ArrTmpG',ArrTmpG
       !write(*,*) 'ArrTmpB',ArrTmpB
       !write(*,*) 'ArrTmpA',ArrTmpA
       IF(IsThere) THEN
          R%AlphaAxis = ArrTmpA
          R%BetaAxis  = ArrTmpB
          R%GammaAxis = ArrTmpG
       ENDIF
       !CALL MondoHalt(PRSE_ERROR,'Well tryed, but Cubic Response has not been implemented yet.')
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
    Prp%Resp%BetaAxis (:) = .TRUE.
    Prp%Resp%GammaAxis(:) = .TRUE.
    !
    ! ...
    !
  END SUBROUTINE PPrp_Init
  !
  !
END MODULE ParseProperties
