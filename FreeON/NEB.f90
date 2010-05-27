!------------------------------------------------------------------------------
!    This code is part of the FreeON suite of programs for linear scaling
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
!    to return derivative works to the FreeON group for review, and possible
!    dissemination in future releases.
!------------------------------------------------------------------------------

#include "MondoConfig.h"

MODULE NEB
  !===============================================================================
  ! Module for calculating reaction (minimum energy) paths between known
  ! reactant and product states.  This module impliments the climbing image
  ! NEB so that the highest energy image will converge to a saddle point.
  !
  ! Module written by Graeme Henkelman and Matt Challacombe
  ! Email: graeme@lanl.gov
  !
  ! NEB References:
  !  H. Jonsson, G. Mills, and K.W. Jacobsen, "Nudged elastic band method for
  !    finding minimum energy paths of transitions," in Classical and Quantum
  !    Dynamics in Condensed Phase Simulations, edited by B.J.Berne, G.
  !    Ciccotti, and D. F. Coker (World Scientific, Singapore, 1998), p.385
  !  G. Henkelman, B.P. Uberuaga, and H. Jonsson, "A climbing-image NEB method
  !    for finding saddle points and minimum energy paths", J. Chem. Phys 113,
  !    9901 (2000).
  !  G. Henkelman, H. Jonsson, "Improved tangent estimate in the nudged elastic
  !    band method for finding minimum energy paths and saddle points", J. Chem.
  !    Phys. 113, 9978 (2000)
  !
  !===============================================================================
  USE InOut
  USE ls_rmsd
  USE PrettyPrint
  USE ControlStructures
  USE Order
  USE MondoLogger

  IMPLICIT NONE

  SAVE

CONTAINS
  !===============================================================================
  ! Initialize the NEB by generating an linear interpolation between intial
  ! and final states.
  !===============================================================================
  SUBROUTINE NEBInit(G)
    TYPE(Geometries)                            :: G
    REAL(DOUBLE),DIMENSION(3,G%Clone(0)%NAtms)  :: ReactionVector
    REAL(DOUBLE)                                :: ImageFraction
    INTEGER                                     :: iCLONE, j, RPBeginClone, RPEndClone

    !Initialize each clone to initial state then interpolate Cartesian coordinates
#if defined(NEB_DEBUG)
    CALL MondoLog(DEBUG_NONE, "NEBInit", "starting...")

    CALL MondoLog(DEBUG_NONE, "NEBInit", "Reactant Clone 0 (in input file units)")
    DO j=1, G%Clone(0)%NAtms
      CALL MondoLog(DEBUG_NONE, "NEBInit", TRIM(G%Clone(0)%AtNam%C(j))//" "// &
        TRIM(DblToChar(G%Clone(0)%Carts%D(1,j)))//" "// &
        TRIM(DblToChar(G%Clone(0)%Carts%D(2,j)))//" "// &
        TRIM(DblToChar(G%Clone(0)%Carts%D(3,j)))//" "// &
        TRIM(IntToChar(G%Clone(0)%CConstrain%I(j))))
    ENDDO

    CALL MondoLog(DEBUG_NONE, "NEBInit", "Product Clone "//TRIM(IntToChar(G%Clones+1))//" (in input file units)")
    DO j=1, G%Clone(G%Clones+1)%NAtms
      CALL MondoLog(DEBUG_NONE, "NEBInit", TRIM(G%Clone(G%Clones+1)%AtNam%C(j))//" "// &
        TRIM(DblToChar(G%Clone(G%Clones+1)%Carts%D(1,j)))//" "// &
        TRIM(DblToChar(G%Clone(G%Clones+1)%Carts%D(2,j)))//" "// &
        TRIM(DblToChar(G%Clone(G%Clones+1)%Carts%D(3,j)))//" "// &
        TRIM(IntToChar(G%Clone(G%Clones+1)%CConstrain%I(j))))
    ENDDO
#endif

    ! Check some things for internal consistency.
    DO j = 1, G%Clone(0)%NAtms
      IF(G%Clone(0)%CConstrain%I(j) /= G%Clone(G%Clones+1)%CConstrain%I(j)) THEN
        CALL MondoLog(DEBUG_NONE, "NEBInit", "constrain on atom "//TRIM(IntToChar(j))//" is different between reactant and product")
        CALL Halt("Constrain mismatch on input")
      ENDIF
    ENDDO

    ! Calculate reaction path vector. We might have alrady read some clone
    ! geometries from input. We will interpolate along the reaction path vector
    ! between clones given in input.
    RPBeginClone = 0
    RPEndClone = 1

    DO WHILE(RPEndClone < G%Clones+1)

      DO iCLONE = RPBeginClone+1, G%Clones+1
        IF(G%Clone(iCLONE)%NAtms > 0) THEN
          ! This clone is already allocated. Construct reaction path.
          RPEndClone = iCLONE
          IF(RPEndClone > RPBeginClone+1) THEN
            CALL MondoLog(DEBUG_NONE, "NEBInit", "interpolating between clones "//TRIM(IntToChar(RPBeginClone))//" and "//TRIM(IntToChar(RPEndClone)))
          ELSE
            CALL MondoLog(DEBUG_NONE, "NEBInit", "adjacent clones "//TRIM(IntToChar(RPBeginClone))//" and " &
              //TRIM(IntToChar(RPEndClone))//" given in input file, no need for interpolation")
          ENDIF

#if defined(NEB_DEBUG)
          CALL MondoLog(DEBUG_NONE, "NEBInit", "Clone "//TRIM(IntToChar(RPBeginClone))//" (in input file units)")
          DO j=1, G%Clone(RPBeginClone)%NAtms
            CALL MondoLog(DEBUG_NONE, "NEBInit", TRIM(G%Clone(RPBeginClone)%AtNam%C(j))//" "// &
              TRIM(DblToChar(G%Clone(RPBeginClone)%Carts%D(1,j)))//" "// &
              TRIM(DblToChar(G%Clone(RPBeginClone)%Carts%D(2,j)))//" "// &
              TRIM(DblToChar(G%Clone(RPBeginClone)%Carts%D(3,j)))//" "// &
              TRIM(IntToChar(G%Clone(RPBeginClone)%CConstrain%I(j))))
          ENDDO

          CALL MondoLog(DEBUG_NONE, "NEBInit", "Clone "//TRIM(IntToChar(RPEndClone))//" (in input file units)")
          DO j=1, G%Clone(RPEndClone)%NAtms
            CALL MondoLog(DEBUG_NONE, "NEBInit", TRIM(G%Clone(RPEndClone)%AtNam%C(j))//" "// &
              TRIM(DblToChar(G%Clone(RPEndClone)%Carts%D(1,j)))//" "// &
              TRIM(DblToChar(G%Clone(RPEndClone)%Carts%D(2,j)))//" "// &
              TRIM(DblToChar(G%Clone(RPEndClone)%Carts%D(3,j)))//" "// &
              TRIM(IntToChar(G%Clone(RPEndClone)%CConstrain%I(j))))
          ENDDO
#endif

          EXIT
        ENDIF
      ENDDO

      ReactionVector = G%Clone(RPEndClone)%Carts%D-G%Clone(RPBeginClone)%Carts%D

#ifdef NEB_DEBUG
      CALL MondoLog(DEBUG_NONE, "NEBInit", "Reaction vector (in input file units)")
      DO j=1, G%Clone(RPBeginClone)%NAtms
        CALL MondoLog(DEBUG_NONE, "NEBInit", "RV["//TRIM(IntToChar(j))//"] = [ "// &
          TRIM(DblToChar(ReactionVector(1,j)))//" "// &
          TRIM(DblToChar(ReactionVector(2,j)))//" "// &
          TRIM(DblToChar(ReactionVector(3,j)))//" ]")
      ENDDO
#endif

      ! Initialize reaction path.
      DO iCLONE = RPBeginClone+1, RPEndClone-1
        ImageFraction = DBLE(iCLONE-RPBeginClone)/DBLE(RPEndClone-RPBeginClone)

        ! Allocate clone.
        G%Clone(iCLONE)%NAtms = G%Clone(0)%NAtms
        G%Clone(iCLONE)%NKind = G%Clone(0)%NKind
        CALL New(G%Clone(iCLONE))

        ! Initialize this clone with the reactant as the template.
        CALL SetEq(G%Clone(iCLONE), G%Clone(0))

        ! Fix constrains.
        G%Clone(iCLONE)%CConstrain%I = G%Clone(0)%CConstrain%I

        ! Linear interpolation of path.
        CALL MondoLog(DEBUG_NONE, "NEBInit", "linear interpolation for clone "//TRIM(IntToChar(iCLONE)))
        G%Clone(iCLONE)%Carts%D = G%Clone(RPBeginClone)%Carts%D + ImageFraction*ReactionVector

        ! Set everything else to 0 in this clone.
        G%Clone(iCLONE)%Velocity%D = Zero
        G%Clone(iCLONE)%Fext%D = Zero
        G%Clone(iCLONE)%Gradients%D = Zero

#ifdef NEB_DEBUG
        CALL MondoLog(DEBUG_NONE, "NEBInit", "Clone "//TRIM(IntToChar(iCLONE))//" (in input file units)")
        DO j=1, G%Clone(iCLONE)%NAtms
          CALL MondoLog(DEBUG_NONE, "NEBInit", TRIM(G%Clone(iCLONE)%AtNam%C(j))//" "// &
            TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(1,j)))//" "// &
            TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(2,j)))//" "// &
            TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(3,j)))//" "// &
            TRIM(IntToChar(G%Clone(iCLONE)%CConstrain%I(j))))
        ENDDO
#endif

      ENDDO

      RPBeginClone = RPEndClone
      RPEndClone = RPBeginClone+1
    ENDDO

    CALL MondoLog(DEBUG_NONE, "NEBInit", "done NEBInit")
  END SUBROUTINE NEBInit

  SUBROUTINE NEBPurify(G,Init_O)
    TYPE(Geometries)                      :: G
    LOGICAL,OPTIONAL                      :: Init_O
    LOGICAL                               :: Init
    INTEGER                               :: I,iCLONE,bCLONE,eCLONE,nCLONE,J
    TYPE(DBL_RNK2), DIMENSION(G%Clones+1) :: GTmp
    INTEGER, DIMENSION(G%Clones+1)        :: I2
    REAL(DOUBLE),DIMENSION(G%Clones+1)    :: R2
    REAL(DOUBLE),DIMENSION(3,3)           :: U
    REAL(DOUBLE),DIMENSION(3)             :: Center1,Center2
    REAL(DOUBLE)                          :: Error
    CHARACTER(LEN=4*DCL)                  :: Mssg

    IF(PRESENT(Init_O))THEN
       Init=.TRUE.
    ELSE
       Init=.FALSE.
    ENDIF

#ifdef NEB_DEBUG
    CALL MondoLog(DEBUG_NONE, "NEB", "Init = "//TRIM(LogicalToChar(Init)))
#endif

    IF(Init)THEN
       bCLONE=G%Clones+1
       eCLONE=G%Clones+1
       ! Check for stupid input
       DO I=1,G%Clone(0)%NAtms
          IF(G%Clone(0)%AtNum%D(I).NE.G%Clone(G%Clones+1)%AtNum%D(I))THEN
             CALL MondoHalt(NEBS_ERROR,'Ordering of Reactant and Product is different!')
          ENDIF
       ENDDO
    ELSE
       bCLONE=1
       eCLONE=G%Clones+1
    ENDIF

    !IF(.NOT. Init) THEN
      CALL MondoLog(DEBUG_NONE, "NEB", "not purifying")
      RETURN
    !ENDIF

#ifdef NEB_DEBUG
    CALL MondoLog(DEBUG_NONE, "NEB", "bCLONE = "//TRIM(IntToChar(bCLONE)))
    CALL MondoLog(DEBUG_NONE, "NEB", "eCLONE = "//TRIM(IntToChar(eCLONE)))
    CALL MondoLog(DEBUG_NONE, "NEB", "CLONE 0 before anything = ")
    DO I = 1, G%Clone(0)%NAtms
      CALL MondoLog(DEBUG_NONE, "NEB", "R["//TRIM(IntToChar(I))//"] = [ "// &
        TRIM(DblToChar(G%Clone(0)%Carts%D(1,I)*AUToAngstroms))//" "// &
        TRIM(DblToChar(G%Clone(0)%Carts%D(2,I)*AUToAngstroms))//" "// &
        TRIM(DblToChar(G%Clone(0)%Carts%D(3,I)*AUToAngstroms))//" ]", "Clone "//TRIM(IntToChar(0)))
    ENDDO
#endif

!!$    ! Scale the coordinates by Z
!!$    DO I=1,G%Clone(0)%NAtms
!!$       G%Clone(0)%Carts%D(:,I)=G%Clone(0)%Carts%D(:,I)*G%Clone(0)%AtNum%D(I)
!!$    ENDDO

#ifdef NEB_DEBUG
    !CALL MondoLog(DEBUG_NONE, "NEB", "CLONE 0 AFTER SCALING = ")
    !DO I = 1, G%Clone(0)%NAtms
    !  CALL MondoLog(DEBUG_NONE, "NEB", "R["//TRIM(IntToChar(I))//"] = [ "// &
    !    TRIM(DblToChar(G%Clone(0)%Carts%D(1,I)*AUToAngstroms))//" "// &
    !    TRIM(DblToChar(G%Clone(0)%Carts%D(2,I)*AUToAngstroms))//" "// &
    !    TRIM(DblToChar(G%Clone(0)%Carts%D(3,I)*AUToAngstroms))//" ]", "Clone "//TRIM(IntToChar(0)))
    !ENDDO
#endif

    ! Constraints over-ride purification
    DO I=1,G%Clone(0)%NAtms
       IF(G%Clone(0)%CConstrain%I(I)/=0)THEN
          ! No RMSD alignment with constraints
          GOTO 101 ! can still re-order based on RMSD
       ENDIF
    ENDDO

    ! Translate and rotate each clone to minimize the rmsd relative to clone zero
    DO iCLONE=bCLONE,eCLONE
#ifdef NEB_DEBUG
      CALL MondoLog(DEBUG_NONE, "NEB", "purifying clone "//TRIM(IntToChar(iclone)))
#endif

!!$       ! Scale the coordinates by Z
!!$       DO I=1,G%Clone(0)%NAtms
!!$          G%Clone(iCLONE)%Carts%D(:,I)=G%Clone(iCLONE)%Carts%D(:,I)*G%Clone(iCLONE)%AtNum%D(I)
!!$       ENDDO

       ! Find the transformation that minimizes the RMS deviation between the
       ! reactants (clone 0), the clones (1-N) and the products (N+1)
       CALL RMSD(G%Clone(0)%NAtms,G%Clone(iCLONE)%Carts%D,G%Clone(0)%Carts%D,  &
            1, U, Center2, Center1, error )! , calc_g, grad)
#ifdef NEB_DEBUG
       CALL MondoLog(DEBUG_NONE, "NEB", "Center   1: "// &
         TRIM(DblToChar(Center1(1)))//" "// &
         TRIM(DblToChar(Center1(2)))//" "// &
         TRIM(DblToChar(Center1(3))))
       CALL MondoLog(DEBUG_NONE, "NEB", "Center   2: "// &
         TRIM(DblToChar(Center2(1)))//" "// &
         TRIM(DblToChar(Center2(2)))//" "// &
         TRIM(DblToChar(Center2(3))))
       CALL MondoLog(DEBUG_NONE, "NEB", "Center 1-2: "// &
         TRIM(DblToChar(Center1(1)-Center2(1)))//" "// &
         TRIM(DblToChar(Center1(2)-Center2(2)))//" "// &
         TRIM(DblToChar(Center1(3)-Center2(3))))
#endif
       IF(Init)THEN
          ! Translate the reactants JUST ONCE to C1
          DO I=1,G%Clone(0)%NAtms
             G%Clone(0)%Carts%D(:,I)=G%Clone(0)%Carts%D(:,I)-Center1
          ENDDO
       ENDIF
#ifdef NEB_DEBUG
       CALL MondoLog(DEBUG_NONE, "NEB", "CLONE 0 AFTER TRANSLATION")
       DO I = 1, G%Clone(0)%NAtms
         CALL MondoLog(DEBUG_NONE, "NEB", "R["//TRIM(IntToChar(I))//"] = [ "// &
           TRIM(DblToChar(G%Clone(0)%Carts%D(1,I)*AUToAngstroms))//" "// &
           TRIM(DblToChar(G%Clone(0)%Carts%D(2,I)*AUToAngstroms))//" "// &
           TRIM(DblToChar(G%Clone(0)%Carts%D(3,I)*AUToAngstroms))//" ]", "Clone "//TRIM(IntToChar(0)))
       ENDDO

       CALL MondoLog(DEBUG_NONE, "NEB", "CLONE "//TRIM(IntToChar(iCLONE))//" BEFORE TRANSLATION")
       DO I = 1, G%Clone(iCLONE)%NAtms
         CALL MondoLog(DEBUG_NONE, "NEB", "R["//TRIM(IntToChar(I))//"] = "// &
           TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(1,I)*AUToAngstroms))//" "// &
           TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(2,I)*AUToAngstroms))//" "// &
           TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(3,I)*AUToAngstroms))//" ]", "Clone "//TRIM(IntToChar(iCLONE)))
       ENDDO
#endif
       ! Translation to C2 ...
       DO I=1,G%Clone(0)%NAtms
          G%Clone(iCLONE)%Carts%D(:,I)=G%Clone(iCLONE)%Carts%D(:,I)-Center2
       ENDDO
       ! ... and rotation
       DO I=1,G%Clone(0)%NAtms
          G%Clone(iCLONE)%Carts%D(:,I)=MATMUL(U,G%Clone(iCLONE)%Carts%D(:,I))
       ENDDO
#ifdef NEB_DEBUG
       CALL MondoLog(DEBUG_NONE, "NEB", "CLONE "//TRIM(IntToChar(iCLONE))//" AFTER TRANSLATION AND ROTATION")
       DO I = 1, G%Clone(iCLONE)%NAtms
         CALL MondoLog(DEBUG_NONE, "NEB", "R["//TRIM(IntToChar(I))//"] = [ "// &
           TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(1,I)*AUToAngstroms))//" "// &
           TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(2,I)*AUToAngstroms))//" "// &
           TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(3,I)*AUToAngstroms))//" ]", "Clone "//TRIM(IntToChar(iCLONE)))
       ENDDO
#endif
    ENDDO

!!$    ! Un-scale the coordinates by Z
!!$    DO I=1,G%Clone(0)%NAtms
!!$       G%Clone(0)%Carts%D(:,I)=G%Clone(0)%Carts%D(:,I)/G%Clone(0)%AtNum%D(I)
!!$    ENDDO
!!$    DO iCLONE=bCLONE,eCLONE
!!$       DO I=1,G%Clone(0)%NAtms
!!$          G%Clone(iCLONE)%Carts%D(:,I)=G%Clone(iCLONE)%Carts%D(:,I)/G%Clone(iCLONE)%AtNum%D(I)
!!$       ENDDO
!!$    ENDDO

101 IF(Init) RETURN

    ! Compute RMSD from first clone
    R2(:)=Zero
    DO iCLONE=bCLONE,eCLONE
       R2(iCLONE)=Zero
       DO I=1,G%Clone(0)%NAtms
          DO J=1,3
             R2(iCLONE)=R2(iCLONE)+(G%Clone(iCLONE)%Carts%D(J,I)-G%Clone(0)%Carts%D(J,I))**2
          ENDDO
       ENDDO
       R2(iCLONE)=SQRT(R2(iCLONE))/G%Clone(0)%NAtms
    ENDDO

    ! Order based on RMSD
    nCLONE=eCLONE-bCLONE+1
    DO I=1,nCLONE
       I2(I)=I
       CALL New(GTmp(I),(/3,G%Clone(0)%NAtms/))
       GTmp(I)%D=G%Clone(I)%Carts%D
    ENDDO
    CALL DblIntSort77(nCLONE,R2,I2,2)

    DO I=bCLONE, eCLONE
      CALL MondoLog(DEBUG_NONE, "NEBPurify", "Clone "//TRIM(IntToChar(I))//": I2 = "//TRIM(IntToChar(I2(I)))//" R2 = "//TRIM(DblToChar(R2(I))))
    ENDDO

    DO I=1,nCLONE
       G%Clone(I)%Carts%D=GTmp(I2(I))%D
    ENDDO
    DO I=1,nCLONE
       CALL Delete(GTmp(I))
    ENDDO

    Mssg='RMSDs = '
    DO I=1,G%Clones
       Mssg=TRIM(Mssg)//' '//TRIM(DblToShrtChar(R2(I)))//','
    ENDDO
    Mssg=TRIM(Mssg)//' '//TRIM(DblToShrtChar(R2(G%Clones+1)))

    CALL MondoLog(DEBUG_NONE, "FreeON", Mssg, "NEBPurify("//TRIM(IntToChar(G%Clone(1)%Confg))//')')

  END SUBROUTINE NEBPurify

  !===============================================================================
  ! Project out the force along the band and add spring forces along the band.
  !===============================================================================
  SUBROUTINE NEBForce(G,O)
    TYPE(Geometries)                                           :: G
    TYPE(Options)                                              :: O
    INTEGER                                                    :: iCLONE, j, UMaxClone
    LOGICAL                                                    :: UPMinus, UPPlus
    REAl(DOUBLE)                                               :: UMin, UMax, UMinus, UPlus, fMax
    REAL(DOUBLE)                                               :: magnitude, RMinusMagnitude, RPlusMagnitude
    REAL(DOUBLE)                                               :: fMagnitude, fSpringMagnitude
    REAL(DOUBLE)                                               :: RMSd, theta, thetaNormMinus, thetaNormPlus
    REAL(DOUBLE), DIMENSION(3, G%Clone(0)%NAtms, 0:G%Clones+1) :: N
    REAL(DOUBLE), DIMENSION(3, G%Clone(0)%NAtms)               :: fParallel, fPerpendicular, RMinus, RPlus, fSpring, fDNEB

    DO iCLONE = 1, G%Clones+1
      RMSd = Zero
      DO j = 1, G%Clone(iCLONE)%NAtms
        RMSd = RMSd &
             + (G%Clone(iCLONE)%Carts%D(1, j) - G%Clone(iCLONE-1)%Carts%D(1, j))**2 &
             + (G%Clone(iCLONE)%Carts%D(2, j) - G%Clone(iCLONE-1)%Carts%D(2, j))**2 &
             + (G%Clone(iCLONE)%Carts%D(3, j) - G%Clone(iCLONE-1)%Carts%D(3, j))**2
      ENDDO
      RMSd = SQRT(RMSd/G%Clone(iCLONE)%NAtms)
      CALL MondoLog(DEBUG_NONE, "NEBForce", "RMSd("//TRIM(IntToChar(iCLONE-1))// &
        " --> "//TRIM(IntToChar(iCLONE))//") = "//TRIM(DblToChar(RMSd*AUToAngstroms))//" A")
    ENDDO

    DO iCLONE = 1, G%Clones+1
      RMSd = Zero
      DO j = 1, G%Clone(iCLONE)%NAtms
        RMSd = RMSd &
             + (G%Clone(iCLONE)%Carts%D(1, j) - G%Clone(0)%Carts%D(1, j))**2 &
             + (G%Clone(iCLONE)%Carts%D(2, j) - G%Clone(0)%Carts%D(2, j))**2 &
             + (G%Clone(iCLONE)%Carts%D(3, j) - G%Clone(0)%Carts%D(3, j))**2
      ENDDO
      RMSd = SQRT(RMSd/G%Clone(iCLONE)%NAtms)
      CALL MondoLog(DEBUG_NONE, "NEBForce", "RMSd(0 --> "// &
        TRIM(IntToChar(iCLONE))//") = "//TRIM(DblToChar(RMSd*AUToAngstroms))//" A")
    ENDDO

    DO iCLONE = 1, G%Clones
      theta = Zero
      thetaNormMinus = Zero
      thetaNormPlus = Zero
      DO j = 1, G%Clone(iCLONE)%NAtms
        theta = theta &
          + (G%Clone(iCLONE-1)%Carts%D(1, j)-G%Clone(iCLONE)%Carts%D(1, j)) &
          * (G%Clone(iCLONE+1)%Carts%D(1, j)-G%Clone(iCLONE)%Carts%D(1, j)) &
          + (G%Clone(iCLONE-1)%Carts%D(2, j)-G%Clone(iCLONE)%Carts%D(2, j)) &
          * (G%Clone(iCLONE+1)%Carts%D(2, j)-G%Clone(iCLONE)%Carts%D(2, j)) &
          + (G%Clone(iCLONE-1)%Carts%D(3, j)-G%Clone(iCLONE)%Carts%D(3, j)) &
          * (G%Clone(iCLONE+1)%Carts%D(3, j)-G%Clone(iCLONE)%Carts%D(3, j))
        thetaNormMinus = thetaNormMinus &
          + (G%Clone(iCLONE-1)%Carts%D(1, j)-G%Clone(iCLONE)%Carts%D(1, j))**2 &
          + (G%Clone(iCLONE-1)%Carts%D(2, j)-G%Clone(iCLONE)%Carts%D(2, j))**2 &
          + (G%Clone(iCLONE-1)%Carts%D(3, j)-G%Clone(iCLONE)%Carts%D(3, j))**2
        thetaNormPlus = thetaNormPlus &
          + (G%Clone(iCLONE+1)%Carts%D(1, j)-G%Clone(iCLONE)%Carts%D(1, j))**2 &
          + (G%Clone(iCLONE+1)%Carts%D(2, j)-G%Clone(iCLONE)%Carts%D(2, j))**2 &
          + (G%Clone(iCLONE+1)%Carts%D(3, j)-G%Clone(iCLONE)%Carts%D(3, j))**2
      ENDDO
      theta = theta/SQRT(thetaNormMinus*thetaNormPlus)
      theta = MIN(1.0D0, theta)
      theta = MAX(-1.0D0, theta)
      CALL MondoLog(DEBUG_NONE, "NEBForce", "theta ("//TRIM(IntToChar(iCLONE-1))//" - "// &
        TRIM(IntToChar(iCLONE))//" - "//TRIM(IntToChar(iCLONE+1))//") = "// &
        TRIM(FltToShrtChar(ACOS(theta)*RadToDeg))//" degrees")
    ENDDO

    ! Set the tangent for the endpoints.
    N(:,:,0) = G%Clone(1)%Carts%D - G%Clone(0)%Carts%D
    N(:,:,G%Clones+1) = G%Clone(G%Clones+1)%Carts%D - G%Clone(G%Clones)%Carts%D

    ! Find the climbing image (in case we want to use it).
    UMax = G%Clone(1)%ETotal
    UMaxClone = 1

    ! Check for climbing image.
    DO iCLONE = 1, G%Clones
      IF(G%Clone(iCLONE)%ETotal > UMax) THEN
        UMax = G%Clone(iCLONE)%ETotal
        UMaxClone = iCLONE
      ENDIF
    ENDDO

    CALL MondoLog(DEBUG_NONE, "NEBForce", "climbing image is Clone "// &
      TRIM(IntToChar(UMaxClone))//" with energy "// &
      TRIM(DblToChar(G%Clone(UMaxClone)%ETotal*au2eV))//" eV")

    ! Loop over the clones.
    DO iCLONE = 1, G%Clones

      ! Are the neighboring images higher in energy?
      IF(G%Clone(iCLONE-1)%ETotal > G%Clone(iCLONE)%ETotal) THEN
        UPMinus = .TRUE.
      ELSE
        UPMinus = .FALSE.
      ENDIF

      IF(G%Clone(iCLONE+1)%ETotal > G%Clone(iCLONE)%ETotal) THEN
        UPPlus = .TRUE.
      ELSE
        UPPlus = .FALSE.
      ENDIF

      ! Get tangent.
      IF(UPMinus .NEQV. UPPlus) THEN
        ! If we are not at an extremum of energy, the tangent is the vector to
        ! the lower energy neighbour
        IF(UPMinus)THEN
          N(:,:,iCLONE) = G%Clone(iCLONE)%Carts%D - G%Clone(iCLONE-1)%Carts%D
        ELSE
          N(:,:,iCLONE) = G%Clone(iCLONE+1)%Carts%D - G%Clone(iCLONE)%Carts%D
        ENDIF
      ELSE
        ! At an extremum of energy, interpolate the tangent linearly with the
        ! energy. This is done so that the direction of the tangent does not
        ! change abruptly, improving the stability of the NEB method and
        ! avoiding artificial kinks in high force regions along the path. See
        ! discussion in Secs. 3 and 4 of JCP 113, 9978.
        CALL MondoLog(DEBUG_NONE, "NEBForce", "UPMinus == UPPlus", "Clone "//TRIM(IntToChar(iCLONE)))

        UMinus = G%Clone(iCLONE-1)%ETotal - G%Clone(iCLONE)%ETotal
        UPlus  = G%Clone(iCLONE+1)%ETotal - G%Clone(iCLONE)%ETotal

        UMin = MIN(ABS(UPlus),ABS(UMinus))
        UMax = MAX(ABS(UPlus),ABS(UMinus))

        CALL MondoLog(DEBUG_NONE, "NEBForce", "UMin = "//TRIM(DblToChar(UMin*au2eV))//" eV", "Clone "//TRIM(IntToChar(iCLONE)))
        CALL MondoLog(DEBUG_NONE, "NEBForce", "UMax = "//TRIM(DblToChar(UMax*au2eV))//" eV", "Clone "//TRIM(IntToChar(iCLONE)))

        IF(UMinus > UPlus) THEN
          N(:,:,iCLONE) = (G%Clone(iCLONE+1)%Carts%D - G%Clone(iCLONE)%Carts%D)*UMin &
                        + (G%Clone(iCLONE)%Carts%D - G%Clone(iCLONE-1)%Carts%D)*UMax
        ELSE
          N(:,:,iCLONE) = (G%Clone(iCLONE+1)%Carts%D - G%Clone(iCLONE)%Carts%D)*UMax &
                        + (G%Clone(iCLONE)%Carts%D - G%Clone(iCLONE-1)%Carts%D)*UMin
        ENDIF
      ENDIF
    ENDDO

    DO iCLONE = 0, G%Clones+1
      CALL MondoLog(DEBUG_NONE, "NEBForce", "Coordinates (in A)", "Clone "//TRIM(IntToChar(iCLONE)))
      DO j = 1, G%Clone(iCLONE)%NAtms
        CALL MondoLog(DEBUG_NONE, "NEBForce", "R["//TRIM(IntToChar(j))//"] = "// &
          TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(1,j)*AUToAngstroms))//" "// &
          TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(2,j)*AUToAngstroms))//" "// &
          TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(3,j)*AUToAngstroms)), &
          "Clone "//TRIM(IntToChar(iCLONE)))
      ENDDO
    ENDDO

    DO iCLONE = 0, G%Clones+1
      CALL MondoLog(DEBUG_NONE, "NEBForce", "Spring tangent (in A)", "Clone "//TRIM(IntToChar(iCLONE)))
      DO j = 1, G%Clone(iCLONE)%NAtms
        CALL MondoLog(DEBUG_NONE, "NEBForce", "N["//TRIM(IntToChar(j))//"] = "// &
          TRIM(DblToChar(N(1,j,iCLONE)*AUToAngstroms))//" "// &
          TRIM(DblToChar(N(2,j,iCLONE)*AUToAngstroms))//" "// &
          TRIM(DblToChar(N(3,j,iCLONE)*AUToAngstroms)), &
          "Clone "//TRIM(IntToChar(iCLONE)))
      ENDDO
    ENDDO

    ! Normalize tangents.
    DO iCLONE = 0, G%Clones+1
      magnitude = Zero
      DO j = 1, G%Clone(iCLONE)%NAtms
        magnitude = magnitude + N(1, j, iCLONE)**2 + N(2, j, iCLONE)**2 + N(3, j, iCLONE)**2
      ENDDO
      IF(magnitude > Zero) THEN
        DO j = 1, G%Clone(iCLONE)%NAtms
          N(:, j, iCLONE) = N(:, j, iCLONE)/SQRT(magnitude)
        ENDDO
      ENDIF
    ENDDO

    DO iCLONE = 0, G%Clones+1
      CALL MondoLog(DEBUG_NONE, "NEBForce", "Spring tangent, normalized", "Clone "//TRIM(IntToChar(iCLONE)))
      DO j = 1, G%Clone(iCLONE)%NAtms
        CALL MondoLog(DEBUG_NONE, "NEBForce", "N["//TRIM(IntToChar(j))//"] = "// &
          TRIM(DblToChar(N(1,j,iCLONE)))//" "// &
          TRIM(DblToChar(N(2,j,iCLONE)))//" "// &
          TRIM(DblToChar(N(3,j,iCLONE))), &
          "Clone "//TRIM(IntToChar(iCLONE)))
      ENDDO
    ENDDO

    DO iCLONE = 1, G%Clones
      CALL MondoLog(DEBUG_NONE, "NEBForce", "Force without NEB (in eV/A)", "Clone "//TRIM(IntToChar(iCLONE)))
      DO j = 1, G%Clone(iCLONE)%NAtms
        CALL MondoLog(DEBUG_NONE, "NEBForce", "F["//TRIM(IntToChar(j))//"] = "// &
          TRIM(DblToChar(-G%Clone(iCLONE)%Gradients%D(1,j)*au2eV/AUToAngstroms))//" "// &
          TRIM(DblToChar(-G%Clone(iCLONE)%Gradients%D(2,j)*au2eV/AUToAngstroms))//" "// &
          TRIM(DblToChar(-G%Clone(iCLONE)%Gradients%D(3,j)*au2eV/AUToAngstroms)), &
          "Clone "//TRIM(IntToChar(iCLONE)))
      ENDDO
    ENDDO

    DO iCLONE = 1, G%Clones
      ! Project out the force along the tangent.
      magnitude = Zero
      DO j = 1, G%Clone(iCLONE)%NAtms
        magnitude = magnitude &
                  - G%Clone(iCLONE)%Gradients%D(1, j)*N(1, j, iCLONE) &
                  - G%Clone(iCLONE)%Gradients%D(2, j)*N(2, j, iCLONE) &
                  - G%Clone(iCLONE)%Gradients%D(3, j)*N(3, j, iCLONE)
      ENDDO

      DO j = 1, G%Clone(iCLONE)%NAtms
        fParallel(1, j) = magnitude * N(1, j, iCLONE)
        fParallel(2, j) = magnitude * N(2, j, iCLONE)
        fParallel(3, j) = magnitude * N(3, j, iCLONE)

        fPerpendicular(1, j) = -G%Clone(iCLONE)%Gradients%D(1, j) - fParallel(1, j)
        fPerpendicular(2, j) = -G%Clone(iCLONE)%Gradients%D(2, j) - fParallel(2, j)
        fPerpendicular(3, j) = -G%Clone(iCLONE)%Gradients%D(3, j) - fParallel(3, j)
      ENDDO

      CALL MondoLog(DEBUG_NONE, "NEBForce", "parallel force component (in eV/A)", "Clone "//TRIM(IntToChar(iCLONE)))
      DO j = 1, G%Clone(iCLONE)%NAtms
        CALL MondoLog(DEBUG_NONE, "NEBForce", "F["//TRIM(IntToChar(j))//"] = "// &
          TRIM(DblToChar(fParallel(1,j)*au2eV/AUToAngstroms))//" "// &
          TRIM(DblToChar(fParallel(2,j)*au2eV/AUToAngstroms))//" "// &
          TRIM(DblToChar(fParallel(3,j)*au2eV/AUToAngstroms)), &
          "Clone "//TRIM(IntToChar(iCLONE)))
      ENDDO

      CALL MondoLog(DEBUG_NONE, "NEBForce", "perpendicular force component (in eV/A)", "Clone "//TRIM(IntToChar(iCLONE)))
      DO j = 1, G%Clone(iCLONE)%NAtms
        CALL MondoLog(DEBUG_NONE, "NEBForce", "F["//TRIM(IntToChar(j))//"] = "// &
          TRIM(DblToChar(fPerpendicular(1,j)*au2eV/AUToAngstroms))//" "// &
          TRIM(DblToChar(fPerpendicular(2,j)*au2eV/AUToAngstroms))//" "// &
          TRIM(DblToChar(fPerpendicular(3,j)*au2eV/AUToAngstroms)), &
          "Clone "//TRIM(IntToChar(iCLONE)))
      ENDDO

      IF(O%NEBClimb .AND. iCLONE == UMaxClone)THEN
        ! If this is the climbing image, the force along the tangent is
        ! inverted.
        CALL MondoLog(DEBUG_NONE, "NEBForce", "no spring force for climbing image", "Clone "//TRIM(IntToChar(iCLONE)))
        G%Clone(iCLONE)%Gradients%D = G%Clone(iCLONE)%Gradients%D + 2.0*fParallel
      ELSE
        ! The non-climbing image.
        !
        ! Calculate the spring force.
        RMinus = G%Clone(iCLONE-1)%Carts%D - G%Clone(iCLONE)%Carts%D
        RPlus  = G%Clone(iCLONE+1)%Carts%D - G%Clone(iCLONE)%Carts%D

        DO j = 1, G%Clone(iCLONE)%NAtms
          RMinusMagnitude = SQRT(RMinus(1, j)**2 + RMinus(2, j)**2 + RMinus(3, j)**2)
          RPlusMagnitude  = SQRT(RPlus(1, j)**2  + RPlus(2, j)**2  + RPlus(3, j)**2)

          fSpring(1, j) = (RPlusMagnitude - RMinusMagnitude)*O%NEBSpring*N(1, j, iCLONE)
          fSpring(2, j) = (RPlusMagnitude - RMinusMagnitude)*O%NEBSpring*N(2, j, iCLONE)
          fSpring(3, j) = (RPlusMagnitude - RMinusMagnitude)*O%NEBSpring*N(3, j, iCLONE)

          ! Double-nudging?
          IF(O%NEBDoubleNudge) THEN
            ! Do double nudging...
            !
            ! JCP 120, 2082 (2004).
          ELSE
            fDNEB(:, j) = Zero
          ENDIF

        ENDDO

        CALL MondoLog(DEBUG_NONE, "NEBForce", "spring force (in eV/A)", "Clone "//TRIM(IntToChar(iCLONE)))
        DO j = 1, G%Clone(iCLONE)%NAtms
          CALL MondoLog(DEBUG_NONE, "NEBForce", "F["//TRIM(IntToChar(j))//"] = "// &
            TRIM(DblToChar(fSpring(1,j)*au2eV/AUToAngstroms))//" "// &
            TRIM(DblToChar(fSpring(2,j)*au2eV/AUToAngstroms))//" "// &
            TRIM(DblToChar(fSpring(3,j)*au2eV/AUToAngstroms)), &
            "Clone "//TRIM(IntToChar(iCLONE)))
        ENDDO

        CALL MondoLog(DEBUG_NONE, "NEBForce", "DNEB spring force (in eV/A)", "Clone "//TRIM(IntToChar(iCLONE)))
        DO j = 1, G%Clone(iCLONE)%NAtms
          CALL MondoLog(DEBUG_NONE, "NEBForce", "F["//TRIM(IntToChar(j))//"] = "// &
            TRIM(DblToChar(fDNEB(1,j)*au2eV/AUToAngstroms))//" "// &
            TRIM(DblToChar(fDNEB(2,j)*au2eV/AUToAngstroms))//" "// &
            TRIM(DblToChar(fDNEB(3,j)*au2eV/AUToAngstroms)), &
            "Clone "//TRIM(IntToChar(iCLONE)))
        ENDDO

        ! Calculate the force magnitudes and compare.
        DO j = 1, G%Clone(iCLONE)%NAtms
          fMagnitude = SQRT(G%Clone(iCLONE)%Gradients%D(1, j)**2 &
                          + G%Clone(iCLONE)%Gradients%D(2, j)**2 &
                          + G%Clone(iCLONE)%Gradients%D(3, j)**2)
          fSpringMagnitude = SQRT(fSpring(1, j)**2 + fSpring(2, j)**2 + fSpring(3, j)**2)
          IF(fMagnitude > Zero) THEN
            IF(fSpringMagnitude/fMagnitude > 1.0D-3) THEN
              CALL MondoLog(DEBUG_NONE, "NEBForce", "fSpring/f["//TRIM(IntToChar(j))// &
                "] = "//TRIM(DblToChar(fSpringMagnitude/fMagnitude))// &
                " --WARNING-- spring force is more than 1/10 % of force", "Clone "//TRIM(IntToChar(iCLONE)))
            ELSE
              CALL MondoLog(DEBUG_NONE, "NEBForce", "fSpring/f["//TRIM(IntToChar(j))// &
                "] = "//TRIM(DblToChar(fSpringMagnitude/fMagnitude)), "Clone "//TRIM(IntToChar(iCLONE)))
            ENDIF
          ELSE
            CALL MondoLog(DEBUG_NONE, "NEBForce", "fSpring/f["//TRIM(IntToChar(j))// &
              "] = N/A (f = 0)", "Clone "//TRIM(IntToChar(iCLONE)))
          ENDIF

        ENDDO

        G%Clone(iCLONE)%Gradients%D = -(fPerpendicular + fSpring)
      ENDIF
    ENDDO

    DO iCLONE = 1, G%Clones
      fMax = Zero
      DO j = 1, G%Clone(iCLONE)%NAtms
        IF(ABS(G%Clone(iCLONE)%Gradients%D(1, j)) > fMax) fMax = ABS(G%Clone(iCLONE)%Gradients%D(1, j))
        IF(ABS(G%Clone(iCLONE)%Gradients%D(2, j)) > fMax) fMax = ABS(G%Clone(iCLONE)%Gradients%D(2, j))
        IF(ABS(G%Clone(iCLONE)%Gradients%D(3, j)) > fMax) fMax = ABS(G%Clone(iCLONE)%Gradients%D(3, j))
      ENDDO
      CALL MondoLog(DEBUG_NONE, "NEBForce", "Force with NEB, fMax = "// &
        TRIM(DblToChar(fMax*au2eV/AUToAngstroms))//" (in eV/A)", "Clone "//TRIM(IntToChar(iCLONE)))
      DO j = 1, G%Clone(iCLONE)%NAtms
        CALL MondoLog(DEBUG_NONE, "NEBForce", "F["//TRIM(IntToChar(j))//"] = "// &
          TRIM(DblToChar(-G%Clone(iCLONE)%Gradients%D(1,j)*au2eV/AUToAngstroms))//" "// &
          TRIM(DblToChar(-G%Clone(iCLONE)%Gradients%D(2,j)*au2eV/AUToAngstroms))//" "// &
          TRIM(DblToChar(-G%Clone(iCLONE)%Gradients%D(3,j)*au2eV/AUToAngstroms)), &
          "Clone "//TRIM(IntToChar(iCLONE)))
      ENDDO
    ENDDO

  END SUBROUTINE NEBForce

  !===============================================================================
  ! Generate a cubic spline along the band and interpolate to find extrema
  !===============================================================================
  SUBROUTINE NEBSpline()
    !----------------------------------------------------------------------------
    ! Not implemented yet
  END SUBROUTINE NEBSpline

END MODULE NEB
