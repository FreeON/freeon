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
  !  H. Jonsson, G. Mills, and K.W. Jacobsen, "Nudged elastic band method
  !    for finding minimum energy paths of transitions," in Classical and
  !    Quantum Dynamics in Condensed Phase Simulations, edited by B.J.Berne,
  !    G. Ciccotti, and D. F. Coker (World Scientific, Singapore, 1998), p.385
  !  G. Henkelman, B.P. Uberuaga, and H. Jonsson, "A climbing-image NEB method
  !     for finding saddle points and minimum energy paths", J. Chem. Phys,
  !     v113, 9901 (2000).
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
    INTEGER                                     :: iCLONE,j
    CHARACTER(LEN=DEFAULT_CHR_LEN)              :: Message

    REAL(DOUBLE), DIMENSION(3, 6, 2), PARAMETER :: HardcodedClone = RESHAPE( (/ &

          0.76197,       -1.00108,        3.14924, &
          0.75124,       -1.15700,        2.16337, &
          0.84932,       -0.00926,        3.14748, &
          0.75830,        0.00000,        0.64446, &
          0.00000,        0.00000,        0.00000, &
          1.51661,        0.00000,        0.00000, &

          1.02993,        0.55781,        3.28112, &
          1.00122,       -0.43566,        3.18940, &
          0.85522,        0.80533,        2.33270, &
          0.75830,        0.00000,        0.64446, &
          0.00000,        0.00000,        0.00000, &
          1.51661,        0.00000,        0.00000  &

    /), (/ 3, 6, 2 /) )
    !----------------------------------------------------------------------------
    !Initialize each clone to initial state then interpolate Cartesian coordinates
#ifdef NEB_DEBUG
    CALL MondoLog(DEBUG_NONE, "NEBInit", "starting...")
#endif

    ! Calculate reaction path vector.
    ReactionVector=G%Clone(G%Clones+1)%Carts%D-G%Clone(0)%Carts%D

#ifdef NEB_DEBUG
    CALL MondoLog(DEBUG_NONE, "NEBInit", "Reaction vector")
    DO j=1, G%Clone(0)%NAtms
      CALL MondoLog(DEBUG_NONE, "NEBInit", "RV["//TRIM(IntToChar(j))//"] = [ "// &
        TRIM(DblToChar(ReactionVector(1,j)))//" "// &
        TRIM(DblToChar(ReactionVector(2,j)))//" "// &
        TRIM(DblToChar(ReactionVector(3,j)))//" ]")
    ENDDO
#endif

#ifdef NEB_DEBUG
    CALL MondoLog(DEBUG_NONE, "NEBInit", "Reactant Clone "//TRIM(IntToChar(0)))
    DO j=1, G%Clone(0)%NAtms
      CALL MondoLog(DEBUG_NONE, "NEBInit", TRIM(G%Clone(0)%AtNam%C(j))//" "// &
        TRIM(DblToChar(G%Clone(0)%Carts%D(1,j)))//" "// &
        TRIM(DblToChar(G%Clone(0)%Carts%D(2,j)))//" "// &
        TRIM(DblToChar(G%Clone(0)%Carts%D(3,j)))//" "// &
        TRIM(IntToChar(G%Clone(0)%CConstrain%I(j))))
    ENDDO
#endif

    ! Check some things for internal consistency.
    DO j = 1, G%Clone(0)%NAtms
      IF(G%Clone(0)%CConstrain%I(j) /= G%Clone(G%Clones+1)%CConstrain%I(j)) THEN
        CALL MondoLog(DEBUG_NONE, "NEBInit", "constrain on atom "//TRIM(IntToChar(j))//" is different between reactant and product")
        CALL Halt("Constrain mismatch on input")
      ENDIF
    ENDDO

    ! Initialize reaction path.
    DO iCLONE=1,G%Clones
       ImageFraction=DBLE(iCLONE)/DBLE(G%Clones+1)
       CALL SetEq(G%Clone(iCLONE),G%Clone(0))

       ! Fix constrain.
       G%Clone(iCLONE)%CConstrain%I = G%Clone(0)%CConstrain%I

       ! Linear interpolation of path.
       !CALL MondoLog(DEBUG_NONE, "NEBInit", "linear interpolation for clone "//TRIM(IntToChar(iCLONE)))
       !G%Clone(iCLONE)%Carts%D=G%Clone(0)%Carts%D+ImageFraction*ReactionVector

       ! Hardcoded path.
       CALL MondoLog(DEBUG_NONE, "NEBInit", "hardcoded configuration for clone "//TRIM(IntToChar(iCLONE)))
       DO j = 1, G%Clone(0)%NAtms
         G%Clone(iCLONE)%Carts%D(:, j) = HardcodedClone(:, j, iCLONE)
       ENDDO

       ! Set everything else to 0 in this clone.
       G%Clone(iCLONE)%Velocity%D = Zero
       G%Clone(iCLONE)%Fext%D = Zero
       G%Clone(iCLONE)%Gradients%D = Zero

#ifdef NEB_DEBUG
       CALL MondoLog(DEBUG_NONE, "NEBInit", "Clone "//TRIM(IntToChar(iCLONE)))
       DO j=1, G%Clone(iCLONE)%NAtms
         CALL MondoLog(DEBUG_NONE, "NEBInit", TRIM(G%Clone(iCLONE)%AtNam%C(j))//" "// &
           TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(1,j)))//" "// &
           TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(2,j)))//" "// &
           TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(3,j)))//" "// &
           TRIM(IntToChar(G%Clone(iCLONE)%CConstrain%I(j))))
       ENDDO
#endif
    ENDDO

    iCLONE=G%Clones+1
#ifdef NEB_DEBUG
    CALL MondoLog(DEBUG_NONE, "NEBInit", "Product Clone "//TRIM(IntToChar(iCLONE)))
    DO j=1, G%Clone(iCLONE)%NAtms
      CALL MondoLog(DEBUG_NONE, "NEBInit", TRIM(G%Clone(iCLONE)%AtNam%C(j))//" "// &
        TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(1,j)))//" "// &
        TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(2,j)))//" "// &
        TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(3,j)))//" "// &
        TRIM(IntToChar(G%Clone(iCLONE)%CConstrain%I(j))))
    ENDDO
    CALL MondoLog(DEBUG_NONE, "NEBInit", "done NEBInit")
#endif
  END SUBROUTINE NEBInit

  SUBROUTINE NEBPurify(G,Init_O,Print_O)
    TYPE(Geometries)                      :: G
    LOGICAL,OPTIONAL                      :: Init_O,Print_O
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

    CALL MondoLog(DEBUG_NONE, "NEB", "not purifying")
    RETURN

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
    CALL MondoLog(DEBUG_NONE, "NEB", "CLONE 0 AFTER SCALING = ")
    DO I = 1, G%Clone(0)%NAtms
      CALL MondoLog(DEBUG_NONE, "NEB", "R["//TRIM(IntToChar(I))//"] = [ "// &
        TRIM(DblToChar(G%Clone(0)%Carts%D(1,I)*AUToAngstroms))//" "// &
        TRIM(DblToChar(G%Clone(0)%Carts%D(2,I)*AUToAngstroms))//" "// &
        TRIM(DblToChar(G%Clone(0)%Carts%D(3,I)*AUToAngstroms))//" ]", "Clone "//TRIM(IntToChar(0)))
    ENDDO
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
!    IF(PRESENT(Print_O))THEN

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
!   ENDIF

    CALL MondoLog(DEBUG_NONE, "FreeON", Mssg, "NEBPurify("//TRIM(IntToChar(G%Clone(1)%Confg))//')')

  END SUBROUTINE NEBPurify
  !===============================================================================
  ! Project out the force along the band and add spring forces along the band.
  !===============================================================================
  SUBROUTINE NEBForce(G,O)
    TYPE(Geometries) :: G
    TYPE(Options)    :: O
    INTEGER          :: I,j,UMaxI,NAtms
    LOGICAL          :: UPm,UPp
    REAl(DOUBLE)     :: UMin,UMax,Um,Up,Rm,Rp,Dist,FProj
    REAL(DOUBLE),DIMENSION(3,G%Clone(0)%NAtms) :: N
    CHARACTER(LEN=DCL) :: Mssg
    !----------------------------------------------------------------------------
    Dist=0
    ! Find the image with the maximum total energy
    UMax=G%Clone(1)%ETotal
    UMaxI=1
    DO I=2,G%Clones
       IF(G%Clone(I)%ETotal>UMax) THEN
          UMaxI=I
          UMax=G%Clone(I)%ETotal
       ENDIF
    ENDDO
    ! Find the tangent to the path at each image
    ! Project out potential forces along the band
    ! Add spring forces along the band

    !GH    write(*,*)'React Crds'
    !GH    write(*,'(3F13.5)') (G%Clone(0)%Carts%D(:,j),j=1,G%Clone(0)%NAtms)

    !GH    write(*,*)'Prod Crds'
    !GH    write(*,'(3F13.5)') (G%Clone(G%Clones+1)%Carts%D(:,j),j=1,G%Clone(0)%NAtms)

    CALL MondoLog(DEBUG_NONE, "NEBForce", &
          "Dist = "//TRIM(FltToShrtChar(Dist)) &
      //', E = '//TRIM(DblToMedmChar(G%Clone(0)%ETotal)), "Reactant")

    DO I=1,G%Clones
       ! Are the neighboring images higher in energy?
       IF(I==1)THEN
          UPm=.FALSE.
       ELSE
          UPm=G%Clone(I-1)%ETotal>G%Clone(I)%ETotal
       ENDIF
       IF(I==G%Clones)THEN
          UPp=.FALSE.
       ELSE
          UPp=G%Clone(I+1)%ETotal>G%Clone(I)%ETotal
       ENDIF
       IF(UPm.NEQV.UPp)THEN
          ! If we are not at an extremum of energy, the tangent is the vector to
          ! the lower energy neighbour
          IF(UPm)THEN
             N=G%Clone(I)%Carts%D-G%Clone(I-1)%Carts%D
          ELSE
             N=G%Clone(I+1)%Carts%D-G%Clone(I)%Carts%D
          ENDIF
       ELSE
          ! At an extremum of energy, interpolate the tangent linearly with the
          ! energy
          Um=G%Clone(I-1)%ETotal-G%Clone(I)%ETotal
          Up=G%Clone(I+1)%ETotal-G%Clone(I)%ETotal
          UMin=MIN(ABS(Up),ABS(Um))
          UMax=MAX(ABS(Up),ABS(Um))
          IF(Um>Up)THEN
             N=(G%Clone(I+1)%Carts%D-G%Clone(I)%Carts%D)*UMin
             N=N+(G%Clone(I)%Carts%D-G%Clone(I-1)%Carts%D)*UMax
          ELSE
             N=(G%Clone(I+1)%Carts%D-G%Clone(I)%Carts%D)*UMax
             N=N+(G%Clone(I)%Carts%D-G%Clone(I-1)%Carts%D)*UMin
          ENDIF
       ENDIF
       N=N/SQRT(SUM(N**2))

       CALL MondoLog(DEBUG_NONE, "NEBForce", "Coordinates", "Clone "//TRIM(IntToChar(I)))
       DO j = 1, G%Clone(I)%NAtms
        CALL MondoLog(DEBUG_NONE, "NEBForce", "R["//TRIM(IntToChar(j))//"] = "// &
          TRIM(DblToChar(G%Clone(I)%Carts%D(1,j)*AUToAngstroms))//" "// &
          TRIM(DblToChar(G%Clone(I)%Carts%D(2,j)*AUToAngstroms))//" "// &
          TRIM(DblToChar(G%Clone(I)%Carts%D(3,j)*AUToAngstroms))//" A", &
          "Clone "//TRIM(IntToChar(I)))
       ENDDO

       CALL MondoLog(DEBUG_NONE, "NEBForce", "Normal", "Clone "//TRIM(IntToChar(I)))
       DO j = 1, G%Clone(I)%NAtms
        CALL MondoLog(DEBUG_NONE, "NEBForce", "N["//TRIM(IntToChar(j))//"] = "// &
          TRIM(DblToChar(N(1,j)*AUToAngstroms))//" "// &
          TRIM(DblToChar(N(2,j)*AUToAngstroms))//" "// &
          TRIM(DblToChar(N(3,j)*AUToAngstroms))//" A", &
          "Clone "//TRIM(IntToChar(I)))
       ENDDO

       CALL MondoLog(DEBUG_NONE, "NEBForce", "Force in", "Clone "//TRIM(IntToChar(I)))
       DO j = 1, G%Clone(I)%NAtms
        CALL MondoLog(DEBUG_NONE, "NEBForce", "F["//TRIM(IntToChar(j))//"] = "// &
          TRIM(DblToChar(G%Clone(I)%Gradients%D(1,j)))//" "// &
          TRIM(DblToChar(G%Clone(I)%Gradients%D(2,j)))//" "// &
          TRIM(DblToChar(G%Clone(I)%Gradients%D(3,j))), &
          "Clone "//TRIM(IntToChar(I)))
       ENDDO

       ! Project out the force along the tangent, unless this is the climbing
       ! image, for which the force along the tangent is inverted
       FProj=SUM(G%Clone(I)%Gradients%D*N)
       CALL MondoLog(DEBUG_NONE, "NEBForce", "Projected Force = "//TRIM(DblToChar(FProj)), "Clone "//TRIM(IntToChar(I)))
       IF(O%NEBClimb.AND.I==UMaxI)THEN
          CALL MondoLog(DEBUG_NONE, "NEBForce", "Climbing Image "//TRIM(IntToChar(I)))
          G%Clone(I)%Gradients%D=G%Clone(I)%Gradients%D-2.0*N*FProj
       ELSE
          G%Clone(I)%Gradients%D=G%Clone(I)%Gradients%D-N*FProj
       ENDIF
       CALL MondoLog(DEBUG_NONE, "NEBForce", "Force", "Clone "//TRIM(IntToChar(I)))
       DO j = 1, G%Clone(I)%NAtms
        CALL MondoLog(DEBUG_NONE, "NEBForce", "F["//TRIM(IntToChar(j))//"] = "// &
          TRIM(DblToChar(G%Clone(I)%Gradients%D(1,j)))//" "// &
          TRIM(DblToChar(G%Clone(I)%Gradients%D(2,j)))//" "// &
          TRIM(DblToChar(G%Clone(I)%Gradients%D(3,j))), &
          "Clone "//TRIM(IntToChar(I)))
       ENDDO

       ! Add spring forces along the band (if we are not the climbing image)
       Rm=SQRT(SUM((G%Clone(I)%Carts%D-G%Clone(I-1)%Carts%D)**2))
       Rp=SQRT(SUM((G%Clone(I)%Carts%D-G%Clone(I+1)%Carts%D)**2))
       CALL MondoLog(DEBUG_NONE, "NEBForce", "Rm = "//TRIM(DblToChar(Rm))//", Rp = "//TRIM(DblToChar(Rp)))
       IF(O%NEBClimb.AND.I==UMaxI)THEN
          ! Do nothing (no springs for the climbing image)
       ELSE
          G%Clone(I)%Gradients%D=G%Clone(I)%Gradients%D+O%NEBSpring*N*(Rp-Rm)
       ENDIF

       CALL MondoLog(DEBUG_NONE, "NEBForce", "Out Force", "Clone "//TRIM(IntToChar(I)))
       DO j = 1, G%Clone(I)%NAtms
        CALL MondoLog(DEBUG_NONE, "NEBForce", "F["//TRIM(IntToChar(j))//"] = "// &
          TRIM(DblToChar(G%Clone(I)%Gradients%D(1,j)))//" "// &
          TRIM(DblToChar(G%Clone(I)%Gradients%D(2,j)))//" "// &
          TRIM(DblToChar(G%Clone(I)%Gradients%D(3,j))), &
          "Clone "//TRIM(IntToChar(I)))
       ENDDO

       ! Write distance, energies and forces
       Dist=Dist+Rm

       CALL MondoLog(DEBUG_NONE, "NEBForce", &
                "Rm = "//TRIM(FltToShrtChar(Rm)) &
            //", Rp = "//TRIM(FltToShrtChar(Rp)) &
            //", Dist = "//TRIM(FltToShrtChar(Dist)) &
            //', E = '//TRIM(DblToMedmChar(G%Clone(I)%ETotal)) &
            //', F = '//TRIM(DblToMedmChar(FProj)), &
            "Clone "//TRIM(IntToChar(I)))
    ENDDO

    Rm=SQRT(SUM((G%Clone(G%Clones+1)%Carts%D-G%Clone(G%Clones)%Carts%D)**2))
    Dist=Dist+Rm
    CALL MondoLog(DEBUG_NONE, "NEBForce", &
             "Dist = "//TRIM(FltToShrtChar(Dist)) &
         //', E = '//TRIM(DblToMedmChar(G%Clone(G%Clones+1)%ETotal)), "Product")
  END SUBROUTINE NEBForce

  !===============================================================================
  ! Generate a cubic spline along the band and interpolate to find extrema
  !===============================================================================
  SUBROUTINE NEBSpline()
    !----------------------------------------------------------------------------
    ! Not implemented yet
  END SUBROUTINE NEBSpline

END MODULE NEB
