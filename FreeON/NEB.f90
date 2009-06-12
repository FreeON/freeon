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

    REAL(DOUBLE), DIMENSION(3, 6, 10), PARAMETER :: HardcodedClone = RESHAPE( (/ &

    -0.0387366602204657,     0.0000393965267929,     0.0410905670798839, &
     0.9530736777049564,     0.0002525697906025,     0.1445955034255581, &
    -0.0752357325259799,    -0.0001496215393734,    -0.9540061478025289, &
     2.6556519212575451,     0.0000281399113110,    -0.0500862325818810, &
     3.0423878849392021,    -0.7590283573696242,     0.4638287628034553, &
     3.0424989088397254,     0.7590187773407788,     0.4638394736799561, &

    -0.0794497850847175,     0.0000277437598462,     0.0620395952368433, &
     0.8888604802485919,     0.0002721687111570,     0.2955608810390860, &
     0.0282124933919739,    -0.0002032686396689,    -0.9280110154374481, &
     2.6040657379359171,     0.0000266791727028,    -0.0791644284574151, &
     3.0689305423806683,    -0.7593321897019759,     0.3646128932842665, &
     3.0690205311195999,     0.7593261701816222,     0.3646167256429593, &

    -0.1178522556007471,     0.0000162244709894,     0.0624815202079251, &
     0.8006330625464896,     0.0002947321882481,     0.4454547972330142, &
     0.1516223428270747,    -0.0002501046220511,    -0.8962243886390246, &
     2.5554102029040431,     0.0000215466708141,    -0.0857827524034297, &
     3.0948807331078116,    -0.7592959963022474,     0.2637122651323487, &
     3.0949459142136417,     0.7592921654072674,     0.2637098448399127, &

    -0.1481236723682279,     0.0000068900831750,     0.0447475195245014, &
     0.6866891940568129,     0.0003122315911442,     0.5856039279406026, &
     0.2903285979244061,    -0.0002864754847743,    -0.8491221071550005, &
     2.5151912642506002,     0.0000135050891977,    -0.0667557960058376, &
     3.1177583294099578,    -0.7591810226142175,     0.1576360147321610, &
     3.1177962867298197,     0.7591786668962897,     0.1576297122637085, &

    -0.1638838343236930,     0.0000015606479687,     0.0153925297309695, &
     0.5575079362564702,     0.0003171628820935,     0.7004606590777634, &
     0.4307653487143945,    -0.0003089456644423,    -0.7826403055898776, &
     2.4925857934647260,     0.0000043345424929,    -0.0242317439001282, &
     3.1313265141579647,    -0.7591389334375256,     0.0500325441935843, &
     3.1313382417208215,     0.7591381316627542,     0.0500248703584339, &

    -0.1638838343236524,    -0.0000015606479686,    -0.0153925297308907, &
     0.4307653487147219,     0.0003089456644423,     0.7826403055897587, &
     0.5575079362562541,    -0.0003171628820937,    -0.7004606590779888, &
     2.4925857934651638,    -0.0000043345424970,     0.0242317439002015, &
     3.1313382417220108,    -0.7591381316629530,    -0.0500248703582935, &
     3.1313265141591455,     0.7591389334377275,    -0.0500325441934447, &

    -0.1481236723680223,    -0.0000068900831748,    -0.0447475195245079, &
     0.2903285979244822,     0.0002864754847735,     0.8491221071550337, &
     0.6866891940569854,    -0.0003122315911435,    -0.5856039279405969, &
     2.5151912642508361,    -0.0000135050892006,     0.0667557960063135, &
     3.1177962867297522,    -0.7591786668962823,    -0.1576297122639786, &
     3.1177583294098885,     0.7591810226142078,    -0.1576360147324317, &

    -0.1178522556008283,    -0.0000162244709897,    -0.0624815202079347, &
     0.1516223428272587,     0.0002501046220514,     0.8962243886389735, &
     0.8006330625463327,    -0.0002947321882482,    -0.4454547972332506, &
     2.5554102029044734,    -0.0000215466708139,     0.0857827524038935, &
     3.0949459142140769,    -0.7592921654074238,    -0.2637098448399897, &
     3.0948807331082504,     0.7592959963024053,    -0.2637122651324246, &

    -0.0794497850847173,    -0.0000277437598465,    -0.0620395952368519, &
     0.0282124933918858,     0.0002032686396687,     0.9280110154374427, &
     0.8888604802485930,    -0.0002721687111571,    -0.2955608810390030, &
     2.6040657379360894,    -0.0000266791727059,     0.0791644284572417, &
     3.0690205311194405,    -0.7593261701817332,    -0.3646167256429189, &
     3.0689305423805067,     0.7593321897020844,    -0.3646128932842275, &

    -0.0387366602206245,    -0.0000393965267927,    -0.0410905670797749, &
    -0.0752357325261671,     0.0001496215393729,     0.9540061478026143, &
     0.9530736777048092,    -0.0002525697906020,    -0.1445955034254412, &
     2.6556519212566352,    -0.0000281399113104,     0.0500862325814744, &
     3.0424989088394474,    -0.7590187773408736,    -0.4638394736798984, &
     3.0423878849389268,     0.7590283573697194,    -0.4638287628033949  &

    /), (/ 3, 6, 10 /) )
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
      CALL MondoLog(DEBUG_NONE, "NEBInit", "R["//TRIM(IntToChar(j))//"] = [ "// &
        TRIM(DblToChar(G%Clone(0)%Carts%D(1,j)))//" "// &
        TRIM(DblToChar(G%Clone(0)%Carts%D(2,j)))//" "// &
        TRIM(DblToChar(G%Clone(0)%Carts%D(3,j)))//" ]")
    ENDDO
#endif

    ! Initialize reaction path.
    DO iCLONE=1,G%Clones
       ImageFraction=DBLE(iCLONE)/DBLE(G%Clones+1)
       CALL SetEq(G%Clone(iCLONE),G%Clone(0))

       ! Linear interpolation of path.
       CALL MondoLog(DEBUG_NONE, "NEBInit", "linear interpolation for clone "//TRIM(IntToChar(iCLONE)))
       G%Clone(iCLONE)%Carts%D=G%Clone(0)%Carts%D+ImageFraction*ReactionVector

       ! Hardcoded path.
       !CALL MondoLog(DEBUG_NONE, "NEBInit", "hardcoded configuration for clone "//TRIM(IntToChar(iCLONE)))
       !DO j = 1, G%Clone(0)%NAtms
       !  G%Clone(iCLONE)%Carts%D(:, j) = HardcodedClone(:, j, iCLONE)*AngstromsToAU
       !ENDDO

       ! Set everything else to 0 in this clone.
       G%Clone(iCLONE)%Velocity%D = Zero
       G%Clone(iCLONE)%Fext%D = Zero
       G%Clone(iCLONE)%Gradients%D = Zero
       G%Clone(iCLONE)%CConstrain%I = 0

#ifdef NEB_DEBUG
       CALL MondoLog(DEBUG_NONE, "NEBInit", "Clone "//TRIM(IntToChar(iCLONE)))
       DO j=1, G%Clone(iCLONE)%NAtms
         CALL MondoLog(DEBUG_NONE, "NEBInit", "R["//TRIM(IntToChar(j))//"] = [ "// &
           TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(1,j)))//" "// &
           TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(2,j)))//" "// &
           TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(3,j)))//" ]")
       ENDDO
#endif
    ENDDO

    iCLONE=G%Clones+1
#ifdef NEB_DEBUG
    CALL MondoLog(DEBUG_NONE, "NEBInit", "Product Clone "//TRIM(IntToChar(iCLONE)))
    DO j=1, G%Clone(iCLONE)%NAtms
      CALL MondoLog(DEBUG_NONE, "NEBInit", "R["//TRIM(IntToChar(j))//"] = [ "// &
        TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(1,j)))//" "// &
        TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(2,j)))//" "// &
        TRIM(DblToChar(G%Clone(iCLONE)%Carts%D(3,j)))//" ]")
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
          ! If we are not at an extrema of energy,
          ! the tangent is the vector to the lower energy neighbour
          IF(UPm)THEN
             N=G%Clone(I)%Carts%D-G%Clone(I-1)%Carts%D
          ELSE
             N=G%Clone(I+1)%Carts%D-G%Clone(I)%Carts%D
          ENDIF
       ELSE
          ! At an extrema of energy,
          ! interpolate the tangent linearly with the energy
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

       !GH       write(*,*)'Crds'
       !GH       write(*,'(3F13.5)') (G%Clone(I)%Carts%D(:,j),j=1,G%Clone(0)%NAtms)
       !GH       write(*,*)'Normal'
       !GH       write(*,'(3F13.5)') (N(:,j),j=1,G%Clone(0)%NAtms)
       !GH       write(*,*)'In Force'
       !GH       write(*,'(3F13.5)') (G%Clone(I)%Gradients%D(:,j),j=1,G%Clone(0)%NAtms)

       ! Project out the force along the tangent, unless this is the
       ! climbing image, for which the force along the tangent is inverted
       FProj=SUM(G%Clone(I)%Gradients%D*N)
       !GH       write(*,*)'Prj Force'
       IF(O%NEBClimb.AND.I==UMaxI)THEN
          !GH          write(*,*)'Climbing Image',I
          G%Clone(I)%Gradients%D=G%Clone(I)%Gradients%D-2.0*N*FProj
       ELSE
          G%Clone(I)%Gradients%D=G%Clone(I)%Gradients%D-N*FProj
       ENDIF
       !GH       write(*,'(3F13.5)') (G%Clone(I)%Gradients%D(:,j),j=1,G%Clone(0)%NAtms)

       ! Add spring forces along the band (if we are not the climbing image)
       Rm=SQRT(SUM((G%Clone(I)%Carts%D-G%Clone(I-1)%Carts%D)**2))
       Rp=SQRT(SUM((G%Clone(I)%Carts%D-G%Clone(I+1)%Carts%D)**2))
       !GH       write(*,*)'Rm,Rp',Rm,Rp
       IF(O%NEBClimb.AND.I==UMaxI)THEN
          ! Do nothing (no springs for the climbing image)
       ELSE
          G%Clone(I)%Gradients%D=G%Clone(I)%Gradients%D+O%NEBSpring*N*(Rp-Rm)
       ENDIF

       !GH       write(*,*)'Out Force'
       !GH       write(*,'(3F13.5)') (G%Clone(I)%Gradients%D(:,j),j=1,G%Clone(0)%NAtms)
       ! Write distance, energies and forces
       Dist=Dist+Rm

       CALL MondoLog(DEBUG_NONE, "NEBForce", &
                "Rm = "//TRIM(FltToShrtChar(Rm)) &
            //", Rp = "//TRIM(FltToShrtChar(Rp)) &
            //", Dist = "//TRIM(FltToShrtChar(Dist)) &
            //', E = '//TRIM(DblToMedmChar(G%Clone(I)%ETotal)) &
            //', F = '//TRIM(DblToMedmChar(FProj)), "Image "//TRIM(IntToChar(I)))
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
