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

    -0.0379684684386368,     0.0000938157183353,     0.0843055783659356, &
     0.9593586619845134,     0.0001247056484773,     0.1103633564968881, &
    -0.1593766828191444,    -0.0001856828884250,    -0.9040954661159468, &
     2.6639259147738930,     0.0000513424528050,    -0.0647317438493449, &
     3.0768002097181673,    -0.7589186205445418,     0.4286750412897766, &
     3.0769003646942430,     0.7589569742538715,     0.4286900897487111, &

    -0.0937559175022203,     0.0000609449922573,     0.1178862604430854, &
     0.8916721664426711,     0.0001472422464066,     0.2670108484107643, &
    -0.0814372110122828,    -0.0002581813358750,    -0.8780429712761841, &
     2.6036472437274876,     0.0000620595554307,    -0.0664420536060491, &
     3.1297000296263628,    -0.7595903101579257,     0.3018980507980705, &
     3.1298136886589778,     0.7596433241210134,     0.3018818912993235, &

    -0.1439871549253217,     0.0000410287472592,     0.1074232619168511, &
     0.8003252538353912,     0.0001792072835162,     0.4226949938813975, &
     0.0465622766603119,    -0.0002907103431003,    -0.8701825411676509, &
     2.5568653608024863,     0.0000585902664696,    -0.0262256429856887, &
     3.1598842918552967,    -0.7594316037271355,     0.1953088798814689, &
     3.1599899717478928,     0.7594792395361162,     0.1952581981232711, &

    -0.1823683779424272,     0.0000247379621522,     0.0713148726816970, &
     0.6762675405830877,     0.0002196667221771,     0.5739130487045684, &
     0.2086768238259256,    -0.0002986118687167,    -0.8444561724936545, &
     2.5340523237421557,     0.0000411061868110,     0.0073397062150339, &
     3.1714698694097723,    -0.7582793475221478,     0.1019525250280853, &
     3.1715418203748329,     0.7583100337399371,     0.1018772217681556, &

    -0.2036571338431797,     0.0000073551072755,     0.0217660556865676, &
     0.5210933525337463,     0.0002605883185721,     0.7035640527560384, &
     0.3818155070835987,    -0.0002855583101745,    -0.7832313639567793, &
     2.5296377881723515,     0.0000130347992749,     0.0064309229381284, &
     3.1753640648476402,    -0.7573127371041237,     0.0273469725005775, &
     3.1753864211314609,     0.7573219385765803,     0.0272615021535959, &

    -0.2036571338431080,    -0.0000073551072719,    -0.0217660556860814, &
     0.3818155070855683,     0.0002855583101707,     0.7832313639559373, &
     0.5210933525325109,    -0.0002605883185758,    -0.7035640527570181, &
     2.5296377881702385,    -0.0000130347992812,    -0.0064309229374554, &
     3.1753864211292484,    -0.7573219385766251,    -0.0272615021530600, &
     3.1753640648454016,     0.7573127371041750,    -0.0273469725000388, &

    -0.1823683779404627,    -0.0000247379621513,    -0.0713148726819314, &
     0.2086768238269042,     0.0002986118687151,     0.8444561724938352, &
     0.6762675405856946,    -0.0002196667221792,    -0.5739130487037164, &
     2.5340523237433081,    -0.0000411061868139,    -0.0073397062147851, &
     3.1715418203751091,    -0.7583100337398314,    -0.1018772217686530, &
     3.1714698694100347,     0.7582793475220512,    -0.1019525250285765, &

    -0.1439871549242878,    -0.0000410287472598,    -0.1074232619166213, &
     0.0465622766623916,     0.0002907103431005,     0.8701825411677275, &
     0.8003252538362471,    -0.0001792072835173,    -0.4226949938820085, &
     2.5568653608064462,    -0.0000585902664731,     0.0262256429861633, &
     3.1599899717503344,    -0.7594792395361484,    -0.1952581981231633, &
     3.1598842918577281,     0.7594316037271713,    -0.1953088798813620, &

    -0.0937559174999328,    -0.0000609449922584,    -0.1178862604432407, &
    -0.0814372110100596,     0.0002581813358763,     0.8780429712760780, &
     0.8916721664450168,    -0.0001472422464081,    -0.2670108484107855, &
     2.6036472437300442,    -0.0000620595554314,     0.0664420536065462, &
     3.1298136886618915,    -0.7596433241209953,    -0.3018818912994259, &
     3.1297000296292712,     0.7595903101579113,    -0.3018980507981722, &

    -0.0379684684412321,    -0.0000938157183345,    -0.0843055783654090, &
    -0.1593766828201745,     0.0001856828884253,     0.9040954661166546, &
     0.9593586619818703,    -0.0001247056484792,    -0.1103633564979257, &
     2.6639259147782028,    -0.0000513424528038,     0.0647317438516321, &
     3.0769003646935897,    -0.7589569742540995,    -0.4286900897497586, &
     3.0768002097175198,     0.7589186205447672,    -0.4286750412908237  &

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
       CALL MondoLog(DEBUG_NONE, "NEBInit", "linear interpolation for clone "//TRIM(IntToChar(iCLONE)))
       G%Clone(iCLONE)%Carts%D=G%Clone(0)%Carts%D+ImageFraction*ReactionVector

       ! Hardcoded path.
       !CALL MondoLog(DEBUG_NONE, "NEBInit", "hardcoded configuration for clone "//TRIM(IntToChar(iCLONE)))
       !DO j = 1, G%Clone(0)%NAtms
       !  G%Clone(iCLONE)%Carts%D(:, j) = HardcodedClone(:, j, iCLONE)
       !ENDDO

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
