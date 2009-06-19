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

    -0.0225464956235612,    -0.0004285713489616,     0.0387454186573339, &
     0.9751275920846620,    -0.0042960384492094,     0.0506019656445768, &
    -0.1621388796463510,     0.0095236342903110,    -0.9471453669480014, &
     2.6938736411393824,    -0.0072318358972327,    -0.0334296623125359, &
     3.0498388619316317,    -0.7575799491160704,     0.5143680388569060, &
     3.0454854227338406,     0.7601845653217942,     0.4935323175486023, &

    -0.1305492259696334,    -0.0082857636535119,     0.0265697116371085, &
     0.8283028711231523,    -0.0292584928958058,     0.2954193765789251, &
     0.0064855967955351,     0.0329530039876786,    -0.9591458276080755, &
     2.5795878063847093,    -0.0058985801234155,     0.0438999292310092, &
     3.1538846728040308,    -0.7542490247663468,     0.3592741207384341, &
     3.1419281383597446,     0.7648725134167307,     0.3247281813295995, &

    -0.1736336720590211,    -0.0078520708498591,     0.0210405054637412, &
     0.7124643721444700,    -0.0247696606808161,     0.4730484787551942, &
     0.1623204675975285,     0.0382669721653310,    -0.9165477049659552, &
     2.5410820825066538,    -0.0137428757477686,     0.0765028837509381, &
     3.1859892099654030,    -0.7560859956125755,     0.2272376462008677, &
     3.1514177048163248,     0.7642791112548130,     0.1835363383651085, &

    -0.1461080880847425,    -0.0028333644132694,     0.0134292927649176, &
     0.6408738335472759,     0.0011940624274350,     0.6198075606252413, &
     0.3754369591777665,     0.0297334082200154,    -0.8352112791223436, &
     2.5325437693113266,    -0.0282364098954122,     0.0875501380912344, &
     3.1315803819179573,    -0.8142504309060065,     0.0895263280182630, &
     3.0453131898924779,     0.8144499960044141,     0.0637888764768322, &

    -0.1473028621193182,     0.0322545410301132,    -0.0198941826005811, &
     0.4803718318184713,     0.0039026802470152,     0.7526539905342892, &
     0.5418594343384163,    -0.0080572787675319,    -0.7374875182427375, &
     2.5423914040085274,    -0.0282762390378007,     0.0310524430331555, &
     3.1206025758993534,    -0.8255101980492303,    -0.0133851199959591, &
     3.0417175202555828,     0.8257055550252480,     0.0000240429094105, &

    -0.1419630214020312,    -0.0094269906069960,     0.0008680941462734, &
     0.5076906512334726,    -0.0086758217502434,     0.7557511515971845, &
     0.5276596232824522,     0.0165263753546337,    -0.7360072065392828, &
     2.5335147818358124,     0.0009233415058229,     0.0348690700492715, &
     3.0765399458916276,    -0.8236778295338349,    -0.0249169566452696, &
     3.0761978116009856,     0.8243118329178392,    -0.0435277812956071, &

    -0.1339390453937901,     0.0129721499871359,    -0.0145904659618596, &
     0.3999647298812341,    -0.0190824064794037,     0.8267760618979123, &
     0.6411350857919261,     0.0060282343684871,    -0.6347691041207585, &
     2.5226678909135596,     0.0022274975188881,    -0.0262084688345510, &
     3.0696635998558346,    -0.8165173399815321,    -0.0871880151669101, &
     3.0801476431522699,     0.8143145873167634,    -0.1029109339812010, &

    -0.1799796184726618,     0.0166641121108839,    -0.0186982697704393, &
     0.1595700834674049,    -0.0092009669993689,     0.9187733795504495, &
     0.7032843437479059,     0.0109462363073876,    -0.4771406801460116, &
     2.5426476559944184,    -0.0152149015293707,    -0.0924593416657353, &
     3.1919881365427090,    -0.7604268303422111,    -0.2050864563875936, &
     3.1621295338892579,     0.7571369065943794,    -0.1902068313046307, &

    -0.1410595033490503,    -0.0045253414975949,    -0.0191928014463507, &
     0.0180928200672215,    -0.0030499514254417,     0.9640745513356053, &
     0.8114115816485612,     0.0048240411196737,    -0.3101856598790826, &
     2.5719087704922985,     0.0061217130811837,    -0.0706893326492596, &
     3.1461136104002869,    -0.7616789115083143,    -0.3354275989087419, &
     3.1731726100405520,     0.7581748147542637,    -0.3193245926171713, &

    -0.0504108072399734,    -0.0022881436498385,    -0.0267860718205956, &
    -0.1531015288650077,    -0.0063306916302273,     0.9637317725843537, &
     0.9459759316799826,     0.0083399030917806,    -0.0733605702538127, &
     2.6620925497342980,     0.0031817943526740,    -0.0049307232880860, &
     3.0840427677407054,    -0.7602395478294008,    -0.4832686375731840, &
     3.0910409054693595,     0.7571648603752852,    -0.4920584801642321  &

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
            //', F = '//TRIM(DblToMedmChar(FProj)), "Clone "//TRIM(IntToChar(I)))
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
