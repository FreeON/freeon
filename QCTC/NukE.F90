!==============================================================================
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe and C.J. Tymczak
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!--  COMPUTE THE NUCLEAR-TOTAL ELECTROSTATIC ENERGY IN O(N Lg N) CPU TIME
!==============================================================================
MODULE NuklarE
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE Thresholding
  USE Globals
  USE AtomPairs
  USE BraBloks
  USE PoleTree 
  USE TreeWalk
#ifdef PERIODIC
  USE PBCFarField
#endif
  IMPLICIT NONE
!----------!
  CONTAINS !
!=============================================================================================
!
!=============================================================================================
    FUNCTION NukE()
       REAL(DOUBLE)                   :: NukE,Tau,NukeCo,NukePole
       REAL(DOUBLE),DIMENSION(1:1)    :: HGBra
       REAL(DOUBLE),DIMENSION(0:0)    :: SPBraC,SPBraS
#ifdef PERIODIC
       INTEGER                        :: NC
       REAL(DOUBLE),DIMENSION(3)      :: QC
#endif
!---------------------------------------------------------------------------------------------
       NukE=Zero 
       DO At=1,NAtoms
          NukeCo    = -GM%AtNum%I(At)*(NuclearExpnt/Pi)**(ThreeHalves)
          NukePole  = -GM%AtNum%I(At)
!         Set MAC and PAC parameters
          Tau       = Thresholds%TwoE
          DP2       = (ABS(NukePole)/Tau)**(Two/DBLE(SPEll+1))
          PoleSwitch= MAX(1.D-12,-LOG(Two*Tau/ABS(NukePole)))
!         Set atomic "primitive"
          Prim%Ell  = 0
          Prim%Zeta = NuclearExpnt
!         Zero Accumulators
          HGBra(1)  = NukeCo
          SPBraC(0) = NukePole
          SPBraS(0) = Zero
          HGKet(1)  = Zero
          SPKetC(0) = Zero
          SPKetS(0) = Zero
#ifdef PERIODIC
          QC(:)=GM%Carts%D(:,At)
          DO NC=1,CSMM1%NCells
!            Set Atomic Coordinates
             Prim%P(:)=QC(:)+CSMM1%CellCarts%D(:,NC)
!            Walk the walk
             CALL VWalk(PoleRoot)
          ENDDO
!         Reset the Atomic Coordinates
          Prim%P(:) = QC(:)
!         Accumulate the atomic contribution
          NukE=NukE+HGBra(1)*HGKet(1)+SPBraC(0)*SPKetC(0)
!         Add in the Far Field, Dipole and Quadripole  Correction
          IF(Dimen > 0) THEN
             NukE = NukE + CTraxFF(Prim,HGBra) &
                         + CTraxQ(Prim,HGBra)
          ENDIF
#else
!         Set Atomic Coordinates
          Prim%P=GM%Carts%D(:,At) 
!         Walk the walk
          CALL VWalk(PoleRoot)
!         Accumulate the atomic contribution
          NukE=NukE+HGBra(1)*HGKet(1)+SPBraC(0)*SPKetC(0)
#endif
       ENDDO
!
       WRITE(*,*)' NukeE = ',NukE
     END FUNCTION NukE
!
END MODULE
