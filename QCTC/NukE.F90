!    COMPUTE THE NUCLEAR-TOTAL ELECTROSTATIC ENERGY IN O(N Lg N) CPU TIME
!    Authors: Matt Challacombe and C.J. Tymczak
!==============================================================================
MODULE NuklarE
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE QCTCThresholds
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
       REAL(DOUBLE)                    :: NukE,NukeCo,NukePole
       REAL(DOUBLE),DIMENSION(1:1)     :: HGBra
       REAL(DOUBLE),DIMENSION(0:0)     :: SPBraC,SPBraS
#ifdef PERIODIC
       INTEGER                         :: NC
       REAL(DOUBLE),DIMENSION(3)       :: QC
#endif
!---------------------------------------------------------------------------------------------
       NukE=Zero 
       DO At=1,NAtoms
          NukeCo  = -GM%AtNum%I(At)*(NuclearExpnt/Pi)**(ThreeHalves)
          NukePole= -GM%AtNum%I(At)
          DP2=(GM%AtNum%I(At)/TauMAC)**(Two/DBLE(SPEll+2))
!         Initialize the  Accumulators
          HGBra(1)  = NukeCo
          HGKet(1)  = Zero
          SPBraC(0) = NukePole
          SPBraS(0) = Zero
          SPKetC(0) = Zero
          SPKetS(0) = Zero 
!         Set atomic "primitive" 
          Prim%Ell  = 0
          Prim%P    = GM%Carts%D(:,At)
          Prim%Zeta = NuclearExpnt
          PBox%BndBox(:,1)=Prim%P
          PBox%BndBox(:,2)=Prim%P
          PBox=ExpandBox(PBox,Extent(0,NuclearExpnt,HGBra,TauPAC))
#ifdef PERIODIC
          QC(:)=GM%Carts%D(:,At)
          DO NC=1,CSMM1%NCells
!            Set Atomic Coordinates
             Prim%P(:)  = QC(:)+CSMM1%CellCarts%D(:,NC)
             PBox%Center= Prim%P
!            Walk the walk
             CALL VWalk(PoleRoot)
          ENDDO
!         Reset the Atomic Coordinates
          Prim%P(:)  = QC(:)
          PBox%Center= Prim%P
!         Accumulate the atomic contribution
          NukE=NukE+HGBra(1)*HGKet(1)+SPBraC(0)*SPKetC(0)
!         Add in the Far Field, Dipole and Quadripole  Correction
          IF(Dimen > 0) THEN
             NukE = NukE + CTraxFF(Prim,HGBra)
          ENDIF
#else
!         Walk the walk
          CALL VWalk(PoleRoot)
!         Accumulate the atomic contribution
          NukE=NukE+HGBra(1)*HGKet(1)+SPBraC(0)*SPKetC(0) 
#endif
       ENDDO
     END FUNCTION NukE
!
END MODULE
