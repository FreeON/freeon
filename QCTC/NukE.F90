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
          NukeCo=-GM%AtNum%I(At)*(NuclearExpnt/Pi)**(ThreeHalves)
          NukePole=-GM%AtNum%I(At)
          Tau=Thresholds%TwoE
          DP2=(GM%AtNum%I(At)/Tau)**(Two/DBLE(SPEll+2))
!         Set atomic "primitive"  
          Prim%P=GM%Carts%D(:,At)
          Prim%Zeta=NuclearExpnt
          PBox%BndBox(:,1)=Prim%P
          PBox%BndBox(:,2)=Prim%P
          PBox=ExpandBox(PBox,Extent(0,NuclearExpnt,ABS(NukePole),TauPAC)
          HGKet(1)=Zero
          SPKetC(0)=Zero       
#ifdef PERIODIC
          QC(:)=GM%Carts%D(:,At)
          DO NC=1,CSMM1%NCells
!            Set Atomic Coordinates
             Prim%P(:)=QC(:)+CSMM1%CellCarts%D(:,NC)
             PBox%Center=Prim%P
!            Walk the walk
             CALL VWalk(PoleRoot)
          ENDDO
!         Reset the Atomic Coordinates
          Prim%P(:) = QC(:)
          PBox%Center=Prim%P
#else
!         Walk the walk
          CALL VWalk(PoleRoot)
#endif
!         Accumulate the atomic contribution
          NukE=NukE+NukeCo*HGKet(1)+NukePole*SPKetC(0) 
     END FUNCTION NukE
!
END MODULE
