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
  IMPLICIT NONE
!----------!
  CONTAINS !
!=============================================================================================
!
!=============================================================================================
    FUNCTION NukE() 
       REAL(DOUBLE) :: NukE,Tau,NukeCo,NukePole
!---------------------------------------------------------------------------------------------
!      Set MAC and PAC parameters
       NukE=Zero
       DO At=1,NAtoms
          NukeCo=-GM%AtNum%I(At)*(NuclearExpnt/Pi)**(ThreeHalves)
          NukePole=-GM%AtNum%I(At)
          Tau=Thresholds%TwoE
          DP2=(GM%AtNum%I(At)/Tau)**(Two/DBLE(SPEll+1))
          PoleSwitch=MAX(1.D-12,-LOG(Two*Tau/GM%AtNum%I(At)))
!         Set atomic "primitive"
          Prim%P=GM%Carts%D(:,At) 
          Prim%Zeta=NuclearExpnt
          HGKet(1)=Zero
          SPKetC(0)=Zero
!         Walk the walk
          CALL VWalk(PoleRoot)
!         Accumulate the atomic contribution
          NukE=NukE+NukeCo*HGKet(1)+NukePole*SPKetC(0)
       ENDDO
     END FUNCTION NukE
!=============================================================================================
!
!=============================================================================================
END MODULE
