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
       REAL(DOUBLE)                   :: NukE
       REAL(DOUBLE),DIMENSION(1)      :: NukeCo,NukePole
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
          DP2=(GM%AtNum%I(At)/TauMAC)**(Two/DBLE(SPEll+2))
!         Set atomic "primitive"  
          Prim%P=GM%Carts%D(:,At)
          Prim%Zeta=NuclearExpnt
          PBox%BndBox(:,1)=Prim%P
          PBox%BndBox(:,2)=Prim%P
          PBox=ExpandBox(PBox,Extent(0,NuclearExpnt,NukeCo,TauPAC))
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
          NukE=NukE+NukeCo(1)*HGKet(1)+NukePole(1)*SPKetC(0) 
       ENDDO
     END FUNCTION NukE
!
END MODULE
