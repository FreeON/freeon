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
    FUNCTION NukE(GM_Loc)
       TYPE(CRDS),  INTENT(IN)  :: GM_Loc
       REAL(DOUBLE)                    :: NukE,NukeCo,NukePole,PExtent
       REAL(DOUBLE),DIMENSION(1:1)     :: HGBra
       REAL(DOUBLE),DIMENSION(0:0)     :: SPBraC,SPBraS
#ifdef PERIODIC
       INTEGER                         :: NC
       REAL(DOUBLE),DIMENSION(3)       :: PTmp
#endif
!---------------------------------------------------------------------------------------------
       NukE=Zero 
       DO At=1,GM_Loc%Natms
	IF(GM_Loc%AtNum%D(At)<105.D0)THEN
!         Initialize |BRA>
          HGBra(1) =-GM_Loc%AtNum%D(At)*(NuclearExpnt/Pi)**(ThreeHalves)
          SPBraC(0)=-GM_Loc%AtNum%D(At)
          Prim%Ell=0
          Prim%P=GM_Loc%Carts%D(:,At)
          Prim%Zeta=NuclearExpnt
!         Set the MAC
!         Set the MAC
          DP2 = (FudgeFactorial(0,SPELL+1)*ABS(GM_Loc%AtNum%D(At))/TauMAC)**(Two/DBLE(SPEll+2))
          DP2 = MIN(1.D10,DP2)
!         Set the PAC
          PExtent=Extent(0,NuclearExpnt,HGBra,TauPAC)
!         Initialize <KET|
          CALL SetKet(Prim,PExtent)
#ifdef PERIODIC
          PTmp=GM_Loc%Carts%D(:,At)
          DO NC=1,CS_IN%NCells
!            Set Atomic Coordinates
             Prim%P=PTmp+CS_IN%CellCarts%D(:,NC)
             PBox%Center= Prim%P
!            Walk the walk
             CALL VWalk(PoleRoot)
          ENDDO
!         Reset the Atomic Coordinates
          Prim%P=PTmp
!         Accumulate the atomic contribution
          NukE=NukE+HGBra(1)*HGKet(1)+SPBraC(0)*SPKetC(0)
!         Add in the Far Field, Dipole and Quadripole  Correction
          IF(GM_Loc%PBC%Dimen > 0) THEN
             NukE = NukE + CTraxFF(Prim,HGBra,GM_Loc)
          ENDIF
#else
!         Walk the walk
          CALL VWalk(PoleRoot)
!         Accumulate the atomic contribution
          NukE=NukE+HGBra(1)*HGKet(1)+SPBraC(0)*SPKetC(0) 
#endif
	ENDIF
       ENDDO
     END FUNCTION NukE
!
END MODULE
