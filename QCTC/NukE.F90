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
  USE PBCFarField
#ifdef PARALLEL
  USE ParallelQCTC
#endif
  IMPLICIT NONE
!----------!
  CONTAINS !
!=============================================================================================
!
!=============================================================================================
    FUNCTION NukE(GMLoc)
       TYPE(CRDS),  INTENT(IN)         :: GMLoc
       REAL(DOUBLE)                    :: NukE,NukeCo,NukePole,PExtent
       REAL(DOUBLE),DIMENSION(1:1)     :: HGBra
       REAL(DOUBLE),DIMENSION(0:0)     :: SPBraC,SPBraS
       INTEGER                         :: NC
       REAL(DOUBLE),DIMENSION(3)       :: PTmp
!---------------------------------------------------------------------------------------------
       NukE=Zero 
       DO At=1,GMLoc%Natms
#ifdef PARALLEL
      IF(At >= BegAtInd%I(MyID) .AND. At <= EndAtInd%I(MyID)) THEN
#endif
          IF(GMLoc%AtNum%D(At)<105.D0)THEN
!            Initialize |BRA>
             HGBra(1) =-GMLoc%AtNum%D(At)*(NuclearExpnt/Pi)**(ThreeHalves)
             SPBraC(0)=-GMLoc%AtNum%D(At)
             Prim%Ell =0
             Prim%P   =GMLoc%Carts%D(:,At)
             Prim%Zeta=NuclearExpnt
!            Set the MAC
             DP2 = (FudgeFactorial(0,SPELL+1)*ABS(GMLoc%AtNum%D(At))/TauMAC)**(Two/DBLE(SPEll+2))
             DP2 = MIN(1.D10,DP2)
!            Set the PAC
#ifdef NewPAC
             PrimBeta  = Prim%Zeta
             PrimWCoef = ABS(GMLoc%AtNum%D(At))
#else
             PExtent=Extent(0,NuclearExpnt,HGBra,TauPAC)
#endif
!            Initialize <KET|
             CALL SetKet(Prim,PExtent)
             PTmp=GMLoc%Carts%D(:,At)
             DO NC=1,CS_IN%NCells
!               Set Atomic Coordinates
                Prim%P=PTmp+CS_IN%CellCarts%D(:,NC) 
                PBox%Center= Prim%P
!               Walk the walk
                CALL VWalk(PoleRoot)
             ENDDO 
!            Reset the Atomic Coordinates
             Prim%P=PTmp
!            Accumulate the atomic contribution
             NukE=NukE+HGBra(1)*HGKet(1)+SPBraC(0)*SPKetC(0)
!            Add in the Far Field, Dipole and Quadripole  Correction 
             IF(GMLoc%PBC%Dimen > 0) THEN
                NukE = NukE + CTraxFF(Prim,HGBra,GMLoc)
             ENDIF
          ENDIF
#ifdef PARALLEL
       ENDIF
#endif
       ENDDO
!
     END FUNCTION NukE
!
END MODULE
