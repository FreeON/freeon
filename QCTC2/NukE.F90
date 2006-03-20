!    COMPUTE THE NUCLEAR-TOTAL ELECTROSTATIC ENERGY IN O(N Lg N) CPU TIME
!    Authors: Matt Challacombe and C.J. Tymczak
!==============================================================================
MODULE NukularE
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
    INTEGER                         :: NC,AT
    TYPE(QCPrim)                    :: QP
    REAL(DOUBLE),DIMENSION(3)       :: PTmp
    !---------------------------------------------------------------------------------------------
    NukE=Zero 
    DO At=1,GMLoc%Natms 
       IF(GMLoc%AtNum%D(At)<105.D0)THEN
          ! Initialize <BRA|
          HGBra(1) =-GMLoc%AtNum%D(At)*(NuclearExpnt/Pi)**(ThreeHalves)
          SPBraC(0)=-GMLoc%AtNum%D(At)
          QP%Prim%Ell = 0
          QP%Prim%P   = GMLoc%Carts%D(:,At)
          QP%Prim%Zeta= NuclearExpnt
#ifdef PAC_DEBUG
          ERRBRA=Zero
#endif
#ifdef MAC_DEBUG
          SPERRORBRAC=Zero
          SPERRORBRAS=Zero
#endif
#ifdef PAC_DEBUG
          ERRBRA(1)=ABS(HGBra(1))
#endif
#ifdef MAC_DEBUG
          SPERRORBRAC(0)=SPBraC(0)
          SPERRORBRAS(0)=Zero
#endif
          ! Set MAC and PAC
          QP%PAC%Zeta=QP%Prim%Zeta
          QP%PAC%Wght=GMLoc%AtNum%D(At)
          QP%MAC%O(0)=GMLoc%AtNum%D(At)
          QP%MAC%Delta=Zero
          QP%IHalf=ABS(HGBra(1)) !??
          ! Initialize the ket
          HGKet(1)=Zero
          SPKetC(0)=Zero
          SPKetS(0)=Zero
          ! Lay down the potential
          PTmp=GMLoc%Carts%D(:,At)
          DO NC=1,CS_IN%NCells
             QP%Prim%P=PTmp+CS_IN%CellCarts%D(:,NC) 
             CALL JWalk2(QP,PoleRoot,Nucular_O=.TRUE.)
          ENDDO
          ! Reset the primitive coordinates
          QP%Prim%P=PTmp
          ! Accumulate the atomic contribution
          NukE=NukE+HGBra(1)*HGKet(1)+SPBraC(0)*SPKetC(0)
         !  Add in the Far Field, Dipole and Quadripole  Correction 
          IF(GMLoc%PBC%Dimen > 0) THEN
             NukE = NukE + CTraxFF(QP%Prim,HGBra,GMLoc)
          ENDIF
       ENDIF
    ENDDO
  END FUNCTION NukE
END MODULE NukularE


!  NukeE =   -3.3319240537953522E+02
!  NukeE =   -4.1961861220141731E+02
!  NukeE =   -7.5267525922534048E+02
!  NukeE =   -8.3887193954951090E+02
!  NukeE =   -9.2510685429043531E+02
!  NukeE =   -1.0115307860803191E+03
!  NukeE =   -1.3442293987235985E+03
!  NukeE =   -1.4306381470191338E+03
!  NukeE =   -1.5169285373755615E+03
