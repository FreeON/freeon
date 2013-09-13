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
  FUNCTION NukE(GMLoc,R,NoWrap_O)
    TYPE(CRDS),  INTENT(IN)         :: GMLoc
    REAL(DOUBLE)                    :: PiZ,NukE,NukeCo,NukePole,PExtent
    REAL(DOUBLE),DIMENSION(1:1)     :: HGBra
    REAL(DOUBLE),DIMENSION(0:0)     :: SPBraC,SPBraS
    INTEGER                         :: NC,AT
    TYPE(QCPrim)                    :: QP
    REAL(DOUBLE),DIMENSION(3)       :: PTmp

    TYPE(HGRho),OPTIONAL            :: R
    REAL(DOUBLE)                  :: XI,XM,CI,CM,ZZ,Zeta,PExt
    INTEGER      :: I,J,K,L,M,N,IAdd,JAdd
    LOGICAL, OPTIONAL         :: NoWrap_O
    LOGICAL                   :: NoWrap
    !---------------------------------------------------------------------------------------------
    IF(PRESENT(NoWrap_O))THEN
       NoWrap=NoWrap_O
    ELSE
       NoWrap=.FALSE.
    ENDIF
    NukE=Zero
    DO At=1,GMLoc%Natms
       IF(GMLoc%AtNum%D(At)<105.D0.AND.GMLoc%AtNum%D(At)>=0D0)THEN
          ! Initialize <BRA|
          PiZ=(Pi/NuclearExpnt)**(ThreeHalves)
          HGBra(1) =-GMLoc%AtNum%D(At)/PiZ
          SPBraC(0)=-GMLoc%AtNum%D(At)
          QP%Prim%Ell=0
          QP%Prim%P=GMLoc%Carts%D(:,At)
          QP%Prim%Zeta=NuclearExpnt
          CALL PWrap(GM,QP%Prim,.NOT.NoWrap)
          !
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
          QP%MAC%O(0)=GMLoc%AtNum%D(At)
          QP%MAC%Delta=Zero
          PExt=Extent(0,NuclearExpnt,HGBra,TauPAC,Potential_O=.TRUE.)
          ! Initialize the ket
          HGKet(1)=Zero
          SPKetC(0)=Zero
          SPKetS(0)=Zero
          ! Lay down the potential
          PTmp=QP%Prim%Pw
          DO NC=1,CS_IN%NCells
             !
             QP%Prim%Pw=PTmp+CS_IN%CellCarts%D(:,NC)
             !
             QP%Box%BndBox(:,1)=QP%Prim%Pw
             QP%Box%BndBox(:,2)=QP%Prim%Pw
             QP%Box=ExpandBox(QP%Box,PExt)
             !
#ifdef MAC_DEBUG
             CALL DBL_VECT_EQ_DBL_SCLR(SPLen+1,SPKetC(0),Zero)
             CALL DBL_VECT_EQ_DBL_SCLR(SPLen+1,SPKetS(0),Zero)
             SPErrorKetS=Zero
             SPErrorKetC=Zero
#endif
             CALL JWalk2(QP,PoleRoot,Nucular_O=.TRUE.)
          ENDDO
          ! Reset the primitive coordinates
          QP%Prim%Pw=PTmp
          ! Accumulate the atomic contribution
          NukE=NukE+HGBra(1)*HGKet(1)+SPBraC(0)*SPKetC(0)
          !  Add in the Far Field, Dipole and Quadripole  Correction

          IF(GMLoc%PBC%Dimen > 0) THEN
             NukE = NukE + CTraxFF(QP%Prim,GMLoc,HGBra,PiZ)
          ENDIF

!!$          WRITE(*,*)'================================================='
!!$          WRITE(*,*)' Nuke = ',CTraxFF(QP%Prim,HGBra,GMLoc)
!!$          Nuke=0
!!$          IF(PRESENT(R))THEN
!!$             XI=GMLoc%Carts%D(1,At)
!!$             CI=-GMLoc%AtNum%D(At)
!!$             DO M=1,R%NExpt
!!$                Zeta=R%Expt%D(M)
!!$                ZZ=(Pi/R%Expt%D(M))**1.5D0
!!$                DO N=1,R%NQ%I(M)
!!$                   IAdd=R%OffQ%I(M)+N
!!$                   JAdd=R%OffR%I(M)+(N-1)+1
!!$                   XM=Rho%Qx%D(IAdd)
!!$                   CM=Rho%Co%D(JAdd)*ZZ
!!$                   DO L=-5000,-2,1
!!$                      NukE=NukE+CI*CM/ABS(XI-(XM+DBLE(L)))
!!$                   ENDDO
!!$                   DO L=5000,2,-1
!!$                      NukE=NukE+CI*CM/ABS(XI-(XM+DBLE(L)))
!!$                   ENDDO
!!$                ENDDO
!!$             ENDDO
!!$          ENDIF
!!$          WRITE(*,*)' Nuke = ',NukE


       ENDIF
    ENDDO
  END FUNCTION NukE
END MODULE NukularE
