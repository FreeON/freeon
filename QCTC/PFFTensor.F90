!
!  Matt Challacombe and CJ Tymczak
!
MODULE PFFTen
  USE Derivedtypes
  USE GlobalScalars
  USE GlobalObjects
  USE ProcessControl
  USE Indexing
  USE Parse
  USE InOut
  USE Macros
  USE Thresholding
  USE BoundingBox
  USE SpecFun
  USE MondoPoles
  USE AtomPairs
  USE PlaneWise
  IMPLICIT NONE
  ! Globals
  REAL(DOUBLE),PARAMETER              :: Small=1.0D-12
  CHARACTER(LEN=7),PARAMETER          :: Cube   ='Cube   '
  CHARACTER(LEN=7),PARAMETER          :: Rect   ='Rect   '
  CHARACTER(LEN=7),PARAMETER          :: NonCube='NonCube'
  CHARACTER(LEN=7),PARAMETER          :: Cube111='Cube111'
  CHARACTER(LEN=7),PARAMETER          :: Cube112='Cube112'
  CHARACTER(LEN=7),PARAMETER          :: Cube122='Cube122'
  CHARACTER(LEN=7),PARAMETER          :: Cube222='Cube222'
  !
  INTEGER,PARAMETER                   :: PlaneWiseEll=0
  INTEGER,PARAMETER                   :: LSwitch=22
  REAL(DOUBLE),PARAMETER              :: Accuracy = 1.D-20
  !
CONTAINS
!========================================================================================
!
!========================================================================================
  SUBROUTINE CalculatePFFT(MaxEll,GM,Args,CS,TenC,TenS)
    TYPE(CRDS)             :: GM
    TYPE(ARGMT)            :: Args
    TYPE(DBL_VECT)         :: TenC,TenS,TenC2,TenS2
    TYPE(CellSet)          :: CS
    INTEGER                :: MaxEll , L,M,MM,LM
    INTEGER                :: Layers,NC
    REAL(DOUBLE)           :: Scale,R1,R2,Radius,Rx,Ry,Rz,ab
    CHARACTER(LEN=7)       :: CellType
!-------------------------------------------------------------------------------------!
!   Determine Cell Type

    CellType=DetCellType(GM,CS%NCells,Scale)
!   Read from Archive or Make the Tensor
!!$    IF(CellType==Cube111) THEN
!!$       CALL GetExactTensor1L3D(2*MaxEll,GM,Scale,TenC,TenS)
!!$    ELSEIF(CellType==Cube222) THEN
!!$       CALL GetExactTensor2L3D(2*MaxEll,GM,Scale,TenC,TenS)
!!$    ELSE
    IF(GM%PBC%Dimen==1) THEN
       CALL MakeTensor1D(2*MaxEll,GM,Args,CS,TenC,TenS)
       !          CALL DirectTensor1D(2*MaxEll,GM,Args,CS,TenC,TenS)
    ELSEIF(GM%PBC%Dimen==2) THEN
       CALL MakeTensor2D(2*MaxEll,GM,Args,CS,TenC,TenS)
    ELSEIF(GM%PBC%Dimen==3) THEN
       CALL MakeTensor3D(2*MaxEll,GM,Args,CS,TenC,TenS)
!!$
!!$       CALL New(TenC2,LSP(2*MaxEll),0)
!!$       CALL New(TenS2,LSP(2*MaxEll),0)
!!$       CALL GetExactTensor1L3D(2*MaxEll,GM,Scale,TenC2,TenS2)
!!$
!!$       DO L=4,80!MaxEll
!!$          DO M=0,L
!!$             LM=LTD(L)+M
!!$             IF(ABS(TenC2%D(LM)).NE.0D0)THEN
!!$                WRITE(*,*)LM,TENC2%D(LM)
!!$                ab=ABS((TenC%D(LM)-TenC2%D(LM))/TenC2%D(LM))
!!$                IF(ab>1D-14)THEN
!!$                   WRITE(*,*)'L = ',L,' M = ',M,' Ten = ',TenC%D(LM),TenC2%D(LM)
!!$                   DO MM=0,L
!!$                      LM=LTD(L)+MM
!!$                      WRITE(*,*)'L = ',L,' M = ',MM,' Ten = ',TenC%D(LM),TenC%D(LM)-TenC2%D(LM)
!!$                   ENDDO
!!$                   STOP 'STOP STOP STOP STOP STOP STOP STOP STOP '
!!$                ENDIF
!!$             ENDIF
!!$          ENDDO
!!$       ENDDO
!!$       STOP 'STOP STOP STOP STOP STOP STOP STOP STOP '

   ENDIF

  END SUBROUTINE CalculatePFFT
!========================================================================================
!
!========================================================================================
  SUBROUTINE GetExactTensor1L3D(MaxL,GM,Scale,TenC,TenS)
    INTEGER              :: MaxL
    TYPE(CRDS)           :: GM
    TYPE(DBL_VECT)       :: TenC,TenS
    INTEGER              :: I,L,M,LM,K,NumberOfRecords
    REAL(DOUBLE)         :: Scale
!-------------------------------------------------------------------------------------!
!
    INCLUDE "PBCTensor/Majik_Kubic_WS1.Inc"
!
    TenC%D = Zero
    TenS%D = Zero
!
!   Load Majik numbers from parameter statements
!
    K=0
    DO L=4,MIN(128,MaxL),2
       DO M=0,L,4
          K=K+1
          LM=LTD(L)+M
          TenC%D(LM)=Majik1(K)/(Scale**(DBLE(L)+One))
          IF(L>80)TenC%D(LM)=0D0
       ENDDO
    ENDDO
!
  END SUBROUTINE GetExactTensor1L3D
!========================================================================================
!
!========================================================================================
  SUBROUTINE GetExactTensor2L3D(MaxL,GM,Scale,TenC,TenS)
    INTEGER              :: MaxL
    TYPE(CRDS)           :: GM
    TYPE(DBL_VECT)       :: TenC,TenS
    INTEGER              :: I,L,M,LM,K,NumberOfRecords
    REAL(DOUBLE)         :: Scale
!-------------------------------------------------------------------------------------!
!
    INCLUDE "PBCTensor/Majik_Kubic_WS2.Inc"
!
    TenC%D = Zero
    TenS%D = Zero
!
!     Load Majik numbers from parameter statements
!
    K=0
    DO L=4,MIN(128,MaxL),2
       DO M=0,L,4
          K=K+1
          LM=LTD(L)+M
          TenC%D(LM)=Majik2(K)/(Scale**(DBLE(L)+One))
       ENDDO
    ENDDO
!
  END SUBROUTINE GetExactTensor2L3D
!========================================================================================
! Calculate the PFFTensor 1D
!========================================================================================
  SUBROUTINE MakeTensor1D(MaxL,GM,Args,CS,TenC,TenS)
    INTEGER                           :: MaxL,LenL
    INTEGER                           :: I,J,L,M,LM,NC
    TYPE(CellSet)                     :: CS, CSMM
    TYPE(DBL_VECT)                    :: TenC,TenS
    TYPE(CRDS)                        :: GM
    TYPE(ARGMT)                       :: Args
!-------------------------------------------------------------------------------------!
!
!   Number of Inner Boxes
!
    LenL=LSP(MaxL)
    TenC%D=Zero
    NC = (CS%NCells-1)/2
!
!   One Dimension: Lef and Right
!
    IF(GM%PBC%AutoW%I(1)==1) THEN
       CALL IrRegular(MaxL, GM%PBC%BoxShape%D(1,1),Zero,Zero)
       TenC%D(0:LenL) = Cpq(0:LenL)
       CALL IrRegular(MaxL,-GM%PBC%BoxShape%D(1,1),Zero,Zero)
       TenC%D(0:LenL) = TenC%D(0:LenL)+Cpq(0:LenL)
    ELSEIF(GM%PBC%AutoW%I(2)==1) THEN
       CALL IrRegular(MaxL,Zero, GM%PBC%BoxShape%D(2,2),Zero)
       TenC%D(0:LenL) = Cpq(0:LenL)
       CALL IrRegular(MaxL,Zero,-GM%PBC%BoxShape%D(2,2),Zero)
       TenC%D(0:LenL) = TenC%D(0:LenL)+Cpq(0:LenL)
    ELSEIF(GM%PBC%AutoW%I(3)==1) THEN
       CALL IrRegular(MaxL,Zero,Zero, GM%PBC%BoxShape%D(3,3))
       TenC%D(0:LenL) = Cpq(0:LenL)
       CALL IrRegular(MaxL,Zero,Zero,-GM%PBC%BoxShape%D(3,3))
       TenC%D(0:LenL) = TenC%D(0:LenL)+Cpq(0:LenL)
    ENDIF

    DO L=1,MaxL
       DO M = 0,L
          LM = LTD(L)+M
          TenC%D(LM) = TenC%D(LM)*RZeta(L+1,NC)
       ENDDO
    ENDDO
!





  END SUBROUTINE MakeTensor1D
!========================================================================================
! Calculate the PFFTensor 2D
!========================================================================================
  SUBROUTINE MakeTensor2D(MaxL,GM,Args,CS,TenC,TenS)
    INTEGER                           :: MaxL
    INTEGER                           :: I,J,K,L,M,LM,NC
    REAL(DOUBLE),DIMENSION(3)         :: PQ
    TYPE(DBL_VECT)                    :: TenC,TenS
    REAL(DOUBLE)                      :: CFac,SFac,BetaSq,Rad,RadSq,ExpFac
    TYPE(CellSet)                     :: CS, CSMM
    REAL(DOUBLE)                      :: Rmin,RMAX,Accuracy,LenScale,SUM
!
    INTEGER                           :: NDiv
    REAL(DOUBLE)                      :: LenMax,Delt,IntFac
    REAL(DOUBLE),DIMENSION(3)         :: Vec
    REAL(DOUBLE),DIMENSION(3,3)       :: RecpLatVec,LatVec
    TYPE(CRDS)                        :: GM
    TYPE(ARGMT)                       :: Args
    COMPLEX(DOUBLE),DIMENSION(0:MaxL,0:MaxL) :: Ylm
    !-------------------------------------------------------------------------------------
    TenC%D=Zero
    TenS%D=Zero
    !
    DO I = 1,3
       DO J = 1,3
          RecpLatVec(I,J)=GM%PBC%InvBoxSh%D(J,I)
          LatVec(I,J)    =GM%PBC%BoxShape%D(I,J)
       ENDDO
    ENDDO
    !---------------------------------------------------------------------------------
    ! THE PLANEWISE LATTICE SUMS, PLANE <<MUST>> BE THE X,Y PLANE
    !---------------------------------------------------------------------------------
    Ylm=0D0
    CALL InnPlane(MaxL,LatVec(:,1),LatVec(:,2),RecpLatVec(:,1),RecpLatVec(:,2),Ylm, &
                  WellSepCells_O=CS) ! Last bit causes subtraction of near-field
    !---------------------------------------------------------------------------------
    DO L = 0,MaxL
       DO M = 0,L
          LM = LTD(L)+M
          TenC%D(LM)=(1D0, 0D0)*Ylm(L,M)*Factorial(L-M)
          TenS%D(LM)=(0D0,-1D0)*Ylm(L,M)*Factorial(L-M)
       ENDDO
    ENDDO
!!$    !---------------------------------------------------------------------------------
!!$    ! SUBTRACT THE RECIPROCAL SPACE PART OF THE NEAR FIELD FROM THE PLANEWISE CELLS
!!$    !---------------------------------------------------------------------------------
!!$    DO NC = 1,CS%NCells-1
!!$       PQ(:) = CS%CellCarts%D(:,NC)
!!$       CALL IrRegular(MaxL,PQ(1),PQ(2),PQ(3))
!!$       DO L = 1,MaxL
!!$          DO M = 0,L
!!$             LM = LTD(L)+M
!!$             TenC%D(LM)=TenC%D(LM)-Cpq(LM)
!!$             TenS%D(LM)=TenS%D(LM)-Spq(LM)
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
!
  END SUBROUTINE MakeTensor2D
!========================================================================================
! Calculate the PFFTensor 3D
!========================================================================================
  SUBROUTINE MakeTensor3D(MaxL,GM,Args,CS,TenC,TenS)
    INTEGER                           :: MaxL
    INTEGER                           :: I,J,K,L,M,LM,NC
    TYPE(DBL_VECT)                    :: TenC,TenS
    TYPE(CellSet)                     :: CS, CSMM
    REAL(DOUBLE),DIMENSION(3)         :: PQ
    REAL(DOUBLE)                      :: CFac,SFac,CFac2,SFac2,BetaSq,Rad,RadSq,ExpFac
    REAL(DOUBLE)                      :: R,Rmin,KMin,KMax,RMax,Beta
    INTEGER                           :: NDiv
    REAL(DOUBLE)                      :: LenMax,Delt,IntFac
    REAL(DOUBLE),DIMENSION(3)         :: Vec
    REAL(DOUBLE),DIMENSION(3,3)       :: RecpLatVec,LatVec
    TYPE(CRDS)                        :: GM
    TYPE(ARGMT)                       :: Args

    COMPLEX(DOUBLE),DIMENSION(0:PlaneWiseEll,0:PlaneWiseEll) :: Ylm
    COMPLEX(DOUBLE) :: CmpxI=(0D0,1D0),RealI=(1D0,0D0)

    INTEGER :: IFF=10
    LOGICAL :: InCell,OutCell
    REAL(DOUBLE) :: X,Y,Z
    !-------------------------------------------------------------------------------------!
    !
    !   Initialize
    !
    TenC%D=Zero
    TenS%D=Zero

    !
    DO I = 1,3
       DO J = 1,3
          RecpLatVec(I,J) = GM%PBC%InvBoxSh%D(J,I)
          LatVec(I,J)     = GM%PBC%BoxShape%D(I,J)
       ENDDO
    ENDDO
    !-------------------------------------------------------------------
    !  THIS NEEDS TO BE TESTED FOR REALLY FLATTENED CELLS:
    !-------------------------------------------------------------------
    KMin=1D10
    RMin=1D10
    DO J=1,3
       R=SQRT(LatVec(1,J)**2+LatVec(2,J)**2+LatVec(3,J)**2)
       RMin=MIN(RMin,R)
    ENDDO
    ! This is the beta that hopefully balances the reciprocal and real space
    ! NBDW lattice sums
    Beta=SQRT(Pi)/RMin
    !-----------------------------------------------------------------------
    ! OK, this bit is important.  RMax is the distance the real real space
    ! lattice sum is carried out over.  HOWEVER, for large L, direct summation
    ! is used, and this MAY BE LARGER than the radius given by the balanced
    ! value.
    !
    ! Here are the balanced values:
    RMax=(One/Beta)*SQRT(ABS(LOG(Accuracy)))
    KMax=    Beta *SQRT(ABS(LOG(Accuracy)))/Pi
    ! And here is the value we use, which takes into account direct summation
    ! in the first loop below (CFac=1).
    RMax=MAX(RMax,Accuracy**(-One/DBLE(LSwitch+1)))
!!$
!!$    RMax=3D0*RMax
!!$    KMax=3D0*KMax
!!$
    !---------------------------------------------------------------------------------
    ! REAL SPACE
    !---------------------------------------------------------------------------------
    BetaSq=Beta**2
    CALL SetCellSets(CSMM,GM%PBC%AutoW%I,GM%PBC%BoxShape%D,'Tensor3D[R]','Radial',RMin_O=RMax)
    !
    DO NC = 1,CSMM%NCells-1 ! Leaves out central cell (CellSet is sorted, so PQ=0 is last)
       PQ(:) = CSMM%CellCarts%D(:,NC)
       IF(.NOT. InCell_CellSet(CS,PQ(1),PQ(2),PQ(3))) THEN
          RadSq = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
          CALL IrRegular(MaxL,PQ(1),PQ(2),PQ(3))
          DO L = 1,MaxL
             IF(L .LE. LSwitch) THEN
                CFac=RegularizedGammaQ(DBLE(L)+5D-1,RadSq)
!!$                CFac = GScript(L,RadSq)
             ELSE
                CFac = One
             ENDIF
             DO M = 0,L
                LM = LTD(L)+M
                TenC%D(LM)=TenC%D(LM)+Cpq(LM)*CFac
                TenS%D(LM)=TenS%D(LM)+Spq(LM)*CFac
             ENDDO
          ENDDO
       ENDIF
    ENDDO

    CALL Delete_CellSet(CSMM)
    !---------------------------------------------------------------------------------
    ! RECIPROCAL SPACE
    !---------------------------------------------------------------------------------
    ExpFac = Pi*Pi/BetaSq
    !
    CALL SetCellSets(CSMM,GM%PBC%AutoW%I,RecpLatVec,'Tensor3D[K]','Radial',RMin_O=KMax)
    !
    DO NC = 1,CSMM%NCells
       PQ(:) = CSMM%CellCarts%D(:,NC)
       Rad   = SQRT(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
       IF(Rad .GT. 1.D-14) THEN
          CALL IrRegular(LSwitch,PQ(1),PQ(2),PQ(3))
          DO L = 1,LSwitch
!!$
!!$          CFac2 = (SQRT(Two*Pi)**3/Two**(L-Half))/DGamma(DBLE(L)+Half) &
!!$               *(Two*Pi*Rad)**(L-2)*EXP(-ExpFac*Rad**2) &
!!$               * CmpxI**(L)*RealI/GM%PBC%CellVolume
!!$
!!$          SFac2 = -(SQRT(Two*Pi)**3/Two**(L-Half))/DGamma(DBLE(L)+Half) &
!!$               *(Two*Pi*Rad)**(L-2)*EXP(-ExpFac*Rad**2) &
!!$               * CmpxI**(L)*CmpxI/GM%PBC%CellVolume

             CFac = FT_FScriptC_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
             SFac = FT_FScriptS_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
             DO M = 0,L
                LM = LTD(L)+M
                TenC%D(LM)=TenC%D(LM)+Cpq(LM)*CFac-Spq(LM)*SFac
                TenS%D(LM)=TenS%D(LM)+Spq(LM)*CFac+Cpq(LM)*SFac
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    !
    CALL Delete_CellSet(CSMM)
!!$    !---------------------------------------------------------------------------------
!!$    ! PLANEWISE LATTICE SUMS (CAPABILITY IS INCOMPLETE, AS SURFACE TERM IS NOT CODED)
!!$    !
!!$    !---------------------------------------------------------------------------------
!!$    Ylm=0D0
!!$    CALL OffPlane(PlaneWiseEll,LatVec(:,1),LatVec(:,2),LatVec(:,3),RecpLatVec(:,1),RecpLatVec(:,2),Ylm)
!!$    CALL InnPlane(PlaneWiseEll,LatVec(:,1),LatVec(:,2),RecpLatVec(:,1),RecpLatVec(:,2),Ylm)
!!$    !---------------------------------------------------------------------------------
!!$    DO L = 0,PlaneWiseEll
!!$       DO M = 0,L
!!$          LM = LTD(L)+M
!!$          TenC%D(LM)=(1D0,0D0)*Ylm(L,M)*Factorial(L-M)
!!$          TenS%D(LM)=(0D0,-1D0)*Ylm(L,M)*Factorial(L-M)
!!$       ENDDO
!!$    ENDDO
!!$    !---------------------------------------------------------------------------------
!!$    ! SUBTRACT THE NEAR FIELD FROM THE PLANEWISE CELLS
!!$    !---------------------------------------------------------------------------------
!!$    DO NC = 1,CS%NCells-1
!!$       PQ(:) = CS%CellCarts%D(:,NC)
!!$       CALL IrRegular(PlaneWiseEll,PQ(1),PQ(2),PQ(3))
!!$       DO L = 1,PlaneWiseEll
!!$          DO M = 0,L
!!$             LM = LTD(L)+M
!!$             TenC%D(LM)=TenC%D(LM)-Cpq(LM)
!!$             TenS%D(LM)=TenS%D(LM)-Spq(LM)
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
    !---------------------------------------------------------------------------------
    ! SUBTRACT THE NEAR FIELD RECIPROCAL COMPONENT FROM THE NIJBOER-DEWETTE CELLS
    !---------------------------------------------------------------------------------
    DO NC = 1,CS%NCells
       PQ(:) = CS%CellCarts%D(:,NC)
       RadSq = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
       IF(RadSq .GT. 1.D-14) THEN
          CALL IrRegular(LSwitch,PQ(1),PQ(2),PQ(3))
          DO L = 1,LSwitch
             CFac=RegularizedGammaP(DBLE(L)+5D-1,RadSq)
             DO M = 0,L
                LM = LTD(L)+M
                TenC%D(LM)=TenC%D(LM)-Cpq(LM)*CFac
                TenS%D(LM)=TenS%D(LM)-Spq(LM)*CFac
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    !---------------------------------------------------------------------------------
    ! ZERO LOW ORDER TERMS
    !---------------------------------------------------------------------------------
    DO L = 0,1
       DO M = 0,L
          LM = LTD(L)+M
          TenC%D(LM)=0D0
          TenS%D(LM)=0D0
       ENDDO
    ENDDO
    !
  END SUBROUTINE MakeTensor3D

#ifdef LEGACY_CRAP
  SUBROUTINE MakeTensor3D_MED_OLD(MaxL,GM,Args,CS,TenC,TenS)
    INTEGER                           :: MaxL
    INTEGER                           :: I,J,K,L,M,LM,NC
    TYPE(DBL_VECT)                    :: TenC,TenS
    TYPE(CellSet)                     :: CS, CSMM
    REAL(DOUBLE),DIMENSION(3)         :: PQ
    REAL(DOUBLE)                      :: CFac,SFac,CFac2,SFac2,BetaSq,Rad,RadSq,ExpFac
    REAL(DOUBLE)                      :: R,Rmin,KMin,KMax,RMax,Beta
    INTEGER                           :: NDiv
    REAL(DOUBLE)                      :: LenMax,Delt,IntFac
    REAL(DOUBLE),DIMENSION(3)         :: Vec
    REAL(DOUBLE),DIMENSION(3,3)       :: RecpLatVec,LatVec
    TYPE(CRDS)                        :: GM
    TYPE(ARGMT)                       :: Args

    COMPLEX(DOUBLE),DIMENSION(0:PlaneWiseEll,0:PlaneWiseEll) :: Ylm
    COMPLEX(DOUBLE) :: CmpxI=(0D0,1D0),RealI=(1D0,0D0)

    INTEGER :: IFF=10
    LOGICAL :: InCell,OutCell
    REAL(DOUBLE) :: X,Y,Z
    !-------------------------------------------------------------------------------------!
    !
    !   Initialize
    !
    TenC%D=Zero
    TenS%D=Zero

    !
    DO I = 1,3
       DO J = 1,3
          RecpLatVec(I,J) = GM%PBC%InvBoxSh%D(J,I)
          LatVec(I,J)     = GM%PBC%BoxShape%D(I,J)
       ENDDO
    ENDDO
    !-------------------------------------------------------------------
    !  THIS NEEDS TO BE TESTED FOR REALLY FLATTENED CELLS:
    !-------------------------------------------------------------------
    KMin=1D10
    RMin=1D10
    DO J=1,3
       R=SQRT(LatVec(1,J)**2+LatVec(2,J)**2+LatVec(3,J)**2)
       RMin=MIN(RMin,R)
    ENDDO
    Beta=SQRT(Pi)/RMin
    RMax=(One/Beta)*SQRT(ABS(LOG(Accuracy)))
    KMax=    Beta *SQRT(ABS(LOG(Accuracy)))/Pi

    !---------------------------------------------------------------------------------
    ! REAL SPACE
    !---------------------------------------------------------------------------------
    BetaSq=Beta**2
    CALL SetCellSets(CSMM,GM%PBC%AutoW%I,GM%PBC%BoxShape%D,'Radial',RMin_O=RMax)
    !
    DO NC = 1,CSMM%NCells-1 ! Leaves out central cell (CellSet is sorted, so PQ=0 is last)
       PQ(:) = CSMM%CellCarts%D(:,NC)
       IF(.NOT. InCell_CellSet(CS,PQ(1),PQ(2),PQ(3))) THEN
          RadSq = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
          CALL IrRegular(MaxL,PQ(1),PQ(2),PQ(3))
          DO L = 1,MaxL
             IF(L .LE. LSwitch) THEN
                CFac=RegularizedGammaQ(DBLE(L)+5D-1,RadSq)
!!$                CFac = GScript(L,RadSq)
             ELSE
                CFac = One
             ENDIF
             DO M = 0,L
                LM = LTD(L)+M
                TenC%D(LM)=TenC%D(LM)+Cpq(LM)*CFac
                TenS%D(LM)=TenS%D(LM)+Spq(LM)*CFac
             ENDDO
          ENDDO
       ENDIF
    ENDDO

    CALL Delete_CellSet(CSMM)
    !---------------------------------------------------------------------------------
    ! RECIPROCAL SPACE
    !---------------------------------------------------------------------------------
    ExpFac = Pi*Pi/BetaSq
    !
    CALL SetCellSets(CSMM,GM%PBC%AutoW%I,RecpLatVec,'Radial',RMin_O=KMax)
    !
    DO NC = 1,CSMM%NCells
       PQ(:) = CSMM%CellCarts%D(:,NC)
       Rad   = SQRT(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
       IF(Rad .GT. 1.D-14) THEN
          CALL IrRegular(MaxL,PQ(1),PQ(2),PQ(3))
          DO L = 1,LSwitch
!!$
!!$          CFac2 = (SQRT(Two*Pi)**3/Two**(L-Half))/DGamma(DBLE(L)+Half) &
!!$               *(Two*Pi*Rad)**(L-2)*EXP(-ExpFac*Rad**2) &
!!$               * CmpxI**(L)*RealI/GM%PBC%CellVolume
!!$
!!$          SFac2 = -(SQRT(Two*Pi)**3/Two**(L-Half))/DGamma(DBLE(L)+Half) &
!!$               *(Two*Pi*Rad)**(L-2)*EXP(-ExpFac*Rad**2) &
!!$               * CmpxI**(L)*CmpxI/GM%PBC%CellVolume

             CFac = FT_FScriptC_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
             SFac = FT_FScriptS_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
             DO M = 0,L
                LM = LTD(L)+M
                TenC%D(LM)=TenC%D(LM)+Cpq(LM)*CFac-Spq(LM)*SFac
                TenS%D(LM)=TenS%D(LM)+Spq(LM)*CFac+Cpq(LM)*SFac
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    !
    CALL Delete_CellSet(CSMM)
!!$    !---------------------------------------------------------------------------------
!!$    ! PLANEWISE LATTICE SUMS (CAPABILITY IS INCOMPLETE, AS SURFACE TERM IS NOT CODED)
!!$    !
!!$    !---------------------------------------------------------------------------------
!!$    Ylm=0D0
!!$    CALL OffPlane(PlaneWiseEll,LatVec(:,1),LatVec(:,2),LatVec(:,3),RecpLatVec(:,1),RecpLatVec(:,2),Ylm)
!!$    CALL InnPlane(PlaneWiseEll,LatVec(:,1),LatVec(:,2),RecpLatVec(:,1),RecpLatVec(:,2),Ylm)
!!$    !---------------------------------------------------------------------------------
!!$    DO L = 0,PlaneWiseEll
!!$       DO M = 0,L
!!$          LM = LTD(L)+M
!!$          TenC%D(LM)=(1D0,0D0)*Ylm(L,M)*Factorial(L-M)
!!$          TenS%D(LM)=(0D0,-1D0)*Ylm(L,M)*Factorial(L-M)
!!$       ENDDO
!!$    ENDDO
!!$    !---------------------------------------------------------------------------------
!!$    ! SUBTRACT THE NEAR FIELD FROM THE PLANEWISE CELLS
!!$    !---------------------------------------------------------------------------------
!!$    DO NC = 1,CS%NCells-1
!!$       PQ(:) = CS%CellCarts%D(:,NC)
!!$       CALL IrRegular(PlaneWiseEll,PQ(1),PQ(2),PQ(3))
!!$       DO L = 1,PlaneWiseEll
!!$          DO M = 0,L
!!$             LM = LTD(L)+M
!!$             TenC%D(LM)=TenC%D(LM)-Cpq(LM)
!!$             TenS%D(LM)=TenS%D(LM)-Spq(LM)
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
    !---------------------------------------------------------------------------------
    ! SUBTRACT THE NEAR FIELD RECIPROCAL COMPONENT FROM THE NIJBOER-DEWETTE CELLS
    !---------------------------------------------------------------------------------
    DO NC = 1,CS%NCells
       PQ(:) = CS%CellCarts%D(:,NC)
       RadSq = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
       IF(RadSq .GT. 1.D-14) THEN
          CALL IrRegular(MaxL,PQ(1),PQ(2),PQ(3))
          DO L = 1,LSwitch
!!$             CFac = FScript(L,RadSq)
             CFac=RegularizedGammaP(DBLE(L)+5D-1,RadSq)
             DO M = 0,L
                LM = LTD(L)+M
                TenC%D(LM)=TenC%D(LM)-Cpq(LM)*CFac
                TenS%D(LM)=TenS%D(LM)-Spq(LM)*CFac
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    !---------------------------------------------------------------------------------
    ! ZERO LOW ORDER TERMS
    !---------------------------------------------------------------------------------
    DO L = 0,1
       DO M = 0,L
          LM = LTD(L)+M
          TenC%D(LM)=0D0
          TenS%D(LM)=0D0
       ENDDO
    ENDDO
    !
  END SUBROUTINE MakeTensor3D_MED_OLD
#endif
!========================================================================================
!   Determine Cell Type
!========================================================================================
  FUNCTION DetCellType(GM,NC,Scale)
    TYPE(CRDS)                   :: GM
    INTEGER                      :: I,NC
    REAL(DOUBLE)                 :: MagA,MagB,MagC,AdotB,AdotC,BdotC,Scale
    CHARACTER(LEN=7)             :: DetCellType
!
!   Determine the cell type
!
    MagA  = Zero
    MagB  = Zero
    MagC  = Zero
    AdotB = Zero
    AdotC = Zero
    BdotC = Zero
    DO I=1,3
       MagA  = MagA+GM%PBC%BoxShape%D(I,1)*GM%PBC%BoxShape%D(I,1)
       MagB  = MagB+GM%PBC%BoxShape%D(I,2)*GM%PBC%BoxShape%D(I,2)
       MagC  = MagC+GM%PBC%BoxShape%D(I,3)*GM%PBC%BoxShape%D(I,3)
       AdotB = AdotB+GM%PBC%BoxShape%D(I,1)*GM%PBC%BoxShape%D(I,2)
       AdotC = AdotC+GM%PBC%BoxShape%D(I,1)*GM%PBC%BoxShape%D(I,3)
       BdotC = BdotC+GM%PBC%BoxShape%D(I,2)*GM%PBC%BoxShape%D(I,3)
    ENDDO
    MagA  = SQRT(MagA)
    MagB  = SQRT(MagB)
    MagC  = SQRT(MagC)
    AdotB = SQRT(AdotB/(MagA*MagB))
    AdotC = SQRT(AdotC/(MagA*MagC))
    BdotC = SQRT(BdotC/(MagB*MagC))
!
    DetCellType=NonCube
    IF(AdotB < Small .AND. AdotC < Small .AND. BdotC < Small) THEN
       IF( ABS(MagA-MagB) < Small .AND. ABS(MagA-MagC) < Small .AND. ABS(MagB-MagC) < Small) THEN
          IF(    NC==27 ) THEN
             DetCellType=Cube111
          ELSEIF(NC==125) THEN
             DetCellType=Cube222
          ELSE
             DetCellType=Cube
          ENDIF
          Scale = (MagA+MagB+MagC)/Three
       ELSE
          DetCellType=Rect
       ENDIF
    ELSE
       DetCellType=NonCube
    ENDIF
!
  END FUNCTION DetCellType
!========================================================================================
!   FT_FSCriptC
!========================================================================================
  FUNCTION FT_FScriptC_3D(L,ExpFac,R)
    INTEGER                    :: L,IFac,Isgn
    REAL(DOUBLE)               :: ExpFac,R,FT_FScriptC_3D,Fac
    REAL(DOUBLE),DIMENSION(41) :: Sfac = (/2.00000000000000000D0 ,1.61199195401646964D0 , &
         1.39399011673806825D0 ,1.24836057075206815D0 ,1.14180213615392840D0 , &
         1.05927697236749128D0 ,9.92824332104422014D-1,9.37767832227969827D-1, &
         8.91149652589504141D-1,8.50992321055495993D-1,8.15915401864279526D-1, &
         7.84921233738737833D-1,7.57267966537096089D-1,7.32390674279586654D-1, &
         7.09850372573553877D-1,6.89299957250757043D-1,6.70460791874183787D-1, &
         6.53106213773379367D-1,6.37049661171763433D-1,6.22135962734544435D-1, &
         6.08234838286273008D-1,5.95235975457893168D-1,5.83045248969098145D-1, &
         5.71581781316674666D-1,5.60775631817133548D-1,5.50565960943842365D-1, &
         5.40899558418448758D-1,5.31729652703953292D-1,5.23014940361368490D-1, &
         5.14718788772630577D-1,5.06808576734283345D-1,4.99255145565447538D-1, &
         4.92032339458264325D-1,4.85116618392624313D-1,4.78486730436807416D-1, &
         4.72123432945047398D-1,4.66009254246335764D-1,4.60128289044821961D-1, &
         4.54466022030378857D-1,4.49009175209461305D-1,4.43745575272018022D-1 /)
    !
    IFac = 1+(-1)**L
    IF(IFac == 0) THEN
       FT_FScriptC_3D = Zero
    ELSE
       Isgn = (-1)**(L/2)
       Fac  = ExpFac/DBLE(2*L-1)
       FT_FScriptC_3D = Isgn*(Sfac(L)*R*EXP(-Fac*R*R))**(2*L-1)
    ENDIF
    !
  END FUNCTION FT_FScriptC_3D
!========================================================================================
! FT_FSCriptS
!========================================================================================
  FUNCTION FT_FScriptS_3D(L,ExpFac,R)
    INTEGER                    :: L,IFac,Isgn
    REAL(DOUBLE)               :: ExpFac,R,FT_FScriptS_3D,Fac
    REAL(DOUBLE),DIMENSION(41) :: Sfac = (/2.00000000000000000D0 ,1.61199195401646964D0 , &
         1.39399011673806825D0 ,1.24836057075206815D0 ,1.14180213615392840D0 , &
         1.05927697236749128D0 ,9.92824332104422014D-1,9.37767832227969827D-1, &
         8.91149652589504141D-1,8.50992321055495993D-1,8.15915401864279526D-1, &
         7.84921233738737833D-1,7.57267966537096089D-1,7.32390674279586654D-1, &
         7.09850372573553877D-1,6.89299957250757043D-1,6.70460791874183787D-1, &
         6.53106213773379367D-1,6.37049661171763433D-1,6.22135962734544435D-1, &
         6.08234838286273008D-1,5.95235975457893168D-1,5.83045248969098145D-1, &
         5.71581781316674666D-1,5.60775631817133548D-1,5.50565960943842365D-1, &
         5.40899558418448758D-1,5.31729652703953292D-1,5.23014940361368490D-1, &
         5.14718788772630577D-1,5.06808576734283345D-1,4.99255145565447538D-1, &
         4.92032339458264325D-1,4.85116618392624313D-1,4.78486730436807416D-1, &
         4.72123432945047398D-1,4.66009254246335764D-1,4.60128289044821961D-1, &
         4.54466022030378857D-1,4.49009175209461305D-1,4.43745575272018022D-1 /)
    !
    IFac = 1-(-1)**L
    IF(IFac == 0) THEN
       FT_FScriptS_3D = Zero
    ELSE
       Isgn = (-1)**((L-1)/2)
       Fac  = ExpFac/DBLE(2*L-1)
       FT_FScriptS_3D = Isgn*(Sfac(L)*R*EXP(-Fac*R*R))**(2*L-1)
    ENDIF
    !
  END FUNCTION FT_FScriptS_3D
!========================================================================================
!   FT_FSCriptC
!========================================================================================
  FUNCTION GScript(L,R)
    INTEGER                    :: L,LS
    REAL(DOUBLE)               :: R,SqrtR,GScript,XSUM
    REAL(DOUBLE),DIMENSION(41) :: Sfac = (/2.0000000000000000D0 ,1.333333333333333333D0 ,&
         7.3029674334022148D-1,5.3412580601946498D-1,4.2897258944050779D-1,&
         3.6130260767122454D-1,3.1338173978619282D-1,2.7736675692301273D-1,&
         2.4916970717917522D-1,2.2642030439700098D-1,2.0763707534036204D-1,&
         1.9184081786310919D-1,1.7835560611619894D-1,1.6669836456773905D-1,&
         1.5651386447776487D-1,1.4753458765116469D-1,1.3955494670757578D-1,&
         1.3241416731800685D-1,1.2598459485542503D-1,1.2016350807025758D-1,&
         1.1486726370711862D-1,1.1002702835105538D-1,1.0558561445605698D-1,&
         1.0149509929347097D-1,9.7715008595692035D-2,9.4210913822889116D-2,&
         9.0953336661706548D-2,8.7916884656620846D-2,8.5079562763772593D-2,&
         8.2422220248006543D-2,7.9928102738676689D-2,7.7582486742601370D-2,&
         7.5372379364895518D-2,7.3286270006162048D-2,7.1313923796257465D-2,&
         6.9446208774383175D-2,6.7674950532187841D-2,6.5992809342867331D-2,&
         6.4393175806982839D-2,6.2870081828999622D-2,6.1418124351705448D-2 /)
    !
    SqrtR      = SQRT(R)
    IF(L == 0) THEN
       GScript = (One-ERF(SqrtR))
    ELSE
       XSUM = SFAC(1)*EXP(-R)
       DO LS = 2,L
          XSUM  = XSUM+(Sfac(LS)*R*EXP(-R/DBLE(LS-1)))**(LS-1)
       ENDDO
       GScript = (One-ERF(SqrtR))+(SqrtR/SqrtPi)*XSUM
    ENDIF
    !
  END FUNCTION GScript
!========================================================================================
!   FT_FSCriptS
!========================================================================================
  FUNCTION FScript(L,R)
    INTEGER                    :: L,LS
    REAL(DOUBLE)               :: R,SqrtR,FScript,XSUM
    REAL(DOUBLE),DIMENSION(41) :: Sfac = (/2.0000000000000000D0 ,1.333333333333333333D0 ,&
         7.3029674334022148D-1,5.3412580601946498D-1,4.2897258944050779D-1,&
         3.6130260767122454D-1,3.1338173978619282D-1,2.7736675692301273D-1,&
         2.4916970717917522D-1,2.2642030439700098D-1,2.0763707534036204D-1,&
         1.9184081786310919D-1,1.7835560611619894D-1,1.6669836456773905D-1,&
         1.5651386447776487D-1,1.4753458765116469D-1,1.3955494670757578D-1,&
         1.3241416731800685D-1,1.2598459485542503D-1,1.2016350807025758D-1,&
         1.1486726370711862D-1,1.1002702835105538D-1,1.0558561445605698D-1,&
         1.0149509929347097D-1,9.7715008595692035D-2,9.4210913822889116D-2,&
         9.0953336661706548D-2,8.7916884656620846D-2,8.5079562763772593D-2,&
         8.2422220248006543D-2,7.9928102738676689D-2,7.7582486742601370D-2,&
         7.5372379364895518D-2,7.3286270006162048D-2,7.1313923796257465D-2,&
         6.9446208774383175D-2,6.7674950532187841D-2,6.5992809342867331D-2,&
         6.4393175806982839D-2,6.2870081828999622D-2,6.1418124351705448D-2 /)
    !
    SqrtR      = SQRT(R)
    IF(L == 0) THEN
       FScript = ERF(SqrtR)
    ELSE
       XSUM = SFAC(1)*EXP(-R)
       DO LS = 2,L
          XSUM  = XSUM+(Sfac(LS)*R*EXP(-R/DBLE(LS-1)))**(LS-1)
       ENDDO
       FScript = ERF(SqrtR)-(SqrtR/SqrtPi)*XSUM
    ENDIF
    !
  END FUNCTION FScript
!========================================================================================
! Rieman Zeta Function
!========================================================================================
  FUNCTION RZeta(N,M)
    INTEGER                    :: N,M,I
    REAL(DOUBLE)               :: RZeta,RSum
    REAL(DOUBLE),DIMENSION(56) :: RZ = (/ 0.00000000000000000D0, 1.64493406684822644D0, &
         1.20205690315959429D0, 1.08232323371113819D0, 1.03692775514336993D0, &
         1.01734306198444914D0, 1.00834927738192283D0, 1.00407735619794434D0, &
         1.00200839282608221D0, 1.00099457512781809D0, 1.00049418860411946D0, &
         1.00024608655330805D0, 1.00012271334757849D0, 1.00006124813505870D0, &
         1.00003058823630702D0, 1.00001528225940865D0, 1.00000763719763790D0, &
         1.00000381729326500D0, 1.00000190821271655D0, 1.00000095396203387D0, &
         1.00000047693298679D0, 1.00000023845050273D0, 1.00000011921992597D0, &
         1.00000005960818905D0, 1.00000002980350351D0, 1.00000001490155483D0, &
         1.00000000745071179D0, 1.00000000372533402D0, 1.00000000186265972D0, &
         1.00000000093132743D0, 1.00000000046566291D0, 1.00000000023283118D0, &
         1.00000000011641550D0, 1.00000000005820772D0, 1.00000000002910385D0, &
         1.00000000001455192D0, 1.00000000000727596D0, 1.00000000000363798D0, &
         1.00000000000181899D0, 1.00000000000090949D0, 1.00000000000045475D0, &
         1.00000000000022737D0, 1.00000000000011369D0, 1.00000000000005684D0, &
         1.00000000000002842D0, 1.00000000000001421D0, 1.00000000000000711D0, &
         1.00000000000000355D0, 1.00000000000000178D0, 1.00000000000000089D0, &
         1.00000000000000044D0, 1.00000000000000022D0, 1.00000000000000011D0, &
         1.00000000000000006D0, 1.00000000000000003D0, 1.00000000000000001D0  /)

    !
    IF(N .LE. 56) THEN
       RZeta = RZ(N)
    ELSE
       RZeta = One
    ENDIF
    !
    RSum = Zero
    DO I=1,M
       RSum = RSum + One/(DBLE(I)**N)
    ENDDO
    RZeta = RZeta - RSum
    !
  END FUNCTION RZeta










  !========================================================================================
  ! Calculate the PFFTensor 3D
  !========================================================================================
  SUBROUTINE DirectTensor3D(MaxL,GM,Args,CS,TenC,TenS)
    INTEGER                           :: MaxL
    INTEGER                           :: I,J,K,L,M,LM,NC,NCell
    INTEGER                           :: LSwitch
    TYPE(DBL_VECT)                    :: TenC,TenS
    TYPE(CellSet)                     :: CS, CSMM
    REAL(DOUBLE),DIMENSION(3)         :: PQ
    REAL(DOUBLE)                      :: CFac,SFac,BetaSq,Rad,RadSq,ExpFac
    REAL(DOUBLE)                      :: Rmin,RMAX,Accuracy,LenScale,SUM,dum,dum2,eum,eum2
    !
    INTEGER                           :: NDiv,IFF,HardP
    REAL(DOUBLE)                      :: LenMax,Delt,IntFac,X,Y,Z
    REAL(DOUBLE),DIMENSION(3)         :: Vec
    REAL(DOUBLE),DIMENSION(3,3)       :: RecpLatVec,LatVec
    TYPE(CRDS)                        :: GM
    TYPE(ARGMT)                       :: Args
    LOGICAL :: OutCell,InCell

    !-------------------------------------------------------------------------------------!
    !
    !   Initialize
    !
    TenC%D=Zero
    TenS%D=Zero

    IFF=100
    HardP=1

    dum=0D0
    dum2=0D0
    eum=0D0
    eum2=0D0

    MaxL=2

    NCell=0

    DO IFF=100,10000,100

    DO I=-IFF,IFF
       DO J=-IFF,IFF
          DO K=-IFF,IFF

             X=I*GM%PBC%BoxShape%D(1,1)+J*GM%PBC%BoxShape%D(1,2)+K*GM%PBC%BoxShape%D(1,3)
             Y=I*GM%PBC%BoxShape%D(2,1)+J*GM%PBC%BoxShape%D(2,2)+K*GM%PBC%BoxShape%D(2,3)
             Z=I*GM%PBC%BoxShape%D(3,1)+J*GM%PBC%BoxShape%D(3,2)+K*GM%PBC%BoxShape%D(3,3)
!             InCell=ABS(I).EQ.0.AND.ABS(J).EQ.0.AND.ABS(K).EQ.0

             InCell=ABS(K).EQ.0
             OutCell=.NOT.InCell

             IF(OutCell)THEN
                CALL IrRegular(MaxL,X,Y,Z,Print_O=.TRUE.)
                DO L = 2,2 !1,MaxL
                   DO M = 0,0
                      LM = LTD(L)+M
                      TenC%D(LM)=TenC%D(LM)+Cpq(LM)
                      TenS%D(LM)=TenS%D(LM)+Spq(LM)
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    WRITE(*,*)IFF,TENC%D(3)

 ENDDO

    DO L = 1,4
       DO M = 0,L
          LM = LTD(L)+M
          WRITE(*,33)L,M,TenC%D(LM)
33        FORMAT(" L = ",I2," M = ",I2,' Sphere = ',1(D12.6,", "))
       ENDDO
    ENDDO

    STOP
    !
  END SUBROUTINE DirectTensor3D
  !========================================================================================
  ! Calculate the PFFTensor 1D
  !========================================================================================
  SUBROUTINE DirectTensor1D(MaxL,GM,Args,CS,TenC,TenS)
    INTEGER                           :: MaxL
    INTEGER                           :: I,J,K,L,M,LM,NC,NCell
    INTEGER                           :: LSwitch
    TYPE(DBL_VECT)                    :: TenC,TenS
    TYPE(CellSet)                     :: CS, CSMM
    REAL(DOUBLE),DIMENSION(3)         :: PQ
    REAL(DOUBLE)                      :: CFac,SFac,BetaSq,Rad,RadSq,ExpFac
    REAL(DOUBLE)                      :: Rmin,RMAX,Accuracy,LenScale,SUM,dum,dum2,eum,eum2
    !
    INTEGER                           :: NDiv,IFF,HardP
    REAL(DOUBLE)                      :: LenMax,Delt,IntFac,X,Y,Z
    REAL(DOUBLE),DIMENSION(3)         :: Vec
    REAL(DOUBLE),DIMENSION(3,3)       :: RecpLatVec,LatVec
    TYPE(CRDS)                        :: GM
    TYPE(ARGMT)                       :: Args
    LOGICAL :: OutCell,InCell

    !-------------------------------------------------------------------------------------!
    !
    !   Initialize
    !
    TenC%D=Zero
    TenS%D=Zero
    IFF=4
    DO I=-IFF,IFF
       J=0
       K=0
       X=I*GM%PBC%BoxShape%D(1,1)+J*GM%PBC%BoxShape%D(1,2)+K*GM%PBC%BoxShape%D(1,3)
       Y=I*GM%PBC%BoxShape%D(2,1)+J*GM%PBC%BoxShape%D(2,2)+K*GM%PBC%BoxShape%D(2,3)
       Z=I*GM%PBC%BoxShape%D(3,1)+J*GM%PBC%BoxShape%D(3,2)+K*GM%PBC%BoxShape%D(3,3)
       IF(ABS(I)>1)THEN
          CALL IrRegular(MaxL,X,Y,Z,Print_O=.TRUE.)
          DO L = 0,MaxL
             DO M = 0,L
                LM = LTD(L)+M
                TenC%D(LM)=TenC%D(LM)+Cpq(LM)
                TenS%D(LM)=TenS%D(LM)+Spq(LM)
             ENDDO
          ENDDO
       ENDIF
    ENDDO
  END SUBROUTINE DirectTensor1D
!
END MODULE PFFTen
