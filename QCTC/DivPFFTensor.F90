MODULE DivPFFTen
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
  USE PFFTen
  IMPLICIT NONE
  ! Globals
CONTAINS
!========================================================================================
!
!========================================================================================
  SUBROUTINE CalculateDivPFFT(MaxEll,GM,Args,CS,dTenC,dTenS)
    TYPE(CRDS)             :: GM
    TYPE(ARGMT)            :: Args
    TYPE(DBL_RNK3)         :: dTenC,dTenS
    TYPE(CellSet)          :: CS
    INTEGER                :: MaxEll
!-------------------------------------------------------------------------------------!
    IF(GM%PBC%Dimen==1) THEN
       CALL MakeDivTensor1D(2*MaxEll,GM,Args,CS,dTenC,dTenS)
    ELSEIF(GM%PBC%Dimen==2) THEN
       CALL MakeDivTensor2D(2*MaxEll,GM,Args,CS,dTenC,dTenS)
    ELSEIF(GM%PBC%Dimen==3) THEN
       CALL MakeDivTensor3D(2*MaxEll,GM,Args,CS,dTenC,dTenS)
    ENDIF
!
  END SUBROUTINE CalculateDivPFFT
!========================================================================================
! Calculate the Derivative PFFTensor 1D
!========================================================================================
  SUBROUTINE MakeDivTensor1D(MaxL,GM,Args,CS,dTenC,dTenS)
    INTEGER                           :: MaxL,LenL
    INTEGER                           :: I,J,L,M,LM,NC,IJ
    TYPE(CellSet)                     :: CS, CSMM
    TYPE(DBL_RNK3)                    :: dTenC,dTenS
    TYPE(CRDS)                        :: GM
    TYPE(ARGMT)                       :: Args
!-------------------------------------------------------------------------------------!
!
!   Number of Inner Boxes
!
    LenL=LSP(MaxL)
    dTenC%D=Zero
    dTenS%D=Zero
    NC = (CS%NCells-1)/2
!
!   One Dimension: Right and Left
!
    IF(GM%PBC%AutoW%I(1)==1) THEN
       IJ = 1
       CALL DIrRegular(MaxL,GM%PBC%BoxShape%D(IJ,IJ),Zero,Zero)
       dTenC%D(0:LenL,IJ,IJ) = -DCpq(0:LenL,IJ)
       CALL DIrRegular(MaxL,-GM%PBC%BoxShape%D(IJ,IJ),Zero,Zero)
       dTenC%D(0:LenL,IJ,IJ) = dTenC%D(0:LenL,IJ,IJ)+DCpq(0:LenL,IJ)
    ELSEIF(GM%PBC%AutoW%I(2)==1) THEN
       IJ = 2
       CALL DIrRegular(MaxL,Zero, GM%PBC%BoxShape%D(IJ,IJ),Zero)
       dTenC%D(0:LenL,IJ,IJ) = -DCpq(0:LenL,IJ)
       CALL DIrRegular(MaxL,Zero,-GM%PBC%BoxShape%D(IJ,IJ),Zero)
       dTenC%D(0:LenL,IJ,IJ) = dTenC%D(0:LenL,IJ,IJ)+DCpq(0:LenL,IJ)
    ELSEIF(GM%PBC%AutoW%I(3)==1) THEN
       IJ = 3
       CALL DIrRegular(MaxL,Zero,Zero, GM%PBC%BoxShape%D(IJ,IJ))
       dTenC%D(0:LenL,IJ,IJ) = -DCpq(0:LenL,IJ)
       CALL DIrRegular(MaxL,Zero,Zero,-GM%PBC%BoxShape%D(IJ,IJ))
       dTenC%D(0:LenL,IJ,IJ) = dTenC%D(0:LenL,IJ,IJ)+DCpq(0:LenL,IJ)
    ENDIF
    DO L=1,MaxL
       DO M = 0,L
          LM = LTD(L)+M
          dTenC%D(LM,IJ,IJ) = dTenC%D(LM,IJ,IJ)*RZeta(L+1,NC)
       ENDDO
    ENDDO
!
  END SUBROUTINE MakeDivTensor1D
!========================================================================================
! Calculate the PFFTensor 2D
!========================================================================================
  SUBROUTINE MakeDivTensor2D(MaxL,GM,Args,CS,dTenC,dTenS)
    INTEGER                           :: MaxL
    INTEGER                           :: I,J,K,L,M,LM,NC,ND
    TYPE(DBL_RNK3)                    :: dTenC,dTenS
    TYPE(CellSet)                     :: CS, CSMM
    REAL(DOUBLE),DIMENSION(3)         :: PQ,FPQ
    REAL(DOUBLE)                      :: CFac,SFac,BetaSq,Rad,RadSq,ExpFac
    REAL(DOUBLE)                      :: Rmin,RMAX,Accuracy,LenScale,SUM,FacC,FacS
!
    INTEGER                           :: NDiv
    REAL(DOUBLE)                      :: LenMax,Delt,IntFac,MCFac,MSFac
!
    REAL(DOUBLE),DIMENSION(3)         :: Vec,DCFac,DSFac
    REAL(DOUBLE),DIMENSION(3,3)       :: RecpLatVec,LatVec,DivVol
    TYPE(CRDS)                        :: GM
    TYPE(ARGMT)                       :: Args
!-------------------------------------------------------------------------------------!
!
!   Initialize
!
    dTenC%D=Zero
    dTenS%D=Zero

    ! Oh god, this routine is total crap.  Need to implement analytic derivatives
    ! based on plane wise summation.  Maybe someone will do that to avoid total
    ! humiliation.

    Accuracy=1D-12


    Rmin     = SQRT(CS%CellCarts%D(1,1)**2+CS%CellCarts%D(2,1)**2+CS%CellCarts%D(3,1)**2)
    DivVol = DivCellVolume(GM%PBC%BoxShape%D,GM%PBC%AutoW%I)
!
!   Get The Lattice and Reciprocal Lattice Vectors
!
    DO I = 1,3
       DO J = 1,3
          RecpLatVec(I,J) = GM%PBC%InvBoxSh%D(J,I)
          LatVec(I,J)     = GM%PBC%BoxShape%D(I,J)
       ENDDO
    ENDDO
!
!   Two Dimension
!
    LenScale = 6.0D0!GM%PBC%CellVolume**(Half)
    Rmax     = Rmin+LenScale*(One/Accuracy)**(One/DBLE(LSwitch))
    BetaSq   = One/(LenScale)**2
!
!   Sum the Real Space
!
    DO
       CALL New_CellSet_Sphere(CSMM,GM%PBC%AutoW%I,LatVec,Rmax)
       IF(CSMM%NCells .GT. 50000) THEN
          Rmax = 0.99*Rmax
          CALL Delete_CellSet(CSMM)
       ELSE
          EXIT
       ENDIF
    ENDDO
    CALL Sort_CellSet(CSMM)
!   Sum over Cells
    DO NC = 1,CSMM%NCells
       PQ(:)  =  CSMM%CellCarts%D(:,NC)
       FPQ(:) = -AtomToFrac(GM,PQ(:))
       IF(.NOT. InCell_CellSet(CS,PQ(1),PQ(2),PQ(3))) THEN
          RadSq = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
          CALL DIrRegular(MaxL,PQ(1),PQ(2),PQ(3))
          DO L = 1,MaxL
             IF(L .LE. LSwitch) THEN
                CFac  =  GScript(L,RadSq)
                DCFac = DGScript(L,BetaSq,PQ(1),PQ(2),PQ(3))
             ELSE
                CFac  = One
                DCFac = Zero
             ENDIF
             DO M = 0,L
                LM = LTD(L)+M
                DO I=1,3
                   DO J=1,3
                      IF(GM%PBC%AutoW%I(I)==1 .AND. GM%PBC%AutoW%I(J)==1) THEN
                         dTenC%D(LM,I,J)=dTenC%D(LM,I,J)+FPQ(J)*(CFac*DCpq(LM,I) + DCFac(I)*Cpq(LM))
                         dTenS%D(LM,I,J)=dTenS%D(LM,I,J)+FPQ(J)*(CFac*DSpq(LM,I) + DCFac(I)*Spq(LM))
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    CALL Delete_CellSet(CSMM)
!
!   Sum the Reciprical Space
!
    ExpFac = Pi*Pi/BetaSq
    Rmax = SQRT(ABS(LOG(Accuracy/(10.D0**(LSwitch)))/ExpFac))


    DO
       CALL New_CellSet_Sphere(CSMM,GM%PBC%AutoW%I,RecpLatVec,Rmax)
       IF(CSMM%NCells .LT. 9) THEN
          Rmax = 1.01D0*Rmax
          CALL Delete_CellSet(CSMM)
       ELSE
          EXIT
       ENDIF
    ENDDO
    CALL Sort_CellSet(CSMM)
!
    NDiv    = 1000
    LenMax  = SQRT(ABS(LOG(1.D-6))/ExpFac)
    Delt    = LenMax/DBLE(NDiv)
    DO ND=-NDiv,NDiv
       Vec(:) = Zero
       IF(GM%PBC%AutoW%I(1)==0) Vec(1) = ND*Delt
       IF(GM%PBC%AutoW%I(2)==0) Vec(2) = ND*Delt
       IF(GM%PBC%AutoW%I(3)==0) Vec(3) = ND*Delt
       DO NC = 1,CSMM%NCells
          PQ(:) = CSMM%CellCarts%D(:,NC)+Vec(:)
          Rad   = SQRT(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
          IF(Rad .GT. 1.D-14) THEN
             CALL DIrRegular(MaxL,PQ(1),PQ(2),PQ(3))
             DO L = 1,LSwitch
                CFac  = Delt*FT_FScriptC_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
                SFac  = Delt*FT_FScriptS_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
                DCFac = Delt*DFT_FScriptC_3D(L,ExpFac,PQ(1),PQ(2),PQ(3))/GM%PBC%CellVolume
                DSFac = Delt*DFT_FScriptS_3D(L,ExpFac,PQ(1),PQ(2),PQ(3))/GM%PBC%CellVolume
                DO M = 0,L
                   LM = LTD(L)+M
                   DO I=1,3
                      DO J=1,3
                         IF(GM%PBC%AutoW%I(I)==1 .AND. GM%PBC%AutoW%I(J)==1) THEN
                            DO K=1,3
                                  dTenC%D(LM,I,J)=dTenC%D(LM,I,J) + PQ(I)*RecpLatVec(K,J)*(CFac*DCpq(LM,K)+ DCFac(K)*Cpq(LM)) &
                                                                  - PQ(I)*RecpLatVec(K,J)*(SFac*DSpq(LM,K)+ DSFac(K)*Spq(LM))
                                  dTenS%D(LM,I,J)=dTenS%D(LM,I,J) + PQ(I)*RecpLatVec(K,J)*(SFac*DCpq(LM,K)+ DSFac(K)*Cpq(LM)) &
                                                                  + PQ(I)*RecpLatVec(K,J)*(CFac*DSpq(LM,K)+ DCFac(K)*Spq(LM))
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                   MCFac = (Cpq(LM)*CFac-Spq(LM)*SFac)/GM%PBC%CellVolume
                   MSFac = (Spq(LM)*CFac+Cpq(LM)*SFac)/GM%PBC%CellVolume
                   DO I=1,3
                      DO J=1,3
                         IF(GM%PBC%AutoW%I(I)==1 .AND. GM%PBC%AutoW%I(J)==1) THEN
                            dTenC%D(LM,I,J)=dTenC%D(LM,I,J) + DivVol(I,J)*MCFac
                            dTenS%D(LM,I,J)=dTenS%D(LM,I,J) + DivVol(I,J)*MSFac
                         ENDIF
                      ENDDO
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    CALL Delete_CellSet(CSMM)
!
!   Add in the k=0 piece for the L=2 multipoles (Zero in 3D)
!
    DO I=1,3
       IF(GM%PBC%AutoW%I(I)==0) THEN
          PQ(:) = Zero
          PQ(I) = 1.D-10
          Rad   = SQRT(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
          L     = 2
          CALL DIrRegular(L,PQ(1),PQ(2),PQ(3))
          CFac  = Half*Delt*FT_FScriptC_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
          SFac  = Half*Delt*FT_FScriptS_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
          DCFac = Half*Delt*DFT_FScriptC_3D(L,ExpFac,PQ(1),PQ(2),PQ(3))/GM%PBC%CellVolume
          DSFac = Half*Delt*DFT_FScriptS_3D(L,ExpFac,PQ(1),PQ(2),PQ(3))/GM%PBC%CellVolume
          DO M = 0,L
             LM = LTD(L)+M
             MCFac = (Cpq(LM)*CFac-Spq(LM)*SFac)/GM%PBC%CellVolume
             MSFac = (Spq(LM)*CFac+Cpq(LM)*SFac)/GM%PBC%CellVolume
             DO J=1,3
                DO K=1,3
                   IF(GM%PBC%AutoW%I(J)==1 .AND. GM%PBC%AutoW%I(K)==1) THEN
                      dTenC%D(LM,J,K)=dTenC%D(LM,J,K) + DivVol(J,K)*MCFac
                      dTenS%D(LM,J,K)=dTenS%D(LM,J,K) + DivVol(J,K)*MSFac
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
          PQ(:) = Zero
          PQ(I) = -1.D-10
          Rad   = SQRT(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
          L     = 2
          CALL DIrRegular(L,PQ(1),PQ(2),PQ(3))
          CFac  = Half*Delt*FT_FScriptC_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
          SFac  = Half*Delt*FT_FScriptS_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
          DCFac = Half*Delt*DFT_FScriptC_3D(L,ExpFac,PQ(1),PQ(2),PQ(3))/GM%PBC%CellVolume
          DSFac = Half*Delt*DFT_FScriptS_3D(L,ExpFac,PQ(1),PQ(2),PQ(3))/GM%PBC%CellVolume
          DO M = 0,L
             LM = LTD(L)+M
             MCFac = (Cpq(LM)*CFac-Spq(LM)*SFac)/GM%PBC%CellVolume
             MSFac = (Spq(LM)*CFac+Cpq(LM)*SFac)/GM%PBC%CellVolume
             DO J=1,3
                DO K=1,3
                   IF(GM%PBC%AutoW%I(J)==1 .AND. GM%PBC%AutoW%I(K)==1) THEN
                      dTenC%D(LM,J,K)=dTenC%D(LM,J,K) + DivVol(J,K)*MCFac
                      dTenS%D(LM,J,K)=dTenS%D(LM,J,K) + DivVol(J,K)*MSFac
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    CALL Delete_CellSet(CSMM)
!
!   Substract the inner boxes
!
    DO NC = 1,CS%NCells
       PQ(:)  =  CS%CellCarts%D(:,NC)
       FPQ(:) = -AtomToFrac(GM,PQ(:))
       RadSq  = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
       IF(RadSq .GT. 1.D-14) THEN
          CALL DIrRegular(MaxL,PQ(1),PQ(2),PQ(3))
          DO L = 1,LSwitch
             CFac  =  FScript(L,RadSq)
             DCFac = DFScript(L,BetaSq,PQ(1),PQ(2),PQ(3))
             DO M = 0,L
                LM = LTD(L)+M
                DO I=1,3
                   DO J=1,3
                      IF(GM%PBC%AutoW%I(I)==1 .AND. GM%PBC%AutoW%I(J)==1) THEN
                         dTenC%D(LM,I,J)=dTenC%D(LM,I,J)-FPQ(J)*(CFac*DCpq(LM,I) + DCFac(I)*Cpq(LM))
                         dTenS%D(LM,I,J)=dTenS%D(LM,I,J)-FPQ(J)*(CFac*DSpq(LM,I) + DCFac(I)*Spq(LM))
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDDO
!
  END SUBROUTINE MakeDivTensor2D

!========================================================================================
! Calculate the Derivative PFFTensor 3D
!========================================================================================
  SUBROUTINE MakeDivTensor3D(MaxL,GM,Args,CS,dTenC,dTenS)
    INTEGER                           :: MaxL,NBDWSummedL
    INTEGER                           :: I,J,K,L,M,LM,NC
    TYPE(DBL_RNK3)                    :: dTenC,dTenS
    TYPE(CellSet)                     :: CS, CSMM
    REAL(DOUBLE),DIMENSION(3)         :: PQ,FPQ
    REAL(DOUBLE)                      :: CFac,SFac,Rad,RadSq,ExpFac
    REAL(DOUBLE)                      :: Beta,BetaSq,R,RMin,KMin,RMax,KMax
    REAL(DOUBLE)                      :: FacC,FacS,MCFac,MSFac
    REAL(DOUBLE),DIMENSION(3)         :: Vec,DCFac,DSFac
    REAL(DOUBLE),DIMENSION(3,3)       :: RecpLatVec,LatVec,DivVol
    TYPE(CRDS)                        :: GM
    TYPE(ARGMT)                       :: Args
    COMPLEX(DOUBLE),DIMENSION(0:PlaneWiseEll,0:PlaneWiseEll,3,3) :: dYlm
    !-------------------------------------------------------------------------------------!
    ! Initialize
    !
    dTenC%D=Zero
    dTenS%D=Zero
    DO I = 1,3
       DO J = 1,3
          RecpLatVec(I,J) = GM%PBC%InvBoxSh%D(J,I)
          LatVec(I,J)     = GM%PBC%BoxShape%D(I,J)
       ENDDO
    ENDDO
    DivVol = DivCellVolume(GM%PBC%BoxShape%D,GM%PBC%AutoW%I)
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
    BetaSq=Beta**2
    !-----------------------------------------------------------------------
    ! This is the highest angular symmetry we will consider for NBDW summation.
    ! If MaxL<LSwitch, we only go as high as MaxL.
    NBDWSummedL=MIN(LSwitch,MaxL)
    !-----------------------------------------------------------------------
    ! OK, this bit is important.  RMax is the distance the real real space
    ! lattice sum is carried out over.  HOWEVER, for large L, direct summation
    ! is used, and this MAY BE LARGER than the radius given by the balanced
    ! value.  Here are the balanced values:
    RMax=(One/Beta)*SQRT(ABS(LOG(Accuracy)))
    KMax=    Beta *SQRT(ABS(LOG(Accuracy)))/Pi
    ! And here is the value we use, which takes into account direct summation
    ! in the first loop below (CFac=1), which typically involves more cells than
    ! nessecary in the balanced NBDW summation.
    RMax=MAX(RMax,Accuracy**(-One/DBLE(NBDWSummedL+1)))
    !---------------------------------------------------------------------------------
    ! REAL SPACE
    !---------------------------------------------------------------------------------
    CALL SetCellSets(CSMM,GM%PBC%AutoW%I,GM%PBC%BoxShape%D,'DivTensor3d[R]','Radial',RMin_O=RMax)
    !
    DO NC = 1,CSMM%NCells-1 ! Leaves out central cell (CellSet is sorted, so PQ=0 is last)
       PQ(:)  =  CSMM%CellCarts%D(:,NC)
       FPQ(:) = -AtomToFrac(GM,PQ(:))
       IF(.NOT. InCell_CellSet(CS,PQ(1),PQ(2),PQ(3))) THEN
          RadSq = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
          CALL DIrRegular(MaxL,PQ(1),PQ(2),PQ(3))
          DO L = 1,MaxL
             IF(L.LE.NBDWSummedL) THEN
                CFac     = GScript(L,RadSq)
                DCFac(:) = DGScript(L,BetaSq,PQ(1),PQ(2),PQ(3))
             ELSE
                CFac     = One
                DCFac(:) = Zero
             ENDIF
             DO M = 0,L
                LM = LTD(L)+M
                DO I=1,3
                   DO J=1,3
                      dTenC%D(LM,I,J)=dTenC%D(LM,I,J)+FPQ(J)*(CFac*DCpq(LM,I) + DCFac(I)*Cpq(LM))
                      dTenS%D(LM,I,J)=dTenS%D(LM,I,J)+FPQ(J)*(CFac*DSpq(LM,I) + DCFac(I)*Spq(LM))
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    CALL Delete_CellSet(CSMM)
    !---------------------------------------------------------------------------------
    ! RECIPROCAL SPACE
    !---------------------------------------------------------------------------------
    ExpFac = Pi*Pi/BetaSq
    CALL SetCellSets(CSMM,GM%PBC%AutoW%I,RecpLatVec,'DivTensor3D[K]','Radial',RMin_O=KMax)
    !
    DO NC = 1,CSMM%NCells-1
       PQ(:)  = CSMM%CellCarts%D(:,NC)
       Rad    = SQRT(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
       CALL DIrRegular(NBDWSummedL,PQ(1),PQ(2),PQ(3))
       DO L = 1,LSwitch
          CFac  = FT_FScriptC_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
          SFac  = FT_FScriptS_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
          DCFac = DFT_FScriptC_3D(L,ExpFac,PQ(1),PQ(2),PQ(3))/GM%PBC%CellVolume
          DSFac = DFT_FScriptS_3D(L,ExpFac,PQ(1),PQ(2),PQ(3))/GM%PBC%CellVolume
          DO M = 0,L
             LM = LTD(L)+M
             DO I=1,3
                DO J=1,3
                   DO K=1,3
                      dTenC%D(LM,I,J)=dTenC%D(LM,I,J) + PQ(I)*RecpLatVec(K,J)*(CFac*DCpq(LM,K)+ DCFac(K)*Cpq(LM)) &
                           - PQ(I)*RecpLatVec(K,J)*(SFac*DSpq(LM,K)+ DSFac(K)*Spq(LM))
                      dTenS%D(LM,I,J)=dTenS%D(LM,I,J) + PQ(I)*RecpLatVec(K,J)*(SFac*DCpq(LM,K)+ DSFac(K)*Cpq(LM)) &
                           + PQ(I)*RecpLatVec(K,J)*(CFac*DSpq(LM,K)+ DCFac(K)*Spq(LM))
                   ENDDO
                ENDDO
             ENDDO
             MCFac = (Cpq(LM)*CFac-Spq(LM)*SFac)/GM%PBC%CellVolume
             MSFac = (Spq(LM)*CFac+Cpq(LM)*SFac)/GM%PBC%CellVolume
             DO I=1,3
                DO J=1,3
                   dTenC%D(LM,I,J)=dTenC%D(LM,I,J) + DivVol(I,J)*MCFac
                   dTenS%D(LM,I,J)=dTenS%D(LM,I,J) + DivVol(I,J)*MSFac
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    CALL Delete_CellSet(CSMM)
    !---------------------------------------------------------------------------------
    ! SUBTRACT THE NEAR FIELD
    !---------------------------------------------------------------------------------
    DO NC = 1,CS%NCells-1
       PQ(:)  =  CS%CellCarts%D(:,NC)
       FPQ(:) = -AtomToFrac(GM,PQ(:))
       RadSq  = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
       CALL DIrRegular(LSwitch,PQ(1),PQ(2),PQ(3))
       DO L = 1,NBDWSummedL
!!$          CFac  = FScript(L,RadSq)
          CFac=RegularizedGammaP(DBLE(L)+5D-1,RadSq)
          DCFac = DFScript(L,BetaSq,PQ(1),PQ(2),PQ(3))
          DO M = 0,L
             LM = LTD(L)+M
             DO I=1,3
                DO J=1,3
                   dTenC%D(LM,I,J)=dTenC%D(LM,I,J)-FPQ(J)*(CFac*DCpq(LM,I) + DCFac(I)*Cpq(LM))
                   dTenS%D(LM,I,J)=dTenS%D(LM,I,J)-FPQ(J)*(CFac*DSpq(LM,I) + DCFac(I)*Spq(LM))
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
!
  END SUBROUTINE MakeDivTensor3D
!!$
!!$
!!$!========================================================================================
!!$! Calculate the Derivative PFFTensor 3D
!!$!========================================================================================
!!$  SUBROUTINE MakeDivTensor3D(MaxL,GM,Args,CS,dTenC,dTenS)
!!$    INTEGER                           :: MaxL
!!$    INTEGER                           :: I,J,K,L,M,LM,NC
!!$    INTEGER                           :: LSwitch
!!$    TYPE(DBL_RNK3)                    :: dTenC,dTenS
!!$    TYPE(CellSet)                     :: CS, CSMM
!!$    REAL(DOUBLE),DIMENSION(3)         :: PQ,FPQ
!!$    REAL(DOUBLE)                      :: CFac,SFac,beta,BetaSq,Rad,RadSq,ExpFac
!!$    REAL(DOUBLE)                      :: Rmin,RMAX,Accuracy,LenScale,FacC,FacS,MCFac,MSFac
!!$    REAL(DOUBLE),DIMENSION(3)         :: Vec,DCFac,DSFac
!!$    REAL(DOUBLE),DIMENSION(3,3)       :: RecpLatVec,LatVec,DivVol
!!$    TYPE(CRDS)                        :: GM
!!$    TYPE(ARGMT)                       :: Args
!!$!-------------------------------------------------------------------------------------!
!!$!
!!$!   Initialize
!!$!
!!$    dTenC%D=Zero
!!$    dTenS%D=Zero
!!$    Accuracy = 1.D-12
!!$    Rmin     = SQRT(CS%CellCarts%D(1,1)**2+CS%CellCarts%D(2,1)**2+CS%CellCarts%D(3,1)**2)
!!$    DivVol = DivCellVolume(GM%PBC%BoxShape%D,GM%PBC%AutoW%I)
!!$!
!!$!   Get The Lattice and Reciprocal Lattice Vectors
!!$!
!!$    DO I = 1,3
!!$       DO J = 1,3
!!$          RecpLatVec(I,J) = GM%PBC%InvBoxSh%D(J,I)
!!$          LatVec(I,J)     = GM%PBC%BoxShape%D(I,J)
!!$       ENDDO
!!$    ENDDO
!!$!
!!$!   Three Dimension
!!$!
!!$    LSwitch  = 10
!!$    LenScale = GM%PBC%CellVolume**(One/Three)
!!$    Rmax     = Rmin+LenScale*(One/Accuracy)**(One/DBLE(LSwitch))
!!$    BetaSq  = One/(LenScale)**2
!!$    LSwitch  = MIN(MaxL,LSwitch)
!!$!
!!$!   Sum the Real Space
!!$!
!!$    DO
!!$       CALL New_CellSet_Sphere(CSMM,GM%PBC%AutoW%I,LatVec,Rmax)
!!$       IF(CSMM%NCells .GT. 600000) THEN
!!$          Rmax = 0.99*Rmax
!!$          CALL Delete_CellSet(CSMM)
!!$       ELSE
!!$          EXIT
!!$       ENDIF
!!$    ENDDO
!!$    CALL Sort_CellSet(CSMM)
!!$!   Sum over Cells
!!$    DO NC = 1,CSMM%NCells
!!$       PQ(:)  =  CSMM%CellCarts%D(:,NC)
!!$       FPQ(:) = -AtomToFrac(GM,PQ(:))
!!$       IF(.NOT. InCell_CellSet(CS,PQ(1),PQ(2),PQ(3))) THEN
!!$          RadSq = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
!!$          CALL DIrRegular(MaxL,PQ(1),PQ(2),PQ(3))
!!$          DO L = 1,MaxL
!!$             IF(L .LE. LSwitch) THEN
!!$                CFac     = GScript(L,RadSq)
!!$                DCFac(:) = DGScript(L,BetaSq,PQ(1),PQ(2),PQ(3))
!!$             ELSE
!!$                CFac     = One
!!$                DCFac(:) = Zero
!!$             ENDIF
!!$             DO M = 0,L
!!$                LM = LTD(L)+M
!!$                DO I=1,3
!!$                   DO J=1,3
!!$                      dTenC%D(LM,I,J)=dTenC%D(LM,I,J)+FPQ(J)*(CFac*DCpq(LM,I) + DCFac(I)*Cpq(LM))
!!$                      dTenS%D(LM,I,J)=dTenS%D(LM,I,J)+FPQ(J)*(CFac*DSpq(LM,I) + DCFac(I)*Spq(LM))
!!$                   ENDDO
!!$                ENDDO
!!$
!!$             ENDDO
!!$          ENDDO
!!$       ENDIF
!!$    ENDDO
!!$    CALL Delete_CellSet(CSMM)
!!$!
!!$!   Sum the Reciprical Space
!!$!
!!$    ExpFac = Pi*Pi/BetaSq
!!$    Rmax   = SQRT(ABS(LOG(Accuracy/(10.D0**(LSwitch)))/ExpFac))
!!$    DO
!!$       CALL New_CellSet_Sphere(CSMM,GM%PBC%AutoW%I,RecpLatVec,Rmax)
!!$       IF(CSMM%NCells .LT. 27) THEN
!!$          Rmax = 1.01D0*Rmax
!!$          CALL Delete_CellSet(CSMM)
!!$       ELSE
!!$          EXIT
!!$       ENDIF
!!$    ENDDO
!!$    CALL Sort_CellSet(CSMM)
!!$!
!!$    DO NC = 1,CSMM%NCells
!!$       PQ(:)  = CSMM%CellCarts%D(:,NC)
!!$       Rad    = SQRT(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
!!$       IF(Rad .GT. 1.D-14) THEN
!!$          CALL DIrRegular(MaxL,PQ(1),PQ(2),PQ(3))
!!$          DO L = 1,LSwitch
!!$             CFac  = FT_FScriptC_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
!!$             SFac  = FT_FScriptS_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
!!$             DCFac = DFT_FScriptC_3D(L,ExpFac,PQ(1),PQ(2),PQ(3))/GM%PBC%CellVolume
!!$             DSFac = DFT_FScriptS_3D(L,ExpFac,PQ(1),PQ(2),PQ(3))/GM%PBC%CellVolume
!!$             DO M = 0,L
!!$                LM = LTD(L)+M
!!$                DO I=1,3
!!$                   DO J=1,3
!!$                      DO K=1,3
!!$                         dTenC%D(LM,I,J)=dTenC%D(LM,I,J) + PQ(I)*RecpLatVec(K,J)*(CFac*DCpq(LM,K)+ DCFac(K)*Cpq(LM)) &
!!$                                                         - PQ(I)*RecpLatVec(K,J)*(SFac*DSpq(LM,K)+ DSFac(K)*Spq(LM))
!!$                         dTenS%D(LM,I,J)=dTenS%D(LM,I,J) + PQ(I)*RecpLatVec(K,J)*(SFac*DCpq(LM,K)+ DSFac(K)*Cpq(LM)) &
!!$                                                         + PQ(I)*RecpLatVec(K,J)*(CFac*DSpq(LM,K)+ DCFac(K)*Spq(LM))
!!$                      ENDDO
!!$                   ENDDO
!!$                ENDDO
!!$                MCFac = (Cpq(LM)*CFac-Spq(LM)*SFac)/GM%PBC%CellVolume
!!$                MSFac = (Spq(LM)*CFac+Cpq(LM)*SFac)/GM%PBC%CellVolume
!!$                DO I=1,3
!!$                   DO J=1,3
!!$                      dTenC%D(LM,I,J)=dTenC%D(LM,I,J) + DivVol(I,J)*MCFac
!!$                      dTenS%D(LM,I,J)=dTenS%D(LM,I,J) + DivVol(I,J)*MSFac
!!$                   ENDDO
!!$                ENDDO
!!$             ENDDO
!!$          ENDDO
!!$       ENDIF
!!$    ENDDO
!!$    CALL Delete_CellSet(CSMM)
!!$!
!!$!   Substract the inner boxes
!!$!
!!$    DO NC = 1,CS%NCells
!!$       PQ(:)  =  CS%CellCarts%D(:,NC)
!!$       FPQ(:) = -AtomToFrac(GM,PQ(:))
!!$       RadSq  = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
!!$       IF(RadSq .GT. 1.D-14) THEN
!!$          CALL DIrRegular(MaxL,PQ(1),PQ(2),PQ(3))
!!$          DO L = 1,LSwitch
!!$             CFac  =  FScript(L,RadSq)
!!$             DCFac = DFScript(L,BetaSq,PQ(1),PQ(2),PQ(3))
!!$             DO M = 0,L
!!$                LM = LTD(L)+M
!!$                DO I=1,3
!!$                   DO J=1,3
!!$                      dTenC%D(LM,I,J)=dTenC%D(LM,I,J)-FPQ(J)*(CFac*DCpq(LM,I) + DCFac(I)*Cpq(LM))
!!$                      dTenS%D(LM,I,J)=dTenS%D(LM,I,J)-FPQ(J)*(CFac*DSpq(LM,I) + DCFac(I)*Spq(LM))
!!$                   ENDDO
!!$                ENDDO
!!$             ENDDO
!!$          ENDDO
!!$       ENDIF
!!$    ENDDO
!!$!
!!$  END SUBROUTINE MakeDivTensor3D
!========================================================================================
! Dirivative of the Irregular Fuction
!========================================================================================
  SUBROUTINE DIrRegular_old(MaxL,Px,Py,Pz)
    INTEGER                    :: MaxL,LMtot
    REAL(DOUBLE)               :: Px,Py,Pz
    REAL(DOUBLE)               :: DDelta,Factor,Rad
!
    DDelta = 1.D-6
    Factor = Half/DDelta
    LMtot=LSP(MaxL)
!   x div
    CALL IrRegular(MaxL,Px+DDelta,Py,Pz)
    DCpq(0:LMtot,1)=Cpq(0:LMtot)
    DSpq(0:LMtot,1)=Spq(0:LMtot)
    CALL IrRegular(MaxL,Px-DDelta,Py,Pz)
    DCpq(0:LMtot,1)=(DCpq(0:LMtot,1)-Cpq(0:LMtot))*Factor
    DSpq(0:LMtot,1)=(DSpq(0:LMtot,1)-Spq(0:LMtot))*Factor
!   y div
    CALL IrRegular(MaxL,Px,Py+DDelta,Pz)
    DCpq(0:LMtot,2)=Cpq(0:LMtot)
    DSpq(0:LMtot,2)=Spq(0:LMtot)
    CALL IrRegular(MaxL,Px,Py-DDelta,Pz)
    DCpq(0:LMtot,2)=(DCpq(0:LMtot,2)-Cpq(0:LMtot))*Factor
    DSpq(0:LMtot,2)=(DSpq(0:LMtot,2)-Spq(0:LMtot))*Factor
!   z div
    CALL IrRegular(MaxL,Px,Py,Pz+DDelta)
    DCpq(0:LMtot,3)=Cpq(0:LMtot)
    DSpq(0:LMtot,3)=Spq(0:LMtot)
    CALL IrRegular(MaxL,Px,Py,Pz-DDelta)
    DCpq(0:LMtot,3)=(DCpq(0:LMtot,3)-Cpq(0:LMtot))*Factor
    DSpq(0:LMtot,3)=(DSpq(0:LMtot,3)-Spq(0:LMtot))*Factor
!
    CALL IrRegular(MaxL,Px,Py,Pz)
!
  END SUBROUTINE DIrRegular_old
!========================================================================================
! Dirivative of the Irregular Fuction
!========================================================================================
  SUBROUTINE DIrRegular(MaxL,Px,Py,Pz)
    INTEGER                    :: MaxL,LM,L,M,LMP
    REAL(DOUBLE)               :: Px,Py,Pz,OvRad,OvRad2,OvRad3,XY,OvXY,OvXY2
    REAL(DOUBLE)               :: DrDx,DrDy,DrDz
    REAL(DOUBLE)               :: DphiDx,DphiDy
    REAL(DOUBLE)               :: DthetaDx,DthetaDy,DthetaDz
    REAL(DOUBLE)               :: COTT
    REAL(DOUBLE)               :: FrX,FrY,FrZ,FthetaX1,FthetaY1,FthetaZ1
    REAL(DOUBLE)               :: FthetaX2,FthetaY2,FthetaZ2
!
    CALL IrRegular(MaxL,Px,Py,Pz)
!
    IF(SQRT((Px*Px+Py*Py)) < 1.D-12) THEN
       DrDz     = One/Pz
       DCpq(1:LSP(MaxL),1:3) = Zero
       DSpq(1:LSP(MaxL),1:3) = Zero
       DO L=0,MaxL
          FrZ      = -DBLE(L+1)*DrDz
          LM  = LTD(L)
          IF(L .NE. 0) THEN
!            x-derivative
             DCpq(LM+1,1) = Zero
             DSpq(LM+1,1) = Half*FrZ*Cpq(LM)
!            y-derivative
             DCpq(LM+1,2) = Half*FrZ*Cpq(LM)
             DSpq(LM+1,2) = Zero
          ENDIF
!         z-derivative
          DCpq(LM,3) = FrZ*Cpq(LM)
          DSpq(LM,3) = Zero
       ENDDO
       RETURN
    ENDIF
!
    OvRad    = One/SQRT(Px*Px+Py*Py+Pz*Pz)
    OvRad2   = OvRad*OvRad
    OvRad3   = OvRad*OvRad2
!
    XY       = SQRT(Px*Px+Py*Py)
    OvXY     = One/XY
    OvXY2    = OvXY*OvXY
!
    DrDx     = Px*OvRad2
    DrDy     = Py*OvRad2
    DrDz     = Pz*OvRad2
!
    DphiDx   =-Py*OvXY2
    DphiDy   = Px*OvXY2
!
    DthetaDx = Px*Pz*OvXY*OvRad2
    DthetaDy = Py*Pz*OvXY*OvRad2
    DthetaDz = -XY*OvRad2
!
    COTT     = Pz*OvXY
!
    DO L=0,MaxL
       FrX      = -DBLE(L+1)*DrDx
       FthetaX1 =  DthetaDx*DBLE(L)*COTT
       FthetaX2 = -DthetaDx*OvXY
       FrY      = -DBLE(L+1)*DrDy
       FthetaY1 =  DthetaDy*DBLE(L)*COTT
       FthetaY2 = -DthetaDy*OvXY
       FrZ      = -DBLE(L+1)*DrDz
       FthetaZ1 =  DthetaDz*DBLE(L)*COTT
       FthetaZ2 = -DthetaDz*OvXY
       DO M = 0,L
          LM  = LTD(L)+M
          LMP = LTD(L-1)+M
!         x-derivative
          DCpq(LM,1) = FrX*Cpq(LM)
          DSpq(LM,1) = FrX*Spq(LM)
          DCpq(LM,1) = DCpq(LM,1) + DBLE(M)*DphiDx*Spq(LM)
          DSpq(LM,1) = DSpq(LM,1) - DBLE(M)*DphiDx*Cpq(LM)
          DCpq(LM,1) = DCpq(LM,1) + FthetaX1*Cpq(LM) + FthetaX2*DBLE(L*L-M*M)*Cpq(LMP)
          DSpq(LM,1) = DSpq(LM,1) + FthetaX1*Spq(LM) + FthetaX2*DBLE(L*L-M*M)*Spq(LMP)
!         y-derivative
          DCpq(LM,2) = FrY*Cpq(LM)
          DSpq(LM,2) = FrY*Spq(LM)
          DCpq(LM,2) = DCpq(LM,2) + DBLE(M)*DphiDy*Spq(LM)
          DSpq(LM,2) = DSpq(LM,2) - DBLE(M)*DphiDy*Cpq(LM)
          DCpq(LM,2) = DCpq(LM,2) + FthetaY1*Cpq(LM) + FthetaY2*DBLE(L*L-M*M)*Cpq(LMP)
          DSpq(LM,2) = DSpq(LM,2) + FthetaY1*Spq(LM) + FthetaY2*DBLE(L*L-M*M)*Spq(LMP)
!         z-derivative
          DCpq(LM,3) = FrZ*Cpq(LM)
          DSpq(LM,3) = FrZ*Spq(LM)
          DCpq(LM,3) = DCpq(LM,3) + FthetaZ1*Cpq(LM) + FthetaZ2*DBLE(L*L-M*M)*Cpq(LMP)
          DSpq(LM,3) = DSpq(LM,3) + FthetaZ1*Spq(LM) + FthetaZ2*DBLE(L*L-M*M)*Spq(LMP)
       ENDDO
    ENDDO
!
  END SUBROUTINE DIrRegular
!========================================================================================
! Dirivative of the GScript Convegence Function
!========================================================================================
  FUNCTION DGScript(MaxL,BetaSq,Px,Py,Pz) RESULT(GDiv)
    INTEGER                    :: MaxL
    REAL(DOUBLE)               :: Px,Py,Pz,BetaSq
    REAL(DOUBLE)               :: DDelta,Factor,RadSq,GRDiv
    REAL(DOUBLE),DIMENSION(3)  :: GDiv
!
    RadSq   = BetaSq*(Px*Px+Py*Py+Pz*Pz)
!   R div
!    GRDiv   = DrGScript(MaxL,RadSq)
    GRDiv   = dRegularizedGammaQ(DBLE(MaxL),RadSq)
!!$
!!$       WRITE(*,*)' RadSq = ',RadSq
!!$       WRITE(*,*)' GRDiv = ',GRDiv
!!$       WRITE(*,*)' GRDIF = ',( GScript(MaxL,RadSq+1D-3)-GScript(MaxL,RadSq-1D-3) )/2.D-3
!!$       WRITE(*,*)' REGd  = ',dRegularizedGammaQ(DBLE(MaxL),RadSq)
!!$

!   x div
    GDiv(1) = Two*Px*BetaSq*GRDiv
!   y div
    GDiv(2) = Two*Py*BetaSq*GRDiv
!   z div
    GDiv(3) = Two*Pz*BetaSq*GRDiv
!
  END FUNCTION DGScript
!========================================================================================
! Dirivative of the FScript Convegence Function
!========================================================================================
  FUNCTION DFScript(MaxL,BetaSq,Px,Py,Pz) RESULT(FDiv)
    INTEGER                    :: MaxL
    REAL(DOUBLE)               :: Px,Py,Pz,BetaSq
    REAL(DOUBLE)               :: DDelta,Factor,RadSq,FRDiv
    REAL(DOUBLE),DIMENSION(3)  :: FDiv
!
    RadSq   = BetaSq*(Px*Px+Py*Py+Pz*Pz)
!   R div
!    FRDiv   = DrFScript(MaxL,RadSq)
    FRDiv=-dRegularizedGammaQ(DBLE(MaxL),RadSq)
!   x div
    FDiv(1) = Two*Px*BetaSq*FRDiv
!   y div
    FDiv(2) = Two*Py*BetaSq*FRDiv
!   z div
    FDiv(3) = Two*Pz*BetaSq*FRDiv
!
  END FUNCTION DFScript
!========================================================================================
!   FT_FSCriptC
!========================================================================================
  FUNCTION DrGScript(L,R)
    INTEGER                    :: L,LS
    REAL(DOUBLE)               :: R,SqrtR,DrGScript,XSUM,DXSUM
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
       DrGScript = -EXP(-R)/(SqrtR*SqrtPi)
    ELSE
       XSUM  =  SFAC(1)*EXP(-R)
       DXSUM = -SFAC(1)*EXP(-R)
       DO LS = 2,L
          XSUM  = XSUM  + (Sfac(LS)*R*EXP(-R/DBLE(LS-1)))**(LS-1)
          DXSUM = DXSUM + ((DBLE(LS-1)-R)/R)*(Sfac(LS)*R*EXP(-R/DBLE(LS-1)))**(LS-1)
       ENDDO
       DrGScript = -EXP(-R)/(SqrtR*SqrtPi)+(Half/(SqrtR*SqrtPi))*XSUM+(SqrtR/SqrtPi)*DXSUM
    ENDIF

  END FUNCTION DrGScript
!========================================================================================
!   FT_FSCriptS
!========================================================================================
  FUNCTION DrFScript(L,R)
    INTEGER                    :: L,LS
    REAL(DOUBLE)               :: R,SqrtR,DrFScript,XSUM,DXSUM
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
       DrFScript = EXP(-R)/(SqrtR*SqrtPi)
    ELSE
       XSUM  =  SFAC(1)*EXP(-R)
       DXSUM = -SFAC(1)*EXP(-R)
       DO LS = 2,L
          XSUM  = XSUM  + (Sfac(LS)*R*EXP(-R/DBLE(LS-1)))**(LS-1)
          DXSUM = DXSUM + ((DBLE(LS-1)-R)/R)*(Sfac(LS)*R*EXP(-R/DBLE(LS-1)))**(LS-1)
       ENDDO
       DrFScript = EXP(-R)/(SqrtR*SqrtPi) - (Half/(SqrtR*SqrtPi))*XSUM-(SqrtR/SqrtPi)*DXSUM
    ENDIF
    !
  END FUNCTION DrFScript
!========================================================================================
! Dirivative of the FT_FScript_C Convegence Function
!========================================================================================
  FUNCTION DFT_FScriptC_3D(MaxL,ExpFac,Px,Py,Pz) RESULT(FT_FCDiv)
    INTEGER                    :: MaxL
    REAL(DOUBLE)               :: Px,Py,Pz,ExpFac
    REAL(DOUBLE)               :: Rad,FTRDiv
    REAL(DOUBLE),DIMENSION(3)  :: FT_FCDiv
!
    Rad =SQRT(Px*Px+Py*Py+Pz*Pz)
!   R div
    FTRDiv = DrFT_FScriptC_3D(MaxL,ExpFac,Rad)
!   x div
    FT_FCDiv(1) = Px*FTRDiv/Rad
!   y div
    FT_FCDiv(2) = Py*FTRDiv/Rad
!   z div
    FT_FCDiv(3) = Pz*FTRDiv/Rad
  END FUNCTION DFT_FScriptC_3D
!========================================================================================
! Dirivative of the FT_FScript_S Convegence Function
!========================================================================================
  FUNCTION DFT_FScriptS_3D(MaxL,ExpFac,Px,Py,Pz) RESULT(FT_FSDiv)
    INTEGER                    :: MaxL
    REAL(DOUBLE)               :: Px,Py,Pz,ExpFac
    REAL(DOUBLE)               :: Rad,FTRDiv
    REAL(DOUBLE),DIMENSION(3)  :: FT_FSDiv
!
    Rad =SQRT(Px*Px+Py*Py+Pz*Pz)
!   R div
    FTRDiv = DrFT_FScriptS_3D(MaxL,ExpFac,Rad)
!   x div
    FT_FSDiv(1) = Px*FTRDiv/Rad
!   y div
    FT_FSDiv(2) = Py*FTRDiv/Rad
!   z div
    FT_FSDiv(3) = Pz*FTRDiv/Rad
  END FUNCTION DFT_FScriptS_3D
!========================================================================================
!   DrFT_FSCriptC
!========================================================================================
  FUNCTION DrFT_FScriptC_3D(L,ExpFac,R)
    INTEGER                    :: L,IFac,Isgn
    REAL(DOUBLE)               :: ExpFac,R,DrFT_FScriptC_3D,Fac
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
       DrFT_FScriptC_3D = Zero
    ELSE
       Isgn = (-1)**(L/2)
       Fac  = ExpFac/DBLE(2*L-1)
       DrFT_FScriptC_3D = (Isgn*(Sfac(L)*R*EXP(-Fac*R*R))**(2*L-1))*(DBLE(2*L-1)/R-Two*ExpFac*R)
    ENDIF
    !
  END FUNCTION DrFT_FScriptC_3D
!========================================================================================
!   DrFT_FSCriptC
!========================================================================================
  FUNCTION DrFT_FScriptS_3D(L,ExpFac,R)
    INTEGER                    :: L,IFac,Isgn
    REAL(DOUBLE)               :: ExpFac,R,DrFT_FScriptS_3D,Fac
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
       DrFT_FScriptS_3D = Zero
    ELSE
       Isgn = (-1)**((L-1)/2)
       Fac  = ExpFac/DBLE(2*L-1)
       DrFT_FScriptS_3D = (Isgn*(Sfac(L)*R*EXP(-Fac*R*R))**(2*L-1))*(DBLE(2*L-1)/R-Two*ExpFac*R)
    ENDIF
    !
  END FUNCTION DrFT_FScriptS_3D
!
!
!
END MODULE DivPFFTen
