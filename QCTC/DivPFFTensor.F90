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
  TYPE(DBL_RNK2)                      :: DCpq,DSpq
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
    CALL New(DCpq,(/LSP(2*MaxEll),3/),(/0,1/))    
    CALL New(DSpq,(/LSP(2*MaxEll),3/),(/0,1/))
    IF(GM%PBC%Dimen==1) THEN
       CALL MakeDivTensor1D(2*MaxEll,GM,Args,CS,dTenC,dTenS)
    ELSEIF(GM%PBC%Dimen==2) THEN
!       CALL MakeDivTensor2D(2*MaxEll,GM,Args,CS,dTenC,dTenS)
    ELSEIF(GM%PBC%Dimen==3) THEN
!       CALL MakeDivTensor3D(2*MaxEll,GM,Args,CS,dTenC,dTenS)
    ENDIF
!
  END SUBROUTINE CalculateDivPFFT
!========================================================================================
! Calculate the Derivative PFFTensor 1D
!========================================================================================
  SUBROUTINE MakeDivTensor1D(MaxL,GM,Args,CS,dTenC,dTenS)
    INTEGER                           :: MaxL
    INTEGER                           :: I,J,L,M,LM,NC,IJ
    TYPE(CellSet)                     :: CS, CSMM
    TYPE(DBL_RNK3)                    :: dTenC,dTenS
    TYPE(CRDS)                        :: GM
    TYPE(ARGMT)                       :: Args
!-------------------------------------------------------------------------------------!
!
!   Number of Inner Boxes
!
    dTenC%D=Zero
    dTenS%D=Zero
    NC = (CS%NCells-1)/2
!
!   One Dimension: Right and Left
!
    IF(GM%PBC%AutoW%I(1)==1) THEN
       IJ = 1
       CALL DIrRegular(MaxL, GM%PBC%BoxShape%D(IJ,IJ),Zero,Zero)
       dTenC%D(:,IJ,IJ) = DCpq%D(:,IJ)
       CALL DIrRegular(MaxL,-GM%PBC%BoxShape%D(IJ,IJ),Zero,Zero)
       dTenC%D(:,IJ,IJ) = dTenC%D(:,IJ,IJ)+DCpq%D(:,IJ)
    ELSEIF(GM%PBC%AutoW%I(2)==1) THEN
       IJ = 2
       CALL DIrRegular(MaxL,Zero, GM%PBC%BoxShape%D(IJ,IJ),Zero)
       dTenC%D(:,IJ,IJ) = DCpq%D(:,IJ)
       CALL DIrRegular(MaxL,Zero,-GM%PBC%BoxShape%D(IJ,IJ),Zero)
       dTenC%D(:,IJ,IJ) = dTenC%D(:,IJ,IJ)+DCpq%D(:,IJ)
    ELSEIF(GM%PBC%AutoW%I(3)==1) THEN 
       IJ = 3
       CALL DIrRegular(MaxL,Zero,Zero, GM%PBC%BoxShape%D(IJ,IJ))
       dTenC%D(:,IJ,IJ) = DCpq%D(:,IJ)
       CALL DIrRegular(MaxL,Zero,Zero,-GM%PBC%BoxShape%D(IJ,IJ))
       dTenC%D(:,IJ,IJ) = dTenC%D(:,IJ,IJ)+DCpq%D(:,IJ)
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
! Calculate the Derivative PFFTensor 3D
!========================================================================================
  SUBROUTINE MakeDivTensor3D(MaxL,GM,Args,CS,dTenC,dTenS)
    INTEGER                           :: MaxL
    INTEGER                           :: I,J,K,L,M,LM,NC
    INTEGER                           :: LSwitch
    TYPE(DBL_RNK3)                    :: dTenC,dTenS
    TYPE(CellSet)                     :: CS, CSMM
    REAL(DOUBLE),DIMENSION(3)         :: PQ,FPQ
    REAL(DOUBLE)                      :: CFac,SFac,BetaSq,Rad,RadSq,ExpFac
    REAL(DOUBLE)                      :: Rmin,RMAX,Accuracy,LenScale,SUM
!
    INTEGER                           :: NDiv
    REAL(DOUBLE)                      :: LenMax,Delt,IntFac
    REAL(DOUBLE),DIMENSION(3)         :: Vec,DCFac,DSFac
    REAL(DOUBLE),DIMENSION(3,3)       :: RecpLatVec,LatVec
    TYPE(CRDS)                        :: GM
    TYPE(ARGMT)                       :: Args
!-------------------------------------------------------------------------------------!
!
!   Initialize 
!
    dTenC%D=Zero
    dTenS%D=Zero
    Accuracy = 1.D-14
    Rmin     = SQRT(CS%CellCarts%D(1,1)**2+CS%CellCarts%D(2,1)**2+CS%CellCarts%D(3,1)**2)
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
!   Three Dimension
!
    LSwitch  = 10
    LenScale = GM%PBC%CellVolume**(One/Three)
    Rmax     = Rmin+LenScale*(One/Accuracy)**(One/DBLE(LSwitch))
    BetaSq   = 0.25D0/(LenScale)**2
!           
!   Sum the Real Space
!
    DO
       CALL New_CellSet_Sphere(CSMM,GM%PBC%AutoW%I,LatVec,Rmax)
       IF(CSMM%NCells .GT. 600000) THEN
          Rmax = 0.99*Rmax
          CALL Delete_CellSet(CSMM)
       ELSE
          EXIT
       ENDIF
    ENDDO
    CALL Sort_CellSet(CSMM)
!
    DO NC = 1,CSMM%NCells
       PQ(:)  = CSMM%CellCarts%D(:,NC)
       FPQ(:) = AtomToFrac(GM,CSMM%CellCarts%D(:,NC))
       IF(.NOT. InCell_CellSet(CS,PQ(1),PQ(2),PQ(3))) THEN
          RadSq = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
          CALL  IrRegular(MaxL,PQ(1),PQ(2),PQ(3))
          CALL DIrRegular(MaxL,PQ(1),PQ(2),PQ(3))
          DO L = 1,MaxL
             IF(L .LE. LSwitch) THEN
                CFac  =  GScript(L,RadSq)
                DCFac = DGScript(L,BetSq,PQ(1),PQ(2),PQ(3))
             ELSE
                GFac  = One
                DGFac = Zero
             ENDIF
             DO M = 0,L
                LM = LTD(L)+M
                DO I=1,3
                   DO J=1,3
                      dTenC%D(LM,I,J)=dTenC%D(LM,I,J)+FPQ(I)*(DCpq(LM,J)*CFac+DCFac(J)*Cpq(LM))
                      dTenS%D(LM,I,J)=dTenS%D(LM,I,J)+FPQ(I)*(DSpq(LM,J)*CFac+DCFac(J)*Spq(LM))
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
    Rmax   = SQRT(ABS(LOG(Accuracy/(10.D0**(LSwitch)))/ExpFac))
    DO
       CALL New_CellSet_Sphere(CSMM,GM%PBC%AutoW%I,RecpLatVec,Rmax)
       IF(CSMM%NCells .LT. 27) THEN
          Rmax = 1.01D0*Rmax
          CALL Delete_CellSet(CSMM)
       ELSE
          EXIT
       ENDIF
    ENDDO
    CALL Sort_CellSet(CSMM)
!
    DO NC = 1,CSMM%NCells
       PQ(:) = CSMM%CellCarts%D(:,NC)
       Rad   = SQRT(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
       IF(Rad .GT. 1.D-14) THEN 
          CALL  IrRegular(MaxL,PQ(1),PQ(2),PQ(3))
          CALL DIrRegular(MaxL,PQ(1),PQ(2),PQ(3))
          DO L = 1,LSwitch
             CFac  = FT_FScriptC_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
             SFac  = FT_FScriptS_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
             DCFac = DFT_FScriptC_3D(L,ExpFac,PQ(1),PQ(2),PQ(3))/GM%PBC%CellVolume
             DSFac = DFT_FScriptS_3D(L,ExpFac,PQ(1),PQ(2),PQ(3))/GM%PBC%CellVolume
             DO M = 0,L
                LM = LTD(L)+M
                DO I=1,3
                   DO J=1,3
                      dTenC%D(LM,I,J)=dTenC%D(LM,I,J)

                      +Cpq(LM)*CFac-Spq(LM)*SFac
                      dTenS%D(LM,I,J)=dTenS%D(LM,I,J)

                      +Spq(LM)*CFac+Cpq(LM)*SFac
                   ENDDO
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
       PQ(:)  = CS%CellCarts%D(:,NC)
       FPQ(:) = AtomToFrac(GM,CSMM%CellCarts%D(:,NC))
       RadSq  = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
       IF(RadSq .GT. 1.D-14) THEN
          CALL IrRegular(MaxL,PQ(1),PQ(2),PQ(3))
          CALL IrRegular(MaxL,PQ(1),PQ(2),PQ(3))
          DO L = 1,LSwitch
             CFac  =  FScript(L,RadSq)
             DCFac = DFScript(L,BetSq,PQ(1),PQ(2),PQ(3))
             DO M = 0,L
                LM = LTD(L)+M
                DO I=1,3
                   DO J=1,3
                      dTenC%D(LM,I,J)=dTenC%D(LM,I,J)-FPQ(I)*(DCpq(LM,J)*CFac+DCFac(J)*Cpq(LM))
                      dTenS%D(LM,I,J)=dTenS%D(LM,I,J)-FPQ(I)*(DSpq(LM,J)*CFac+DCFac(J)*Spq(LM))
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDDO
!
  END SUBROUTINE MakeDivTensor3D

!========================================================================================
! Dirivative of the Irregular Fuction
!========================================================================================
  SUBROUTINE DIrRegular(MaxL,Px,Py,Pz)
    INTEGER                    :: MaxL
    REAL(DOUBLE)               :: Px,Py,Pz
    REAL(DOUBLE)               :: DDelta,Factor
!
    DDelta = 1.d-8
    Factor = Half/DDelta
!   x div
    CALL IrRegular(MaxL,Px+DDelta,Py,Pz)
    DCpq%D(:,1)=Factor*Cpq(:)
    CALL IrRegular(MaxL,Px-DDelta,Py,Pz)
    DCpq%D(:,1)=DCpq%D(:,1)-Factor*Cpq(:)
!   y div
    CALL IrRegular(MaxL,Px,Py+DDelta,Pz)
    DCpq%D(:,2)=Factor*Cpq(:)
    CALL IrRegular(MaxL,Px,Py-DDelta,Pz)
    DCpq%D(:,2)=DCpq%D(:,2)-Factor*Cpq(:)
!   z div
    CALL IrRegular(MaxL,Px,Py,Pz+DDelta)
    DCpq%D(:,3)=Factor*Cpq(:)
    CALL IrRegular(MaxL,Px,Py,Pz-DDelta)
    DCpq%D(:,3)=DCpq%D(:,2)-Factor*Cpq(:)
!
  END SUBROUTINE DIrRegular
!========================================================================================
! Dirivative of the GScript Convegence Function
!========================================================================================
  FUNCTION DGScript(MaxL,BetSq,Px,Py,Pz) RESULT(GDiv)
    INTEGER                    :: MaxL
    REAL(DOUBLE)               :: Px,Py,Pz,BetSq
    REAL(DOUBLE)               :: DDelta,Factor,RadSq,GRDiv
    REAL(DOUBLE),DIMENSION(3)  :: GDiv
!
    DDelta  = 1.D-8
    Factor = Half/DDelta
    RadSq  = BetaSq*(Px*Px+Py*Py+Pz*Pz)
!   R div
    GRDiv = Factor*(GScript(MaxL,RadSq+DDelta)-GScript(MaxL,RadSq-DDelta))
!   x div
    GDiv(1) = Two*Px*BetSq*GRDiv
!   y div
    GDiv(1) = Two*Py*BetSq*GRDiv
!   z div
    GDiv(1) = Two*Pz*BetSq*GRDiv
  END FUNCTION DGScript
!========================================================================================
! Dirivative of the FScript Convegence Function
!========================================================================================
  FUNCTION DFScript(MaxL,BetSq,Px,Py,Pz) RESULT(FDiv)
    INTEGER                    :: MaxL
    REAL(DOUBLE)               :: Px,Py,Pz,BetSq
    REAL(DOUBLE)               :: DDelta,Factor,RadSq,FRDiv
    REAL(DOUBLE),DIMENSION(3)  :: FDiv
!
    DDelta  = 1.D-8
    Factor = Half/DDelta
    RadSq  = BetaSq*(Px*Px+Py*Py+Pz*Pz)
!   R div
    FRDiv = Factor*(FScript(MaxL,RadSq+DDelta)-FScript(MaxL,RadSq-DDelta))
!   x div
    FDiv(1) = Two*Px*BetSq*FRDiv
!   y div
    FDiv(1) = Two*Py*BetSq*FRDiv
!   z div
    FDiv(1) = Two*Pz*BetSq*FRDiv
  END FUNCTION DFScript
!
!
!
END MODULE PFFTen

