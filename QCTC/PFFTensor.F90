MODULE PFFTen
#ifdef PERIODIC  
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
  IMPLICIT NONE
!====================================================================================
! Globals
!====================================================================================
  TYPE(CellSet)                       :: CSMM
  TYPE(DBL_VECT)                      :: TensorC,TensorS
  REAL(DOUBLE), DIMENSION(0:FFLen)    :: RhoC,RhoS
  REAL(DOUBLE), DIMENSION(0:FFLen)    :: PFFBraC,PFFBraS,PFFKetC,PFFKetS,FarFS,FarFC
!
  REAL(DOUBLE),PARAMETER              :: Small=1.0D-12
  CHARACTER(LEN=7),PARAMETER          :: Cube   ='Cube   '
  CHARACTER(LEN=7),PARAMETER          :: Rect   ='Rect   '
  CHARACTER(LEN=7),PARAMETER          :: NonCube='NonCube'
  CHARACTER(LEN=7),PARAMETER          :: Cube111='Cube111'
  CHARACTER(LEN=7),PARAMETER          :: Cube112='Cube112'
  CHARACTER(LEN=7),PARAMETER          :: Cube122='Cube122'
  CHARACTER(LEN=7),PARAMETER          :: Cube222='Cube112'

  CONTAINS
!========================================================================================
! 
!========================================================================================
    SUBROUTINE CalculatePFFT(MaxEll,GM,Args)
      TYPE(CRDS)             :: GM
      TYPE(ARGMT)            :: Args
!-------------------------------------------------------------------------------- 
      INTEGER                :: MaxEll
      INTEGER                :: Layers
      REAL(DOUBLE)           :: Scale,R1,R2,Radius
      CHARACTER(LEN=7)       :: CellType
!
      IF(GM%PBC%Dimen==0) THEN 
!
         CALL New_CellSet_Cube(CS_IN,GM%PBC%AutoW,GM%PBC%BoxShape,(/0,0,0/))    
         CALL Put(1.D-12,'CS_IN%Radius'//CurBase//CurGeom)
         CALL Put_CellSet(CS_IN,'CS_IN'//CurBase//CurGeom)
!
         MaxEll = 0
         CALL New(TensorC,LSP(2*MaxEll),0)
         CALL New(TensorS,LSP(2*MaxEll),0)
         TensorC%D = Zero
         TensorS%D = Zero 
         CALL Put(MaxEll ,'MaxEll')
         CALL Put(TensorC,'PFFTensorC'//CurGeom)
         CALL Put(TensorS,'PFFTensorS'//CurGeom)
!
         RETURN
      ENDIF
! 
!     Make The Tensor
!
      IF(GM%PBC%PFFOvRide) THEN
!
!        Initialize
!
         CALL New(TensorC,LSP(2*MaxEll),0)
         CALL New(TensorS,LSP(2*MaxEll),0)
         TensorC%D = Zero
         TensorS%D = Zero 
!             
!        Make and Store the Inner Box
!
         Layers = GM%PBC%PFFMaxLay
         Radius = DBLE(Layers)*MaxBoxDim(GM)
         CALL New_CellSet_Cube(CS_IN,GM%PBC%AutoW,GM%PBC%BoxShape,(/Layers,Layers,Layers/))
         CALL Sort_CellSet(CS_IN)
!
         CALL Put(Radius,'CS_IN%Radius'//CurBase//CurGeom)
         CALL Put_CellSet(CS_IN,'CS_IN'//CurBase//CurGeom)
!
!        Determine Cell Type
!
         CellType=DetCellType(GM,CS_IN%NCells,Scale)
!     
!        Read from Archive or Make the Tensor
!
         IF(    CellType==Cube111) THEN
            CALL GetExactTensor1L3D(2*MaxEll,GM,Scale)
         ELSEIF(CellType==Cube222) THEN
            CALL GetExactTensor2L3D(2*MaxEll,GM,Scale)
         ELSE
            IF(GM%PBC%Dimen==1) THEN
               CALL MakeTensor1D(2*MaxEll,GM,Args)
            ELSEIF(GM%PBC%Dimen==2) THEN
               CALL MakeTensor2D(2*MaxEll,GM,Args)
            ELSEIF(GM%PBC%Dimen==3) THEN
               CALL MakeTensor3D(2*MaxEll,GM,Args)
            ENDIF
         ENDIF
!
!        Store the Tensor
!
         CALL Put(MaxEll ,'MaxEll'//CurGeom)
         CALL Put(TensorC,'PFFTensorC'//CurGeom)
         CALL Put(TensorS,'PFFTensorS'//CurGeom)
         RETURN
      ELSE
!
!        Initialize
!
         CALL New(TensorC,LSP(2*FFELL),0)
         CALL New(TensorS,LSP(2*FFELL),0)
         TensorC%D = Zero
         TensorS%D = Zero 
!
!        Determin The Minuim Radius
!
         Layers = GM%PBC%PFFMaxLay
         R1     = Two*MaxAtomDist(GM) + SQRT(AtomPairDistanceThreshold)
         R2     = DBLE(Layers)*MaxBoxDim(GM)*(One+1.D-12)
         Radius = MAX(R1,R2)
         GM%PBC%PFFMaxLay = INT(Radius/MaxBoxDim(GM)+Half)
!             
!        Make and Store the Inner Box
!
         CALL New_CellSet_Sphere(CS_IN,GM%PBC%AutoW,GM%PBC%BoxShape,Radius)
         CALL Sort_CellSet(CS_IN)
!
         CALL Put(Radius,'CS_IN%Radius'//CurBase//CurGeom)
         CALL Put_CellSet(CS_IN,'CS_IN'//CurBase//CurGeom)
!
!        Determine Cell Type
!
         CellType=DetCellType(GM,CS_IN%NCells,Scale)
!     
!        Read from Archive or Make the Tensor
!
         IF(    CellType==Cube111) THEN
            CALL GetExactTensor1L3D(2*FFELL,GM,Scale)
         ELSEIF(CellType==Cube222) THEN
            CALL GetExactTensor2L3D(2*FFELL,GM,Scale)
         ELSE
            IF(GM%PBC%Dimen==1) THEN
               CALL MakeTensor1D(2*FFELL,GM,Args)
            ELSEIF(GM%PBC%Dimen==2) THEN
               CALL MakeTensor2D(2*FFELL,GM,Args)
            ELSEIF(GM%PBC%Dimen==3) THEN
               CALL MakeTensor3D(2*FFELL,GM,Args)
            ENDIF
         ENDIF
!
!        Store the Tensor
!
         CALL Put(FFELL  ,'MaxEll'//CurGeom)
         CALL Put(TensorC,'PFFTensorC'//CurGeom)
         CALL Put(TensorS,'PFFTensorS'//CurGeom)
         RETURN
      ENDIF
!
    END SUBROUTINE CalculatePFFT
!========================================================================================
! 
!========================================================================================
    SUBROUTINE GetExactTensor1L3D(MaxL,GM,Scale)
      INTEGER              :: MaxL
      TYPE(CRDS)           :: GM
      INTEGER              :: I,L,M,LM,K,NumberOfRecords
      REAL(DOUBLE)         :: Scale
!
      INCLUDE "PBCTensor/Majik_Kubic_WS1.Inc"			
!
      TensorC%D = Zero
      TensorS%D = Zero 
!
!     Load Majik numbers from parameter statements	
!
      K=0
      DO L=4,MIN(128,MaxL),2
         DO M=0,L,4
            K=K+1	
            LM=LTD(L)+M   
            TensorC%D(LM)=Majik1(K)/(Scale**(DBLE(L)+One))
         ENDDO
      ENDDO
!
     END SUBROUTINE GetExactTensor1L3D
!========================================================================================
! 
!========================================================================================
    SUBROUTINE GetExactTensor2L3D(MaxL,GM,Scale)
      INTEGER              :: MaxL
      TYPE(CRDS)           :: GM
      INTEGER              :: I,L,M,LM,K,NumberOfRecords
      REAL(DOUBLE)         :: Scale
!
      INCLUDE "PBCTensor/Majik_Kubic_WS2.Inc"			
!
      TensorC%D = Zero
      TensorS%D = Zero 
!
!     Load Majik numbers from parameter statements	
!
      K=0
      DO L=4,MIN(128,MaxL),2
         DO M=0,L,4
            K=K+1	
            LM=LTD(L)+M   
            TensorC%D(LM)=Majik2(K)/(Scale**(DBLE(L)+One))
         ENDDO
      ENDDO
!
    END SUBROUTINE GetExactTensor2L3D
!========================================================================================
! Calculate the PFFTensor 1D
!========================================================================================
    SUBROUTINE MakeTensor1D(MaxL,GM,Args)
      INTEGER                           :: MaxL
      INTEGER                           :: I,J,L,M,LM,NC
      REAL(DOUBLE),DIMENSION(3,3)       :: RecpLatVec,LatVec
      TYPE(CRDS)                        :: GM
      TYPE(ARGMT)                       :: Args
!
!     Get The Lattice and Reciprocal Lattice Vectors
!
      DO I = 1,3
         DO J = 1,3
            RecpLatVec(I,J) = GM%PBC%InvBoxSh(J,I)
            LatVec(I,J)     = GM%PBC%BoxShape(I,J)
         ENDDO
      ENDDO
!
!     Number of Inner Boxes
!
      NC = (CS_IN%NCells-1)/2
!
!     One Dimension: Right
!
      IF(GM%PBC%AutoW(1)) THEN
         CALL IrRegular(MaxL,LatVec(1,1),Zero,Zero)
      ELSEIF(GM%PBC%AutoW(2)) THEN
         CALL IrRegular(MaxL,Zero,LatVec(2,2),Zero)
      ELSEIF(GM%PBC%AutoW(3)) THEN      
         CALL IrRegular(MaxL,Zero,Zero,LatVec(3,3))
      ENDIF
      DO L=1,MaxL
         DO M = 0,L
            LM = LTD(L)+M
            TensorC%D(LM) = Cpq(LM)*RZeta(L+1,NC)
         ENDDO
      ENDDO
!
!     One Dimension: Left
!
      IF(GM%PBC%AutoW(1)) THEN
         CALL IrRegular(MaxL,-LatVec(1,1),Zero,Zero)
      ELSEIF(GM%PBC%AutoW(2)) THEN
         CALL IrRegular(MaxL,Zero,-LatVec(2,2),Zero)
      ELSEIF(GM%PBC%AutoW(3)) THEN      
         CALL IrRegular(MaxL,Zero,Zero,-LatVec(3,3))
      ENDIF
      NC = (CS_IN%NCells-1)/2
      DO L=1,MaxL
         DO M = 0,L
            LM = LTD(L)+M
            TensorC%D(LM) =  TensorC%D(LM) + Cpq(LM)*RZeta(L+1,NC)
         ENDDO
      ENDDO
!
    END SUBROUTINE MakeTensor1D
!========================================================================================
! Calculate the PFFTensor 2D
!========================================================================================
    SUBROUTINE MakeTensor2D(MaxL,GM,Args)
      INTEGER                           :: MaxL
      INTEGER                           :: I,J,K,L,M,LM,NC
      INTEGER                           :: LSwitch
      REAL(DOUBLE),DIMENSION(3)         :: PQ
      REAL(DOUBLE)                      :: CFac,SFac,BetaSq,Rad,RadSq,ExpFac
      REAL(DOUBLE)                      :: Rmin,RMAX,Accuracy,LenScale,SUM
!
      INTEGER                           :: NDiv
      REAL(DOUBLE)                      :: LenMax,Delt,IntFac
      REAL(DOUBLE),DIMENSION(3)         :: Vec
      REAL(DOUBLE),DIMENSION(3,3)       :: RecpLatVec,LatVec
      TYPE(CRDS)                        :: GM
      TYPE(ARGMT)                       :: Args
!
!     Initialize 
!
      Accuracy = 1.D-14
      Rmin     = SQRT(CS_IN%CellCarts%D(1,1)**2+CS_IN%CellCarts%D(2,1)**2+CS_IN%CellCarts%D(3,1)**2)
!
!     Get The Lattice and Reciprocal Lattice Vectors
!
      DO I = 1,3
         DO J = 1,3
            RecpLatVec(I,J) = GM%PBC%InvBoxSh(J,I)
            LatVec(I,J)     = GM%PBC%BoxShape(I,J)
         ENDDO
      ENDDO

!
!     Two Dimension
!
      LSwitch  = 10
      LenScale = GM%PBC%CellVolume**(Half)
      Rmax     = Rmin+LenScale*(One/Accuracy)**(One/DBLE(LSwitch))
      BetaSq   = 0.25D0/(LenScale)**2
!           
!     Sum the Real Space
!
      DO
         CALL New_CellSet_Sphere(CSMM,GM%PBC%AutoW,LatVec,Rmin)
         IF(CSMM%NCells .GT. 500000) THEN
            Rmax = 0.99*Rmax
            CALL Delete_CellSet(CSMM)
         ELSE
            EXIT
         ENDIF
      ENDDO
!
      DO NC = 1,CSMM%NCells
         PQ(:) = CSMM%CellCarts%D(:,NC)
         IF(.NOT. InCell_CellSet(CS_IN,PQ(1),PQ(2),PQ(3))) THEN
            RadSq = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
            CALL IrRegular(MaxL,PQ(1),PQ(2),PQ(3))
            DO L = 1,MaxL
               IF(L .LE. LSwitch) THEN
                  CFac = GScript(L,RadSq)
               ELSE
                  CFac = One
               ENDIF
               DO M = 0,L
                  LM = LTD(L)+M
                  TensorC%D(LM)=TensorC%D(LM)+Cpq(LM)*CFac
                  TensorS%D(LM)=TensorS%D(LM)+Spq(LM)*CFac
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      CALL Delete_CellSet(CSMM)
!
!     Sum the Reciprical Space 
!
      ExpFac = Pi*Pi/BetaSq
      Rmax = SQRT(ABS(LOG(Accuracy/(10.D0**(LSwitch)))/ExpFac))
      DO
         CALL New_CellSet_Sphere(CSMM,GM%PBC%AutoW,RecpLatVec,Rmax)
         IF(CSMM%NCells .LT. 9) THEN
            Rmax = 1.01D0*Rmax
         ELSE
            EXIT
         ENDIF
      ENDDO
!
      NDiv    = 1000
      LenMax  = SQRT(ABS(LOG(1.D-6))/ExpFac)
      Delt    = LenMax/DBLE(NDiv)
      DO I=-NDiv,NDiv
         Vec(:) = Zero
         IF(.NOT.GM%PBC%AutoW(1)) Vec(1) = I*Delt
         IF(.NOT.GM%PBC%AutoW(2)) Vec(2) = I*Delt
         IF(.NOT.GM%PBC%AutoW(3)) Vec(3) = I*Delt
         DO NC = 1,CSMM%NCells
            PQ(:) = CSMM%CellCarts%D(:,NC)+Vec(:)
            Rad   = SQRT(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
            IF(Rad .GT. 1.D-14) THEN
               CALL IrRegular(MaxL,PQ(1),PQ(2),PQ(3))
               DO L = 1,LSwitch
                  CFac = Delt*FT_FScriptC_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
                  SFac = Delt*FT_FScriptS_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
                  DO M = 0,L
                     LM = LTD(L)+M
                     TensorC%D(LM)=TensorC%D(LM)+Cpq(LM)*CFac-Spq(LM)*SFac
                     TensorS%D(LM)=TensorS%D(LM)+Spq(LM)*CFac+Cpq(LM)*SFac
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO
!
!        Add in the k=0 piece for the L=2 multipoles (Zero in 3D)
!
      DO I=1,3
         IF(.NOT. GM%PBC%AutoW(I)) THEN
            PQ(:) = Zero
            PQ(I) = 1.D-10
            Rad   = SQRT(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
            L     = 2
            CALL IrRegular(L,PQ(1),PQ(2),PQ(3))
            CFac = Half*Delt*FT_FScriptC_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
            SFac = Half*Delt*FT_FScriptS_3D(L,ExpFac,Rad)/GM%PBC%CellVolume           
            DO M = 0,L
               LM = LTD(L)+M
               TensorC%D(LM)=TensorC%D(LM)+Cpq(LM)*CFac-Spq(LM)*SFac
               TensorS%D(LM)=TensorS%D(LM)+Spq(LM)*CFac+Cpq(LM)*SFac
            ENDDO
            PQ(:) = Zero
            PQ(I) = -1.D-10
            Rad   = SQRT(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
            L     = 2
            CALL IrRegular(L,PQ(1),PQ(2),PQ(3))
            CFac = Half*Delt*FT_FScriptC_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
            SFac = Half*Delt*FT_FScriptS_3D(L,ExpFac,Rad)/GM%PBC%CellVolume           
            DO M = 0,L
               LM = LTD(L)+M
               TensorC%D(LM)=TensorC%D(LM)+Cpq(LM)*CFac-Spq(LM)*SFac
               TensorS%D(LM)=TensorS%D(LM)+Spq(LM)*CFac+Cpq(LM)*SFac
            ENDDO
         ENDIF
      ENDDO
      CALL Delete_CellSet(CSMM)
!
!        Substract the inner boxes
!
      DO NC = 1,CS_IN%NCells
         PQ(:) = CS_IN%CellCarts%D(:,NC)
         RadSq = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
         IF(RadSq .GT. 1.D-14) THEN
            CALL IrRegular(MaxL,PQ(1),PQ(2),PQ(3))
            DO L = 1,LSwitch
               CFac = FScript(L,RadSq)
               DO M = 0,L
                  LM = LTD(L)+M
                  TensorC%D(LM)=TensorC%D(LM)-Cpq(LM)*CFac
                  TensorS%D(LM)=TensorS%D(LM)-Spq(LM)*CFac
               ENDDO
            ENDDO
         ENDIF
      ENDDO
!
    END SUBROUTINE MakeTensor2D
!========================================================================================
! Calculate the PFFTensor 3D
!========================================================================================
    SUBROUTINE MakeTensor3D(MaxL,GM,Args)
      INTEGER                           :: MaxL
      INTEGER                           :: I,J,K,L,M,LM,NC
      INTEGER                           :: LSwitch
      REAL(DOUBLE),DIMENSION(3)         :: PQ
      REAL(DOUBLE)                      :: CFac,SFac,BetaSq,Rad,RadSq,ExpFac
      REAL(DOUBLE)                      :: Rmin,RMAX,Accuracy,LenScale,SUM
!
      INTEGER                           :: NDiv
      REAL(DOUBLE)                      :: LenMax,Delt,IntFac
      REAL(DOUBLE),DIMENSION(3)         :: Vec
      REAL(DOUBLE),DIMENSION(3,3)       :: RecpLatVec,LatVec
      TYPE(CRDS)                        :: GM
      TYPE(ARGMT)                       :: Args
!
!     Initialize 
!
      Accuracy = 1.D-14
      Rmin     = SQRT(CS_IN%CellCarts%D(1,1)**2+CS_IN%CellCarts%D(2,1)**2+CS_IN%CellCarts%D(3,1)**2)
!
!     Get The Lattice and Reciprocal Lattice Vectors
!
      DO I = 1,3
         DO J = 1,3
            RecpLatVec(I,J) = GM%PBC%InvBoxSh(J,I)
            LatVec(I,J)     = GM%PBC%BoxShape(I,J)
         ENDDO
      ENDDO
!
!     Three Dimension
!
      LSwitch  = 10
      LenScale = GM%PBC%CellVolume**(One/Three)
      Rmax     = Rmin+LenScale*(One/Accuracy)**(One/DBLE(LSwitch))
      BetaSq   = 0.25D0/(LenScale)**2
!           
!     Sum the Real Space
!
      DO
         CALL New_CellSet_Sphere(CSMM,GM%PBC%AutoW,LatVec,Rmax)
         IF(CSMM%NCells .GT. 500000) THEN
            Rmax = 0.99*Rmax
            CALL Delete_CellSet(CSMM)
         ELSE
            EXIT
         ENDIF
      ENDDO
      CALL Sort_CellSet(CSMM)
!
      DO NC = 1,CSMM%NCells
         PQ(:) = CSMM%CellCarts%D(:,NC)
         IF(.NOT. InCell_CellSet(CS_IN,PQ(1),PQ(2),PQ(3))) THEN
            RadSq = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
            CALL IrRegular(MaxL,PQ(1),PQ(2),PQ(3))
            DO L = 1,MaxL
               IF(L .LE. LSwitch) THEN
                  CFac = GScript(L,RadSq)
               ELSE
                  CFac = One
               ENDIF
               DO M = 0,L
                  LM = LTD(L)+M
                  TensorC%D(LM)=TensorC%D(LM)+Cpq(LM)*CFac
                  TensorS%D(LM)=TensorS%D(LM)+Spq(LM)*CFac
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      CALL Delete_CellSet(CSMM)
!
!     Sum the Reciprical Space 
!
      ExpFac = Pi*Pi/BetaSq
      Rmax   = SQRT(ABS(LOG(Accuracy/(10.D0**(LSwitch)))/ExpFac))
      DO
         CALL New_CellSet_Sphere(CSMM,GM%PBC%AutoW,RecpLatVec,Rmax)
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
            CALL IrRegular(MaxL,PQ(1),PQ(2),PQ(3))
            DO L = 1,LSwitch
               CFac = FT_FScriptC_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
               SFac = FT_FScriptS_3D(L,ExpFac,Rad)/GM%PBC%CellVolume
               DO M = 0,L
                  LM = LTD(L)+M
                  TensorC%D(LM)=TensorC%D(LM)+Cpq(LM)*CFac-Spq(LM)*SFac
                  TensorS%D(LM)=TensorS%D(LM)+Spq(LM)*CFac+Cpq(LM)*SFac
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      CALL Delete_CellSet(CSMM)
!
!     Substract the inner boxes
!
      DO NC = 1,CS_IN%NCells
         PQ(:) = CS_IN%CellCarts%D(:,NC)
         RadSq = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
         IF(RadSq .GT. 1.D-14) THEN
            CALL IrRegular(MaxL,PQ(1),PQ(2),PQ(3))
            DO L = 1,LSwitch
               CFac = FScript(L,RadSq)
               DO M = 0,L
                  LM = LTD(L)+M
                  TensorC%D(LM)=TensorC%D(LM)-Cpq(LM)*CFac
                  TensorS%D(LM)=TensorS%D(LM)-Spq(LM)*CFac
               ENDDO
            ENDDO
         ENDIF
      ENDDO
!
    END SUBROUTINE MakeTensor3D
!========================================================================================
!   Determine Cell Type
!========================================================================================
    FUNCTION DetCellType(GM,NC,Scale)
      TYPE(CRDS)                   :: GM
      INTEGER                      :: I,NC
      REAL(DOUBLE)                 :: MagA,MagB,MagC,AdotB,AdotC,BdotC,Scale
      CHARACTER(LEN=7)             :: DetCellType
!
!     Determine the cell type
!
      MagA  = Zero
      MagB  = Zero
      MagC  = Zero
      AdotB = Zero
      AdotC = Zero
      BdotC = Zero
      DO I=1,3
         MagA  = MagA+GM%PBC%BoxShape(I,1)*GM%PBC%BoxShape(I,1)
         MagB  = MagB+GM%PBC%BoxShape(I,2)*GM%PBC%BoxShape(I,2)
         MagC  = MagC+GM%PBC%BoxShape(I,3)*GM%PBC%BoxShape(I,3)
         AdotB = AdotB+GM%PBC%BoxShape(I,1)*GM%PBC%BoxShape(I,2)
         AdotC = AdotC+GM%PBC%BoxShape(I,1)*GM%PBC%BoxShape(I,3)
         BdotC = BdotC+GM%PBC%BoxShape(I,2)*GM%PBC%BoxShape(I,3)
      ENDDO      
      MagA  = SQRT(MagA)
      MagB  = SQRT(MagB)      
      MagC  = SQRT(MagC)
      AdotB = SQRT(AdotB/(MagA*MagB))
      AdotC = SQRT(AdotC/(MagA*MagC))
      BdotC = SQRT(BdotC/(MagB*MagC))
!
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
#endif
  END MODULE PFFTen
!!$!========================================================================================
!!$! Calculate the PFF
!!$!========================================================================================
!!$    SUBROUTINE PFFTensor(MaxEll,GM,Args)
!!$      INTEGER              :: MaxEll
!!$      LOGICAL              :: HaveTensor
!!$      TYPE(CRDS)           :: GM
!!$      TYPE(ARGMT)          :: Args
!!$!
!!$      TensorC%D = Zero
!!$      TensorS%D = Zero
!!$      IF(GM%PBC%Dimen==0) RETURN
!!$!
!!$      IF(.NOT. GetTensor(MaxEll,GM,Args)) THEN
!!$         IF(GM%PBC%Dimen==1) THEN
!!$            CALL MakeTensor1D(2*MaxEll,GM,Args)
!!$         ELSEIF(GM%PBC%Dimen==2) THEN
!!$            CALL MakeTensor2D(2*MaxEll,GM,Args)
!!$         ELSEIF(GM%PBC%Dimen==3) THEN
!!$            CALL MakeTensor3D(2*MaxEll,GM,Args)
!!$         ENDIF
!!$         CALL PutTensor(MaxEll,GM,Args)
!!$      ENDIF
!!$!
!!$    END SUBROUTINE PFFTensor
!!$!========================================================================================
!!$! 
!!$!========================================================================================
!!$    FUNCTION GetTensor(MaxL,GM,Args)
!!$      INTEGER              :: MaxL,L,M,LM
!!$      LOGICAL              :: GetTensor
!!$      CHARACTER(LEN=120)   :: FileName
!!$      CHARACTER(LEN=1)     :: AWX,AWY,AWZ
!!$      TYPE(CRDS)           :: GM
!!$      TYPE(ARGMT)          :: Args
!!$!
!!$      AWX = '0'
!!$      AWY = '0'
!!$      AWZ = '0'
!!$      IF(GM%PBC%AutoW(1)) AWX = '1'
!!$      IF(GM%PBC%AutoW(2)) AWY = '1'      
!!$      IF(GM%PBC%AutoW(3)) AWZ = '1'
!!$      FileName= TRIM(MONDO_SCRATCH) // TRIM(Args%C%C(1)) // "_Geom#" // TRIM(CurGeom) //    &
!!$                "_AW"   // TRIM(AWX) // TRIM(AWY) // TRIM(AWZ) //    &
!!$                "_NC"   // TRIM(IntToChar(CS_IN%NCells))   //        &    
!!$                "_LM"   // TRIM(IntToChar(MaxL)) //                  &
!!$                ".PFFT"
!!$!
!!$      INQUIRE(FILE=FileName,EXIST=GetTensor)
!!$      IF(GetTensor) THEN
!!$         OPEN(UNIT=77,FILE=FileName,FORM='UNFORMATTED',STATUS='OLD')
!!$         DO L = 1,MaxL
!!$            DO M = 0,L
!!$               LM = LTD(L)+M
!!$               READ(77) TensorC%D(LM),TensorS%D(LM)
!!$            ENDDO
!!$         ENDDO
!!$         CLOSE(77)
!!$      ENDIF
!!$!
!!$    END FUNCTION GetTensor
!!$!========================================================================================
!!$! 
!!$!========================================================================================
!!$    SUBROUTINE PutTensor(MaxL,GM,Args)
!!$      INTEGER             :: MaxL,IMin,JMin,KMin,L,M,LM
!!$      CHARACTER(LEN=120)  :: FileName
!!$      CHARACTER(LEN=1)    :: AWX,AWY,AWZ
!!$      LOGICAL             :: Exists
!!$      TYPE(CRDS)          :: GM
!!$      TYPE(ARGMT)         :: Args
!!$!
!!$      AWX = '0'
!!$      AWY = '0'
!!$      AWZ = '0'
!!$      IF(GM%PBC%AutoW(1)) AWX = '1'
!!$      IF(GM%PBC%AutoW(2)) AWY = '1'      
!!$      IF(GM%PBC%AutoW(3)) AWZ = '1'
!!$      FileName= TRIM(MONDO_SCRATCH) // TRIM(Args%C%C(1)) // "_Geom#" // TRIM(CurGeom) //    &
!!$                "_AW"   // TRIM(AWX) // TRIM(AWY) // TRIM(AWZ) //    &
!!$                "_NC"   // TRIM(IntToChar(CS_IN%NCells))   //        &    
!!$                "_LM"   // TRIM(IntToChar(MaxL)) //                  &
!!$                ".PFFT"
!!$!
!!$      INQUIRE(FILE=FileName,EXIST=Exists)
!!$      IF(Exists)THEN
!!$         RETURN
!!$      ELSE
!!$         OPEN(FILE=Filename,Unit=77,FORM='UNFORMATTED')
!!$      ENDIF
!!$      DO L = 1,MaxL
!!$         DO M = 0,L
!!$            LM = LTD(L)+M
!!$            WRITE(77) TensorC%D(LM),TensorS%D(LM)
!!$         ENDDO
!!$      ENDDO
!!$      CLOSE(77)
!!$!
!!$    END SUBROUTINE PutTensor
