!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!--  Matt Challacombe and  C. J. Tymczak
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!    Compute the Atom Pair and the Cells for PBC
MODULE AtomPairs
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE GlobalObjects
  USE BoundingBox
  USE Thresholding
  USE CellSets
  USE LinAlg
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Global variables
!-------------------------------------------------------------------------------
  TYPE(CellSet)                :: CS_OUT,CS_IN
CONTAINS
!-------------------------------------------------------------------------------
! Test to see if the atomic overlap is above tolerance and create pair
!-------------------------------------------------------------------------------
  FUNCTION SetAtomPair(GM,BS,I,J,Pair,Box_O)
    LOGICAL                   :: SetAtomPair
    INTEGER                   :: I,J,K
    TYPE(AtomPair)            :: Pair
    TYPE(CRDS)                :: GM
    TYPE(BSet)                :: BS
    TYPE(BBox),OPTIONAL       :: Box_O
    Pair%AB2  = MinImageDist(GM,GM%Carts%D(1:3,I),GM%Carts%D(1:3,J))
    Pair%KA   = GM%AtTyp%I(I)
    Pair%KB   = GM%AtTyp%I(J)
    Pair%NA   = BS%BFKnd%I(Pair%KA)
    Pair%NB   = BS%BFKnd%I(Pair%KB)
    Pair%A(:) = GM%Carts%D(:,I) 
    Pair%B(:) = GM%Carts%D(:,J) 
    IF(I==J) THEN 
       Pair%SameAtom = .TRUE.
    ELSE
       Pair%SameAtom = .FALSE.
    ENDIF
!
    SetAtomPair=TestAtomPair(Pair,Box_O)
!
  END FUNCTION SetAtomPair
!-------------------------------------------------------------------------------
! Convert a Block Matrix to a Vector
!-------------------------------------------------------------------------------
  FUNCTION BlockToVect(M,N,A) RESULT (B)
    INTEGER                      :: M,N,I,J,K
    REAL(DOUBLE), DIMENSION(M,N) :: A
    REAL(DOUBLE), DIMENSION(M*N) :: B
    K=0
    DO J=1,N
       DO I=1,M
          K=K+1
          B(K)=A(I,J)
       ENDDO
    ENDDO
  END FUNCTION BlockToVect
!-------------------------------------------------------------------------------
! Convert a Block Matrix to a Vector
!-------------------------------------------------------------------------------
  FUNCTION VectToBlock(M,N,A) RESULT (B)
    INTEGER                      :: M,N,I,J,K
    REAL(DOUBLE), DIMENSION(M*N) :: A
    REAL(DOUBLE), DIMENSION(M,N) :: B
    K=0
    DO J=1,N
       DO I=1,M
          K=K+1
          B(I,J) = A(K)
       ENDDO
    ENDDO
  END FUNCTION VectToBlock
!========================================================================================
! From HG Coefs, calculate the Dipole
!========================================================================================
  FUNCTION CalculateDiPole(LQ,Expt,RX,RY,RZ,Coef)
    INTEGER                   :: LQ
    REAL(DOUBLE)              :: RX,RY,RZ,Expt,PiExpt
    REAL(DOUBLE),DIMENSION(3) :: CalculateDiPole
    REAL(DOUBLE),DIMENSION(:) :: Coef
!
    PiExpt = (Pi/Expt)**1.5D0
    SELECT CASE(LQ)
    CASE (0)
       CalculateDiPole(1) = -PiExpt*Coef(1)*RX
       CalculateDiPole(2) = -PiExpt*Coef(1)*RY
       CalculateDiPole(3) = -PiExpt*Coef(1)*RZ
    CASE(1:)
       CalculateDiPole(1) = -PiExpt*(Coef(2)+Coef(1)*RX)
       CalculateDiPole(2) = -PiExpt*(Coef(3)+Coef(1)*RY)
       CalculateDiPole(3) = -PiExpt*(Coef(4)+Coef(1)*RZ)
    END SELECT
!
  END FUNCTION CalculateDiPole
!========================================================================================
! From HG Coefs, calculate the Quadrupole
!========================================================================================
  FUNCTION CalculateQuadruPole(LQ,Expt,RX,RY,RZ,Coef)
    INTEGER                   :: LQ
    REAL(DOUBLE)              :: RX,RY,RZ,PiExpt,Expt,RX2,RY2,RZ2
    REAL(DOUBLE),DIMENSION(6) :: CalculateQuadruPole
    REAL(DOUBLE),DIMENSION(:) :: Coef
!
    PiExpt = (Pi/Expt)**1.5D0
    RX2 = Half/Expt + RX*RX 
    RY2 = Half/Expt + RY*RY    
    RZ2 = Half/Expt + RZ*RZ 
    SELECT CASE(LQ)
    CASE(0)
       CalculateQuadruPole(1) = PiExpt*Coef(1)*RX2
       CalculateQuadruPole(2) = PiExpt*Coef(1)*RY2
       CalculateQuadruPole(3) = PiExpt*Coef(1)*RZ2
       CalculateQuadruPole(4) = PiExpt*Coef(1)*RX*RY
       CalculateQuadruPole(5) = PiExpt*Coef(1)*RX*RZ
       CalculateQuadruPole(6) = PiExpt*Coef(1)*RY*RZ
    CASE(1)
       CalculateQuadruPole(1) = PiExpt*(Coef(1)*RX2 + Two*Coef(2)*RX)
       CalculateQuadruPole(2) = PiExpt*(Coef(1)*RY2 + Two*Coef(3)*RY)
       CalculateQuadruPole(3) = PiExpt*(Coef(1)*RZ2 + Two*Coef(4)*RZ)
       CalculateQuadruPole(4) = PiExpt*(Coef(1)*RX*RY + Coef(2)*RY + Coef(3)*RX)
       CalculateQuadruPole(5) = PiExpt*(Coef(1)*RX*RZ + Coef(2)*RZ + Coef(4)*RX)
       CalculateQuadruPole(6) = PiExpt*(Coef(1)*RY*RZ + Coef(3)*RZ + Coef(4)*RY)
    CASE(2:)
       CalculateQuadruPole(1) = PiExpt*(Coef(1)*RX2 + Two*Coef(2)*RX + Two*Coef(5))
       CalculateQuadruPole(2) = PiExpt*(Coef(1)*RY2 + Two*Coef(3)*RY + Two*Coef(7))
       CalculateQuadruPole(3) = PiExpt*(Coef(1)*RZ2 + Two*Coef(4)*RZ + Two*Coef(10))
       CalculateQuadruPole(4) = PiExpt*(Coef(1)*RX*RY + Coef(2)*RY + Coef(3)*RX + Coef(6))
       CalculateQuadruPole(5) = PiExpt*(Coef(1)*RX*RZ + Coef(2)*RZ + Coef(4)*RX + Coef(8))
       CalculateQuadruPole(6) = PiExpt*(Coef(1)*RY*RZ + Coef(3)*RZ + Coef(4)*RY + Coef(9))
    END SELECT
!  
  END FUNCTION CalculateQuadruPole
!-------------------------------------------------------------------------------
! Set up the lattice vecters to sum over 
! Also, Added a Factor that Takes into Account the Box Size and Shape
!-------------------------------------------------------------------------------
  SUBROUTINE SetCellNumber(GM,CS,Rad_O)
    TYPE(CRDS)                     :: GM
    TYPE(CellSet)                  :: CS
    REAL(DOUBLE)                   :: Radius
    REAL(DOUBLE), OPTIONAL         :: Rad_O
    INTEGER                        :: IL
!
!   Calculate Radius
!
    IF(PRESENT(Rad_O)) THEN
       Radius = Rad_O
    ELSE
       Radius = (One+1.D-14)*MaxBoxDim(GM) + SQRT(AtomPairDistanceThreshold)
    ENDIF
!   Determine PFFOverde and Make CS
    IF(GM%PBC%PFFOvRide) THEN
       IL = GM%PBC%PFFMaxLay
       CALL New_CellSet_Cube(CS,GM%PBC%AutoW%I,GM%PBC%BoxShape%D,(/IL,IL,IL/))
    ELSE
       CALL New_CellSet_Sphere(CS,GM%PBC%AutoW%I,GM%PBC%BoxShape%D,Radius)
    ENDIF
!
    CALL Sort_CellSet(CS)
    CS%Radius  = SQRT(CS%CellCarts%D(1,1)**2 +CS%CellCarts%D(2,1)**2 +CS%CellCarts%D(3,1)**2)
!
  END SUBROUTINE SetCellNumber
!-------------------------------------------------------------------------------
! Convert from Atomic Coordinates  to Fractional Coordinates
!-------------------------------------------------------------------------------
  FUNCTION AtomToFrac(GM,VecA) RESULT(VecF)
    TYPE(CRDS)                 :: GM
    REAL(DOUBLE),DIMENSION(3)  :: VecA,VecF

    VecF(1) = VecA(1)*GM%PBC%InvBoxSh%D(1,1) + VecA(2)*GM%PBC%InvBoxSh%D(1,2) + VecA(3)*GM%PBC%InvBoxSh%D(1,3)
    VecF(2) = VecA(1)*GM%PBC%InvBoxSh%D(2,1) + VecA(2)*GM%PBC%InvBoxSh%D(2,2) + VecA(3)*GM%PBC%InvBoxSh%D(2,3)
    VecF(3) = VecA(1)*GM%PBC%InvBoxSh%D(3,1) + VecA(2)*GM%PBC%InvBoxSh%D(3,2) + VecA(3)*GM%PBC%InvBoxSh%D(3,3)

  END FUNCTION AtomToFrac
!-------------------------------------------------------------------------------
! Convert from Fractional Coordinates to Atomic Coordinates
!-------------------------------------------------------------------------------
  FUNCTION FracToAtom(GM,VecF) RESULT(VecA)
    TYPE(CRDS)                 :: GM
    REAL(DOUBLE),DIMENSION(3)  :: VecA,VecF

    VecA(1) = VecF(1)*GM%PBC%BoxShape%D(1,1) + VecF(2)*GM%PBC%BoxShape%D(1,2) + VecF(3)*GM%PBC%BoxShape%D(1,3)
    VecA(2) = VecF(1)*GM%PBC%BoxShape%D(2,1) + VecF(2)*GM%PBC%BoxShape%D(2,2) + VecF(3)*GM%PBC%BoxShape%D(2,3)
    VecA(3) = VecF(1)*GM%PBC%BoxShape%D(3,1) + VecF(2)*GM%PBC%BoxShape%D(3,2) + VecF(3)*GM%PBC%BoxShape%D(3,3)

  END FUNCTION FracToAtom
!-------------------------------------------------------------------------------
! In Perodic Systems, Cyclically put positions back in the Fracional Box (in Frac Coord)
!-------------------------------------------------------------------------------
  SUBROUTINE FracCyclic(GM,VecF)
    TYPE(CRDS)                 :: GM        
    REAL(DOUBLE),DIMENSION(3)  :: VecF
    INTEGER                    :: I,N     
    REAL(DOUBLE),PARAMETER     :: MDelt = 1.D-12

    DO I=1,3
       IF(GM%PBC%AutoW%I(I)==1) THEN
          DO 
             IF(VecF(I) > (One+MDelt)) THEN
                VecF(I) = VecF(I) - One
             ELSE
                EXIT
             ENDIF
          ENDDO
          DO 
             IF(VecF(I) < (Zero-MDelt)) THEN
                VecF(I) = VecF(I) + One
             ELSE
                EXIT
             ENDIF
          ENDDO
       ENDIF
    ENDDO
  END SUBROUTINE FracCyclic
!-------------------------------------------------------------------------------
! In Perodic Systems, Cyclically put positions back in the Atomic Box (in Atom Coord)
!-------------------------------------------------------------------------------
  SUBROUTINE AtomCyclic(GM,VecA)
    TYPE(CRDS)                 :: GM        
    REAL(DOUBLE),DIMENSION(3)  :: VecA,VecF    
!
    VecF = AtomToFrac(GM,VecA)
    CALL FracCyclic(GM,VecF)
    VecA = FracToAtom(GM,VecF)
!
  END SUBROUTINE AtomCyclic
!===================================================================================
!
!===================================================================================
  SUBROUTINE  WrapAtoms(G)
    TYPE(CRDS)     :: G
    INTEGER        :: I
!   Check for no wrapping at all
    IF(G%PBC%Dimen==0) RETURN
!   Wrap Fractional coordinates
    DO I=1,G%NAtms
       CALL FracCyclic(G,G%BoxCarts%D(:,I))
    ENDDO
!   Wrap Atomic  coordinates 
    DO I=1,G%NAtms
       CALL AtomCyclic(G,G%Carts%D(:,I))
    ENDDO
!
  END SUBROUTINE WrapAtoms
!===================================================================================
!
!===================================================================================
  SUBROUTINE  CalFracCarts(GM)
    TYPE(CRDS)                 :: GM
    INTEGER                    :: I
    DO I=1,GM%NAtms
!      Compute fracts from Carts:
       GM%BoxCarts%D(:,I) = AtomToFrac(GM,GM%Carts%D(:,I))
    ENDDO
  END SUBROUTINE CalFracCarts
!===================================================================================
!
!===================================================================================
  SUBROUTINE  CalAtomCarts(GM)
    TYPE(CRDS)                 :: GM
    INTEGER                    :: I
    DO I=1,GM%NAtms
!      Only wrap carts
       GM%Carts%D(:,I)   = FracToAtom(GM,GM%BoxCarts%D(:,I))
    ENDDO
  END SUBROUTINE CalAtomCarts
!===================================================================================
!
!===================================================================================
  SUBROUTINE  Translate(GM,ATvec)
    TYPE(CRDS)                 :: GM
    REAL(DOUBLE),DIMENSION(3)  :: ATvec,FTvec
    INTEGER                    :: I
!
    FTvec(:) = AtomToFrac(GM,ATvec(:))
!
!   Tranaslate The Atoms
!
    DO I=1,GM%NAtms
       GM%Carts%D(:,I)      = GM%Carts%D(:,I)    + ATvec(:)
       GM%AbCarts%D(:,I)    = GM%AbCarts%D(:,I)  + ATvec(:)
       GM%BoxCarts%D(:,I)   = GM%BoxCarts%D(:,I) + FTvec(:)
    ENDDO
!
  END SUBROUTINE Translate
!-------------------------------------------------------------------------------
! Test to See if in the Box (Fractional Coordinates)
!-------------------------------------------------------------------------------
  FUNCTION InFracBox(GM,VecF)
    TYPE(CRDS)                 :: GM     
    LOGICAL                    :: InFracBox
    INTEGER                    :: I
    REAL(DOUBLE),DIMENSION(3)  :: VecF

    InFracBox = .TRUE.
    DO I=1,3
       IF(GM%PBC%AutoW%I(I)==1) THEN
          IF(VecF(I) < zero .OR. VecF(I) > one) THEN
             InFracBox = .FALSE.
             RETURN
          ENDIF
       ENDIF
    ENDDO

  END FUNCTION InFracBox
!-------------------------------------------------------------------------------
! Test to See if in the Box (Atomic Coordinates)
!-------------------------------------------------------------------------------
  FUNCTION InAtomBox(GM,VecA)
    TYPE(CRDS)                 :: GM   
    LOGICAL                    :: InAtomBox
    REAL(DOUBLE),DIMENSION(3)  :: VecA,VecF     

    VecF = AtomToFrac(GM,VecA)
    InAtomBox = InFracBox(GM,VecF)

  END FUNCTION InAtomBox
!-------------------------------------------------------------------------------
! Calculate the Minmiun Image Distance Between Atoms I and J 
!-------------------------------------------------------------------------------
  FUNCTION MinImageDist(GM,Carts1,Carts2)
    TYPE(CRDS)                :: GM
    INTEGER                   :: K,N1,N2,N3
    REAL(DOUBLE)              :: MinImageDist
    INTEGER,DIMENSION(3)      :: NRgn
    REAL(DOUBLE),DIMENSION(3) :: Carts1,Carts2
    REAL(DOUBLE),DIMENSION(3) :: VecF,VecA
!  
    NRgn = 0
    DO K=1,3
       IF(GM%PBC%AutoW%I(K)==1) NRgn(K) = 1
    ENDDO
!
    MinImageDist = 1.D16
    DO N1 = -NRgn(1),NRgn(1)
       DO N2 = -NRgn(2),NRgn(2)
          DO N3 = -NRgn(3),NRgn(3)
             VecF         = (/DBLE(N1),DBLE(N2),DBLE(N3)/)
             VecA         = FracToAtom(GM,VecF)+Carts1-Carts2
             MinImageDist = MIN(MinImageDist,VecA(1)**2+VecA(2)**2+VecA(3)**2)
          ENDDO
       ENDDO
    ENDDO

  END FUNCTION MinImageDist
!-------------------------------------------------------------------------------
! Calculate the Minmiun Image Distance Between Atoms I and J  via Box Size
!-------------------------------------------------------------------------------
  FUNCTION MaxBoxDim(GM)
    TYPE(CRDS)                  :: GM
    INTEGER                     :: I
    REAL(DOUBLE)                :: MaxBoxDim
    REAL(DOUBLE),DIMENSION(3)   :: A,B,C
    !
    MaxBoxDim = Zero
    A(:) = GM%PBC%BoxShape%D(:,1)+GM%PBC%BoxShape%D(:,2)+GM%PBC%BoxShape%D(:,3)
    B(:) = GM%PBC%BoxShape%D(:,1)+GM%PBC%BoxShape%D(:,2)-GM%PBC%BoxShape%D(:,3)
    C(:) = GM%PBC%BoxShape%D(:,1)-GM%PBC%BoxShape%D(:,2)-GM%PBC%BoxShape%D(:,3)

    DO I=1,3
       IF(GM%PBC%AutoW%I(I)==0) THEN
         A(I) = Zero
         B(I) = Zero
         C(I) = Zero
       ENDIF
    ENDDO
    MaxBoxDim = SQRT(A(1)**2+A(2)**2+A(3)**2)
    MaxBoxDim = MAX(MaxBoxDim,SQRT(B(1)**2+B(2)**2+B(3)**2))
    MaxBoxDim = MAX(MaxBoxDim,SQRT(C(1)**2+C(2)**2+C(3)**2))
!
  END FUNCTION MaxBoxDim
!-------------------------------------------------------------------------------
! Calculate the Minmiun Image Distance Between Atoms I and J via Atom Positions
!-------------------------------------------------------------------------------
  FUNCTION MaxAtomDist(GM)
    TYPE(CRDS)                  :: GM
    INTEGER                     :: I,At
    REAL(DOUBLE)                :: MaxAtomDist,Dist
!
    MaxAtomDist = Zero
    DO At=1,GM%NAtms
       Dist = Zero
       DO I=1,3
          IF(GM%PBC%AutoW%I(I)==1) THEN
             Dist = Dist+(GM%Carts%D(I,At)-GM%PBC%CellCenter%D(I))**2
          ENDIF
       ENDDO
       MaxAtomDist = MAX(MaxAtomDist,SQRT(Dist))
    ENDDO
    MaxAtomDist = Two*MaxAtomDist
!
  END FUNCTION MaxAtomDist
!-------------------------------------------------------------------------------
!  Do the Latice Force Contraction
!-------------------------------------------------------------------------------
  FUNCTION LaticeForce(GM,nlm,F) RESULT(LatF)
    TYPE(CRDS)                     :: GM
    INTEGER                        :: I,J
    REAL(DOUBLE),DIMENSION(3)      :: nlm,F
    REAL(DOUBLE),DIMENSION(3,3)    :: LatF
!
    LatF=Zero
    DO I=1,3
       DO J=1,3
          IF(GM%PBC%AutoW%I(I)==1 .AND. GM%PBC%AutoW%I(J)==1) THEN
             LatF(I,J)= nlm(J)*F(I)
          ENDIF
       ENDDO
    ENDDO
!
  END FUNCTION LaticeForce
!-------------------------------------------------------------------------------
! Make GM Periodic
!-------------------------------------------------------------------------------
  SUBROUTINE PBCInfoFromNewCarts(PBC)
    TYPE(PBCInfo)             :: PBC
!
!   This routine rebuilds PBC data based on PBC%BoxShape
!
    PBC%CellVolume   = CellVolume(PBC%BoxShape%D,PBC%AutoW%I)
    PBC%CellCenter%D = CellCenter(PBC%BoxShape%D,PBC%AutoW%I)
    PBC%InvBoxSh%D   = InverseMatrix(PBC%BoxShape%D)
!
  END SUBROUTINE PBCInfoFromNewCarts
!-------------------------------------------------------------------------------
  FUNCTION CellCenter(BoxShape,AutoW)
    REAL(DOUBLE),DIMENSION(3)   :: CellCenter
    REAL(DOUBLE),DIMENSION(3,3) :: BoxShape
    INTEGER,DIMENSION(3)        :: AutoW
    INTEGER                     :: I,J
    !
!   Find the center of the cell     
    DO I=1,3
       CellCenter(I)=Zero
       IF(AutoW(I)==1)THEN           
          DO J=1,3
             IF(AutoW(J)==1) THEN
                CellCenter(I)=CellCenter(I)+Half*BoxShape(I,J)
             ENDIF
          ENDDO
       ENDIF
    ENDDO
  END FUNCTION CellCenter
!-------------------------------------------------------------------------------
  FUNCTION CellVolume(BoxShape,AutoW)
    REAL(DOUBLE)                :: CellVolume
    REAL(DOUBLE),DIMENSION(3,3) :: BoxShape
    REAL(DOUBLE),DIMENSION(3)   :: AB
    INTEGER,DIMENSION(3)        :: AutoW
    INTEGER                     :: I,J,D
!
    D=0
    DO I=1,3
       D = D + AutoW(I)
    ENDDO
!
    IF(D==1) THEN
       IF(AutoW(1)==1) CellVolume=BoxShape(1,1)
       IF(AutoW(2)==1) CellVolume=BoxShape(2,2)
       IF(AutoW(3)==1) CellVolume=BoxShape(3,3)
    ELSEIF(D==2) THEN
       IF(AutoW(1)==0) THEN
          CellVolume = BoxShape(2,2)*BoxShape(3,3)-BoxShape(3,2)*BoxShape(2,3)
       ENDIF  
       IF(AutoW(2)==0) THEN
          CellVolume = BoxShape(1,1)*BoxShape(3,3)-BoxShape(1,3)*BoxShape(3,1)
       ENDIF  
       IF(AutoW(3)==0) THEN
          CellVolume = BoxShape(1,1)*BoxShape(2,2)-BoxShape(1,2)*BoxShape(2,1)
       ENDIF  
    ELSEIF(D==3) THEN
       CALL CROSS_PRODUCT(BoxShape(:,1),BoxShape(:,2),AB)
       CellVolume = DOT_PRODUCT(AB,BoxShape(:,3))
    ENDIF
!
  END FUNCTION CellVolume
!-------------------------------------------------------------------------------
  FUNCTION DivCellVolume(BoxShape,AutoW) RESULT(DivCV)
    REAL(DOUBLE),DIMENSION(3,3) :: BoxShape,DivCV,Temp
    INTEGER,DIMENSION(3)        :: AutoW
    REAL(DOUBLE)                :: Phase
    INTEGER                     :: I,J,D
!
    D=0
    DO I=1,3
       D = D + AutoW(I)
    ENDDO
!
    DivCV = Zero
    IF(D==1) THEN
       IF(AutoW(1)==1) DivCV(1,1) = One
       IF(AutoW(2)==1) DivCV(2,2) = One
       IF(AutoW(3)==1) DivCV(3,3) = One
    ELSEIF(D==2) THEN
       IF(AutoW(1)==0) THEN
          DivCV(2,2) =  BoxShape(3,3)
          DivCV(2,3) = -BoxShape(3,2)
          DivCV(3,2) = -BoxShape(2,3)
          DivCV(3,3) =  BoxShape(2,2)
       ENDIF
       IF(AutoW(2)==0) THEN
          DivCV(1,1) =  BoxShape(3,3)
          DivCV(1,3) = -BoxShape(3,1)
          DivCV(3,1) = -BoxShape(1,3)
          DivCV(3,3) =  BoxShape(1,1)
       ENDIF
       IF(AutoW(3)==0) THEN
          DivCV(1,1) =  BoxShape(2,2)
          DivCV(1,2) = -BoxShape(2,1)
          DivCV(2,1) = -BoxShape(1,2)
          DivCV(2,2) =  BoxShape(1,1)
       ENDIF
    ELSEIF(D==3) THEN
       DO I=1,3
          DO J=1,3
             Temp = BoxShape
             Temp(:,J) = Zero
             Temp(I,J) = One
             DivCV(I,J) = CellVolume(Temp,AutoW)             
          ENDDO
       ENDDO
    ENDIF
!
    Phase = CellVolume(BoxShape,AutoW)/ABS(CellVolume(BoxShape,AutoW))
    DivCV = Phase*DivCV
!
  END FUNCTION DivCellVolume

!--------------------------------------------------------------------------------
  SUBROUTINE BoxParsToCart(Vec,BoxShape)
    REAL(DOUBLE),DIMENSION(6)   :: Vec
    REAL(DOUBLE),DIMENSION(3,3) :: BoxShape
    BoxShape=Zero
    BoxShape(1,1)=Vec(1)
    BoxShape(1,2)=Vec(2)*COS(Vec(6))
    BoxShape(2,2)=Vec(2)*SIN(Vec(6))
    BoxShape(1,3)=Vec(3)*COS(Vec(5))
    BoxShape(2,3)=(Vec(2)*Vec(3)*COS(Vec(4)) &
                  -BoxShape(1,2)*BoxShape(1,3))/BoxShape(2,2) 
    BoxShape(3,3)=SQRT(Vec(3)**2-BoxShape(1,3)**2-BoxShape(2,3)**2) 
  END SUBROUTINE BoxParsToCart
!--------------------------------------------------------------------
  SUBROUTINE CalcBoxPars(Vec,BoxShape)
    REAL(DOUBLE),DIMENSION(6)  :: Vec
    REAL(DOUBLE),DIMENSION(3)  :: VecA,VecB,VecC,VecCr
    REAL(DOUBLE),DIMENSION(:,:):: BoxShape
    !
    Vec(1)=SQRT(DOT_PRODUCT(BoxShape(1:3,1),BoxShape(1:3,1)))
    Vec(2)=SQRT(DOT_PRODUCT(BoxShape(1:3,2),BoxShape(1:3,2)))
    Vec(3)=SQRT(DOT_PRODUCT(BoxShape(1:3,3),BoxShape(1:3,3)))
    VecA=BoxShape(1:3,1)/Vec(1)
    VecB=BoxShape(1:3,2)/Vec(2)
    VecC=BoxShape(1:3,3)/Vec(3)
    Vec(4)=DOT_PRODUCT(VecB,VecC)
    Vec(5)=DOT_PRODUCT(VecC,VecA)
    Vec(6)=DOT_PRODUCT(VecA,VecB)
    Vec(4)=ACOS(Vec(4)) 
    Vec(5)=ACOS(Vec(5)) 
    Vec(6)=ACOS(Vec(6)) 
  END SUBROUTINE CalcBoxPars
!--------------------------------------------------------------------
  SUBROUTINE MakeGMPeriodic(GM,WP_O)
    TYPE(CRDS)                     :: GM
    INTEGER                        :: K
    LOGICAL,OPTIONAL,DIMENSION(3)  :: WP_O
    LOGICAL,DIMENSION(3)           :: WP
!
!   Calculate the Cell Volume
!
    IF(PRESENT(WP_O)) THEN
       WP = WP_O
    ELSE
       WP = (/.true.,.true.,.true./)
    ENDIF
!   Calculate the Volume and Dipole term
    IF(WP(1)) THEN
       GM%PBC%CellVolume = CellVolume(GM%PBC%BoxShape%D,GM%PBC%AutoW%I)
!      Calculate the Dipole and Quadripole Factors
       IF(GM%PBC%Dimen < 2) THEN
          GM%PBC%DipoleFAC = Zero
          GM%PBC%QupoleFAC = Zero
       ELSEIF(GM%PBC%Dimen ==2) THEN
          GM%PBC%DipoleFAC = (Four*Pi/GM%PBC%CellVolume)*(One/(GM%PBC%Epsilon+One))
          GM%PBC%QupoleFAC =  Zero
          IF(ABS(GM%PBC%DipoleFAC) .LT. 1.D-14) GM%PBC%DipoleFAC = Zero
       ELSEIF(GM%PBC%Dimen ==3) THEN
          GM%PBC%DipoleFAC = -(Four*Pi/GM%PBC%CellVolume)*(One/Three - One/(Two*GM%PBC%Epsilon+One))
          GM%PBC%QupoleFAC =  (Two*Pi/GM%PBC%CellVolume)*(One/Three - One/(Two*GM%PBC%Epsilon+One))
          IF(ABS(GM%PBC%DipoleFAC) .LT. 1.D-14) GM%PBC%DipoleFAC = Zero
          IF(ABS(GM%PBC%QupoleFAC) .LT. 1.D-14) GM%PBC%QupoleFAC = Zero
       ENDIF
    ENDIF
!   Calculate The Inverse of BoxShape  InvBoxSh = [BoxShape]^(-1)
    IF(WP(2)) THEN
       GM%PBC%InvBoxSh%D = InverseMatrix(GM%PBC%BoxShape%D)
    ENDIF
!   Calculate Atom Positions
    IF(WP(3)) THEN
       CALL CalAtomCarts(GM)
    ENDIF
!
  END SUBROUTINE MakeGMPeriodic
!-------------------------------------------------------------------------------
!   Calculate Inverse Matrix
!-------------------------------------------------------------------------------
  FUNCTION InverseMatrix(Mat) RESULT(InvMat)
    REAL(DOUBLE)                    :: Det,Norm
    REAL(DOUBLE),DIMENSION(3,3)     :: Mat,InvMat
!
    Det  = Mat(1,1)*Mat(2,2)*Mat(3,3) + Mat(1,2)*Mat(2,3)*Mat(3,1) + Mat(1,3)*Mat(2,1)*Mat(3,2) &
         - Mat(1,3)*Mat(2,2)*Mat(3,1) - Mat(1,1)*Mat(2,3)*Mat(3,2) - Mat(1,2)*Mat(2,1)*Mat(3,3)
    Norm = One/Det
    InvMat(1,1) = Norm*(Mat(2,2)*Mat(3,3) - Mat(2,3)*Mat(3,2))
    InvMat(1,2) = Norm*(Mat(1,3)*Mat(3,2) - Mat(1,2)*Mat(3,3))
    InvMat(1,3) = Norm*(Mat(1,2)*Mat(2,3) - Mat(1,3)*Mat(2,2))
    InvMat(2,1) = Norm*(Mat(2,3)*Mat(3,1) - Mat(2,1)*Mat(3,3))
    InvMat(2,2) = Norm*(Mat(1,1)*Mat(3,3) - Mat(1,3)*Mat(3,1))
    InvMat(2,3) = Norm*(Mat(1,3)*Mat(2,1) - Mat(1,1)*Mat(2,3))
    InvMat(3,1) = Norm*(Mat(2,1)*Mat(3,2) - Mat(2,2)*Mat(3,1))
    InvMat(3,2) = Norm*(Mat(1,2)*Mat(3,1) - Mat(1,1)*Mat(3,2))
    InvMat(3,3) = Norm*(Mat(1,1)*Mat(2,2) - Mat(1,2)*Mat(2,1))
!
  END FUNCTION InverseMatrix
!----------------------------------------------------------------
   SUBROUTINE ConvertToXYZRef(XYZ,RefXYZ,PBCDim,BoxShape_O)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ,RefXYZ
     REAL(DOUBLE),DIMENSION(3)   :: VectA,VectB,VectC,VectAux
     REAL(DOUBLE),DIMENSION(3,27):: TrPos
     REAL(DOUBLE),DIMENSION(27)  :: Dist 
     REAL(DOUBLE)                :: Dist0,DistCrit
     INTEGER                     :: PBCDim,NA,NB,NC
     INTEGER                     :: NatmsLoc,I,J,III,IA,IB,IC,ICrit
     REAL(DOUBLE),DIMENSION(3,3),OPTIONAL :: BoxShape_O
     !
     IF(PBCDim==0) RETURN
     IF(PRESENT(BoxShape_O)) THEN
       NatmsLoc=SIZE(XYZ,2)
       DO I=1,3
         VectA(I)=BoxShape_O(I,1)
         VectB(I)=BoxShape_O(I,2)
         VectC(I)=BoxShape_O(I,3)
       ENDDO
     ELSE
       NatmsLoc=SIZE(XYZ,2)-3
       DO I=1,3
         VectA(I)=XYZ(I,NatmsLoc+1)
         VectB(I)=XYZ(I,NatmsLoc+2)
         VectC(I)=XYZ(I,NatmsLoc+3)
       ENDDO
     ENDIF
     !
     NA=0 ; NB=0 ; NC=0
     IF(PBCDim>0) NA=1
     IF(PBCDim>1) NB=1
     IF(PBCDim>2) NC=1
     VectAux=VectA+VectB+VectC
     Dist0=DOT_PRODUCT(VectAux,VectAux)
     DO I=1,NatmsLoc
       III=0
       DistCrit=Dist0
       DO IA=-NA,NA   
         DO IB=-NB,NB   
           DO IC=-NC,NC   
             III=III+1
             TrPos(1:3,III)=XYZ(1:3,I)+IA*VectA+IB*VectB+IC*VectC
             VectAux=TrPos(1:3,III)-RefXYZ(1:3,I)
             Dist(III)=DOT_PRODUCT(VectAux,VectAux)
             IF(Dist(III)<DistCrit) THEN
               ICrit=III
               DistCrit=Dist(III)
             ENDIF
           ENDDO
         ENDDO
       ENDDO
       XYZ(1:3,I)=TrPos(1:3,ICrit) 
     ENDDO
     !
   END SUBROUTINE ConvertToXYZRef
END MODULE AtomPairs

