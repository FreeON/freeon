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
    Pair%AB2 = MinImageDist(GM,GM%Carts%D(1:3,I),GM%Carts%D(1:3,J))
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
!
!   Calculate Radius
!
    IF(PRESENT(Rad_O)) THEN
       Radius = Rad_O
    ELSE
       Radius = Two*MaxAtomDist(GM) + SQRT(AtomPairDistanceThreshold)
    ENDIF
    CALL New_CellSet_Sphere(CS,GM%PBC%AutoW,GM%PBC%BoxShape,Radius)
    CALL Sort_CellSet(CS)
!
  END SUBROUTINE SetCellNumber
!-------------------------------------------------------------------------------
! Convert from Atomic Coordinates  to Fractional Coordinates
!-------------------------------------------------------------------------------
  FUNCTION AtomToFrac(GM,VecA) RESULT(VecF)
    TYPE(CRDS)                 :: GM
    REAL(DOUBLE),DIMENSION(3)  :: VecA,VecF

    VecF(1) = VecA(1)*GM%PBC%InvBoxSh(1,1) + VecA(2)*GM%PBC%InvBoxSh(1,2) + VecA(3)*GM%PBC%InvBoxSh(1,3)
    VecF(2) = VecA(1)*GM%PBC%InvBoxSh(2,1) + VecA(2)*GM%PBC%InvBoxSh(2,2) + VecA(3)*GM%PBC%InvBoxSh(2,3)
    VecF(3) = VecA(1)*GM%PBC%InvBoxSh(3,1) + VecA(2)*GM%PBC%InvBoxSh(3,2) + VecA(3)*GM%PBC%InvBoxSh(3,3)

  END FUNCTION AtomToFrac
!-------------------------------------------------------------------------------
! Convert from Fractional Coordinates to Atomic Coordinates
!-------------------------------------------------------------------------------
  FUNCTION FracToAtom(GM,VecF) RESULT(VecA)
    TYPE(CRDS)                 :: GM
    REAL(DOUBLE),DIMENSION(3)  :: VecA,VecF

    VecA(1) = VecF(1)*GM%PBC%BoxShape(1,1) + VecF(2)*GM%PBC%BoxShape(1,2) + VecF(3)*GM%PBC%BoxShape(1,3)
    VecA(2) = VecF(1)*GM%PBC%BoxShape(2,1) + VecF(2)*GM%PBC%BoxShape(2,2) + VecF(3)*GM%PBC%BoxShape(2,3)
    VecA(3) = VecF(1)*GM%PBC%BoxShape(3,1) + VecF(2)*GM%PBC%BoxShape(3,2) + VecF(3)*GM%PBC%BoxShape(3,3)

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
       IF(GM%PBC%AutoW(I)) THEN
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
    TYPE(CRDS)                         :: G
    REAL(DOUBLE),DIMENSION(3,G%NAtms) :: BoxCarts
    INTEGER        :: I
    !--------------------------------------------------------------------------------!
    ! First account for unwrapped dimensions
    G%Carts%D=G%AbCarts%D
    ! Check for no wrapping at all
    IF(G%PBC%Dimen==0)RETURN
    ! Generate fractional coordinates from unwrapped coordinates
    DO I=1,G%NAtms
       BoxCarts(:,I)=AtomToFrac(G,G%AbCarts%D(:,I))
    ENDDO
    ! Wrap fractional coordinates
    IF(G%PBC%AtomW) THEN
       DO I=1,G%NAtms
          CALL FracCyclic(G,BoxCarts(:,I))
       ENDDO
    ENDIF
    ! Here are the wrapped coordinates
    DO I=1,G%NAtms
       G%Carts%D(:,I)=FracToAtom(G,BoxCarts(:,I))
    ENDDO
  END SUBROUTINE WrapAtoms
!===================================================================================
!
!===================================================================================
  SUBROUTINE  CalFracCarts(GM)
    TYPE(CRDS)                 :: GM
    INTEGER                    :: I
    DO I=1,GM%NAtms
       ! Only compute fracts from unwrapped ab carts:
       GM%BoxCarts%D(:,I) = AtomToFrac(GM,GM%Carts%D(:,I))
       ! This appears to be an incorrect algorithm why would we wrap AbBoxCarts????
       GM%AbBoxCarts%D(:,I) = AtomToFrac(GM,GM%AbCarts%D(:,I))
    ENDDO
  END SUBROUTINE CalFracCarts
!===================================================================================
!
!===================================================================================
  SUBROUTINE  CalAtomCarts(GM)
    TYPE(CRDS)                 :: GM
    INTEGER                    :: I
!   Generate the Atomic Coordinates
    DO I=1,GM%NAtms
       ! Only wrap carts
       GM%Carts%D(:,I)   = FracToAtom(GM,GM%BoxCarts%D(:,I))
       ! This appears to be an incorrect algorithm: Why would we wrap AbCarts????
       GM%AbCarts%D(:,I) = FracToAtom(GM,GM%AbBoxCarts%D(:,I))
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
       GM%Carts%D(:,I)      = GM%Carts%D(:,I)      + ATvec(:)
       GM%AbCarts%D(:,I)    = GM%AbCarts%D(:,I)    + ATvec(:)
       GM%BoxCarts%D(:,I)   = GM%BoxCarts%D(:,I)   + FTvec(:)
       GM%AbBoxCarts%D(:,I) = GM%AbBoxCarts%D(:,I) + FTvec(:)
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
       IF(GM%PBC%AutoW(I)) THEN
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
       IF(GM%PBC%AutoW(K)) NRgn(K) = 1
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
    A(:) = GM%PBC%BoxShape(:,1)+GM%PBC%BoxShape(:,2)+GM%PBC%BoxShape(:,3)
    B(:) = GM%PBC%BoxShape(:,1)+GM%PBC%BoxShape(:,2)-GM%PBC%BoxShape(:,3)
    C(:) = GM%PBC%BoxShape(:,1)-GM%PBC%BoxShape(:,2)-GM%PBC%BoxShape(:,3)

    DO I=1,3
       IF(.NOT. GM%PBC%AutoW(I)) THEN
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
          IF(GM%PBC%AutoW(I)) THEN
             Dist = Dist+(GM%Carts%D(I,At)-GM%PBC%CellCenter(I))**2
          ENDIF
       ENDDO
       MaxAtomDist = MAX(MaxAtomDist,SQRT(Dist))
    ENDDO
    MaxAtomDist = Two*MaxAtomDist
!
  END FUNCTION MaxAtomDist
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
          IF(GM%PBC%AutoW(I).AND.GM%PBC%AutoW(J)) THEN
             LatF(I,J)= nlm(I)*F(J)
          ENDIF
       ENDDO
    ENDDO
!
  END FUNCTION LaticeForce
!-------------------------------------------------------------------------------
! Make GM Periodic
!-------------------------------------------------------------------------------
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
       GM%PBC%CellVolume = One
       DO K=1,3
          IF(GM%PBC%AutoW(K)) THEN
             GM%PBC%CellVolume = GM%PBC%CellVolume*GM%PBC%BoxShape(K,K)
          ENDIF
       ENDDO
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
       GM%PBC%InvBoxSh = InverseMatrix(GM%PBC%BoxShape)
    ENDIF
!   Calculate Atom Positions
    IF(WP(3)) THEN
       CALL CalAtomCarts(GM)
    ENDIF
!
  END SUBROUTINE MakeGMPeriodic
!
END MODULE AtomPairs

