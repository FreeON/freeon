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
#ifdef PERIODIC
  USE CellSets
#endif
  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Global variables
!-------------------------------------------------------------------------------
#ifdef PERIODIC
  TYPE(CellSet)                :: CS_OUT,CS_IN
  TYPE(CellSet)                :: CS_Kxc,CS_Grid
#endif
CONTAINS
!-------------------------------------------------------------------------------
! Test to see if the atomic overlap is above tolerance and create pair
!-------------------------------------------------------------------------------
  FUNCTION SetAtomPair(GM,BS,I,J,Pair)
    LOGICAL                   :: SetAtomPair
    INTEGER                   :: I,J,K
    TYPE(AtomPair)            :: Pair
    TYPE(CRDS)                :: GM
    TYPE(BSet)                :: BS
#ifdef PERIODIC
    Pair%AB2 = MinImageDist(GM,GM%Carts%D(1:3,I),GM%Carts%D(1:3,J))
#else
    Pair%AB2 = (GM%Carts%D(1,I)-GM%Carts%D(1,J))**2 &
             + (GM%Carts%D(2,I)-GM%Carts%D(2,J))**2 &
             + (GM%Carts%D(3,I)-GM%Carts%D(3,J))**2 
#endif
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
    SetAtomPair=TestAtomPair(Pair)
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
#ifdef PERIODIC
!-------------------------------------------------------------------------------
! Set up the lattice vecters to sum over 
! Also, Added a Factor that Takes into Account the Box Size and Shape
!-------------------------------------------------------------------------------
  SUBROUTINE SetCellNumber(GM)
    TYPE(CRDS)                     :: GM
    INTEGER                        :: I       
    REAL(DOUBLE)                   :: Radius,Radd
    REAL(DOUBLE),DIMENSION(3,3)    :: ABC
!
    ABC(:,1) = GM%PBC%BoxShape(:,1)+GM%PBC%BoxShape(:,2)+GM%PBC%BoxShape(:,3)
    ABC(:,2) = GM%PBC%BoxShape(:,1)+GM%PBC%BoxShape(:,2)-GM%PBC%BoxShape(:,3)
    ABC(:,3) = GM%PBC%BoxShape(:,1)-GM%PBC%BoxShape(:,2)-GM%PBC%BoxShape(:,3)
    Radd = Zero
    DO I=1,3
       IF(GM%PBC%AutoW(I)) THEN
          Radd = MAX(Radd,SQRT(ABC(1,I)**2+ABC(2,I)**2+ABC(3,I)**2))
       ENDIF
    ENDDO
!
    Radius = Radd+SQRT(AtomPairDistanceThreshold)
    CALL New_CellSet_Sphere(CS_OUT,GM%PBC%AutoW,GM%PBC%BoxShape,Radius)
!
  END SUBROUTINE SetCellNumber
!-------------------------------------------------------------------------------
! Set Up the Cell Set for the HiCu Grid
!-------------------------------------------------------------------------------
  SUBROUTINE SetGridCell(GM,ExtraEll_O)
    TYPE(CRDS)                     :: GM
    INTEGER,OPTIONAL               :: ExtraEll_O
    INTEGER                        :: I,ExtraEll       
    REAL(DOUBLE)                   :: Radd,BoxExtent,Coef
    REAL(DOUBLE),DIMENSION(3,3)    :: ABC
!
    IF(GM%PBC%Dimen == 0) THEN
       CALL New_CellSet(CS_Grid,1,3)
       CS_Grid%CellCarts%D = Zero
       RETURN
    ENDIF
!
    IF(Present(ExtraEll_O)) THEN
       ExtraEll = ExtraEll_O
    ELSE
       ExtraEll = 0
    ENDIF
!
    ABC(:,1) = GM%PBC%BoxShape(:,1)+GM%PBC%BoxShape(:,2)+GM%PBC%BoxShape(:,3)
    ABC(:,2) = GM%PBC%BoxShape(:,1)+GM%PBC%BoxShape(:,2)-GM%PBC%BoxShape(:,3)
    ABC(:,3) = GM%PBC%BoxShape(:,1)-GM%PBC%BoxShape(:,2)-GM%PBC%BoxShape(:,3)
    Radd = Zero
    DO I=1,3 
       IF(GM%PBC%AutoW(I)) THEN
          Radd = MAX(Radd,SQRT(ABC(1,I)**2+ABC(2,I)**2+ABC(3,I)**2))
       ENDIF
    ENDDO
!
    Coef = One
    BoxExtent=AtomPairExtent(MaxRadialAngSym+ExtraEll,Two*MinRadialExponent,Thresholds%Dist,Coef)
    CALL New_CellSet_Sphere(CS_Grid,GM%PBC%AutoW,GM%PBC%BoxShape,Radd+BoxExtent)
!
  END SUBROUTINE SetGridCell
!-------------------------------------------------------------------------------
! Set Up the Double Cell Set Sum for Kxc
!-------------------------------------------------------------------------------
  SUBROUTINE SetKxcCell(GM,Ell,ZetaA,ZetaB,ExtraEll_O)
    TYPE(CRDS)                     :: GM
    TYPE(CellSet)                  :: CSTemp
    INTEGER,OPTIONAL               :: ExtraEll_O
    INTEGER                        :: I,J,K,Ell,NC1,NC2,Icount,ExtraEll
    REAL(DOUBLE)                   :: Radius,Radd,A0,B0,C0,AtmExtent,BoxExtent
    REAL(DOUBLE)                   :: ZetaA,ZetaB,ZA,ZB,Xi,Zeta,R12,R22,Coef
    REAL(DOUBLE),DIMENSION(3)      :: R1,R2
    REAL(DOUBLE),DIMENSION(3,3)    :: ABC
!
    IF(GM%PBC%Dimen == 0) THEN
       CALL New_CellSet(CS_Kxc,1,6)
       CS_Kxc%CellCarts%D=Zero
       RETURN
    ENDIF
!
    IF(Present(ExtraEll_O)) THEN
       ExtraEll = ExtraEll_O
    ELSE
       ExtraEll = 0
    ENDIF
!
    ABC(:,1) = GM%PBC%BoxShape(:,1)+GM%PBC%BoxShape(:,2)+GM%PBC%BoxShape(:,3)
    ABC(:,2) = GM%PBC%BoxShape(:,1)+GM%PBC%BoxShape(:,2)-GM%PBC%BoxShape(:,3)
    ABC(:,3) = GM%PBC%BoxShape(:,1)-GM%PBC%BoxShape(:,2)-GM%PBC%BoxShape(:,3)
    Radd = Zero
    DO I=1,3
       IF(GM%PBC%AutoW(I)) THEN
          Radd = MAX(Radd,SQRT(ABC(1,I)**2+ABC(2,I)**2+ABC(3,I)**2))
       ENDIF
    ENDDO
!
    Zeta = ZetaA+ZetaB
    Xi   = ZetaA*ZetaB/Zeta
    Coef = One
    AtmExtent = AtomPairExtent(Ell+ExtraEll,Xi,Thresholds%Dist,Coef)
    BoxExtent = AtomPairExtent(Ell+ExtraEll,Zeta,Thresholds%Dist,Coef)
    ZA = ZetaA/Zeta
    ZB = ZetaB/Zeta
!
    CALL New_CellSet_Sphere(CSTemp,GM%PBC%AutoW,GM%PBC%BoxShape,Radd+AtmExtent)
!
    ICount = 0
    DO NC1 = 1,CSTemp%NCells
       DO NC2 = 1,CSTemp%NCells
          R1(1) = ZA*CSTemp%CellCarts%D(1,NC1)+ZB*CSTemp%CellCarts%D(1,NC2)
          R1(2) = ZA*CSTemp%CellCarts%D(2,NC1)+ZB*CSTemp%CellCarts%D(2,NC2)
          R1(3) = ZA*CSTemp%CellCarts%D(3,NC1)+ZB*CSTemp%CellCarts%D(3,NC2)
          R2(1) = CSTemp%CellCarts%D(1,NC1)-CSTemp%CellCarts%D(1,NC2)
          R2(2) = CSTemp%CellCarts%D(2,NC1)-CSTemp%CellCarts%D(2,NC2)
          R2(3) = CSTemp%CellCarts%D(3,NC1)-CSTemp%CellCarts%D(3,NC2)
          R12  = SQRT(R1(1)**2+R1(2)**2+R1(3)**2)
          R22  = SQRT(R2(1)**2+R2(2)**2+R2(3)**2)
          IF(R12 < BoxExtent+Radd) THEN
             IF(R22 < AtmExtent) THEN
                ICount = ICount+1
             ENDIF
          ENDIF
       ENDDO
    ENDDO
!
    CALL New_CellSet(CS_Kxc,ICount,6)
!
    ICount = 0
    DO NC1 = 1,CSTemp%NCells
       DO NC2 = 1,CSTemp%NCells
          R1(1) = ZA*CSTemp%CellCarts%D(1,NC1)+ZB*CSTemp%CellCarts%D(1,NC2)
          R1(2) = ZA*CSTemp%CellCarts%D(2,NC1)+ZB*CSTemp%CellCarts%D(2,NC2)
          R1(3) = ZA*CSTemp%CellCarts%D(3,NC1)+ZB*CSTemp%CellCarts%D(3,NC2)
          R2(1) = CSTemp%CellCarts%D(1,NC1)-CSTemp%CellCarts%D(1,NC2)
          R2(2) = CSTemp%CellCarts%D(2,NC1)-CSTemp%CellCarts%D(2,NC2)
          R2(3) = CSTemp%CellCarts%D(3,NC1)-CSTemp%CellCarts%D(3,NC2)
          R12  = SQRT(R1(1)**2+R1(2)**2+R1(3)**2)
          R22  = SQRT(R2(1)**2+R2(2)**2+R2(3)**2)
          IF(R12 < BoxExtent+Radd) THEN
             IF(R22 < AtmExtent) THEN
                ICount = ICount+1
                CS_Kxc%CellCarts%D(1,ICount) = CSTemp%CellCarts%D(1,NC1)
                CS_Kxc%CellCarts%D(2,ICount) = CSTemp%CellCarts%D(2,NC1)
                CS_Kxc%CellCarts%D(3,ICount) = CSTemp%CellCarts%D(3,NC1)
                CS_Kxc%CellCarts%D(4,ICount) = CSTemp%CellCarts%D(1,NC2)
                CS_Kxc%CellCarts%D(5,ICount) = CSTemp%CellCarts%D(2,NC2)
                CS_Kxc%CellCarts%D(6,ICount) = CSTemp%CellCarts%D(3,NC2)
             ENDIF
          ENDIF
       ENDDO
    ENDDO
!
  END SUBROUTINE SetKxcCell
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

    DO I=1,3
       IF(GM%PBC%AutoW(I)) THEN
          DO 
             IF(VecF(I) > One) THEN
                VecF(I) = VecF(I) - One
             ELSE
                EXIT
             ENDIF
          ENDDO
          DO 
             IF(VecF(I) < Zero) THEN
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

    VecF = AtomToFrac(GM,VecA)
    CALL FracCyclic(GM,VecF)
    VecA = FracToAtom(GM,VecF)

  END SUBROUTINE AtomCyclic
!===================================================================================
!
!===================================================================================
  SUBROUTINE  CalFracCarts(GM)
    TYPE(CRDS)                 :: GM
    INTEGER                    :: I
!
!   Generate the Fractioanl Coordinates
!
    DO I=1,GM%NAtms
       GM%BoxCarts%D(:,I) = AtomToFrac(GM,GM%Carts%D(:,I))
       GM%BoxVects%D(:,I) = AtomToFrac(GM,GM%Vects%D(:,I))
    ENDDO
!
  END SUBROUTINE CalFracCarts
!===================================================================================
!
!===================================================================================
  SUBROUTINE  CalAtomCarts(GM)
    TYPE(CRDS)                 :: GM
    INTEGER                    :: I
!
!   Generate the Atomic Coordinates
!
    DO I=1,GM%NAtms
       GM%Carts%D(:,I)   = FracToAtom(GM,GM%BoxCarts%D(:,I))
       GM%Vects%D(:,I)   = FracToAtom(GM,GM%BoxVects%D(:,I))
    ENDDO
!
  END SUBROUTINE CalAtomCarts
!===================================================================================
!
!===================================================================================
  SUBROUTINE  WrapAtoms(GM)
    TYPE(CRDS)     :: GM
    INTEGER        :: I
!
    CALL CalFracCarts(GM)
    IF(GM%PBC%AtomW) THEN
       DO I=1,GM%NAtms
          CALL FracCyclic(GM,GM%BoxCarts%D(:,I))
       ENDDO
    ENDIF
    CALL CalAtomCarts(GM)
!
  END SUBROUTINE WrapAtoms
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
       GM%Carts%D(:,I)    = GM%Carts%D(:,I) + ATvec(:)
       GM%BoxCarts%D(:,I) = GM%BoxCarts%D(:,I) + FTvec(:)
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
!
#endif
!
END MODULE AtomPairs
!!$!===================================================================================
!!$!
!!$!===================================================================================
!!$      SUBROUTINE VecFormatToAngFormat(GM,A,B,C,Alpha,Beta,Gamma)
!!$        TYPE(CRDS)                  :: GM
!!$        REAL(DOUBLE)                :: A,B,C,Alpha,Beta,Gamma
!!$        REAL(DOUBLE),PARAMETER      :: DegToRad =  1.745329251994329576923D-2
!!$!
!!$        A = SQRT(GM%PBC%BoxShape(1,1)**2 + GM%PBC%BoxShape(2,1)**2+ GM%PBC%BoxShape(3,1)**2)
!!$        B = SQRT(GM%PBC%BoxShape(1,2)**2 + GM%PBC%BoxShape(2,2)**2+ GM%PBC%BoxShape(3,2)**2)
!!$        C = SQRT(GM%PBC%BoxShape(1,3)**2 + GM%PBC%BoxShape(2,3)**2+ GM%PBC%BoxShape(3,3)**2)
!!$!
!!$        Gamma = ACOS((GM%PBC%BoxShape(1,1)*GM%PBC%BoxShape(1,2))/(A*B))/DegToRad
!!$        Beta  = ACOS((GM%PBC%BoxShape(1,1)*GM%PBC%BoxShape(1,3))/(A*C))/DegToRad
!!$        Alpha = GM%PBC%BoxShape(1,2)*GM%PBC%BoxShape(1,3)+GM%PBC%BoxShape(2,2)*GM%PBC%BoxShape(2,3)   
!!$        Alpha = ACOS(Alpha/(B*C))/DegToRad
!!$!
!!$      END SUBROUTINE VecFormatToAngFormat
