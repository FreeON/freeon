!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe and  C. J. Tymczak
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!    Compute the Atom Pair and the Cells for PBC
!
MODULE AtomPairs
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE GlobalObjects
  USE BoundingBox
#ifdef PERIODIC
  USE CellSets
#endif
  USE Thresholding
  IMPLICIT NONE
!---------------------------------------------------------------------------
! Global variables
!---------------------------------------------------------------------------
#ifdef PERIODIC
!  LOGICAL,DIMENSION(3)        :: PBCXYZ     
!  REAL(DOUBLE),DIMENSION(3,3) :: LATVEC     
  TYPE(CellSet)               :: CS,CSMM1,CSMM2
#endif
CONTAINS
!--------------------------------------------------------------------------
! Test to see if the atomic overlap is above tolerance and create pair
!--------------------------------------------------------------------------
  FUNCTION SetAtomPair(GM,BS,I,J,Pair)
    LOGICAL                   :: SetAtomPair
    INTEGER                   :: I,J
    TYPE(AtomPair)            :: Pair
    TYPE(CRDS)                :: GM
    TYPE(BSet)                :: BS
#ifdef PERIODIC
    REAL(DOUBLE)              :: LX
    REAL(DOUBLE)              :: LY
    REAL(DOUBLE)              :: LZ
    LX = ABS(GM%BoxCarts%D(1,I)-GM%BoxCarts%D(1,J))
    LY = ABS(GM%BoxCarts%D(2,I)-GM%BoxCarts%D(2,J))
    LZ = ABS(GM%BoxCarts%D(3,I)-GM%BoxCarts%D(3,J))
    IF(GM%AutoW(1) .AND. LX > half) LX = one-LX
    IF(GM%AutoW(2) .AND. LY > half) LY = one-LY    
    IF(GM%AutoW(3) .AND. LZ > half) LZ = one-LZ
    LX = LX*GM%BoxShape%D(1,1)
    LY = LX*GM%BoxShape%D(1,2)+LY*GM%BoxShape%D(2,2)
    LZ = LX*GM%BoxShape%D(1,3)+LY*GM%BoxShape%D(2,3)+LZ*GM%BoxShape%D(3,3)
    Pair%AB2 = LX*LX+LY*LY+LZ*LZ
#else
    Pair%AB2 = (GM%Carts%D(1,I)-GM%Carts%D(1,J))**2 &
             + (GM%Carts%D(2,I)-GM%Carts%D(2,J))**2 &
             + (GM%Carts%D(3,I)-GM%Carts%D(3,J))**2 
#endif
    SetAtomPair=TestAtomPair(Pair)
    IF(SetAtomPair) THEN
       Pair%A(:) = GM%Carts%D(:,I) 
       Pair%B(:) = GM%Carts%D(:,J) 
       Pair%KA   = GM%AtTyp%I(I)
       Pair%KB   = GM%AtTyp%I(J)
       Pair%NA   = BS%BFKnd%I(Pair%KA)
       Pair%NB   = BS%BFKnd%I(Pair%KB)
    ENDIF
    IF(I==J) THEN 
       Pair%SameAtom = .TRUE.
    ELSE
       Pair%SameAtom = .FALSE.
    ENDIF
  END FUNCTION SetAtomPair
!--------------------------------------------------------------------------
! Convert a Block Matrix to a Vector
!--------------------------------------------------------------------------
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
!--------------------------------------------------------------------------
! Convert a Block Matrix to a Vector
!--------------------------------------------------------------------------
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
#ifdef PERIODIC
!--------------------------------------------------------------------------
! Set up the lattice vecters to sum over 
! Also, Added a Factor that Takes into Account the Box Size and Shape
!--------------------------------------------------------------------------
  SUBROUTINE SetCellNumber(GM)
    TYPE(CRDS)                     :: GM
    INTEGER                        :: I
    REAL(DOUBLE)                   :: Radius,Radd,A0,B0,C0
!
    A0 = Zero
    B0 = Zero
    C0 = Zero
    DO I=1,3
       A0 = A0 + (GM%BoxShape%D(I,1)+GM%BoxShape%D(I,2)+GM%BoxShape%D(I,3))**2
       B0 = B0 + (GM%BoxShape%D(I,1)+GM%BoxShape%D(I,2)-GM%BoxShape%D(I,3))**2
       C0 = C0 + (GM%BoxShape%D(I,1)-GM%BoxShape%D(I,2)-GM%BoxShape%D(I,3))**2
    ENDDO
    Radd = MAX(SQRT(A0),SQRT(B0))
    Radd = MAX(Radd,SQRT(C0))
!
    Radius = Radd+SQRT(AtomPairDistanceThreshold)
    CALL New_CellSet_Sphere(CS,GM%AutoW,GM%BoxShape%D,Radius)
!
    WRITE(*,*) 'THE NUMBER OF CELLS  in CS = ',CS%NCells  
!
  END SUBROUTINE SetCellNumber
!----------------------------------------------------------------------------
!     Convert from Atomic Coordinates  to Fractional Coordinates
!----------------------------------------------------------------------------
  FUNCTION AtomToFrac(GM,VecA) RESULT(VecF)
    TYPE(CRDS)                 :: GM
    REAL(DOUBLE),DIMENSION(3)  :: VecA,VecF
!
    VecF(1) = VecA(1)*GM%InvBoxSh%D(1,1)
    VecF(2) = VecA(1)*GM%InvBoxSh%D(1,2)+VecA(2)*GM%InvBoxSh%D(2,2)
    VecF(3) = VecA(1)*GM%InvBoxSh%D(1,3)+VecA(2)*GM%InvBoxSh%D(2,3) &
         + VecA(3)*GM%InvBoxSh%D(3,3)
!
  END FUNCTION AtomToFrac
!-----------------------------------------------------------------------------
!     Convert from Fractional Coordinates to Atomic Coordinates
!-----------------------------------------------------------------------------
  FUNCTION FracToAtom(GM,VecF) RESULT(VecA)
    TYPE(CRDS)                 :: GM
    REAL(DOUBLE),DIMENSION(3)  :: VecA,VecF
!
    VecA(1) = VecF(1)*GM%BoxShape%D(1,1)
    VecA(2) = VecF(1)*GM%BoxShape%D(1,2)+VecF(2)*GM%BoxShape%D(2,2)
    VecA(3) = VecF(1)*GM%BoxShape%D(1,3)+VecF(2)*GM%BoxShape%D(2,3) &
         + VecF(3)*GM%BoxShape%D(3,3)
!
  END FUNCTION FracToAtom
!-----------------------------------------------------------------------------
!     In Perodic Systems, Cyclically put positions back in the Fracional Box
!-----------------------------------------------------------------------------
  SUBROUTINE FracCyclic(GM,VecF)
    TYPE(CRDS)                 :: GM        
    REAL(DOUBLE),DIMENSION(3)  :: VecF
    INTEGER                    :: I,N      
!
    DO I=1,3
       IF(GM%AutoW(I)) THEN
          IF(VecF(I) > one) THEN
             N = INT(VecF(I))
             VecF(I) = VecF(I)-DBLE(N)
          ENDIF
          IF(VecF(I) < zero) THEN
             N = INT(-VecF(I)+one)
             VecF(I) = VecF(I)+DBLE(N)
          ENDIF
       ENDIF
    ENDDO
  END SUBROUTINE FracCyclic
!------------------------------------------------------------------------------
!     Test to See if in the Box (Fractional Coordinates)
!------------------------------------------------------------------------------
  FUNCTION InFracBox(GM,VecF)
    TYPE(CRDS)                 :: GM     
    LOGICAL                    :: InFracBox
    INTEGER                    :: I
    REAL(DOUBLE),DIMENSION(3)  :: VecF
!
    InFracBox = .TRUE.
    DO I=1,3
       IF(GM%AutoW(I)) THEN
          IF(VecF(I) < zero .OR. VecF(I) > one) THEN
             InFracBox = .FALSE.
             RETURN
          ENDIF
       ENDIF
    ENDDO
!
  END FUNCTION InFracBox
#endif
!
END MODULE AtomPairs
