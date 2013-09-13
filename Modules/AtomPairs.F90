!------------------------------------------------------------------------------
!    This code is part of the FreeON suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the FreeON group for review, and possible
!    dissemination in future releases.
!------------------------------------------------------------------------------
!    Compute the Atom Pair and the Cells for PBC
!--  Matt Challacombe and  C. J. Tymczak
MODULE AtomPairs
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE GlobalObjects
  USE BoundingBox
  USE Thresholding
  USE PBC
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

!!========================================================================
!!  OMG WTF IS THIS KRAP?? AND WHY IS IT HERE?
   Pair%AB2 = MinImageDist(GM,GM%Carts%D(1:3,I),GM%Carts%D(1:3,J))

!!$!!
!!$!!========================================================================

    Pair%KA   = GM%AtTyp%I(I)
    Pair%KB   = GM%AtTyp%I(J)
    Pair%NA   = BS%BFKnd%I(Pair%KA)
    Pair%NB   = BS%BFKnd%I(Pair%KB)
    Pair%A(:) = GM%Carts%D(:,I)
    Pair%B(:) = GM%Carts%D(:,J)


!!$    WRITE(*,*)' ------------------------------------------------'
!!$    WRITE(*,*)' A = ',Pair%A(:)
!!$    WRITE(*,*)' B = ',Pair%B(:)
!!$    WRITE(*,*)' AB2 = ',Pair%AB2
!!$
    IF(I==J) THEN
       Pair%SameAtom = .TRUE.
    ELSE
       Pair%SameAtom = .FALSE.
    ENDIF
!
    SetAtomPair=TestAtomPair(Pair,Box_O)


    Pair%AB2  = DOT_PRODUCT(GM%Carts%D(1:3,I)-GM%Carts%D(1:3,J),GM%Carts%D(1:3,I)-GM%Carts%D(1:3,J))

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
! Calculate the Minmiun Image Distance Between Atoms I and J
!-------------------------------------------------------------------------------
  FUNCTION MinImageDist(GM,Carts1,Carts2)
    TYPE(CRDS)                :: GM
    INTEGER                   :: K,N1,N2,N3
    REAL(DOUBLE)              :: MinImageDist
    INTEGER,DIMENSION(3)      :: NRgn
    REAL(DOUBLE),DIMENSION(3) :: Carts1,Carts2
    REAL(DOUBLE),DIMENSION(3) :: VecF,VecA

    NRgn = 0
    DO K=1,3
       IF(GM%PBC%AutoW%I(K)==1) NRgn(K) = 1
    ENDDO

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
!!$!-------------------------------------------------------------------------------
!!$! Make GM Periodic
!!$!-------------------------------------------------------------------------------
!!$  SUBROUTINE PBCInfoFromNewCarts(PBC,WP_O)
!!$    TYPE(PBCInfo)             :: PBC
!!$    LOGICAL,DIMENSION(3)      :: WP
!!$    LOGICAL,DIMENSION(3),OPTIONAL      :: WP_O
!!$!
!!$!   This routine rebuilds PBC data based on PBC%BoxShape
!!$!
!!$    PBC%CellVolume   = CellVolume(PBC%BoxShape%D,PBC%AutoW%I)
!!$    PBC%CellCenter%D = CellCenter(PBC%BoxShape%D,PBC%AutoW%I)
!!$    PBC%InvBoxSh%D   = InverseMatrix(PBC%BoxShape%D)
!!$!
!!$    IF(PRESENT(WP_O)) THEN
!!$       WP = WP_O
!!$    ELSE
!!$       WP = (/.true.,.true.,.true./)
!!$    ENDIF
!!$!
!!$    IF(WP(1)) THEN
!!$       PBC%CellVolume = CellVolume(PBC%BoxShape%D,PBC%AutoW%I)
!!$!      Calculate the Dipole and Quadripole Factors
!!$       IF(PBC%Dimen < 2) THEN
!!$          PBC%DipoleFAC = Zero
!!$          PBC%QupoleFAC = Zero
!!$       ELSEIF(PBC%Dimen ==2) THEN
!!$          PBC%DipoleFAC = (Four*Pi/PBC%CellVolume)*(One/(PBC%Epsilon+One))
!!$          PBC%QupoleFAC =  Zero
!!$          IF(ABS(PBC%DipoleFAC) .LT. 1.D-14) PBC%DipoleFAC = Zero
!!$       ELSEIF(PBC%Dimen ==3) THEN
!!$          PBC%DipoleFAC = -(Four*Pi/PBC%CellVolume)*(One/Three - One/(Two*PBC%Epsilon+One))
!!$          PBC%QupoleFAC =  ( Two*Pi/PBC%CellVolume)*(One/Three - One/(Two*PBC%Epsilon+One))
!!$          IF(ABS(PBC%DipoleFAC) .LT. 1.D-14) PBC%DipoleFAC = Zero
!!$          IF(ABS(PBC%QupoleFAC) .LT. 1.D-14) PBC%QupoleFAC = Zero
!!$       ENDIF
!!$    ENDIF
!!$!   Calculate The Inverse of BoxShape  InvBoxSh = [BoxShape]^(-1)
!!$    IF(WP(2)) THEN
!!$       PBC%InvBoxSh%D = InverseMatrix(PBC%BoxShape%D)
!!$    ENDIF
!!$!
!!$  END SUBROUTINE PBCInfoFromNewCarts

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
