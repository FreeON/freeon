!    SPARSE UPDATE OF THE BFGS HESSIAN 
!    Author: Matt Challacombe
!------------------------------------------------------------------------------   
PROGRAM BFGSHess
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalCharacters
   USE InOut
   USE PrettyPrint
   USE MemMan
   USE Parse
   USE Macros
   USE SetXYZ
   USE LinAlg
   USE Functionals
   USE MatFunk
#ifdef PARALLEL
   USE MondoMPI
#endif
   USE BFGSMisc
   IMPLICIT NONE
   TYPE(ARGMT)                :: Args
   CHARACTER(LEN=8),PARAMETER :: Prog='BFGSHess'
   TYPE(AtomPair)             :: Pair
   TYPE(BSet)                 :: BS
   TYPE(CRDS)                 :: GM0,GM1
   TYPE(BCSR)                 :: B0,B1,sDX,sDXB
   TYPE(DBL_VECT)             :: G0,G1,X0,X1,DX,DG,DXB
   TYPE(DBL_RNK2)             :: X, G, DnsB,BInv
   INTEGER                    :: N3,I,J,K,IGeom,M1,M2,M3,I1,I2
   REAL(DOUBLE)               :: D1,D2,F,AdEl,MinEl,A, &
                                 AB2,AB2Max,GradEDotDeltaX
   REAL(DOUBLE),DIMENSION(6)  :: AA
#ifdef PERIODIC
    REAL(DOUBLE),DIMENSION(3) :: VecF,VecA,NRgn
#endif
!----------------------------------------------------------------------------
   CALL StartUp(Args,Prog)
   IGeom=Args%I%I(3)
!  Override default blocking...
   N3=3*NAtoms
   NBasF=N3

   MaxAtms=NAtoms+1
   MaxBlks=NAtoms**2+1
   MaxNon0=N3**2+1

   PrintFlags%Mat=DEBUG_MATRICES
   BSiz%I=3
   OffS%I(1)=1
   DO I=2,NAtoms
      OffS%I(I)=OffS%I(I-1)+3
   ENDDO      
!  Allocations (probably a lot more than nessesary)
   CALL New(G0,N3)
   CALL New(G1,N3)
   CALL New(X0,N3)
   CALL New(X1,N3)
   CALL New(DG,N3)
   CALL New(DX,N3)
   CALL New(DXB,N3)
   CALL New(B0)
   CALL New(B1)
   IF(IGeom==1)THEN
!     Bo=I
      CALL SetToI(B1)
      CALL Put(1.D-10,'GradEDotDeltaX')
      GO TO 10
   ENDIF
   CALL New(sDX)
   CALL New(sDXB)
!  Get position and gradient info
   CALL Get(BS,Tag_O=CurBase)
   CALL Get(GM0,Tag_O=PrvGeom)
!   WRITE(*,*)' BFGSHess : OLDGEOM # ',PrvGeom
!   CALL PPrint(GM0,Unit_O=6,PrintGeom_O='XYZ')
   CALL Get(GM1,Tag_O=CurGeom)
!   WRITE(*,*)' BFGSHess : NEWGEOM # ',CurGeom
!   CALL PPrint(GM1,Unit_O=6,PrintGeom_O='XYZ')
   CALL Get(G0,'GradE',Tag_O=PrvGeom)
!   WRITE(*,*)' BFGSHess : OLDGRAD # ',PrvGeom
   CALL Get(G1,'GradE',Tag_O=CurGeom)
!   WRITE(*,*)' BFGSHess : NEWGRAD # ',CurGeom
#ifdef PERIODIC
    NRgn=0
    DO K=1,3
      IF(GM1%PBC%AutoW(K))NRgn(K)=1
    ENDDO
#endif
   K=0
   DO I=1,NAtoms
      I1=(I-1)*3+1
      I2=I*3
      X0%D(I1:I2)=GM0%Carts%D(:,I)
      X1%D(I1:I2)=GM1%Carts%D(:,I)
#ifdef PERIODIC
!     Find the minimum image distance
      AB2MAX = 1D8
      DO M1=-NRgn(1),NRgn(1)
         DO M2=-NRgn(2),NRgn(2)
            DO M3=-NRgn(3),NRgn(3)  
               VecF=(/DBLE(M1),DBLE(M2),DBLE(M3)/)
               VecA=ABS(FracToAtom(GM1,VecF)+X0%D(I1:I2)-X1%D(I1:I2))
               AB2=VecA(1)**2+VecA(2)**2+VecA(3)**2
               IF(AB2<AB2MAX)THEN
                  DX%D(I1:I2)=VecA
                  AB2MAX=AB2
                ENDIF
            ENDDO
         ENDDO
      ENDDO
#else
      DX%D(I1:I2)=X1%D(I1:I2)-X0%D(I1:I2)
#endif
   ENDDO
!WRITE(*,*)' DX = ',DX%D
   DG%D=G1%D-G0%D
   GradEDotDeltaX=DOT_PRODUCT(G1%D,DX%D)/DBLE(N3)
   CALL Put(GradEDotDeltaX,'GradEDotDeltaX')
!  Get previous steps BFGS approx to Hessian
   CALL Get(B0,TrixFile('B',Args,NoTags_O=.TRUE.))
!  Compute intermediate scalars and vectors
   D1=DOT_PRODUCT(DX%D,DG%D)             !  D1=[/\G^T].[/\X]
   CALL Set_BCSR_EQ_VECT(sDX,DX)         !  sDX=SparseM[DX]
   CALL Multiply(B0,sDX,sDXB)            !  sDXB=[/\X^T].[B_I]
   CALL Set_VECT_EQ_BCSR(DXB,sDXB)       !  DXB=DenseV[sDXB]
   D2=DOT_PRODUCT(DX%D,DXB%D)            !  D2=[/\X^T].[B_I].[/\X]
   MinEl=MinDiagEl(B0,DXB,DG,D1,D2)      !  Min diag element of [B_(I+1)]
   Thresholds%Trix=1.D-1*MinEl
!  Sparse BFGS update 
   CALL ReUpB(GM1,BS,B0,DG,DXB,D1,D2,B1) ! [B_I+1]=BFGS([B_I],[/\X],[/\G])
10 CONTINUE
!   CALL PPrint(B1,'B1',Unit_O=6)
!  Filter out small blocks
   CALL Filter(B0,B1)
   CALL Put(B0,TrixFile('B',Args,NoTags_O=.TRUE.))
   CALL Delete(G0)
   CALL Delete(G1)
   CALL Delete(X0)
   CALL Delete(X1)
   CALL Delete(DG)
   CALL Delete(DX)
   CALL Delete(DXB)
   CALL Delete(B0)
   CALL Delete(B1)
   CALL Shutdown(Prog)
END PROGRAM 
