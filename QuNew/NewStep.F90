!    ULTIMATELY USE THE AINV FACTORIZED INVERSE HESSIAN TO TAKE A NEWTON STEP
!    FOR NOW, JUST DO EIGENVALUES TO INVERT THE BFGS HESSIAN 
!    Autho:  Matt Challacombe
!------------------------------------------------------------------------------
PROGRAM NewStep
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
   IMPLICIT NONE
   TYPE(ARGMT)                :: Args
   CHARACTER(LEN=7),PARAMETER :: Prog='NewStep'
   TYPE(AtomPair)             :: Pair
   TYPE(BSet)                 :: BS
   TYPE(CRDS)                 :: GM
   TYPE(BCSR)                 :: Z,ZT,sHInvG,B
   TYPE(DBL_VECT)             :: G1,X1,X2
   TYPE(DBL_RNK2)             :: X, G, DnsB,BInv
   INTEGER                    :: N3,I,J,K,IGeom,I1
   REAL(DOUBLE)               :: D1,D2,F,AdEl,MinEl,A, &
                                 StepSz,GradEDotNewStep,GMAX,GRMS
   REAL(DOUBLE),DIMENSION(18) :: NGrad
!----------------------------------------------------------------------------
   CALL StartUp(Args,Prog)
   IGeom=Args%I%I(3)
!  Override default blocking...
   N3=3*NAtoms
   NBasF=N3
   PrintFlags%Mat=DEBUG_MATRICES
   CALL SetDSYEVWork(N3) 
   BSiz%I=3
   OffS%I(1)=1
   DO I=2,NAtoms
      OffS%I(I)=OffS%I(I-1)+3
   ENDDO      
!  Allocations (probably a lot more than nessesary)
   CALL New(G1,N3)
   CALL New(X1,N3)
   CALL New(X2,N3)!
   CALL Get(GM,Tag_O=PrvGeom)
!   WRITE(*,*)' NewStep: OLD GEOMETRY # ',PrvGeom
!   CALL PPrint(GM,Unit_O=6,PrintGeom_O='XYZ')
   CALL Get(G1,'GradE',Tag_O=PrvGeom)
 !  WRITE(*,*)' NewStep: NEW GRADIENT # ',CurGeom
 !  CALL PPrint(G1,'GradE',Unit_O=6)
   K=0
   DO I=1,NAtoms
      DO J=1,3
         K=K+1
         X1%D(K)=GM%Carts%D(J,I)
      ENDDO
   ENDDO
   CALL Get(B,TrixFile('B',Args,NoTags_O=.TRUE.))
!   CALL PPrint(B,'B',Unit_O=6)
   CALL New(DnsB,(/N3,N3/))
   CALL New(BInv,(/N3,N3/))
   CALL SetEq(DnsB,B)
   CALL FunkOnSqMat(N3,Inverse,DnsB%D,BInv%D)
   X2%D=-MATMUL(BInv%D,G1%D)
!   CALL PPrint(X2,'DescentDir',Unit_O=6)
   GradEDotNewStep=DOT_PRODUCT(X2%D,G1%D)
   CALL Put(GradEDotNewStep,'GradEDotNewStep')
   CALL Get(StepSz,'StepSize')
   X2%D=X1%D+StepSz*X2%D
   K=0
   DO I=1,NAtoms
      DO J=1,3
         K=K+1
         GM%Carts%D(J,I)=X2%D(K)
      ENDDO
   ENDDO   
   CALL Put(GM,CurGeom)   
!   WRITE(*,*)' NewStep: NEW GEOMETRY # ',CurGeom
!   CALL PPrint(GM,Unit_O=6,PrintGeom_O='XYZ')
   GRMS=SQRT(DOT_PRODUCT(G1%D,G1%D))/DBLE(N3)
   GMAX=Zero
   DO I=1,N3;GMAX=MAX(GMAX,ABS(G1%D(I)));ENDDO
   CALL Put(GRMS,'RMSGrad',CurGeom)
   CALL Put(GMAX,'MaxGrad',CurGeom)
!   WRITE(*,*)NxtGeom,'GradStats ',GRMS,GMAX

   CALL Shutdown(Prog)
END PROGRAM 
