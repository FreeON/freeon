!    COMPUTE THE EFFECTIVE CORE POTENTIAL MATRIX U
!    Author: Matt Challacombe 
!------------------------------------------------------------
PROGRAM MakeU
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters  
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE AtomPairs
  USE ECPBlock
  USE BraBloks
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef PARALLEL
  TYPE(DBCSR)                :: U,T1
#else
  TYPE(BCSR)                 :: U,T1
#endif
#ifdef PERIODIC 
  INTEGER                    :: NC
  REAL(DOUBLE),DIMENSION(3)  :: B
#endif
  TYPE(AtomPair)             :: Pair
  TYPE(BSET)                 :: BS
  TYPE(CRDS)                 :: GM
  TYPE(DBL_RNK4)             :: MD
  TYPE(ARGMT)                :: Args
  INTEGER                    :: P,R,AtA,AtB,NN                         
  CHARACTER(LEN=5),PARAMETER :: Prog='MakeU'
  !--------------------------------------- 
  ! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
  !----------------------------------------------
  ! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
  !---------------------------------------------- 
  ! Allocations 
  CALL NewBraBlok(BS)
  CALL New(U)
  !-----------------------------------------------
  ! Initialize the matrix and associated indecies
  P=1; R=1; U%RowPt%I(1)=1
  CALL SetEq(U%MTrix,Zero)
  !-----------------------------------------------
  ! Main loops
#ifdef PARALLEL
  U%NAtms=0
  DO AtA=Beg%I(MyId),End%I(MyId)
     U%NAtms=U%NAtms+1  
#else
     U%NAtms=NAtoms
     DO AtA=1,NAtoms
#endif
        DO AtB=1,NAtoms
           IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
              NN = Pair%NA*Pair%NB
#ifdef PERIODIC
              B = Pair%B
              DO NC = 1,CS_OUT%NCells
                 Pair%B = B+CS_OUT%CellCarts%D(:,NC)
                 Pair%AB2  = (Pair%A(1)-Pair%B(1))**2 &
                      + (Pair%A(2)-Pair%B(2))**2 &
                      + (Pair%A(3)-Pair%B(3))**2
                 IF(TestAtomPair(Pair)) THEN
                    U%MTrix%D(R:R+NN-1)=U%MTrix%D(R:R+NN-1)+UBlok(BS,GM,Pair)
                 ENDIF
              ENDDO
#else
              U%MTrix%D(R:R+NN-1)=SBlok(BS,Pair)
#endif
              U%ColPt%I(P)=AtB
              U%BlkPt%I(P)=R
              R=R+NN 
              P=P+1 
#ifdef PARALLEL           
              U%RowPt%I(U%NAtms+1)=P
              IF(R>MaxNon0Node.OR.P>MaxBlksNode)THEN
                 WRITE(*,*)' MyId = ',MyId,MaxBlksNode,MaxNon0Node
                 WRITE(*,*)' MyId = ',MyId,' P = ',P,' R = ',R
                 CALL Halt(' DBCSR dimensions blown in MakeS ')
              ENDIF
#else  
              U%RowPt%I(AtA+1)=P        
              IF(R>MaxNon0.OR.P>MaxBlks) &
                   CALL Halt(' BCSR dimensions blown in MakeS ')
#endif
           ENDIF
        ENDDO
     ENDDO
     U%NBlks=P-1
     U%NNon0=R-1
     !------------------------------------------------------------
     ! Put U to disk
     Thresholds%Trix = Thresholds%Trix*1.D-2
     CALL Filter(T1,U)
     Thresholds%Trix = Thresholds%Trix*1.D2
     CALL Put(T1,TrixFile('U',Args))
     !     CALL PPrint( T1,'U',Unit_O=6)
     !-----------------------------------------------------------
     ! Printing
     CALL PChkSum(T1,'U',Prog)
     CALL PPrint( T1,'U')
     CALL Plot(   T1,'U')
     !------------------------------------------------------------
     ! Tidy up
     CALL Delete(U)
     CALL Delete(T1)
     CALL Delete(BS)
     CALL Delete(GM)
     CALL DeleteBraBlok()
     ! didn't count flops, any accumulation is residual
     ! from matrix routines
     PerfMon%FLOP=Zero 
     CALL ShutDown(Prog)
   END PROGRAM MakeU
