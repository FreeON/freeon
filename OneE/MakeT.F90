!    COMPUTE THE KINETIC ENERGY MATRIX
!    Author: Matt Challacombe and C.J. Tymczak
!----------------------------------------------------------------------
PROGRAM MakeT
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
  USE TBlock
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef PARALLEL
  TYPE(DBCSR)         :: T,T2
#else
  TYPE(BCSR)          :: T,T2
#endif
#ifdef PERIODIC 
  INTEGER                   :: NC
  REAL(DOUBLE),DIMENSION(3) :: B
#endif
  TYPE(AtomPair)      :: Pair
!
  TYPE(BSET)          :: BS
  TYPE(CRDS)          :: GM
  TYPE(DBL_RNK4)      :: MD
!
  TYPE(ARGMT)         :: Args
  INTEGER             :: P,R,AtA,AtB,NN                          
  CHARACTER(LEN=5),PARAMETER :: Prog='MakeT'
!--------------------------------------- 
! Start up macro
!
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
!----------------------------------------------
! Get basis set and geometry
!
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
!---------------------------------------------- 
! Allocations 
!
  CALL New(MD,(/3,BS%NASym+1,BS%NASym+1,2*BS%NASym+2/),(/1,-1,-1,-1/))
  CALL New(T)
#ifdef PERIODIC
!-----------------------------------------------
! Calculate the Number of Cells
!
  CALL SetCellNumber(GM)
  CALL PPrint(CS_OUT,'outer sum',Prog)
#endif
!-----------------------------------------------
! Initialize the matrix and associated indecies
!
  P=1; R=1; T%RowPt%I(1)=1
  CALL SetEq(T%MTrix,Zero)
!-----------------------------------------------
! Main loops
!
#ifdef PARALLEL
  T%NAtms=0
  DO AtA=Beg%I(MyId),End%I(MyId)
     T%NAtms=T%NAtms+1  
#else
  T%NAtms=NAtoms
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
                 T%MTrix%D(R:R+NN-1)=T%MTrix%D(R:R+NN-1)+TBlok(BS,MD,Pair)
              ENDIF
           ENDDO
#else
           T%MTrix%D(R:R+NN-1)=TBlok(BS,MD,Pair)
#endif
           T%ColPt%I(P)=AtB
           T%BlkPt%I(P)=R
           R=R+NN 
           P=P+1 
#ifdef PARALLEL           
           T%RowPt%I(T%NAtms+1)=P
           IF(R>MaxNon0Node.OR.P>MaxBlksNode)THEN
              WRITE(*,*)' MyId = ',MyId,MaxBlksNode,MaxNon0Node
              WRITE(*,*)' MyId = ',MyId,' P = ',P,' R = ',R
              CALL Halt(' DBCSR dimensions blown in MakeT ')
           ENDIF
#else  
           T%RowPt%I(AtA+1)=P        
           IF(R>MaxNon0.OR.P>MaxBlks) &
                CALL Halt(' BCSR dimensions blown in MakeT ')
#endif
        ENDIF
     ENDDO
  ENDDO
  T%NBlks=P-1
  T%NNon0=R-1
!------------------------------------------------------------
! Put T to disk
!  
  CALL Filter(T2,T)
  CALL Put(T2,TrixFile('T',Args))
!-----------------------------------------------------------
! Printing
!
  CALL PChkSum(T2,'T',Prog)
  CALL PPrint( T2,'T')
  CALL Plot(   T2,'T')
!------------------------------------------------------------
! Tidy up
! 
  CALL Delete(T)
  CALL Delete(T2)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(MD)
! didn't count flops, any accumulation is residual
! from matrix routines
  PerfMon%FLOP=Zero 
  CALL ShutDown(Prog)
!
END PROGRAM MakeT
