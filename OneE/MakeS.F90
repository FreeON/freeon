!    COMPUTE THE OVERLAP MATRIX S
!    Authors: Matt Challacombe and C.J. Tymczak
!------------------------------------------------------------
PROGRAM MakeS
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
  USE OverlapBlock
  USE BraBloks
#ifdef PARALLEL
  USE MondoMPI
#endif
#ifdef PARALLEL_DEVELOPMENT
  USE FastMatrices
#endif
  IMPLICIT NONE
#ifdef PARALLEL
#ifdef PARALLEL_DEVELOPMENT
  TYPE(FastMat),POINTER :: S
  TYPE(BCSR)            :: T1
#else
  TYPE(DBCSR)         :: S,T1
#endif
#else
  TYPE(BCSR)          :: S,T1
#endif
#ifdef PERIODIC 
  INTEGER                   :: NC
  REAL(DOUBLE),DIMENSION(3) :: B
#endif
  TYPE(AtomPair)      :: Pair

  TYPE(BSET)          :: BS
  TYPE(CRDS)          :: GM
  TYPE(DBL_RNK4)      :: MD

  TYPE(ARGMT)         :: Args
  INTEGER             :: P,R,AtA,AtB,NN                         
  CHARACTER(LEN=5),PARAMETER :: Prog='MakeS'
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
#ifdef PARALLEL_DEVELOPMENT
  CALL New_FASTMAT(S,0,(/0,0/))
#else
  CALL New(S)
#endif
#ifdef PERIODIC
!-----------------------------------------------
! Calculate the Number of Cells

  CALL SetCellNumber(GM)
  CALL PPrint(CS_OUT,'outer sum',Prog)
#endif
!-----------------------------------------------
! Initialize the matrix and associated indecies

#ifdef PARALLEL_DEVELOPMENT
#else
  P=1; R=1; S%RowPt%I(1)=1
  CALL SetEq(S%MTrix,Zero)
#endif
!-----------------------------------------------
! Main loops

#ifdef PARALLEL
#ifdef PARALLEL_DEVELOPMENT
  DO AtA=1,NAtoms
#else
  S%NAtms=0
  DO AtA=Beg%I(MyId),End%I(MyId)
     S%NAtms=S%NAtms+1  
#endif
#else
  S%NAtms=NAtoms
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
#ifdef PARALLEL_DEVELOPMENT
                 CALL AddFASTMATBlok(S,AtA,AtB,SBlok(BS,Pair))
#else
                 S%MTrix%D(R:R+NN-1)=S%MTrix%D(R:R+NN-1)+SBlok(BS,Pair)
#endif
              ENDIF
           ENDDO
#else
#ifdef PARALLEL_DEVELOPMENT
           CALL AddFASTMATBlok(S,AtA,AtB,SBlok(BS,Pair))
#else
           S%MTrix%D(R:R+NN-1)=SBlok(BS,Pair)
#endif
#endif

#ifdef PARALLEL_DEVELOPMENT
! do nothing
#else
           S%ColPt%I(P)=AtB
           S%BlkPt%I(P)=R
           R=R+NN 
           P=P+1 
#ifdef PARALLEL           
           S%RowPt%I(S%NAtms+1)=P
           IF(R>MaxNon0Node.OR.P>MaxBlksNode)THEN
              WRITE(*,*)' MyId = ',MyId,MaxBlksNode,MaxNon0Node
              WRITE(*,*)' MyId = ',MyId,' P = ',P,' R = ',R
              CALL Halt(' DBCSR dimensions blown in MakeS ')
           ENDIF
#else  
           S%RowPt%I(AtA+1)=P        
           IF(R>MaxNon0.OR.P>MaxBlks) &
                CALL Halt(' BCSR dimensions blown in MakeS ')
#endif
#endif
        ENDIF
     ENDDO
  ENDDO
#ifdef PARALLEL_DEVELOPMENT
! do nothing
#else
  S%NBlks=P-1
  S%NNon0=R-1
#endif
!------------------------------------------------------------
! Put S to disk

#ifdef PARALLEL_DEVELOPMENT
  CALL Redistribute_FASTMAT(S)
  CALL AlignNodes()
  CALL Multiply_FASTMAT_SCALAR(S,One/DBLE(NPrc))
!  >>>>>>>>>> FilterOut is BROKEN! <<<<<<<<<<<<<
!  CALL FilterOut_FASTMAT(S)
  CALL AlignNodes()
  CALL Set_BCSR_EQ_DFASTMAT(T1,S)
!  CALL PPrint(T1,'S',Unit_O=6)
!  CALL PChkSum(T1,'S',Prog)
#else
  CALL Filter(T1,S)
  !! the following line gives a problem!!
! CALL PPrint(T1,'S',Unit_O=6)
#endif
  CALL Put(T1,TrixFile('S',Args))
!-----------------------------------------------------------
! Printing

  CALL PChkSum(T1,'S',Prog)
  CALL PPrint( T1,'S')
  CALL Plot(   T1,'S')
!------------------------------------------------------------
! Tidy up

#ifdef PARALLEL_DEVELOPMENT
  CALL Delete_FastMat1(S)
#else
  CALL Delete(S)
#endif
  CALL Delete(T1)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL DeleteBraBlok()
! didn't count flops, any accumulation is residual
! from matrix routines
  PerfMon%FLOP=Zero 
  CALL ShutDown(Prog)
END PROGRAM MakeS
