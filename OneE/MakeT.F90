!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe and  C. J. Tymczak
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!    COMPUTE THE Kintic Energy Matrix
!
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
  TYPE(DBCSR)         :: TK,T2
#else
  TYPE(BCSR)          :: TK,T2
#endif
#ifdef PERIODIC 
  LOGICAL             :: NotMakeBlock
  INTEGER             :: NC
  REAL(DOUBLE)        :: Bx,By,Bz
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
  CALL StartUp(Args,Prog)
!----------------------------------------------
! Get basis set and geometry
!
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
!---------------------------------------------- 
! Allocations 
!
  CALL New(MD,(/3,BS%NASym+1,BS%NASym+1,2*BS%NASym+2/),(/1,-1,-1,0/))
  CALL New(TK)
#ifdef PERIODIC
!-----------------------------------------------
! Calculate the Number of Cells
!
  CALL SetCellNumber(GM)
#endif
!-----------------------------------------------
! Initialize the matrix and associated indecies
!
  P=1; R=1; TK%RowPt%I(1)=1
  CALL SetEq(TK%MTrix,Zero)
!-----------------------------------------------
! Main loops
!
#ifdef PARALLEL
  TK%NAtms=0
  DO AtA=Beg%I(MyId),End%I(MyId)
     TK%NAtms=TK%NAtms+1  
#else
  TK%NAtms=NAtoms
  DO AtA=1,NAtoms
#endif
     DO AtB=1,NAtoms
        IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
           NN = Pair%NA*Pair%NB
#ifdef PERIODIC
           Bx = Pair%B(1)
           By = Pair%B(2)           
           Bz = Pair%B(3)
           NotMakeBlock = .TRUE.
           DO NC = 1,CS%NCells
              Pair%B(1) = Bx+CS%CellCarts%D(1,NC)
              Pair%B(2) = By+CS%CellCarts%D(2,NC)
              Pair%B(3) = Bz+CS%CellCarts%D(3,NC)
              Pair%AB2  = (Pair%A(1)-Pair%B(1))**2 &
                        + (Pair%A(2)-Pair%B(2))**2 &
                        + (Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair)) THEN
                 TK%MTrix%D(R:R+NN-1)=TK%MTrix%D(R:R+NN-1)+TBlok(BS,MD,Pair)
                 NotMakeBlock = .FALSE.
              ENDIF
           ENDDO
#else
           TK%MTrix%D(R:R+NN-1)=TBlok(BS,MD,Pair)
#endif
           TK%ColPt%I(P)=AtB
           TK%BlkPt%I(P)=R
           R=R+NN 
           P=P+1 
#ifdef PARALLEL           
           TK%RowPt%I(TK%NAtms+1)=P
           IF(R>MaxNon0Node.OR.P>MaxBlksNode)THEN
              WRITE(*,*)' MyId = ',MyId,MaxBlksNode,MaxNon0Node
              WRITE(*,*)' MyId = ',MyId,' P = ',P,' R = ',R
              CALL Halt(' DBCSR dimensions blown in MakeT ')
           ENDIF
#else  
           TK%RowPt%I(AtA+1)=P        
           IF(R>MaxNon0.OR.P>MaxBlks) &
                CALL Halt(' BCSR dimensions blown in MakeT ')
#endif
#ifdef PERIODIC
           IF(NotMakeBlock) &
                CALL Halt(' Making a Zero Block in MakeT ')
#endif
        ENDIF
     ENDDO
  ENDDO
  TK%NBlks=P-1
  TK%NNon0=R-1
!------------------------------------------------------------
! Put T to disk
!  
  CALL Filter(T2,TK)
  CALL Put(T2,TrixFile('T',Args))
!
! BlksName_O='nti',Non0Name_O='ntm')
!-----------------------------------------------------------
! Printing
!
  CALL PChkSum(T2,'TK',Prog)
  CALL PPrint( T2,'TK')
  CALL Plot(   T2,'TK')
!------------------------------------------------------------
! Tidy up
! 
  CALL Delete(TK)
  CALL Delete(T2)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(MD)
  CALL ShutDown(Prog)
!
END PROGRAM MakeT
