!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 1999, The University of California
!
!    COMPUTE THE DERIVATIVE OF THE OVERLAP MATRIX S
!
PROGRAM GradS
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
  USE GradSBlock
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef PARALLEL
  TYPE(DBCSR)         :: T1,W,C,F,P
#else
  TYPE(BCSR)          :: T1,W,C,F,P
  TYPE(TIME)          :: M1,M2
#endif
#ifdef PERIODIC 
  INTEGER             :: NCells,NX,NY,NZ,I
  LOGICAL             :: InitBlock
#endif
  TYPE(BSET)          :: BS
  TYPE(CRDS)          :: GM
  TYPE(ARGMT)         :: Args
  TYPE(DBL_RNK4)      :: MD
  TYPE(DBL_RNK2)      :: TrWS
  TYPE(AtomPair)                          :: Pair
  REAL(DOUBLE)        :: Ax,Ay,Az,Bx,By,Bz,    &
                         ABx,ABy,ABz,AB2,Test, &
                         BlockSize

  INTEGER             :: Q,R,AtA,AtB,NBFA,NBFB,NN,KA,KB, &
                         IL,MA,J,JP,K,KP,PP,MB,QQ, &
                         IStrtW,IStopW,IStrtW1,IStopW1, &
                         JP1,J1,PP1,DiffAtm,LP,Top
  CHARACTER(LEN=5),PARAMETER :: Prog='GradS'
!--------------------------------------- 
! Start up macro
!
  CALL StartUp(Args,Prog)
!--------------------------------------- 
! Get basis set and geometry
!
  CALL Get(BS)
  CALL Get(GM)
!-------------------------------------- 
! Allocations 
! 
  CALL New(MD,(/3,BS%NASym+1,BS%NASym+1,2*BS%NASym+2/),(/1,-1,-1,0/))
  CALL New(W)
  CALL New(C)
  CALL New(F)
  CALL New(P)
  CALL New(TrWS,(/3,NAtoms/))
!-----------------------------------------------
! Initialize the matrix and associated indecies
!
  W%RowPt%I(1)=1;C%RowPt%I(1)=1
  CALL SetEq(W%MTrix,Zero)
  CALL SetEq(C%MTrix,Zero)
!--------------------------------------
! First compute the W matrix
!
!  CALL Elapsed_TIME(M1,'Init',Proc_O='MakeW')
  CALL Get(F,TrixFile('F',Args,0))
  CALL Get(P,TrixFile('D',Args,1))
  CALL Multiply(F,P,T1)
  CALL Filter(C,T1)
  CALL Multiply(C,P,T1)
  CALL Filter(W,T1)
  CALL Delete(F)
  CALL Delete(P)
  CALL Delete(C)
!  CALL Elapsed_TIME(M1,'Accum',Proc_O='MakeW')
!  CALL PPrint(M1,'MakeW')
!-------------------------------------- 
! Main loops
!
!  CALL Elapsed_TIME(M2,'Init',Proc_O='Loops')
  TrWS%D=Zero
  DO AtA=1,NAtoms
     MA=BSiz%I(AtA)
     IStrtW=W%RowPt%I(AtA)
     IStopW=W%RowPt%I(AtA+1)-1
     DO JP=IStrtW,IStopW
        AtB=W%ColPt%I(JP)
        PP=W%BlkPt%I(JP)
        MB=BSiz%I(AtB)
        IF(SetAtomPair(GM,BS,AtB,AtA,Pair)) THEN
           TrWS%D(:,AtA)=TrWS%D(:,AtA) + &
                GradSBlok(BS,MD,AtB,AtA,AtA,W%MTrix%D(PP:PP+MA*MB-1),Pair)
           IF (AtA.NE.AtB) THEN
              TrWS%D(:,AtB)=TrWS%D(:,AtB) + &
                   GradSBlok(BS,MD,AtB,AtA,AtB,W%MTrix%D(PP:PP+MA*MB-1),Pair)
           ENDIF
        ENDIF
     ENDDO
  ENDDO
!  CALL Elapsed_TIME(M2,'Accum',Proc_O='Loops')
!  CALL PPrint(M2,'Loops')

!-----------------------------------------------------------
! Printing
!
!  PrintFlags%Fmt=DBG_DBL_STYLE
!  PrintFlags%Key=DBG_MAXIMUM 
  CALL PChkSum(T1,'W',Prog)
!  CALL PPrint( T1,'W')
  CALL Plot(   T1,'W')

!------------------------------------------------------------
! Tidy up
! 
  CALL Delete(T1)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(MD)
  CALL ShutDown(Prog)
!
END PROGRAM GradS
