!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 1999, The University of California
!
!    COMPUTES A GUESS NON-ORTHOGONAL DENSITY MATRIX EITHER FROM
!    A SUPERPOSITION OF DIAGONAL, ATOMIC LEWIS STRUCTURE BLOCKS (Action=Direct)
!    OR FROM A ORTHOGONAL DENSITY MATRIX COORESPONDING TO A DIFFERENT GEOMETRY (Action=Switch) 
!
PROGRAM P2Use
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE PBlokGuess
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
  TYPE(ARGMT)                   :: Args
  TYPE(BSET)                    :: BS
  TYPE(CRDS)                    :: GM
#ifdef PARALLEL
  TYPE(DBCSR) &
#else
  TYPE(BCSR)  & 
#endif
                                :: P,T1,T2
  TYPE(DBL_RNK2)                :: BlkP
  REAL(DOUBLE)                  :: TrP,Fact
  INTEGER                       :: I,J,AtA,Q,R,KA,NBFA
  CHARACTER(LEN=2)              :: Cycl
  CHARACTER(LEN=5),PARAMETER    :: Prog='P2Use'
  LOGICAL                       :: Present
!--------------------------------------- 
! Start up macro
!
  CALL StartUp(Args,Prog)
  Cycl=IntToChar(Args%I%I(1))
!--------------------------------------- 
! Get basis set and geometry
!
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
!-------------------------------------- 
! Allocations 
!
  CALL New(P)
  CALL New(BlkP,(/MaxBlkSize**2,NAtoms/))
!
  IF(SCFActn=='Direct')THEN

!    WRITE(*,*)' In P2Use, using SuperPos Guess'
!   Compute a diagonal guess as the superposition of 
!   atomic lewis structure occupancies 
    DO I=1,NAtoms
       CALL FillPBlok(BSiz%I(I),GM%AtNum%I(I),BlkP%D(:,I))
    ENDDO
    CALL SetToI(P,BlkP)
!   Check for the correct elctron count
    TrP=Trace(P)
    IF(ABS(TrP-DBLE(NEl/Two))>1.D-10) &
       CALL Halt(' In P2Use, TrP = '//TRIM(DblToChar(TrP)))
!   IO for the orthogonal P 
    CALL Put(P,TrixFile('OrthoD',Args,0)) 
    CALL PChkSum(P,'OrthoP['//TRIM(Cycl)//']',Prog)
    CALL PPrint( P,'OrthoP['//TRIM(Cycl)//']')
    CALL Plot(   P,'OrthoP_'//TRIM(Cycl))
  ELSEIF(SCFActn=='Switch')THEN
!    WRITE(*,*)' In P2Use, using previous density matrix from file: '
!    WRITE(*,*)' =  ',TrixFile('OrthoD',Args,-1)
!    CALL Get(P,TrixFile('OrthoD',Args,-1))

    WRITE(*,*)' In P2Use, using previous density matrix from info file '
    CALL Get(P,'CurrentOrthoD',CheckPoint_O=.TRUE.)   
    CALL Put(P,TrixFile('OrthoD',Args,0)) 
  ELSEIF(SCFActn=='Restart')THEN
!    WRITE(*,*)' In P2Use, using previous density matrix from info file '
    CALL Get(P,'CurrentOrthoD',CheckPoint_O=.TRUE.)   
    CALL Put(P,TrixFile('OrthoD',Args,0)) 
  ENDIF    
!
  CALL Delete(BlkP)
!---------------------------------------------
!  Orthog to non-orthogonal xformation
!
  CALL New(T1)
  CALL New(T2)
!
   INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
   IF(Present)THEN     
      CALL Get(T1,TrixFile('X',Args))   ! T1=S^(-1/2)
      CALL Multiply(T1,P,T2)            ! T2=S^(-1/2).DiagonalGuess
      CALL Multiply(T2,T1,P)            ! P=S^(-1/2).DiagonalGuess.S^(-1/2)
      CALL Filter(T1,P)                 ! T1=Filter[S^(-1/2).DiagonalGuess.S^(-1/2)]
   ELSE
      CALL Get(T1,TrixFile('Z',Args))   ! T1=Z
      CALL Multiply(T1,P,T2)            ! T2=Z.DiagonalGuess
      CALL Get(T1,TrixFile('ZT',Args))  ! T1=Z^T
      CALL Multiply(T2,T1,P)            ! P=Z.DiagonalGuess.Z^T
      CALL Filter(T1,P)                 ! T1=Filter[Z.DiagonalGuess.Z^T]
   ENDIF     
!--------------------------------------------------------------------
!  IO for the non-orthogonal P 
!
   CALL Put(T1,TrixFile('D',Args,0),        &
            BlksName_O='ndi'//TRIM(Cycl),   &
            Non0Name_O='ndm'//TRIM(Cycl)) 
!----------------------------------------------
   CALL PChkSum(T1,'P['//TRIM(Cycl)//']',Prog)
   CALL PPrint( T1,'P['//TRIM(Cycl)//']')
   CALL Plot(   T1,'P_'//TRIM(Cycl))
!---------------------------------
!  Tidy up ...
!
   CALL Delete(GM)
   CALL Delete(BS)
   CALL Delete(P)
   CALL Delete(T1)
   CALL Delete(T2)
   CALL ShutDown(Prog)   
!------------------------
END PROGRAM 

      
