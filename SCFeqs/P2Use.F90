!    COMPUTES AN ORTHOGONAL GUESS DENSITY MATRIX EITHER FROM
!    A SUPERPOSITION OF DIAGONAL, ATOMIC LEWIS STRUCTURE BLOCKS 
!    OR FROM A ORTHOGONAL DENSITY MATRIX COORESPONDING TO A DIFFERENT 
!    (HOPEFULLY CLOSE) GEOMETRY OR LEVEL OF ACCURACY 
!    Author: Matt Challacombe
!-------------------------------------------------------------------------
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
                                :: P,T0,T1,T2
  TYPE(DBL_RNK2)                :: BlkP
  REAL(DOUBLE)                  :: TrP,Fact,ECount,DensityDev
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
  IF(SCFActn=='Extrapolate')THEN
     CALL New(T0)
     CALL New(T1)
     CALL New(T2)
     CALL Get(T0,TrixFile('S',Args,Stats_O=Previous))
     CALL Get(T1,TrixFile('S',Args,Stats_O=Current))
     CALL Multiply(T1,-One)
!    dS=Sp-Sc 
     CALL Add(T1,T0,T2)        
!    Get previous density matrix 
     CALL Get(P,TrixFile('D',Args,-1))     
!    P.dS
     CALL Multiply(P,T2,T1)
!    dP=-P.dS.P
     CALL Multiply(T1,P,T0)
!    P'=P+dP Can probably come up with an extrapolation to higher order...
!     CALL Multiply(T0,Zero)
     CALL Add(T0,P,T1)
!    Check for normalization in corrected DP
     CALL Filter(P,T1)
     CALL Get(T0,TrixFile('S',Args,Stats_O=Current))
!    The following is clumsy.  Should have better diagonal only renormalization
!    followed by a purify....
     ECount=Trace(T0,P)
     DensityDev=DBLE(Nel)-Two*ECount
!    Renormalize corrected density matrix 
     CALL Multiply(P,DBLE(NEl)/(Two*ECount))
  ELSE
     IF(SCFActn=='Restart')THEN
       CALL Get(P,'CurrentOrthoD',CheckPoint_O=.TRUE.)   
       CALL Put(P,TrixFile('OrthoD',Args,0)) 
     ELSEIF(SCFActn=='Project')THEN
!      Get previous geometries orthogonal density matrix 
       CALL Get(P,TrixFile('OrthoD',Args,-1))     
     ELSE
!      Compute a diagonal guess as the superposition of atomic lewis 
!      structure occupancies--works only for minimal (STO) basis sets
       DO I=1,NAtoms
          CALL FillPBlok(BSiz%I(I),GM%AtNum%I(I),BlkP%D(:,I))
       ENDDO
       CALL SetToI(P,BlkP)
!      Check for the correct elctron count
       TrP=Trace(P)
       IF(ABS(TrP-DBLE(NEl/Two))>1.D-10) &
          CALL Halt(' In P2Use, TrP = '//TRIM(DblToChar(TrP)))
       CALL Delete(BlkP)
     ENDIF
     CALL New(T0)
     CALL New(T1)
     CALL New(T2)
     INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
     IF(Present)THEN     
         CALL Get(T1,TrixFile('X',Args))   ! T1=S^(-1/2)
         CALL Multiply(T1,P,T2)            ! T2=S^(-1/2).DiagonalGuess
         CALL Multiply(T2,T1,T0)           ! P=S^(-1/2).DiagonalGuess.S^(-1/2)
         CALL Filter(P,T0)                 ! T1=Filter[S^(-1/2).DiagonalGuess.S^(-1/2)]
      ELSE
         CALL Get(T1,TrixFile('Z',Args))   ! T1=Z
         CALL Multiply(T1,P,T2)            ! T2=Z.DiagonalGuess
         CALL Get(T1,TrixFile('ZT',Args))  ! T1=Z^T
         CALL Multiply(T2,T1,T0)           ! P=Z.DiagonalGuess.Z^T
         CALL Filter(P,T0)                 ! T1=Filter[Z.DiagonalGuess.Z^T]
      ENDIF     
  ENDIF    
! IO for the non-orthogonal P 
  CALL Put(P,TrixFile('D',Args,0))
!----------------------------------------------
   CALL PChkSum(P,'P['//TRIM(Cycl)//']',Prog)
   CALL PPrint( P,'P['//TRIM(Cycl)//']')
   CALL Plot(   P,'P_'//TRIM(Cycl))
!---------------------------------
!  Tidy up ...
!
   CALL Delete(GM)
   CALL Delete(BS)
   CALL Delete(P)
   CALL Delete(T0)
   CALL Delete(T1)
   CALL Delete(T2)
   CALL ShutDown(Prog)   
!------------------------
END PROGRAM 

      
