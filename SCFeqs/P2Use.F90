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
                                :: S,P,T0,T1,T2
  TYPE(DBL_RNK2)                :: BlkP
  REAL(DOUBLE)                  :: Scale,TrP,Fact,ECount, &
                                   MxDelP,DensityDev
  INTEGER                       :: I,J,AtA,Q,R,KA,NBFA, &
                                   NPur,PcntPNon0
  CHARACTER(LEN=2)              :: Cycl
  LOGICAL                       :: Present
  CHARACTER(LEN=DEFAULT_CHR_LEN):: Mssg,BName,RestartHDF
  CHARACTER(LEN=5),PARAMETER    :: Prog='P2Use'
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
  CALL New(T0)
  CALL New(T1)
  CALL New(T2)
!
  IF(SCFActn=='Extrapolate')THEN
     CALL Get(T0,TrixFile('S',Args,Stats_O=Previous))
     CALL Get(S,TrixFile('S',Args,Stats_O=Current))
     CALL Multiply(S,-One)
!    dS=Sp-Sc 
     CALL Add(S,T0,T2)        
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
     CALL Multiply(S,-One)
     NPur=0
     DO I=1,20
        NPur=NPur+1
!       P.S
        CALL Multiply(P,S,T1)
!       P.S.P.S
        CALL Multiply(T1,T1,T0)
!       P.S.P
        CALL Multiply(T1,P,T2)
!       P.S.P.S.P
        CALL Multiply(T0,P,T1)
!       3*P.S.P
        CALL Multiply(T2,Three)
!       -2.P.S.P.S.P
        CALL Multiply(T1,-Two)
!       P[i+1]=3*Pi.S.Pi-2.Pi.S.Pi.S.Pi
        CALL Add(T2,T1,T0)
!       DeltaP
        CALL Multiply(P,-One)
        CALL Add(P,T0,T1)
        MxDelP=MAX(T1)
        CALL Filter(P,T0)
        ECount=Trace(S,P)
        DensityDev=DBLE(Nel)-Two*ECount
#ifdef PARALLEL
        IF(MyId==ROOT)THEN
#endif
          IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
              CALL OpenASCII(OutFile,Out)
              CALL PrintProtectL(Out)
              PcntPNon0=INT(1.D2*DBLE(P%NNon0)/DBLE(NBasF*NBasF))
              Mssg=ProcessName(Prog,'Pure '//TRIM(IntToChar(NPur)))            &
                             //'Tr(P)-NEl= '//TRIM(DblToShrtChar(DensityDev))  &
                             //', %Non0= '//TRIM(IntToChar(PcntPNon0))         &
                             //', MAX(/P) = '//TRIM(DblToShrtChar(MxDelP)) 
            WRITE(*,*)TRIM(Mssg)
            WRITE(Out,*)TRIM(Mssg)
            CALL PrintProtectR(Out)
            CLOSE(UNIT=Out,STATUS='KEEP')
         ENDIF
#ifdef PARALLEL
      ENDIF
#endif
        IF(DensityDev<1.D-5)GOTO 999
     ENDDO
     CALL Warn('In P2Use, failed to converge McWeeny purification.'//RTRN    &
              //'   Still missing '//TRIM(DblToShrtChar(DensityDev))//' electrons.'//RTRN &
              //'   Will try projection instead of interpolation.')
     SCFActn='Project'      
     GOTO 1
999 CONTINUE
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
       IF(PrintFlags%Key>DEBUG_MINIMUM)THEN
          CALL OpenASCII(OutFile,Out)
          CALL PrintProtectL(Out)
          Mssg=ProcessName(Prog)//TRIM(IntToChar(NPur))//' purification steps taken.'
          IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
             WRITE(*,*)TRIM(Mssg)
          ENDIF
          WRITE(Out,*)TRIM(Mssg)
          Mssg=ProcessName(Prog)//'MAX(/\P) = '//TRIM(DblToShrtChar(MxDelP))//', ' &
                //'|Tr(P)-NEl| = '//TRIM(DblToShrtChar(DensityDev))//' .'
          IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
             WRITE(*,*)TRIM(Mssg)
          ENDIF
          WRITE(Out,*)TRIM(Mssg)
          CALL PrintProtectR(Out)
          CLOSE(UNIT=Out,STATUS='KEEP')
       ENDIF
#ifdef PARALLEL
   ENDIF
#endif
  ENDIF
1 CONTINUE
!
  IF(SCFActn=='Restart')THEN
     CALL Get(RestartHDF,'OldInfo')
     CALL CloseHDF()
     CALL OpenHDF(RestartHDF)
     CALL Get(P,'CurrentDM',CheckPoint_O=.TRUE.)   
     CALL CloseHDF()
     CALL OpenHDF(InfFile)
  ELSEIF(SCFActn=='Project')THEN
!    Get previous geometries orthogonal density matrix 
     CALL Get(P,TrixFile('OrthoD',Args,-1))     
  ELSEIF(SCFActn=='DensitySuperposition')THEN
!    Make a diagonal guess
     CALL Get(BName,'bsetname',CurBase)
     IF(INDEX(BName,'STO')/=0)THEN
!         Compute a diagonal guess as the superposition of atomic lewis 
!         structure occupancies--works only for minimal (STO) basis sets
          CALL New(BlkP,(/MaxBlkSize**2,NAtoms/))
          DO I=1,NAtoms
             CALL FillPBlok(BSiz%I(I),GM%AtNum%I(I),BlkP%D(:,I))
          ENDDO
          CALL SetToI(P,BlkP)
!         Check for the correct elctron count
          TrP=Trace(P)
          IF(ABS(TrP-DBLE(NEl/Two))>1.D-10) &
             CALL Halt(' In P2Use, TrP = '//TRIM(DblToChar(TrP)))
          CALL Delete(BlkP)
        ELSE
!         Evenly wheigted diagonal guess (better than core ;)
          CALL SetToI(P)
          Scale=Half*DBLE(NEl)/DBLE(NBasF)
          CALL Multiply(P,Scale)
        ENDIF
  ENDIF
! Create density matrix in AO representation
  IF(SCFActn/='Extrapolate'.AND.SCFActn/='Restart')THEN
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
  CALL PChkSum(P,'P['//TRIM(Cycl)//']',Prog)
  CALL PPrint( P,'P['//TRIM(Cycl)//']')
  CALL Plot(   P,'P_'//TRIM(Cycl))
! Tidy up ...
  CALL Delete(GM)
  CALL Delete(BS)
  CALL Delete(P)
  CALL Delete(T0)
  CALL Delete(T1)
  CALL Delete(T2)
  CALL ShutDown(Prog)   
END PROGRAM 

      
