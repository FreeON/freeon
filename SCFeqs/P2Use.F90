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
                                :: S,P,X,T0,T1,T2
  TYPE(INT_VECT)                :: Stat
  TYPE(DBL_RNK2)                :: BlkP
  REAL(DOUBLE)                  :: Scale,TrP,Fact,ECount, &
       DeltaP,OldDeltaP,DensityDev
  INTEGER                       :: I,J,JP,AtA,Q,R,T,KA,NBFA, &
       NPur,PcntPNon0,OldFileID
  CHARACTER(LEN=2)              :: Cycl
  LOGICAL                       :: Present
  CHARACTER(LEN=DEFAULT_CHR_LEN):: Mssg,BName
  CHARACTER(LEN=5),PARAMETER    :: Prog='P2Use'
!------------------------------------------------------------------------------- 
  ! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
  Cycl=IntToChar(Args%I%I(1))
  ! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)

  ! Allocations 
  CALL New(P)
  CALL New(T0)
  CALL New(T1)
  CALL New(T2)

  IF(SCFActn=='Extrapolate')THEN
     CALL Halt(' Extrapolation turned off, need non-orthogonal SP2 or TS4... ')
  ELSEIF(SCFActn=='Restart')THEN  
     ! Close current group and HDF
     CALL CloseHDFGroup(H5GroupID)
     CALL CloseHDF(HDFFileID)
     ! Open old group and HDF
     OldFileID=OpenHDF(Restart)
     HDF_CurrentID=OpenHDF(Restart)
     ! Get old basis set stuff
     CALL New(Stat,3)
     CALL Get(Stat,'current_state')
     SCFCycl=TRIM(IntToChar(Stat%I(1)))
     CurBase=TRIM(IntToChar(Stat%I(2)))
     CurGeom=TRIM(IntToChar(Stat%I(3)))
     ! Open the old group
     HDF_CurrentID=OpenHDFGroup(OldFileID,"Clone #"//TRIM(IntToChar(MyClone)))
     ! Get the old basis set and geometry for indexing the DM
     CALL Get(BS,CurBase)
     CALL Get(GM,CurGeom)
     ! Compute a new sparse matrix blocking scheme for the old BS
     CALL BlockBuild(GM,BS,BSiz,OffS)
#ifdef PARALLEL
     CALL BCast(BSiz)
     CALL BCast(OffS)
#endif
     ! Find the current orthogonal density matrix
     CALL Get(P,'CurrentDM',CheckPoint_O=.TRUE.)
!     CALL PPrint(P,'CURRENTDM',Unit_O=6)
     ! Close it up 
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(OldFileID)
     ! Reopen current group and HDF
     HDFFileID=OpenHDF(H5File)
     H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
     HDF_CurrentID=H5GroupID
  ELSEIF(SCFActn=='Project')THEN
     ! Get previous geometries orthogonal density matrix 
     CALL Get(P,TrixFile('OrthoD',Args,-1))     
#ifdef PARALLEL
     IF(MyId==ROOT)THEN
#endif
        INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
#ifdef PARALLEL
     ENDIF
     CALL BCast(Present)
#endif
  ELSEIF(SCFActn=='DensitySuperposition')THEN
     ! Make a diagonal guess
     CALL Get(BName,'bsetname',CurBase)
     IF(INDEX(BName,'STO')/=0)THEN
        ! Compute a diagonal guess as the superposition of atomic lewis 
        ! structure occupancies--works only for minimal (STO) basis sets
        CALL New(BlkP,(/MaxBlkSize**2,NAtoms/))
        DO I=1,NAtoms
           IF(GM%AtNum%D(I) < 105.D0) THEN
              CALL FillPBlok(BSiz%I(I),INT(GM%AtNum%D(I)),BlkP%D(:,I))
           ENDIF
        ENDDO
        CALL SetToI(P,BlkP)
        ! Check for the correct elctron count
        TrP=Trace(P)
        IF(ABS(TrP-DBLE(NEl/Two))>1.D-10) &
             CALL Warn(' In P2Use, TrP = '//TRIM(DblToChar(TrP)))
        CALL Delete(BlkP)
     ELSE
        CALL SetToI(P)
        CALL Warn('Attempting to use density superpostion with a non STO basis set. Going for scaled I.')
        CALL Multiply(P,DBLE(NEl)/DBLE(2*NBasF))
        TrP=Trace(P)
        IF(ABS(TrP-DBLE(NEl/Two))>1.D-10) &
             CALL Warn(' In P2Use, TrP = '//TRIM(DblToChar(TrP)))
     ENDIF
  ELSEIF(SCFActn=='GuessEqCore')THEN
     ! Guess == Core
     CALL New(BlkP,(/MaxBlkSize**2,NAtoms/))
     DO I=1,NAtoms
        BlkP%D(:,I)=Zero
     ENDDO
     CALL SetToI(P,BlkP)
     CALL Delete(BlkP)
  ELSE
     CALL Halt(' Unknown option '//TRIM(SCFActn))
  ENDIF
  IF(SCFActn/='GuessEqCore')THEN
#ifdef PARALLEL
        IF(MyId==ROOT)THEN
#endif
           INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
#ifdef PARALLEL
        ENDIF
        CALL BCast(Present)
#endif
     IF(Present)THEN     
        CALL Get(T1,TrixFile('X',Args))   ! T1=S_new^(-1/2)
        CALL Multiply(T1,P,T2)            ! T2=S_new^(-1/2).P_old
        CALL Multiply(T2,T1,T0)           ! P_new_AO=S_new^(-1/2).P_old.S_new^(-1/2)
        CALL Filter(P,T0)                 ! T1=Filter[P_new_AO]
     ELSE
        CALL Get(T1,TrixFile('Z',Args))   ! T1=Z_new
        CALL Multiply(T1,P,T2)            ! T2=Z.P_old
        CALL Get(T1,TrixFile('ZT',Args))  ! T1=Z^T
        CALL Multiply(T2,T1,T0)           ! P_new_AO=Z.P_old.Z^T
        CALL Filter(P,T0)                 ! T1=Filter[P_new_AO]
     ENDIF
     CALL Delete(T0)
     CALL Delete(T1)
     CALL Delete(T2)
  ENDIF
  ! IO for the non-orthogonal P 
  CALL Put(P,TrixFile('D',Args,0))
  CALL PChkSum(P,'P['//TRIM(Cycl)//']',Prog)
  CALL PPrint( P,'P['//TRIM(Cycl)//']')
!  CALL PPrint( P,'P['//TRIM(Cycl)//']',Unit_O=6)
  CALL Plot(   P,'P_'//TRIM(Cycl))
  ! Tidy up ...
  CALL Delete(GM)
  CALL Delete(BS)
  CALL Delete(P)
  CALL ShutDown(Prog)   
END PROGRAM P2Use

      
