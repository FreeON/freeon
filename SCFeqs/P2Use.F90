
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
  TYPE(BSET)                    :: BS,OldBS
  TYPE(CRDS)                    :: GM,OldGM
#ifdef PARALLEL
  TYPE(DBCSR) &
#else
  TYPE(BCSR)  & 
#endif
                                :: S,P,X,T0,T1,T2,S1,P0,dP,dS
  TYPE(INT_VECT)                :: Stat
  TYPE(DBL_RNK2)                :: BlkP
  REAL(DOUBLE)                  :: Scale,TrP,Fact,ECount, &
       DeltaP,OldDeltaP,DensityDev,dN,MaxGDIff,GDIff,OldN,M
  INTEGER                       :: I,J,JP,AtA,Q,R,T,KA,NBFA, &
       NPur,PcntPNon0,OldFileID,ICart
  CHARACTER(LEN=2)              :: Cycl
  LOGICAL                       :: Present,SameBasis,SameGeom
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
  IF(SCFActn=='Extrapolate'.OR.SCFActn=='Restart')THEN
     SameBasis=.TRUE.
     SameBasis=.TRUE.
     IF(SCFActn=='Restart')THEN
        ! Old basis and geometry are the same?
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
        CALL Get(OldBS,CurBase)
        CALL Get(OldGM,CurGeom)
        IF(OldBS%BName==BS%BName.AND.OldBS%NBasF==BS%NBasF)THEN
           SameBasis=.TRUE.
        ELSE
           SameBasis=.FALSE.
        ENDIF
        MaxGDiff=Zero
        DO ICart=1,GM%Natms
           GDiff=ABS(GM%Carts%D(1,ICart)-OldGM%Carts%D(1,ICart)) + &
                 ABS(GM%Carts%D(2,ICart)-OldGM%Carts%D(2,ICart)) + &
                 ABS(GM%Carts%D(3,ICart)-OldGM%Carts%D(3,ICart)) 
           MaxGDiff=MAX(MaxGDiff,GDiff)
        ENDDO
        IF(MaxGDiff<1D-6)THEN
           SameGeom=.TRUE.
        ELSE
           SameGeom=.FALSE.
        ENDIF
        IF(.NOT.SameGeom.AND..NOT.SameBasis)THEN ! Don't know how to do this
           CALL Halt(' Restart with different basis and geometry in P2Use not supported by P2Use ')
        ELSEIF(.NOT.SameGeom.AND.SameBasis)THEN  ! Restart with new geometry
           CALL Halt(' Restart with same basis different geometry not yet enabled in P2Use ')
           ! Get the old AO-DM
           CALL Get(P0,'CurrentDM',CheckPoint_O=.TRUE.)
           ! Get the old overlap matrix
           CALL Get(T0,'CurrentS',CheckPoint_O=.TRUE.)
        ELSEIF(.NOT.SameBasis.AND.SameGeom)THEN  ! Restart with basis set switch...
           ! Overwrite the new with the old
           CALL Delete(BS)
           CALL Delete(GM)
           CALL Get(BS,CurBase)
           CALL Get(GM,CurGeom)
           ! Compute a new sparse matrix blocking scheme for the old BS
           CALL BlockBuild(GM,BS,BSiz,OffS)
#ifdef PARALLEL
           CALL BCast(BSiz)
           CALL BCast(OffS)
#endif           ! Get the old AO-DM in the old basis 
           CALL Get(P,'CurrentDM',CheckPoint_O=.TRUE.)
           !CALL PChkSum(P,'CurrentDM',Prog)
           !CALL PPrint(P,'CURRENTDM',Unit_O=6)
        ELSE                                     ! Simple restart 
           ! Get the old AO-DM
           CALL Get(P,'CurrentDM',CheckPoint_O=.TRUE.)
           ! IO for the non-orthogonal P 
           CALL Put(P,TrixFile('D',Args,0))
           CALL PChkSum(P,'P['//TRIM(Cycl)//']',Prog)
           CALL PPrint( P,'P['//TRIM(Cycl)//']')
           !  CALL PPrint( P,'P['//TRIM(Cycl)//']',Unit_O=6)
           CALL Plot(   P,'P_'//TRIM(Cycl))
           CALL CloseHDFGroup(HDF_CurrentID)
           CALL CloseHDF(OldFileID)
           ! Reopen current group and HDF
           HDFFileID=OpenHDF(H5File)
           H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
           HDF_CurrentID=H5GroupID
           ! go down 
           CALL ShutDown(Prog)   
        ENDIF
        CALL CloseHDFGroup(HDF_CurrentID)
        CALL CloseHDF(OldFileID)
        ! Reopen current group and HDF
        HDFFileID=OpenHDF(H5File)
        H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
        HDF_CurrentID=H5GroupID
     ENDIF
     ! If extrapolating or restarting with different geometry, use DMPT to get a starting DM
     IF(SCFActn=='Extrapolate'.OR.(SCFActn=='Restart'.AND.SameBasis))THEN
        CALL Delete(P)
        CALL New(dS)
        CALL New(dP)
        CALL SetToI(dP)
        CALL Multiply(dP,Zero)
        IF(SCFActn=='Extrapolate')THEN
           ! Get previous non-orthogonal density matrix 
           CALL Get(P0,TrixFile('D',Args,-1))     
           CALL Get(T0,TrixFile('S',Args,Stats_O=Previous))
        ENDIF
        CALL Get(S1,TrixFile('S',Args,Stats_O=Current))
        CALL Multiply(T0,-One)
        ! dS=S1-S0
        CALL Add(S1,T0,dS)        
        OldN=BIG_DBL
        DO I=1,30
           ! T0=P0*[dS*P0+S1*dP]
           CALL Multiply(dS,P0,T0)
           CALL Multiply(S1,dP,T1)
           CALL Add(T0,T1,T2)
           CALL Multiply(P0,T2,T0)
           ! T1=dP*S1*[P0+dP]
           CALL Add(P0,dP,T1)
           CALL Multiply(S1,T1,T2)
           CALL Multiply(dP,T2,T1)
           ! T2=T0+T1=P0*[dS*P0+S1*dP]+dP*S1*[P0+dP]
           CALL Add(T1,T0,T2)
           ! P1~P0+dP; dN~Tr(P1.S1)
           CALL Add(P0,dP,T1)
#ifdef PARALLEL
           CALL Multiply(T1,S1,T0)
           dN=Two*Trace(T0)-NEl
           IF(MyId==ROOT)THEN
#else
           dN=Two*Trace(T1,S1)-NEl
#endif
!              IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
                 Mssg=ProcessName(Prog,'AO-DMX '//TRIM(IntToChar(I)))//'dN='//TRIM(DblToShrtChar(dN))
                 CALL OpenASCII(OutFile,Out)
                 CALL PrintProtectL(Out)
                 WRITE(*,*)TRIM(Mssg)
                 WRITE(Out,*)TRIM(Mssg)
                 CALL PrintProtectR(Out)
                 CLOSE(UNIT=Out,STATUS='KEEP')
!              ENDIF
#ifdef PARALLEL
           ENDIF
#endif
           IF(I>6.AND.ABS(dN)>ABS(OldN))EXIT
           OldN=dN
           IF(MOD(I,2)==0.AND.I<6)THEN
              CALL Multiply(dP,Two)
              CALL Multiply(T2,-One)
              CALL Add(dP,T2,T1)
              CALL Filter(dP,T1)
           ELSEIF(I<6)THEN
              CALL Filter(dP,T2)
           ELSEIF(dN<0)THEN
              CALL Multiply(dP,Two)
              CALL Multiply(T2,-One)
              CALL Add(dP,T2,T1)
              CALL Filter(dP,T1)
           ELSE
              CALL Filter(dP,T2)
           ENDIF
        ENDDO
        CALL Add(dp,P0,T0)
        CALL Put(T0,TrixFile('D',Args,0))     
        CALL Put(T0,'CurrentDM',CheckPoint_O=.TRUE.)
        CALL ShutDown(Prog)
     ENDIF
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
