
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
  USE DenMatMethods
  IMPLICIT NONE
  TYPE(ARGMT)                   :: Args
  TYPE(BSET)                    :: BS,BS_old
  TYPE(CRDS)                    :: GM,GM_old
#ifdef PARALLEL
  TYPE(DBCSR) &
#else
  TYPE(BCSR)  & 
#endif
                                :: P,P0,X,S,S0,S1,dS,Tmp1,Tmp2
  TYPE(INT_VECT)                :: Stat
  TYPE(DBL_RNK2)                :: BlkP
  REAL(DOUBLE)                  :: MaxDS
  REAL(DOUBLE)                  :: Scale,Fact,ECount,RelNErr, DeltaP,OldDeltaP, & 
                                   DensityDev,dN,MaxGDIff,GDIff,OldN,M,PNon0s,PSMin,PSMax, &
                                   Ipot_Error,Norm_Error,Lam,TError0
  INTEGER                       :: I,J,JP,AtA,Q,R,T,KA,NBFA,NPur,PcntPNon0,Qstep, & 
                                   OldFileID,ICart,N,NStep
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
  ! Do what needs to be done CASE by CASE
  SELECT CASE(SCFActn)
  ! P=0
  CASE('GuessEqCore')
     CALL New(P)
     CALL New(BlkP,(/MaxBlkSize**2,NAtoms/))
     DO I=1,NAtoms
        BlkP%D(:,I)=Zero
     ENDDO
     CALL SetToI(P,BlkP)
     CALL Delete(BlkP)
     ! IO for the non-orthogonal P 
     CALL Put(P,TrixFile('D',Args,0))
     CALL PChkSum(P,'P['//TRIM(Cycl)//']',Prog)
     CALL PPrint( P,'P['//TRIM(Cycl)//']')
     CALL Plot(   P,'P_'//TRIM(Cycl))
     CALL Delete(P)
  ! Density SuperPosition 
  CASE('DensitySuperposition')
     CALL New(P)
     CALL New(Tmp1)
     CALL New(Tmp2)
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
#ifdef PARALLEL
  IF(MyID == ROOT) THEN
#endif
        CALL Warn('Attempting to use density superpostion with a non STO basis set. Going for scaled I.')
#ifdef PARALLEL
  ENDIF
#endif
        CALL Multiply(P,DBLE(NEl)/DBLE(2*NBasF))
        TrP=Trace(P)
        IF(ABS(TrP-DBLE(NEl/Two))>1.D-10) CALL Warn(' In P2Use, TrP = '//TRIM(DblToChar(TrP)))
     ENDIF
#ifdef PARALLEL    
     IF(MyId==ROOT)THEN
#endif
        INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
#ifdef PARALLEL
     ENDIF
     CALL BCast(Present)
#endif
     INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
     IF(Present)THEN     
        CALL Get(X,TrixFile('X',Args))    ! T1=S_new^(-1/2)
        CALL Multiply(X,P,Tmp1)           ! T2=S_new^(-1/2).P_old
        CALL Multiply(Tmp1,X,Tmp2)        ! P_new_AO=S_new^(-1/2).P_old.S_new^(-1/2)
        CALL Filter(P,Tmp2)               ! T1=Filter[P_new_AO]
     ELSE
        CALL Get(X,TrixFile('Z',Args))    ! T1=Z_new
        CALL Multiply(X,P,Tmp1)           ! T2=Z.P_old
        CALL Get(X,TrixFile('ZT',Args))   ! T1=Z^T
        CALL Multiply(Tmp1,X,Tmp2)        ! P_new_AO=Z.P_old.Z^T
        CALL Filter(P,Tmp2)               ! T1=Filter[P_new_AO]
     ENDIF
     ! IO for the non-orthogonal P 
     CALL Put(P,TrixFile('D',Args,0))
     CALL PChkSum(P,'P['//TRIM(Cycl)//']',Prog)
     CALL PPrint( P,'P['//TRIM(Cycl)//']')
     CALL Plot(   P,'P_'//TRIM(Cycl))
     CALL Delete(P)
     CALL Delete(X)
     CALL Delete(Tmp1)
     CALL Delete(Tmp2)
  ! Restarting without Geometry of BasisSet Change
  CASE('Restart')
     ! Close Current Group
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
     ! Get the old AO-DM
     CALL Get(P,'CurrentDM',CheckPoint_O=.TRUE.)
     ! IO for the non-orthogonal P 
     CALL Put(P,TrixFile('D',Args,0))
     CALL PChkSum(P,'P['//TRIM(Cycl)//']',Prog)
     CALL PPrint( P,'P['//TRIM(Cycl)//']')
     CALL Plot(   P,'P_'//TRIM(Cycl))
     ! Close Old group
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(OldFileID)
     ! Reopen current group and HDF
     HDFFileID=OpenHDF(H5File)
     H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
     HDF_CurrentID=H5GroupID
     CALL Delete(P)
  ! Restarting with BasisSet Change
  CASE('RestartBasisSwitch')
     ! Close Current Group
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
     ! Get the Old Basis Set and Geometry
     CALL Get(BS_old,CurBase)
     CALL Get(GM_old,CurGeom)      
     ! Compute a new sparse matrix blocking scheme for the old BS
     CALL BlockBuild(GM_old,BS_old,BSiz,OffS)
     NBasF=BS_old%NBasF
#ifdef PARALLEL
     CALL BCast(BSiz)
     CALL BCast(OffS)
     CALL BCast(NBasF)
#endif           
     ! Get the old AO-DM
     CALL Get(P,'CurrentDM',CheckPoint_O=.TRUE.)
     ! IO for the non-orthogonal P 
     CALL Put(P,TrixFile('D',Args,0))
     CALL PChkSum(P,'P['//TRIM(Cycl)//']',Prog)
     CALL PPrint( P,'P['//TRIM(Cycl)//']')
     CALL Plot(   P,'P_'//TRIM(Cycl))
     ! Close Old group
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(OldFileID)
     ! Reopen current group and HDF
     HDFFileID=OpenHDF(H5File)
     H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
     HDF_CurrentID=H5GroupID
     CALL Delete(P)
  ! Geometry Change 
  CASE('Extrapolate')
     ! Get Marices
     IF(PrvGeom==CurGeom) THEN
        CALL Get(S0,TrixFile('S',Args,Stats_O=(/Current(1),Current(2),Current(3)-1/)))
        CALL Get(S1,TrixFile('S',Args,Stats_O=Current))
        ! Close Current Group
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
        ! Get the old AO-DM
        CALL Get(P0,'CurrentDM',CheckPoint_O=.TRUE.) 
        ! Close Old group
        CALL CloseHDFGroup(HDF_CurrentID)
        CALL CloseHDF(OldFileID)
        ! Reopen current group and HDF
        HDFFileID=OpenHDF(H5File)
        H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
        HDF_CurrentID=H5GroupID
     ELSE
        CALL Get(S0,TrixFile('S',Args,Stats_O=Previous))
        CALL Get(S1,TrixFile('S',Args,Stats_O=Current))
        CALL Get(P0,TrixFile('D',Args,-1))
     ENDIF
     ! Allocate
     CALL New(P)
     CALL New(S)
     CALL New(dS)
     CALL New(Tmp1)
     CALL New(Tmp2)
     ! Compute dS
     CALL Multiply(S1,-One)
     CALL Add(S1,S0,dS)
     CALL Multiply(S1,-One)
     MaxDS = SQRT(Trace(dS,dS))
     ! Initialize P
     CALL SetEq(P,P0)
     ! Initial Trace Error
     TError0  = ABS(Trace(P,S1)-Half*DBLE(NEl))
     !Number of Steps
     NStep=INT(ABS(Trace(P,S1)-Half*DBLE(NEl)))/4+1
     ! Purify 
     IF(MaxDS > 1.D-10) THEN
        DO N=1,NStep
           ! Initialize S
           Lam  = DBLE(N)/DBLE(NStep)
           CALL SetEq(Tmp1,S0)
           CALL SetEq(Tmp2,S1)
           CALL Multiply(Tmp1,One-Lam)
           CALL Multiply(Tmp2,Lam)
           CALL Add(Tmp1,Tmp2,S)
           TrP        = Trace(P,S)
           Norm_Error = TrP-Half*DBLE(NEl)
           Ipot_Error = One
!
           Qstep=0
           DO I=1,20
              IF((ABS(Ipot_Error) > 0.5D0 .OR. ABS(Norm_Error) > 0.5D0) .AND. I > 1) THEN
                 IF(Norm_Error > Zero) THEN
                    CALL AOSP2(P,S,Tmp1,Tmp2,.TRUE.)
                 ELSE
                    CALL AOSP2(P,S,Tmp1,Tmp2,.FALSE.)
                 ENDIF
              ELSE
                 IF(I > 1) Qstep = Qstep+1
                 IF(Norm_Error > Zero) THEN
                    CALL AOSP2(P,S,Tmp1,Tmp2,.TRUE.)
                    CALL AOSP2(P,S,Tmp1,Tmp2,.FALSE.)
                 ELSE
                    CALL AOSP2(P,S,Tmp1,Tmp2,.FALSE.)
                    CALL AOSP2(P,S,Tmp1,Tmp2,.TRUE.)
                 ENDIF
              ENDIF
              Norm_Error = TrP-Half*DBLE(NEl)
              Ipot_Error = TrP2-Trace(P)
#ifdef PARALLEL
              IF(MyId==ROOT)THEN
#endif
                 PNon0s=100.D0*DBLE(P%NNon0)/DBLE(NBasF*NBasF)
                 Mssg=ProcessName(Prog,'AO-DMX '//TRIM(IntToChar(I))) &
                      //'dN='//TRIM(DblToShrtChar(Norm_Error)) &
                      //', %Non0='//TRIM(DblToShrtChar(PNon0s))                  
                 CALL OpenASCII(OutFile,Out)
                 CALL PrintProtectL(Out)
                 WRITE(*,*)TRIM(Mssg)
                 WRITE(Out,*)TRIM(Mssg)
                 CALL PrintProtectR(Out)
                 CLOSE(UNIT=Out,STATUS='KEEP')
#ifdef PARALLEL          
              ENDIF
#endif
              IF(ABS(Ipot_Error) < 1.0D-10 .AND. ABS(Norm_Error) < 1.0D-10) EXIT
              IF(ABS(Norm_Error) > 100.D0*TError0) EXIT
              IF(Qstep > 4) EXIT
           ENDDO
        ENDDO
!
        IF(ABS(Norm_Error) > TError0) THEN
           CALL Warn("Using Old Density Matrix: Norm Error to Large")
           CALL Filter(P,P0)
        ENDIF
     ELSE
        CALL Warn("Using Old Density Matrix: MaxDS too small")
        CALL Filter(P,P0)
     ENDIF
     ! Save back to be sure.
     CALL Put(P,TrixFile('D',Args,0))
     CALL Put(P,'CurrentDM',CheckPoint_O=.TRUE.) 
     ! Clean Up
     CALL Delete(P)
     CALL Delete(S)
     CALL Delete(dS) 
     CALL Delete(Tmp1)
     CALL Delete(Tmp2)
  CASE('Project')
     CALL Halt(' Not Implimented '//TRIM(SCFActn))
  CASE DEFAULT
     CALL Halt(' Unknown option '//TRIM(SCFActn))
  END SELECT
  ! Tidy up ...
  CALL Delete(GM)
  CALL Delete(BS)
  CALL ShutDown(Prog)   
!
END PROGRAM P2Use
