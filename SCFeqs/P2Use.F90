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
  TYPE(BCSR) :: P_BCSR,F_BCSR
  TYPE(DBCSR) &
#else
  TYPE(BCSR)  & 
#endif
                                :: F,P,P0,X,S,S0,S1,Tmp1,Tmp2

  TYPE(INT_VECT)                :: Stat
  TYPE(DBL_RNK2)                :: BlkP
  REAL(DOUBLE)                  :: MaxDS,NoiseLevel 
  REAL(DOUBLE)                  :: Scale,Fact,ECount,RelNErr, DeltaP,OldDeltaP, & 
                                   DensityDev,dN,MaxGDIff,GDIff,OldN,M,PNon0s,PSMin,PSMax, &
                                   Ipot_Error,Norm_Error,Lam,DLam,TError0,SFac,Dum,Fmin,Fmax
  INTEGER                       :: I,J,JP,AtA,Q,R,T,KA,NBFA,NPur,PcntPNon0,Qstep, & 
                                   OldFileID,ICart,N,NStep,iGEO,DMPOrder,NSMat,MM,ICycle,Cycle
  CHARACTER(LEN=2)              :: Cycl
  LOGICAL                       :: Present,DoingMD,ConvergeAOSP,ConvergeAll,AOSPExit
  CHARACTER(LEN=DEFAULT_CHR_LEN):: Mssg,BName,FileName,DMFile
  CHARACTER(LEN=8)              :: MDGeuss
  CHARACTER(LEN=5),PARAMETER    :: Prog='P2Use'
  !------------------------------------------------------------------------------- 
  ! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
  Cycl=IntToChar(Args%I%I(1))
  !
  NSMat=1
  IF(NAlph.NE.NBeta)NSMat=2
!write(*,*) 'Enter NSMat'
!read(5,*)NSMat
  !
  ! Select spin factor for R/U/G theory.
  SFac=2D0
  IF(NSMat.GT.1)SFac=1D0
  !
  ! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
  ! Do what needs to be done CASE by CASE
  SELECT CASE(SCFActn)
  ! P=0
  CASE('GuessEqCore')
     CALL New(P,NSMat_O=NSMat)
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
     CALL New(P,NSMat_O=NSMat)
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
              ! Rescale for charged molecules.
              IF(TotCh.NE.Zero)BlkP%D(:,I)=NEl*BlkP%D(:,I)/(NEl+TotCh)
           ENDIF
        ENDDO
        !
        IF(NSMat.EQ.1)THEN
           !Set the P
           CALL SetToI(P,BlkP)
        ELSEIF(NSMat.EQ.2)THEN
           !Set the P_alpha
           CALL DSCAL(MaxBlkSize**2*NAtoms,2D0*DBLE(NAlph)/DBLE(NEl)  ,BlkP%D(1,1),1)
           CALL SetToI(P,BlkP,Expert_O=1)
           !Set the P_beta
           CALL DSCAL(MaxBlkSize**2*NAtoms,    DBLE(NBeta)/DBLE(NAlph),BlkP%D(1,1),1)
           CALL SetToI(P,BlkP,Expert_O=2)
        ELSEIF(NSMat.EQ.4)THEN
           !Set the P_alpha
           CALL DSCAL(MaxBlkSize**2*NAtoms,2D0*DBLE(NAlph)/DBLE(NEl)  ,BlkP%D(1,1),1)
           CALL SetToI(P,BlkP,Expert_O=1)
           !Set the off diag bloks P_alpha,beta P_beta,alpha
           Dum=-0.1D0
           CALL DSCAL(MaxBlkSize**2*NAtoms,            Dum/DBLE(NAlph),BlkP%D(1,1),1)
           CALL SetToI(P,BlkP,Expert_O=2)
           CALL SetToI(P,BlkP,Expert_O=3)
           !Set the P_beta
           CALL DSCAL(MaxBlkSize**2*NAtoms,    DBLE(NBeta)/Dum        ,BlkP%D(1,1),1)
           CALL SetToI(P,BlkP,Expert_O=4)
!           BlkP%D=0d0
!           CALL SetToI(P,BlkP,Expert_O=2)
!           CALL SetToI(P,BlkP,Expert_O=3)
        ELSE
           CALL Halt('P2Use: Something wrong when computing P_guess!')
        ENDIF
        !
        ! Check for the correct elctron count
        TrP=Trace(P)
        IF(ABS(TrP-DBLE(NEl/SFac))>1.D-10) &
             & CALL Warn(' In P2Use, TrP = '//TRIM(DblToChar(TrP)))
        CALL Delete(BlkP)
     ELSE
        CALL SetToI(P)
#ifdef PARALLEL
        IF(MyID == ROOT) &
#endif
             CALL Warn('Attempting to use density superpostion with a non STO basis set. Going for scaled I.')
        !
        IF(NSMat.EQ.1)THEN
           !Set the P
           CALL Multiply(P,DBLE(NEl  )/(Two*DBLE(NBasF)))
        ELSEIF(NSMat.EQ.2)THEN
           !Set the P_alpha
           CALL Multiply(P,DBLE(NAlph)/DBLE(NBasF),Expert_O=1)
           !Set the P_beta
           CALL Multiply(P,DBLE(NBeta)/DBLE(NBasF),Expert_O=2)
        ELSEIF(NSMat.EQ.4)THEN
           !Set the P_alpha
           CALL Multiply(P,DBLE(NAlph)/DBLE(NBasF),Expert_O=1)
           !Set the P_beta
           CALL Multiply(P,DBLE(NBeta)/DBLE(NBasF),Expert_O=4)
           !Set the P_ab and P_ba
           CALL Multiply(P,DBLE(NEl  )/DBLE(NBasF),Expert_O=2)
           CALL Multiply(P,DBLE(NEl  )/DBLE(NBasF),Expert_O=3)
        ELSE
           CALL Halt('P2Use: Something wrong when computing P_guess!')
        ENDIF
        !
        TrP=Trace(P)
        IF(ABS(TrP-DBLE(NEl/SFac))>1.D-10) CALL Warn(' In P2Use, TrP = '//TRIM(DblToChar(TrP)))
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
     CALL New(P,NSMat_O=NSMat)
     CALL New(S)
     CALL New(Tmp1)
     CALL New(Tmp2)
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
     ! Get the old AO-dM
     CALL Get(P,'CurrentDM',CheckPoint_O=.TRUE.)
     CALL Get(S,TrixFile('S',Args))
     CALL Multiply(P,S,Tmp1)
     TError0 = ABS(SFac*Trace(Tmp1)-DBLE(NEl))/DBLE(NEl)
     IF(TError0>1D-3) THEN
        !Print out the checksum for P, S and P*S.
        CALL PChkSum(P   ,'CurrentDM',Prog)
        CALL PChkSum(S   ,'Overlap  ',Prog)
        CALL PChkSum(Tmp1,'P*S      ',Prog)
        CALL Halt(' Possible geometry, density matrix mismatch '//Rtrn// &
                  ' Relative error in Tr(S.P)-Nel = '//DblToMedmChar(TError0)//Rtrn// &
                  ' NEl = '//DblToMedmChar(DBLE(NEl))//Rtrn// &
                  ' Try Guess=(Restart,Reparse) with previous geometry ' )
     ENDIF

!     CALL AOSP2(P,S,Tmp1,Tmp2,.TRUE.)
!     CALL AOSP2(P,S,Tmp1,Tmp2,.FALSE.)
!     CALL PChkSum(P,'P['//TRIM(Cycl)//']',Prog,Unit_O=6)

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
     ! Clean Up
     CALL Delete(P)
!     CALL Delete(S)
!     CALL Delete(Tmp1)
!     CALL Delete(Tmp2)
  !  Restarting with BasisSet Change
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
!
! Density Matrix Verlet
!
  CASE('DMVerlet0')
     iGEO = Args%I%I(3)
     CALL New(P)
     CALL New(Tmp1)
     CALL New(Tmp2)
     IF(iGEO .LE. 2) CALL Halt('P2Use:PMVerlet0: No PreviousDensity Matrix Defined') 
!    Save P(p-1,n) as P(p-1,0), where p<3
     IF(iGEO==3) THEN
        DO I=1,2
           FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(I))//'_C#'//TRIM(IntToChar(MyClone))//'.DOsave'
           CALL Get(Tmp1,FileName)
           FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(I))//'_C#'//TRIM(IntToChar(MyClone))//'.DOPsave'  
           CALL Put(Tmp1,FileName)  
        ENDDO
     ENDIF
!    P(p,0) = 2P(p-1,n)-P(p-2,0)
!    Get P(p-1,n)
     FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO-1))//'_C#'//TRIM(IntToChar(MyClone))//'.DOsave'
     CALL Get(Tmp1,FileName)
     CALL Multiply(Tmp1, 2.0D0)
     CALL SetEq(P,Tmp1)
!    Get P(p-2,0)
     FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO-2))//'_C#'//TRIM(IntToChar(MyClone))//'.DOPsave'
     CALL Get(Tmp1,FileName)
     CALL Multiply(Tmp1,-1.0D0)
     CALL Add(P,Tmp1,Tmp2)
     CALL SetEq(P,Tmp2)
!    Save P(p,0)
     FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO  ))//'_C#'//TRIM(IntToChar(MyClone))//'.DOPsave'
     CALL Put(P,FileName)
!    Purify P
#ifdef PARALLEL
     CALL SetEq(P_BCSR,P)
     CALL SpectralBounds(P_BCSR,Fmin,Fmax)
     CALL Delete(P_BCSR)
#else
     CALL SpectralBounds(P,Fmin,Fmax)
#endif
     CALL Add(P,-Fmin)
     CALL Multiply(P,One/(Fmax-Fmin))
     MM = 0
     DO I=1,40
        CALL SP2(P,Tmp1,Tmp2,Half*DBLE(NEl),MM)
        IF(ABS(TrP -Half*DBLE(NEl))< 1.0D-8) THEN
           IF(ABS(TrP2-Half*DBLE(NEl)) < 1.D-8) EXIT
        ENDIF
     ENDDO
     WRITE(*,*) "P2Use: Trace(P) = ",TrP,TrP2
!    Convert to AO Rep
     INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
     IF(Present)THEN
        CALL Get(Tmp1,TrixFile('X',Args))   ! Z=S^(-1/2)
        CALL Multiply(Tmp1,P,Tmp2)
        CALL Multiply(Tmp2,Tmp1,P)
     ELSE
        CALL Get(Tmp1,TrixFile('Z',Args))   ! Z=S^(-L)
        CALL Multiply(Tmp1,P,Tmp2)
        CALL Get(Tmp1,TrixFile('ZT',Args))
        CALL Multiply(Tmp2,Tmp1,P)
     ENDIF
     CALL Filter(Tmp1,P)  
!    Put to Disk
     CALL Put(Tmp1,'CurrentDM',CheckPoint_O=.TRUE.)
     CALL Put(Tmp1,TrixFile('D',Args,0))
     CALL PChkSum(Tmp1,'P[0]',Prog)
!    Clean Up
     CALL Delete(P)
     CALL Delete(Tmp1)
     CALL Delete(Tmp2)
  CASE('DMVerlet1')
     iGEO = Args%I%I(3)
     CALL New(P)
     CALL New(Tmp1)
     CALL New(Tmp2)
     IF(iGEO .LE. 2) CALL Halt('P2Use:PMVerlet0: No PreviousDensity Matrix Defined') 
!
!    Compute Pnew in ortho space: This is specifically for MD
!    Save P(p-1,n) as P(p-1,0), where p<5
     IF(iGEO==5) THEN
        DO I=1,4
           FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(I))//'_C#'//TRIM(IntToChar(MyClone))//'.DOsave'
           CALL Get(Tmp1,FileName)
           FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(I))//'_C#'//TRIM(IntToChar(MyClone))//'.DOPsave'  
           CALL Put(Tmp1,FileName)  
        ENDDO
     ENDIF
!    P(p,0) = 4*P(p-1,n)-6*P(p-2,n)+4*P(p-3,n)-P(p-4,0)
!    Get P(p-1,n)
     FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO-1))//'_C#'//TRIM(IntToChar(MyClone))//'.DOsave'
     CALL Get(Tmp1,FileName)
     CALL Multiply(Tmp1, 4.0D0)
     CALL SetEq(P,Tmp1)
!    Get P(p-2,n)
     FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO-2))//'_C#'//TRIM(IntToChar(MyClone))//'.DOsave'
     CALL Get(Tmp1,FileName)
     CALL Multiply(Tmp1,-6.0D0)
     CALL Add(P,Tmp1,Tmp2)
     CALL SetEq(P,Tmp2)
!    Get P(p-3,n)
     FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO-3))//'_C#'//TRIM(IntToChar(MyClone))//'.DOsave'
     CALL Get(Tmp1,FileName)
     CALL Multiply(Tmp1, 4.0D0)
     CALL Add(P,Tmp1,Tmp2)
     CALL SetEq(P,Tmp2)
!    Get P(p-4,0)
     FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO-4))//'_C#'//TRIM(IntToChar(MyClone))//'.DOPsave'
     CALL Get(Tmp1,FileName)
     CALL Multiply(Tmp1,-1.0D0)
     CALL Add(P,Tmp1,Tmp2)
     CALL SetEq(P,Tmp2)
!    Save P(p,0)
     FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO  ))//'_C#'//TRIM(IntToChar(MyClone))//'.DOPsave'
     CALL Put(P,FileName)
!    Purify P
#ifdef PARALLEL
     CALL SetEq(P_BCSR,P)
     CALL SpectralBounds(P_BCSR,Fmin,Fmax)
     CALL Delete(P_BCSR)
#else
     CALL SpectralBounds(P,Fmin,Fmax)
#endif
     CALL Add(P,-Fmin)
     CALL Multiply(P,One/(Fmax-Fmin))
     MM = 0
     DO I=1,40
        CALL SP2(P,Tmp1,Tmp2,Half*DBLE(NEl),MM)
        IF(ABS(TrP -Half*DBLE(NEl))< 1.0D-8) THEN
           IF(ABS(TrP2-Half*DBLE(NEl)) < 1.D-8) EXIT
        ENDIF
     ENDDO
     WRITE(*,*) "P2Use: Trace(P) = ",TrP,TrP2
!    Convert to AO Rep
     INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
     IF(Present)THEN
        CALL Get(Tmp1,TrixFile('X',Args))   ! Z=S^(-1/2)
        CALL Multiply(Tmp1,P,Tmp2)
        CALL Multiply(Tmp2,Tmp1,P)
     ELSE
        CALL Get(Tmp1,TrixFile('Z',Args))   ! Z=S^(-L)
        CALL Multiply(Tmp1,P,Tmp2)
        CALL Get(Tmp1,TrixFile('ZT',Args))
        CALL Multiply(Tmp2,Tmp1,P)
     ENDIF
     CALL Filter(Tmp1,P)  
!    Put to Disk
     CALL Put(Tmp1,'CurrentDM',CheckPoint_O=.TRUE.)
     CALL Put(Tmp1,TrixFile('D',Args,0))
     CALL PChkSum(Tmp1,'P[0]',Prog)
!    Clean Up
     CALL Delete(P)
     CALL Delete(Tmp1)
     CALL Delete(Tmp2)
  CASE('FMVerlet0')
     iGEO = Args%I%I(3)
     CALL New(F)
     CALL New(P)
     CALL New(Tmp1)
     CALL New(Tmp2)
     IF(iGEO .LE. 2) CALL Halt('P2Use:FMVerlet0: No Previous Fock Matrix Defined') 
!
!    Compute Fnew in ortho space: This is specifically for MD
!
!    Save F(p-1,n) as F(p-1,0), when p==3
     IF(iGEO==3) THEN
        DO I=1,2
           FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(I))//'_C#'//TRIM(IntToChar(MyClone))//'.FOsave'
           CALL Get(Tmp1,FileName)
           FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(I))//'_C#'//TRIM(IntToChar(MyClone))//'.FOPsave'  
           CALL Put(Tmp1,FileName)  
        ENDDO
     ENDIF
!    F(p,0) = 2*F(p-1,n)-F(p-2,0)
!    Get F(p-1,n)
     FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO-1))//'_C#'//TRIM(IntToChar(MyClone))//'.FOsave'
     CALL Get(Tmp1,FileName)
     CALL Multiply(Tmp1, 2.0D0)
     CALL SetEq(F,Tmp1)
!    Get F(p-2,0)
     FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO-2))//'_C#'//TRIM(IntToChar(MyClone))//'.FOPsave'
     CALL Get(Tmp1,FileName)
     CALL Multiply(Tmp1,-1.0D0)
     CALL Add(F,Tmp1,Tmp2)
     CALL SetEq(F,Tmp2)
!    Save F(p,0)
     FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO  ))//'_C#'//TRIM(IntToChar(MyClone))//'.FOPsave'
     CALL Put(F,FileName)
!
!    Solve for the Density Matrix using SP2
!
!    Guess P from F
#ifdef PARALLEL
     CALL SetEq(F_BCSR,F)
     CALL FockGuess(F_BCSR,P,Half*DBLE(NEl),1)
     CALL Delete(F_BCSR)
#else
     CALL FockGuess(F,P,Half*DBLE(NEl),1)
#endif
!    Solve for P
     DO I=1,40
        CALL SP2(P,Tmp1,Tmp2,Half*DBLE(NEl),MM)
        IF(ABS(TrP -Half*DBLE(NEl))< 1.0D-8) THEN
           IF(ABS(TrP2-Half*DBLE(NEl)) < 1.D-8) EXIT
        ENDIF
     ENDDO
     WRITE(*,*) "P2Use: Trace(P) = ",TrP
!    Convert to AO Rep
     INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
     IF(Present)THEN
        CALL Get(Tmp1,TrixFile('X',Args))   ! Z=S^(-1/2)
        CALL Multiply(Tmp1,P,Tmp2)
        CALL Multiply(Tmp2,Tmp1,P)
     ELSE
        CALL Get(Tmp1,TrixFile('Z',Args))   ! Z=S^(-L)
        CALL Multiply(Tmp1,P,Tmp2)
        CALL Get(Tmp1,TrixFile('ZT',Args))
        CALL Multiply(Tmp2,Tmp1,P)
     ENDIF
     CALL Filter(Tmp1,P)  
!    Put to Disk
     CALL Put(Tmp1,'CurrentDM',CheckPoint_O=.TRUE.)
     CALL Put(Tmp1,TrixFile('D',Args,0))
     CALL PChkSum(Tmp1,'P[0]',Prog)
!    Clean Up
     CALL Delete(F)
     CALL Delete(P)
     CALL Delete(Tmp1)
     CALL Delete(Tmp2)
!
!
!
  CASE('FMVerlet1')
     iGEO = Args%I%I(3)
     CALL New(F)
     CALL New(P)
     CALL New(Tmp1)
     CALL New(Tmp2)
     IF(iGEO .LE. 4) CALL Halt('P2Use:FMVerlet0: No Previous Fock Matrix Defined') 
!
!    Compute Fnew in ortho space: This is specifically for MD
!
!    Save F(p-1,n) as F(p-1,0), when p==3
     IF(iGEO==5) THEN
        DO I=1,4
           FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(I))//'_C#'//TRIM(IntToChar(MyClone))//'.FOsave'
           CALL Get(Tmp1,FileName)
           FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(I))//'_C#'//TRIM(IntToChar(MyClone))//'.FOPsave'  
           CALL Put(Tmp1,FileName)  
        ENDDO
     ENDIF
!    F(p,0) = 4*P(p-1,n)-6*P(p-2,n)+4*P(p-3,n)-P(p-4,0)
!    Get F(p-1,n)
     FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO-1))//'_C#'//TRIM(IntToChar(MyClone))//'.FOsave'
     CALL Get(Tmp1,FileName)
     CALL Multiply(Tmp1, 4.0D0)
     CALL SetEq(F,Tmp1)
!    Get F(p-2,0)
     FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO-2))//'_C#'//TRIM(IntToChar(MyClone))//'.FOsave'
     CALL Get(Tmp1,FileName)
     CALL Multiply(Tmp1,-6.0D0)
     CALL Add(F,Tmp1,Tmp2)
     CALL SetEq(F,Tmp2)
!    Get F(p-3,0)
     FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO-3))//'_C#'//TRIM(IntToChar(MyClone))//'.FOsave'
     CALL Get(Tmp1,FileName)
     CALL Multiply(Tmp1, 4.0D0)
     CALL Add(F,Tmp1,Tmp2)
     CALL SetEq(F,Tmp2)
!    Get F(p-4,0)
     FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO-4))//'_C#'//TRIM(IntToChar(MyClone))//'.FOPsave'
     CALL Get(Tmp1,FileName)
     CALL Multiply(Tmp1,-1.0D0)
     CALL Add(F,Tmp1,Tmp2)
     CALL SetEq(F,Tmp2)
!    Save F(p,0)
     FileName = TRIM(SCRName)//'_G#'//TRIM(IntToChar(iGEO  ))//'_C#'//TRIM(IntToChar(MyClone))//'.FOPsave'
     CALL Put(F,FileName)
!
!
!    Solve for the Density Matrix using SP2
!
!    Guess P from F
#ifdef PARALLEL
     CALL SetEq(F_BCSR,F)
     CALL FockGuess(F_BCSR,P,Half*DBLE(NEl),1)
     CALL Delete(F_BCSR)
#else
     CALL FockGuess(F,P,Half*DBLE(NEl),1)
#endif
!    Solve for P
     DO I=1,40
        CALL SP2(P,Tmp1,Tmp2,Half*DBLE(NEl),MM)
        IF(ABS(TrP -Half*DBLE(NEl))< 1.0D-8) THEN
           IF(ABS(TrP2-Half*DBLE(NEl)) < 1.D-8) EXIT
        ENDIF
     ENDDO
     WRITE(*,*) "P2Use: Trace(P) = ",TrP
!    Convert to AO Rep
     INQUIRE(FILE=TrixFile('X',Args),EXIST=Present)
     IF(Present)THEN
        CALL Get(Tmp1,TrixFile('X',Args))   ! Z=S^(-1/2)
        CALL Multiply(Tmp1,P,Tmp2)
        CALL Multiply(Tmp2,Tmp1,P)
     ELSE
        CALL Get(Tmp1,TrixFile('Z',Args))   ! Z=S^(-L)
        CALL Multiply(Tmp1,P,Tmp2)
        CALL Get(Tmp1,TrixFile('ZT',Args))
        CALL Multiply(Tmp2,Tmp1,P)
     ENDIF
     CALL Filter(Tmp1,P)  
!    Put to Disk
     CALL Put(Tmp1,'CurrentDM',CheckPoint_O=.TRUE.)
     CALL Put(Tmp1,TrixFile('D',Args,0))
     CALL PChkSum(Tmp1,'P[0]',Prog)
!    Clean Up
     CALL Delete(F)
     CALL Delete(P)
     CALL Delete(Tmp1)
     CALL Delete(Tmp2)
!
!    Geometry Change 
!
  CASE('Extrapolate')
     ! Allocate
     CALL New(P)
     CALL New(P0)
     CALL New(S)
     CALL New(S0)
     CALL New(S1)
     CALL New(Tmp1)
     CALL New(Tmp2)
     ! Get Marices

     DO ICycle=1,1000
       DMFile=TRIM(SCRName)//'_Geom#'//TRIM(IntToChar(Current(3)-1)) &
                           //'_Base#'//TRIM(IntToChar(Current(2))) &
                           //'_Cycl#'//TRIM(IntToChar(ICycle)) &
                           //'_Clone#'//TRIM(IntToChar(MyClone)) &
                           //'.D'
!       WRITE(*,*)' Looking for ',TRIM(DMFile)
        INQUIRE(FILE=DMFile,EXIST=Present)
!       WRITE(*,*)' Found it? ',Present
        IF(.NOT.Present)THEN
           Cycle=ICycle-1		
           DMFile=TRIM(SCRName)//'_Geom#'//TRIM(IntToChar(Current(3)-1)) &
                               //'_Base#'//TRIM(IntToChar(Current(2))) &
                               //'_Cycl#'//TRIM(IntToChar(Cycle)) &
                               //'_Clone#'//TRIM(IntToChar(MyClone)) &
                               //'.D'
           CALL Logger('On extrapolation, P2Use is openning DM '//TRIM(DMFile),.TRUE.)
           EXIT
        ENDIF
     ENDDO	

     IF(Cycle<0)THEN
	WRITE(*,*)' Assuming this is a restart!  If there is no restart, its going to die...'
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
        CALL Get(S1,TrixFile('S',Args,Stats_O=Current ))
        DMFile=TRIM(SCRName)//'_Geom#'//TRIM(IntToChar(Current(3)-1)) &
                            //'_Base#'//TRIM(IntToChar(Current(2))) &
                            //'_Cycl#'//TRIM(IntToChar(Cycle)) &
                            //'_Clone#'//TRIM(IntToChar(MyClone)) &
                            //'.D'
        CALL Get(P0,DMFile)
        CALL Get(DoingMD ,'DoingMD')
        IF(DoingMD) THEN
           DMPorder=0
           CALL Get(MDGeuss ,"MDGeuss")
           IF(MDGeuss=='DMProj0') DMPorder=0
           IF(MDGeuss=='DMProj1') DMPorder=1
           IF(MDGeuss=='DMProj2') DMPorder=2
           IF(MDGeuss=='DMProj3') DMPorder=3
           IF(MDGeuss=='DMProj4') DMPorder=4
           iGEO     = Args%I%I(3)
           DMPOrder = MIN(MAX(iGEO-2,0),DMPOrder)
           CALL DMPProj(iGEO,DMPOrder,P0,Tmp1,Tmp2)
        ELSE
           IF(DMPorder > 0) THEN
              CALL Warn('P2Use:DMPOrder is only implimented for MD')
           ENDIF
        ENDIF
     ENDIF
     ! Initial Trace Error
#ifdef PARALLEL
     CALL Multiply(P0,S0,Tmp1)
     TError0 = ABS(Trace(Tmp1)-DBLE(NEl)/SFac)
     IF(MyID.EQ.ROOT) THEN 
        IF(PrintFlags%Key==DEBUG_MAXIMUM) THEN
           WRITE(*,*) "Trace Error: Tr[P0,S0] = ",TError0
           CALL OpenASCII(OutFile,Out)
           CALL PrintProtectL(Out)
           WRITE(oUT,*) "Trace Error: Tr[P0,S1] = ",TError0 
           CALL PrintProtectR(Out)
           CLOSE(UNIT=Out,STATUS='KEEP')
        ENDIF
     ENDIF
     CALL Multiply(P0,S1,Tmp1)
     TError0 = ABS(Trace(Tmp1)-DBLE(NEl)/SFac)
     IF(MyID.EQ.ROOT)THEN 
        IF(PrintFlags%Key==DEBUG_MAXIMUM) THEN 
           WRITE(*,*) "Trace Error: Tr[P0,S1] = ",TError0 
           CALL OpenASCII(OutFile,Out)
           CALL PrintProtectL(Out)
           WRITE(oUT,*) "Trace Error: Tr[P0,S1] = ",TError0 
           CALL PrintProtectR(Out)
           CLOSE(UNIT=Out,STATUS='KEEP')
        ENDIF
     ENDIF
#else
     TError0 = ABS(Trace(P0,S1)-DBLE(NEl)/SFac)
     IF(PrintFlags%Key==DEBUG_MAXIMUM) THEN
        TrP=ABS(Trace(P0,S0)-DBLE(NEl)/SFac)
        WRITE(*,*) "Trace Error: Tr[P0,S0] = ",TrP
        TrP=ABS(Trace(P0,S1)-DBLE(NEl)/SFac)
        WRITE(*,*) "Trace Error: Tr[P0,S1] = ",TrP
        CALL OpenASCII(OutFile,Out)
        CALL PrintProtectL(Out)
        TrP=ABS(Trace(P0,S0)-DBLE(NEl)/SFac)
        WRITE(Out,*) "Trace Error: Tr[P0,S0] = ",TrP
        TrP=ABS(Trace(P0,S1)-DBLE(NEl)/SFac)
        WRITE(Out,*) "Trace Error: Tr[P0,S1] = ",TrP
        CALL PrintProtectR(Out)
        CLOSE(UNIT=Out,STATUS='KEEP')
     ENDIF
#endif
!    Initialize
     NStep       = 0
     Lam         = Zero
     DLam        = One
     ConvergeAll = .FALSE.
!    Purify 
     DO 
        NStep = NStep + 1 
        Lam  = Lam + DLam
!       Set Up S and P        
        CALL SetEq(Tmp1,S0)
        CALL SetEq(Tmp2,S1)
        CALL Multiply(Tmp1,One-Lam)
        CALL Multiply(Tmp2,Lam)
        CALL Add(Tmp1,Tmp2,S)
        CALL SetEq(P,P0) 
!
#ifdef PARALLEL 
        IF(MyId==ROOT)THEN
#endif
           Mssg=ProcessName(Prog,'AO-DMX ') &
                //' Nstep  = '//TRIM(IntToChar(NStep)) &
                //' Lambda = '//TRIM(DblToMedmChar(Lam))      
           CALL OpenASCII(OutFile,Out)
           CALL PrintProtectL(Out)
           IF(PrintFlags%Key>=DEBUG_MEDIUM) THEN
              WRITE(*,*)TRIM(Mssg)
              WRITE(Out,*)TRIM(Mssg)
           ENDIF
           CALL PrintProtectR(Out)
           CLOSE(UNIT=Out,STATUS='KEEP')
#ifdef PARALLEL          
        ENDIF
#endif
!
#ifdef PARALLEL
        CALL Multiply(P,S,Tmp1)
        TrP        = Trace(Tmp1)
#else
        TrP        = Trace(P,S)
#endif
        Norm_Error = TrP-DBLE(NEl)/SFac
        Ipot_Error = One
!
        Qstep        = 0
        ConvergeAOSP = .FALSE.
        AOSPExit     = .FALSE.
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
           Norm_Error = TrP-DBLE(NEl)/SFac
           Ipot_Error = TrP2-Trace(P)
#ifdef PARALLEL
           IF(MyId==ROOT)THEN
#endif
              PNon0s=100.D0*DBLE(P%NNon0)/DBLE(NBasF*NBasF)   
              Mssg=ProcessName(Prog,'AO-DMX '//TRIM(IntToChar(NStep))//":"//TRIM(IntToChar(I))) &
                   //' dN='//TRIM(DblToShrtChar(Norm_Error))       &
                   //', %Non0='//TRIM(DblToShrtChar(PNon0s))              
              CALL OpenASCII(OutFile,Out)
              CALL PrintProtectL(Out)
              IF(PrintFlags%Key>=DEBUG_MEDIUM) THEN
                 WRITE(*,*)TRIM(Mssg)
                 WRITE(Out,*)TRIM(Mssg)
              ENDIF
              CALL PrintProtectR(Out)
              CLOSE(UNIT=Out,STATUS='KEEP')
#ifdef PARALLEL          
           ENDIF
#endif
!          Logic
           IF(ABS(Ipot_Error) < 1.0D-10 .AND. ABS(Norm_Error) < 1.0D-10) THEN
              AOSPExit     = .TRUE.
              ConvergeAOSP = .TRUE.
              IF(Lam > (One-1.D-14)) ConvergeAll = .TRUE.
           ELSEIF(Qstep > 5) THEN
              AOSPExit  = .TRUE.
              IF(ABS(Ipot_Error) < 1.0D-4 .AND. ABS(Norm_Error) < 1.0D-4) ConvergeAOSP = .TRUE.
              IF(ConvergeAOSP) THEN
                 IF(Lam > (One-1.D-14)) ConvergeAll = .TRUE.
              ENDIF
           ELSEIF(ABS(Norm_Error) > 100.D0*TError0) THEN
              AOSPExit     = .TRUE.
           ENDIF
           IF(AOSPExit) EXIT
        ENDDO
!       Logic
        IF(ConvergeAll) THEN
           CALL SetEq(P0,P)
           EXIT
        ELSE
           IF(DLam < 1.D-2 .OR. NStep > 50) THEN
              CALL Halt('P2Use:Density Matrix Extrapolater Failed to converge in '//TRIM(IntToChar(NStep))//' steps')
           ENDIF
           IF(.NOT. ConvergeAOSP) THEN
              Lam  = Lam - DLam
              DLam = Half*DLam
           ELSE
              CALL SetEq(P0,P)
           ENDIF
        ENDIF
     ENDDO
     ! Save back to be sure.
     CALL PChkSum(P,'P[0]',Prog)
     CALL Put(P,TrixFile('D',Args,0))
     CALL Put(P,'CurrentDM',CheckPoint_O=.TRUE.) 
     ! Clean Up
     CALL Delete(P)
     CALL Delete(P0)
     CALL Delete(S)
     CALL Delete(S0)
     CALL Delete(S1)
     CALL Delete(Tmp1)
     CALL Delete(Tmp2)
  CASE DEFAULT
     CALL Halt(' Unknown option '//TRIM(SCFActn))
  END SELECT
  ! Tidy up ...
  CALL Delete(GM)
  CALL Delete(BS)
  CALL ShutDown(Prog)   
!
END PROGRAM P2Use

