!    FAST O(N lg N) COMPUTATION OF THE COULOMB MATRIX
!    Authors:  Matt Challacombe and CJ Tymczak
!===============================================================================
PROGRAM QCTC
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE ProcessControl
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE AtomPairs
  USE BraBloks
  USE QCTCThresholds
  USE PoleTree
  USE Globals
  USE PBCFarField
  USE JGen
  USE NukularE
  USE Density
  USE Clock
  USE TreeWalk
  IMPLICIT NONE
  TYPE(BCSR)                     :: J,  P
  TYPE(BCSR)                     :: T1,T2
  TYPE(BCSR)                     :: Dmat,D1,D2
  TYPE(INT_VECT)                 :: Stat
  TYPE(HGLL),POINTER             :: RhoHead
  REAL(DOUBLE)                   :: E_Nuc_Tot,Etot,SdvErrorJ,MaxErrorJ,JMExact,E_PFF,D_PFF
  REAL(DOUBLE)                   :: QCTC_TotalTime_Start
  TYPE(TIME)                     :: TimeMakeJ,TimeMakeTree,TimeNukE
  CHARACTER(LEN=4),PARAMETER     :: Prog='QCTC'
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
  INTEGER                        :: I,K,NLink,OldFileID
  !
!!$  LOGICAL                        :: NoWrap=.TRUE.  ! WRAPPING IS OFF
  LOGICAL                        :: NoWrap=.FALSE. ! WRAPPING IS ON
  !------------------------------------------------------------------------------- 
  QCTC_TotalTime_Start=MTimer()
  ! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.TRUE.)
  ! Chose a density matrix
  IF(SCFActn=='BasisSetSwitch')THEN
     ! Get the previous information
     CALL Get(BS,PrvBase)
     CALL Get(GM,CurGeom)
     CALL SetThresholds(PrvBase)
     CALL Get(BSiz,'atsiz',PrvBase)
     CALL Get(OffS,'atoff',PrvBase)
     CALL Get(NBasF,'nbasf',PrvBase)
     CALL Get(Dmat,TrixFile('D',Args,-1))
  ELSEIF(SCFActn=='Restart'.OR. SCFActn=='RestartBasisSwitch')THEN
     ! Get the current geometry from the current HDF first
     CALL Get(GM,CurGeom)
     ! then close current group and HDF
     CALL CloseHDFGroup(H5GroupID)
     CALL CloseHDF(HDFFileID)
     ! and open the old group and HDF
     HDF_CurrentID=OpenHDF(Restart)
     OldFileID=HDF_CurrentID
     CALL New(Stat,3)
     CALL Get(Stat,'current_state')
     HDF_CurrentID=OpenHDFGroup(HDF_CurrentID,"Clone #"//TRIM(IntToChar(MyClone)))
     ! Get old basis set stuff
     SCFCycl=TRIM(IntToChar(Stat%I(1)))
     CurBase=TRIM(IntToChar(Stat%I(2)))
     CurGeom=TRIM(IntToChar(Stat%I(3)))
     CALL Get(BS,CurBase)
     ! Compute a sparse matrix blocking scheme for the old BS
     CALL BlockBuild(GM,BS,BSiz,OffS)
     NBasF=BS%NBasF
     ! Close the old hdf up 
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(OldFileID)
     ! Reopen current group and HDF
     HDFFileID=OpenHDF(H5File)
     H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
     HDF_CurrentID=H5GroupID
     CALL Get(Dmat,TrixFile('D',Args,0))
  ELSE
     ! Get the current information
     CALL Get(BS,CurBase)
     CALL Get(GM,CurGeom)
     IF(SCFActn=='InkFok')THEN
        CALL Get(D1,TrixFile('D',Args,-1))
        CALL Get(D2,TrixFile('D',Args,0))
        CALL Multiply(D1,-One)
        CALL Add(D1,D2,DMat)
        IF(HasHF(ModelChem))THEN                
           CALL Filter(D1,DMat)
           CALL Put(D1,TrixFile('DeltaD',Args,0))
        ENDIF
        CALL Delete(D1)
        CALL Delete(D2)
     ELSEIF(SCFActn=='ForceEvaluation')THEN
        CALL Get(Dmat,TrixFile('D',Args,1))
     ELSEIF(SCFActn=='StartResponse')THEN
        CALL Halt('MakeRho: SCFActn cannot be equal to <StartResponse>')
     ELSEIF(SCFActn=='DensityPrime')THEN
        CALL Get(Dmat,TrixFile('DPrime'//TRIM(Args%C%C(3)),Args,0))
     ELSEIF(SCFActn/='Core')THEN
        ! Default
        CALL Get(Dmat,TrixFile('D',Args,0))
     ENDIF
  ENDIF
 

!!  CALL Get(P,TrixFile('D',Args,0))
!!$!  CALL PPrint(GM,Unit_O=6)
!!$!  CALL PPrint(GM%PBC,Unit_O=6)

  ! Allocate some memory for bra HG shenanigans 
  CALL NewBraBlok(BS)
  ! Set thresholds local to QCTC (for PAC and MAC)
  CALL SetLocalThresholds(Thresholds%TwoE)
  ! RhoHead is the start of a linked density list
  ALLOCATE(RhoHead)
  RhoHead%LNum=0
!!$
  ! Here, the LL is filled out
  CALL MakeRhoList(GM,BS,DMat,NLink,RhoHead,'QCTC',NoWrap_O=NoWrap)
  ! Add in the nuclear charges only in certain cases
  IF(SCFActn/='InkFok'.AND.SCFActn/='StartResponse'.AND.SCFActn/='DensityPrime')THEN
     CALL AddNukes(GM,RhoHead,NoWrap)
     NLink=NLink+GM%NAtms
  ENDIF
  ! Load density into arrays and delete the linked list
  CALL Collate(GM,RhoHead,Rho,'QCTC',RhoPoles,NLink)
!
  WRITE(*,*)' PUNTED IN DELETE OF DENSITY !!! SEE FOLLOWING COMMENT OUT :::: '
  !  CALL DeleteHGLL(RhoHead)
  ! Allocate and compute multipole moments of the density
  !
  MaxPFFFEll=GM%PBC%PFFMaxEll
  ! Local expansion order of the multipoles to use in the tree
  MaxPoleEll=MIN(2*BS%NASym+6,MaxPFFFEll)
  IF(MaxPoleEll<2*(BS%NASym+1)) &
     CALL Halt('Bombed in QCTC. Please set PFFMaxEll larger ')
  ! Find the total energy from past calculations
  CALL Get(Etot,'Etot')
  !
  IF(NoWrap)THEN
     Mssg=ProcessName('QCTC','No wrap')
  ELSE
     Mssg=ProcessName('QCTC','Wrapping on')
  ENDIF
  Mssg=TRIM(Mssg)//' Gaussians in density = '//TRIM(IntToChar(NLink))//', Cluster Size = '//TRIM(IntToChar(ClusterSize)) &
           //', MaxPFFEll = '//TRIM(IntToChar(MaxPFFFEll))
  WRITE(*,*)TRIM(Mssg)
  ! Initialize addressing for tensor contraction loops
  CALL TensorIndexingSetUp()
  ! Setup global arrays for computation of multipole tensors ...
  CALL MultipoleSetUp()
  ! Initialize some counters
  ! This preliminary density shit is very sloppy...
  MaxTier=0
  RhoLevel=0
  PoleNodes=0
  ! Initialize the root node   
  CALL NewPoleNode(PoleRoot,0)
  ! Initialize the auxiliary density arrays
  CALL InitRhoAux
  ! Here we set the max multipole expansion relative
  ! to the max angular symmetry
  NInts=0
  NPrim=0
  NFarAv=0
  NNearAv=0
  ! Build the global PoleTree representation of the total density
  CALL RhoToPoleTree
  ! Set up the crystal field     
  CALL PBCFarFuckingFieldSetUp(GM,Rho,'QCTC',MaxPFFFEll,ETot)
  ! Delete the auxiliary density arrays
  CALL DeleteRhoAux
  ! Delete the Density
!  CALL Delete(Rho)
  ! Allocate J
  CALL New(J)
  ! Compute the Coulomb matrix J in O(N Lg N)
  CALL Elapsed_Time(TimeMakeJ,'Init')

  JWalk_Time=0D0
  Integral_Time=0D0
  Multipole_Time=0D0

  ALLOCATE(NNearCount(1:CS_IN%NCells))
  NNearCount=0D0

!  WRITE(*,*)' Poles(1) = ',RhoPoles%DPole%D
!  NewPole=0D0

  CALL MakeJ(J,NoWrap_O=NoWrap)

  K=0
  DO I=1,CS_IN%NCells
!     WRITE(*,*)I,NNearCount(I)
     IF(NNearCount(I)==0D0)K=K+1
  ENDDO
  WRITE(*,*)' % of NoPAC = ',DBLE(K)/DBLE(CS_IN%NCells)

  WRITE(*,11)' Decompos_Time = ',Decompose_Time
  WRITE(*,11)' TreeMake_Time = ',TreeMake_Time
  WRITE(*,11)' JWalking_Time = ',JWalk_Time
  WRITE(*,11)' Integral_Time = ',Integral_Time
  WRITE(*,11)' Multipol_Time = ',Multipole_Time
  WRITE(*,11)' Total J Time  = ',Decompose_Time+TreeMake_Time+JWalk_Time+Multipole_Time+Integral_Time
  WRITE(*,11)' Total JWalks  = ',DBLE(NPrim)
  WRITE(*,11)' Av  Ints/Prim = ',DBLE(NInts)/DBLE(NPrim)
  WRITE(*,11)' Av  # NF/Prim = ',DBLE(NNearAv)/DBLE(NPrim)
  WRITE(*,11)' Av  # FF/Prim = ',DBLE(NFarAv)/DBLE(NPrim)
  WRITE(*,11)' Time per INode= ',Integral_Time/DBLE(NNearAv)
  WRITE(*,11)' Time per MNode= ',Multipole_Time/DBLE(NFarAv)
11 FORMAT(A20,D12.6)
!!$
   CALL PChkSum(J,'J',Prog,Unit_O=6)

  CALL Elapsed_TIME(TimeMakeJ,'Accum')
  IF(SCFActn=='InkFok')THEN
     ! Add in correction if incremental J build
     CALL New(T1)
     CALL New(T2)
     CALL Get(T1,TrixFile('J',Args,-1))
     CALL Add(T1,J,T2)
     CALL Filter(T1,T2)
     CALL Delete(T2)
  ELSE
     CALL Filter(T1,J)
  ENDIF
  ! Put J to disk
  IF(SCFActn=='FockPrimeBuild'.OR.SCFActn=='StartResponse')THEN
     CALL Put(T1,TrixFile('JPrime'//TRIM(Args%C%C(3)),Args,0))
  ELSE
     CALL Put(T1,TrixFile('J',Args,0))
  ENDIF
  ! Compute the nuclear-total electrostatic energy in O(N Lg N)
  IF(SCFActn=='InkFok')THEN
     CALL Get(E_Nuc_Tot,'E_NuclearTotal')
     CALL Elapsed_Time(TimeNukE,'Init')
     E_Nuc_Tot=E_Nuc_Tot+NukE(GM,Rho,NoWrap)
     CALL Elapsed_Time(TimeNukE,'Accum')
  ELSE     
     CALL Elapsed_Time(TimeNukE,'Init')
     E_Nuc_Tot=NukE(GM,Rho,NoWrap)
     CALL Elapsed_Time(TimeNukE,'Accum')
  ENDIF

!  WRITE(*,*)' E_EL-TOT = ',Trace(DMat,T1)    
!  WRITE(*,*)' E_NUC-TOT = ',E_Nuc_Tot
  Mssg=ProcessName(Prog, ' ')
  Mssg=TRIM(Mssg)//' Coulomb Energy      = <'//TRIM(DblToChar(E_Nuc_Tot+Trace(DMat,T1)))//'>'
  WRITE(*,*)TRIM(Mssg)
  
  CALL Put(E_Nuc_Tot,'E_NuclearTotal',StatsToChar(Current))
  !-------------------------------------------------------------------------------
  ! Printing



  IF(SCFActn=='FockPrimeBuild'.OR.SCFActn=='StartResponse')THEN
     CALL PChkSum(T1,'J'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']',Prog)
     CALL PPrint( T1,'J'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']')
     CALL Plot(   T1,'J'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']')
  ELSE
     CALL PChkSum(T1,'J['//TRIM(SCFCycl)//']',Prog)
     CALL PPrint( T1,'J['//TRIM(SCFCycl)//']')!,Unit_O=6)
     CALL Plot(   T1,'J['//TRIM(SCFCycl)//']')
  ENDIF
  ! Print Periodic Info
!  CALL Print_Periodic(GM,Prog)
  CALL Delete(J)
  CALL Delete(T1)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(Args)
  CALL Delete(RhoPoles)
  ! didn't count flops, any accumulation is residual from matrix routines
  PerfMon%FLOP=Zero 
  ! Shutdown 
!!$  WRITE(*,11)' QCTC Total Time = ',MTimer()-QCTC_TotalTime_Start
!!$  CALL OpenASCII("Times.dat",111)
!!$  WRITE(111,12)DBLE(NAtoms),DBLE(NInts)/DBLE(NPrim),DBLE(NNearAv)/DBLE(NPrim),DBLE(NFarAv)/DBLE(NPrim), &
!!$               Decompose_Time,TreeMake_Time,JWalk_Time,Integral_Time,Multipole_Time,MTimer()-QCTC_TotalTime_Start
!!$12 FORMAT(100(" ",D12.6))
!!$  CLOSE(111)
!!$
  Mssg=ProcessName('QCTC','Timing')
  Mssg=TRIM(Mssg)//' Total='//TRIM(DblToMedmChar(MTimer()-QCTC_TotalTime_Start)) & 
                //'; Bisect='//TRIM(DblToShrtChar(Decompose_Time))//', Tree='//TRIM(DblToShrtChar(TreeMake_Time)) 
  WRITE(*,*)TRIM(Mssg)
  Mssg=ProcessName('QCTC','Timing')
  Mssg=TRIM(Mssg)//' Walk='//TRIM(DblToShrtChar(JWalk_Time))//', Ints='//TRIM(DblToShrtChar(Integral_Time)) &
      //', Mults='//TRIM(DblToShrtChar(Multipole_Time)) 
  WRITE(*,*)TRIM(Mssg)

  CALL ShutDown(Prog)
END PROGRAM QCTC
