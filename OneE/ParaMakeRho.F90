! COMPUTE THE DENSITY IN A HGTF BASIS FROM THE DENSITY MATRIX
! IN PARALLEL
! C. K. Gan, modify MakeRho.f90
! deleted MP stuff

PROGRAM ParaMakeRho
#ifdef PARALLEL
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
  USE BraBloks
  USE RhoBlok
  USE RhoTools
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef PARALLEL
  TYPE(DBCSR)                    :: Dmat,D1,D2
  INTEGER                        :: LocalAtom
  REAL(DOUBLE)                   :: TotRSumE,TotRSumN
#else
  TYPE(BCSR)                     :: Dmat,D1,D2
#endif
#ifdef PERIODIC
  INTEGER                        :: NC
  REAL(DOUBLE),DIMENSION(3)      :: B
#endif
  TYPE(AtomPair)                 :: Pair
  TYPE(BSET)                     :: BS
  TYPE(CRDS)                     :: GM,GM_MM
  TYPE(DBL_RNK4)                 :: MD
  TYPE(ARGMT)                    :: Args
  TYPE(HGRho)                    :: Rho
  TYPE(HGRho_new)                :: RhoA
  TYPE(INT_VECT)                 :: Stat
  INTEGER                        :: P,R,AtA,AtB,NN,iSwitch,IC1,IC2
  INTEGER                        :: NDist,NCoef,I,J,Pbeg,Pend,NDist_old,NDist_new,NumAtoms
  INTEGER                        :: N1,N2,QMOffSetQ,QMOffSetR,PcntDist,OldFileID
  REAL(DOUBLE)                   :: DistThresh,RSumE,RSumN,RSumMM,RelRhoErr, &
       QMCharge,dQMCharge,MMCharge,dMMCharge,PcntCharge
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg1,Mssg2,RestartHDF
  CHARACTER(LEN=11),PARAMETER    :: Prog='ParaMakeRho'
  !-------------------------------------------------------------------------------
  ! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
#ifdef PERIODIC
#ifdef PARALLEL_CLONES
#else
  ! Get the Outer Cell Set
  CALL Get(CS_OUT,'CS_OUT',Tag_O=CurBase)
#endif
#endif


  IF(SCFActn=='BasisSetSwitch')THEN
     !  Get the previous information
     CALL Get(BS,PrvBase)
     CALL Get(GM,CurGeom)
     CALL SetThresholds(PrvBase)
     CALL Get(BSiz,'atsiz',PrvBase)
     CALL Get(OffS,'atoff',PrvBase)
     CALL Get(NBasF,'nbasf',PrvBase)
     CALL Get(Dmat,TrixFile('D',Args,-1))
  ELSEIF(SCFActn=='Restart')THEN
     CALL Get(GM,CurGeom)
#ifdef PARALLEL_CLONES
     ! Close current group and HDF
     CALL CloseHDFGroup(H5GroupID)
     CALL CloseHDF(HDFFileID)
     ! Open old group and HDF
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
     CALL BCast(BSiz)
     CALL BCast(OffS)
     ! Close the old hdf up 
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(OldFileID)
     ! Reopen current group and HDF
     HDFFileID=OpenHDF(H5File)
     H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
     HDF_CurrentID=H5GroupID
#else
     ! Get the old information
     CALL Get(RestartHDF,'OldInfo')
     CALL CloseHDF(HDF_CurrentID)
     HDF_CurrentID=OpenHDF(RestartHDF)
     CALL New(Stat,3)
     CALL Get(Stat,'current')
     SCFCycl=TRIM(IntToChar(Stat%I(1)))
     CurBase=TRIM(IntToChar(Stat%I(2)))
     CurGeom=TRIM(IntToChar(Stat%I(3)))
     CALL Get(BS,CurBase)
     CALL Get(GM,CurGeom)
     CALL CloseHDF(HDF_CurrentID)
     HDF_CurrentID=OpenHDF(InfFile)     
#endif
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
     ELSEIF(SCFActn/='Core')THEN
        CALL Get(Dmat,TrixFile('D',Args,0))
     ENDIF
     ENDIF
!  Allocations and precalculations
   CALL NewBraBlok(BS)  
!--------------------------------------------------------------
! Main loops: First pass calculates the size.
!             Second pass calculates the density
!-------------------------------------------------------------
!  Initailize 
   NDist = 0
   NCoef = 0
!  Loop over atoms and count primatives
#ifdef PARALLEL
   DO LocalAtom = 1, Dmat%NAtms
      AtA  = Beg%I(MyID)+(LocalAtom-1)
      Pbeg = Dmat%RowPt%I(LocalAtom)
      Pend = Dmat%RowPt%I(LocalAtom+1)-1
#else
   DO AtA=1,NAtoms
      Pbeg = Dmat%RowPt%I(AtA)
      Pend = Dmat%RowPt%I(AtA+1)-1
#endif
      DO P = Pbeg,Pend
         AtB = Dmat%ColPt%I(P)
         IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
#ifdef PERIODIC
            B = Pair%B
            DO NC = 1,CS_OUT%NCells
               Pair%B = B+CS_OUT%CellCarts%D(:,NC)
               Pair%AB2  = (Pair%A(1)-Pair%B(1))**2 &
                    + (Pair%A(2)-Pair%B(2))**2 &
                    + (Pair%A(3)-Pair%B(3))**2
               IF(TestAtomPair(Pair)) THEN
                  CALL PrimCount(BS,Pair,NDist,NCoef)
               ENDIF
            ENDDO
#else
            CALL PrimCount(BS,Pair,NDist,NCoef)
#endif
         ENDIF
      ENDDO
   ENDDO
!  Allocate the Density
   CALL New_HGRho_new(RhoA,(/NDist,NCoef/))
!  Initailize  Counters
   NDist        = 0
   NCoef        = 0
!  Loop over atoms and calculate the electronic density
#ifdef PARALLEL
   DO LocalAtom = 1, Dmat%NAtms
      AtA  = Beg%I(MyID)+(LocalAtom-1)
      Pbeg = Dmat%RowPt%I(LocalAtom)
      Pend = Dmat%RowPt%I(LocalAtom+1)-1
#else
   DO AtA=1,NAtoms
      Pbeg = Dmat%RowPt%I(AtA)
      Pend = Dmat%RowPt%I(AtA+1)-1
#endif
      DO P = Pbeg,Pend
         AtB = Dmat%ColPt%I(P)
         R   = Dmat%BlkPt%I(P)
         IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN 
#ifdef PERIODIC
            B = Pair%B
            DO NC = 1,CS_OUT%NCells
               Pair%B = B+CS_OUT%CellCarts%D(:,NC)
               Pair%AB2  = (Pair%A(1)-Pair%B(1))**2 &
                    + (Pair%A(2)-Pair%B(2))**2 &
                    + (Pair%A(3)-Pair%B(3))**2
               IF(TestAtomPair(Pair)) THEN
                  NN = Pair%NA*Pair%NB
                  CALL RhoBlk(BS,Dmat%MTrix%D(R:R+NN-1),Pair,NDist,NCoef,RhoA)
               ENDIF
            ENDDO
#else
            NN = Pair%NA*Pair%NB
            CALL RhoBlk(BS,Dmat%MTrix%D(R:R+NN-1),Pair,NDist,NCoef,RhoA)
#endif
         ENDIF
      ENDDO
   ENDDO
!  Add in the density for the nuclear centers
   IF(SCFActn/='InkFok') THEN
      CALL AddDist(RhoA,GM,NuclearExpnt,Beg%I(MyID),End%I(MyID))
   ENDIF
!  Prune negligible distributions from the electronic density
   NDist_old = RhoA%NDist
   CALL Prune_Rho_new(Thresholds%Dist,RhoA) 
   NDist_new = RhoA%NDist
!  Compute integrated electron and nuclear densities
#ifdef PERIODIC
   CALL Fold_Rho_new(GM,RhoA)
#endif
!  Compute the Electron and Nuclear Densities
   NumAtoms = End%I(MyID)-Beg%I(MyID)+1
   RSumE    = Integrate_HGRho_new(RhoA,1,RhoA%NDist-NumAtoms)
   RSumN    = Integrate_HGRho_new(RhoA,RhoA%NDist-NumAtoms+1,RhoA%NDist)
   !   WRITE(*,*) 'MyID = ', MyID, ' RSumE = ', RSumE, ' RSumN = ',RSumN
   TotRSumE = AllReduce(RSumE)
   TotRSumN = AllReduce(RSumN)
   RSumE = TotRSumE
   RSumN = TotRSumN
!
  QMCharge=Zero
  dQMCharge=Zero
  MMCharge=Zero
  dMMCharge=Zero  
!
  QMCharge=Two*(RSumE+RSumN)
  dQMCharge=QMCharge+GM%TotCh 
!
  PcntDist=FLOOR(1.D2*DBLE(NDist_new)/DBLE(NDist_old))
!
  Mssg1=TRIM(Mssg1)//' dNel = '//TRIM(DblToShrtChar(dQMCharge))//', kept '  &
          //TRIM(IntToChar(PcntDist))//'% of distributions.'
  RelRhoErr=ABS(dQMCharge)/DBLE(NEl)
! Check error
  IF(RelRhoErr>Thresholds%Dist*5.D3.AND.SCFActn/='NumForceEvaluation')THEN
     IF(SCFActn=='Restart')THEN
        CALL Warn(ProcessName(Prog)//'relative error in density = '//TRIM(DblToShrtChar(RelRhoErr)))
     ELSE
        WRITE(*,*) 'MyID = ', MyID, ' prepare to crash!!!! Halt is to be called'
        CALL Halt(ProcessName(Prog)//'relative error in density = '//TRIM(DblToShrtChar(RelRhoErr)) &
             //'. Distribution threshold = '//TRIM(DblToShrtChar(Thresholds%Dist))      &
             //'. Total charge lost = '//TRIM(DblToShrtChar(dQMCharge+dMMCharge)))
     ENDIF
  ENDIF

  CALL ConvertToOldRho(Rho,RhoA)

! Put Rho and MPs to disk
  IF(SCFActn=='ForceEvaluation')THEN
#ifdef PARALLEL
     CALL Put_HGRho(Rho,'Rho'//IntToChar(MyID),Args,1)
#else
     CALL Put_HGRho(Rho,'Rho',Args,1) 
#endif

  ELSEIF(SCFActn=='InkFok')THEN
     CALL Put_HGRho(Rho,'DeltaRho',Args,0)
  ELSE

#ifdef PARALLEL
     CALL Put_HGRho(Rho,'Rho'//IntToChar(MyID),Args,0)
#else
     CALL Put_HGRho(Rho,'Rho',Args,0)
#endif
  ENDIF
! Tidy up

     CALL Delete(GM)
     CALL Delete(Dmat)
     CALL Delete(BS)
     CALL DeleteBraBlok()

  CALL Delete_HGRho(Rho)
  CALL Delete_HGRho_New(RhoA)
  CALL ShutDown(Prog)
#endif
END PROGRAM ParaMakeRho
