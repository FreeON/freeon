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
  TYPE(DBCSR)               :: Dmat,D1,D2
  INTEGER                   :: LocalAtom
  REAL(DOUBLE)              :: TotRSumE,TotRSumN
#else
  TYPE(BCSR)                :: Dmat,D1,D2
#endif

#ifdef PERIODIC
  INTEGER                   :: NC
  REAL(DOUBLE),DIMENSION(3) :: B
#endif

  TYPE(AtomPair)            :: Pair
  TYPE(BSET)                :: BS
  TYPE(CRDS)                :: GM,GM_MM
  TYPE(DBL_RNK4)            :: MD
  TYPE(ARGMT)               :: Args
  TYPE(HGRho)               :: Rho,Rho2
!!  TYPE(CMPoles)             :: MP,PrvMP
  TYPE(INT_VECT)            :: Stat
  INTEGER                   :: P,R,AtA,AtB,NN,iSwitch,IC1,IC2
  INTEGER                   :: NExpt,NDist,NCoef,I,J,Iq,Ir,Pbeg,Pend
  INTEGER                   :: N1,N2,QMOffSetQ,QMOffSetR,PcntDist
  LOGICAL                   :: First
  REAL(DOUBLE)              :: DistThresh,RSumE,RSumN,RSumMM,RelRhoErr, &
                               QMCharge,dQMCharge,MMCharge,dMMCharge,PcntCharge
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg1,Mssg2,RestartHDF
  CHARACTER(LEN=11),PARAMETER :: Prog='ParaMakeRho'
  REAL(DOUBLE),DIMENSION(10,10) :: DD
  REAL(DOUBLE),DIMENSION(3) :: Center
!-------------------------------------------------------------------------------
  ! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)

     IF(SCFActn=='BasisSetSwitch')THEN
        !  Get the previous information
        CALL Get(BS,PrvBase)
        CALL Get(GM,CurGeom)
        CALL SetThresholds(PrvBase)
        CALL Get(NExpt,'nexpt',PrvBase)
        CALL New_HGRho(Rho,(/NExpt,0,0/))
        CALL Get(Rho%Lndx ,'lndex',PrvBase)
        CALL Get(Rho%Expt,'dexpt',PrvBase)
        CALL Get(BSiz,'atsiz',PrvBase)
        CALL Get(OffS,'atoff',PrvBase)
        CALL Get(NBasF,'nbasf',PrvBase)
        CALL Get(Dmat,TrixFile('D',Args,-1))
     ELSEIF(SCFActn=='Restart')THEN
        ! Get the old information
        CALL Get(RestartHDF,'OldInfo')
        CALL CloseHDF()
        CALL OpenHDF(RestartHDF)
        CALL New(Stat,3)
        CALL Get(Stat,'current')
        SCFCycl=TRIM(IntToChar(Stat%I(1)))
        CurBase=TRIM(IntToChar(Stat%I(2)))
        CurGeom=TRIM(IntToChar(Stat%I(3)))
        CALL Get(BS,CurBase)
        CALL Get(GM,CurGeom)
        CALL Get(NExpt,'nexpt',CurBase)
        CALL New_HGRho(Rho,(/NExpt,0,0/))
        CALL Get(Rho%Expt,'dexpt',CurBase)
        CALL Get(Rho%Lndx ,'lndex',CurBase)
        CALL CloseHDF()
        CALL OpenHDF(InfFile)     
        CALL Get(Dmat,TrixFile('D',Args,0))
     ELSE
        ! Get the current information
        CALL Get(BS,CurBase)
        CALL Get(GM,CurGeom)
        CALL Get(NExpt,'nexpt',CurBase)
        CALL New_HGRho(Rho,(/NExpt,0,0/))
        CALL Get(Rho%Expt,'dexpt',CurBase)
        CALL Get(Rho%Lndx ,'lndex',CurBase)
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
     ! Allocations and precalculations
     CALL NewBraBlok(BS)  
     CALL New(MD,(/3,BS%NASym,BS%NASym,2*BS%NASym/),(/1,0,0,0/))

#ifdef PERIODIC
!    Get the Outer Cell Set
     CALL Get_CellSet(CS_OUT,'CS_OUT'//CurBase//CurGeom)
     CALL PPrint(CS_OUT,'outer sum',Prog)
#endif

     !--------------------------------------------------------------
     ! Main loops: First pass calculates the size.
     !             Second pass calculates the density
     !-------------------------------------------------------------
     ! Initailize  NQ
     Rho%NQ%I = 0
     IF(SCFActn/='InkFok')Rho%NQ%I(Rho%NExpt)=NAtoms
     ! Loop over atoms and count primatives
#ifdef PARALLEL
     DO LocalAtom = 1, Dmat%NAtms
        AtA = Beg%I(MyID)+(LocalAtom-1)
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
                    CALL PrimCount(BS,Pair,Rho)
                 ENDIF
              ENDDO
#else
              CALL PrimCount(BS,Pair,Rho)
#endif
           ENDIF
        ENDDO
     ENDDO
     !

! Calculate NDist and NCoef from NQ and Lndx
  NDist = CalNDist(Rho)
  NCoef = CalNCoef(Rho)
  Rho%OffQ%I=CalOffQ(Rho)
  Rho%OffR%I=CalOffR(Rho)
  CALL New_HGRho(Rho,(/NExpt,NDist,NCoef/))
  ! Initailize  RhoCo and First
  First = .TRUE.
  Rho%Co%D=zero
  ! Loop over atoms and calculate the electronic density

#ifdef PARALLEL
     DO LocalAtom = 1, Dmat%NAtms
        AtA = Beg%I(MyID)+(LocalAtom-1)
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
                    CALL RhoBlk(BS,MD,Dmat%MTrix%D(R:R+NN-1),Pair,First,Rho)
                 ENDIF
              ENDDO
#else
              NN = Pair%NA*Pair%NB
              CALL RhoBlk(BS,MD,Dmat%MTrix%D(R:R+NN-1),Pair,First,Rho)
#endif
           ENDIF
        ENDDO
     ENDDO
     ! Add in the density for the nuclear centers
     IF(SCFActn/='InkFok') CALL AddNukes(GM,Rho)

#if defined(PERIODIC) && defined(WRAPDIST)
  CALL Fold_Rho(GM,Rho)
#endif


  ! Prune negligible distributions from the electronic density
  CALL Prune_Rho(Thresholds%Dist,Rho,Rho2) 
  ! Compute integrated electron and nuclear densities

  !! CALL Integrate_HGRho(Rho2,RSumE,RSumN,RSumMM,0,GM%Natms)
  CALL ParaIntegrate_HGRho(Rho2,RSumE,RSumN,RSumMM,0,GM%Natms)
#ifdef PARALLEL
  !! WRITE(*,*) 'MyID = ', MyID, ' RSumE = ', RSumE, ' RSumN = ',RSumN, ', RSumMM = ', RSumMM
  TotRSumE = AllReduce(RSumE)
  TotRSumN = AllReduce(RSumN)
  IF(MyID == 0) THEN
    WRITE(*,*) 'MyID = ', MyID, 'TotRSumE = ', TotRSumE,', TotRSumN = ',TotRSumN
  ENDIF
  RSumE = TotRSumE
  RSumN = TotRSumN
#endif

  QMCharge=Zero
  dQMCharge=Zero
  MMCharge=Zero
  dMMCharge=Zero  

     QMCharge=Two*(RSumE+RSumN)
     dQMCharge=QMCharge+GM%TotCh 
     PcntDist=FLOOR(1.D2*DBLE(Rho2%NDist)/DBLE(Rho%NDist))

     Mssg1=TRIM(Mssg1)//' dNel = '//TRIM(DblToShrtChar(dQMCharge))//', kept '  &
          //TRIM(IntToChar(PcntDist))//'% of distributions.'
     RelRhoErr=ABS(dQMCharge)/DBLE(NEl)

  ! Output pruning and multipole stats
  IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
     !! I=OpenPU()
     !! WRITE(*,*)TRIM(Mssg1)
     !! WRITE(*,*)TRIM(Mssg2)
     !! WRITE(I,*)TRIM(Mssg1)
     !! WRITE(I,*)TRIM(Mssg2)
     !! CALL ClosePU(I)
  ELSEIF(PrintFlags%Key==DEBUG_MEDIUM)THEN
     !! I=OpenPU()
     !! WRITE(I,*)TRIM(Mssg1)
     !! WRITE(I,*)TRIM(Mssg2)
     !! CALL ClosePU(I)
  ENDIF
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
  ! Put Rho and MPs to disk
  IF(SCFActn=='ForceEvaluation')THEN
#ifdef PARALLEL
     CALL Put_HGRho(Rho2,'Rho'//IntToChar(MyID),Args,1)
#else
     CALL Put_HGRho(Rho2,'Rho',Args,1) 
#endif

  ELSEIF(SCFActn=='InkFok')THEN
     CALL Put_HGRho(Rho2,'DeltaRho',Args,0)
  ELSE

#ifdef PARALLEL
     CALL Put_HGRho(Rho2,'Rho'//IntToChar(MyID),Args,0)
#else
     CALL Put_HGRho(Rho2,'Rho',Args,0)
#endif
  ENDIF
  !! CALL PChkSum(Rho2,'Rho',Prog)
! CALL Print_HGRho(Rho2,'Rho2',Unit_O=6)
  ! Tidy up

     CALL Delete(GM)
     CALL Delete(Dmat)
     CALL Delete(BS)
     CALL Delete(MD)
     CALL DeleteBraBlok()

  CALL Delete_HGRho(Rho)
  CALL Delete_HGRho(Rho2)
  CALL ShutDown(Prog)
#endif
END PROGRAM ParaMakeRho
