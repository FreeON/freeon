!    COMPUTE THE DENSITY IN A HGTF BASIS FROM THE DENSITY MATRIX
!    BASED ON AHMADI AND ALMLOF, CPL 246 p.364 (1995) 
!    Authors: Matt Challacombe and C.J. Tymczak
!----------------------------------------------------------------
PROGRAM MakeRho
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
  IMPLICIT NONE
  TYPE(BCSR)                :: Dmat,D1,D2
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
  TYPE(CMPoles)             :: MP,PrvMP
  TYPE(INT_VECT)            :: Stat
  INTEGER                   :: P,R,AtA,AtB,NN,iSwitch,IC1,IC2
  INTEGER                   :: NExpt,NDist,NCoef,I,J,Iq,Ir,Pbeg,Pend
  INTEGER                   :: N1,N2,QMOffSetQ,QMOffSetR,PcntDist
  LOGICAL                   :: First
  REAL(DOUBLE)              :: DistThresh,RSumE,RSumN,RSumMM,RelRhoErr, &
                               QMCharge,dQMCharge,MMCharge,dMMCharge,PcntCharge
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg1,Mssg2,RestartHDF
  CHARACTER(LEN=7),PARAMETER :: Prog='MakeRho'
  REAL(DOUBLE),DIMENSION(10,10) :: DD
  REAL(DOUBLE),DIMENSION(3) :: Center
!---------------------------------------------------------------------------------------
  ! Start up macro
  CALL StartUp(Args,Prog)
#ifdef MMech
  IF(HasMM())THEN
     CALL Get(GM_MM,Tag_O='GM_MM'//CurGeom)
  ENDIF
  IF(HasQM())THEN
#endif
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
     ! Calculate the Number of Cells
     CALL SetCellNumber(GM)
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
     DO AtA=1,NAtoms
        Pbeg = Dmat%RowPt%I(AtA)
        Pend = Dmat%RowPt%I(AtA+1)-1
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
#ifdef MMech
  ELSEIF(MMOnly())THEN
     Nexpt=1
     CALL New_HGRho(Rho,(/Nexpt,0,0/))
     Rho%NQ%I=0
     Rho%Lndx%I=0
     Rho%Expt%D=NuclearExpnt
  ELSE
     CALL Halt(' Logic problem in MakeRho ')
  ENDIF !!!! Mechanics
  IF(HasMM()) THEN
     ! In case of MM, append MM atoms to nuclear-exponent list 
     ! and update counters and offsets
     QMOffSetQ=CalNDist(Rho)
     QMOffSetR=CalNCoef(Rho)
     Rho%NQ%I(Rho%NExpt)=QMOffSetQ+GM_MM%Natms
  ENDIF
#endif
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
#ifdef MMech
  IF(HasQM())THEN
#endif
     DO AtA=1,NAtoms
        Pbeg = Dmat%RowPt%I(AtA)
        Pend = Dmat%RowPt%I(AtA+1)-1
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
#ifdef MMech
  ENDIF
  IF(MMOnly())THEN
     CALL AddNukes(GM_MM,Rho,(/0,0/))
  ELSEIF(HasMM())THEN
     CALL AddNukes(GM_MM,Rho,(/QMOffSetQ,QMOffSetR/))
  ENDIF
#endif
#ifdef PERIODIC 
  !  Fold the Distributions back into the Cell
  CALL Fold_Rho(GM,Rho)
#endif
  ! Prune negligible distributions from the electronic density
  CALL Prune_Rho(Thresholds%Dist,Rho,Rho2) 
  ! Compute integrated electron and nuclear densities
#ifdef MMech
  N1=0
  N2=0
  IF(HasMM())N1=GM_MM%Natms
  IF(HasQM())N2=GM%Natms
  CALL Integrate_HGRho(Rho2,RSumE,RSumN,RSumMM,N1,N2)
#else
  CALL Integrate_HGRho(Rho2,RSumE,RSumN,RSumMM,0,GM%Natms)
#endif
  ! Calculate dipole and quadrupole moments
  CALL New(MP)
  ! Find center of system to expand multipole moments about
  ! This is arbitrary for non-periodic systems...
#ifdef PERIODIC
#ifdef MMech
  IF(HasMM()) THEN
     CALL Warn(' PBCs and QM/MM not yet hacked together...')
     Center(:) = GM_MM%PBC%CellCenter(:)
  ELSE
#endif
     Center(:) = GM%PBC%CellCenter(:)
#ifdef MMech
  ENDIF
#endif
#else
#ifdef MMech
  IF(HasMM()) THEN
     Center(:) = Half*(GM_MM%BndBox%D(:,2)+GM_MM%BndBox%D(:,1))
  ELSE
#endif
     Center(:) = Half*(GM%BndBox%D(:,2)+GM%BndBox%D(:,1))
#ifdef MMech
  ENDIF
#endif
#endif
  CALL CalRhoPoles(MP,Center,Rho2)
  ! Format output for pruning and multipole stats
#ifdef MMech
  IF(HasQM()) THEN
#endif
     IF(SCFActn=='InkFok')THEN
        Mssg1=ProcessName(Prog,'InkFok')
        Mssg2=Mssg1
     ELSE
        Mssg1=ProcessName(Prog,'Pruned Rho')
        Mssg2=ProcessName(Prog,'Moments')
     ENDIF
#ifdef MMech
  ENDIF
#endif
  QMCharge=Zero
  dQMCharge=Zero
  MMCharge=Zero
  dMMCharge=Zero  
#ifdef MMech
  IF(HasQM())THEN
#endif
     QMCharge=Two*(RSumE-RSumN)
     IF(ABS(QMCharge-GM%TotCh)>1.D-2)THEN
        CALL Halt(' Wrong charge state in MakeRho: '//Rtrn &
	          //' Integrated electron population = '   &
                  //TRIM(DblToShrtChar(RSumE))//Rtrn       &
	          //', Integrated nuclear population = '   &
                  //TRIM(DblToShrtChar(RSumN))//Rtrn       &
                  //TRIM(DblToShrtChar(GM%TotCh)) )
     ENDIF
     dQMCharge=Two*RSumE-NEl
     PcntDist=FLOOR(1.D2*DBLE(Rho2%NDist)/DBLE(Rho%NDist))
#ifdef MMech
  ENDIF
  IF(HasMM())THEN
     MMCharge=Two*RSumMM
     dMMCharge=MMCharge+GM_MM%TotCh
  ENDIF
  IF(MMOnly())THEN
     Mssg1=TRIM(Mssg1)//' dMMChg = ' //TRIM(DblToShrtChar(dMMCharge))
     RelRhoErr=ABS(MMCharge)/DBLE(GM_MM%Natms)
  ELSEIF(QMOnly())THEN
#endif
     Mssg1=TRIM(Mssg1)//' dNel = '//TRIM(DblToShrtChar(dQMCharge))//', kept '  &
          //TRIM(IntToChar(PcntDist))//'% of distributions.'
     RelRhoErr=ABS(dQMCharge)/DBLE(NEl)
#ifdef MMech
  ELSE
     Mssg1=TRIM(Mssg1)//' QMChg = '//TRIM(DblToShrtChar(QMCharge))//', '       &
          //'MMChg = '//TRIM(DblToShrtChar(MMCharge))//', '       &
          //'dNel = '//TRIM(DblToShrtChar(dQMCharge))//', kept '  &
          //TRIM(IntToChar(PcntDist))//'% of distributions.'
     RelRhoErr=ABS(QMCharge+MMCharge)/DBLE(NEl) 
  ENDIF
#endif
  Mssg2=TRIM(Mssg2)                                      &
        //' <r> = ('//TRIM(DblToShrtChar(MP%DPole%D(1))) &
        //', '//TRIM(DblToShrtChar(MP%DPole%D(2)))       &
        //', '//TRIM(DblToShrtChar(MP%DPole%D(3)))       &
        //'), <r^2> = '//TRIM(DblToShrtChar(             &
           MP%QPole%D(1)+MP%QPole%D(2)+MP%QPole%D(3)))
  ! Output pruning and multipole stats
  IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
     I=OpenPU()
     WRITE(*,*)TRIM(Mssg1)
     WRITE(*,*)TRIM(Mssg2)
     WRITE(I,*)TRIM(Mssg1)
     WRITE(I,*)TRIM(Mssg2)
     CALL ClosePU(I)
  ELSEIF(PrintFlags%Key==DEBUG_MEDIUM)THEN
     I=OpenPU()
     WRITE(I,*)TRIM(Mssg1)
     WRITE(I,*)TRIM(Mssg2)
     CALL ClosePU(I)
  ENDIF
  ! Check error
  IF(RelRhoErr>Thresholds%Dist*5.D3.AND.SCFActn/='NumForceEvaluation') &
       CALL Halt('In MakeRho, relative error in density = '//TRIM(DblToShrtChar(RelRhoErr)) &
       //'. Distribution threshold = '//TRIM(DblToShrtChar(Thresholds%Dist))      &
       //'. Total charge lost = '//TRIM(DblToShrtChar(dQMCharge+dMMCharge)))
  ! Put Rho and MPs to disk
  IF(SCFActn=='ForceEvaluation')THEN
     CALL Put_HGRho(Rho2,'Rho',Args,1) 
     CALL Put(MP,NxtCycl)
  ELSEIF(SCFActn=='InkFok')THEN
     CALL Put_HGRho(Rho2,'DeltaRho',Args,0)
     CALL Put(MP,'Delta'//TRIM(SCFCycl))
  ELSE
     CALL Put_HGRho(Rho2,'Rho',Args,0) 
     CALL Put(MP,IntToChar(Current(1))) 
  ENDIF
  CALL PChkSum(Rho2,'Rho',Prog)
  ! Tidy up
#ifdef MMech
  IF(HasMM()) THEN
     CALL Delete(GM_MM)
  ENDIF
  IF(HasQM()) THEN
#endif
     CALL Delete(GM)
     CALL Delete(Dmat)
     CALL Delete(BS)
     CALL Delete(MD)
     CALL DeleteBraBlok()
#ifdef MMech
  ENDIF
#endif
  CALL Delete_HGRho(Rho)
  CALL Delete_HGRho(Rho2)
  CALL ShutDown(Prog)
END PROGRAM MakeRho

