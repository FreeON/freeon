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
  TYPE(BCSR)                      :: Dmat,D1,D2,S
#ifdef PERIODIC 
  INTEGER                         :: NC
  REAL(DOUBLE),DIMENSION(3)       :: B
#endif
  TYPE(AtomPair)                  :: Pair
  TYPE(BSET)                      :: BS
  TYPE(CRDS)                      :: GM,GM_MM
  TYPE(ARGMT)                     :: Args
  TYPE(HGRho)                     :: Rho
  TYPE(HGRho_new)                 :: RhoA
  TYPE(CMPoles)                   :: MP,PrvMP
  TYPE(INT_VECT)                  :: Stat
  INTEGER                         :: P,R,AtA,AtB,NN,iSwitch,IC1,IC2
  INTEGER                         :: NExpt,NDist,NCoef,I,J,K,Iq,Ir,Pbeg,Pend,NDist_old,NDist_new
  INTEGER                         :: N1,N2,QMOffSetQ,QMOffSetR,PcntDist,OldFileID
  REAL(DOUBLE)                    :: DistThresh,RSumE,RSumN,RSumMM,RSum_TPS,RelRhoErr, &
                                     QMCharge,dQMCharge,MMCharge,dMMCharge,PcntCharge
  CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Mssg1,Mssg2,RestartHDF
  CHARACTER(LEN=7),PARAMETER      :: Prog='MakeRho'
!---------------------------------------------------------------------------------------
! Start up macro
  CALL StartUp(Args,Prog)
#ifdef PERIODIC
#ifdef PARALLEL_CLONES
#else
!    Get the Outer Cell Set
     CALL Get(CS_OUT,'CS_OUT',Tag_O=CurBase)
#endif
#endif
#ifdef MMech
!---------------------------------------------------------------
  IF(HasMM())THEN
     CALL Get(GM_MM,Tag_O='GM_MM'//CurGeom)
  ENDIF
!---------------------------------------------------------------
  IF(HasQM())THEN
#endif
     IF(SCFActn=='BasisSetSwitch')THEN
!       Get the previous information
        CALL Get(BS,PrvBase)
        CALL Get(GM,CurGeom)
        CALL SetThresholds(PrvBase)
        CALL Get(BSiz,'atsiz',PrvBase)
        CALL Get(OffS,'atoff',PrvBase)
        CALL Get(NBasF,'nbasf',PrvBase)
        CALL Get(Dmat,TrixFile('D',Args,-1))
     ELSEIF(SCFActn=='Restart')THEN
#ifdef PARALLEL_CLONES
        ! Get the current geometry from the current HDF first
        CALL Get(GM,CurGeom)
        ! ... then close current group and HDF
        CALL CloseHDFGroup(H5GroupID)
        CALL CloseHDF(HDFFileID)
        ! ... and open the old group and HDF
        HDF_CurrentID=OpenHDF(Restart)
        OldFileID=HDF_CurrentID
        CALL New(Stat,3)
        CALL Get(Stat,'current_state')
        HDF_CurrentID=OpenHDFGroup(HDF_CurrentID,"Clone #"//TRIM(IntToChar(MyClone)))
!       Get old basis set stuff
        SCFCycl=TRIM(IntToChar(Stat%I(1)))
        CurBase=TRIM(IntToChar(Stat%I(2)))
        CurGeom=TRIM(IntToChar(Stat%I(3)))
        CALL Get(BS,CurBase)
        ! Compute a sparse matrix blocking scheme for the old BS
        CALL BlockBuild(GM,BS,BSiz,OffS)
!       Close the old hdf up 
        CALL CloseHDFGroup(HDF_CurrentID)
        CALL CloseHDF(OldFileID)
!       Reopen current group and HDF
        HDFFileID=OpenHDF(H5File)
        H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
        HDF_CurrentID=H5GroupID
#else
!       Get the old information
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
!       Get the current information
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
! Allocations and precalculations
!***********************************************************************************
!     CALL Get(S,TrixFile('S',Args))
!     RSum_TPS = Trace(Dmat,S)
!     WRITE(*,*) 'Trace[P_AO*S]  = ', RSum_TPS
!***********************************************************************************
     CALL NewBraBlok(BS)  
!--------------------------------------------------------------
! Main loops: First pass calculates the size.
!             Second pass calculates the density
!-------------------------------------------------------------
!    Initailize      
     NDist = 0
     NCoef = 0
!    Loop over atoms and count primatives
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
                    CALL PrimCount(BS,Pair,NDist,NCoef)
                 ENDIF
              ENDDO
#else
              CALL PrimCount(BS,Pair,NDist,NCoef)
#endif
           ENDIF
        ENDDO
     ENDDO
!    Allocate the Density
     CALL New_HGRho_new(RhoA,(/NDist,NCoef/))
!    Initailize  Counters
     NDist        = 0
     NCoef        = 0
     DO AtA=1,NAtoms
        Pbeg = Dmat%RowPt%I(AtA)
        Pend = Dmat%RowPt%I(AtA+1)-1
        DO P=Pbeg,Pend
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
     IF(SCFActn/='InkFok') THEN
        CALL AddDist(RhoA,GM,NuclearExpnt,1,GM%NAtms)
     ENDIF
#ifdef MMech
  ELSE
     CALL New_HGRho_new(RhoA,(/0,0/))
  ENDIF
! In case of MM, append MM atoms to nuclear-exponent list 
  IF(HasMM()) THEN
     CALL AddDist(RhoA,GM_MM,NuclearExpnt,1,GM_MM%NAtms)
  ENDIF
#endif
!***********************************************************************************
  RSumE  =  Integrate_HGRho_new(RhoA,1,RhoA%NDist-GM%NAtms)
!  WRITE(*,*) 'Integrate[Rho] = ',RSumE
!***********************************************************************************
! Prune negligible distributions from the electronic density
  NDist_old = RhoA%NDist
  CALL Prune_Rho_new(Thresholds%Dist,RhoA)
  NDist_new = RhoA%NDist
#ifdef PERIODIC
! Fold distributions back into the box
  CALL Fold_Rho_new(GM,RhoA)
#endif
#ifdef MMech
! Compute integrated electron and nuclear densities
  IF(HasQM()) THEN
     IF(HasMM()) THEN
        RSumE  =  Integrate_HGRho_new(RhoA,1                                ,RhoA%NDist-GM%NAtms-GM_MM%NAtms)
        RSumN  =  Integrate_HGRho_new(RhoA,RhoA%NDist-GM%NAtms-GM_MM%NAtms+1,RhoA%NDist-GM_MM%NAtms         )
     ELSE
        RSumE  =  Integrate_HGRho_new(RhoA,1                    ,RhoA%NDist-GM%NAtms)
        RSumN  =  Integrate_HGRho_new(RhoA,RhoA%NDist-GM%NAtms+1,RhoA%NDist         )
     ENDIF
  ENDIF
  IF(HasMM()) THEN
     RSumMM =  Integrate_HGRho_new(RhoA,RhoA%NDist-GM_MM%NAtms+1,RhoA%NDist)
  ENDIF
! Calculate dipole and quadrupole moments
  CALL New(MP)
  IF(HasMM()) THEN
     CALL CalRhoPoles_new(MP,RhoA,GM_MM)
  ELSE
     CALL CalRhoPoles_new(MP,RhoA,GM)
  ENDIF
#else
! Compute integrated electron and nuclear densities
  RSumE  =  Integrate_HGRho_new(RhoA,1                    ,RhoA%NDist-GM%NAtms)
  RSumN  =  Integrate_HGRho_new(RhoA,RhoA%NDist-GM%NAtms+1,RhoA%NDist         )
! Calculate dipole and quadrupole moments
  CALL New(MP)
  CALL CalRhoPoles_new(MP,RhoA,GM)
#endif
!
!  Convert to the old format
!
  CALL ConvertToOldRho(Rho,RhoA)
! What follow is a complete Mess: I did not touch this - CJ
! Format output for pruning and multipole stats
#ifdef MMech
  IF(HasQM()) THEN
     IF(SCFActn=='InkFok')THEN
        Mssg1=ProcessName(Prog,'InkFok')
        Mssg2=Mssg1
     ELSE
        Mssg1=ProcessName(Prog,'Pruned Rho')
        Mssg2=ProcessName(Prog,'Moments')
     ENDIF
  ELSE
        Mssg1=ProcessName(Prog,'Pruned Rho')
        Mssg2=ProcessName(Prog,'Moments')
  ENDIF
#else
     IF(SCFActn=='InkFok')THEN
        Mssg1=ProcessName(Prog,'InkFok')
        Mssg2=Mssg1
     ELSE
        Mssg1=ProcessName(Prog,'Pruned Rho')
        Mssg2=ProcessName(Prog,'Moments')
     ENDIF
#endif
  QMCharge=Zero
  dQMCharge=Zero
  MMCharge=Zero
  dMMCharge=Zero  
#ifdef MMech
  IF(HasQM())THEN
#endif
     QMCharge=Two*(RSumE+RSumN)
     dQMCharge=QMCharge+GM%TotCh 
     PcntDist=FLOOR(1.D2*DBLE(NDist_new)/DBLE(NDist_old))
#ifdef MMech
  ENDIF
  IF(HasMM())THEN
     MMCharge=Two*RSumMM
     dMMCharge=MMCharge+GM_MM%TotCh
  ENDIF
  IF(MMOnly())THEN
     Mssg1=TRIM(Mssg1)//' dMMChg = ' //TRIM(DblToShrtChar(dMMCharge))
     RelRhoErr=ABS(dMMCharge)/DBLE(GM_MM%Natms)
  ELSEIF(QMOnly())THEN
#endif
     Mssg1=TRIM(Mssg1)//' dNel = '//TRIM(DblToShrtChar(dQMCharge))//', kept '  &
          //TRIM(IntToChar(PcntDist))//'% of distributions.'
     RelRhoErr=ABS(dQMCharge)/DBLE(NEl)
#ifdef MMech
  ELSE
     Mssg1=TRIM(Mssg1)//' QMChg = '//TRIM(DblToShrtChar(QMCharge))//', '   &
          //'MMChg = '//TRIM(DblToShrtChar(MMCharge))//', '                &
          //'dNel = '//TRIM(DblToShrtChar(dQMCharge))//','//Rtrn           &
          //'                         kept '//TRIM(IntToChar(PcntDist))//'% of distributions.'
     RelRhoErr=ABS(dQMCharge+dMMCharge)/DBLE(NEl) 
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
  IF(RelRhoErr>Thresholds%Dist*5.D3.AND.SCFActn/='NumForceEvaluation')THEN
        CALL Warn(ProcessName(Prog)//'relative error in density = '//TRIM(DblToShrtChar(RelRhoErr)) &
             //'. Distribution threshold = '//TRIM(DblToShrtChar(Thresholds%Dist))      &
             //'. Total charge lost = '//TRIM(DblToShrtChar(dQMCharge+dMMCharge)))
  ENDIF
!------------------------------------------------------------------------------------
! Put Rho and MPs to disk
#ifdef MMech
  IF(SCFActn=='ForceEvaluation')THEN
     IF(MMOnly()) THEN
        CALL Put_HGRho(Rho,'Rho',Args,Current(1)) 
#ifdef PARALLEL_CLONES
       CALL Put(MP)  
#else
       CALL Put(MP,CurGeom)
#endif
     ELSE
        CALL Put_HGRho(Rho,'Rho',Args,1) 
#ifdef PARALLEL_CLONES
        CALL Put(MP)   
#else
        CALL Put(MP,NxtCycl)
#endif
     ENDIF
  ELSEIF(SCFActn=='InkFok')THEN
     CALL Put_HGRho(Rho,'DeltaRho',Args,0)
     CALL Put(MP,'Delta'//TRIM(SCFCycl))
  ELSE
     IF(MMOnly()) THEN
        CALL Put_HGRho(Rho,'Rho',Args,0)
#ifdef PARALLEL_CLONES
        CALL Put(MP)   
#else
        CALL Put(MP,CurGeom) 
#endif
     ELSE
        CALL Put_HGRho(Rho,'Rho',Args,0) 
#ifdef PARALLEL_CLONES
        CALL Put(MP)   
#else
        CALL Put(MP,IntToChar(Current(1))) 
#endif
     ENDIF
  ENDIF
#else
  IF(SCFActn=='ForceEvaluation')THEN
     CALL Put_HGRho(Rho,'Rho',Args,1)
#ifdef PARALLEL_CLONES
       CALL Put(MP)  
#else
     CALL Put(MP,NxtCycle)
#endif
  ELSEIF(SCFActn=='InkFok')THEN
     CALL Put_HGRho(Rho,'DeltaRho',Args,0)
     CALL Put(MP,'Delta'//TRIM(SCFCycl))
  ELSE
     CALL Put_HGRho(Rho,'Rho',Args,0) 
#ifdef PARALLEL_CLONES
       CALL Put(MP)  
#else
     CALL Put(MP,IntToChar(Current(1))) 
#endif
  ENDIF
#endif
  CALL PChkSum(Rho,'Rho',Prog)
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
     CALL DeleteBraBlok()
#ifdef MMech
  ENDIF
#endif
  CALL Delete_HGRho(Rho)
  CALL Delete_HGRho_new(RhoA)
  CALL ShutDown(Prog)
END PROGRAM MakeRho

