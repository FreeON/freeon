!    COMPUTE THE DENSITY IN A HGTF BASIS FROM THE DENSITY MATRIX
!    BASED ON AHMADI AND ALMLOF, CPL 246 p.364 (1995) 
!    Authors: Matt Challacombe and C.J. Tymczak
!-------------------------------------------------------------------
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
!#ifdef PARALLEL
!  USE MondoMPI
!#endif
  IMPLICIT NONE
!#ifdef PARALLEL
!  TYPE(DBCSR)               :: Dmat,D1,D2
!#else
  TYPE(BCSR)                :: Dmat,D1,D2
!#endif
#ifdef PERIODIC 
  INTEGER                   :: NC
  REAL(DOUBLE),DIMENSION(3) :: B
#endif
  TYPE(AtomPair)            :: Pair
  TYPE(BSET)                :: BS
  TYPE(CRDS)                :: GM
  TYPE(DBL_RNK4)            :: MD
  TYPE(ARGMT)               :: Args
  TYPE(HGRho)               :: Rho,Rho2
  TYPE(CMPoles)             :: MP,PrvMP
  TYPE(INT_VECT)            :: Stat
  INTEGER                   :: P,R,AtA,AtB,NN,iSwitch,IC1,IC2
  INTEGER                   :: NExpt,NDist,NCoef,I,J,Iq,Ir,Pbeg,Pend
  LOGICAL                   :: First
  REAL(DOUBLE)              :: DistThresh,RSumE,RSumN,RelRhoErr
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg1,Mssg2,RestartHDF
  CHARACTER(LEN=7),PARAMETER :: Prog='MakeRho'
  REAL(DOUBLE),DIMENSION(10,10) :: DD
!----------------------------------------------
! Start up macro
!----------------------------------------------
  CALL StartUp(Args,Prog)
!----------------------------------------------
  IF(SCFActn=='BasisSetSwitch')THEN
!    Get the previous information
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
!    Get the old information
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
!     WRITE(*,*)TRIM(TrixFile('D',Args,0))
!     CALL PChkSum(DMat,'RestartDM',Prog,Unit_O=6)
  ELSE
!    Get the current information
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
!        CALL PPrint(DMat,'D',Unit_O=6)
     ELSEIF(SCFActn=='ForceEvaluation')THEN
        CALL Get(Dmat,TrixFile('D',Args,1))
     ELSEIF(SCFActn/='Core')THEN
        CALL Get(Dmat,TrixFile('D',Args,0))
     ENDIF
  ENDIF
! Allocations and precalculations
  CALL NewBraBlok(BS)  
  CALL New(MD,(/3,BS%NASym,BS%NASym,2*BS%NASym/),(/1,0,0,0/))
  CALL New(MP)
#ifdef PERIODIC
! Calculate the Number of Cells
  CALL SetCellNumber(GM)
  CALL PPrint(CS_OUT,'outer sum',Prog)
#endif
!---------------------------------------------------
! Main loops: First pass calculates the size.
!             Second pass calculates the density
!---------------------------------------------------
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
! Calculate NDist and NCoef from NQ and Lndx
  NDist = CalNDist(Rho)
  NCoef = CalNCoef(Rho)
! Initailize  OffQ,OffR and RhoCo
  Rho%OffQ%I=CalOffQ(Rho)
  Rho%OffR%I=CalOffR(Rho)
! Re-allocate the density
  CALL New_HGRho(Rho,(/NExpt,NDist,NCoef/))
! Initailize  RhoCo and First
  First = .TRUE.
  Rho%Co%D=zero
!-----------------------------------------------------
! Loop over atoms and calculate the electronic 
! density
!
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
!--------------------------------------------------------
! Add in the density for the nuclear centers
!
  IF(SCFActn/='InkFok')CALL AddNukes(GM,Rho)
!--------------------------------------------------------
! Prune negligible distributions from density
  CALL Prune_Rho(Thresholds%Dist,Rho,Rho2) 
! Compute integrated electron and nuclear densities
  CALL Integrate_HGRho(Rho2,RSumE,RSumN)
! Calculate dipole and quadrupole moments
  CALL CalRhoPoles(MP,GM,Rho2)
! Format output for pruning and multipole stats
  IF(SCFActn=='InkFok')THEN
     Mssg1=ProcessName(Prog,'InkFok')
     Mssg2=Mssg1
  ELSE
     Mssg1=ProcessName(Prog,'Pruned Rho')
     Mssg2=ProcessName(Prog,'Moments')
  ENDIF
  RelRhoErr=Two*ABS(RSumE+RSumN)/DBLE(NEl)
  Mssg1=TRIM(Mssg1)//'dNel = '                       &
       //TRIM(DblToShrtChar(Two*ABS(RSumE+RSumN)))   &
       //', '//TRIM(IntToChar(FLOOR(1.D2*DBLE(Rho2%NDist)/DBLE(Rho%NDist)))) &
       //'% of distributions retained'
  Mssg2=TRIM(Mssg2)                                     &
        //'<r> = ('//TRIM(DblToShrtChar(MP%DPole%D(1))) &
        //', '//TRIM(DblToShrtChar(MP%DPole%D(2)))      &
        //', '//TRIM(DblToShrtChar(MP%DPole%D(3)))      &
        //'), <r^2> = '//TRIM(DblToShrtChar(            &
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
! Check for error
  IF(RelRhoErr>Thresholds%Dist*1.D4.AND.SCFActn/='NumForceEvaluation') &
       CALL Halt('In MakeRho, missing '//TRIM(DblToShrtChar(Two*ABS(RSumE+RSumN)))   &
       //' electrons after pruning.')
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
  CALL Delete(Dmat)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(MD)
  CALL DeleteBraBlok()
  CALL Delete_HGRho(Rho)
  CALL Delete_HGRho(Rho2)
  CALL ShutDown(Prog)
END PROGRAM MakeRho

