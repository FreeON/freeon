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
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef PARALLEL
  TYPE(DBCSR)               :: Dmat
#else
  TYPE(BCSR)                :: Dmat
#endif
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
  TYPE(CMPoles)             :: MP
  INTEGER                   :: P,R,AtA,AtB,NN,iSwitch,IC1,IC2
  INTEGER                   :: NExpt,NDist,NCoef,I,J,Iq,Ir,Pbeg,Pend
  LOGICAL                   :: First
  REAL(DOUBLE)              :: DistThresh,RSumE,RSumN,RelRhoErr
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
  CHARACTER(LEN=7),PARAMETER :: Prog='MakeRho'
!----------------------------------------------
! Start up macro
!----------------------------------------------
  CALL StartUp(Args,Prog)
!----------------------------------------------
  IF(SCFActn=='BasisSetSwitch')THEN
!    Get the previous information
     CALL Get(BS,PrvBase)
     CALL Get(GM,CurGeom)
     CALL Get(NExpt,'nexpt',PrvBase)
     CALL New_HGRho(Rho,(/NExpt,0,0/))
     CALL Get(Rho%Lndx ,'lndex',PrvBase)
     CALL Get(Rho%Expt,'dexpt',PrvBase)
     CALL Get(BSiz,'atsiz',PrvBase)
     CALL Get(OffS,'atoff',PrvBase)
     CALL Get(NBasF,'nbasf',PrvBase)
     CALL Get(Dmat,TrixFile('D',Args,-1))
  ELSE
!    Get the current information
     CALL Get(BS,CurBase)
     CALL Get(GM,CurGeom)
     CALL Get(NExpt,'nexpt',CurBase)
     CALL New_HGRho(Rho,(/NExpt,0,0/))
     CALL Get(Rho%Expt,'dexpt',CurBase)
     CALL Get(Rho%Lndx ,'lndex',CurBase)
     IF(SCFActn=='InkFok')THEN
        CALL Get(Dmat,TrixFile('DeltaD',Args,0))
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
  IF(SCFActn=='Core')THEN
!    Re-allocate the density
     CALL New_HGRho(Rho,(/NExpt,NAtoms,NAtoms/))
!    Initailize  NQ
     Rho%NQ%I            = 0
     Rho%NQ%I(Rho%NExpt) = NAtoms
!    Initailize  OffQ,OffR and RhoCo
     Rho%OffQ%I=CalOffQ(Rho)
     Rho%OffR%I=CalOffR(Rho)
     Rho%Co%D=Zero
!    Add in the density for the nuclear centers
     CALL AddNukes(GM,Rho)
  ELSE
!    Initailize  NQ
     Rho%NQ%I = 0
     IF(SCFActn/='InkFok')Rho%NQ%I(Rho%NExpt)=NAtoms
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
                    CALL PrimCount(BS,Pair,Rho)
                 ENDIF
              ENDDO
#else
              CALL PrimCount(BS,Pair,Rho)
#endif
           ENDIF
        ENDDO
     ENDDO
!    Calculate NDist and NCoef from NQ and Lndx
     NDist = CalNDist(Rho)
     NCoef = CalNCoef(Rho)
!    Initailize  OffQ,OffR and RhoCo
     Rho%OffQ%I=CalOffQ(Rho)
     Rho%OffR%I=CalOffR(Rho)
!    Re-allocate the density
     CALL New_HGRho(Rho,(/NExpt,NDist,NCoef/))
!    Initailize  RhoCo and First
     First = .TRUE.
     Rho%Co%D=zero
!-----------------------------------------------------
!    Loop over atoms and calculate the electronic 
!    density
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
!---------------------------------------------------
!    Add in the density for the nuclear centers
!
     IF(SCFActn/='InkFok')CALL AddNukes(GM,Rho)
  ENDIF
#ifdef PERIODIC
!-----------------------------------------------------------
!  Fold the Distributions back into the Cell
!
!  CALL Fold_Rho(GM,Rho)
#endif
!------------------------------------------------------------
!  Remove distribution which do not contibute significantly 
!  to the density
!
  CALL Prune_Rho(Thresholds%Dist,Rho,Rho2) 
!------------------------------------------------------------
! Halt if we lost to many electrons
!
  CALL Integrate_HGRho(Rho2,RSumE,RSumN)
  RelRhoErr=Two*ABS(RSumE+RSumN)/DBLE(NEl)
  Mssg=ProcessName(Prog,'Pruned Rho')//'dNel = '    &
       //TRIM(DblToShrtChar(Two*ABS(RSumE+RSumN)))   &
       //', '//TRIM(IntToChar(FLOOR(1.D2*DBLE(Rho2%NDist)/DBLE(Rho%NDist)))) &
       //'% of distributions retained'
  IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
     I=OpenPU()
     WRITE(*,*)TRIM(Mssg)
     WRITE(I,*)TRIM(Mssg)
     CALL ClosePU(I)
  ELSEIF(PrintFlags%Key==DEBUG_MEDIUM)THEN
     I=OpenPU()
     WRITE(I,*)TRIM(Mssg)
     CALL ClosePU(I)
  ENDIF
  IF(RelRhoErr>Thresholds%Dist*1.D3) &
  CALL Halt('In MakeRho, missing '//TRIM(DblToShrtChar(Two*ABS(RSumE+RSumN)))   &
           //' electrons after pruning.')
!------------------------------------------------------------
!  Calculate the DiPole and QuadruPole moments of Rho
!
  CALL CalRhoPoles(MP,GM,RHo2)
  IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
     CALL PPrint(MP,Prog,Unit_O=6)  
     CALL PPrint(MP,Prog)  
  ELSEIF(PrintFlags%Key==DEBUG_MEDIUM)THEN
     CALL PPrint(MP,Prog)  
  ENDIF
!------------------------------------------------------------
! Put Rho and MP to disk
! 
  IF(SCFActn=='ForceEvaluation')THEN
     CALL Put_HGRho(Rho2,'Rho',Args,1) 
     CALL Put(MP,NxtCycl)
 ELSE 
     CALL Put_HGRho(Rho2,'Rho',Args,0)
     CALL Put(MP,SCFCycl)
  ENDIF
!------------------------------------------------------------
! Printing
!
  CALL PChkSum(Rho2,'Rho',Prog)
!
!  PrintFlags%Fmt=DEBUG_MMASTYLE
!  CALL PPrint(Rho,'Rho',Unit_O=6)
!  CALL PPrint(Rho2,'Rho2',Unit_O=6)
!---------------------------------------------------
! Tidy up
! 
  IF(SCFActn/='Core') &
  CALL Delete(Dmat)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(MD)
  CALL DeleteBraBlok()
  CALL Delete_HGRho(Rho)
  CALL Delete_HGRho(Rho2)
  CALL ShutDown(Prog)
!
END PROGRAM MakeRho

