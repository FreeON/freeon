!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe and  C. J. Tymczak
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!    COMPUTE THE DENSITY FROM THE DENSITY MATRIX
!
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
  USE RhoBlok
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef PARALLEL
  TYPE(DBCSR)         :: Dmat
#else
  TYPE(BCSR)          :: Dmat
#endif
#ifdef PERIODIC 
  INTEGER             :: NC
  REAL(DOUBLE)        :: Bx,By,Bz
#endif
  TYPE(AtomPair)      :: Pair
  TYPE(BSET)          :: BS
  TYPE(CRDS)          :: GM
  TYPE(DBL_RNK4)      :: MD
  TYPE(ARGMT)         :: Args
  TYPE(HGRho)         :: Rho
  INTEGER             :: P,R,AtA,AtB,NN,iSwitch
  INTEGER             :: NExpt,NDist,NCoef,I,J,Iq,Ir,Pbeg,Pend
  LOGICAL             :: First
!
  CHARACTER(LEN=7),PARAMETER :: Prog='MakeRho'
!----------------------------------------------
! Start up macro
!----------------------------------------------
  CALL StartUp(Args,Prog)
!----------------------------------------------
! Get basis set and geometry
!----------------------------------------------
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
  CALL NewBraKetBlok(BS)
!----------------------------------------------
! Set up the appropriate action
!----------------------------------------------
  IF(Args%C%C(2)=='Core') THEN
     iSwitch=0
  ELSEIF(Args%C%C(2)=='Direct') THEN
     iSwitch=1
  ELSEIF(Args%C%C(2)=='Switch') THEN
     iSwitch=1
  ELSEIF(Args%C%C(2)=='InkFok') THEN
     iSwitch=2
  ELSE
     CALL MondoHalt(-100,' Inappropriate Action in MakeRho:'//Args%C%C(2))
  ENDIF
!----------------------------------------------
! Get the Density Matrix
!----------------------------------------------
  IF(iSwitch==1) THEN
     CALL Get(Dmat,TrixFile('D',Args,0))
  ELSEIF(iSwitch==2) THEN
     CALL Get(Dmat,TrixFile('DeltaD',Args,0))
  ENDIF
!---------------------------------------------- 
! Allocations 
!
  CALL New(MD,(/3,BS%NASym,BS%NASym,2*BS%NASym/),(/1,0,0,0/))
#ifdef PERIODIC
!-----------------------------------------------
! Calculate the Number of Cells
!-----------------------------------------------
  CALL SetCellNumber(GM)
#endif
!---------------------------------------------------
! Get the Exponents and Angular Symmetries
!
  CALL Get(NExpt,'nexpt',Tag_O=CurBase)
  CALL New_HGRho(Rho,(/NExpt,0,0/),.TRUE.)
  CALL Get(Rho%Expt,'dexpt',Tag_O=CurBase)
  CALL Get(Rho%Lndx ,'lndex',Tag_O=CurBase)
!---------------------------------------------------
! Main loops: First pass calculates the size.
!             Second pass calculates the density
!---------------------------------------------------
  IF(iSwitch==0) THEN
!---------------------------------------------------
!    Re-allocate the density
!
     CALL New_HGRho(Rho,(/NExpt,NAtoms,NAtoms/),.TRUE.)
!-----------------------------------------------------
!    Initailize  NQ
!
     Rho%NQ%I            = 0
     Rho%NQ%I(Rho%NExpt) = NAtoms
!-----------------------------------------------------
!    Initailize  OffQ,OffR and RhoCo
!
     Rho%OffQ%I=CalOffQ(Rho)
     Rho%OffR%I=CalOffR(Rho)
     Rho%Co%D=Zero
!---------------------------------------------------
!    Add in the density for the nuclear centers
!
     CALL AddNukes(GM,Rho)
!---------------------------------------------------
!    Calculate Integral estimates for Rho
!
     CALL RhoEst(Rho)
  ELSE
!-----------------------------------------------------
!    Initailize  NQ
!
     Rho%NQ%I            = 0
     IF(iSwitch == 1) Rho%NQ%I(Rho%NExpt) = NAtoms
!----------------------------------------------------
!    Loop over atoms and count primatives
!
     DO AtA=1,NAtoms
        Pbeg = Dmat%RowPt%I(AtA)
        Pend = Dmat%RowPt%I(AtA+1)-1
        DO P = Pbeg,Pend
           AtB = Dmat%ColPt%I(P)
           IF(AtB >= AtA) THEN
              IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
#ifdef PERIODIC                 
                 Bx = Pair%B(1)
                 By = Pair%B(2)           
                 Bz = Pair%B(3)
                 DO NC = 1,CS%NCells
                    Pair%B(1) = Bx+CS%CellCarts%D(1,NC)
                    Pair%B(2) = By+CS%CellCarts%D(2,NC)
                    Pair%B(3) = Bz+CS%CellCarts%D(3,NC)
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
           ENDIF
        ENDDO
     ENDDO
!----------------------------------------------------
!    Calculate NDist and NCoef from NQ and Lndx
!
     NDist = CalNDist(Rho)
     NCoef = CalNCoef(Rho)
!-----------------------------------------------------
!    Initailize  OffQ,OffR and RhoCo
!
     Rho%OffQ%I=CalOffQ(Rho)
     Rho%OffR%I=CalOffR(Rho)
!-----------------------------------------------------
!    Re-allocate the density
!
     CALL New_HGRho(Rho,(/NExpt,NDist,NCoef/),.TRUE.)
!-----------------------------------------------------
!    Initailize  RhoCo and First
!
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
           IF(AtB >= AtA) THEN
              IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
#ifdef PERIODIC                 
                 Bx = Pair%B(1)
                 By = Pair%B(2)           
                 Bz = Pair%B(3)
                 DO NC = 1,CS%NCells
                    Pair%B(1) = Bx+CS%CellCarts%D(1,NC)
                    Pair%B(2) = By+CS%CellCarts%D(2,NC)
                    Pair%B(3) = Bz+CS%CellCarts%D(3,NC)
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
           ENDIF
        ENDDO
     ENDDO 
!---------------------------------------------------
!    Add in the density for the nuclear centers
!
     IF(iSwitch == 1) CALL AddNukes(GM,Rho)
!---------------------------------------------------
!    Calculate Integral estimates for Rho
!
     CALL RhoEst(Rho)
  ENDIF
!------------------------------------------------------------
! Put Rho to disk
! 
  CALL Put_HGRho(Rho,'Rho',Args,0,.TRUE.)
!------------------------------------------------------------
! Printing
!
  CALL PChkSum(Rho,'Rho',Prog)
  CALL PPrint(Rho,'Rho')
!---------------------------------------------------
! Tidy up
! 
  IF( iSwitch /= 0) CALL Delete(Dmat)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(MD)
  CALL Delete_HGRho(Rho)
  CALL ShutDown(Prog)
!
END PROGRAM MakeRho

