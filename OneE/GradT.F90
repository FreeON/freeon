!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe and  Michael R. Salazar
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!    COMPUTE THE GRADIENT
!
PROGRAM GradT
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
  USE BraKetBloks
  USE GradTBlock
#ifdef PARALLEL
  USE MondoMPI
  TYPE(DBCSR)         :: T1,P
#else
  TYPE(BCSR)          :: T1,P
#endif
#ifdef PERIODIC 
  LOGICAL             :: NotMakeBlock
  INTEGER             :: NC,NLay,MMLow,MMHig,MaxL
  REAL(DOUBLE)        :: Bx,By,Bz
#endif
  TYPE(AtomPair)      :: Pair
  TYPE(BSET)          :: BS
  TYPE(CRDS)          :: GM
  TYPE(DBL_RNK4)      :: MD
  TYPE(ARGMT)         :: Args
  INTEGER             :: PP,R,AtA,AtB,NN,iSwitch,IStrtP,IStopP,LP,JP,MB,MA
  REAL(DOUBLE)        :: E_Nuc_Tot,E_Nuc_TotX,E_Nuc_TotY,E_Nuc_TotZ,VBlokX,VBlokY,VBlokZ
  TYPE(HGRho)         :: Rho
  TYPE(DBL_RNK2)      :: TrPdT
  CHARACTER(LEN=5),PARAMETER :: Prog='GradT'
!---------------------------------------------- 
! Start up macro
!
  CALL StartUp(Args,Prog)
!----------------------------------------------
! Get basis set and geometry
!
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
!---------------------------------------------- 
! Allocations 
!
  CALL New(MD,(/3,BS%NASym+4,BS%NASym+4,2*BS%NASym+4/),(/1,-2,-2,0/))
  CALL New(P)
  CALL Get(P,TrixFile('D',Args,1))
  CALL New(TrPdT,(/3,NAtoms/))
  CALL NewGradBraKetBlok(BS)
!-----------------------------------------------
! Initialize the matrix and associated indecies
!
  P%RowPt%I(1)=1
#ifdef PERIODIC
!-----------------------------------------------
! Calculate the Number of Cells
!-----------------------------------------------
  CALL SetCellNumber(GM)
!-----------------------------------------------
! Set Up the Multipoles  
!-----------------------------------------------
  MaxL = 32
  CALL MMSetup(MaxL,GM,Rho)
!-----------------------------------------------
! Calculate the Multipoles for Rho
!-----------------------------------------------
  CALL CalMMRho(Rho)
#endif
!-----------------------------------------------
!  Main loops over P
!
  TrPdT%D=Zero
#ifdef PARALLEL
  DO AtA=Beg%I(MyId),End%I(MyId)
#else
  DO AtA=1,NAtoms
#endif
     MA=BSiz%I(AtA)
     IStrtP=P%RowPt%I(AtA)
     IStopP=P%RowPt%I(AtA+1)-1
     DO JP=IStrtP,IStopP
        AtB=P%ColPt%I(JP)
        PP=P%BlkPt%I(JP)
        MB=BSiz%I(AtB)
        IF(SetAtomPair(GM,BS,AtB,AtA,Pair)) THEN
#ifdef PERIODIC
           Bx = Pair%B(1)
           By = Pair%B(2)           
           Bz = Pair%B(3)
           NotMakeBlock = .TRUE.
           DO NC = 1,CS%NCells
              Pair%B(1) = Bx+CS%CellCarts%D(1,NC)
              Pair%B(2) = By+CS%CellCarts%D(2,NC)
              Pair%B(3) = Bz+CS%CellCarts%D(3,NC)
              Pair%AB2  = (Pair%A(1)-Pair%B(1))**2 &
                        + (Pair%A(2)-Pair%B(2))**2 &
                        + (Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair)) THEN
                 TrPdT%D(:,AtA)=TrPdT%D(:,AtA) + &
                      GradTBlok(BS,MD,AtB,AtA,AtA,P%MTrix%D(PP:PP+MA*MB-1),Pair)
                 IF (AtA.NE.AtB) THEN
                    TrPdT%D(:,AtB)=TrPdT%D(:,AtB) + &
                         GradTBlok(BS,MD,AtB,AtA,AtB,P%MTrix%D(PP:PP+MA*MB-1),Pair)
                 ENDIF
              ENDIF
           ENDDO
#else
           TrPdT%D(:,AtA)=TrPdT%D(:,AtA) + &
                GradTBlok(BS,MD,AtB,AtA,AtA,P%MTrix%D(PP:PP+MA*MB-1),Pair)
           IF (AtA.NE.AtB) THEN
              TrPdT%D(:,AtB)=TrPdT%D(:,AtB) + &
                   GradTBlok(BS,MD,AtB,AtA,AtB,P%MTrix%D(PP:PP+MA*MB-1),Pair)
           ENDIF
        ENDIF
#endif
     ENDDO
  ENDDO
  DO AtA=1,NAtoms
     WRITE(*,*) TrPdT%D(1,AtA),TrPdT%D(2,AtA),TrPdT%D(3,AtA)
  ENDDO
!------------------------------------------------------------
! Put P to disk
!  
  CALL Filter(T1,P)
!  CALL Put(T1,TrixFile('W',Args),BlksName_O='nsi',Non0Name_O='nsm')
!-----------------------------------------------------------
! Printing
!
  CALL PChkSum(T1,'P',Prog)
  CALL PPrint( T1,'P')
  CALL Plot(   T1,'P')

!---------------------------------------------------
! Tidy up
! 
  CALL Delete(T1)
  CALL Delete(P)
  CALL DeleteGradBraKetBlok()
  CALL Delete_HGRho(Rho)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(MD)
  CALL ShutDown(Prog)
!
END PROGRAM GradT
