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
PROGRAM GradS
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
  USE GradSBlock
#ifdef PARALLEL
  USE MondoMPI
  TYPE(DBCSR)         :: T1,W,C,F,P
#else
  TYPE(BCSR)          :: T1,W,C,F,P
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
  TYPE(DBL_RNK2)      :: TrWdS
  CHARACTER(LEN=5),PARAMETER :: Prog='GradS'
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
  CALL New(MD,(/3,BS%NASym+2,BS%NASym+2,2*BS%NASym+2/),(/1,-1,-1,0/))
  CALL New(P)
  CALL New(TrWdS,(/3,NAtoms/))
  CALL New(W)
  CALL New(C)
  CALL New(F)
  CALL NewGradBraKetBlok(BS)
!-----------------------------------------------
! Initialize the matrix and associated indecies
!
  W%RowPt%I(1)=1;C%RowPt%I(1)=1
  CALL SetEq(W%MTrix,Zero)
  CALL SetEq(C%MTrix,Zero)
!--------------------------------------
! First compute the W matrix
!
  CALL Get(F,TrixFile('F',Args,0))
  CALL Get(P,TrixFile('D',Args,1))

  CALL Multiply(P,F,C)
  CALL Multiply(C,P,T1)
  CALL Filter(W,T1)

  CALL Delete(F)
  CALL Delete(C)
  CALL Delete(P)
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
!  Main loops over W
!
  TrWdS%D=Zero
#ifdef PARALLEL
  DO AtA=Beg%I(MyId),End%I(MyId)
#else
  DO AtA=1,NAtoms
#endif
     MA=BSiz%I(AtA)
     IStrtP=W%RowPt%I(AtA)
     IStopP=W%RowPt%I(AtA+1)-1
     DO JP=IStrtP,IStopP
        AtB=W%ColPt%I(JP)
        PP=W%BlkPt%I(JP)
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
                 TrWdS%D(:,AtA)=TrWdS%D(:,AtA) + &
                      GradSBlok(BS,MD,AtB,AtA,AtA,W%MTrix%D(PP:PP+MA*MB-1),Pair)
                 IF (AtA.NE.AtB) THEN
                    TrWdS%D(:,AtB)=TrWdS%D(:,AtB) + &
                         GradSBlok(BS,MD,AtB,AtA,AtB,W%MTrix%D(PP:PP+MA*MB-1),Pair)
                 ENDIF
              ENDIF
           ENDDO
#else
           TrWdS%D(:,AtA)=TrWdS%D(:,AtA) + &
                GradSBlok(BS,MD,AtB,AtA,AtA,W%MTrix%D(PP:PP+MA*MB-1),Pair)
           IF (AtA.NE.AtB) THEN
              TrWdS%D(:,AtB)=TrWdS%D(:,AtB) + &
                   GradSBlok(BS,MD,AtB,AtA,AtB,W%MTrix%D(PP:PP+MA*MB-1),Pair)
           ENDIF
        ENDIF
#endif
     ENDDO
  ENDDO
  DO AtA=1,NAtoms
     WRITE(*,*) TrWdS%D(1,AtA),TrWdS%D(2,AtA),TrWdS%D(3,AtA)
  ENDDO
!------------------------------------------------------------
! Put W to disk
!  
  CALL Filter(T1,W)
!  CALL Put(T1,TrixFile('W',Args),BlksName_O='nsi',Non0Name_O='nsm')
!-----------------------------------------------------------
! Printing
!
  CALL PChkSum(T1,'W',Prog)
  CALL PPrint( T1,'W')
  CALL Plot(   T1,'W')

!---------------------------------------------------
! Tidy up
! 
  CALL Delete(T1)
  CALL Delete(W)
  CALL DeleteGradBraKetBlok()
  CALL Delete_HGRho(Rho)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(MD)
  CALL ShutDown(Prog)
!
END PROGRAM GradS
