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
PROGRAM GradQCTC
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE Globals
  USE AtomPairs
  USE BraKetBloks
  USE PoleTree
  USE GradJGen
  USE GradNukE
#ifdef PARALLEL
  USE MondoMPI
  TYPE(DBCSR)         :: J,T1,P
#else
  TYPE(BCSR)          :: J,T1,P
#endif
  TYPE(AtomPair)      :: Pair
  TYPE(ARGMT)         :: Args
  INTEGER             :: PP,R,AtA,AtB,NN,iSwitch,IStrtP,IStopP,LP,JP,MB,MA
  TYPE(DBL_RNK2)      :: TrPdJ,dE_Nuc
  CHARACTER(LEN=8),PARAMETER :: Prog='GradQCTC'
!-------------------------------------------------------------------------------- 
! Start up macro
  CALL StartUp(Args,Prog)
! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
!  CALL Get(Rho,'Rho',Args,0)
!  CALL Print_HGRho(Rho,'Rho',Unit_o=6)
!  CALL DELETE(Rho)
! Allocations
  CALL New(MD,(/3,BS%NASym+2,BS%NASym+2,2*BS%NASym+2/),(/1,-1,-1,0/))
  CALL NewGradBraKetBlok(BS)
  CALL New(P)
  CALL Get(P,TrixFile('D',Args,1))
  CALL New(TrPdJ,(/3,NAtoms/))
  CALL New(dE_Nuc,(/3,NAtoms/))
! Setup global arrays for computation of multipole tensors
  CALL MultipoleSetUp(SPell2)
! Build the global PoleTree representation of the total density
  CALL RhoToPoleTree(Args)
! Set the electrostatic background 
  CALL FarField(PoleRoot)
! Compute the Derivative of the Coulomb matrix in O(N Lg N)
!-----------------------------------------------
! Initialize the matrix and associated indecies
!
  P%RowPt%I(1)=1
!-----------------------------------------------
!  Main loops over P
!
  TrPdJ%D=Zero
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
                 TrPdJ%D(:,AtA)=TrPdJ%D(:,AtA) + &
                      GradJBlock(BS,MD,AtB,AtA,AtA,P%MTrix%D(PP:PP+MA*MB-1),Pair,PoleRoot)
                 IF (AtA.NE.AtB) THEN
                    TrPdJ%D(:,AtB)=TrPdJ%D(:,AtB) + &
                         GradJBlock(BS,MD,AtB,AtA,AtB,P%MTrix%D(PP:PP+MA*MB-1),Pair,PoleRoot)
                 ENDIF
              ENDIF
           ENDDO
#else
           TrPdJ%D(:,AtA)=TrPdJ%D(:,AtA) + &
                GradJBlock(BS,MD,AtB,AtA,AtA,P%MTrix%D(PP:PP+MA*MB-1),Pair,PoleRoot)
           IF (AtA.NE.AtB) THEN
              TrPdJ%D(:,AtB)=TrPdJ%D(:,AtB) + &
                   GradJBlock(BS,MD,AtB,AtA,AtB,P%MTrix%D(PP:PP+MA*MB-1),Pair,PoleRoot)
           ENDIF
        ENDIF
#endif
     ENDDO
  ENDDO
!------------------------------------------------------------
! Calculate the derivative of Vnn and Vne
  dE_Nuc%D=Zero
  CALL GradNukE(dE_Nuc)
!
  DO AtA=1,NAtoms
!     WRITE(*,*) AtA,TrPdJ%D(1,AtA),TrPdJ%D(2,AtA),TrPdJ%D(3,AtA),'ELECTRONIC'
!     WRITE(*,*) AtA,dE_Nuc%D(1,AtA),dE_Nuc%D(2,AtA),dE_Nuc%D(3,AtA),'NUCLEAR'
     WRITE(*,*) TrPdJ%D(1,AtA)+dE_Nuc%D(1,AtA),TrPdJ%D(2,AtA)+dE_Nuc%D(2,AtA), &
          TrPdJ%D(3,AtA)+dE_Nuc%D(3,AtA)
  ENDDO
!------------------------------------------------------------
! Put P to disk
!  
  CALL Filter(T1,P)
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
END PROGRAM GradQCTC

