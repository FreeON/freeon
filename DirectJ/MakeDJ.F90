!==============================================================================
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe and  C. J. Tymczak
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!    COMPUTE THE COULOMB MATRIX
!==============================================================================
PROGRAM MakeDJ
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
  USE JBlock
  USE Multipoles
#ifdef PARALLEL
  USE MondoMPI
  TYPE(DBCSR)         :: Jmat,T1
#else
  TYPE(BCSR)          :: Jmat,T1
  TYPE(BCSR)          :: Jold,Jnew
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
  INTEGER             :: P,R,AtA,AtB,NN,iSwitch
  REAL(DOUBLE)        :: E_Nuc_Tot
  TYPE(HGRho)         :: Rho                       
  CHARACTER(LEN=6),PARAMETER :: Prog='MakeDJ'
!---------------------------------------------- 
! Start up macro
  CALL StartUp(Args,Prog)
! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
! Get the Density 
  CALL Get_HGRho(Rho,'Rho',Args,0,.TRUE.)
! Allocations 
  CALL New(MD,(/3,BS%NASym,BS%NASym,2*BS%NASym/),(/1,0,0,0/))
  CALL NewBraKetBlok(BS)
  CALL New(Jmat)
!
#ifdef PERIODIC
!-----------------------------------------------
! Calculate the Number of Cells
  CALL SetCellNumber(GM)
! Set Up the Multipoles  
  MaxL = 16
  CALL MMSetup(MaxL,GM,Rho)
! Calculate the Multipoles for Rho
  CALL CalMMRho(Rho)
! Calculate the Multipoles for Rho
  CALL CalQRho(Rho)
#endif
!----------------------------------------------
! Set up the appropriate action
!----------------------------------------------
  IF(Args%C%C(2)=='Core') THEN
     iSwitch = 0
  ELSEIF(Args%C%C(2)=='Switch') THEN
     iSwitch = 1
  ELSEIF(Args%C%C(2)=='Direct') THEN
     iSwitch = 1
  ELSEIF(Args%C%C(2)=='InkFok') THEN
     iSwitch = 2
  ELSE
     CALL MondoHalt(-100,' Inappropriate Action in MakeDJ:'//Args%C%C(2))
  ENDIF
!-----------------------------------------------
! Initialize the matrix and associated indecies
!
  P=1; R=1; Jmat%RowPt%I(1)=1
  CALL SetEq(Jmat%MTrix,Zero)
!-----------------------------------------------
!  Main loops
!
  Jmat%NAtms=NAtoms
  DO AtA=1,NAtoms
     DO AtB=1,NAtoms
        IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
           NN = Pair%NA*Pair%NB
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
                 Jmat%MTrix%D(R:R+NN-1)=Jmat%MTrix%D(R:R+NN-1) &
                                       +Two*JBlok(BS,MD,Pair,Rho)

                 NotMakeBlock = .FALSE.
              ENDIF
           ENDDO
#else
           Jmat%MTrix%D(R:R+NN-1)=Two*JBlok(BS,MD,Pair,Rho)
#endif
           Jmat%ColPt%I(P)=AtB
           Jmat%BlkPt%I(P)=R
           R=R+NN 
           P=P+1 
           Jmat%RowPt%I(AtA+1)=P        
           IF(R>MaxNon0.OR.P>MaxBlks) &
                CALL Halt(' BCSR dimensions blown in MakeDJ ')
#ifdef PERIODIC
           IF(NotMakeBlock) &
                CALL Halt(' Making a Zero Block in MakeDJ ')
#endif
        ENDIF
     ENDDO
  ENDDO
  Jmat%NBlks=P-1
  Jmat%NNon0=R-1
!------------------------------------------------------------
! Put J to disk
!
  CALL Filter(T1,Jmat)
  IF(iSwitch == 0) THEN
     CALL Put(T1,TrixFile('V',Args))
  ELSE
     CALL Put(T1,TrixFile('J',Args,0))
  ENDIF
!------------------------------------------------------------
! Calculate the Nuc-Nuc + Nuc-Elec Piece
!
  E_Nuc_Tot = Zero
  DO AtA = 1,NAtoms 
     E_Nuc_Tot = E_Nuc_Tot+Two*VBlok(AtA,BS,MD,Rho)
  ENDDO
  WRITE(*,*)' NukE = ',E_Nuc_Tot
  CALL Put(E_Nuc_Tot,'enn+ene',Tag_O=SCFCycl)
!------------------------------------------------------------
! Printing
!
  IF(iSwitch==0) THEN
     CALL PChkSum(JMat,'Vne',Prog)
     CALL PPrint(Jmat,'Vne')
  ELSEIF(iSwitch==1) THEN
     CALL PChkSum(JMat,'Vte['//TRIM(SCFCycl)//']',Prog)
     CALL PPrint( JMat,'Vte['//TRIM(SCFCycl)//']')
     CALL Plot(   JMat,'Vte['//TRIM(SCFCycl)//']')
  ELSEIF(iSwitch==2) THEN
     CALL PChkSum(JMat,'DeltaVee',Prog)
     CALL PPrint(Jmat,'DeltaVee')
  ENDIF
!---------------------------------------------------
! Tidy up
! 
  CALL Delete(Jmat)
  CALL Delete(T1)
  CALL DeleteBraKetBlok()
  CALL Delete_HGRho(Rho)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(MD)
  CALL ShutDown(Prog)
!
END PROGRAM MakeDJ
