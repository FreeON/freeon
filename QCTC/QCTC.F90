!==============================================================================
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!    FAST O(N lg N) COMPUTATION OF THE COULOMB MATRIX
!==============================================================================
PROGRAM QCTC
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
  USE JGen
  USE NukE
#ifdef PARALLEL
  USE MondoMPI
  TYPE(DBCSR)         :: J,T1
#else
  TYPE(BCSR)          :: J,T1
#endif
  TYPE(AtomPair)      :: Pair
  TYPE(ARGMT)         :: Args
  INTEGER             :: P,R,AtA,AtB,NN,iSwitch
  REAL(DOUBLE)        :: E_Nuc_Tot
  CHARACTER(LEN=4),PARAMETER :: Prog='QCTC'
!-------------------------------------------------------------------------------- 
! Start up macro
  CALL StartUp(Args,Prog)
! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
! Allocations 
  CALL New(MD,(/3,BS%NASym,BS%NASym,2*BS%NASym/),(/1,0,0,0/))
  CALL NewBraKetBlok(BS)
  CALL New(J)
! Setup global arrays for computation of multipole tensors
  CALL MultipoleSetUp(SPell2)
! Build the global PoleTree representation of the total density
  CALL RhoToPoleTree(Args)
! Set the electrostatic background 
  CALL FarField(PoleRoot)
! Compute the Coulomb matrix J in O(N Lg N)
  CALL MakeJ(J)
! Put J to disk
  CALL Filter(T1,J)
  CALL Put(T1,TrixFile('J',Args,0))
! Compute the nuclear-total electrostatic energy
  E_Nuc_Tot=NukE()
  CALL Put(E_Nuc_Tot,'enn+ene',Tag_O=SCFCycl)
! Printing
  CALL PChkSum(T1,'J['//TRIM(SCFCycl)//']',Prog,Unit_O=6)
  CALL PChkSum(T1,'J['//TRIM(SCFCycl)//']',Prog)
  CALL PPrint( T1,'J['//TRIM(SCFCycl)//']')
  CALL Plot(   T1,'J['//TRIM(SCFCycl)//']')
! Tidy up
  CALL Delete(J)
  CALL Delete(T1)
  CALL Delete(BS)
  CALL Delete(GM)
! Shutdown 
  CALL ShutDown(Prog)
END PROGRAM QCTC
