!------------------------------------------------------------------------------
!--  This code is part of the MondoSCF suite of programs for linear scaling 
!    electronic structure theory and ab initio molecular dynamics.
!
!--  Copyright (c) 2001, the Regents of the University of California.  
!    This SOFTWARE has been authored by an employee or employees of the 
!    University of California, operator of the Los Alamos National Laboratory 
!    under Contract No. W-7405-ENG-36 with the U.S. Department of Energy.  
!    The U.S. Government has rights to use, reproduce, and distribute this 
!    SOFTWARE.  The public may copy, distribute, prepare derivative works 
!    and publicly display this SOFTWARE without charge, provided that this 
!    Notice and any statement of authorship are reproduced on all copies.  
!    Neither the Government nor the University makes any warranty, express 
!    or implied, or assumes any liability or responsibility for the use of 
!    this SOFTWARE.  If SOFTWARE is modified to produce derivative works, 
!    such modified SOFTWARE should be clearly marked, so as not to confuse 
!    it with the version available from LANL.  The return of derivative works
!    to the primary author for integration and general release is encouraged. 
!    The first publication realized with the use of MondoSCF shall be
!    considered a joint work.  Publication of the results will appear
!    under the joint authorship of the researchers nominated by their
!    respective institutions. In future publications of work performed
!    with MondoSCF, the use of the software shall be properly acknowledged,
!    e.g. in the form "These calculations have been performed using MondoSCF, 
!    a suite of programs for linear scaling electronic structure theory and
!    ab initio molecular dynamics", and given appropriate citation.  
!------------------------------------------------------------------------------
!    Author:  Matt Challacombe
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
  USE BraBloks
  USE PoleTree
  USE JGen
  USE NuklarE
#ifdef PARALLEL
  USE MondoMPI
  TYPE(DBCSR)         :: J,T1
#else
  TYPE(BCSR)          :: J,T1
#endif
  TYPE(AtomPair)      :: Pair
  TYPE(ARGMT)         :: Args
  INTEGER             :: P,AtA,AtB,NN,iSwitch
  REAL(DOUBLE)        :: E_Nuc_Tot
  CHARACTER(LEN=4),PARAMETER :: Prog='QCTC'
!-------------------------------------------------------------------------------- 
! Start up macro
  CALL StartUp(Args,Prog)
! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
! Allocations 
  CALL NewBraBlok(BS)
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
  IF(Args%C%C(2)=='Core')THEN
     CALL Put(T1,TrixFile('V',Args))
  ELSE
     CALL Put(T1,TrixFile('J',Args,0))
  ENDIF
! Compute the nuclear-total electrostatic energy
  E_Nuc_Tot=NukE()
!  WRITE(*,*)' E_Nuc_Tot = ',E_Nuc_Tot
  CALL Put(E_Nuc_Tot,'enn+ene',Tag_O=SCFCycl)
!---------------------------------------------------------------
! Printing
!  CALL PPrint(E_Nuc_Tot,'NukE['//TRIM(SCFCycl)//']')
  IF(Args%C%C(2)=='Core')THEN
     CALL PChkSum(T1,'V',Prog)
     CALL PPrint( T1,'V')
     CALL Plot(   T1,'V')
  ELSE
     CALL PChkSum(T1,'J['//TRIM(SCFCycl)//']',Prog)
     CALL PPrint( T1,'J['//TRIM(SCFCycl)//']')
     CALL Plot(   T1,'J['//TRIM(SCFCycl)//']')
  ENDIF
!---------------------------------------------------------------
! Tidy up
  CALL Delete(J)
  CALL Delete(T1)
  CALL Delete(BS)
  CALL Delete(GM)
! Shutdown 
  CALL ShutDown(Prog)
END PROGRAM QCTC
