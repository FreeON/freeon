!    FAST O(N lg N) COMPUTATION OF THE COULOMB MATRIX
!    Authors:  Matt Challacombe and CJ Tymczak
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
  USE QCTCThresholds
  USE PoleTree
#ifdef PERIODIC
  USE PBCFarField
#endif
  USE JGen
  USE NuklarE
#ifdef PARALLEL
  USE MondoMPI
  TYPE(DBCSR)                :: J,S,T1,P
#else
  TYPE(BCSR)                 :: J,S,T1,P
#endif
  REAL(DOUBLE)               :: E_Nuc_Tot
  TYPE(TIME)                 :: TimeMakeJ
  CHARACTER(LEN=4),PARAMETER :: Prog='QCTC'
!-------------------------------------------------------------------------------- 
! Start up macro
  CALL StartUp(Args,Prog)
! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
! Allocations 
  CALL NewBraBlok(BS)
! Get the Density for Poletree
  CALL Get(Rho,'Rho',Args,0)
! Set thresholds local to QCTC (for PAC and MAC)
  CALL SetLocalThresholds(Thresholds%TwoE)
! Initialize the auxiliary density arrays
  CALL InitRhoAux
! Setup global arrays for computation of multipole tensors
  CALL MultipoleSetUp(FFEll2)
! Build the global PoleTree representation of the total density
  CALL RhoToPoleTree
#ifdef PERIODIC
! Calculate the Number of Cells
  CALL SetCellNumber(GM)
! Set the electrostatic background 
  CALL PBCFarFieldSetUp(FFEll,PoleRoot)
#endif
! Delete the auxiliary density arrays
  CALL DeleteRhoAux
! Delete the Density
  CALL Delete(Rho)
! Allocate J
  CALL New(J)
! Compute the Coulomb matrix J in O(N Lg N)
  CALL Elapsed_Time(TimeMakeJ,'Init')
  CALL MakeJ(J)
  CALL Elapsed_TIME(TimeMakeJ,'Accum')
! Put J to disk
!
! now skiping filtration until after orthogonal Fock is xformed
!  CALL Filter(T1,J)

  IF(Args%C%C(2)=='Core')THEN
     CALL Put(J,TrixFile('V',Args))
  ELSE
     CALL Put(J,TrixFile('J',Args,0))
  ENDIF
! Compute the nuclear-total electrostatic energy
  E_Nuc_Tot=NukE()
  CALL Put(E_Nuc_Tot,'enn+ene',Tag_O=SCFCycl)
!---------------------------------------------------------------
! Printing
  CALL PPrint(E_Nuc_Tot,'NukE['//TRIM(SCFCycl)//']')
  IF(Args%C%C(2)=='Core')THEN
     CALL PChkSum(J,'V',Prog)
     CALL PPrint( J,'V')
     CALL Plot(   J,'V')
  ELSE
     CALL PChkSum(J,'J['//TRIM(SCFCycl)//']',Prog)
     CALL PPrint( J,'J['//TRIM(SCFCycl)//']')
     CALL Plot(   J,'J['//TRIM(SCFCycl)//']')
  ENDIF
! Tidy up
  CALL Delete(J)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(Args)
! Shutdown 
  CALL ShutDown(Prog)
END PROGRAM QCTC
