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
!    Author: Matt Challacombe
!    FAST O(N Lg N) COMPUTATION OF GRADIENTS OF THE COULOMB ENERGY 
!    WRT TO NUCLEAR COORDINATES
!==============================================================================
PROGRAM JForce
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
  USE BlokTrPdJ
  USE NuklarE
#ifdef PARALLEL
  USE MondoMPI
  TYPE(DBCSR)         :: P
#else
  TYPE(BCSR)          :: P
#endif
  TYPE(AtomPair)      :: Pair
  TYPE(ARGMT)         :: Args
  INTEGER             :: Q,AtA,AtB,A1,A2,iSwitch
  REAL(DOUBLE)        :: E_Nuc_Tot
  TYPE(DBL_VECT)      :: Frc,JFrc
  REAL(DOUBLE),DIMENSION(3) :: FrcAB
  CHARACTER(LEN=6),PARAMETER :: Prog='JForce'
!-------------------------------------------------------------------------------- 
! Start up macro
  CALL StartUp(Args,Prog)
! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
! Allocations 
  CALL NewBraBlok(BS,Gradients_O=.TRUE.)
  CALL New(Frc,3*NAtoms)
  CALL New(JFrc,3*NAtoms)
  CALL Get(P,TrixFile('D',Args,1))
! Setup global arrays for computation of multipole tensors
  CALL MultipoleSetUp(SPell2)
! Build the global PoleTree representation of the total density
  CALL RhoToPoleTree(Args)
! Set the electrostatic background 
  CALL FarField(PoleRoot)
!--------------------------------------------------------------
! Compute the Coulomb contribution to the force in O(N Lg N)
!
  JFrc%D=Zero
  DO AtA=1,NAtoms
     MA=BSiz%I(AtA)
     A1=3*(AtA-1)+1
     A2=3*AtA
     JFrc%D(A1:A2)=dNukE(AtA)
     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1
        AtB=P%ColPt%I(JP)
        IF(SetAtomPair(GM,BS,AtA,AtB,Pair))THEN
           Q=P%BlkPt%I(JP)
           JFrc%D(A1:A2)=JFrc%D(A1:A2)+TrPdJ(Pair,P%MTrix%D(Q:),AtA,AtB)
        ENDIF
     ENDDO
  ENDDO
  JFrc%D=Two*JFrc%D
!  CALL PPrint(JFrc,'JForce',Unit_O=6)
!---------------------------------------------------------------
! Update forces
  CALL Get(Frc,'GradE',Tag_O=CurGeom)
  Frc%D=Frc%D+JFrc%D
  CALL Put(Frc,'GradE',Tag_O=CurGeom)
  JFrcChk=SQRT(DOT_PRODUCT(JFrc%D,JFrc%D))
!  WRITE(*,*)' JFrcChk = ',JFrcChk
!---------------------------------------------------------------
! Tidy up
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(P)
  CALL Delete(Frc)
  CALL Delete(JFrc)
  CALL DeleteBraBlok(Gradients_O=.TRUE.)
  CALL ShutDown(Prog)
END PROGRAM JForce
