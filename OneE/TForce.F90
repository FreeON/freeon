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
!----------------------------------------------------------------------------------
!    Author: Matt Challacombe
!    COMPUTE THE FORCE CORESPONDING TO THE DERIVATIVE OF THE 
!    KINETIC ENERGY MATRIX, TForce=2*Tr{P.dT}
!----------------------------------------------------------------------------------
PROGRAM TForce
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
  USE BlokTrPdT
#ifdef PARALLEL
  USE MondoMPI
  TYPE(DBCSR)         :: P
#else
  TYPE(BCSR)          :: P
#endif
#ifdef PERIODIC 
  INTEGER             :: NC,NLay,MMLow,MMHig,MaxL
  REAL(DOUBLE)        :: Bx,By,Bz
#endif
  TYPE(AtomPair)      :: Pair
  TYPE(BSET)          :: BS
  TYPE(CRDS)          :: GM
  TYPE(ARGMT)         :: Args
  INTEGER             :: Q,R,AtA,AtB,NN,iSwitch,IStrtP,IStopP,LP,JP,MB,MA,A1,A2
  TYPE(HGRho)         :: Rho
  TYPE(DBL_VECT)      :: TFrc,Frc
  CHARACTER(LEN=6),PARAMETER :: Prog='TForce'
!------------------------------------------------------------------------------------- 
! Start up macro
!
  CALL StartUp(Args,Prog)
!----------------------------------------------
! Get basis set and geometry
!
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
  CALL Get(P,TrixFile('D',Args,0))
  CALL New(Frc,3*NAtoms)
  CALL New(TFrc,3*NAtoms)
  TFrc%D=Zero
!---------------------------------------------------------------------------
! TForce=2*Tr{P.dT} (Extra 2 to account for symmetry of T in the trace)
  DO AtA=1,NAtoms
     A1=3*(AtA-1)+1
     A2=3*AtA
     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1
        AtB=P%ColPt%I(JP)
        IF(SetAtomPair(GM,BS,AtA,AtB,Pair))THEN
           Q=P%BlkPt%I(JP)
           TFrc%D(A1:A2)=TFrc%D(A1:A2)+Four*TrPdT(BS,Pair,P%MTrix%D(Q:),AtA,AtB)
        ENDIF
     ENDDO
  ENDDO
!  CALL PPrint(TFrc,'TForce',Unit_O=6)
!---------------------------------------------------------------
! Update forces
  CALL Get(Frc,'GradE',Tag_O=CurGeom)
  Frc%D=Frc%D+TFrc%D
  CALL Put(Frc,'GradE',Tag_O=CurGeom)
  TFrcChk=SQRT(DOT_PRODUCT(TFrc%D,TFrc%D))
!  WRITE(*,*)' TFrcChk = ',TFrcChk
!------------------------------------------------------------------------------
! Tidy up 
  CALL Delete(P)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(TFrc)
  CALL ShutDown(Prog)
END PROGRAM TForce
