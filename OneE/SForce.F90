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
!    Author: Matt Challacombe and CJ Tymczak
!    COMPUTE THE FORCE DUE TO CHANGES IN THE DENSITY MATRIX:
!    dP.(2T+J+K)=-2*dS.P.F.P (Early work by McWeeny, need a cite...)
!------------------------------------------------------------------------------
PROGRAM SForce
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
  USE BraBloks
  USE BlokTrWdS
#ifdef PARALLEL
  USE MondoMPI
  TYPE(DBCSR)         :: T1,F,P
#else
  TYPE(BCSR)          :: T1,F,P
#endif
#ifdef PERIODIC 
  INTEGER             :: NC
  REAL(DOUBLE)        :: Bx,By,Bz
#endif
  TYPE(AtomPair)      :: Pair
  TYPE(BSET)          :: BS
  TYPE(CRDS)          :: GM
  TYPE(ARGMT)         :: Args
  INTEGER             :: Q,R,AtA,AtB,NN,iSwitch,IStrtP,IStopP,LP,JP,MB,MA,A1,A2
  TYPE(HGRho)         :: Rho
  TYPE(DBL_VECT)      :: SFrc
  CHARACTER(LEN=6),PARAMETER :: Prog='SForce'
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
  CALL New(P)
  CALL New(F)
  CALL New(T1)
!--------------------------------------
! Compute W=P.F.P
  CALL Get(P,TrixFile('D',Args,0))
!
  CALL Get(F,TrixFile('F',Args,0))
!  CALL Get(F,TrixFile('T',Args))
  PrintFlags%Fmt=DEBUG_DOUBLES
  CALL Multiply(P,F,T1)       ! T1:=P.F
  CALL Multiply(T1,P,F)       ! F:=P.F.P
  CALL Filter(P,F)            ! P=Filter[P.F.P]      
  CALL Delete(F)
  CALL Delete(T1)
#ifdef PERIODIC
!-----------------------------------------------
! Calculate the Number of Cells
!
  CALL SetCellNumber(GM)
#endif
!  CALL PPrint(P,'W',Unit_O=6)
!---------------------------------------------------------------------------
! SForce=-2*Tr{P.F.P.dS} (Extra 2 to account for symmetry of S in the trace)
  CALL New(SFrc,3*NAtoms)
  SFrc%D=Zero
  CALL NewBraBlok(BS,Gradients_O=.TRUE.)
  DO AtA=1,NAtoms
     A1=3*(AtA-1)+1
     A2=3*AtA
     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1
        AtB=P%ColPt%I(JP)
        IF(SetAtomPair(GM,BS,AtA,AtB,Pair))THEN
           Q=P%BlkPt%I(JP)
#ifdef PERIODIC
           Bx = Pair%B(1)
           By = Pair%B(2)           
           Bz = Pair%B(3)
           DO NC=1,CS%NCells
              Pair%B(1)=Bx+CS%CellCarts%D(1,NC)
              Pair%B(2)=By+CS%CellCarts%D(2,NC)
              Pair%B(3)=Bz+CS%CellCarts%D(3,NC)
              Pair%AB2=(Pair%A(1)-Pair%B(1))**2+(Pair%A(2)-Pair%B(2))**2+(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair)) THEN
                 SFrc%D(A1:A2)=SFrc%D(A1:A2)-Four*TrWdS(BS,Pair,P%MTrix%D(Q:))
              ENDIF
           ENDDO
#else
           SFrc%D(A1:A2)=SFrc%D(A1:A2)-Four*TrWdS(BS,Pair,P%MTrix%D(Q:))
#endif
        ENDIF
     ENDDO
  ENDDO
!------------------------------------------------------------------------------
! Put SForce to info file
  CALL PPrint(SFrc,'SForce')
  CALL Put(SFrc,'GradE',Tag_O=CurGeom)
  SFrcChk=SQRT(DOT_PRODUCT(SFrc%D,SFrc%D))
!------------------------------------------------------------------------------
! Tidy up 
  CALL Delete(P)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(SFrc)
  CALL DeleteBraBlok(Gradients_O=.TRUE.)
  CALL ShutDown(Prog)
END PROGRAM SForce
