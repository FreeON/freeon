!------------------------------------------------------------------------------
!--  This source code is part of the MondoSCF suite of programs for 
!    linear scaling electronic structure theory and ab initio molecular 
!    dynamics.
!
!--  Matt Challacombe
!    Los Alamos National Laboratory
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
!    FICTIOUS, PURELY QUADRATIC ENERGY AND GRADIENTS FOR DEBUGING QuNew
!    Author(s):  Matt Challacombe
!------------------------------------------------------------------------------
PROGRAM QuGrad
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalCharacters
   USE InOut
   USE PrettyPrint
   USE MemMan
   USE Parse
   USE Macros
   USE SetXYZ
   USE LinAlg
   USE Functionals
   USE MatFunk
#ifdef PARALLEL
   USE MondoMPI
#endif
   IMPLICIT NONE
   TYPE(ARGMT)                :: Args
   CHARACTER(LEN=6),PARAMETER :: Prog='QuGrad'
   TYPE(BSet)                 :: BS
   TYPE(CRDS)                 :: GM
   TYPE(DBL_VECT)             :: G,X
   INTEGER                    :: N3,I,J,K,IGeom
   REAL(DOUBLE)               :: A
!----------------------------------------------------------------------------
   CALL StartUp(Args,Prog)
   IGeom=Args%I%I(3)
!  Override default blocking...
   N3=3*NAtoms
   NBasF=N3
   PrintFlags%Mat=DEBUG_MATRICES
   CALL New(BSiz,NAtoms)
   CALL New(OffS,NAtoms)
   BSiz%I=3
   OffS%I(1)=1
   DO I=2,NAtoms
      OffS%I(I)=OffS%I(I-1)+3
   ENDDO      
!  Allocations
   CALL New(G,N3)
   CALL New(X,N3)
!
   CALL Get(GM,Tag_O=CurGeom)
   K=0
   DO I=1,NAtoms
      DO J=1,3
         K=K+1
         X1%D(K)=GM1%Carts%D(J,I)
      ENDDO
   ENDDO   
!
   F=Zero
   K=0
   DO I=1,NAtoms
      A=10**(DBLE(NAtoms)/2.D0-DBLE(I))
      DO J=1,3
         K=K+1
         Co=A*J
         F=F+Co*X%D(K,IGeom)**2
         G%D(K)=2.0D0*A*X%D(K)
      ENDDO
   ENDDO
   CALL Put(G,'GradE',Tag_O=CurGeom)
   CALL Shutdown(Prog)
!
END PROGRAM
