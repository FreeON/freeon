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
!    DRIVER ROUTINES FOR DYNAMICS, OPTIMIZATION AND TS METHODS
!    Author:  Matt Challacombe
!------------------------------------------------------------------------------
MODULE DrvFrcs
  USE DerivedTypes
  USE GlobalScalars
#ifdef NAG
  USE F90_UNIX_PROC
#endif
  USE PrettyPrint
  USE SCFLocals
  USE Overlay
  USE ParsingKeys
  USE Macros
  USE Functionals
  USE DrvSCFs
  IMPLICIT NONE 
  CONTAINS
!========================================================================================
!
!========================================================================================
     SUBROUTINE NewStep(Ctrl,IChk)
        TYPE(SCFControls)  :: Ctrl
        INTEGER            :: IChk,ICyc,IGeo
!----------------------------------------------------------------------------------------
        ICyc=Ctrl%Current(1)
        IGeo=Ctrl%Current(3)
        CtrlVect=SetCtrlVect(Ctrl,'QuNew')
        IF(IChk==1)THEN
!          Current geometry is good, compute gradients and update Hessian
           CALL Invoke('SForce',     CtrlVect)
           CALL Invoke('TForce',     CtrlVect)
           CALL Invoke('JForce',  CtrlVect)
           IF(HasDFT(Ctrl%Model(Ctrl%Current(2)))) &
              CALL Invoke('XCForce',     CtrlVect)       
           IF(HasHF(Ctrl%Model(Ctrl%Current(2))))THEN
!              CALL Invoke('XForce',     CtrlVect)       
              CALL MondoHalt(-999,'ONX Gradients not in yet, soon... ')
           ENDIF
           CALL Invoke('BFGSHess',   CtrlVect)
        ENDIF
        CALL Invoke('NewStep',    CtrlVect)
   END SUBROUTINE
!========================================================================================
!
!========================================================================================
END MODULE




