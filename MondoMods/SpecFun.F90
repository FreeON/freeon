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
!    SPECIAL (AND NOT SO-SPECIAL) FUNCTIONS
!    Author: Matt Challacombe
!------------------------------------------------------------------------------
!
!
!==============================================================================
MODULE GammaF0
   USE DerivedTypes
   IMPLICIT NONE
   INTEGER :: IF0
   INCLUDE "GammaDimensions.Inc"
   REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F0_0,F0_1,F0_2,F0_3,F0_4
   INCLUDE 'F0.Inc'
END MODULE 

MODULE GammaF1
   USE DerivedTypes
   IMPLICIT NONE
   INTEGER :: IF1
   INCLUDE 'GammaDimensions.Inc'
   REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F1_0,F1_1,F1_2,F1_3,F1_4
   INCLUDE 'F1.Inc'
END MODULE 

MODULE GammaF2
   USE DerivedTypes
   IMPLICIT NONE
   INTEGER :: IF2
   INCLUDE 'GammaDimensions.Inc'
   REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F2_0,F2_1,F2_2,F2_3,F2_4
   INCLUDE 'F2.Inc'
END MODULE 

MODULE GammaF3
   USE DerivedTypes
   IMPLICIT NONE
   INTEGER :: IF3
   INCLUDE 'GammaDimensions.Inc'
   REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F3_0,F3_1,F3_2,F3_3,F3_4
   INCLUDE 'F3.Inc'
END MODULE 

MODULE GammaF4
   USE DerivedTypes
   IMPLICIT NONE
   INTEGER :: IF4
   INCLUDE 'GammaDimensions.Inc'
   REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F4_0,F4_1,F4_2,F4_3,F4_4
   INCLUDE 'F4.Inc'
END MODULE 

MODULE GammaF5
   USE DerivedTypes
   IMPLICIT NONE
   INTEGER :: IF5
   INCLUDE 'GammaDimensions.Inc'
   REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F5_0,F5_1,F5_2,F5_3,F5_4
   INCLUDE 'F5.Inc'
END MODULE 

MODULE GammaF6
   USE DerivedTypes
   IMPLICIT NONE
   INTEGER :: IF6
   INCLUDE 'GammaDimensions.Inc'
   REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F6_0,F6_1,F6_2,F6_3,F6_4
   INCLUDE 'F6.Inc'
END MODULE 

MODULE GammaF7
   USE DerivedTypes
   IMPLICIT NONE
   INTEGER :: IF7
   INCLUDE 'GammaDimensions.Inc'
   REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F7_0,F7_1,F7_2,F7_3,F7_4
   INCLUDE 'F7.Inc'
END MODULE 

MODULE GammaF8
   USE DerivedTypes
   IMPLICIT NONE
   INTEGER :: IF8
   INCLUDE 'GammaDimensions.Inc'
   REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F8_0,F8_1,F8_2,F8_3,F8_4
   INCLUDE 'F8.Inc'
END MODULE 

MODULE GammaF9
   USE DerivedTypes
   IMPLICIT NONE
   INTEGER :: IF9
   INCLUDE 'GammaDimensions.Inc'
   REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F9_0,F9_1,F9_2,F9_3,F9_4
   INCLUDE 'F9.Inc'
END MODULE 

MODULE GammaF10
   USE DerivedTypes
   IMPLICIT NONE
   INTEGER :: IF10
   INCLUDE 'GammaDimensions.Inc'
   REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F10_0,F10_1,F10_2,F10_3,F10_4
   INCLUDE 'F10.Inc'
END MODULE 

MODULE GammaF11
   USE DerivedTypes
   IMPLICIT NONE
   INTEGER :: IF11
   INCLUDE 'GammaDimensions.Inc'
   REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F11_0,F11_1,F11_2,F11_3,F11_4
   INCLUDE 'F11.Inc'
END MODULE 

MODULE InvExp
   USE DerivedTypes
   IMPLICIT NONE
   INCLUDE "InverseExp.Inc"
   CONTAINS
!========================================================================================
!     Compute the inverse exponential function, EXP(-X)
!========================================================================================
      FUNCTION EXPInv(X)
         REAL(DOUBLE) :: EXPInv,X
         INTEGER      :: J ,I
         IF(X.GE.Exp_Switch)THEN
            EXPinv=0.0D0
         ELSE
            J=AINT(X*Exp_Grid)
            EXPinv=Exp_0(J)+X*(Exp_1(J)+X*(Exp_2(J)+X*(Exp_3(J)+X*Exp_4(J))))
         ENDIF
      END FUNCTION EXPinv
END MODULE

MODULE ErfFunk
   USE DerivedTypes
   IMPLICIT NONE
   INCLUDE "ErrorFunction.Inc"
   CONTAINS
!========================================================================================
!     Compute the erf function
!========================================================================================
      FUNCTION ERF(W)
         REAL(DOUBLE) :: ERF,W,X,Sgn
         INTEGER      :: J
         IF(W<0.0D0)THEN
            Sgn=-1.0D0
         ELSE
            Sgn=1.0D0 
         ENDIF        
         X=Sgn*W
         IF(X>Erf_Switch)THEN
            Erf=1.0D0
         ELSE
            J=AINT(X*Erf_Grid)
            Erf=Erf_0(J)+X*(Erf_1(J)+X*(Erf_2(J)+X*(Erf_3(J)+X*Erf_4(J))))
         ENDIF
         Erf=Sgn*Erf
      END FUNCTION ERF
!========================================================================================
!     The complimentary error function 
!========================================================================================
      FUNCTION ERFC(X)
         REAL(DOUBLE) :: ERFC,X
         ERFC=1.0D0-ERF(X)
      END FUNCTION ERFC
END MODULE

MODULE GammaFunctions
   USE GammaF0
   USE GammaF1
   USE GammaF2
   USE GammaF3
   USE GammaF4
   USE GammaF5
   USE GammaF6
   USE GammaF7
   USE GammaF8
   USE GammaF9
   USE GammaF10
   USE GammaF11
   USE ProcessControl
   IMPLICIT NONE
   INCLUDE "GammaGrid.Inc"
   CONTAINS
      FUNCTION GammaF(M,T)   
         INTEGER,      INTENT(IN) :: M
         REAL(DOUBLE), INTENT(IN) :: T
         INTEGER                  :: J
         REAL(DOUBLE)             :: GammaF
!-------------------------------------------------------------------------         
         IF(T>Gamma_Switch) &
            CALL Halt(' GammaF does not cover multipole approximation ')
         J=AINT(T*Gamma_Grid)
         SELECT CASE(M)
         CASE(0)
            GammaF=F0_0(J)+T*(F0_1(J)+T*(F0_2(J)+T*(F0_3(J)+T*F0_4(J))))
         CASE(1)
            GammaF=F1_0(J)+T*(F1_1(J)+T*(F1_2(J)+T*(F1_3(J)+T*F1_4(J))))
         CASE(2)
            GammaF=F2_0(J)+T*(F2_1(J)+T*(F2_2(J)+T*(F2_3(J)+T*F2_4(J))))
         CASE(3)
            GammaF=F3_0(J)+T*(F3_1(J)+T*(F3_2(J)+T*(F3_3(J)+T*F3_4(J))))
         CASE(4)
            GammaF=F4_0(J)+T*(F4_1(J)+T*(F4_2(J)+T*(F4_3(J)+T*F4_4(J))))
         CASE(5)
            GammaF=F5_0(J)+T*(F5_1(J)+T*(F5_2(J)+T*(F5_3(J)+T*F5_4(J))))
         CASE(6)
            GammaF=F6_0(J)+T*(F6_1(J)+T*(F6_2(J)+T*(F6_3(J)+T*F6_4(J))))
         CASE(7)
            GammaF=F7_0(J)+T*(F7_1(J)+T*(F7_2(J)+T*(F7_3(J)+T*F7_4(J))))
         CASE(8)
            GammaF=F8_0(J)+T*(F8_1(J)+T*(F8_2(J)+T*(F8_3(J)+T*F8_4(J))))
         CASE(9)
            GammaF=F9_0(J)+T*(F9_1(J)+T*(F9_2(J)+T*(F9_3(J)+T*F9_4(J))))
         CASE(10)
            GammaF=F10_0(J)+T*(F10_1(J)+T*(F10_2(J)+T*(F10_3(J)+T*F10_4(J))))
         CASE(11)
            GammaF=F11_0(J)+T*(F11_1(J)+T*(F11_2(J)+T*(F11_3(J)+T*F11_4(J))))
         CASE(12:)
            CALL Halt('Only M<=11 has been implimented in GammaF ')
         END SELECT
      END FUNCTION GammaF
END MODULE

MODULE SpecFun
   USE GammaFunctions
   USE ErfFunk
   USE InvExp
END MODULE
