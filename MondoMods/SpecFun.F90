!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
MODULE SpecFun
   USE DerivedTypes
   IMPLICIT NONE
   INTEGER :: I
   INCLUDE 'MMA/Functions/Erf.Inc'
   INCLUDE 'MMA/Functions/Exp.Inc'
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
!========================================================================================
!     Compute the erf function
!========================================================================================
      FUNCTION ERF(W)
         REAL(DOUBLE) :: ERF,W,X,Sgn
         INTEGER      :: J
!         REAL(DOUBLE),EXTERNAL :: DERF
         IF(W<0.0D0)THEN
            Sgn=-1.0D0
         ELSE
            Sgn=1.0D0 
         ENDIF        
         X=Sgn*W
         IF(X.GE.Erf_Switch)THEN
            Erf=1.0D0
         ELSE
            J=AINT(X*Erf_Grid)
            Erf=Erf_0(J)+X*(Erf_1(J)+X*(Erf_2(J)+X*Erf_3(J))) 
         ENDIF
         Erf=Sgn*Erf
      END FUNCTION ERF
!========================================================================================
!     The complimentary error function 
!========================================================================================
      FUNCTION ERFC(X)
         REAL(DOUBLE)           :: ERFC,X
         ERFC=1.0D0-ERF(X)
      END FUNCTION ERFC
END MODULE
