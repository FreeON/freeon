!    TESTS PRETABULATED SPECIAL FUNCTION APPROXIMATIONS
!    Author: Matt Challacombe
!------------------------------------------------------------------------------
PROGRAM FTest
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE Parse
  USE MemMan
  USE GammaFunctions
  INTEGER,PARAMETER             :: LMx=200
  REAL(DOUBLE),DIMENSION(0:LMx) :: F
  REAL(DOUBLE),DIMENSION(0:20)  :: MaxErr,MaxErr2,TestF
  REAL(DOUBLE)                  :: DT,T,T2,T4,ET,TwoT,FJ,Test
  INTEGER                       :: NSamples,K,L,J
  REAL(DOUBLE),EXTERNAL         :: DErf,DErfC
!---------------------------------------------------------------------
  INCLUDE "InverseExp.Inc" 
  INCLUDE "ErrorFunction.Inc"
  NSamples=100000
!---------------------------------------------------------------------
! TEST GAMMA FUNCTIONS
  T=0.D0
  MaxErr=-1.D10
  DT=Gamma_Switch/DBLE(NSamples+1)
  DO K=1,NSamples
     T=T+DT
!    Downward recursion: F_{j}(T)=(2*T*F_{j+1}+Exp[-T])/(2*j+1)
     TwoT=Two*T
     ET=EXP(-T)
     FJ=Zero
     DO L=LMx,0,-1
        F(L)=FJ
        FJ=(TwoT*F(L)+ET)/(Two*DBLE(L)-One)
     ENDDO
     DO J=0,11
        MaxErr(J)=MAX(MaxErr(J),ABS(GammaF(J,T)-F(J)))
     ENDDO
  ENDDO
  WRITE(*,*)' Maximum absolute error in gamma function interpolation:'
  DO J=0,11
     WRITE(*,*)J,MaxErr(J)
  ENDDO
!=================================================================================
  T=0.D0
  MaxErr=-1.D10
  DT=Exp_Switch/DBLE(NSamples+1)
  DO K=1,NSamples
     T=T+DT
     IF(T>Exp_Switch)STOP '2: Bad logic in FTest'
     J=AINT(T*Exp_Grid)
     Test=Exp_0(J)+T*(Exp_1(J)+T*(Exp_2(J)+T*(+Exp_3(J)+T*Exp_4(J))))
     MaxErr(0)=MAX(MaxErr(0),ABS(Test-EXP(-T)))
  ENDDO
  WRITE(*,*)' Maximum absolute error in inverse exponential interpolation = ',MaxErr(0)
!=================================================================================
!=================================================================================
  T=0.D0
  MaxErr=-1.D10
  DT=Erf_Switch/DBLE(NSamples+1)
  DO K=1,NSamples
     T=T+DT
     IF(T>Erf_Switch)STOP '3: Bad logic in FTest'
     J=AINT(T*Erf_Grid)
     Test=Erf_0(J)+T*(Erf_1(J)+T*(Erf_2(J)+T*(+Erf_3(J)+T*Erf_4(J))))
     MaxErr(0)=MAX(MaxErr(0),ABS(Test-DERF(T)))
  ENDDO
  WRITE(*,*)' Maximum absolute error in error function interpolation = ',MaxErr(0)
!=================================================================================
END PROGRAM
