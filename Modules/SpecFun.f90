!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
!    SPECIAL (AND NOT SO-SPECIAL) FUNCTIONS
!    Authors: Matt Challacombe and C.J. Tymczak
!==============================================================================
MODULE GammaAssymp
  USE DerivedTypes
  IMPLICIT NONE
  INCLUDE "GammaAssymptotics.Inc"
END MODULE

MODULE SpecFunMesh
  USE DerivedTypes
  IMPLICIT NONE
  INCLUDE "GammaGrid.Inc"
  INCLUDE "GammaDimensions.Inc"
END MODULE SpecFunMesh

MODULE GammaF0
  USE DerivedTypes
  USE SpecFunMesh
  USE GammaAssymp
  IMPLICIT NONE
  INTEGER :: IF0
  REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F0_0,F0_1,F0_2,F0_3,F0_4
  INCLUDE 'F0.Inc'
END MODULE GammaF0

MODULE GammaF1
  USE DerivedTypes
  USE SpecFunMesh
  USE GammaAssymp
  IMPLICIT NONE
  INTEGER :: IF1
  REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F1_0,F1_1,F1_2,F1_3,F1_4
  INCLUDE 'F1.Inc'
END MODULE GammaF1

MODULE GammaF2
  USE DerivedTypes
  USE SpecFunMesh
  USE GammaAssymp
  IMPLICIT NONE
  INTEGER :: IF2
  REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F2_0,F2_1,F2_2,F2_3,F2_4
  INCLUDE 'F2.Inc'
END MODULE GammaF2

MODULE GammaF3
  USE DerivedTypes
  USE SpecFunMesh
  USE GammaAssymp
  IMPLICIT NONE
  INTEGER :: IF3
  REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F3_0,F3_1,F3_2,F3_3,F3_4
  INCLUDE 'F3.Inc'
END MODULE GammaF3

MODULE GammaF4
  USE DerivedTypes
  USE SpecFunMesh
  USE GammaAssymp
  IMPLICIT NONE
  INTEGER :: IF4
  REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F4_0,F4_1,F4_2,F4_3,F4_4
  INCLUDE 'F4.Inc'
END MODULE GammaF4

MODULE GammaF5
  USE DerivedTypes
  USE SpecFunMesh
  USE GammaAssymp
  IMPLICIT NONE
  INTEGER :: IF5
  REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F5_0,F5_1,F5_2,F5_3,F5_4
  INCLUDE 'F5.Inc'
END MODULE GammaF5

MODULE GammaF6
  USE DerivedTypes
  USE SpecFunMesh
  USE GammaAssymp
  IMPLICIT NONE
  INTEGER :: IF6
  REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F6_0,F6_1,F6_2,F6_3,F6_4
  INCLUDE 'F6.Inc'
END MODULE GammaF6

MODULE GammaF7
  USE DerivedTypes
  USE SpecFunMesh
  USE GammaAssymp
  IMPLICIT NONE
  INTEGER :: IF7
  REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F7_0,F7_1,F7_2,F7_3,F7_4
  INCLUDE 'F7.Inc'
END MODULE GammaF7

MODULE GammaF8
  USE DerivedTypes
  USE SpecFunMesh
  USE GammaAssymp
  IMPLICIT NONE
  INTEGER :: IF8
  REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F8_0,F8_1,F8_2,F8_3,F8_4
  INCLUDE 'F8.Inc'
END MODULE GammaF8

MODULE GammaF9
  USE DerivedTypes
  USE SpecFunMesh
  USE GammaAssymp
  IMPLICIT NONE
  INTEGER :: IF9
  REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F9_0,F9_1,F9_2,F9_3,F9_4
  INCLUDE 'F9.Inc'
END MODULE GammaF9

MODULE GammaF10
  USE DerivedTypes
  USE SpecFunMesh
  USE GammaAssymp
  IMPLICIT NONE
  INTEGER :: IF10
  REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F10_0,F10_1,F10_2,F10_3,F10_4
  INCLUDE 'F10.Inc'
END MODULE GammaF10

MODULE GammaF11
  USE DerivedTypes
  USE SpecFunMesh
  USE GammaAssymp
  IMPLICIT NONE
  INTEGER :: IF11
  REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F11_0,F11_1,F11_2,F11_3,F11_4
  INCLUDE 'F11.Inc'
END MODULE GammaF11

MODULE GammaF12
  USE DerivedTypes
  USE SpecFunMesh
  USE GammaAssymp
  IMPLICIT NONE
  INTEGER :: IF12
  REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F12_0,F12_1,F12_2,F12_3,F12_4
  INCLUDE 'F12.Inc'
END MODULE GammaF12

MODULE GammaF13
  USE DerivedTypes
  USE SpecFunMesh
  USE GammaAssymp
  IMPLICIT NONE
  INTEGER :: IF13
  REAL(DOUBLE),DIMENSION(0:Gamma_Mesh) :: F13_0,F13_1,F13_2,F13_3,F13_4
  INCLUDE 'F13.Inc'
END MODULE GammaF13


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
END MODULE InvExp

MODULE ErfFunk
  USE DerivedTypes
  USE ProcessControl
  USE Parse
  IMPLICIT NONE
  INCLUDE "ErrorFunction.Inc"
CONTAINS
  !========================================================================================
  !     Compute the erf function
  !========================================================================================
  FUNCTION Erf(W)
    REAL(DOUBLE) :: Erf,W,X,Sgn
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
  END FUNCTION Erf
  !========================================================================================
  !     The complimentary error function
  !========================================================================================
  FUNCTION ERFC(X)
    REAL(DOUBLE) :: ERFC,X
    ERFC=1.0D0-Erf(X)
  END FUNCTION ERFC
  !========================================================================================
  !     The Funtion ==> -ProductLog[-1,-x] from [0,1/E]
  !========================================================================================
  FUNCTION ProductLog1(X)
    INTEGER                  :: I
    REAL(DOUBLE)             :: X,ProductLog1,Error
    REAL(DOUBLE),PARAMETER   :: Xmax = 0.3678794411714423216D0
    !
    IF(X > Xmax)CALL Halt('In ProductLog1[-1,-x] :: |x| > 1/E')
    ProductLog1 = One
    DO I=1,100
       Error       = ProductLog1
       ProductLog1 = -LOG(X/ProductLog1)
       Error = ABS((ProductLog1-Error)/ProductLog1)
       IF(Error<1D-12)GOTO 99
    ENDDO
    CALL Halt('Failed to converge ProductLog1: X = '  &
         //TRIM(DblToChar(X))//' ProductLog1 = ' &
         //TRIM(DblToChar(ProductLog1))//'.')
99  RETURN
  END FUNCTION ProductLog1
END MODULE ErfFunk
!-------------------------------------------------------------------------------------
! Incomplete Associated Gamma Functions I. Shavitt ""The Gaussian function in
! calculations of statistical mechanics and quantum mechanics"
! Meth. Comp. Phys. 2, p.1 (1963)
!-------------------------------------------------------------------------------------
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
  USE GammaF12
  USE GammaF13
  USE ErfFunk
  USE ProcessControl
  IMPLICIT NONE
CONTAINS
  FUNCTION GammaF(M,T)
    INTEGER,      INTENT(IN) :: M
    REAL(DOUBLE), INTENT(IN) :: T
    INTEGER                  :: J
    REAL(DOUBLE)             :: GammaF
    !-------------------------------------------------------------------------
    IF(T<Gamma_Switch)THEN
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
       CASE(12)
          GammaF=F12_0(J)+T*(F12_1(J)+T*(F12_2(J)+T*(F12_3(J)+T*F12_4(J))))
       CASE(13)
          GammaF=F13_0(J)+T*(F13_1(J)+T*(F13_2(J)+T*(F13_3(J)+T*F13_4(J))))
       CASE(14:)
          CALL Halt('Only M<=13 has been implimented in GammaF ')
       END SELECT
    ELSE
       GammaF=GammAss(M)*T**(-DBLE(M)-Half)
    ENDIF
  END FUNCTION GammaF

  FUNCTION GammaHalf(L,X)
    INTEGER                    :: L,LS
    REAL(DOUBLE)               :: X,SqrtX,XSUM,GammaHalf
    REAL(DOUBLE),DIMENSION(64) :: Sfac = (/2.00000000000000000D0,1.33333333333333333D0,&
         7.3029674334022148D-1,5.3412580601946498D-1,4.2897258944050779D-1,&
         3.6130260767122454D-1,3.1338173978619282D-1,2.7736675692301273D-1,&
         2.4916970717917522D-1,2.2642030439700098D-1,2.0763707534036204D-1,&
         1.9184081786310919D-1,1.7835560611619894D-1,1.6669836456773905D-1,&
         1.5651386447776487D-1,1.4753458765116469D-1,1.3955494670757578D-1,&
         1.3241416731800685D-1,1.2598459485542503D-1,1.2016350807025758D-1,&
         1.1486726370711862D-1,1.1002702835105538D-1,1.0558561445605698D-1,&
         1.0149509929347097D-1,9.7715008595692035D-2,9.4210913822889116D-2,&
         9.0953336661706548D-2,8.7916884656620846D-2,8.5079562763772593D-2,&
         8.2422220248006543D-2,7.9928102738676689D-2,7.7582486742601370D-2,&
         7.5372379364895518D-2,7.3286270006162048D-2,7.1313923796257465D-2,&
         6.9446208774383175D-2,6.7674950532187841D-2,6.5992809342867331D-2,&
         6.4393175806982839D-2,6.2870081828999622D-2,6.1418124351705448D-2,&
         6.0032399758875509D-2,5.8708447239769698D-2,5.7442199714798083D-2,&
         5.6229941167024131D-2,5.5068269422098354D-2,5.3954063579712824D-2,&
         5.2884455430455434D-2,5.1856804299020200D-2,5.0868674842786776D-2,&
         4.9917817407506679D-2,4.9002150602140566D-2,4.8119745805096224D-2,&
         4.7268813356069430D-2,4.6447690222872336D-2,4.5654828962240815D-2,&
         4.4888787818609763D-2,4.4148221826018521D-2,4.3431874796297674D-2,&
         4.2738572092017457D-2,4.2067214095777623D-2,4.1416770298643776D-2,&
         4.0786273940179798D-2,4.0174817140833709D-2 /)
    IF(L > 64) CALL MondoHalt(-100,"L > 64 in GammaHalf")
    SqrtX      = SQRT(X)
    IF(L == 0) THEN
       GammaHalf = (One-Erf(SqrtX))
    ELSE
       XSUM = SFAC(1)*EXP(-X)
       DO LS = 2,L
          XSUM  = XSUM+(Sfac(LS)*X*EXP(-X/DBLE(LS-1)))**(LS-1)
       ENDDO
       GammaHalf= (One-Erf(SqrtX))+(SqrtX/SqrtPi)*XSUM
    ENDIF
  END FUNCTION GammaHalf
  !
  !
  FUNCTION GammaOne(L,X)
    INTEGER                    :: L,LS
    REAL(DOUBLE)               :: X,GammaOne
    REAL(DOUBLE),DIMENSION(64) :: Sfac = (/1.000000000000000000D0,7.07106781186547524D-1,&
         5.50321208149104447D-1,4.51801001804922416D-1,3.83851949637377488D-1,&
         3.34024188266401231D-1,2.95856661262446280D-1,2.65650069930254085D-1,&
         2.41128504100169337D-1,2.20812521320600886D-1,2.03697567972983717D-1,&
         1.89076949234936271D-1,1.76438688859779138D-1,1.65402572778695283D-1,&
         1.55680192257920798D-1,1.47048729104116938D-1,1.39333258314618537D-1,&
         1.32394499854830953D-1,1.26120154361686513D-1,1.20418654182541255D-1,&
         1.15214577818197263D-1,1.10445232321394289D-1,1.06058070184465669D-1,&
         1.02008711931231563D-1,9.82594146986036343D-2,9.47778735343361583D-2,&
         9.15362739042765808D-2,8.85105359766346556D-2,8.56797068131988230D-2,&
         8.30254677164270479D-2,8.05317320243545763D-2,7.81843145303305836D-2,&
         7.59706580569001366D-2,7.38796059641125900D-2,7.19012118235820398D-2,&
         7.00265793537853832D-2,6.82477271415646394D-2,6.65574737794576092D-2,&
         6.49493399083641569D-2,6.34174643290225174D-2,6.19565318774671916D-2,&
         6.05617111816940929D-2,5.92286007537456755D-2,5.79531821419947414D-2,&
         5.67317790867945570D-2,5.55610217998274493D-2,5.44378156318983940D-2,&
         5.33593135121798543D-2,5.23228916391817847D-2,5.13261279840590129D-2,&
         5.03667832334909196D-2,4.94427838548302931D-2,4.85522070125600978D-2,&
         4.76932671039560510D-2,4.68643037145489063D-2,4.60637708215811051D-2,&
         4.52902270970256948D-2,4.45423271815889401D-2,4.38188138180314242D-2,&
         4.31185107465918583D-2,4.24403162776759877D-2,4.17831974676047934D-2,&
         4.11461848323713144D-2,4.05283675422583969D-2 /)
    IF(L > 64) CALL MondoHalt(-100,"L > 64 in GammaOne")
    GammaOne = EXP(-X)
    DO LS = 1,L
       GammaOne = GammaOne+(Sfac(LS)*X*EXP(-X/DBLE(LS)))**LS
    ENDDO
  END FUNCTION GammaOne
END MODULE GammaFunctions

MODULE SpecFun
  USE GammaFunctions
  USE ErfFunk
  USE InvExp
END MODULE SpecFun
