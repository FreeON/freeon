!    SPECIAL (AND NOT SO-SPECIAL) FUNCTIONS
!    Authors: Matt Challacombe and C.J. Tymczak
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
   USE ErfFunk
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
           GammaHalf = (One-ERF(SqrtX))
        ELSE
           XSUM = SFAC(1)*EXP(-X)
           DO LS = 2,L
              XSUM  = XSUM+(Sfac(LS)*X*EXP(-X/DBLE(LS-1)))**(LS-1)
           ENDDO
           GammaHalf= (One-ERF(SqrtX))+(SqrtX/SqrtPi)*XSUM
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
END MODULE

MODULE SpecFun
   USE GammaFunctions
   USE ErfFunk
   USE InvExp
END MODULE
