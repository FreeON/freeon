      PROGRAM Convert
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER Mesh,I
      INCLUDE "Mesh.Inc"
      INTEGER Ltot
      PARAMETER (Ltot=12)
      REAL*8 Table(0:3,0:Mesh)
      REAL*8 GammaA(0:Ltot)
      INCLUDE "Gamma_Asymptotics.Inc"

      CALL ConvertG0()
      CALL ConvertG1()
      CALL ConvertG2()
      CALL ConvertG3()
      CALL ConvertG4()
      CALL ConvertG5()
      CALL ConvertG6()
      CALL ConvertG7()
      CALL ConvertG8()
      CALL ConvertG9()
      CALL ConvertG10()
      CALL ConvertG11()
      CALL ConvertG12()
      CALL ConvertE()

      GammaA(0) =F0Asymp
      GammaA(1) =F1Asymp
      GammaA(2) =F2Asymp
      GammaA(3) =F3Asymp
      GammaA(4) =F4Asymp
      GammaA(5) =F5Asymp
      GammaA(6) =F6Asymp
      GammaA(7) =F7Asymp
      GammaA(8) =F8Asymp
      GammaA(9) =F9Asymp
      GammaA(10)=F10Asymp
      GammaA(11)=F11Asymp
      GammaA(12)=F12Asymp

      CALL HDump(Mesh,Switch,Grid)
      CALL ADump(Ltot,GammaA)

      STOP
 9000 FORMAT(I3,1X,E18.10,1X,E18.10,1X,E18.10,1X,E18.10)
      END
