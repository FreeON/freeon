      PROGRAM Convert
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "Mesh.Inc"
      REAL*8 Table(0:3,0:Mesh)

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

c      CALL GetGammaTable(0,4*(Mesh+1),Table)
c      DO I=0,Mesh
c      WRITE(77,9000) I,Table(0,I),Table(1,I),Table(2,I),Table(3,I)
c      ENDDO
      STOP
 9000 FORMAT(I3,1X,E18.10,1X,E18.10,1X,E18.10,1X,E18.10)
      END
