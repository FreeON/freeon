      SUBROUTINE ConvertG12()
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "Mesh.Inc"
      REAL*8 F12_0(0:Mesh),F12_1(0:Mesh),F12_2(0:Mesh),F12_3(0:Mesh)
      REAL*8 Table(0:3,0:Mesh)
      INCLUDE "Gamma_12.Inc"
      INCLUDE "Gamma_Asymptotics.Inc"
     
      DO I=0,Mesh
         Table(0,I)=F12_0(I)
         Table(1,I)=F12_1(I)
         Table(2,I)=F12_2(I)
         Table(3,I)=F12_3(I)
      ENDDO

      CALL GDump(12,4*(Mesh+1),F12Asymp,Table)

      RETURN
      END
