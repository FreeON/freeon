      SUBROUTINE ConvertG8()
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER Mesh,I
      INCLUDE "Mesh.Inc"
      REAL*8 F8_0(0:Mesh),F8_1(0:Mesh),F8_2(0:Mesh),F8_3(0:Mesh)
      REAL*8 Table(0:3,0:Mesh)
      INCLUDE "Gamma_8.Inc"
      INCLUDE "Gamma_Asymptotics.Inc"
     
      DO I=0,Mesh
         Table(0,I)=F8_0(I)
         Table(1,I)=F8_1(I)
         Table(2,I)=F8_2(I)
         Table(3,I)=F8_3(I)
      ENDDO

      CALL GDump(8,4*(Mesh+1),Table)

      RETURN
      END
