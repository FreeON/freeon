      SUBROUTINE ConvertG7()
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER Mesh,I
      INCLUDE "Mesh.Inc"
      REAL*8 F7_0(0:Mesh),F7_1(0:Mesh),F7_2(0:Mesh),F7_3(0:Mesh)
      REAL*8 Table(0:3,0:Mesh)
      INCLUDE "Gamma_7.Inc"
      INCLUDE "Gamma_Asymptotics.Inc"
     
      DO I=0,Mesh
         Table(0,I)=F7_0(I)
         Table(1,I)=F7_1(I)
         Table(2,I)=F7_2(I)
         Table(3,I)=F7_3(I)
      ENDDO

      CALL GDump(7,4*(Mesh+1),Table)

      RETURN
      END
