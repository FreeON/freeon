      SUBROUTINE ConvertG9()
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER Mesh,I
      INCLUDE "Mesh.Inc"
      REAL*8 F9_0(0:Mesh),F9_1(0:Mesh),F9_2(0:Mesh),F9_3(0:Mesh)
      REAL*8 Table(0:3,0:Mesh)
      INCLUDE "Gamma_9.Inc"
      INCLUDE "Gamma_Asymptotics.Inc"
     
      DO I=0,Mesh
         Table(0,I)=F9_0(I)
         Table(1,I)=F9_1(I)
         Table(2,I)=F9_2(I)
         Table(3,I)=F9_3(I)
      ENDDO

      CALL GDump(9,4*(Mesh+1),Table)

      RETURN
      END
