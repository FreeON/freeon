      SUBROUTINE ConvertG3()
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER Mesh,I
      INCLUDE "Mesh.Inc"
      REAL*8 F3_0(0:Mesh),F3_1(0:Mesh),F3_2(0:Mesh),F3_3(0:Mesh)
      REAL*8 Table(0:3,0:Mesh)
      INCLUDE "Gamma_3.Inc"
      INCLUDE "Gamma_Asymptotics.Inc"
     
      DO I=0,Mesh
         Table(0,I)=F3_0(I)
         Table(1,I)=F3_1(I)
         Table(2,I)=F3_2(I)
         Table(3,I)=F3_3(I)
      ENDDO

      CALL GDump(3,4*(Mesh+1),Table)

      RETURN
      END
