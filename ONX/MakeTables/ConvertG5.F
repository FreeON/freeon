      SUBROUTINE ConvertG5()
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER Mesh,I
      INCLUDE "Mesh.Inc"
      REAL*8 F5_0(0:Mesh),F5_1(0:Mesh),F5_2(0:Mesh),F5_3(0:Mesh)
      REAL*8 Table(0:3,0:Mesh)
      INCLUDE "Gamma_5.Inc"
      INCLUDE "Gamma_Asymptotics.Inc"
     
      DO I=0,Mesh
         Table(0,I)=F5_0(I)
         Table(1,I)=F5_1(I)
         Table(2,I)=F5_2(I)
         Table(3,I)=F5_3(I)
      ENDDO

      CALL GDump(5,4*(Mesh+1),Table)

      RETURN
      END
