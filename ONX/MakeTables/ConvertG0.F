      SUBROUTINE ConvertG0()
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER Mesh,I
      INCLUDE "Mesh.Inc"
      REAL*8 F0_0(0:Mesh),F0_1(0:Mesh),F0_2(0:Mesh),F0_3(0:Mesh)
      REAL*8 Table(0:3,0:Mesh)
      INCLUDE "Gamma_0.Inc"
      INCLUDE "Gamma_Asymptotics.Inc"

     
      DO I=0,Mesh
         Table(0,I)=F0_0(I)
         Table(1,I)=F0_1(I)
         Table(2,I)=F0_2(I)
         Table(3,I)=F0_3(I)
      ENDDO

      CALL GDump(0,4*(Mesh+1),Table)

      RETURN
      END
