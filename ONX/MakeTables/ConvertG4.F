      SUBROUTINE ConvertG4()
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "Mesh.Inc"
      REAL*8 F4_0(0:Mesh),F4_1(0:Mesh),F4_2(0:Mesh),F4_3(0:Mesh)
      REAL*8 Table(0:3,0:Mesh)
      INCLUDE "Gamma_4.Inc"
      INCLUDE "Gamma_Asymptotics.Inc"
     
      DO I=0,Mesh
         Table(0,I)=F4_0(I)
         Table(1,I)=F4_1(I)
         Table(2,I)=F4_2(I)
         Table(3,I)=F4_3(I)
      ENDDO

      CALL GDump(4,4*(Mesh+1),Table)

      RETURN
      END
