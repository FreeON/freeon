      SUBROUTINE ConvertG1()
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "Mesh.Inc"
      REAL*8 F1_0(0:Mesh),F1_1(0:Mesh),F1_2(0:Mesh),F1_3(0:Mesh)
      REAL*8 Table(0:3,0:Mesh)
      INCLUDE "Gamma_1.Inc"
      INCLUDE "Gamma_Asymptotics.Inc"
     
      DO I=0,Mesh
         Table(0,I)=F1_0(I)
         Table(1,I)=F1_1(I)
         Table(2,I)=F1_2(I)
         Table(3,I)=F1_3(I)
      ENDDO

      CALL GDump(1,4*(Mesh+1),Table)

      RETURN
      END
