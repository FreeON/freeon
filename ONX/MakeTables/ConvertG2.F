      SUBROUTINE ConvertG2()
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "Mesh.Inc"
      REAL*8 F2_0(0:Mesh),F2_1(0:Mesh),F2_2(0:Mesh),F2_3(0:Mesh)
      REAL*8 Table(0:3,0:Mesh)
      INCLUDE "Gamma_2.Inc"
      INCLUDE "Gamma_Asymptotics.Inc"
     
      DO I=0,Mesh
         Table(0,I)=F2_0(I)
         Table(1,I)=F2_1(I)
         Table(2,I)=F2_2(I)
         Table(3,I)=F2_3(I)
      ENDDO

      CALL GDump(2,4*(Mesh+1),F2Asymp,Table)

      RETURN
      END
