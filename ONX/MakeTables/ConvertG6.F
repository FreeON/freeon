      SUBROUTINE ConvertG6()
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "Mesh.Inc"
      REAL*8 F6_0(0:Mesh),F6_1(0:Mesh),F6_2(0:Mesh),F6_3(0:Mesh)
      REAL*8 Table(0:3,0:Mesh)
      INCLUDE "Gamma_6.Inc"
      INCLUDE "Gamma_Asymptotics.Inc"
     
      DO I=0,Mesh
         Table(0,I)=F6_0(I)
         Table(1,I)=F6_1(I)
         Table(2,I)=F6_2(I)
         Table(3,I)=F6_3(I)
      ENDDO

      CALL GDump(6,4*(Mesh+1),F6Asymp,Table)

      RETURN
      END
