      SUBROUTINE ConvertG11()
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER Mesh,I
      INCLUDE "Mesh.Inc"
      REAL*8 F11_0(0:Mesh),F11_1(0:Mesh),F11_2(0:Mesh),F11_3(0:Mesh)
      REAL*8 Table(0:3,0:Mesh)
      INCLUDE "Gamma_11.Inc"
      INCLUDE "Gamma_Asymptotics.Inc"
     
      DO I=0,Mesh
         Table(0,I)=F11_0(I)
         Table(1,I)=F11_1(I)
         Table(2,I)=F11_2(I)
         Table(3,I)=F11_3(I)
      ENDDO

      CALL GDump(11,4*(Mesh+1),Table)

      RETURN
      END
