      SUBROUTINE ConvertG10()
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "Mesh.Inc"
      REAL*8 F10_0(0:Mesh),F10_1(0:Mesh),F10_2(0:Mesh),F10_3(0:Mesh)
      REAL*8 Table(0:3,0:Mesh)
      INCLUDE "Gamma_10.Inc"
      INCLUDE "Gamma_Asymptotics.Inc"
     
      DO I=0,Mesh
         Table(0,I)=F10_0(I)
         Table(1,I)=F10_1(I)
         Table(2,I)=F10_2(I)
         Table(3,I)=F10_3(I)
      ENDDO

      CALL GDump(10,4*(Mesh+1),F10Asymp,Table)

      RETURN
      END
