      SUBROUTINE ConvertE()
      IMPLICIT REAL*8 (a-h,o-z)
      INTEGER Mesh,I
      INCLUDE "Mesh.Inc"
      REAL*8 EX_0(0:Mesh),EX_1(0:Mesh),EX_2(0:Mesh),EX_3(0:Mesh)
      REAL*8 Table(0:3,0:Mesh)
      INCLUDE "Exponential.Inc"

      DO I=0,Mesh
         Table(0,I)=EX_0(I)
         Table(1,I)=EX_1(I)
         Table(2,I)=EX_2(I)
         Table(3,I)=EX_3(I)
      ENDDO

      CALL EDump(4*(Mesh+1),Table)

      RETURN
      END
