   SUBROUTINE KetHRR61(LB,HRR) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB
      REAL(DOUBLE) :: HRR(1:LB,10,1)
      !=|5,1)
      !=|6,1)
      HRR(1:LB,6,1)=1.73205080756888D0*HRR(1:LB,6,1)
      !=|7,1)
      !=|8,1)
      HRR(1:LB,8,1)=1.73205080756888D0*HRR(1:LB,8,1)
      !=|9,1)
      HRR(1:LB,9,1)=1.73205080756888D0*HRR(1:LB,9,1)
      !=|10,1)
END SUBROUTINE KetHRR61