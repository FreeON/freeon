   SUBROUTINE KetHRR21(LB,HRR) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB
      REAL(DOUBLE) :: HRR(1:LB,5,1)
      !=|1,1)
      HRR(1:LB,1,1)=HRR(1:LB,5,1)
END SUBROUTINE KetHRR21