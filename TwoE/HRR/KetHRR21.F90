   SUBROUTINE KetHRR21(LB,HRR) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER :: LB
      REAL(DOUBLE) :: HRR(1:LB,11,4)
      !=|1,1)
      HRR(1:LB,1,1)=HRR(1:LB,5,1)
END SUBROUTINE KetHRR21