   SUBROUTINE KetHRR41(LB,HRR) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER :: LB
      REAL(DOUBLE) :: HRR(1:LB,20,4)
      !=|6,1)
      HRR(1:LB,6,1)=1.73205080756888D0*HRR(1:LB,6,1)
      !=|8,1)
      HRR(1:LB,8,1)=1.73205080756888D0*HRR(1:LB,8,1)
      !=|9,1)
      HRR(1:LB,9,1)=1.73205080756888D0*HRR(1:LB,9,1)
END SUBROUTINE KetHRR41