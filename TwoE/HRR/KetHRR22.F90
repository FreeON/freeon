   SUBROUTINE KetHRR22(LB,HRR) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER :: LB
      REAL(DOUBLE) :: HRR(1:LB,34,10)
      !=|2,2)
      HRR(1:LB,2,2)=CDx*HRR(1:LB,2,1)+  & 
                        HRR(1:LB,5,1)
      !=|3,2)
      HRR(1:LB,3,2)=CDx*HRR(1:LB,3,1)+  & 
                        HRR(1:LB,6,1)
      !=|4,2)
      HRR(1:LB,4,2)=CDx*HRR(1:LB,4,1)+  & 
                        HRR(1:LB,8,1)
      !=|2,3)
      HRR(1:LB,2,3)=CDy*HRR(1:LB,2,1)+  & 
                        HRR(1:LB,6,1)
      !=|3,3)
      HRR(1:LB,3,3)=CDy*HRR(1:LB,3,1)+  & 
                        HRR(1:LB,7,1)
      !=|4,3)
      HRR(1:LB,4,3)=CDy*HRR(1:LB,4,1)+  & 
                        HRR(1:LB,9,1)
      !=|2,4)
      HRR(1:LB,2,4)=CDz*HRR(1:LB,2,1)+  & 
                        HRR(1:LB,8,1)
      !=|3,4)
      HRR(1:LB,3,4)=CDz*HRR(1:LB,3,1)+  & 
                        HRR(1:LB,9,1)
      !=|4,4)
      HRR(1:LB,4,4)=CDz*HRR(1:LB,4,1)+  & 
                        HRR(1:LB,10,1)
      !=|1,2)
      HRR(1:LB,1,2)=CDx*HRR(1:LB,15,1)+  & 
                        HRR(1:LB,16,1)
      !=|1,3)
      HRR(1:LB,1,3)=CDy*HRR(1:LB,15,1)+  & 
                        HRR(1:LB,17,1)
      !=|1,4)
      HRR(1:LB,1,4)=CDz*HRR(1:LB,15,1)+  & 
                        HRR(1:LB,18,1)
      !=|2,1)
      HRR(1:LB,2,1)=HRR(1:LB,12,1)
      !=|3,1)
      HRR(1:LB,3,1)=HRR(1:LB,13,1)
      !=|4,1)
      HRR(1:LB,4,1)=HRR(1:LB,14,1)
      !=|1,1)
      HRR(1:LB,1,1)=HRR(1:LB,11,1)
END SUBROUTINE KetHRR22