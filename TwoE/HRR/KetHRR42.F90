   SUBROUTINE KetHRR42(LB,HRR) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER :: LB
      REAL(DOUBLE) :: HRR(1:LB,45,10)
      !=|5,2)
      HRR(1:LB,5,2)=CDx*HRR(1:LB,5,1)+  & 
                        HRR(1:LB,11,1)
      !=|6,2)
      HRR(1:LB,6,2)=1.73205080756888D0*CDx*HRR(1:LB,6,1)+  & 
                        1.73205080756888D0*HRR(1:LB,12,1)
      !=|7,2)
      HRR(1:LB,7,2)=CDx*HRR(1:LB,7,1)+  & 
                        HRR(1:LB,13,1)
      !=|8,2)
      HRR(1:LB,8,2)=1.73205080756888D0*CDx*HRR(1:LB,8,1)+  & 
                        1.73205080756888D0*HRR(1:LB,15,1)
      !=|9,2)
      HRR(1:LB,9,2)=1.73205080756888D0*CDx*HRR(1:LB,9,1)+  & 
                        1.73205080756888D0*HRR(1:LB,16,1)
      !=|10,2)
      HRR(1:LB,10,2)=CDx*HRR(1:LB,10,1)+  & 
                        HRR(1:LB,18,1)
      !=|5,3)
      HRR(1:LB,5,3)=CDy*HRR(1:LB,5,1)+  & 
                        HRR(1:LB,12,1)
      !=|6,3)
      HRR(1:LB,6,3)=1.73205080756888D0*CDy*HRR(1:LB,6,1)+  & 
                        1.73205080756888D0*HRR(1:LB,13,1)
      !=|7,3)
      HRR(1:LB,7,3)=CDy*HRR(1:LB,7,1)+  & 
                        HRR(1:LB,14,1)
      !=|8,3)
      HRR(1:LB,8,3)=1.73205080756888D0*CDy*HRR(1:LB,8,1)+  & 
                        1.73205080756888D0*HRR(1:LB,16,1)
      !=|9,3)
      HRR(1:LB,9,3)=1.73205080756888D0*CDy*HRR(1:LB,9,1)+  & 
                        1.73205080756888D0*HRR(1:LB,17,1)
      !=|10,3)
      HRR(1:LB,10,3)=CDy*HRR(1:LB,10,1)+  & 
                        HRR(1:LB,19,1)
      !=|5,4)
      HRR(1:LB,5,4)=CDz*HRR(1:LB,5,1)+  & 
                        HRR(1:LB,15,1)
      !=|6,4)
      HRR(1:LB,6,4)=1.73205080756888D0*CDz*HRR(1:LB,6,1)+  & 
                        1.73205080756888D0*HRR(1:LB,16,1)
      !=|7,4)
      HRR(1:LB,7,4)=CDz*HRR(1:LB,7,1)+  & 
                        HRR(1:LB,17,1)
      !=|8,4)
      HRR(1:LB,8,4)=1.73205080756888D0*CDz*HRR(1:LB,8,1)+  & 
                        1.73205080756888D0*HRR(1:LB,18,1)
      !=|9,4)
      HRR(1:LB,9,4)=1.73205080756888D0*CDz*HRR(1:LB,9,1)+  & 
                        1.73205080756888D0*HRR(1:LB,19,1)
      !=|10,4)
      HRR(1:LB,10,4)=CDz*HRR(1:LB,10,1)+  & 
                        HRR(1:LB,20,1)
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
      HRR(1:LB,1,2)=CDx*HRR(1:LB,1,1)+  & 
                        HRR(1:LB,2,1)
      !=|1,3)
      HRR(1:LB,1,3)=CDy*HRR(1:LB,1,1)+  & 
                        HRR(1:LB,3,1)
      !=|1,4)
      HRR(1:LB,1,4)=CDz*HRR(1:LB,1,1)+  & 
                        HRR(1:LB,4,1)
      !=|5,1)
      HRR(1:LB,5,1)=HRR(1:LB,25,1)
      !=|6,1)
      HRR(1:LB,6,1)=1.73205080756888D0*HRR(1:LB,26,1)
      !=|7,1)
      HRR(1:LB,7,1)=HRR(1:LB,27,1)
      !=|8,1)
      HRR(1:LB,8,1)=1.73205080756888D0*HRR(1:LB,28,1)
      !=|9,1)
      HRR(1:LB,9,1)=1.73205080756888D0*HRR(1:LB,29,1)
      !=|10,1)
      HRR(1:LB,10,1)=HRR(1:LB,30,1)
      !=|2,1)
      HRR(1:LB,2,1)=HRR(1:LB,22,1)
      !=|3,1)
      HRR(1:LB,3,1)=HRR(1:LB,23,1)
      !=|4,1)
      HRR(1:LB,4,1)=HRR(1:LB,24,1)
      !=|1,1)
      HRR(1:LB,1,1)=HRR(1:LB,21,1)
END SUBROUTINE KetHRR42