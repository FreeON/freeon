   SUBROUTINE KetHRR153(LB,HRR) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB
      REAL(DOUBLE) :: HRR(1:LB,56,4)
      !=|21,2)
      HRR(1:LB,21,2)=CDx*HRR(1:LB,21,1)+  & 
                        HRR(1:LB,36,1)
      !=|22,2)
      HRR(1:LB,22,2)=CDx*HRR(1:LB,22,1)+  & 
                        HRR(1:LB,37,1)
      !=|23,2)
      HRR(1:LB,23,2)=CDx*HRR(1:LB,23,1)+  & 
                        HRR(1:LB,38,1)
      !=|24,2)
      HRR(1:LB,24,2)=CDx*HRR(1:LB,24,1)+  & 
                        HRR(1:LB,39,1)
      !=|25,2)
      HRR(1:LB,25,2)=CDx*HRR(1:LB,25,1)+  & 
                        HRR(1:LB,40,1)
      !=|26,2)
      HRR(1:LB,26,2)=CDx*HRR(1:LB,26,1)+  & 
                        HRR(1:LB,42,1)
      !=|27,2)
      HRR(1:LB,27,2)=CDx*HRR(1:LB,27,1)+  & 
                        HRR(1:LB,43,1)
      !=|28,2)
      HRR(1:LB,28,2)=CDx*HRR(1:LB,28,1)+  & 
                        HRR(1:LB,44,1)
      !=|29,2)
      HRR(1:LB,29,2)=CDx*HRR(1:LB,29,1)+  & 
                        HRR(1:LB,45,1)
      !=|30,2)
      HRR(1:LB,30,2)=CDx*HRR(1:LB,30,1)+  & 
                        HRR(1:LB,47,1)
      !=|31,2)
      HRR(1:LB,31,2)=CDx*HRR(1:LB,31,1)+  & 
                        HRR(1:LB,48,1)
      !=|32,2)
      HRR(1:LB,32,2)=CDx*HRR(1:LB,32,1)+  & 
                        HRR(1:LB,49,1)
      !=|33,2)
      HRR(1:LB,33,2)=CDx*HRR(1:LB,33,1)+  & 
                        HRR(1:LB,51,1)
      !=|34,2)
      HRR(1:LB,34,2)=CDx*HRR(1:LB,34,1)+  & 
                        HRR(1:LB,52,1)
      !=|35,2)
      HRR(1:LB,35,2)=CDx*HRR(1:LB,35,1)+  & 
                        HRR(1:LB,54,1)
      !=|21,3)
      HRR(1:LB,21,3)=CDy*HRR(1:LB,21,1)+  & 
                        HRR(1:LB,37,1)
      !=|22,3)
      HRR(1:LB,22,3)=CDy*HRR(1:LB,22,1)+  & 
                        HRR(1:LB,38,1)
      !=|23,3)
      HRR(1:LB,23,3)=CDy*HRR(1:LB,23,1)+  & 
                        HRR(1:LB,39,1)
      !=|24,3)
      HRR(1:LB,24,3)=CDy*HRR(1:LB,24,1)+  & 
                        HRR(1:LB,40,1)
      !=|25,3)
      HRR(1:LB,25,3)=CDy*HRR(1:LB,25,1)+  & 
                        HRR(1:LB,41,1)
      !=|26,3)
      HRR(1:LB,26,3)=CDy*HRR(1:LB,26,1)+  & 
                        HRR(1:LB,43,1)
      !=|27,3)
      HRR(1:LB,27,3)=CDy*HRR(1:LB,27,1)+  & 
                        HRR(1:LB,44,1)
      !=|28,3)
      HRR(1:LB,28,3)=CDy*HRR(1:LB,28,1)+  & 
                        HRR(1:LB,45,1)
      !=|29,3)
      HRR(1:LB,29,3)=CDy*HRR(1:LB,29,1)+  & 
                        HRR(1:LB,46,1)
      !=|30,3)
      HRR(1:LB,30,3)=CDy*HRR(1:LB,30,1)+  & 
                        HRR(1:LB,48,1)
      !=|31,3)
      HRR(1:LB,31,3)=CDy*HRR(1:LB,31,1)+  & 
                        HRR(1:LB,49,1)
      !=|32,3)
      HRR(1:LB,32,3)=CDy*HRR(1:LB,32,1)+  & 
                        HRR(1:LB,50,1)
      !=|33,3)
      HRR(1:LB,33,3)=CDy*HRR(1:LB,33,1)+  & 
                        HRR(1:LB,52,1)
      !=|34,3)
      HRR(1:LB,34,3)=CDy*HRR(1:LB,34,1)+  & 
                        HRR(1:LB,53,1)
      !=|35,3)
      HRR(1:LB,35,3)=CDy*HRR(1:LB,35,1)+  & 
                        HRR(1:LB,55,1)
      !=|21,4)
      HRR(1:LB,21,4)=CDz*HRR(1:LB,21,1)+  & 
                        HRR(1:LB,42,1)
      !=|22,4)
      HRR(1:LB,22,4)=CDz*HRR(1:LB,22,1)+  & 
                        HRR(1:LB,43,1)
      !=|23,4)
      HRR(1:LB,23,4)=CDz*HRR(1:LB,23,1)+  & 
                        HRR(1:LB,44,1)
      !=|24,4)
      HRR(1:LB,24,4)=CDz*HRR(1:LB,24,1)+  & 
                        HRR(1:LB,45,1)
      !=|25,4)
      HRR(1:LB,25,4)=CDz*HRR(1:LB,25,1)+  & 
                        HRR(1:LB,46,1)
      !=|26,4)
      HRR(1:LB,26,4)=CDz*HRR(1:LB,26,1)+  & 
                        HRR(1:LB,47,1)
      !=|27,4)
      HRR(1:LB,27,4)=CDz*HRR(1:LB,27,1)+  & 
                        HRR(1:LB,48,1)
      !=|28,4)
      HRR(1:LB,28,4)=CDz*HRR(1:LB,28,1)+  & 
                        HRR(1:LB,49,1)
      !=|29,4)
      HRR(1:LB,29,4)=CDz*HRR(1:LB,29,1)+  & 
                        HRR(1:LB,50,1)
      !=|30,4)
      HRR(1:LB,30,4)=CDz*HRR(1:LB,30,1)+  & 
                        HRR(1:LB,51,1)
      !=|31,4)
      HRR(1:LB,31,4)=CDz*HRR(1:LB,31,1)+  & 
                        HRR(1:LB,52,1)
      !=|32,4)
      HRR(1:LB,32,4)=CDz*HRR(1:LB,32,1)+  & 
                        HRR(1:LB,53,1)
      !=|33,4)
      HRR(1:LB,33,4)=CDz*HRR(1:LB,33,1)+  & 
                        HRR(1:LB,54,1)
      !=|34,4)
      HRR(1:LB,34,4)=CDz*HRR(1:LB,34,1)+  & 
                        HRR(1:LB,55,1)
      !=|35,4)
      HRR(1:LB,35,4)=CDz*HRR(1:LB,35,1)+  & 
                        HRR(1:LB,56,1)
      !=|11,2)
      HRR(1:LB,11,2)=CDx*HRR(1:LB,11,1)+  & 
                        HRR(1:LB,21,1)
      !=|12,2)
      HRR(1:LB,12,2)=CDx*HRR(1:LB,12,1)+  & 
                        HRR(1:LB,22,1)
      !=|13,2)
      HRR(1:LB,13,2)=CDx*HRR(1:LB,13,1)+  & 
                        HRR(1:LB,23,1)
      !=|14,2)
      HRR(1:LB,14,2)=CDx*HRR(1:LB,14,1)+  & 
                        HRR(1:LB,24,1)
      !=|15,2)
      HRR(1:LB,15,2)=CDx*HRR(1:LB,15,1)+  & 
                        HRR(1:LB,26,1)
      !=|16,2)
      HRR(1:LB,16,2)=CDx*HRR(1:LB,16,1)+  & 
                        HRR(1:LB,27,1)
      !=|17,2)
      HRR(1:LB,17,2)=CDx*HRR(1:LB,17,1)+  & 
                        HRR(1:LB,28,1)
      !=|18,2)
      HRR(1:LB,18,2)=CDx*HRR(1:LB,18,1)+  & 
                        HRR(1:LB,30,1)
      !=|19,2)
      HRR(1:LB,19,2)=CDx*HRR(1:LB,19,1)+  & 
                        HRR(1:LB,31,1)
      !=|20,2)
      HRR(1:LB,20,2)=CDx*HRR(1:LB,20,1)+  & 
                        HRR(1:LB,33,1)
      !=|11,3)
      HRR(1:LB,11,3)=CDy*HRR(1:LB,11,1)+  & 
                        HRR(1:LB,22,1)
      !=|12,3)
      HRR(1:LB,12,3)=CDy*HRR(1:LB,12,1)+  & 
                        HRR(1:LB,23,1)
      !=|13,3)
      HRR(1:LB,13,3)=CDy*HRR(1:LB,13,1)+  & 
                        HRR(1:LB,24,1)
      !=|14,3)
      HRR(1:LB,14,3)=CDy*HRR(1:LB,14,1)+  & 
                        HRR(1:LB,25,1)
      !=|15,3)
      HRR(1:LB,15,3)=CDy*HRR(1:LB,15,1)+  & 
                        HRR(1:LB,27,1)
      !=|16,3)
      HRR(1:LB,16,3)=CDy*HRR(1:LB,16,1)+  & 
                        HRR(1:LB,28,1)
      !=|17,3)
      HRR(1:LB,17,3)=CDy*HRR(1:LB,17,1)+  & 
                        HRR(1:LB,29,1)
      !=|18,3)
      HRR(1:LB,18,3)=CDy*HRR(1:LB,18,1)+  & 
                        HRR(1:LB,31,1)
      !=|19,3)
      HRR(1:LB,19,3)=CDy*HRR(1:LB,19,1)+  & 
                        HRR(1:LB,32,1)
      !=|20,3)
      HRR(1:LB,20,3)=CDy*HRR(1:LB,20,1)+  & 
                        HRR(1:LB,34,1)
      !=|11,4)
      HRR(1:LB,11,4)=CDz*HRR(1:LB,11,1)+  & 
                        HRR(1:LB,26,1)
      !=|12,4)
      HRR(1:LB,12,4)=CDz*HRR(1:LB,12,1)+  & 
                        HRR(1:LB,27,1)
      !=|13,4)
      HRR(1:LB,13,4)=CDz*HRR(1:LB,13,1)+  & 
                        HRR(1:LB,28,1)
      !=|14,4)
      HRR(1:LB,14,4)=CDz*HRR(1:LB,14,1)+  & 
                        HRR(1:LB,29,1)
      !=|15,4)
      HRR(1:LB,15,4)=CDz*HRR(1:LB,15,1)+  & 
                        HRR(1:LB,30,1)
      !=|16,4)
      HRR(1:LB,16,4)=CDz*HRR(1:LB,16,1)+  & 
                        HRR(1:LB,31,1)
      !=|17,4)
      HRR(1:LB,17,4)=CDz*HRR(1:LB,17,1)+  & 
                        HRR(1:LB,32,1)
      !=|18,4)
      HRR(1:LB,18,4)=CDz*HRR(1:LB,18,1)+  & 
                        HRR(1:LB,33,1)
      !=|19,4)
      HRR(1:LB,19,4)=CDz*HRR(1:LB,19,1)+  & 
                        HRR(1:LB,34,1)
      !=|20,4)
      HRR(1:LB,20,4)=CDz*HRR(1:LB,20,1)+  & 
                        HRR(1:LB,35,1)
      !=|5,2)
      HRR(1:LB,5,2)=CDx*HRR(1:LB,5,1)+  & 
                        HRR(1:LB,11,1)
      !=|6,2)
      HRR(1:LB,6,2)=CDx*HRR(1:LB,6,1)+  & 
                        HRR(1:LB,12,1)
      !=|7,2)
      HRR(1:LB,7,2)=CDx*HRR(1:LB,7,1)+  & 
                        HRR(1:LB,13,1)
      !=|8,2)
      HRR(1:LB,8,2)=CDx*HRR(1:LB,8,1)+  & 
                        HRR(1:LB,15,1)
      !=|9,2)
      HRR(1:LB,9,2)=CDx*HRR(1:LB,9,1)+  & 
                        HRR(1:LB,16,1)
      !=|10,2)
      HRR(1:LB,10,2)=CDx*HRR(1:LB,10,1)+  & 
                        HRR(1:LB,18,1)
      !=|5,3)
      HRR(1:LB,5,3)=CDy*HRR(1:LB,5,1)+  & 
                        HRR(1:LB,12,1)
      !=|6,3)
      HRR(1:LB,6,3)=CDy*HRR(1:LB,6,1)+  & 
                        HRR(1:LB,13,1)
      !=|7,3)
      HRR(1:LB,7,3)=CDy*HRR(1:LB,7,1)+  & 
                        HRR(1:LB,14,1)
      !=|8,3)
      HRR(1:LB,8,3)=CDy*HRR(1:LB,8,1)+  & 
                        HRR(1:LB,16,1)
      !=|9,3)
      HRR(1:LB,9,3)=CDy*HRR(1:LB,9,1)+  & 
                        HRR(1:LB,17,1)
      !=|10,3)
      HRR(1:LB,10,3)=CDy*HRR(1:LB,10,1)+  & 
                        HRR(1:LB,19,1)
      !=|5,4)
      HRR(1:LB,5,4)=CDz*HRR(1:LB,5,1)+  & 
                        HRR(1:LB,15,1)
      !=|6,4)
      HRR(1:LB,6,4)=CDz*HRR(1:LB,6,1)+  & 
                        HRR(1:LB,16,1)
      !=|7,4)
      HRR(1:LB,7,4)=CDz*HRR(1:LB,7,1)+  & 
                        HRR(1:LB,17,1)
      !=|8,4)
      HRR(1:LB,8,4)=CDz*HRR(1:LB,8,1)+  & 
                        HRR(1:LB,18,1)
      !=|9,4)
      HRR(1:LB,9,4)=CDz*HRR(1:LB,9,1)+  & 
                        HRR(1:LB,19,1)
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
END SUBROUTINE KetHRR153