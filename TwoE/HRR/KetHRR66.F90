   SUBROUTINE KetHRR66(LB,HRR) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB
      REAL(DOUBLE) :: HRR(1:LB,35,10)
      !=|5,5)
      HRR(1:LB,5,5)=CDx*(CDx*HRR(1:LB,5,1)+  & 
                        2.D0*HRR(1:LB,11,1))+  & 
                        HRR(1:LB,21,1)
      !=|6,5)
      HRR(1:LB,6,5)=CDx*(CDx*HRR(1:LB,6,1)+  & 
                        2.D0*HRR(1:LB,12,1))+  & 
                        HRR(1:LB,22,1)
      !=|7,5)
      HRR(1:LB,7,5)=CDx*(CDx*HRR(1:LB,7,1)+  & 
                        2.D0*HRR(1:LB,13,1))+  & 
                        HRR(1:LB,23,1)
      !=|8,5)
      HRR(1:LB,8,5)=CDx*(CDx*HRR(1:LB,8,1)+  & 
                        2.D0*HRR(1:LB,15,1))+  & 
                        HRR(1:LB,26,1)
      !=|9,5)
      HRR(1:LB,9,5)=CDx*(CDx*HRR(1:LB,9,1)+  & 
                        2.D0*HRR(1:LB,16,1))+  & 
                        HRR(1:LB,27,1)
      !=|10,5)
      HRR(1:LB,10,5)=CDx*(CDx*HRR(1:LB,10,1)+  & 
                        2.D0*HRR(1:LB,18,1))+  & 
                        HRR(1:LB,30,1)
      !=|5,6)
      HRR(1:LB,5,6)=CDy*HRR(1:LB,11,1)+  & 
                        CDx*(CDy*HRR(1:LB,5,1)+  & 
                        HRR(1:LB,12,1))+  & 
                        HRR(1:LB,22,1)
      !=|6,6)
      HRR(1:LB,6,6)=CDy*HRR(1:LB,12,1)+  & 
                        CDx*(CDy*HRR(1:LB,6,1)+  & 
                        HRR(1:LB,13,1))+  & 
                        HRR(1:LB,23,1)
      !=|7,6)
      HRR(1:LB,7,6)=CDy*HRR(1:LB,13,1)+  & 
                        CDx*(CDy*HRR(1:LB,7,1)+  & 
                        HRR(1:LB,14,1))+  & 
                        HRR(1:LB,24,1)
      !=|8,6)
      HRR(1:LB,8,6)=CDy*HRR(1:LB,15,1)+  & 
                        CDx*(CDy*HRR(1:LB,8,1)+  & 
                        HRR(1:LB,16,1))+  & 
                        HRR(1:LB,27,1)
      !=|9,6)
      HRR(1:LB,9,6)=CDy*HRR(1:LB,16,1)+  & 
                        CDx*(CDy*HRR(1:LB,9,1)+  & 
                        HRR(1:LB,17,1))+  & 
                        HRR(1:LB,28,1)
      !=|10,6)
      HRR(1:LB,10,6)=CDy*HRR(1:LB,18,1)+  & 
                        CDx*(CDy*HRR(1:LB,10,1)+  & 
                        HRR(1:LB,19,1))+  & 
                        HRR(1:LB,31,1)
      !=|5,7)
      HRR(1:LB,5,7)=CDy*(CDy*HRR(1:LB,5,1)+  & 
                        2.D0*HRR(1:LB,12,1))+  & 
                        HRR(1:LB,23,1)
      !=|6,7)
      HRR(1:LB,6,7)=CDy*(CDy*HRR(1:LB,6,1)+  & 
                        2.D0*HRR(1:LB,13,1))+  & 
                        HRR(1:LB,24,1)
      !=|7,7)
      HRR(1:LB,7,7)=CDy*(CDy*HRR(1:LB,7,1)+  & 
                        2.D0*HRR(1:LB,14,1))+  & 
                        HRR(1:LB,25,1)
      !=|8,7)
      HRR(1:LB,8,7)=CDy*(CDy*HRR(1:LB,8,1)+  & 
                        2.D0*HRR(1:LB,16,1))+  & 
                        HRR(1:LB,28,1)
      !=|9,7)
      HRR(1:LB,9,7)=CDy*(CDy*HRR(1:LB,9,1)+  & 
                        2.D0*HRR(1:LB,17,1))+  & 
                        HRR(1:LB,29,1)
      !=|10,7)
      HRR(1:LB,10,7)=CDy*(CDy*HRR(1:LB,10,1)+  & 
                        2.D0*HRR(1:LB,19,1))+  & 
                        HRR(1:LB,32,1)
      !=|5,8)
      HRR(1:LB,5,8)=CDz*HRR(1:LB,11,1)+  & 
                        CDx*(CDz*HRR(1:LB,5,1)+  & 
                        HRR(1:LB,15,1))+  & 
                        HRR(1:LB,26,1)
      !=|6,8)
      HRR(1:LB,6,8)=CDz*HRR(1:LB,12,1)+  & 
                        CDx*(CDz*HRR(1:LB,6,1)+  & 
                        HRR(1:LB,16,1))+  & 
                        HRR(1:LB,27,1)
      !=|7,8)
      HRR(1:LB,7,8)=CDz*HRR(1:LB,13,1)+  & 
                        CDx*(CDz*HRR(1:LB,7,1)+  & 
                        HRR(1:LB,17,1))+  & 
                        HRR(1:LB,28,1)
      !=|8,8)
      HRR(1:LB,8,8)=CDz*HRR(1:LB,15,1)+  & 
                        CDx*(CDz*HRR(1:LB,8,1)+  & 
                        HRR(1:LB,18,1))+  & 
                        HRR(1:LB,30,1)
      !=|9,8)
      HRR(1:LB,9,8)=CDz*HRR(1:LB,16,1)+  & 
                        CDx*(CDz*HRR(1:LB,9,1)+  & 
                        HRR(1:LB,19,1))+  & 
                        HRR(1:LB,31,1)
      !=|10,8)
      HRR(1:LB,10,8)=CDz*HRR(1:LB,18,1)+  & 
                        CDx*(CDz*HRR(1:LB,10,1)+  & 
                        HRR(1:LB,20,1))+  & 
                        HRR(1:LB,33,1)
      !=|5,9)
      HRR(1:LB,5,9)=CDz*HRR(1:LB,12,1)+  & 
                        CDy*(CDz*HRR(1:LB,5,1)+  & 
                        HRR(1:LB,15,1))+  & 
                        HRR(1:LB,27,1)
      !=|6,9)
      HRR(1:LB,6,9)=CDz*HRR(1:LB,13,1)+  & 
                        CDy*(CDz*HRR(1:LB,6,1)+  & 
                        HRR(1:LB,16,1))+  & 
                        HRR(1:LB,28,1)
      !=|7,9)
      HRR(1:LB,7,9)=CDz*HRR(1:LB,14,1)+  & 
                        CDy*(CDz*HRR(1:LB,7,1)+  & 
                        HRR(1:LB,17,1))+  & 
                        HRR(1:LB,29,1)
      !=|8,9)
      HRR(1:LB,8,9)=CDz*HRR(1:LB,16,1)+  & 
                        CDy*(CDz*HRR(1:LB,8,1)+  & 
                        HRR(1:LB,18,1))+  & 
                        HRR(1:LB,31,1)
      !=|9,9)
      HRR(1:LB,9,9)=CDz*HRR(1:LB,17,1)+  & 
                        CDy*(CDz*HRR(1:LB,9,1)+  & 
                        HRR(1:LB,19,1))+  & 
                        HRR(1:LB,32,1)
      !=|10,9)
      HRR(1:LB,10,9)=CDz*HRR(1:LB,19,1)+  & 
                        CDy*(CDz*HRR(1:LB,10,1)+  & 
                        HRR(1:LB,20,1))+  & 
                        HRR(1:LB,34,1)
      !=|5,10)
      HRR(1:LB,5,10)=CDz*(CDz*HRR(1:LB,5,1)+  & 
                        2.D0*HRR(1:LB,15,1))+  & 
                        HRR(1:LB,30,1)
      !=|6,10)
      HRR(1:LB,6,10)=CDz*(CDz*HRR(1:LB,6,1)+  & 
                        2.D0*HRR(1:LB,16,1))+  & 
                        HRR(1:LB,31,1)
      !=|7,10)
      HRR(1:LB,7,10)=CDz*(CDz*HRR(1:LB,7,1)+  & 
                        2.D0*HRR(1:LB,17,1))+  & 
                        HRR(1:LB,32,1)
      !=|8,10)
      HRR(1:LB,8,10)=CDz*(CDz*HRR(1:LB,8,1)+  & 
                        2.D0*HRR(1:LB,18,1))+  & 
                        HRR(1:LB,33,1)
      !=|9,10)
      HRR(1:LB,9,10)=CDz*(CDz*HRR(1:LB,9,1)+  & 
                        2.D0*HRR(1:LB,19,1))+  & 
                        HRR(1:LB,34,1)
      !=|10,10)
      HRR(1:LB,10,10)=CDz*(CDz*HRR(1:LB,10,1)+  & 
                        2.D0*HRR(1:LB,20,1))+  & 
                        HRR(1:LB,35,1)
      !=|2,5)
      HRR(1:LB,2,5)=CDx*(CDx*HRR(1:LB,2,1)+  & 
                        2.D0*HRR(1:LB,5,1))+  & 
                        HRR(1:LB,11,1)
      !=|3,5)
      HRR(1:LB,3,5)=CDx*(CDx*HRR(1:LB,3,1)+  & 
                        2.D0*HRR(1:LB,6,1))+  & 
                        HRR(1:LB,12,1)
      !=|4,5)
      HRR(1:LB,4,5)=CDx*(CDx*HRR(1:LB,4,1)+  & 
                        2.D0*HRR(1:LB,8,1))+  & 
                        HRR(1:LB,15,1)
      !=|2,6)
      HRR(1:LB,2,6)=CDy*HRR(1:LB,5,1)+  & 
                        CDx*(CDy*HRR(1:LB,2,1)+  & 
                        HRR(1:LB,6,1))+  & 
                        HRR(1:LB,12,1)
      !=|3,6)
      HRR(1:LB,3,6)=CDy*HRR(1:LB,6,1)+  & 
                        CDx*(CDy*HRR(1:LB,3,1)+  & 
                        HRR(1:LB,7,1))+  & 
                        HRR(1:LB,13,1)
      !=|4,6)
      HRR(1:LB,4,6)=CDy*HRR(1:LB,8,1)+  & 
                        CDx*(CDy*HRR(1:LB,4,1)+  & 
                        HRR(1:LB,9,1))+  & 
                        HRR(1:LB,16,1)
      !=|2,7)
      HRR(1:LB,2,7)=CDy*(CDy*HRR(1:LB,2,1)+  & 
                        2.D0*HRR(1:LB,6,1))+  & 
                        HRR(1:LB,13,1)
      !=|3,7)
      HRR(1:LB,3,7)=CDy*(CDy*HRR(1:LB,3,1)+  & 
                        2.D0*HRR(1:LB,7,1))+  & 
                        HRR(1:LB,14,1)
      !=|4,7)
      HRR(1:LB,4,7)=CDy*(CDy*HRR(1:LB,4,1)+  & 
                        2.D0*HRR(1:LB,9,1))+  & 
                        HRR(1:LB,17,1)
      !=|2,8)
      HRR(1:LB,2,8)=CDz*HRR(1:LB,5,1)+  & 
                        CDx*(CDz*HRR(1:LB,2,1)+  & 
                        HRR(1:LB,8,1))+  & 
                        HRR(1:LB,15,1)
      !=|3,8)
      HRR(1:LB,3,8)=CDz*HRR(1:LB,6,1)+  & 
                        CDx*(CDz*HRR(1:LB,3,1)+  & 
                        HRR(1:LB,9,1))+  & 
                        HRR(1:LB,16,1)
      !=|4,8)
      HRR(1:LB,4,8)=CDz*HRR(1:LB,8,1)+  & 
                        CDx*(CDz*HRR(1:LB,4,1)+  & 
                        HRR(1:LB,10,1))+  & 
                        HRR(1:LB,18,1)
      !=|2,9)
      HRR(1:LB,2,9)=CDz*HRR(1:LB,6,1)+  & 
                        CDy*(CDz*HRR(1:LB,2,1)+  & 
                        HRR(1:LB,8,1))+  & 
                        HRR(1:LB,16,1)
      !=|3,9)
      HRR(1:LB,3,9)=CDz*HRR(1:LB,7,1)+  & 
                        CDy*(CDz*HRR(1:LB,3,1)+  & 
                        HRR(1:LB,9,1))+  & 
                        HRR(1:LB,17,1)
      !=|4,9)
      HRR(1:LB,4,9)=CDz*HRR(1:LB,9,1)+  & 
                        CDy*(CDz*HRR(1:LB,4,1)+  & 
                        HRR(1:LB,10,1))+  & 
                        HRR(1:LB,19,1)
      !=|2,10)
      HRR(1:LB,2,10)=CDz*(CDz*HRR(1:LB,2,1)+  & 
                        2.D0*HRR(1:LB,8,1))+  & 
                        HRR(1:LB,18,1)
      !=|3,10)
      HRR(1:LB,3,10)=CDz*(CDz*HRR(1:LB,3,1)+  & 
                        2.D0*HRR(1:LB,9,1))+  & 
                        HRR(1:LB,19,1)
      !=|4,10)
      HRR(1:LB,4,10)=CDz*(CDz*HRR(1:LB,4,1)+  & 
                        2.D0*HRR(1:LB,10,1))+  & 
                        HRR(1:LB,20,1)
      !=|1,5)
      HRR(1:LB,1,5)=CDx*(CDx*HRR(1:LB,1,1)+  & 
                        2.D0*HRR(1:LB,2,1))+  & 
                        HRR(1:LB,5,1)
      !=|1,6)
      HRR(1:LB,1,6)=CDy*HRR(1:LB,2,1)+  & 
                        CDx*(CDy*HRR(1:LB,1,1)+  & 
                        HRR(1:LB,3,1))+  & 
                        HRR(1:LB,6,1)
      !=|1,7)
      HRR(1:LB,1,7)=CDy*(CDy*HRR(1:LB,1,1)+  & 
                        2.D0*HRR(1:LB,3,1))+  & 
                        HRR(1:LB,7,1)
      !=|1,8)
      HRR(1:LB,1,8)=CDz*HRR(1:LB,2,1)+  & 
                        CDx*(CDz*HRR(1:LB,1,1)+  & 
                        HRR(1:LB,4,1))+  & 
                        HRR(1:LB,8,1)
      !=|1,9)
      HRR(1:LB,1,9)=CDz*HRR(1:LB,3,1)+  & 
                        CDy*(CDz*HRR(1:LB,1,1)+  & 
                        HRR(1:LB,4,1))+  & 
                        HRR(1:LB,9,1)
      !=|1,10)
      HRR(1:LB,1,10)=CDz*(CDz*HRR(1:LB,1,1)+  & 
                        2.D0*HRR(1:LB,4,1))+  & 
                        HRR(1:LB,10,1)
END SUBROUTINE KetHRR66