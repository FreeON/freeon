   SUBROUTINE VRRs0f0(LB,LK,VRR0,VRR1) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2E*VRR0(1,2)
      V(2)=r1x2E*ZxZpE*VRR1(1,2)
      V(3)=r1x2E*VRR0(1,3)
      V(4)=r1x2E*ZxZpE*VRR1(1,3)
      V(5)=-V(4)
      V(6)=-V(2)
      V(7)=r1x2E*VRR0(1,4)
      V(8)=r1x2E*ZxZpE*VRR1(1,4)
      V(9)=-V(8)
      VRR0(1,11)=2.D0*V(1)-2.D0*V(2)+QCx*VRR0(1,5)+WQx*VRR1(1,5)
      VRR0(1,12)=V(3)+V(5)+QCx*VRR0(1,6)+WQx*VRR1(1,6)
      VRR0(1,13)=V(1)+V(6)+QCy*VRR0(1,6)+WQy*VRR1(1,6)
      VRR0(1,14)=2.D0*V(3)-2.D0*V(4)+QCy*VRR0(1,7)+WQy*VRR1(1,7)
      VRR0(1,15)=V(7)+V(9)+QCx*VRR0(1,8)+WQx*VRR1(1,8)
      VRR0(1,16)=QCx*VRR0(1,9)+WQx*VRR1(1,9)
      VRR0(1,17)=V(7)+V(9)+QCy*VRR0(1,9)+WQy*VRR1(1,9)
      VRR0(1,18)=V(1)+V(6)+QCz*VRR0(1,8)+WQz*VRR1(1,8)
      VRR0(1,19)=V(3)+V(5)+QCz*VRR0(1,9)+WQz*VRR1(1,9)
      VRR0(1,20)=2.D0*V(7)-2.D0*V(8)+QCz*VRR0(1,10)+WQz*VRR1(1,10)
END SUBROUTINE VRRs0f0