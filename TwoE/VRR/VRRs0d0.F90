   SUBROUTINE VRRs0d0(LB,LK,VRR0,VRR1) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2E*VRR0(1,1)
      V(2)=r1x2E*ZxZpE*VRR1(1,1)
      V(3)=-V(2)
      VRR0(1,5)=V(1)+V(3)+QCx*VRR0(1,2)+WQx*VRR1(1,2)
      VRR0(1,6)=QCx*VRR0(1,3)+WQx*VRR1(1,3)
      VRR0(1,7)=V(1)+V(3)+QCy*VRR0(1,3)+WQy*VRR1(1,3)
      VRR0(1,8)=QCx*VRR0(1,4)+WQx*VRR1(1,4)
      VRR0(1,9)=QCy*VRR0(1,4)+WQy*VRR1(1,4)
      VRR0(1,10)=V(1)+V(3)+QCz*VRR0(1,4)+WQz*VRR1(1,4)
END SUBROUTINE VRRs0d0