   SUBROUTINE VRRd0s0(LB,LK,VRR0,VRR1) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2Z*VRR0(1,1)
      V(2)=ExZpE*r1x2Z*VRR1(1,1)
      V(3)=-V(2)
      VRR0(5,1)=V(1)+V(3)+PAx*VRR0(2,1)+WPx*VRR1(2,1)
      VRR0(6,1)=PAx*VRR0(3,1)+WPx*VRR1(3,1)
      VRR0(7,1)=V(1)+V(3)+PAy*VRR0(3,1)+WPy*VRR1(3,1)
      VRR0(8,1)=PAx*VRR0(4,1)+WPx*VRR1(4,1)
      VRR0(9,1)=PAy*VRR0(4,1)+WPy*VRR1(4,1)
      VRR0(10,1)=V(1)+V(3)+PAz*VRR0(4,1)+WPz*VRR1(4,1)
END SUBROUTINE VRRd0s0