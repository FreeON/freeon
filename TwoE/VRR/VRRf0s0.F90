   SUBROUTINE VRRf0s0(LB,LK,VRR0,VRR1) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      V(1)=r1x2Z*VRR0(2,1)
      V(2)=ExZpE*r1x2Z*VRR1(2,1)
      V(3)=r1x2Z*VRR0(3,1)
      V(4)=ExZpE*r1x2Z*VRR1(3,1)
      V(5)=-V(4)
      V(6)=-V(2)
      V(7)=r1x2Z*VRR0(4,1)
      V(8)=ExZpE*r1x2Z*VRR1(4,1)
      V(9)=-V(8)
      VRR0(11,1)=2.D0*V(1)-2.D0*V(2)+PAx*VRR0(5,1)+WPx*VRR1(5,1)
      VRR0(12,1)=V(3)+V(5)+PAx*VRR0(6,1)+WPx*VRR1(6,1)
      VRR0(13,1)=V(1)+V(6)+PAy*VRR0(6,1)+WPy*VRR1(6,1)
      VRR0(14,1)=2.D0*V(3)-2.D0*V(4)+PAy*VRR0(7,1)+WPy*VRR1(7,1)
      VRR0(15,1)=V(7)+V(9)+PAx*VRR0(8,1)+WPx*VRR1(8,1)
      VRR0(16,1)=PAx*VRR0(9,1)+WPx*VRR1(9,1)
      VRR0(17,1)=V(7)+V(9)+PAy*VRR0(9,1)+WPy*VRR1(9,1)
      VRR0(18,1)=V(1)+V(6)+PAz*VRR0(8,1)+WPz*VRR1(8,1)
      VRR0(19,1)=V(3)+V(5)+PAz*VRR0(9,1)+WPz*VRR1(9,1)
      VRR0(20,1)=2.D0*V(7)-2.D0*V(8)+PAz*VRR0(10,1)+WPz*VRR1(10,1)
END SUBROUTINE VRRf0s0