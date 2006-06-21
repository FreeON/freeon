    Tmp1       = One/RhoXX(1)
    dRhoExc(0) = -1.861051472698200D0*RhoXX(1)**(4.D0/3.D0)
    dRhoExc(1) =  1.333333333333333D0*dRhoExc(0)*Tmp1
    dRhoExc(2) =  0.333333333333333D0*dRhoExc(1)*Tmp1
    dRhoExc(3) = -0.666666666666666D0*dRhoExc(2)*Tmp1
    dRhoExc(4) = -1.666666666666666D0*dRhoExc(3)*Tmp1
    dRhoExc(5) = -2.666666666666666D0*dRhoExc(4)*Tmp1
    dRhoExc(6) = -3.666666666666666D0*dRhoExc(5)*Tmp1
    dRhoExc(7) = -4.666666666666666D0*dRhoExc(6)*Tmp1
    dRhoExc(8) = -5.666666666666666D0*dRhoExc(7)*Tmp1
    dRhoExc(9) = -6.666666666666666D0*dRhoExc(8)*Tmp1
! 
    dRhoVxc(0) = Half*dRhoExc(1)
    dRhoVxc(1) = Half*dRhoExc(2)
    dRhoVxc(2) = Half*dRhoExc(3)
    dRhoVxc(3) = Half*dRhoExc(4)
    dRhoVxc(4) = Half*dRhoExc(5)
    dRhoVxc(5) = Half*dRhoExc(6)
    dRhoVxc(6) = Half*dRhoExc(6)
    dRhoVxc(7) = Half*dRhoExc(6)
    dRhoVxc(8) = Half*dRhoExc(6)
    dRhoVxc(9) = -3.833333333333333D0*dRhoExc(8)*Tmp1
!
