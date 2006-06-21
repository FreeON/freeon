    Tmp1       = One/RhoXX(1)
    dRhoExc(0) = -1.861051472698200D0*RhoXX(1)**(4.D0/3.D0)
    dRhoExc(1) =  1.333333333333333D0*dRhoExc(0)*Tmp1
    dRhoExc(2) =  0.333333333333333D0*dRhoExc(1)*Tmp1
    dRhoExc(3) = -0.666666666666666D0*dRhoExc(2)*Tmp1
! 
    dRhoVxc(0) = Half*dRhoExc(1)
    dRhoVxc(1) = Half*dRhoExc(2)
    dRhoVxc(2) = Half*dRhoExc(3)
    dRhoVxc(3) = -0.833333333333333D0*dRhoExc(3)*Tmp1
!
