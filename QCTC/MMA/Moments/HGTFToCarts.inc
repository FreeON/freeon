!
!   Parameters
!
    REAL(DOUBLE),PARAMETER           :: N2 = 1.D0/2.D0
    REAL(DOUBLE),PARAMETER           :: N4 = 1.D0/4.D0
    REAL(DOUBLE),PARAMETER           :: N6 = 1.D0/6.D0
    REAL(DOUBLE),PARAMETER           :: N8 = 1.D0/8.D0
!
!   Factors
!
    PiZeta = (Pi/Zeta)**(threehalves)
    LCode  = 100*LQ+LP
!
    SELECT CASE(LCode)
    CASE (0:99)
       QTemp(1) = Coefs(1)*PiZeta + QTemp(1)
!
    CASE (100)
       QTemp(1) = Coefs(1)*PiZeta + QTemp(1)
       QTemp(4) = Coefs(1)*PiZeta*RZ + QTemp(4)
       QTemp(3) = Coefs(1)*PiZeta*RY + QTemp(3)
       QTemp(2) = Coefs(1)*PiZeta*RX + QTemp(2)
    CASE (101:199)
       QTemp(1) = Coefs(1)*PiZeta + QTemp(1)
       QTemp(4) = PiZeta*(Coefs(4) + Coefs(1)*RZ) + QTemp(4)
       QTemp(3) = PiZeta*(Coefs(3) + Coefs(1)*RY) + QTemp(3)
       QTemp(2) = PiZeta*(Coefs(2) + Coefs(1)*RX) + QTemp(2)
!
    CASE (200)
       QTemp(1)  = Coefs(1)*PiZeta + QTemp(1)
       QTemp(4)  = Coefs(1)*PiZeta*RZ + QTemp(4)
       QTemp(10) = QTemp(10) + PiZeta*(N2*Coefs(1)*RZ**2 + (N4*Coefs(1))/Zeta)
       QTemp(3)  = Coefs(1)*PiZeta*RY + QTemp(3)
       QTemp(9)  = Coefs(1)*PiZeta*RY*RZ + QTemp(9)
       QTemp(7)  = QTemp(7) + PiZeta*(N2*Coefs(1)*RY**2 + (N4*Coefs(1))/Zeta)
       QTemp(2)  = Coefs(1)*PiZeta*RX + QTemp(2)
       QTemp(8)  = Coefs(1)*PiZeta*RX*RZ + QTemp(8)
       QTemp(6)  = Coefs(1)*PiZeta*RX*RY + QTemp(6)
       QTemp(5)  = QTemp(5) + PiZeta*(N2*Coefs(1)*RX**2 + (N4*Coefs(1))/Zeta)
    CASE (201)
       QTemp(1)  = Coefs(1)*PiZeta + QTemp(1)
       QTemp(4)  = PiZeta*(Coefs(4) + Coefs(1)*RZ) + QTemp(4)
       QTemp(10) = QTemp(10) + PiZeta*(Coefs(4)*RZ + N2*Coefs(1)*RZ**2 + (N4*Coefs(1))/Zeta)
       QTemp(3)  = PiZeta*(Coefs(3) + Coefs(1)*RY) + QTemp(3)
       QTemp(9)  = PiZeta*(Coefs(4)*RY + Coefs(3)*RZ + Coefs(1)*RY*RZ) + QTemp(9)
       QTemp(7)  = QTemp(7) + PiZeta*(Coefs(3)*RY + N2*Coefs(1)*RY**2 + (N4*Coefs(1))/Zeta)
       QTemp(2)  = PiZeta*(Coefs(2) + Coefs(1)*RX) + QTemp(2)
       QTemp(8)  = PiZeta*(Coefs(4)*RX + Coefs(2)*RZ + Coefs(1)*RX*RZ) + QTemp(8)
       QTemp(6)  = PiZeta*(Coefs(3)*RX + Coefs(2)*RY + Coefs(1)*RX*RY) + QTemp(6)
       QTemp(5)  = QTemp(5) + PiZeta*(Coefs(2)*RX + N2*Coefs(1)*RX**2 + (N4*Coefs(1))/Zeta)
    CASE (202:299)
       QTemp(1)  = Coefs(1)*PiZeta + QTemp(1)
       QTemp(4)  = PiZeta*(Coefs(4) + Coefs(1)*RZ) + QTemp(4)
       QTemp(10) = QTemp(10) + PiZeta*(Coefs(10) + Coefs(4)*RZ + N2*Coefs(1)*RZ**2 + (N4*Coefs(1))/Zeta)
       QTemp(3)  = PiZeta*(Coefs(3) + Coefs(1)*RY) + QTemp(3)
       QTemp(9)  = PiZeta*(Coefs(9) + Coefs(4)*RY + Coefs(3)*RZ + Coefs(1)*RY*RZ) + QTemp(9)
       QTemp(7)  = QTemp(7) + PiZeta*(Coefs(7) + Coefs(3)*RY + N2*Coefs(1)*RY**2 + (N4*Coefs(1))/Zeta)
       QTemp(2)  = PiZeta*(Coefs(2) + Coefs(1)*RX) + QTemp(2)
       QTemp(8)  = PiZeta*(Coefs(8) + Coefs(4)*RX + Coefs(2)*RZ + Coefs(1)*RX*RZ) + QTemp(8)
       QTemp(6)  = PiZeta*(Coefs(6) + Coefs(3)*RX + Coefs(2)*RY + Coefs(1)*RX*RY) + QTemp(6)
       QTemp(5)  = QTemp(5) + PiZeta*(Coefs(5) + Coefs(2)*RX + N2*Coefs(1)*RX**2 + (N4*Coefs(1))/Zeta)
!
!    CASE (300)
!    CASE (301)
!    CASE (302)
!    CASE (303:399)
!
!    CASE (400)
!    CASE (401)
!    CASE (402)
!    CASE (403)
!    CASE (404:499)
!
    CASE DEFAULT
       WRITE(*,*) 'LQ = ',LQ,' Bad indexing in HGTFToCarts' 
       STOP
    END SELECT

    IF(LP == 0) THEN 

    ELSEIF(LP == 1) THEN 

    ELSEIF(LP >= 2) THEN 

    ENDIF 
