    IF(LP == 0) THEN 
        QCarts(0,0,0) = Coefs(1)*PiZeta
        QCarts(0,0,1) = -Coefs(1)*PiZeta*RZ
        QCarts(0,0,2) = PiZeta*(5.d-1*Coefs(1)*RZ**2 + (2.5d-1*Coefs(1))/Zeta)
        QCarts(0,1,0) = -Coefs(1)*PiZeta*RY
        QCarts(0,1,1) = Coefs(1)*PiZeta*RY*RZ
        QCarts(0,2,0) = PiZeta*(5.d-1*Coefs(1)*RY**2 + (2.5d-1*Coefs(1))/Zeta)
        QCarts(1,0,0) = -Coefs(1)*PiZeta*RX
        QCarts(1,0,1) = Coefs(1)*PiZeta*RX*RZ
        QCarts(1,1,0) = Coefs(1)*PiZeta*RX*RY
        QCarts(2,0,0) = PiZeta*(5.d-1*Coefs(1)*RX**2 + (2.5d-1*Coefs(1))/Zeta)
    ELSEIF(LP == 1) THEN 
        QCarts(0,0,0) = Coefs(1)*PiZeta
        QCarts(0,0,1) = PiZeta*(-Coefs(4) - Coefs(1)*RZ)
        QCarts(0,0,2) = PiZeta*(Coefs(4)*RZ + 5.d-1*Coefs(1)*RZ**2 + (2.5d-1*Coefs(1))/Zeta)
        QCarts(0,1,0) = PiZeta*(-Coefs(3) - Coefs(1)*RY)
        QCarts(0,1,1) = PiZeta*(Coefs(9) + Coefs(4)*RY + Coefs(3)*RZ + Coefs(1)*RY*RZ)
        QCarts(0,2,0) = PiZeta*(Coefs(3)*RY + 5.d-1*Coefs(1)*RY**2 + (2.5d-1*Coefs(1))/Zeta)
        QCarts(1,0,0) = PiZeta*(-Coefs(2) - Coefs(1)*RX)
        QCarts(1,0,1) = PiZeta*(Coefs(8) + Coefs(4)*RX + Coefs(2)*RZ + Coefs(1)*RX*RZ)
        QCarts(1,1,0) = PiZeta*(Coefs(6) + Coefs(3)*RX + Coefs(2)*RY + Coefs(1)*RX*RY)
        QCarts(2,0,0) = PiZeta*(Coefs(2)*RX + 5.d-1*Coefs(1)*RX**2 + (2.5d-1*Coefs(1))/Zeta)
    ELSEIF(LP >= 2) THEN 
        QCarts(0,0,0) = Coefs(1)*PiZeta
        QCarts(0,0,1) = PiZeta*(-Coefs(4) - Coefs(1)*RZ)
        QCarts(0,0,2) = PiZeta*(Coefs(10) + Coefs(4)*RZ + 5.d-1*Coefs(1)*RZ**2 + (2.5d-1*Coefs(1))/Zeta)
        QCarts(0,1,0) = PiZeta*(-Coefs(3) - Coefs(1)*RY)
        QCarts(0,1,1) = PiZeta*(Coefs(9) + Coefs(4)*RY + Coefs(3)*RZ + Coefs(1)*RY*RZ)
        QCarts(0,2,0) = PiZeta*(Coefs(7) + Coefs(3)*RY + 5.d-1*Coefs(1)*RY**2 + (2.5d-1*Coefs(1))/Zeta)
        QCarts(1,0,0) = PiZeta*(-Coefs(2) - Coefs(1)*RX)
        QCarts(1,0,1) = PiZeta*(Coefs(8) + Coefs(4)*RX + Coefs(2)*RZ + Coefs(1)*RX*RZ)
        QCarts(1,1,0) = PiZeta*(Coefs(6) + Coefs(3)*RX + Coefs(2)*RY + Coefs(1)*RX*RY)
        QCarts(2,0,0) = PiZeta*(Coefs(5) + Coefs(2)*RX + 5.d-1*Coefs(1)*RX**2 + (2.5d-1*Coefs(1))/Zeta)
    ENDIF 
