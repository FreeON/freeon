    IF(LP == 0) THEN 
        QCarts(0,0,0) = Coefs(1)*PiZeta
        QCarts(0,0,1) = -Coefs(1)*PiZeta*RZ
        QCarts(0,1,0) = -Coefs(1)*PiZeta*RY
        QCarts(1,0,0) = -Coefs(1)*PiZeta*RX
    ELSEIF(LP >= 1) THEN 
        QCarts(0,0,0) = Coefs(1)*PiZeta
        QCarts(0,0,1) = PiZeta*(-Coefs(4) - Coefs(1)*RZ)
        QCarts(0,1,0) = PiZeta*(-Coefs(3) - Coefs(1)*RY)
        QCarts(1,0,0) = PiZeta*(-Coefs(2) - Coefs(1)*RX)
    ENDIF 
