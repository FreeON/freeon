!
    REAL(DOUBLE),PARAMETER           :: N05    = 0.5000D0
    REAL(DOUBLE),PARAMETER           :: N025   = 0.2500D0
    REAL(DOUBLE),PARAMETER           :: N0375  = 0.3750D0
    REAL(DOUBLE),PARAMETER           :: N0125  = 0.1250D0
    REAL(DOUBLE),PARAMETER           :: N00625 = 0.0625D0
!
    IF(LQ == 0) THEN 
       RC(0) = PiZeta*Coefs(1)
       RS(0) = Zero
    ELSEIF(LQ == 1) THEN 
       RC(0) = PiZeta*Coefs(1)
       RS(0) = Zero
       RC(1) = PiZeta*Coefs(4)
       RS(1) = Zero
       RC(2) = -N05*PiZeta*Coefs(3)
       RS(2) = -N05*PiZeta*Coefs(2)
    ELSEIF(LQ == 2) THEN 
       RC(0) = PiZeta*Coefs(1)
       RS(0) = Zero
       RC(1) = PiZeta*Coefs(4)
       RS(1) = Zero
       RC(2) = -N05*PiZeta*Coefs(3)
       RS(2) = -N05*PiZeta*Coefs(2)
       RC(3) = PiZeta*(-N05*Coefs(5) - N05*Coefs(7) + Coefs(10))
       RS(3) = Zero
       RC(4) = -N05*PiZeta*Coefs(9)
       RS(4) = -N05*PiZeta*Coefs(8)
       RC(5) = PiZeta*(-N025*Coefs(5) + N025*Coefs(7))
       RS(5) = N025*PiZeta*Coefs(6)
    ELSEIF(LQ == 3) THEN 
       RC(0) = PiZeta*Coefs(1)
       RS(0) = Zero
       RC(1) = PiZeta*Coefs(4)
       RS(1) = Zero
       RC(2) = -N05*PiZeta*Coefs(3)
       RS(2) = -N05*PiZeta*Coefs(2)
       RC(3) = PiZeta*(-N05*Coefs(5) - N05*Coefs(7) + Coefs(10))
       RS(3) = Zero
       RC(4) = -N05*PiZeta*Coefs(9)
       RS(4) = -N05*PiZeta*Coefs(8)
       RC(5) = PiZeta*(-N025*Coefs(5) + N025*Coefs(7))
       RS(5) = N025*PiZeta*Coefs(6)
       RC(6) = PiZeta*(-N05*Coefs(15) - N05*Coefs(17) + Coefs(20))
       RS(6) = Zero
       RC(7) = PiZeta*(N0125*Coefs(12) + N0375*Coefs(14) - N05*Coefs(19))
       RS(7) = PiZeta*(N0375*Coefs(11) + N0125*Coefs(13) - N05*Coefs(18))
       RC(8) = PiZeta*(-N025*Coefs(15) + N025*Coefs(17))
       RS(8) = N025*PiZeta*Coefs(16)
       RC(9) = PiZeta*(N0125*Coefs(12) - N0125*Coefs(14))
       RS(9) = PiZeta*(N0125*Coefs(11) - N0125*Coefs(13))
    ELSEIF(LQ == 4) THEN 
       RC(0) = PiZeta*Coefs(1)
       RS(0) = Zero
       RC(1) = PiZeta*Coefs(4)
       RS(1) = Zero
       RC(2) = -N05*PiZeta*Coefs(3)
       RS(2) = -N05*PiZeta*Coefs(2)
       RC(3) = PiZeta*(-N05*Coefs(5) - N05*Coefs(7) + Coefs(10))
       RS(3) = Zero
       RC(4) = -N05*PiZeta*Coefs(9)
       RS(4) = -N05*PiZeta*Coefs(8)
       RC(5) = PiZeta*(-N025*Coefs(5) + N025*Coefs(7))
       RS(5) = N025*PiZeta*Coefs(6)
       RC(6) = PiZeta*(-N05*Coefs(15) - N05*Coefs(17) + Coefs(20))
       RS(6) = Zero
       RC(7) = PiZeta*(N0125*Coefs(12) + N0375*Coefs(14) - N05*Coefs(19))
       RS(7) = PiZeta*(N0375*Coefs(11) + N0125*Coefs(13) - N05*Coefs(18))
       RC(8) = PiZeta*(-N025*Coefs(15) + N025*Coefs(17))
       RS(8) = N025*PiZeta*Coefs(16)
       RC(9) = PiZeta*(N0125*Coefs(12) - N0125*Coefs(14))
       RS(9) = PiZeta*(N0125*Coefs(11) - N0125*Coefs(13))
       RC(10) = PiZeta*(N0375*Coefs(21) + N0125*Coefs(23) + N0375*Coefs(25) - &
                   N05*Coefs(30) - N05*Coefs(32) + Coefs(35))
       RS(10) = Zero
       RC(11) = PiZeta*(N0125*Coefs(27) + N0375*Coefs(29) - N05*Coefs(34))
       RS(11) = PiZeta*(N0375*Coefs(26) + N0125*Coefs(28) - N05*Coefs(33))
       RC(12) = PiZeta*(N025*Coefs(21) - N025*Coefs(25) - N025*Coefs(30) + N025*Coefs(32))
       RS(12) = PiZeta*(-N0125*Coefs(22) - N0125*Coefs(24) + N025*Coefs(31))
       RC(13) = PiZeta*(N0125*Coefs(27) - N0125*Coefs(29))
       RS(13) = PiZeta*(N0125*Coefs(26) - N0125*Coefs(28))
       RC(14) = PiZeta*(N00625*Coefs(21) - N00625*Coefs(23) + N00625*Coefs(25))
       RS(14) = PiZeta*(-N00625*Coefs(22) + N00625*Coefs(24))
    ELSE 
       WRITE(*,*) 'LQ = ',LQ,' Bad indexing in HGToSP' 
       STOP 
    ENDIF 
