!
!   Parameters
!
    REAL(DOUBLE),PARAMETER           :: N05      = 0.500000D0
    REAL(DOUBLE),PARAMETER           :: N0375    = 0.375000D0
    REAL(DOUBLE),PARAMETER           :: N03125   = 0.312500D0
    REAL(DOUBLE),PARAMETER           :: N025     = 0.250000D0
    REAL(DOUBLE),PARAMETER           :: N0234375 = 0.234375D0
    REAL(DOUBLE),PARAMETER           :: N015625  = 0.156250D0
    REAL(DOUBLE),PARAMETER           :: N0125    = 0.125000D0
    REAL(DOUBLE),PARAMETER           :: N009375  = 0.093750D0
    REAL(DOUBLE),PARAMETER           :: N0078125 = 0.078125D0
    REAL(DOUBLE),PARAMETER           :: N00625   = 0.062500D0
    REAL(DOUBLE),PARAMETER           :: N0046875 = 0.046875D0
    REAL(DOUBLE),PARAMETER           :: N003125  = 0.031250D0
    REAL(DOUBLE),PARAMETER           :: N0015625 = 0.015625D0
!
!   Factors
!
    PiZeta = (Pi/Zeta)**(threehalves)
!
    SELECT CASE (LQ)
    CASE (0)  
       RC(0) = PiZeta*Coefs(1)
       RS(0) = Zero
    CASE (1)  
       RC(0) = PiZeta*Coefs(1)
       RS(0) = Zero
       RC(1) = PiZeta*Coefs(4)
       RS(1) = Zero
       RC(2) = -N05*PiZeta*Coefs(3)
       RS(2) = -N05*PiZeta*Coefs(2)
    CASE (2)  
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
    CASE (3) 
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
    CASE (4)  
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
    CASE (5)
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
       RC(15) = PiZeta*(N0375*Coefs(42) + N0125*Coefs(44) + N0375*Coefs(46) - &
                        N05*Coefs(51) - N05*Coefs(53) + Coefs(56))
       RS(15) = Zero
       RC(16) = PiZeta*(-N00625*Coefs(37) - N00625*Coefs(39) - N03125*Coefs(41) + &
                         N0125*Coefs(48) + N0375*Coefs(50) - N05*Coefs(55))
       RS(16) = PiZeta*(-N03125*Coefs(36) - N00625*Coefs(38) - N00625*Coefs(40) + &
                         N0375*Coefs(47) + N0125*Coefs(49) - N05*Coefs(54))
       RC(17) = PiZeta*(N025*Coefs(42) - N025*Coefs(46) - N025*Coefs(51) + N025*Coefs(53))
       RS(17) = PiZeta*(-(N0125*Coefs(43)) - N0125*Coefs(45) + N025*Coefs(52))
       RC(18) = PiZeta*(-N009375*Coefs(37) - N003125*Coefs(39) + N015625*Coefs(41) + &
                         N0125*Coefs(48) - N0125*Coefs(50))
       RS(18) = PiZeta*(-N015625*Coefs(36) + N003125*Coefs(38) + N009375*Coefs(40) + &
                          N0125*Coefs(47) - N0125*Coefs(49))
       RC(19) = PiZeta*(N00625*Coefs(42) - N00625*Coefs(44) + N00625*Coefs(46))
       RS(19) = PiZeta*(-(N00625*Coefs(43)) + N00625*Coefs(45))
       RC(20) = PiZeta*(-(N003125*Coefs(37)) + N003125*Coefs(39) - N003125*Coefs(41))
       RS(20) = PiZeta*(-(N003125*Coefs(36)) + N003125*Coefs(38) - N003125*Coefs(40))
    CASE (6)
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
       RC(15) = PiZeta*(N0375*Coefs(42) + N0125*Coefs(44) + N0375*Coefs(46) - &
                        N05*Coefs(51) - N05*Coefs(53) + Coefs(56))
       RS(15) = Zero
       RC(16) = PiZeta*(-N00625*Coefs(37) - N00625*Coefs(39) - N03125*Coefs(41) + &
                         N0125*Coefs(48) + N0375*Coefs(50) - N05*Coefs(55))
       RS(16) = PiZeta*(-N03125*Coefs(36) - N00625*Coefs(38) - N00625*Coefs(40) + &
                         N0375*Coefs(47) + N0125*Coefs(49) - N05*Coefs(54))
       RC(17) = PiZeta*(N025*Coefs(42) - N025*Coefs(46) - N025*Coefs(51) + N025*Coefs(53))
       RS(17) = PiZeta*(-(N0125*Coefs(43)) - N0125*Coefs(45) + N025*Coefs(52))
       RC(18) = PiZeta*(-N009375*Coefs(37) - N003125*Coefs(39) + N015625*Coefs(41) + &
                         N0125*Coefs(48) - N0125*Coefs(50))
       RS(18) = PiZeta*(-N015625*Coefs(36) + N003125*Coefs(38) + N009375*Coefs(40) + &
                          N0125*Coefs(47) - N0125*Coefs(49))
       RC(19) = PiZeta*(N00625*Coefs(42) - N00625*Coefs(44) + N00625*Coefs(46))
       RS(19) = PiZeta*(-(N00625*Coefs(43)) + N00625*Coefs(45))
       RC(20) = PiZeta*(-(N003125*Coefs(37)) + N003125*Coefs(39) - N003125*Coefs(41))
       RS(20) = PiZeta*(-(N003125*Coefs(36)) + N003125*Coefs(38) - N003125*Coefs(40))
       RC(21) = PiZeta*(-N03125*Coefs(57) - N00625*Coefs(59) - N00625*Coefs(61) - &
                         N03125*Coefs(63) + N0375*Coefs(70) + N0125*Coefs(72) + &
                         N0375*Coefs(74) - N05*Coefs(79) - N05*Coefs(81) + Coefs(84))
       RS(21) = Zero
       RC(22) = PiZeta*(-N00625*Coefs(65) - N00625*Coefs(67) - N03125*Coefs(69) + & 
                         N0125*Coefs(76) + N0375*Coefs(78) - N05*Coefs(83))
       RS(22) = PiZeta*(-N03125*Coefs(64) - N00625*Coefs(66) - N00625*Coefs(68) + &
                         N0375*Coefs(75) + N0125*Coefs(77) - N05*Coefs(82))
       RC(23) = PiZeta*(-N0234375*Coefs(57) - N0015625*Coefs(59) + N0015625*Coefs(61) + &
                         N0234375*Coefs(63) + N025*Coefs(70) - N025*Coefs(74) - &
                         N025*Coefs(79) + N025*Coefs(81))
       RS(23) = PiZeta*(N0078125*Coefs(58) + N0046875*Coefs(60) + N0078125*Coefs(62) - &
                        N0125*Coefs(71) - N0125*Coefs(73) + N025*Coefs(80))
       RC(24) = PiZeta*(-N009375*Coefs(65) - N003125*Coefs(67) + N015625*Coefs(69) + & 
                         N0125*Coefs(76) - N0125*Coefs(78))
       RS(24) = PiZeta*(-N015625*Coefs(64) + N003125*Coefs(66) + N009375*Coefs(68) + &
                         N0125*Coefs(75) - N0125*Coefs(77))
       RC(25) = PiZeta*(-N009375*Coefs(57) + N003125*Coefs(59) + N003125*Coefs(61) - &
                         N009375*Coefs(63) + N00625*Coefs(70) - N00625*Coefs(72) + &
                         N00625*Coefs(74))
       RS(25) = PiZeta*(N00625*Coefs(58) - N00625*Coefs(62) - N00625*Coefs(71) + &
                        N00625*Coefs(73))
       RC(26) = PiZeta*(-N003125*Coefs(65) + N003125*Coefs(67) - N003125*Coefs(69))
       RS(26) = PiZeta*(-N003125*Coefs(64) + N003125*Coefs(66) - N003125*Coefs(68))
       RC(27) = PiZeta*(-N0015625*Coefs(57) + N0015625*Coefs(59) - N0015625*Coefs(61) + &
                         N0015625*Coefs(63))
       RS(27) = PiZeta*(N0015625*Coefs(58) - N0015625*Coefs(60) + N0015625*Coefs(62))
    CASE DEFAULT
       WRITE(*,*) 'LQ = ',LQ,' Bad indexing in HGToSP' 
       STOP  
    END SELECT 
