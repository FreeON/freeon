   SUBROUTINE BraHRR1010(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(11,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(11)+  & 
        3.D0*HRR(21))+  & 
        3.D0*HRR(36))+  & 
        HRR(57)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(12,11|
      INTGRL(OffSet)=ABx*(ABx*(2.23606797749979D0*ABx*HRR(12)+  & 
        6.70820393249937D0*HRR(22))+  & 
        6.70820393249937D0*HRR(37))+  & 
        2.23606797749979D0*HRR(58)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(13,11|
      INTGRL(OffSet)=ABx*(ABx*(2.23606797749979D0*ABx*HRR(13)+  & 
        6.70820393249937D0*HRR(23))+  & 
        6.70820393249937D0*HRR(38))+  & 
        2.23606797749979D0*HRR(59)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=(14,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(14)+  & 
        3.D0*HRR(24))+  & 
        3.D0*HRR(39))+  & 
        HRR(60)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=(15,11|
      INTGRL(OffSet)=ABx*(ABx*(2.23606797749979D0*ABx*HRR(15)+  & 
        6.70820393249937D0*HRR(26))+  & 
        6.70820393249937D0*HRR(42))+  & 
        2.23606797749979D0*HRR(64)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=(16,11|
      INTGRL(OffSet)=ABx*(ABx*(3.87298334620742D0*ABx*HRR(16)+  & 
        1.16189500386223D1*HRR(27))+  & 
        1.16189500386223D1*HRR(43))+  & 
        3.87298334620742D0*HRR(65)
      OffSet=(OA+6)*LDA+(OB+0)*LDB+CDOffSet !=(17,11|
      INTGRL(OffSet)=ABx*(ABx*(2.23606797749979D0*ABx*HRR(17)+  & 
        6.70820393249937D0*HRR(28))+  & 
        6.70820393249937D0*HRR(44))+  & 
        2.23606797749979D0*HRR(66)
      OffSet=(OA+7)*LDA+(OB+0)*LDB+CDOffSet !=(18,11|
      INTGRL(OffSet)=ABx*(ABx*(2.23606797749979D0*ABx*HRR(18)+  & 
        6.70820393249937D0*HRR(30))+  & 
        6.70820393249937D0*HRR(47))+  & 
        2.23606797749979D0*HRR(70)
      OffSet=(OA+8)*LDA+(OB+0)*LDB+CDOffSet !=(19,11|
      INTGRL(OffSet)=ABx*(ABx*(2.23606797749979D0*ABx*HRR(19)+  & 
        6.70820393249937D0*HRR(31))+  & 
        6.70820393249937D0*HRR(48))+  & 
        2.23606797749979D0*HRR(71)
      OffSet=(OA+9)*LDA+(OB+0)*LDB+CDOffSet !=(20,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(20)+  & 
        3.D0*HRR(33))+  & 
        3.D0*HRR(51))+  & 
        HRR(75)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(11,12|
      INTGRL(OffSet)=2.23606797749979D0*ABy*HRR(36)+  & 
        ABx*(4.47213595499958D0*ABy*HRR(21)+  & 
        ABx*(2.23606797749979D0*ABy*HRR(11)+  & 
        2.23606797749979D0*HRR(22))+  & 
        4.47213595499958D0*HRR(37))+  & 
        2.23606797749979D0*HRR(58)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=(12,12|
      INTGRL(OffSet)=5.D0*ABy*HRR(37)+  & 
        ABx*(1.D1*ABy*HRR(22)+  & 
        ABx*(5.D0*ABy*HRR(12)+  & 
        5.D0*HRR(23))+  & 
        1.D1*HRR(38))+  & 
        5.D0*HRR(59)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=(13,12|
      INTGRL(OffSet)=5.D0*ABy*HRR(38)+  & 
        ABx*(1.D1*ABy*HRR(23)+  & 
        ABx*(5.D0*ABy*HRR(13)+  & 
        5.D0*HRR(24))+  & 
        1.D1*HRR(39))+  & 
        5.D0*HRR(60)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=(14,12|
      INTGRL(OffSet)=2.23606797749979D0*ABy*HRR(39)+  & 
        ABx*(4.47213595499958D0*ABy*HRR(24)+  & 
        ABx*(2.23606797749979D0*ABy*HRR(14)+  & 
        2.23606797749979D0*HRR(25))+  & 
        4.47213595499958D0*HRR(40))+  & 
        2.23606797749979D0*HRR(61)
      OffSet=(OA+4)*LDA+(OB+1)*LDB+CDOffSet !=(15,12|
      INTGRL(OffSet)=5.D0*ABy*HRR(42)+  & 
        ABx*(1.D1*ABy*HRR(26)+  & 
        ABx*(5.D0*ABy*HRR(15)+  & 
        5.D0*HRR(27))+  & 
        1.D1*HRR(43))+  & 
        5.D0*HRR(65)
      OffSet=(OA+5)*LDA+(OB+1)*LDB+CDOffSet !=(16,12|
      INTGRL(OffSet)=8.66025403784439D0*ABy*HRR(43)+  & 
        ABx*(1.73205080756888D1*ABy*HRR(27)+  & 
        ABx*(8.66025403784439D0*ABy*HRR(16)+  & 
        8.66025403784439D0*HRR(28))+  & 
        1.73205080756888D1*HRR(44))+  & 
        8.66025403784439D0*HRR(66)
      OffSet=(OA+6)*LDA+(OB+1)*LDB+CDOffSet !=(17,12|
      INTGRL(OffSet)=5.D0*ABy*HRR(44)+  & 
        ABx*(1.D1*ABy*HRR(28)+  & 
        ABx*(5.D0*ABy*HRR(17)+  & 
        5.D0*HRR(29))+  & 
        1.D1*HRR(45))+  & 
        5.D0*HRR(67)
      OffSet=(OA+7)*LDA+(OB+1)*LDB+CDOffSet !=(18,12|
      INTGRL(OffSet)=5.D0*ABy*HRR(47)+  & 
        ABx*(1.D1*ABy*HRR(30)+  & 
        ABx*(5.D0*ABy*HRR(18)+  & 
        5.D0*HRR(31))+  & 
        1.D1*HRR(48))+  & 
        5.D0*HRR(71)
      OffSet=(OA+8)*LDA+(OB+1)*LDB+CDOffSet !=(19,12|
      INTGRL(OffSet)=5.D0*ABy*HRR(48)+  & 
        ABx*(1.D1*ABy*HRR(31)+  & 
        ABx*(5.D0*ABy*HRR(19)+  & 
        5.D0*HRR(32))+  & 
        1.D1*HRR(49))+  & 
        5.D0*HRR(72)
      OffSet=(OA+9)*LDA+(OB+1)*LDB+CDOffSet !=(20,12|
      INTGRL(OffSet)=2.23606797749979D0*ABy*HRR(51)+  & 
        ABx*(4.47213595499958D0*ABy*HRR(33)+  & 
        ABx*(2.23606797749979D0*ABy*HRR(20)+  & 
        2.23606797749979D0*HRR(34))+  & 
        4.47213595499958D0*HRR(52))+  & 
        2.23606797749979D0*HRR(76)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(11,13|
      INTGRL(OffSet)=ABy*(2.23606797749979D0*ABy*HRR(21)+  & 
        4.47213595499958D0*HRR(37))+  & 
        ABx*(ABy*(2.23606797749979D0*ABy*HRR(11)+  & 
        4.47213595499958D0*HRR(22))+  & 
        2.23606797749979D0*HRR(38))+  & 
        2.23606797749979D0*HRR(59)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=(12,13|
      INTGRL(OffSet)=ABy*(5.D0*ABy*HRR(22)+  & 
        1.D1*HRR(38))+  & 
        ABx*(ABy*(5.D0*ABy*HRR(12)+  & 
        1.D1*HRR(23))+  & 
        5.D0*HRR(39))+  & 
        5.D0*HRR(60)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=(13,13|
      INTGRL(OffSet)=ABy*(5.D0*ABy*HRR(23)+  & 
        1.D1*HRR(39))+  & 
        ABx*(ABy*(5.D0*ABy*HRR(13)+  & 
        1.D1*HRR(24))+  & 
        5.D0*HRR(40))+  & 
        5.D0*HRR(61)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=(14,13|
      INTGRL(OffSet)=ABy*(2.23606797749979D0*ABy*HRR(24)+  & 
        4.47213595499958D0*HRR(40))+  & 
        ABx*(ABy*(2.23606797749979D0*ABy*HRR(14)+  & 
        4.47213595499958D0*HRR(25))+  & 
        2.23606797749979D0*HRR(41))+  & 
        2.23606797749979D0*HRR(62)
      OffSet=(OA+4)*LDA+(OB+2)*LDB+CDOffSet !=(15,13|
      INTGRL(OffSet)=ABy*(5.D0*ABy*HRR(26)+  & 
        1.D1*HRR(43))+  & 
        ABx*(ABy*(5.D0*ABy*HRR(15)+  & 
        1.D1*HRR(27))+  & 
        5.D0*HRR(44))+  & 
        5.D0*HRR(66)
      OffSet=(OA+5)*LDA+(OB+2)*LDB+CDOffSet !=(16,13|
      INTGRL(OffSet)=ABy*(8.66025403784439D0*ABy*HRR(27)+  & 
        1.73205080756888D1*HRR(44))+  & 
        ABx*(ABy*(8.66025403784439D0*ABy*HRR(16)+  & 
        1.73205080756888D1*HRR(28))+  & 
        8.66025403784439D0*HRR(45))+  & 
        8.66025403784439D0*HRR(67)
      OffSet=(OA+6)*LDA+(OB+2)*LDB+CDOffSet !=(17,13|
      INTGRL(OffSet)=ABy*(5.D0*ABy*HRR(28)+  & 
        1.D1*HRR(45))+  & 
        ABx*(ABy*(5.D0*ABy*HRR(17)+  & 
        1.D1*HRR(29))+  & 
        5.D0*HRR(46))+  & 
        5.D0*HRR(68)
      OffSet=(OA+7)*LDA+(OB+2)*LDB+CDOffSet !=(18,13|
      INTGRL(OffSet)=ABy*(5.D0*ABy*HRR(30)+  & 
        1.D1*HRR(48))+  & 
        ABx*(ABy*(5.D0*ABy*HRR(18)+  & 
        1.D1*HRR(31))+  & 
        5.D0*HRR(49))+  & 
        5.D0*HRR(72)
      OffSet=(OA+8)*LDA+(OB+2)*LDB+CDOffSet !=(19,13|
      INTGRL(OffSet)=ABy*(5.D0*ABy*HRR(31)+  & 
        1.D1*HRR(49))+  & 
        ABx*(ABy*(5.D0*ABy*HRR(19)+  & 
        1.D1*HRR(32))+  & 
        5.D0*HRR(50))+  & 
        5.D0*HRR(73)
      OffSet=(OA+9)*LDA+(OB+2)*LDB+CDOffSet !=(20,13|
      INTGRL(OffSet)=ABy*(2.23606797749979D0*ABy*HRR(33)+  & 
        4.47213595499958D0*HRR(52))+  & 
        ABx*(ABy*(2.23606797749979D0*ABy*HRR(20)+  & 
        4.47213595499958D0*HRR(34))+  & 
        2.23606797749979D0*HRR(53))+  & 
        2.23606797749979D0*HRR(77)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=(11,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(11)+  & 
        3.D0*HRR(22))+  & 
        3.D0*HRR(38))+  & 
        HRR(60)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=(12,14|
      INTGRL(OffSet)=ABy*(ABy*(2.23606797749979D0*ABy*HRR(12)+  & 
        6.70820393249937D0*HRR(23))+  & 
        6.70820393249937D0*HRR(39))+  & 
        2.23606797749979D0*HRR(61)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=(13,14|
      INTGRL(OffSet)=ABy*(ABy*(2.23606797749979D0*ABy*HRR(13)+  & 
        6.70820393249937D0*HRR(24))+  & 
        6.70820393249937D0*HRR(40))+  & 
        2.23606797749979D0*HRR(62)
      OffSet=(OA+3)*LDA+(OB+3)*LDB+CDOffSet !=(14,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(14)+  & 
        3.D0*HRR(25))+  & 
        3.D0*HRR(41))+  & 
        HRR(63)
      OffSet=(OA+4)*LDA+(OB+3)*LDB+CDOffSet !=(15,14|
      INTGRL(OffSet)=ABy*(ABy*(2.23606797749979D0*ABy*HRR(15)+  & 
        6.70820393249937D0*HRR(27))+  & 
        6.70820393249937D0*HRR(44))+  & 
        2.23606797749979D0*HRR(67)
      OffSet=(OA+5)*LDA+(OB+3)*LDB+CDOffSet !=(16,14|
      INTGRL(OffSet)=ABy*(ABy*(3.87298334620742D0*ABy*HRR(16)+  & 
        1.16189500386223D1*HRR(28))+  & 
        1.16189500386223D1*HRR(45))+  & 
        3.87298334620742D0*HRR(68)
      OffSet=(OA+6)*LDA+(OB+3)*LDB+CDOffSet !=(17,14|
      INTGRL(OffSet)=ABy*(ABy*(2.23606797749979D0*ABy*HRR(17)+  & 
        6.70820393249937D0*HRR(29))+  & 
        6.70820393249937D0*HRR(46))+  & 
        2.23606797749979D0*HRR(69)
      OffSet=(OA+7)*LDA+(OB+3)*LDB+CDOffSet !=(18,14|
      INTGRL(OffSet)=ABy*(ABy*(2.23606797749979D0*ABy*HRR(18)+  & 
        6.70820393249937D0*HRR(31))+  & 
        6.70820393249937D0*HRR(49))+  & 
        2.23606797749979D0*HRR(73)
      OffSet=(OA+8)*LDA+(OB+3)*LDB+CDOffSet !=(19,14|
      INTGRL(OffSet)=ABy*(ABy*(2.23606797749979D0*ABy*HRR(19)+  & 
        6.70820393249937D0*HRR(32))+  & 
        6.70820393249937D0*HRR(50))+  & 
        2.23606797749979D0*HRR(74)
      OffSet=(OA+9)*LDA+(OB+3)*LDB+CDOffSet !=(20,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(20)+  & 
        3.D0*HRR(34))+  & 
        3.D0*HRR(53))+  & 
        HRR(78)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=(11,15|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(36)+  & 
        ABx*(4.47213595499958D0*ABz*HRR(21)+  & 
        ABx*(2.23606797749979D0*ABz*HRR(11)+  & 
        2.23606797749979D0*HRR(26))+  & 
        4.47213595499958D0*HRR(42))+  & 
        2.23606797749979D0*HRR(64)
      OffSet=(OA+1)*LDA+(OB+4)*LDB+CDOffSet !=(12,15|
      INTGRL(OffSet)=5.D0*ABz*HRR(37)+  & 
        ABx*(1.D1*ABz*HRR(22)+  & 
        ABx*(5.D0*ABz*HRR(12)+  & 
        5.D0*HRR(27))+  & 
        1.D1*HRR(43))+  & 
        5.D0*HRR(65)
      OffSet=(OA+2)*LDA+(OB+4)*LDB+CDOffSet !=(13,15|
      INTGRL(OffSet)=5.D0*ABz*HRR(38)+  & 
        ABx*(1.D1*ABz*HRR(23)+  & 
        ABx*(5.D0*ABz*HRR(13)+  & 
        5.D0*HRR(28))+  & 
        1.D1*HRR(44))+  & 
        5.D0*HRR(66)
      OffSet=(OA+3)*LDA+(OB+4)*LDB+CDOffSet !=(14,15|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(39)+  & 
        ABx*(4.47213595499958D0*ABz*HRR(24)+  & 
        ABx*(2.23606797749979D0*ABz*HRR(14)+  & 
        2.23606797749979D0*HRR(29))+  & 
        4.47213595499958D0*HRR(45))+  & 
        2.23606797749979D0*HRR(67)
      OffSet=(OA+4)*LDA+(OB+4)*LDB+CDOffSet !=(15,15|
      INTGRL(OffSet)=5.D0*ABz*HRR(42)+  & 
        ABx*(1.D1*ABz*HRR(26)+  & 
        ABx*(5.D0*ABz*HRR(15)+  & 
        5.D0*HRR(30))+  & 
        1.D1*HRR(47))+  & 
        5.D0*HRR(70)
      OffSet=(OA+5)*LDA+(OB+4)*LDB+CDOffSet !=(16,15|
      INTGRL(OffSet)=8.66025403784439D0*ABz*HRR(43)+  & 
        ABx*(1.73205080756888D1*ABz*HRR(27)+  & 
        ABx*(8.66025403784439D0*ABz*HRR(16)+  & 
        8.66025403784439D0*HRR(31))+  & 
        1.73205080756888D1*HRR(48))+  & 
        8.66025403784439D0*HRR(71)
      OffSet=(OA+6)*LDA+(OB+4)*LDB+CDOffSet !=(17,15|
      INTGRL(OffSet)=5.D0*ABz*HRR(44)+  & 
        ABx*(1.D1*ABz*HRR(28)+  & 
        ABx*(5.D0*ABz*HRR(17)+  & 
        5.D0*HRR(32))+  & 
        1.D1*HRR(49))+  & 
        5.D0*HRR(72)
      OffSet=(OA+7)*LDA+(OB+4)*LDB+CDOffSet !=(18,15|
      INTGRL(OffSet)=5.D0*ABz*HRR(47)+  & 
        ABx*(1.D1*ABz*HRR(30)+  & 
        ABx*(5.D0*ABz*HRR(18)+  & 
        5.D0*HRR(33))+  & 
        1.D1*HRR(51))+  & 
        5.D0*HRR(75)
      OffSet=(OA+8)*LDA+(OB+4)*LDB+CDOffSet !=(19,15|
      INTGRL(OffSet)=5.D0*ABz*HRR(48)+  & 
        ABx*(1.D1*ABz*HRR(31)+  & 
        ABx*(5.D0*ABz*HRR(19)+  & 
        5.D0*HRR(34))+  & 
        1.D1*HRR(52))+  & 
        5.D0*HRR(76)
      OffSet=(OA+9)*LDA+(OB+4)*LDB+CDOffSet !=(20,15|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(51)+  & 
        ABx*(4.47213595499958D0*ABz*HRR(33)+  & 
        ABx*(2.23606797749979D0*ABz*HRR(20)+  & 
        2.23606797749979D0*HRR(35))+  & 
        4.47213595499958D0*HRR(54))+  & 
        2.23606797749979D0*HRR(79)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=(11,16|
      INTGRL(OffSet)=3.87298334620742D0*ABz*HRR(37)+  & 
        ABy*(3.87298334620742D0*ABz*HRR(21)+  & 
        3.87298334620742D0*HRR(42))+  & 
        ABx*(3.87298334620742D0*ABz*HRR(22)+  & 
        ABy*(3.87298334620742D0*ABz*HRR(11)+  & 
        3.87298334620742D0*HRR(26))+  & 
        3.87298334620742D0*HRR(43))+  & 
        3.87298334620742D0*HRR(65)
      OffSet=(OA+1)*LDA+(OB+5)*LDB+CDOffSet !=(12,16|
      INTGRL(OffSet)=8.66025403784439D0*ABz*HRR(38)+  & 
        ABy*(8.66025403784439D0*ABz*HRR(22)+  & 
        8.66025403784439D0*HRR(43))+  & 
        ABx*(8.66025403784439D0*ABz*HRR(23)+  & 
        ABy*(8.66025403784439D0*ABz*HRR(12)+  & 
        8.66025403784439D0*HRR(27))+  & 
        8.66025403784439D0*HRR(44))+  & 
        8.66025403784439D0*HRR(66)
      OffSet=(OA+2)*LDA+(OB+5)*LDB+CDOffSet !=(13,16|
      INTGRL(OffSet)=8.66025403784439D0*ABz*HRR(39)+  & 
        ABy*(8.66025403784439D0*ABz*HRR(23)+  & 
        8.66025403784439D0*HRR(44))+  & 
        ABx*(8.66025403784439D0*ABz*HRR(24)+  & 
        ABy*(8.66025403784439D0*ABz*HRR(13)+  & 
        8.66025403784439D0*HRR(28))+  & 
        8.66025403784439D0*HRR(45))+  & 
        8.66025403784439D0*HRR(67)
      OffSet=(OA+3)*LDA+(OB+5)*LDB+CDOffSet !=(14,16|
      INTGRL(OffSet)=3.87298334620742D0*ABz*HRR(40)+  & 
        ABy*(3.87298334620742D0*ABz*HRR(24)+  & 
        3.87298334620742D0*HRR(45))+  & 
        ABx*(3.87298334620742D0*ABz*HRR(25)+  & 
        ABy*(3.87298334620742D0*ABz*HRR(14)+  & 
        3.87298334620742D0*HRR(29))+  & 
        3.87298334620742D0*HRR(46))+  & 
        3.87298334620742D0*HRR(68)
      OffSet=(OA+4)*LDA+(OB+5)*LDB+CDOffSet !=(15,16|
      INTGRL(OffSet)=8.66025403784439D0*ABz*HRR(43)+  & 
        ABy*(8.66025403784439D0*ABz*HRR(26)+  & 
        8.66025403784439D0*HRR(47))+  & 
        ABx*(8.66025403784439D0*ABz*HRR(27)+  & 
        ABy*(8.66025403784439D0*ABz*HRR(15)+  & 
        8.66025403784439D0*HRR(30))+  & 
        8.66025403784439D0*HRR(48))+  & 
        8.66025403784439D0*HRR(71)
      OffSet=(OA+5)*LDA+(OB+5)*LDB+CDOffSet !=(16,16|
      INTGRL(OffSet)=1.5D1*ABz*HRR(44)+  & 
        ABy*(1.5D1*ABz*HRR(27)+  & 
        1.5D1*HRR(48))+  & 
        ABx*(1.5D1*ABz*HRR(28)+  & 
        ABy*(1.5D1*ABz*HRR(16)+  & 
        1.5D1*HRR(31))+  & 
        1.5D1*HRR(49))+  & 
        1.5D1*HRR(72)
      OffSet=(OA+6)*LDA+(OB+5)*LDB+CDOffSet !=(17,16|
      INTGRL(OffSet)=8.66025403784439D0*ABz*HRR(45)+  & 
        ABy*(8.66025403784439D0*ABz*HRR(28)+  & 
        8.66025403784439D0*HRR(49))+  & 
        ABx*(8.66025403784439D0*ABz*HRR(29)+  & 
        ABy*(8.66025403784439D0*ABz*HRR(17)+  & 
        8.66025403784439D0*HRR(32))+  & 
        8.66025403784439D0*HRR(50))+  & 
        8.66025403784439D0*HRR(73)
      OffSet=(OA+7)*LDA+(OB+5)*LDB+CDOffSet !=(18,16|
      INTGRL(OffSet)=8.66025403784439D0*ABz*HRR(48)+  & 
        ABy*(8.66025403784439D0*ABz*HRR(30)+  & 
        8.66025403784439D0*HRR(51))+  & 
        ABx*(8.66025403784439D0*ABz*HRR(31)+  & 
        ABy*(8.66025403784439D0*ABz*HRR(18)+  & 
        8.66025403784439D0*HRR(33))+  & 
        8.66025403784439D0*HRR(52))+  & 
        8.66025403784439D0*HRR(76)
      OffSet=(OA+8)*LDA+(OB+5)*LDB+CDOffSet !=(19,16|
      INTGRL(OffSet)=8.66025403784439D0*ABz*HRR(49)+  & 
        ABy*(8.66025403784439D0*ABz*HRR(31)+  & 
        8.66025403784439D0*HRR(52))+  & 
        ABx*(8.66025403784439D0*ABz*HRR(32)+  & 
        ABy*(8.66025403784439D0*ABz*HRR(19)+  & 
        8.66025403784439D0*HRR(34))+  & 
        8.66025403784439D0*HRR(53))+  & 
        8.66025403784439D0*HRR(77)
      OffSet=(OA+9)*LDA+(OB+5)*LDB+CDOffSet !=(20,16|
      INTGRL(OffSet)=3.87298334620742D0*ABz*HRR(52)+  & 
        ABy*(3.87298334620742D0*ABz*HRR(33)+  & 
        3.87298334620742D0*HRR(54))+  & 
        ABx*(3.87298334620742D0*ABz*HRR(34)+  & 
        ABy*(3.87298334620742D0*ABz*HRR(20)+  & 
        3.87298334620742D0*HRR(35))+  & 
        3.87298334620742D0*HRR(55))+  & 
        3.87298334620742D0*HRR(80)
      OffSet=(OA+0)*LDA+(OB+6)*LDB+CDOffSet !=(11,17|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(38)+  & 
        ABy*(4.47213595499958D0*ABz*HRR(22)+  & 
        ABy*(2.23606797749979D0*ABz*HRR(11)+  & 
        2.23606797749979D0*HRR(26))+  & 
        4.47213595499958D0*HRR(43))+  & 
        2.23606797749979D0*HRR(66)
      OffSet=(OA+1)*LDA+(OB+6)*LDB+CDOffSet !=(12,17|
      INTGRL(OffSet)=5.D0*ABz*HRR(39)+  & 
        ABy*(1.D1*ABz*HRR(23)+  & 
        ABy*(5.D0*ABz*HRR(12)+  & 
        5.D0*HRR(27))+  & 
        1.D1*HRR(44))+  & 
        5.D0*HRR(67)
      OffSet=(OA+2)*LDA+(OB+6)*LDB+CDOffSet !=(13,17|
      INTGRL(OffSet)=5.D0*ABz*HRR(40)+  & 
        ABy*(1.D1*ABz*HRR(24)+  & 
        ABy*(5.D0*ABz*HRR(13)+  & 
        5.D0*HRR(28))+  & 
        1.D1*HRR(45))+  & 
        5.D0*HRR(68)
      OffSet=(OA+3)*LDA+(OB+6)*LDB+CDOffSet !=(14,17|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(41)+  & 
        ABy*(4.47213595499958D0*ABz*HRR(25)+  & 
        ABy*(2.23606797749979D0*ABz*HRR(14)+  & 
        2.23606797749979D0*HRR(29))+  & 
        4.47213595499958D0*HRR(46))+  & 
        2.23606797749979D0*HRR(69)
      OffSet=(OA+4)*LDA+(OB+6)*LDB+CDOffSet !=(15,17|
      INTGRL(OffSet)=5.D0*ABz*HRR(44)+  & 
        ABy*(1.D1*ABz*HRR(27)+  & 
        ABy*(5.D0*ABz*HRR(15)+  & 
        5.D0*HRR(30))+  & 
        1.D1*HRR(48))+  & 
        5.D0*HRR(72)
      OffSet=(OA+5)*LDA+(OB+6)*LDB+CDOffSet !=(16,17|
      INTGRL(OffSet)=8.66025403784439D0*ABz*HRR(45)+  & 
        ABy*(1.73205080756888D1*ABz*HRR(28)+  & 
        ABy*(8.66025403784439D0*ABz*HRR(16)+  & 
        8.66025403784439D0*HRR(31))+  & 
        1.73205080756888D1*HRR(49))+  & 
        8.66025403784439D0*HRR(73)
      OffSet=(OA+6)*LDA+(OB+6)*LDB+CDOffSet !=(17,17|
      INTGRL(OffSet)=5.D0*ABz*HRR(46)+  & 
        ABy*(1.D1*ABz*HRR(29)+  & 
        ABy*(5.D0*ABz*HRR(17)+  & 
        5.D0*HRR(32))+  & 
        1.D1*HRR(50))+  & 
        5.D0*HRR(74)
      OffSet=(OA+7)*LDA+(OB+6)*LDB+CDOffSet !=(18,17|
      INTGRL(OffSet)=5.D0*ABz*HRR(49)+  & 
        ABy*(1.D1*ABz*HRR(31)+  & 
        ABy*(5.D0*ABz*HRR(18)+  & 
        5.D0*HRR(33))+  & 
        1.D1*HRR(52))+  & 
        5.D0*HRR(77)
      OffSet=(OA+8)*LDA+(OB+6)*LDB+CDOffSet !=(19,17|
      INTGRL(OffSet)=5.D0*ABz*HRR(50)+  & 
        ABy*(1.D1*ABz*HRR(32)+  & 
        ABy*(5.D0*ABz*HRR(19)+  & 
        5.D0*HRR(34))+  & 
        1.D1*HRR(53))+  & 
        5.D0*HRR(78)
      OffSet=(OA+9)*LDA+(OB+6)*LDB+CDOffSet !=(20,17|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(53)+  & 
        ABy*(4.47213595499958D0*ABz*HRR(34)+  & 
        ABy*(2.23606797749979D0*ABz*HRR(20)+  & 
        2.23606797749979D0*HRR(35))+  & 
        4.47213595499958D0*HRR(55))+  & 
        2.23606797749979D0*HRR(81)
      OffSet=(OA+0)*LDA+(OB+7)*LDB+CDOffSet !=(11,18|
      INTGRL(OffSet)=ABz*(2.23606797749979D0*ABz*HRR(21)+  & 
        4.47213595499958D0*HRR(42))+  & 
        ABx*(ABz*(2.23606797749979D0*ABz*HRR(11)+  & 
        4.47213595499958D0*HRR(26))+  & 
        2.23606797749979D0*HRR(47))+  & 
        2.23606797749979D0*HRR(70)
      OffSet=(OA+1)*LDA+(OB+7)*LDB+CDOffSet !=(12,18|
      INTGRL(OffSet)=ABz*(5.D0*ABz*HRR(22)+  & 
        1.D1*HRR(43))+  & 
        ABx*(ABz*(5.D0*ABz*HRR(12)+  & 
        1.D1*HRR(27))+  & 
        5.D0*HRR(48))+  & 
        5.D0*HRR(71)
      OffSet=(OA+2)*LDA+(OB+7)*LDB+CDOffSet !=(13,18|
      INTGRL(OffSet)=ABz*(5.D0*ABz*HRR(23)+  & 
        1.D1*HRR(44))+  & 
        ABx*(ABz*(5.D0*ABz*HRR(13)+  & 
        1.D1*HRR(28))+  & 
        5.D0*HRR(49))+  & 
        5.D0*HRR(72)
      OffSet=(OA+3)*LDA+(OB+7)*LDB+CDOffSet !=(14,18|
      INTGRL(OffSet)=ABz*(2.23606797749979D0*ABz*HRR(24)+  & 
        4.47213595499958D0*HRR(45))+  & 
        ABx*(ABz*(2.23606797749979D0*ABz*HRR(14)+  & 
        4.47213595499958D0*HRR(29))+  & 
        2.23606797749979D0*HRR(50))+  & 
        2.23606797749979D0*HRR(73)
      OffSet=(OA+4)*LDA+(OB+7)*LDB+CDOffSet !=(15,18|
      INTGRL(OffSet)=ABz*(5.D0*ABz*HRR(26)+  & 
        1.D1*HRR(47))+  & 
        ABx*(ABz*(5.D0*ABz*HRR(15)+  & 
        1.D1*HRR(30))+  & 
        5.D0*HRR(51))+  & 
        5.D0*HRR(75)
      OffSet=(OA+5)*LDA+(OB+7)*LDB+CDOffSet !=(16,18|
      INTGRL(OffSet)=ABz*(8.66025403784439D0*ABz*HRR(27)+  & 
        1.73205080756888D1*HRR(48))+  & 
        ABx*(ABz*(8.66025403784439D0*ABz*HRR(16)+  & 
        1.73205080756888D1*HRR(31))+  & 
        8.66025403784439D0*HRR(52))+  & 
        8.66025403784439D0*HRR(76)
      OffSet=(OA+6)*LDA+(OB+7)*LDB+CDOffSet !=(17,18|
      INTGRL(OffSet)=ABz*(5.D0*ABz*HRR(28)+  & 
        1.D1*HRR(49))+  & 
        ABx*(ABz*(5.D0*ABz*HRR(17)+  & 
        1.D1*HRR(32))+  & 
        5.D0*HRR(53))+  & 
        5.D0*HRR(77)
      OffSet=(OA+7)*LDA+(OB+7)*LDB+CDOffSet !=(18,18|
      INTGRL(OffSet)=ABz*(5.D0*ABz*HRR(30)+  & 
        1.D1*HRR(51))+  & 
        ABx*(ABz*(5.D0*ABz*HRR(18)+  & 
        1.D1*HRR(33))+  & 
        5.D0*HRR(54))+  & 
        5.D0*HRR(79)
      OffSet=(OA+8)*LDA+(OB+7)*LDB+CDOffSet !=(19,18|
      INTGRL(OffSet)=ABz*(5.D0*ABz*HRR(31)+  & 
        1.D1*HRR(52))+  & 
        ABx*(ABz*(5.D0*ABz*HRR(19)+  & 
        1.D1*HRR(34))+  & 
        5.D0*HRR(55))+  & 
        5.D0*HRR(80)
      OffSet=(OA+9)*LDA+(OB+7)*LDB+CDOffSet !=(20,18|
      INTGRL(OffSet)=ABz*(2.23606797749979D0*ABz*HRR(33)+  & 
        4.47213595499958D0*HRR(54))+  & 
        ABx*(ABz*(2.23606797749979D0*ABz*HRR(20)+  & 
        4.47213595499958D0*HRR(35))+  & 
        2.23606797749979D0*HRR(56))+  & 
        2.23606797749979D0*HRR(82)
      OffSet=(OA+0)*LDA+(OB+8)*LDB+CDOffSet !=(11,19|
      INTGRL(OffSet)=ABz*(2.23606797749979D0*ABz*HRR(22)+  & 
        4.47213595499958D0*HRR(43))+  & 
        ABy*(ABz*(2.23606797749979D0*ABz*HRR(11)+  & 
        4.47213595499958D0*HRR(26))+  & 
        2.23606797749979D0*HRR(47))+  & 
        2.23606797749979D0*HRR(71)
      OffSet=(OA+1)*LDA+(OB+8)*LDB+CDOffSet !=(12,19|
      INTGRL(OffSet)=ABz*(5.D0*ABz*HRR(23)+  & 
        1.D1*HRR(44))+  & 
        ABy*(ABz*(5.D0*ABz*HRR(12)+  & 
        1.D1*HRR(27))+  & 
        5.D0*HRR(48))+  & 
        5.D0*HRR(72)
      OffSet=(OA+2)*LDA+(OB+8)*LDB+CDOffSet !=(13,19|
      INTGRL(OffSet)=ABz*(5.D0*ABz*HRR(24)+  & 
        1.D1*HRR(45))+  & 
        ABy*(ABz*(5.D0*ABz*HRR(13)+  & 
        1.D1*HRR(28))+  & 
        5.D0*HRR(49))+  & 
        5.D0*HRR(73)
      OffSet=(OA+3)*LDA+(OB+8)*LDB+CDOffSet !=(14,19|
      INTGRL(OffSet)=ABz*(2.23606797749979D0*ABz*HRR(25)+  & 
        4.47213595499958D0*HRR(46))+  & 
        ABy*(ABz*(2.23606797749979D0*ABz*HRR(14)+  & 
        4.47213595499958D0*HRR(29))+  & 
        2.23606797749979D0*HRR(50))+  & 
        2.23606797749979D0*HRR(74)
      OffSet=(OA+4)*LDA+(OB+8)*LDB+CDOffSet !=(15,19|
      INTGRL(OffSet)=ABz*(5.D0*ABz*HRR(27)+  & 
        1.D1*HRR(48))+  & 
        ABy*(ABz*(5.D0*ABz*HRR(15)+  & 
        1.D1*HRR(30))+  & 
        5.D0*HRR(51))+  & 
        5.D0*HRR(76)
      OffSet=(OA+5)*LDA+(OB+8)*LDB+CDOffSet !=(16,19|
      INTGRL(OffSet)=ABz*(8.66025403784439D0*ABz*HRR(28)+  & 
        1.73205080756888D1*HRR(49))+  & 
        ABy*(ABz*(8.66025403784439D0*ABz*HRR(16)+  & 
        1.73205080756888D1*HRR(31))+  & 
        8.66025403784439D0*HRR(52))+  & 
        8.66025403784439D0*HRR(77)
      OffSet=(OA+6)*LDA+(OB+8)*LDB+CDOffSet !=(17,19|
      INTGRL(OffSet)=ABz*(5.D0*ABz*HRR(29)+  & 
        1.D1*HRR(50))+  & 
        ABy*(ABz*(5.D0*ABz*HRR(17)+  & 
        1.D1*HRR(32))+  & 
        5.D0*HRR(53))+  & 
        5.D0*HRR(78)
      OffSet=(OA+7)*LDA+(OB+8)*LDB+CDOffSet !=(18,19|
      INTGRL(OffSet)=ABz*(5.D0*ABz*HRR(31)+  & 
        1.D1*HRR(52))+  & 
        ABy*(ABz*(5.D0*ABz*HRR(18)+  & 
        1.D1*HRR(33))+  & 
        5.D0*HRR(54))+  & 
        5.D0*HRR(80)
      OffSet=(OA+8)*LDA+(OB+8)*LDB+CDOffSet !=(19,19|
      INTGRL(OffSet)=ABz*(5.D0*ABz*HRR(32)+  & 
        1.D1*HRR(53))+  & 
        ABy*(ABz*(5.D0*ABz*HRR(19)+  & 
        1.D1*HRR(34))+  & 
        5.D0*HRR(55))+  & 
        5.D0*HRR(81)
      OffSet=(OA+9)*LDA+(OB+8)*LDB+CDOffSet !=(20,19|
      INTGRL(OffSet)=ABz*(2.23606797749979D0*ABz*HRR(34)+  & 
        4.47213595499958D0*HRR(55))+  & 
        ABy*(ABz*(2.23606797749979D0*ABz*HRR(20)+  & 
        4.47213595499958D0*HRR(35))+  & 
        2.23606797749979D0*HRR(56))+  & 
        2.23606797749979D0*HRR(83)
      OffSet=(OA+0)*LDA+(OB+9)*LDB+CDOffSet !=(11,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(11)+  & 
        3.D0*HRR(26))+  & 
        3.D0*HRR(47))+  & 
        HRR(75)
      OffSet=(OA+1)*LDA+(OB+9)*LDB+CDOffSet !=(12,20|
      INTGRL(OffSet)=ABz*(ABz*(2.23606797749979D0*ABz*HRR(12)+  & 
        6.70820393249937D0*HRR(27))+  & 
        6.70820393249937D0*HRR(48))+  & 
        2.23606797749979D0*HRR(76)
      OffSet=(OA+2)*LDA+(OB+9)*LDB+CDOffSet !=(13,20|
      INTGRL(OffSet)=ABz*(ABz*(2.23606797749979D0*ABz*HRR(13)+  & 
        6.70820393249937D0*HRR(28))+  & 
        6.70820393249937D0*HRR(49))+  & 
        2.23606797749979D0*HRR(77)
      OffSet=(OA+3)*LDA+(OB+9)*LDB+CDOffSet !=(14,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(14)+  & 
        3.D0*HRR(29))+  & 
        3.D0*HRR(50))+  & 
        HRR(78)
      OffSet=(OA+4)*LDA+(OB+9)*LDB+CDOffSet !=(15,20|
      INTGRL(OffSet)=ABz*(ABz*(2.23606797749979D0*ABz*HRR(15)+  & 
        6.70820393249937D0*HRR(30))+  & 
        6.70820393249937D0*HRR(51))+  & 
        2.23606797749979D0*HRR(79)
      OffSet=(OA+5)*LDA+(OB+9)*LDB+CDOffSet !=(16,20|
      INTGRL(OffSet)=ABz*(ABz*(3.87298334620742D0*ABz*HRR(16)+  & 
        1.16189500386223D1*HRR(31))+  & 
        1.16189500386223D1*HRR(52))+  & 
        3.87298334620742D0*HRR(80)
      OffSet=(OA+6)*LDA+(OB+9)*LDB+CDOffSet !=(17,20|
      INTGRL(OffSet)=ABz*(ABz*(2.23606797749979D0*ABz*HRR(17)+  & 
        6.70820393249937D0*HRR(32))+  & 
        6.70820393249937D0*HRR(53))+  & 
        2.23606797749979D0*HRR(81)
      OffSet=(OA+7)*LDA+(OB+9)*LDB+CDOffSet !=(18,20|
      INTGRL(OffSet)=ABz*(ABz*(2.23606797749979D0*ABz*HRR(18)+  & 
        6.70820393249937D0*HRR(33))+  & 
        6.70820393249937D0*HRR(54))+  & 
        2.23606797749979D0*HRR(82)
      OffSet=(OA+8)*LDA+(OB+9)*LDB+CDOffSet !=(19,20|
      INTGRL(OffSet)=ABz*(ABz*(2.23606797749979D0*ABz*HRR(19)+  & 
        6.70820393249937D0*HRR(34))+  & 
        6.70820393249937D0*HRR(55))+  & 
        2.23606797749979D0*HRR(83)
      OffSet=(OA+9)*LDA+(OB+9)*LDB+CDOffSet !=(20,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(20)+  & 
        3.D0*HRR(35))+  & 
        3.D0*HRR(56))+  & 
        HRR(84)
END SUBROUTINE BraHRR1010