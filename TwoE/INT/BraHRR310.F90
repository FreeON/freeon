   SUBROUTINE BraHRR310(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(2,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(2)+  & 
        3.D0*HRR(5))+  & 
        3.D0*HRR(11))+  & 
        HRR(21)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(3,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(3)+  & 
        3.D0*HRR(6))+  & 
        3.D0*HRR(12))+  & 
        HRR(22)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(4,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(4)+  & 
        3.D0*HRR(8))+  & 
        3.D0*HRR(15))+  & 
        HRR(26)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(2,12|
      INTGRL(OffSet)=2.23606797749979D0*ABy*HRR(11)+  & 
        ABx*(4.47213595499958D0*ABy*HRR(5)+  & 
        ABx*(2.23606797749979D0*ABy*HRR(2)+  & 
        2.23606797749979D0*HRR(6))+  & 
        4.47213595499958D0*HRR(12))+  & 
        2.23606797749979D0*HRR(22)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=(3,12|
      INTGRL(OffSet)=2.23606797749979D0*ABy*HRR(12)+  & 
        ABx*(4.47213595499958D0*ABy*HRR(6)+  & 
        ABx*(2.23606797749979D0*ABy*HRR(3)+  & 
        2.23606797749979D0*HRR(7))+  & 
        4.47213595499958D0*HRR(13))+  & 
        2.23606797749979D0*HRR(23)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=(4,12|
      INTGRL(OffSet)=2.23606797749979D0*ABy*HRR(15)+  & 
        ABx*(4.47213595499958D0*ABy*HRR(8)+  & 
        ABx*(2.23606797749979D0*ABy*HRR(4)+  & 
        2.23606797749979D0*HRR(9))+  & 
        4.47213595499958D0*HRR(16))+  & 
        2.23606797749979D0*HRR(27)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(2,13|
      INTGRL(OffSet)=ABy*(2.23606797749979D0*ABy*HRR(5)+  & 
        4.47213595499958D0*HRR(12))+  & 
        ABx*(ABy*(2.23606797749979D0*ABy*HRR(2)+  & 
        4.47213595499958D0*HRR(6))+  & 
        2.23606797749979D0*HRR(13))+  & 
        2.23606797749979D0*HRR(23)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=(3,13|
      INTGRL(OffSet)=ABy*(2.23606797749979D0*ABy*HRR(6)+  & 
        4.47213595499958D0*HRR(13))+  & 
        ABx*(ABy*(2.23606797749979D0*ABy*HRR(3)+  & 
        4.47213595499958D0*HRR(7))+  & 
        2.23606797749979D0*HRR(14))+  & 
        2.23606797749979D0*HRR(24)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=(4,13|
      INTGRL(OffSet)=ABy*(2.23606797749979D0*ABy*HRR(8)+  & 
        4.47213595499958D0*HRR(16))+  & 
        ABx*(ABy*(2.23606797749979D0*ABy*HRR(4)+  & 
        4.47213595499958D0*HRR(9))+  & 
        2.23606797749979D0*HRR(17))+  & 
        2.23606797749979D0*HRR(28)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=(2,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(2)+  & 
        3.D0*HRR(6))+  & 
        3.D0*HRR(13))+  & 
        HRR(24)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=(3,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(3)+  & 
        3.D0*HRR(7))+  & 
        3.D0*HRR(14))+  & 
        HRR(25)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=(4,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(4)+  & 
        3.D0*HRR(9))+  & 
        3.D0*HRR(17))+  & 
        HRR(29)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=(2,15|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(11)+  & 
        ABx*(4.47213595499958D0*ABz*HRR(5)+  & 
        ABx*(2.23606797749979D0*ABz*HRR(2)+  & 
        2.23606797749979D0*HRR(8))+  & 
        4.47213595499958D0*HRR(15))+  & 
        2.23606797749979D0*HRR(26)
      OffSet=(OA+1)*LDA+(OB+4)*LDB+CDOffSet !=(3,15|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(12)+  & 
        ABx*(4.47213595499958D0*ABz*HRR(6)+  & 
        ABx*(2.23606797749979D0*ABz*HRR(3)+  & 
        2.23606797749979D0*HRR(9))+  & 
        4.47213595499958D0*HRR(16))+  & 
        2.23606797749979D0*HRR(27)
      OffSet=(OA+2)*LDA+(OB+4)*LDB+CDOffSet !=(4,15|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(15)+  & 
        ABx*(4.47213595499958D0*ABz*HRR(8)+  & 
        ABx*(2.23606797749979D0*ABz*HRR(4)+  & 
        2.23606797749979D0*HRR(10))+  & 
        4.47213595499958D0*HRR(18))+  & 
        2.23606797749979D0*HRR(30)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=(2,16|
      INTGRL(OffSet)=3.87298334620742D0*ABz*HRR(12)+  & 
        ABy*(3.87298334620742D0*ABz*HRR(5)+  & 
        3.87298334620742D0*HRR(15))+  & 
        ABx*(3.87298334620742D0*ABz*HRR(6)+  & 
        ABy*(3.87298334620742D0*ABz*HRR(2)+  & 
        3.87298334620742D0*HRR(8))+  & 
        3.87298334620742D0*HRR(16))+  & 
        3.87298334620742D0*HRR(27)
      OffSet=(OA+1)*LDA+(OB+5)*LDB+CDOffSet !=(3,16|
      INTGRL(OffSet)=3.87298334620742D0*ABz*HRR(13)+  & 
        ABy*(3.87298334620742D0*ABz*HRR(6)+  & 
        3.87298334620742D0*HRR(16))+  & 
        ABx*(3.87298334620742D0*ABz*HRR(7)+  & 
        ABy*(3.87298334620742D0*ABz*HRR(3)+  & 
        3.87298334620742D0*HRR(9))+  & 
        3.87298334620742D0*HRR(17))+  & 
        3.87298334620742D0*HRR(28)
      OffSet=(OA+2)*LDA+(OB+5)*LDB+CDOffSet !=(4,16|
      INTGRL(OffSet)=3.87298334620742D0*ABz*HRR(16)+  & 
        ABy*(3.87298334620742D0*ABz*HRR(8)+  & 
        3.87298334620742D0*HRR(18))+  & 
        ABx*(3.87298334620742D0*ABz*HRR(9)+  & 
        ABy*(3.87298334620742D0*ABz*HRR(4)+  & 
        3.87298334620742D0*HRR(10))+  & 
        3.87298334620742D0*HRR(19))+  & 
        3.87298334620742D0*HRR(31)
      OffSet=(OA+0)*LDA+(OB+6)*LDB+CDOffSet !=(2,17|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(13)+  & 
        ABy*(4.47213595499958D0*ABz*HRR(6)+  & 
        ABy*(2.23606797749979D0*ABz*HRR(2)+  & 
        2.23606797749979D0*HRR(8))+  & 
        4.47213595499958D0*HRR(16))+  & 
        2.23606797749979D0*HRR(28)
      OffSet=(OA+1)*LDA+(OB+6)*LDB+CDOffSet !=(3,17|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(14)+  & 
        ABy*(4.47213595499958D0*ABz*HRR(7)+  & 
        ABy*(2.23606797749979D0*ABz*HRR(3)+  & 
        2.23606797749979D0*HRR(9))+  & 
        4.47213595499958D0*HRR(17))+  & 
        2.23606797749979D0*HRR(29)
      OffSet=(OA+2)*LDA+(OB+6)*LDB+CDOffSet !=(4,17|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(17)+  & 
        ABy*(4.47213595499958D0*ABz*HRR(9)+  & 
        ABy*(2.23606797749979D0*ABz*HRR(4)+  & 
        2.23606797749979D0*HRR(10))+  & 
        4.47213595499958D0*HRR(19))+  & 
        2.23606797749979D0*HRR(32)
      OffSet=(OA+0)*LDA+(OB+7)*LDB+CDOffSet !=(2,18|
      INTGRL(OffSet)=ABz*(2.23606797749979D0*ABz*HRR(5)+  & 
        4.47213595499958D0*HRR(15))+  & 
        ABx*(ABz*(2.23606797749979D0*ABz*HRR(2)+  & 
        4.47213595499958D0*HRR(8))+  & 
        2.23606797749979D0*HRR(18))+  & 
        2.23606797749979D0*HRR(30)
      OffSet=(OA+1)*LDA+(OB+7)*LDB+CDOffSet !=(3,18|
      INTGRL(OffSet)=ABz*(2.23606797749979D0*ABz*HRR(6)+  & 
        4.47213595499958D0*HRR(16))+  & 
        ABx*(ABz*(2.23606797749979D0*ABz*HRR(3)+  & 
        4.47213595499958D0*HRR(9))+  & 
        2.23606797749979D0*HRR(19))+  & 
        2.23606797749979D0*HRR(31)
      OffSet=(OA+2)*LDA+(OB+7)*LDB+CDOffSet !=(4,18|
      INTGRL(OffSet)=ABz*(2.23606797749979D0*ABz*HRR(8)+  & 
        4.47213595499958D0*HRR(18))+  & 
        ABx*(ABz*(2.23606797749979D0*ABz*HRR(4)+  & 
        4.47213595499958D0*HRR(10))+  & 
        2.23606797749979D0*HRR(20))+  & 
        2.23606797749979D0*HRR(33)
      OffSet=(OA+0)*LDA+(OB+8)*LDB+CDOffSet !=(2,19|
      INTGRL(OffSet)=ABz*(2.23606797749979D0*ABz*HRR(6)+  & 
        4.47213595499958D0*HRR(16))+  & 
        ABy*(ABz*(2.23606797749979D0*ABz*HRR(2)+  & 
        4.47213595499958D0*HRR(8))+  & 
        2.23606797749979D0*HRR(18))+  & 
        2.23606797749979D0*HRR(31)
      OffSet=(OA+1)*LDA+(OB+8)*LDB+CDOffSet !=(3,19|
      INTGRL(OffSet)=ABz*(2.23606797749979D0*ABz*HRR(7)+  & 
        4.47213595499958D0*HRR(17))+  & 
        ABy*(ABz*(2.23606797749979D0*ABz*HRR(3)+  & 
        4.47213595499958D0*HRR(9))+  & 
        2.23606797749979D0*HRR(19))+  & 
        2.23606797749979D0*HRR(32)
      OffSet=(OA+2)*LDA+(OB+8)*LDB+CDOffSet !=(4,19|
      INTGRL(OffSet)=ABz*(2.23606797749979D0*ABz*HRR(9)+  & 
        4.47213595499958D0*HRR(19))+  & 
        ABy*(ABz*(2.23606797749979D0*ABz*HRR(4)+  & 
        4.47213595499958D0*HRR(10))+  & 
        2.23606797749979D0*HRR(20))+  & 
        2.23606797749979D0*HRR(34)
      OffSet=(OA+0)*LDA+(OB+9)*LDB+CDOffSet !=(2,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(2)+  & 
        3.D0*HRR(8))+  & 
        3.D0*HRR(18))+  & 
        HRR(33)
      OffSet=(OA+1)*LDA+(OB+9)*LDB+CDOffSet !=(3,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(3)+  & 
        3.D0*HRR(9))+  & 
        3.D0*HRR(19))+  & 
        HRR(34)
      OffSet=(OA+2)*LDA+(OB+9)*LDB+CDOffSet !=(4,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(4)+  & 
        3.D0*HRR(10))+  & 
        3.D0*HRR(20))+  & 
        HRR(35)
END SUBROUTINE BraHRR310