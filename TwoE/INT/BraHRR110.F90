   SUBROUTINE BraHRR110(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(1,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(1)+  & 
        3.D0*HRR(2))+  & 
        3.D0*HRR(5))+  & 
        HRR(11)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(1,12|
      INTGRL(OffSet)=2.23606797749979D0*ABy*HRR(5)+  & 
        ABx*(4.47213595499958D0*ABy*HRR(2)+  & 
        ABx*(2.23606797749979D0*ABy*HRR(1)+  & 
        2.23606797749979D0*HRR(3))+  & 
        4.47213595499958D0*HRR(6))+  & 
        2.23606797749979D0*HRR(12)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(1,13|
      INTGRL(OffSet)=ABy*(2.23606797749979D0*ABy*HRR(2)+  & 
        4.47213595499958D0*HRR(6))+  & 
        ABx*(ABy*(2.23606797749979D0*ABy*HRR(1)+  & 
        4.47213595499958D0*HRR(3))+  & 
        2.23606797749979D0*HRR(7))+  & 
        2.23606797749979D0*HRR(13)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=(1,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(1)+  & 
        3.D0*HRR(3))+  & 
        3.D0*HRR(7))+  & 
        HRR(14)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=(1,15|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(5)+  & 
        ABx*(4.47213595499958D0*ABz*HRR(2)+  & 
        ABx*(2.23606797749979D0*ABz*HRR(1)+  & 
        2.23606797749979D0*HRR(4))+  & 
        4.47213595499958D0*HRR(8))+  & 
        2.23606797749979D0*HRR(15)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=(1,16|
      INTGRL(OffSet)=3.87298334620742D0*ABz*HRR(6)+  & 
        ABy*(3.87298334620742D0*ABz*HRR(2)+  & 
        3.87298334620742D0*HRR(8))+  & 
        ABx*(3.87298334620742D0*ABz*HRR(3)+  & 
        ABy*(3.87298334620742D0*ABz*HRR(1)+  & 
        3.87298334620742D0*HRR(4))+  & 
        3.87298334620742D0*HRR(9))+  & 
        3.87298334620742D0*HRR(16)
      OffSet=(OA+0)*LDA+(OB+6)*LDB+CDOffSet !=(1,17|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(7)+  & 
        ABy*(4.47213595499958D0*ABz*HRR(3)+  & 
        ABy*(2.23606797749979D0*ABz*HRR(1)+  & 
        2.23606797749979D0*HRR(4))+  & 
        4.47213595499958D0*HRR(9))+  & 
        2.23606797749979D0*HRR(17)
      OffSet=(OA+0)*LDA+(OB+7)*LDB+CDOffSet !=(1,18|
      INTGRL(OffSet)=ABz*(2.23606797749979D0*ABz*HRR(2)+  & 
        4.47213595499958D0*HRR(8))+  & 
        ABx*(ABz*(2.23606797749979D0*ABz*HRR(1)+  & 
        4.47213595499958D0*HRR(4))+  & 
        2.23606797749979D0*HRR(10))+  & 
        2.23606797749979D0*HRR(18)
      OffSet=(OA+0)*LDA+(OB+8)*LDB+CDOffSet !=(1,19|
      INTGRL(OffSet)=ABz*(2.23606797749979D0*ABz*HRR(3)+  & 
        4.47213595499958D0*HRR(9))+  & 
        ABy*(ABz*(2.23606797749979D0*ABz*HRR(1)+  & 
        4.47213595499958D0*HRR(4))+  & 
        2.23606797749979D0*HRR(10))+  & 
        2.23606797749979D0*HRR(19)
      OffSet=(OA+0)*LDA+(OB+9)*LDB+CDOffSet !=(1,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(1)+  & 
        3.D0*HRR(4))+  & 
        3.D0*HRR(10))+  & 
        HRR(20)
END SUBROUTINE BraHRR110