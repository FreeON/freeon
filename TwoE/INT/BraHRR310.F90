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
        HRR(21)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(3,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(3)+  & 
        3.D0*HRR(6))+  & 
        3.D0*HRR(12))+  & 
        HRR(22)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(4,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(4)+  & 
        3.D0*HRR(8))+  & 
        3.D0*HRR(15))+  & 
        HRR(26)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(2,12|
      INTGRL(OffSet)=ABy*HRR(11)+  & 
        ABx*(2.D0*ABy*HRR(5)+  & 
        ABx*(ABy*HRR(2)+  & 
        HRR(6))+  & 
        2.D0*HRR(12))+  & 
        HRR(22)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=(3,12|
      INTGRL(OffSet)=ABy*HRR(12)+  & 
        ABx*(2.D0*ABy*HRR(6)+  & 
        ABx*(ABy*HRR(3)+  & 
        HRR(7))+  & 
        2.D0*HRR(13))+  & 
        HRR(23)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=(4,12|
      INTGRL(OffSet)=ABy*HRR(15)+  & 
        ABx*(2.D0*ABy*HRR(8)+  & 
        ABx*(ABy*HRR(4)+  & 
        HRR(9))+  & 
        2.D0*HRR(16))+  & 
        HRR(27)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(2,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(5)+  & 
        2.D0*HRR(12))+  & 
        ABx*(ABy*(ABy*HRR(2)+  & 
        2.D0*HRR(6))+  & 
        HRR(13))+  & 
        HRR(23)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=(3,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(6)+  & 
        2.D0*HRR(13))+  & 
        ABx*(ABy*(ABy*HRR(3)+  & 
        2.D0*HRR(7))+  & 
        HRR(14))+  & 
        HRR(24)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=(4,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(8)+  & 
        2.D0*HRR(16))+  & 
        ABx*(ABy*(ABy*HRR(4)+  & 
        2.D0*HRR(9))+  & 
        HRR(17))+  & 
        HRR(28)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=(2,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(2)+  & 
        3.D0*HRR(6))+  & 
        3.D0*HRR(13))+  & 
        HRR(24)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=(3,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(3)+  & 
        3.D0*HRR(7))+  & 
        3.D0*HRR(14))+  & 
        HRR(25)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=(4,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(4)+  & 
        3.D0*HRR(9))+  & 
        3.D0*HRR(17))+  & 
        HRR(29)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=(2,15|
      INTGRL(OffSet)=ABz*HRR(11)+  & 
        ABx*(2.D0*ABz*HRR(5)+  & 
        ABx*(ABz*HRR(2)+  & 
        HRR(8))+  & 
        2.D0*HRR(15))+  & 
        HRR(26)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+4)*LDB+CDOffSet !=(3,15|
      INTGRL(OffSet)=ABz*HRR(12)+  & 
        ABx*(2.D0*ABz*HRR(6)+  & 
        ABx*(ABz*HRR(3)+  & 
        HRR(9))+  & 
        2.D0*HRR(16))+  & 
        HRR(27)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+4)*LDB+CDOffSet !=(4,15|
      INTGRL(OffSet)=ABz*HRR(15)+  & 
        ABx*(2.D0*ABz*HRR(8)+  & 
        ABx*(ABz*HRR(4)+  & 
        HRR(10))+  & 
        2.D0*HRR(18))+  & 
        HRR(30)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=(2,16|
      INTGRL(OffSet)=ABz*HRR(12)+  & 
        ABy*(ABz*HRR(5)+  & 
        HRR(15))+  & 
        ABx*(ABz*HRR(6)+  & 
        ABy*(ABz*HRR(2)+  & 
        HRR(8))+  & 
        HRR(16))+  & 
        HRR(27)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+5)*LDB+CDOffSet !=(3,16|
      INTGRL(OffSet)=ABz*HRR(13)+  & 
        ABy*(ABz*HRR(6)+  & 
        HRR(16))+  & 
        ABx*(ABz*HRR(7)+  & 
        ABy*(ABz*HRR(3)+  & 
        HRR(9))+  & 
        HRR(17))+  & 
        HRR(28)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+5)*LDB+CDOffSet !=(4,16|
      INTGRL(OffSet)=ABz*HRR(16)+  & 
        ABy*(ABz*HRR(8)+  & 
        HRR(18))+  & 
        ABx*(ABz*HRR(9)+  & 
        ABy*(ABz*HRR(4)+  & 
        HRR(10))+  & 
        HRR(19))+  & 
        HRR(31)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+6)*LDB+CDOffSet !=(2,17|
      INTGRL(OffSet)=ABz*HRR(13)+  & 
        ABy*(2.D0*ABz*HRR(6)+  & 
        ABy*(ABz*HRR(2)+  & 
        HRR(8))+  & 
        2.D0*HRR(16))+  & 
        HRR(28)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+6)*LDB+CDOffSet !=(3,17|
      INTGRL(OffSet)=ABz*HRR(14)+  & 
        ABy*(2.D0*ABz*HRR(7)+  & 
        ABy*(ABz*HRR(3)+  & 
        HRR(9))+  & 
        2.D0*HRR(17))+  & 
        HRR(29)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+6)*LDB+CDOffSet !=(4,17|
      INTGRL(OffSet)=ABz*HRR(17)+  & 
        ABy*(2.D0*ABz*HRR(9)+  & 
        ABy*(ABz*HRR(4)+  & 
        HRR(10))+  & 
        2.D0*HRR(19))+  & 
        HRR(32)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+7)*LDB+CDOffSet !=(2,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(5)+  & 
        2.D0*HRR(15))+  & 
        ABx*(ABz*(ABz*HRR(2)+  & 
        2.D0*HRR(8))+  & 
        HRR(18))+  & 
        HRR(30)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+7)*LDB+CDOffSet !=(3,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(6)+  & 
        2.D0*HRR(16))+  & 
        ABx*(ABz*(ABz*HRR(3)+  & 
        2.D0*HRR(9))+  & 
        HRR(19))+  & 
        HRR(31)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+7)*LDB+CDOffSet !=(4,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(8)+  & 
        2.D0*HRR(18))+  & 
        ABx*(ABz*(ABz*HRR(4)+  & 
        2.D0*HRR(10))+  & 
        HRR(20))+  & 
        HRR(33)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+8)*LDB+CDOffSet !=(2,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(6)+  & 
        2.D0*HRR(16))+  & 
        ABy*(ABz*(ABz*HRR(2)+  & 
        2.D0*HRR(8))+  & 
        HRR(18))+  & 
        HRR(31)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+8)*LDB+CDOffSet !=(3,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(7)+  & 
        2.D0*HRR(17))+  & 
        ABy*(ABz*(ABz*HRR(3)+  & 
        2.D0*HRR(9))+  & 
        HRR(19))+  & 
        HRR(32)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+8)*LDB+CDOffSet !=(4,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(9)+  & 
        2.D0*HRR(19))+  & 
        ABy*(ABz*(ABz*HRR(4)+  & 
        2.D0*HRR(10))+  & 
        HRR(20))+  & 
        HRR(34)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+9)*LDB+CDOffSet !=(2,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(2)+  & 
        3.D0*HRR(8))+  & 
        3.D0*HRR(18))+  & 
        HRR(33)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+9)*LDB+CDOffSet !=(3,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(3)+  & 
        3.D0*HRR(9))+  & 
        3.D0*HRR(19))+  & 
        HRR(34)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+9)*LDB+CDOffSet !=(4,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(4)+  & 
        3.D0*HRR(10))+  & 
        3.D0*HRR(20))+  & 
        HRR(35)+  & 
        INTGRL(OffSet)
END SUBROUTINE BraHRR310