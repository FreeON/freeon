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
        HRR(11)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(1,12|
      INTGRL(OffSet)=ABy*HRR(5)+  & 
        ABx*(2.D0*ABy*HRR(2)+  & 
        ABx*(ABy*HRR(1)+  & 
        HRR(3))+  & 
        2.D0*HRR(6))+  & 
        HRR(12)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(1,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(2)+  & 
        2.D0*HRR(6))+  & 
        ABx*(ABy*(ABy*HRR(1)+  & 
        2.D0*HRR(3))+  & 
        HRR(7))+  & 
        HRR(13)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=(1,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(1)+  & 
        3.D0*HRR(3))+  & 
        3.D0*HRR(7))+  & 
        HRR(14)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=(1,15|
      INTGRL(OffSet)=ABz*HRR(5)+  & 
        ABx*(2.D0*ABz*HRR(2)+  & 
        ABx*(ABz*HRR(1)+  & 
        HRR(4))+  & 
        2.D0*HRR(8))+  & 
        HRR(15)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=(1,16|
      INTGRL(OffSet)=ABz*HRR(6)+  & 
        ABy*(ABz*HRR(2)+  & 
        HRR(8))+  & 
        ABx*(ABz*HRR(3)+  & 
        ABy*(ABz*HRR(1)+  & 
        HRR(4))+  & 
        HRR(9))+  & 
        HRR(16)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+6)*LDB+CDOffSet !=(1,17|
      INTGRL(OffSet)=ABz*HRR(7)+  & 
        ABy*(2.D0*ABz*HRR(3)+  & 
        ABy*(ABz*HRR(1)+  & 
        HRR(4))+  & 
        2.D0*HRR(9))+  & 
        HRR(17)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+7)*LDB+CDOffSet !=(1,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(2)+  & 
        2.D0*HRR(8))+  & 
        ABx*(ABz*(ABz*HRR(1)+  & 
        2.D0*HRR(4))+  & 
        HRR(10))+  & 
        HRR(18)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+8)*LDB+CDOffSet !=(1,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(3)+  & 
        2.D0*HRR(9))+  & 
        ABy*(ABz*(ABz*HRR(1)+  & 
        2.D0*HRR(4))+  & 
        HRR(10))+  & 
        HRR(19)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+9)*LDB+CDOffSet !=(1,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(1)+  & 
        3.D0*HRR(4))+  & 
        3.D0*HRR(10))+  & 
        HRR(20)+  & 
        INTGRL(OffSet)
END SUBROUTINE BraHRR110