   SUBROUTINE BraHRR36(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(2,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(2)+  & 
        2.D0*HRR(5))+  & 
        HRR(11)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(3,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(3)+  & 
        2.D0*HRR(6))+  & 
        HRR(12)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(4,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(4)+  & 
        2.D0*HRR(8))+  & 
        HRR(15)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(2,6|
      INTGRL(OffSet)=ABy*HRR(5)+  & 
        ABx*(ABy*HRR(2)+  & 
        HRR(6))+  & 
        HRR(12)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=(3,6|
      INTGRL(OffSet)=ABy*HRR(6)+  & 
        ABx*(ABy*HRR(3)+  & 
        HRR(7))+  & 
        HRR(13)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=(4,6|
      INTGRL(OffSet)=ABy*HRR(8)+  & 
        ABx*(ABy*HRR(4)+  & 
        HRR(9))+  & 
        HRR(16)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(2,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(2)+  & 
        2.D0*HRR(6))+  & 
        HRR(13)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=(3,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(3)+  & 
        2.D0*HRR(7))+  & 
        HRR(14)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=(4,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(4)+  & 
        2.D0*HRR(9))+  & 
        HRR(17)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=(2,8|
      INTGRL(OffSet)=ABz*HRR(5)+  & 
        ABx*(ABz*HRR(2)+  & 
        HRR(8))+  & 
        HRR(15)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=(3,8|
      INTGRL(OffSet)=ABz*HRR(6)+  & 
        ABx*(ABz*HRR(3)+  & 
        HRR(9))+  & 
        HRR(16)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=(4,8|
      INTGRL(OffSet)=ABz*HRR(8)+  & 
        ABx*(ABz*HRR(4)+  & 
        HRR(10))+  & 
        HRR(18)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=(2,9|
      INTGRL(OffSet)=ABz*HRR(6)+  & 
        ABy*(ABz*HRR(2)+  & 
        HRR(8))+  & 
        HRR(16)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+4)*LDB+CDOffSet !=(3,9|
      INTGRL(OffSet)=ABz*HRR(7)+  & 
        ABy*(ABz*HRR(3)+  & 
        HRR(9))+  & 
        HRR(17)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+4)*LDB+CDOffSet !=(4,9|
      INTGRL(OffSet)=ABz*HRR(9)+  & 
        ABy*(ABz*HRR(4)+  & 
        HRR(10))+  & 
        HRR(19)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=(2,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(2)+  & 
        2.D0*HRR(8))+  & 
        HRR(18)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+5)*LDB+CDOffSet !=(3,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(3)+  & 
        2.D0*HRR(9))+  & 
        HRR(19)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+5)*LDB+CDOffSet !=(4,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(4)+  & 
        2.D0*HRR(10))+  & 
        HRR(20)+  & 
        INTGRL(OffSet)
END SUBROUTINE BraHRR36