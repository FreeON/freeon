   SUBROUTINE BraHRR103(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(11,2|
      INTGRL(OffSet)=ABx*HRR(11)+  & 
        HRR(21)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(12,2|
      INTGRL(OffSet)=2.23606797749979D0*ABx*HRR(12)+  & 
        2.23606797749979D0*HRR(22)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(13,2|
      INTGRL(OffSet)=2.23606797749979D0*ABx*HRR(13)+  & 
        2.23606797749979D0*HRR(23)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=(14,2|
      INTGRL(OffSet)=ABx*HRR(14)+  & 
        HRR(24)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=(15,2|
      INTGRL(OffSet)=2.23606797749979D0*ABx*HRR(15)+  & 
        2.23606797749979D0*HRR(26)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=(16,2|
      INTGRL(OffSet)=3.87298334620742D0*ABx*HRR(16)+  & 
        3.87298334620742D0*HRR(27)
      OffSet=(OA+6)*LDA+(OB+0)*LDB+CDOffSet !=(17,2|
      INTGRL(OffSet)=2.23606797749979D0*ABx*HRR(17)+  & 
        2.23606797749979D0*HRR(28)
      OffSet=(OA+7)*LDA+(OB+0)*LDB+CDOffSet !=(18,2|
      INTGRL(OffSet)=2.23606797749979D0*ABx*HRR(18)+  & 
        2.23606797749979D0*HRR(30)
      OffSet=(OA+8)*LDA+(OB+0)*LDB+CDOffSet !=(19,2|
      INTGRL(OffSet)=2.23606797749979D0*ABx*HRR(19)+  & 
        2.23606797749979D0*HRR(31)
      OffSet=(OA+9)*LDA+(OB+0)*LDB+CDOffSet !=(20,2|
      INTGRL(OffSet)=ABx*HRR(20)+  & 
        HRR(33)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(11,3|
      INTGRL(OffSet)=ABy*HRR(11)+  & 
        HRR(22)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=(12,3|
      INTGRL(OffSet)=2.23606797749979D0*ABy*HRR(12)+  & 
        2.23606797749979D0*HRR(23)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=(13,3|
      INTGRL(OffSet)=2.23606797749979D0*ABy*HRR(13)+  & 
        2.23606797749979D0*HRR(24)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=(14,3|
      INTGRL(OffSet)=ABy*HRR(14)+  & 
        HRR(25)
      OffSet=(OA+4)*LDA+(OB+1)*LDB+CDOffSet !=(15,3|
      INTGRL(OffSet)=2.23606797749979D0*ABy*HRR(15)+  & 
        2.23606797749979D0*HRR(27)
      OffSet=(OA+5)*LDA+(OB+1)*LDB+CDOffSet !=(16,3|
      INTGRL(OffSet)=3.87298334620742D0*ABy*HRR(16)+  & 
        3.87298334620742D0*HRR(28)
      OffSet=(OA+6)*LDA+(OB+1)*LDB+CDOffSet !=(17,3|
      INTGRL(OffSet)=2.23606797749979D0*ABy*HRR(17)+  & 
        2.23606797749979D0*HRR(29)
      OffSet=(OA+7)*LDA+(OB+1)*LDB+CDOffSet !=(18,3|
      INTGRL(OffSet)=2.23606797749979D0*ABy*HRR(18)+  & 
        2.23606797749979D0*HRR(31)
      OffSet=(OA+8)*LDA+(OB+1)*LDB+CDOffSet !=(19,3|
      INTGRL(OffSet)=2.23606797749979D0*ABy*HRR(19)+  & 
        2.23606797749979D0*HRR(32)
      OffSet=(OA+9)*LDA+(OB+1)*LDB+CDOffSet !=(20,3|
      INTGRL(OffSet)=ABy*HRR(20)+  & 
        HRR(34)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(11,4|
      INTGRL(OffSet)=ABz*HRR(11)+  & 
        HRR(26)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=(12,4|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(12)+  & 
        2.23606797749979D0*HRR(27)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=(13,4|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(13)+  & 
        2.23606797749979D0*HRR(28)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=(14,4|
      INTGRL(OffSet)=ABz*HRR(14)+  & 
        HRR(29)
      OffSet=(OA+4)*LDA+(OB+2)*LDB+CDOffSet !=(15,4|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(15)+  & 
        2.23606797749979D0*HRR(30)
      OffSet=(OA+5)*LDA+(OB+2)*LDB+CDOffSet !=(16,4|
      INTGRL(OffSet)=3.87298334620742D0*ABz*HRR(16)+  & 
        3.87298334620742D0*HRR(31)
      OffSet=(OA+6)*LDA+(OB+2)*LDB+CDOffSet !=(17,4|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(17)+  & 
        2.23606797749979D0*HRR(32)
      OffSet=(OA+7)*LDA+(OB+2)*LDB+CDOffSet !=(18,4|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(18)+  & 
        2.23606797749979D0*HRR(33)
      OffSet=(OA+8)*LDA+(OB+2)*LDB+CDOffSet !=(19,4|
      INTGRL(OffSet)=2.23606797749979D0*ABz*HRR(19)+  & 
        2.23606797749979D0*HRR(34)
      OffSet=(OA+9)*LDA+(OB+2)*LDB+CDOffSet !=(20,4|
      INTGRL(OffSet)=ABz*HRR(20)+  & 
        HRR(35)
END SUBROUTINE BraHRR103