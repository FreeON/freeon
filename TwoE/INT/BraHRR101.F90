   SUBROUTINE BraHRR101(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(11,1|
      INTGRL(OffSet)=HRR(11)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(12,1|
      INTGRL(OffSet)=2.23606797749979D0*HRR(12)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(13,1|
      INTGRL(OffSet)=2.23606797749979D0*HRR(13)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=(14,1|
      INTGRL(OffSet)=HRR(14)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=(15,1|
      INTGRL(OffSet)=2.23606797749979D0*HRR(15)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=(16,1|
      INTGRL(OffSet)=3.87298334620742D0*HRR(16)
      OffSet=(OA+6)*LDA+(OB+0)*LDB+CDOffSet !=(17,1|
      INTGRL(OffSet)=2.23606797749979D0*HRR(17)
      OffSet=(OA+7)*LDA+(OB+0)*LDB+CDOffSet !=(18,1|
      INTGRL(OffSet)=2.23606797749979D0*HRR(18)
      OffSet=(OA+8)*LDA+(OB+0)*LDB+CDOffSet !=(19,1|
      INTGRL(OffSet)=2.23606797749979D0*HRR(19)
      OffSet=(OA+9)*LDA+(OB+0)*LDB+CDOffSet !=(20,1|
      INTGRL(OffSet)=HRR(20)
END SUBROUTINE BraHRR101