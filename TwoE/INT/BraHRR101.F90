   SUBROUTINE BraHRR101(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(11,1|
      INTGRL(OffSet)=HRR(11)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(12,1|
      INTGRL(OffSet)=HRR(12)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(13,1|
      INTGRL(OffSet)=HRR(13)+  & 
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=(14,1|
      INTGRL(OffSet)=HRR(14)+  & 
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=(15,1|
      INTGRL(OffSet)=HRR(15)+  & 
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=(16,1|
      INTGRL(OffSet)=HRR(16)+  & 
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+0)*LDB+CDOffSet !=(17,1|
      INTGRL(OffSet)=HRR(17)+  & 
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+0)*LDB+CDOffSet !=(18,1|
      INTGRL(OffSet)=HRR(18)+  & 
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+0)*LDB+CDOffSet !=(19,1|
      INTGRL(OffSet)=HRR(19)+  & 
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+0)*LDB+CDOffSet !=(20,1|
      INTGRL(OffSet)=HRR(20)+  & 
        INTGRL(OffSet)
END SUBROUTINE BraHRR101