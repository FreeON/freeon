   SUBROUTINE BraHRR31(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(2,1|
      INTGRL(OffSet)=HRR(2)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(3,1|
      INTGRL(OffSet)=HRR(3)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(4,1|
      INTGRL(OffSet)=HRR(4)
END SUBROUTINE BraHRR31