   SUBROUTINE BraHRR61(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(5,1|
      INTGRL(OffSet)=HRR(5)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(6,1|
      INTGRL(OffSet)=HRR(6)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(7,1|
      INTGRL(OffSet)=HRR(7)+  & 
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=(8,1|
      INTGRL(OffSet)=HRR(8)+  & 
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=(9,1|
      INTGRL(OffSet)=HRR(9)+  & 
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=(10,1|
      INTGRL(OffSet)=HRR(10)+  & 
        INTGRL(OffSet)
END SUBROUTINE BraHRR61