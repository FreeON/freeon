   SUBROUTINE BraHRR61(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(5,1|
      INTGRL(OffSet)=HRR(5)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(6,1|
      INTGRL(OffSet)=1.73205080756888D0*HRR(6)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(7,1|
      INTGRL(OffSet)=HRR(7)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=(8,1|
      INTGRL(OffSet)=1.73205080756888D0*HRR(8)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=(9,1|
      INTGRL(OffSet)=1.73205080756888D0*HRR(9)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=(10,1|
      INTGRL(OffSet)=HRR(10)
END SUBROUTINE BraHRR61