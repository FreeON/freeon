   SUBROUTINE BraHRR12(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(1,1|
      INTGRL(OffSet)=HRR(5)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(1,2|
      INTGRL(OffSet)=ABx*HRR(1)+  & 
        HRR(2)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(1,3|
      INTGRL(OffSet)=ABy*HRR(1)+  & 
        HRR(3)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=(1,4|
      INTGRL(OffSet)=ABz*HRR(1)+  & 
        HRR(4)
END SUBROUTINE BraHRR12