   SUBROUTINE BraHRR23(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(1,2|
      INTGRL(OffSet)=ABx*HRR(11)+  & 
        HRR(12)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(1,3|
      INTGRL(OffSet)=ABy*HRR(11)+  & 
        HRR(13)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(1,4|
      INTGRL(OffSet)=ABz*HRR(11)+  & 
        HRR(14)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(2,2|
      INTGRL(OffSet)=ABx*HRR(2)+  & 
        HRR(5)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(3,2|
      INTGRL(OffSet)=ABx*HRR(3)+  & 
        HRR(6)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=(4,2|
      INTGRL(OffSet)=ABx*HRR(4)+  & 
        HRR(8)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=(2,3|
      INTGRL(OffSet)=ABy*HRR(2)+  & 
        HRR(6)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=(3,3|
      INTGRL(OffSet)=ABy*HRR(3)+  & 
        HRR(7)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=(4,3|
      INTGRL(OffSet)=ABy*HRR(4)+  & 
        HRR(9)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=(2,4|
      INTGRL(OffSet)=ABz*HRR(2)+  & 
        HRR(8)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=(3,4|
      INTGRL(OffSet)=ABz*HRR(3)+  & 
        HRR(9)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=(4,4|
      INTGRL(OffSet)=ABz*HRR(4)+  & 
        HRR(10)
END SUBROUTINE BraHRR23