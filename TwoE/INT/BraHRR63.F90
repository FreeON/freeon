   SUBROUTINE BraHRR63(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(5,2|
      INTGRL(OffSet)=ABx*HRR(5)+  & 
        HRR(11)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(6,2|
      INTGRL(OffSet)=1.73205080756888D0*ABx*HRR(6)+  & 
        1.73205080756888D0*HRR(12)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(7,2|
      INTGRL(OffSet)=ABx*HRR(7)+  & 
        HRR(13)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=(8,2|
      INTGRL(OffSet)=1.73205080756888D0*ABx*HRR(8)+  & 
        1.73205080756888D0*HRR(15)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=(9,2|
      INTGRL(OffSet)=1.73205080756888D0*ABx*HRR(9)+  & 
        1.73205080756888D0*HRR(16)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=(10,2|
      INTGRL(OffSet)=ABx*HRR(10)+  & 
        HRR(18)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(5,3|
      INTGRL(OffSet)=ABy*HRR(5)+  & 
        HRR(12)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=(6,3|
      INTGRL(OffSet)=1.73205080756888D0*ABy*HRR(6)+  & 
        1.73205080756888D0*HRR(13)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=(7,3|
      INTGRL(OffSet)=ABy*HRR(7)+  & 
        HRR(14)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=(8,3|
      INTGRL(OffSet)=1.73205080756888D0*ABy*HRR(8)+  & 
        1.73205080756888D0*HRR(16)
      OffSet=(OA+4)*LDA+(OB+1)*LDB+CDOffSet !=(9,3|
      INTGRL(OffSet)=1.73205080756888D0*ABy*HRR(9)+  & 
        1.73205080756888D0*HRR(17)
      OffSet=(OA+5)*LDA+(OB+1)*LDB+CDOffSet !=(10,3|
      INTGRL(OffSet)=ABy*HRR(10)+  & 
        HRR(19)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(5,4|
      INTGRL(OffSet)=ABz*HRR(5)+  & 
        HRR(15)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=(6,4|
      INTGRL(OffSet)=1.73205080756888D0*ABz*HRR(6)+  & 
        1.73205080756888D0*HRR(16)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=(7,4|
      INTGRL(OffSet)=ABz*HRR(7)+  & 
        HRR(17)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=(8,4|
      INTGRL(OffSet)=1.73205080756888D0*ABz*HRR(8)+  & 
        1.73205080756888D0*HRR(18)
      OffSet=(OA+4)*LDA+(OB+2)*LDB+CDOffSet !=(9,4|
      INTGRL(OffSet)=1.73205080756888D0*ABz*HRR(9)+  & 
        1.73205080756888D0*HRR(19)
      OffSet=(OA+5)*LDA+(OB+2)*LDB+CDOffSet !=(10,4|
      INTGRL(OffSet)=ABz*HRR(10)+  & 
        HRR(20)
END SUBROUTINE BraHRR63