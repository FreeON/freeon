   SUBROUTINE BraHRR16(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(1,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(1)+  & 
        2.D0*HRR(2))+  & 
        HRR(5)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(1,6|
      INTGRL(OffSet)=1.73205080756888D0*ABy*HRR(2)+  & 
        ABx*(1.73205080756888D0*ABy*HRR(1)+  & 
        1.73205080756888D0*HRR(3))+  & 
        1.73205080756888D0*HRR(6)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(1,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(1)+  & 
        2.D0*HRR(3))+  & 
        HRR(7)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=(1,8|
      INTGRL(OffSet)=1.73205080756888D0*ABz*HRR(2)+  & 
        ABx*(1.73205080756888D0*ABz*HRR(1)+  & 
        1.73205080756888D0*HRR(4))+  & 
        1.73205080756888D0*HRR(8)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=(1,9|
      INTGRL(OffSet)=1.73205080756888D0*ABz*HRR(3)+  & 
        ABy*(1.73205080756888D0*ABz*HRR(1)+  & 
        1.73205080756888D0*HRR(4))+  & 
        1.73205080756888D0*HRR(9)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=(1,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(1)+  & 
        2.D0*HRR(4))+  & 
        HRR(10)
END SUBROUTINE BraHRR16