   SUBROUTINE BraHRR66(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(5,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(5)+  & 
        2.D0*HRR(11))+  & 
        HRR(21)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(6,5|
      INTGRL(OffSet)=ABx*(1.73205080756888D0*ABx*HRR(6)+  & 
        3.46410161513775D0*HRR(12))+  & 
        1.73205080756888D0*HRR(22)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(7,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(7)+  & 
        2.D0*HRR(13))+  & 
        HRR(23)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=(8,5|
      INTGRL(OffSet)=ABx*(1.73205080756888D0*ABx*HRR(8)+  & 
        3.46410161513775D0*HRR(15))+  & 
        1.73205080756888D0*HRR(26)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=(9,5|
      INTGRL(OffSet)=ABx*(1.73205080756888D0*ABx*HRR(9)+  & 
        3.46410161513775D0*HRR(16))+  & 
        1.73205080756888D0*HRR(27)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=(10,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(10)+  & 
        2.D0*HRR(18))+  & 
        HRR(30)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(5,6|
      INTGRL(OffSet)=1.73205080756888D0*ABy*HRR(11)+  & 
        ABx*(1.73205080756888D0*ABy*HRR(5)+  & 
        1.73205080756888D0*HRR(12))+  & 
        1.73205080756888D0*HRR(22)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=(6,6|
      INTGRL(OffSet)=3.D0*ABy*HRR(12)+  & 
        ABx*(3.D0*ABy*HRR(6)+  & 
        3.D0*HRR(13))+  & 
        3.D0*HRR(23)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=(7,6|
      INTGRL(OffSet)=1.73205080756888D0*ABy*HRR(13)+  & 
        ABx*(1.73205080756888D0*ABy*HRR(7)+  & 
        1.73205080756888D0*HRR(14))+  & 
        1.73205080756888D0*HRR(24)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=(8,6|
      INTGRL(OffSet)=3.D0*ABy*HRR(15)+  & 
        ABx*(3.D0*ABy*HRR(8)+  & 
        3.D0*HRR(16))+  & 
        3.D0*HRR(27)
      OffSet=(OA+4)*LDA+(OB+1)*LDB+CDOffSet !=(9,6|
      INTGRL(OffSet)=3.D0*ABy*HRR(16)+  & 
        ABx*(3.D0*ABy*HRR(9)+  & 
        3.D0*HRR(17))+  & 
        3.D0*HRR(28)
      OffSet=(OA+5)*LDA+(OB+1)*LDB+CDOffSet !=(10,6|
      INTGRL(OffSet)=1.73205080756888D0*ABy*HRR(18)+  & 
        ABx*(1.73205080756888D0*ABy*HRR(10)+  & 
        1.73205080756888D0*HRR(19))+  & 
        1.73205080756888D0*HRR(31)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(5,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(5)+  & 
        2.D0*HRR(12))+  & 
        HRR(23)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=(6,7|
      INTGRL(OffSet)=ABy*(1.73205080756888D0*ABy*HRR(6)+  & 
        3.46410161513775D0*HRR(13))+  & 
        1.73205080756888D0*HRR(24)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=(7,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(7)+  & 
        2.D0*HRR(14))+  & 
        HRR(25)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=(8,7|
      INTGRL(OffSet)=ABy*(1.73205080756888D0*ABy*HRR(8)+  & 
        3.46410161513775D0*HRR(16))+  & 
        1.73205080756888D0*HRR(28)
      OffSet=(OA+4)*LDA+(OB+2)*LDB+CDOffSet !=(9,7|
      INTGRL(OffSet)=ABy*(1.73205080756888D0*ABy*HRR(9)+  & 
        3.46410161513775D0*HRR(17))+  & 
        1.73205080756888D0*HRR(29)
      OffSet=(OA+5)*LDA+(OB+2)*LDB+CDOffSet !=(10,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(10)+  & 
        2.D0*HRR(19))+  & 
        HRR(32)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=(5,8|
      INTGRL(OffSet)=1.73205080756888D0*ABz*HRR(11)+  & 
        ABx*(1.73205080756888D0*ABz*HRR(5)+  & 
        1.73205080756888D0*HRR(15))+  & 
        1.73205080756888D0*HRR(26)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=(6,8|
      INTGRL(OffSet)=3.D0*ABz*HRR(12)+  & 
        ABx*(3.D0*ABz*HRR(6)+  & 
        3.D0*HRR(16))+  & 
        3.D0*HRR(27)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=(7,8|
      INTGRL(OffSet)=1.73205080756888D0*ABz*HRR(13)+  & 
        ABx*(1.73205080756888D0*ABz*HRR(7)+  & 
        1.73205080756888D0*HRR(17))+  & 
        1.73205080756888D0*HRR(28)
      OffSet=(OA+3)*LDA+(OB+3)*LDB+CDOffSet !=(8,8|
      INTGRL(OffSet)=3.D0*ABz*HRR(15)+  & 
        ABx*(3.D0*ABz*HRR(8)+  & 
        3.D0*HRR(18))+  & 
        3.D0*HRR(30)
      OffSet=(OA+4)*LDA+(OB+3)*LDB+CDOffSet !=(9,8|
      INTGRL(OffSet)=3.D0*ABz*HRR(16)+  & 
        ABx*(3.D0*ABz*HRR(9)+  & 
        3.D0*HRR(19))+  & 
        3.D0*HRR(31)
      OffSet=(OA+5)*LDA+(OB+3)*LDB+CDOffSet !=(10,8|
      INTGRL(OffSet)=1.73205080756888D0*ABz*HRR(18)+  & 
        ABx*(1.73205080756888D0*ABz*HRR(10)+  & 
        1.73205080756888D0*HRR(20))+  & 
        1.73205080756888D0*HRR(33)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=(5,9|
      INTGRL(OffSet)=1.73205080756888D0*ABz*HRR(12)+  & 
        ABy*(1.73205080756888D0*ABz*HRR(5)+  & 
        1.73205080756888D0*HRR(15))+  & 
        1.73205080756888D0*HRR(27)
      OffSet=(OA+1)*LDA+(OB+4)*LDB+CDOffSet !=(6,9|
      INTGRL(OffSet)=3.D0*ABz*HRR(13)+  & 
        ABy*(3.D0*ABz*HRR(6)+  & 
        3.D0*HRR(16))+  & 
        3.D0*HRR(28)
      OffSet=(OA+2)*LDA+(OB+4)*LDB+CDOffSet !=(7,9|
      INTGRL(OffSet)=1.73205080756888D0*ABz*HRR(14)+  & 
        ABy*(1.73205080756888D0*ABz*HRR(7)+  & 
        1.73205080756888D0*HRR(17))+  & 
        1.73205080756888D0*HRR(29)
      OffSet=(OA+3)*LDA+(OB+4)*LDB+CDOffSet !=(8,9|
      INTGRL(OffSet)=3.D0*ABz*HRR(16)+  & 
        ABy*(3.D0*ABz*HRR(8)+  & 
        3.D0*HRR(18))+  & 
        3.D0*HRR(31)
      OffSet=(OA+4)*LDA+(OB+4)*LDB+CDOffSet !=(9,9|
      INTGRL(OffSet)=3.D0*ABz*HRR(17)+  & 
        ABy*(3.D0*ABz*HRR(9)+  & 
        3.D0*HRR(19))+  & 
        3.D0*HRR(32)
      OffSet=(OA+5)*LDA+(OB+4)*LDB+CDOffSet !=(10,9|
      INTGRL(OffSet)=1.73205080756888D0*ABz*HRR(19)+  & 
        ABy*(1.73205080756888D0*ABz*HRR(10)+  & 
        1.73205080756888D0*HRR(20))+  & 
        1.73205080756888D0*HRR(34)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=(5,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(5)+  & 
        2.D0*HRR(15))+  & 
        HRR(30)
      OffSet=(OA+1)*LDA+(OB+5)*LDB+CDOffSet !=(6,10|
      INTGRL(OffSet)=ABz*(1.73205080756888D0*ABz*HRR(6)+  & 
        3.46410161513775D0*HRR(16))+  & 
        1.73205080756888D0*HRR(31)
      OffSet=(OA+2)*LDA+(OB+5)*LDB+CDOffSet !=(7,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(7)+  & 
        2.D0*HRR(17))+  & 
        HRR(32)
      OffSet=(OA+3)*LDA+(OB+5)*LDB+CDOffSet !=(8,10|
      INTGRL(OffSet)=ABz*(1.73205080756888D0*ABz*HRR(8)+  & 
        3.46410161513775D0*HRR(18))+  & 
        1.73205080756888D0*HRR(33)
      OffSet=(OA+4)*LDA+(OB+5)*LDB+CDOffSet !=(9,10|
      INTGRL(OffSet)=ABz*(1.73205080756888D0*ABz*HRR(9)+  & 
        3.46410161513775D0*HRR(19))+  & 
        1.73205080756888D0*HRR(34)
      OffSet=(OA+5)*LDA+(OB+5)*LDB+CDOffSet !=(10,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(10)+  & 
        2.D0*HRR(20))+  & 
        HRR(35)
END SUBROUTINE BraHRR66