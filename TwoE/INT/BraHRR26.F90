   SUBROUTINE BraHRR26(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(1,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(21)+  & 
        2.D0*HRR(22))+  & 
        HRR(25)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(1,6|
      INTGRL(OffSet)=1.73205080756888D0*ABy*HRR(22)+  & 
        ABx*(1.73205080756888D0*ABy*HRR(21)+  & 
        1.73205080756888D0*HRR(23))+  & 
        1.73205080756888D0*HRR(26)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(1,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(21)+  & 
        2.D0*HRR(23))+  & 
        HRR(27)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=(1,8|
      INTGRL(OffSet)=1.73205080756888D0*ABz*HRR(22)+  & 
        ABx*(1.73205080756888D0*ABz*HRR(21)+  & 
        1.73205080756888D0*HRR(24))+  & 
        1.73205080756888D0*HRR(28)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=(1,9|
      INTGRL(OffSet)=1.73205080756888D0*ABz*HRR(23)+  & 
        ABy*(1.73205080756888D0*ABz*HRR(21)+  & 
        1.73205080756888D0*HRR(24))+  & 
        1.73205080756888D0*HRR(29)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=(1,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(21)+  & 
        2.D0*HRR(24))+  & 
        HRR(30)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(2,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(2)+  & 
        2.D0*HRR(5))+  & 
        HRR(11)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(3,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(3)+  & 
        2.D0*HRR(6))+  & 
        HRR(12)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=(4,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(4)+  & 
        2.D0*HRR(8))+  & 
        HRR(15)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=(2,6|
      INTGRL(OffSet)=1.73205080756888D0*ABy*HRR(5)+  & 
        ABx*(1.73205080756888D0*ABy*HRR(2)+  & 
        1.73205080756888D0*HRR(6))+  & 
        1.73205080756888D0*HRR(12)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=(3,6|
      INTGRL(OffSet)=1.73205080756888D0*ABy*HRR(6)+  & 
        ABx*(1.73205080756888D0*ABy*HRR(3)+  & 
        1.73205080756888D0*HRR(7))+  & 
        1.73205080756888D0*HRR(13)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=(4,6|
      INTGRL(OffSet)=1.73205080756888D0*ABy*HRR(8)+  & 
        ABx*(1.73205080756888D0*ABy*HRR(4)+  & 
        1.73205080756888D0*HRR(9))+  & 
        1.73205080756888D0*HRR(16)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=(2,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(2)+  & 
        2.D0*HRR(6))+  & 
        HRR(13)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=(3,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(3)+  & 
        2.D0*HRR(7))+  & 
        HRR(14)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=(4,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(4)+  & 
        2.D0*HRR(9))+  & 
        HRR(17)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=(2,8|
      INTGRL(OffSet)=1.73205080756888D0*ABz*HRR(5)+  & 
        ABx*(1.73205080756888D0*ABz*HRR(2)+  & 
        1.73205080756888D0*HRR(8))+  & 
        1.73205080756888D0*HRR(15)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=(3,8|
      INTGRL(OffSet)=1.73205080756888D0*ABz*HRR(6)+  & 
        ABx*(1.73205080756888D0*ABz*HRR(3)+  & 
        1.73205080756888D0*HRR(9))+  & 
        1.73205080756888D0*HRR(16)
      OffSet=(OA+3)*LDA+(OB+3)*LDB+CDOffSet !=(4,8|
      INTGRL(OffSet)=1.73205080756888D0*ABz*HRR(8)+  & 
        ABx*(1.73205080756888D0*ABz*HRR(4)+  & 
        1.73205080756888D0*HRR(10))+  & 
        1.73205080756888D0*HRR(18)
      OffSet=(OA+1)*LDA+(OB+4)*LDB+CDOffSet !=(2,9|
      INTGRL(OffSet)=1.73205080756888D0*ABz*HRR(6)+  & 
        ABy*(1.73205080756888D0*ABz*HRR(2)+  & 
        1.73205080756888D0*HRR(8))+  & 
        1.73205080756888D0*HRR(16)
      OffSet=(OA+2)*LDA+(OB+4)*LDB+CDOffSet !=(3,9|
      INTGRL(OffSet)=1.73205080756888D0*ABz*HRR(7)+  & 
        ABy*(1.73205080756888D0*ABz*HRR(3)+  & 
        1.73205080756888D0*HRR(9))+  & 
        1.73205080756888D0*HRR(17)
      OffSet=(OA+3)*LDA+(OB+4)*LDB+CDOffSet !=(4,9|
      INTGRL(OffSet)=1.73205080756888D0*ABz*HRR(9)+  & 
        ABy*(1.73205080756888D0*ABz*HRR(4)+  & 
        1.73205080756888D0*HRR(10))+  & 
        1.73205080756888D0*HRR(19)
      OffSet=(OA+1)*LDA+(OB+5)*LDB+CDOffSet !=(2,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(2)+  & 
        2.D0*HRR(8))+  & 
        HRR(18)
      OffSet=(OA+2)*LDA+(OB+5)*LDB+CDOffSet !=(3,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(3)+  & 
        2.D0*HRR(9))+  & 
        HRR(19)
      OffSet=(OA+3)*LDA+(OB+5)*LDB+CDOffSet !=(4,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(4)+  & 
        2.D0*HRR(10))+  & 
        HRR(20)
END SUBROUTINE BraHRR26