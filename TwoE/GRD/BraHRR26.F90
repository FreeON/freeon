    SUBROUTINE BraHRR26ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR,HRRA,HRRB,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*),HRRA(*),HRRB(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(1_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*(ABx*HRRA(2)+&
                         2.D0*HRRA(5))+&
                         HRRA(11)
      !=(1,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(2)+&
                         ABx*(-2.D0*HRR(1)+&
                         ABx*(ABx*HRRB(1)+&
                         3.D0*HRRB(2))+&
                         3.D0*HRRB(5))+&
                         HRRB(11)
      !=(1_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*(ABx*HRRA(3)+&
                         2.D0*HRRA(6))+&
                         HRRA(12)
      !=(1,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(5)+&
                         ABx*(2.D0*ABy*HRRB(2)+&
                         ABx*(ABy*HRRB(1)+&
                         HRRB(3))+&
                         2.D0*HRRB(6))+&
                         HRRB(12)
      !=(1_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*(ABx*HRRA(4)+&
                         2.D0*HRRA(8))+&
                         HRRA(15)
      !=(1,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(5)+&
                         ABx*(2.D0*ABz*HRRB(2)+&
                         ABx*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         2.D0*HRRB(8))+&
                         HRRB(15)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(1_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(5)+&
                         ABx*(ABy*HRRA(2)+&
                         HRRA(6))+&
                         HRRA(12)
      !=(1,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRB(5))+&
                         ABx*(2.D0*ABy*HRRB(2)+&
                         ABx*(ABy*HRRB(1)+&
                         HRRB(3))+&
                         2.D0*HRRB(6))+&
                         HRRB(12)
      !=(1_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(6)+&
                         ABx*(ABy*HRRA(3)+&
                         HRRA(7))+&
                         HRRA(13)
      !=(1,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(2)+&
                         ABy*(ABy*HRRB(2)+&
                         2.D0*HRRB(6))+&
                         ABx*(-HRR(1)+&
                         ABy*(ABy*HRRB(1)+&
                         2.D0*HRRB(3))+&
                         HRRB(7))+&
                         HRRB(13)
      !=(1_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(8)+&
                         ABx*(ABy*HRRA(4)+&
                         HRRA(9))+&
                         HRRA(16)
      !=(1,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         ABx*(ABz*HRRB(3)+&
                         ABy*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         HRRB(9))+&
                         HRRB(16)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(1_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*(ABy*HRRA(2)+&
                         2.D0*HRRA(6))+&
                         HRRA(13)
      !=(1,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(2)+&
                         2.D0*HRRB(6))+&
                         ABx*(ABy*(ABy*HRRB(1)+&
                         2.D0*HRRB(3))+&
                         HRRB(7))+&
                         HRRB(13)
      !=(1_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*(ABy*HRRA(3)+&
                         2.D0*HRRA(7))+&
                         HRRA(14)
      !=(1,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(3)+&
                         ABy*(-2.D0*HRR(1)+&
                         ABy*(ABy*HRRB(1)+&
                         3.D0*HRRB(3))+&
                         3.D0*HRRB(7))+&
                         HRRB(14)
      !=(1_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*(ABy*HRRA(4)+&
                         2.D0*HRRA(9))+&
                         HRRA(17)
      !=(1,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(7)+&
                         ABy*(2.D0*ABz*HRRB(3)+&
                         ABy*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         2.D0*HRRB(9))+&
                         HRRB(17)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(1_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(5)+&
                         ABx*(ABz*HRRA(2)+&
                         HRRA(8))+&
                         HRRA(15)
      !=(1,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRB(5))+&
                         ABx*(2.D0*ABz*HRRB(2)+&
                         ABx*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         2.D0*HRRB(8))+&
                         HRRB(15)
      !=(1_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(6)+&
                         ABx*(ABz*HRRA(3)+&
                         HRRA(9))+&
                         HRRA(16)
      !=(1,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         ABx*(ABz*HRRB(3)+&
                         ABy*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         HRRB(9))+&
                         HRRB(16)
      !=(1_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(8)+&
                         ABx*(ABz*HRRA(4)+&
                         HRRA(10))+&
                         HRRA(18)
      !=(1,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(2)+&
                         ABz*(ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         ABx*(-HRR(1)+&
                         ABz*(ABz*HRRB(1)+&
                         2.D0*HRRB(4))+&
                         HRRB(10))+&
                         HRRB(18)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(1_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(6)+&
                         ABy*(ABz*HRRA(2)+&
                         HRRA(8))+&
                         HRRA(16)
      !=(1,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         ABx*(ABz*HRRB(3)+&
                         ABy*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         HRRB(9))+&
                         HRRB(16)
      !=(1_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(7)+&
                         ABy*(ABz*HRRA(3)+&
                         HRRA(9))+&
                         HRRA(17)
      !=(1,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRB(7))+&
                         ABy*(2.D0*ABz*HRRB(3)+&
                         ABy*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         2.D0*HRRB(9))+&
                         HRRB(17)
      !=(1_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(9)+&
                         ABy*(ABz*HRRA(4)+&
                         HRRA(10))+&
                         HRRA(19)
      !=(1,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(3)+&
                         ABz*(ABz*HRRB(3)+&
                         2.D0*HRRB(9))+&
                         ABy*(-HRR(1)+&
                         ABz*(ABz*HRRB(1)+&
                         2.D0*HRRB(4))+&
                         HRRB(10))+&
                         HRRB(19)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(1_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*(ABz*HRRA(2)+&
                         2.D0*HRRA(8))+&
                         HRRA(18)
      !=(1,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         ABx*(ABz*(ABz*HRRB(1)+&
                         2.D0*HRRB(4))+&
                         HRRB(10))+&
                         HRRB(18)
      !=(1_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*(ABz*HRRA(3)+&
                         2.D0*HRRA(9))+&
                         HRRA(19)
      !=(1,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(3)+&
                         2.D0*HRRB(9))+&
                         ABy*(ABz*(ABz*HRRB(1)+&
                         2.D0*HRRB(4))+&
                         HRRB(10))+&
                         HRRB(19)
      !=(1_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*(ABz*HRRA(4)+&
                         2.D0*HRRA(10))+&
                         HRRA(20)
      !=(1,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(4)+&
                         ABz*(-2.D0*HRR(1)+&
                         ABz*(ABz*HRRB(1)+&
                         3.D0*HRRB(4))+&
                         3.D0*HRRB(10))+&
                         HRRB(20)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(2_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(5)+&
                         ABx*(-2.D0*HRR(2)+&
                         ABx*(-HRR(1)+&
                         HRRA(5))+&
                         2.D0*HRRA(11))+&
                         HRRA(21)
      !=(2,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(5)+&
                         ABx*(-2.D0*HRR(2)+&
                         ABx*(ABx*HRRB(2)+&
                         3.D0*HRRB(5))+&
                         3.D0*HRRB(11))+&
                         HRRB(21)
      !=(2_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*(ABx*HRRA(6)+&
                         2.D0*HRRA(12))+&
                         HRRA(22)
      !=(2,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(11)+&
                         ABx*(2.D0*ABy*HRRB(5)+&
                         ABx*(ABy*HRRB(2)+&
                         HRRB(6))+&
                         2.D0*HRRB(12))+&
                         HRRB(22)
      !=(2_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*(ABx*HRRA(8)+&
                         2.D0*HRRA(15))+&
                         HRRA(26)
      !=(2,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(11)+&
                         ABx*(2.D0*ABz*HRRB(5)+&
                         ABx*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         2.D0*HRRB(15))+&
                         HRRB(26)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(3_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*(ABx*HRRA(6)+&
                         2.D0*HRRA(12))+&
                         HRRA(22)
      !=(3,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(6)+&
                         ABx*(-2.D0*HRR(3)+&
                         ABx*(ABx*HRRB(3)+&
                         3.D0*HRRB(6))+&
                         3.D0*HRRB(12))+&
                         HRRB(22)
      !=(3_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(5)+&
                         ABx*(-2.D0*HRR(2)+&
                         ABx*(-HRR(1)+&
                         HRRA(7))+&
                         2.D0*HRRA(13))+&
                         HRRA(23)
      !=(3,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(12)+&
                         ABx*(2.D0*ABy*HRRB(6)+&
                         ABx*(ABy*HRRB(3)+&
                         HRRB(7))+&
                         2.D0*HRRB(13))+&
                         HRRB(23)
      !=(3_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*(ABx*HRRA(9)+&
                         2.D0*HRRA(16))+&
                         HRRA(27)
      !=(3,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(12)+&
                         ABx*(2.D0*ABz*HRRB(6)+&
                         ABx*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         2.D0*HRRB(16))+&
                         HRRB(27)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(4_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*(ABx*HRRA(8)+&
                         2.D0*HRRA(15))+&
                         HRRA(26)
      !=(4,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(8)+&
                         ABx*(-2.D0*HRR(4)+&
                         ABx*(ABx*HRRB(4)+&
                         3.D0*HRRB(8))+&
                         3.D0*HRRB(15))+&
                         HRRB(26)
      !=(4_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*(ABx*HRRA(9)+&
                         2.D0*HRRA(16))+&
                         HRRA(27)
      !=(4,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(15)+&
                         ABx*(2.D0*ABy*HRRB(8)+&
                         ABx*(ABy*HRRB(4)+&
                         HRRB(9))+&
                         2.D0*HRRB(16))+&
                         HRRB(27)
      !=(4_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(5)+&
                         ABx*(-2.D0*HRR(2)+&
                         ABx*(-HRR(1)+&
                         HRRA(10))+&
                         2.D0*HRRA(18))+&
                         HRRA(30)
      !=(4,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(15)+&
                         ABx*(2.D0*ABz*HRRB(8)+&
                         ABx*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         2.D0*HRRB(18))+&
                         HRRB(30)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(2_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(6)+&
                         ABy*(-HRR(2)+&
                         HRRA(11))+&
                         ABx*(-HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(5))+&
                         HRRA(12))+&
                         HRRA(22)
      !=(2,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(6)+&
                         ABy*(-HRR(2)+&
                         HRRB(11))+&
                         ABx*(2.D0*ABy*HRRB(5)+&
                         ABx*(ABy*HRRB(2)+&
                         HRRB(6))+&
                         2.D0*HRRB(12))+&
                         HRRB(22)
      !=(2_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(12)+&
                         ABx*(ABy*HRRA(6)+&
                         HRRA(13))+&
                         HRRA(23)
      !=(2,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(5)+&
                         ABy*(ABy*HRRB(5)+&
                         2.D0*HRRB(12))+&
                         ABx*(-HRR(2)+&
                         ABy*(ABy*HRRB(2)+&
                         2.D0*HRRB(6))+&
                         HRRB(13))+&
                         HRRB(23)
      !=(2_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(15)+&
                         ABx*(ABy*HRRA(8)+&
                         HRRA(16))+&
                         HRRA(27)
      !=(2,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(12)+&
                         ABy*(ABz*HRRB(5)+&
                         HRRB(15))+&
                         ABx*(ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         HRRB(16))+&
                         HRRB(27)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(3_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(12)+&
                         ABx*(ABy*HRRA(6)+&
                         HRRA(13))+&
                         HRRA(23)
      !=(3,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(7)+&
                         ABy*(-HRR(3)+&
                         HRRB(12))+&
                         ABx*(2.D0*ABy*HRRB(6)+&
                         ABx*(ABy*HRRB(3)+&
                         HRRB(7))+&
                         2.D0*HRRB(13))+&
                         HRRB(23)
      !=(3_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(6)+&
                         ABy*(-HRR(2)+&
                         HRRA(13))+&
                         ABx*(-HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(7))+&
                         HRRA(14))+&
                         HRRA(24)
      !=(3,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(6)+&
                         ABy*(ABy*HRRB(6)+&
                         2.D0*HRRB(13))+&
                         ABx*(-HRR(3)+&
                         ABy*(ABy*HRRB(3)+&
                         2.D0*HRRB(7))+&
                         HRRB(14))+&
                         HRRB(24)
      !=(3_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(16)+&
                         ABx*(ABy*HRRA(9)+&
                         HRRA(17))+&
                         HRRA(28)
      !=(3,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(13)+&
                         ABy*(ABz*HRRB(6)+&
                         HRRB(16))+&
                         ABx*(ABz*HRRB(7)+&
                         ABy*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         HRRB(17))+&
                         HRRB(28)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(4_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(15)+&
                         ABx*(ABy*HRRA(8)+&
                         HRRA(16))+&
                         HRRA(27)
      !=(4,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(9)+&
                         ABy*(-HRR(4)+&
                         HRRB(15))+&
                         ABx*(2.D0*ABy*HRRB(8)+&
                         ABx*(ABy*HRRB(4)+&
                         HRRB(9))+&
                         2.D0*HRRB(16))+&
                         HRRB(27)
      !=(4_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(16)+&
                         ABx*(ABy*HRRA(9)+&
                         HRRA(17))+&
                         HRRA(28)
      !=(4,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(8)+&
                         ABy*(ABy*HRRB(8)+&
                         2.D0*HRRB(16))+&
                         ABx*(-HRR(4)+&
                         ABy*(ABy*HRRB(4)+&
                         2.D0*HRRB(9))+&
                         HRRB(17))+&
                         HRRB(28)
      !=(4_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(6)+&
                         ABy*(-HRR(2)+&
                         HRRA(18))+&
                         ABx*(-HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(10))+&
                         HRRA(19))+&
                         HRRA(31)
      !=(4,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(16)+&
                         ABy*(ABz*HRRB(8)+&
                         HRRB(18))+&
                         ABx*(ABz*HRRB(9)+&
                         ABy*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         HRRB(19))+&
                         HRRB(31)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(2_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(7)+&
                         ABy*(-2.D0*HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(5))+&
                         2.D0*HRRA(12))+&
                         HRRA(23)
      !=(2,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(5)+&
                         2.D0*HRRB(12))+&
                         ABx*(ABy*(ABy*HRRB(2)+&
                         2.D0*HRRB(6))+&
                         HRRB(13))+&
                         HRRB(23)
      !=(2_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*(ABy*HRRA(6)+&
                         2.D0*HRRA(13))+&
                         HRRA(24)
      !=(2,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(6)+&
                         ABy*(-2.D0*HRR(2)+&
                         ABy*(ABy*HRRB(2)+&
                         3.D0*HRRB(6))+&
                         3.D0*HRRB(13))+&
                         HRRB(24)
      !=(2_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*(ABy*HRRA(8)+&
                         2.D0*HRRA(16))+&
                         HRRA(28)
      !=(2,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(13)+&
                         ABy*(2.D0*ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         2.D0*HRRB(16))+&
                         HRRB(28)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(3_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*(ABy*HRRA(6)+&
                         2.D0*HRRA(13))+&
                         HRRA(24)
      !=(3,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(6)+&
                         2.D0*HRRB(13))+&
                         ABx*(ABy*(ABy*HRRB(3)+&
                         2.D0*HRRB(7))+&
                         HRRB(14))+&
                         HRRB(24)
      !=(3_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(7)+&
                         ABy*(-2.D0*HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(7))+&
                         2.D0*HRRA(14))+&
                         HRRA(25)
      !=(3,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(7)+&
                         ABy*(-2.D0*HRR(3)+&
                         ABy*(ABy*HRRB(3)+&
                         3.D0*HRRB(7))+&
                         3.D0*HRRB(14))+&
                         HRRB(25)
      !=(3_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*(ABy*HRRA(9)+&
                         2.D0*HRRA(17))+&
                         HRRA(29)
      !=(3,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(14)+&
                         ABy*(2.D0*ABz*HRRB(7)+&
                         ABy*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         2.D0*HRRB(17))+&
                         HRRB(29)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(4_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*(ABy*HRRA(8)+&
                         2.D0*HRRA(16))+&
                         HRRA(28)
      !=(4,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(8)+&
                         2.D0*HRRB(16))+&
                         ABx*(ABy*(ABy*HRRB(4)+&
                         2.D0*HRRB(9))+&
                         HRRB(17))+&
                         HRRB(28)
      !=(4_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*(ABy*HRRA(9)+&
                         2.D0*HRRA(17))+&
                         HRRA(29)
      !=(4,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(9)+&
                         ABy*(-2.D0*HRR(4)+&
                         ABy*(ABy*HRRB(4)+&
                         3.D0*HRRB(9))+&
                         3.D0*HRRB(17))+&
                         HRRB(29)
      !=(4_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(7)+&
                         ABy*(-2.D0*HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(10))+&
                         2.D0*HRRA(19))+&
                         HRRA(32)
      !=(4,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(17)+&
                         ABy*(2.D0*ABz*HRRB(9)+&
                         ABy*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         2.D0*HRRB(19))+&
                         HRRB(32)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(2_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRA(11))+&
                         ABx*(-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(5))+&
                         HRRA(15))+&
                         HRRA(26)
      !=(2,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRB(11))+&
                         ABx*(2.D0*ABz*HRRB(5)+&
                         ABx*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         2.D0*HRRB(15))+&
                         HRRB(26)
      !=(2_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(12)+&
                         ABx*(ABz*HRRA(6)+&
                         HRRA(16))+&
                         HRRA(27)
      !=(2,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(12)+&
                         ABy*(ABz*HRRB(5)+&
                         HRRB(15))+&
                         ABx*(ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         HRRB(16))+&
                         HRRB(27)
      !=(2_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(15)+&
                         ABx*(ABz*HRRA(8)+&
                         HRRA(18))+&
                         HRRA(30)
      !=(2,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(5)+&
                         ABz*(ABz*HRRB(5)+&
                         2.D0*HRRB(15))+&
                         ABx*(-HRR(2)+&
                         ABz*(ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         HRRB(18))+&
                         HRRB(30)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(3_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(12)+&
                         ABx*(ABz*HRRA(6)+&
                         HRRA(16))+&
                         HRRA(27)
      !=(3,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRB(12))+&
                         ABx*(2.D0*ABz*HRRB(6)+&
                         ABx*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         2.D0*HRRB(16))+&
                         HRRB(27)
      !=(3_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRA(13))+&
                         ABx*(-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(7))+&
                         HRRA(17))+&
                         HRRA(28)
      !=(3,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(13)+&
                         ABy*(ABz*HRRB(6)+&
                         HRRB(16))+&
                         ABx*(ABz*HRRB(7)+&
                         ABy*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         HRRB(17))+&
                         HRRB(28)
      !=(3_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(16)+&
                         ABx*(ABz*HRRA(9)+&
                         HRRA(19))+&
                         HRRA(31)
      !=(3,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(6)+&
                         ABz*(ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         ABx*(-HRR(3)+&
                         ABz*(ABz*HRRB(3)+&
                         2.D0*HRRB(9))+&
                         HRRB(19))+&
                         HRRB(31)
      OffSet=(OA+3)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(4_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(15)+&
                         ABx*(ABz*HRRA(8)+&
                         HRRA(18))+&
                         HRRA(30)
      !=(4,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(10)+&
                         ABz*(-HRR(4)+&
                         HRRB(15))+&
                         ABx*(2.D0*ABz*HRRB(8)+&
                         ABx*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         2.D0*HRRB(18))+&
                         HRRB(30)
      !=(4_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(16)+&
                         ABx*(ABz*HRRA(9)+&
                         HRRA(19))+&
                         HRRA(31)
      !=(4,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(16)+&
                         ABy*(ABz*HRRB(8)+&
                         HRRB(18))+&
                         ABx*(ABz*HRRB(9)+&
                         ABy*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         HRRB(19))+&
                         HRRB(31)
      !=(4_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRA(18))+&
                         ABx*(-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(10))+&
                         HRRA(20))+&
                         HRRA(33)
      !=(4,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(8)+&
                         ABz*(ABz*HRRB(8)+&
                         2.D0*HRRB(18))+&
                         ABx*(-HRR(4)+&
                         ABz*(ABz*HRRB(4)+&
                         2.D0*HRRB(10))+&
                         HRRB(20))+&
                         HRRB(33)
      OffSet=(OA+1)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(2_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRA(12))+&
                         ABy*(-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(5))+&
                         HRRA(15))+&
                         HRRA(27)
      !=(2,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(12)+&
                         ABy*(ABz*HRRB(5)+&
                         HRRB(15))+&
                         ABx*(ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         HRRB(16))+&
                         HRRB(27)
      !=(2_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(13)+&
                         ABy*(ABz*HRRA(6)+&
                         HRRA(16))+&
                         HRRA(28)
      !=(2,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRB(13))+&
                         ABy*(2.D0*ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         2.D0*HRRB(16))+&
                         HRRB(28)
      !=(2_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(16)+&
                         ABy*(ABz*HRRA(8)+&
                         HRRA(18))+&
                         HRRA(31)
      !=(2,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(6)+&
                         ABz*(ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         ABy*(-HRR(2)+&
                         ABz*(ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         HRRB(18))+&
                         HRRB(31)
      OffSet=(OA+2)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(3_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(13)+&
                         ABy*(ABz*HRRA(6)+&
                         HRRA(16))+&
                         HRRA(28)
      !=(3,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(13)+&
                         ABy*(ABz*HRRB(6)+&
                         HRRB(16))+&
                         ABx*(ABz*HRRB(7)+&
                         ABy*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         HRRB(17))+&
                         HRRB(28)
      !=(3_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRA(14))+&
                         ABy*(-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(7))+&
                         HRRA(17))+&
                         HRRA(29)
      !=(3,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRB(14))+&
                         ABy*(2.D0*ABz*HRRB(7)+&
                         ABy*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         2.D0*HRRB(17))+&
                         HRRB(29)
      !=(3_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(17)+&
                         ABy*(ABz*HRRA(9)+&
                         HRRA(19))+&
                         HRRA(32)
      !=(3,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(7)+&
                         ABz*(ABz*HRRB(7)+&
                         2.D0*HRRB(17))+&
                         ABy*(-HRR(3)+&
                         ABz*(ABz*HRRB(3)+&
                         2.D0*HRRB(9))+&
                         HRRB(19))+&
                         HRRB(32)
      OffSet=(OA+3)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(4_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(16)+&
                         ABy*(ABz*HRRA(8)+&
                         HRRA(18))+&
                         HRRA(31)
      !=(4,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(16)+&
                         ABy*(ABz*HRRB(8)+&
                         HRRB(18))+&
                         ABx*(ABz*HRRB(9)+&
                         ABy*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         HRRB(19))+&
                         HRRB(31)
      !=(4_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(17)+&
                         ABy*(ABz*HRRA(9)+&
                         HRRA(19))+&
                         HRRA(32)
      !=(4,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(10)+&
                         ABz*(-HRR(4)+&
                         HRRB(17))+&
                         ABy*(2.D0*ABz*HRRB(9)+&
                         ABy*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         2.D0*HRRB(19))+&
                         HRRB(32)
      !=(4_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRA(19))+&
                         ABy*(-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(10))+&
                         HRRA(20))+&
                         HRRA(34)
      !=(4,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(9)+&
                         ABz*(ABz*HRRB(9)+&
                         2.D0*HRRB(19))+&
                         ABy*(-HRR(4)+&
                         ABz*(ABz*HRRB(4)+&
                         2.D0*HRRB(10))+&
                         HRRB(20))+&
                         HRRB(34)
      OffSet=(OA+1)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(2_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(5))+&
                         2.D0*HRRA(15))+&
                         HRRA(30)
      !=(2,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(5)+&
                         2.D0*HRRB(15))+&
                         ABx*(ABz*(ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         HRRB(18))+&
                         HRRB(30)
      !=(2_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*(ABz*HRRA(6)+&
                         2.D0*HRRA(16))+&
                         HRRA(31)
      !=(2,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         ABy*(ABz*(ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         HRRB(18))+&
                         HRRB(31)
      !=(2_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*(ABz*HRRA(8)+&
                         2.D0*HRRA(18))+&
                         HRRA(33)
      !=(2,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(8)+&
                         ABz*(-2.D0*HRR(2)+&
                         ABz*(ABz*HRRB(2)+&
                         3.D0*HRRB(8))+&
                         3.D0*HRRB(18))+&
                         HRRB(33)
      OffSet=(OA+2)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(3_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*(ABz*HRRA(6)+&
                         2.D0*HRRA(16))+&
                         HRRA(31)
      !=(3,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         ABx*(ABz*(ABz*HRRB(3)+&
                         2.D0*HRRB(9))+&
                         HRRB(19))+&
                         HRRB(31)
      !=(3_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(7))+&
                         2.D0*HRRA(17))+&
                         HRRA(32)
      !=(3,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(7)+&
                         2.D0*HRRB(17))+&
                         ABy*(ABz*(ABz*HRRB(3)+&
                         2.D0*HRRB(9))+&
                         HRRB(19))+&
                         HRRB(32)
      !=(3_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*(ABz*HRRA(9)+&
                         2.D0*HRRA(19))+&
                         HRRA(34)
      !=(3,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(9)+&
                         ABz*(-2.D0*HRR(3)+&
                         ABz*(ABz*HRRB(3)+&
                         3.D0*HRRB(9))+&
                         3.D0*HRRB(19))+&
                         HRRB(34)
      OffSet=(OA+3)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(4_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*(ABz*HRRA(8)+&
                         2.D0*HRRA(18))+&
                         HRRA(33)
      !=(4,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(8)+&
                         2.D0*HRRB(18))+&
                         ABx*(ABz*(ABz*HRRB(4)+&
                         2.D0*HRRB(10))+&
                         HRRB(20))+&
                         HRRB(33)
      !=(4_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*(ABz*HRRA(9)+&
                         2.D0*HRRA(19))+&
                         HRRA(34)
      !=(4,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(9)+&
                         2.D0*HRRB(19))+&
                         ABy*(ABz*(ABz*HRRB(4)+&
                         2.D0*HRRB(10))+&
                         HRRB(20))+&
                         HRRB(34)
      !=(4_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(10))+&
                         2.D0*HRRA(20))+&
                         HRRA(35)
      !=(4,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         ABz*(ABz*HRRB(4)+&
                         3.D0*HRRB(10))+&
                         3.D0*HRRB(20))+&
                         HRRB(35)
    END SUBROUTINE BraHRR26ab
    SUBROUTINE BraHRR26cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,Cart,HRR,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,Cart,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(1,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(1)+&
                                2.D0*HRR(2))+&
                                HRR(5)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(1,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(2)+&
                                ABx*(ABy*HRR(1)+&
                                HRR(3))+&
                                HRR(6)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(1,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(1)+&
                                2.D0*HRR(3))+&
                                HRR(7)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(1,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(2)+&
                                ABx*(ABz*HRR(1)+&
                                HRR(4))+&
                                HRR(8)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(1,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(3)+&
                                ABy*(ABz*HRR(1)+&
                                HRR(4))+&
                                HRR(9)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(1,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(1)+&
                                2.D0*HRR(4))+&
                                HRR(10)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(2,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(2)+&
                                2.D0*HRR(5))+&
                                HRR(11)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(3,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(3)+&
                                2.D0*HRR(6))+&
                                HRR(12)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(4,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(4)+&
                                2.D0*HRR(8))+&
                                HRR(15)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(2,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(5)+&
                                ABx*(ABy*HRR(2)+&
                                HRR(6))+&
                                HRR(12)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(3,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(6)+&
                                ABx*(ABy*HRR(3)+&
                                HRR(7))+&
                                HRR(13)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(4,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(8)+&
                                ABx*(ABy*HRR(4)+&
                                HRR(9))+&
                                HRR(16)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(2,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(2)+&
                                2.D0*HRR(6))+&
                                HRR(13)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(3,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(3)+&
                                2.D0*HRR(7))+&
                                HRR(14)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(4,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(4)+&
                                2.D0*HRR(9))+&
                                HRR(17)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(2,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(5)+&
                                ABx*(ABz*HRR(2)+&
                                HRR(8))+&
                                HRR(15)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(3,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(6)+&
                                ABx*(ABz*HRR(3)+&
                                HRR(9))+&
                                HRR(16)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(4,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(8)+&
                                ABx*(ABz*HRR(4)+&
                                HRR(10))+&
                                HRR(18)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(2,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(6)+&
                                ABy*(ABz*HRR(2)+&
                                HRR(8))+&
                                HRR(16)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(3,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(7)+&
                                ABy*(ABz*HRR(3)+&
                                HRR(9))+&
                                HRR(17)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(4,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(9)+&
                                ABy*(ABz*HRR(4)+&
                                HRR(10))+&
                                HRR(19)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(2,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(2)+&
                                2.D0*HRR(8))+&
                                HRR(18)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(3,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(3)+&
                                2.D0*HRR(9))+&
                                HRR(19)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(4,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(4)+&
                                2.D0*HRR(10))+&
                                HRR(20)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
    END SUBROUTINE BraHRR26cd
