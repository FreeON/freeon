    SUBROUTINE BraHRR110ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR,HRRA,HRRB,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*),HRRA(*),HRRB(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(1_x,11|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*(ABx*(ABx*HRRA(2)+&
                         3.D0*HRRA(5))+&
                         3.D0*HRRA(11))+&
                         HRRA(21)
      !=(1,11_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-3.D0*HRR(5)+&
                         ABx*(-6.D0*HRR(2)+&
                         ABx*(-3.D0*HRR(1)+&
                         ABx*(ABx*HRRB(1)+&
                         4.D0*HRRB(2))+&
                         6.D0*HRRB(5))+&
                         4.D0*HRRB(11))+&
                         HRRB(21)
      !=(1_y,11|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*(ABx*(ABx*HRRA(3)+&
                         3.D0*HRRA(6))+&
                         3.D0*HRRA(12))+&
                         HRRA(22)
      !=(1,11_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(11)+&
                         ABx*(3.D0*ABy*HRRB(5)+&
                         ABx*(3.D0*ABy*HRRB(2)+&
                         ABx*(ABy*HRRB(1)+&
                         HRRB(3))+&
                         3.D0*HRRB(6))+&
                         3.D0*HRRB(12))+&
                         HRRB(22)
      !=(1_z,11|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*(ABx*(ABx*HRRA(4)+&
                         3.D0*HRRA(8))+&
                         3.D0*HRRA(15))+&
                         HRRA(26)
      !=(1,11_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(11)+&
                         ABx*(3.D0*ABz*HRRB(5)+&
                         ABx*(3.D0*ABz*HRRB(2)+&
                         ABx*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         3.D0*HRRB(8))+&
                         3.D0*HRRB(15))+&
                         HRRB(26)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(1_x,12|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(11)+&
                         ABx*(2.D0*ABy*HRRA(5)+&
                         ABx*(ABy*HRRA(2)+&
                         HRRA(6))+&
                         2.D0*HRRA(12))+&
                         HRRA(22)
      !=(1,12_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(6)+&
                         ABy*(-2.D0*HRR(2)+&
                         HRRB(11))+&
                         ABx*(-2.D0*HRR(3)+&
                         ABy*(-2.D0*HRR(1)+&
                         3.D0*HRRB(5))+&
                         ABx*(3.D0*ABy*HRRB(2)+&
                         ABx*(ABy*HRRB(1)+&
                         HRRB(3))+&
                         3.D0*HRRB(6))+&
                         3.D0*HRRB(12))+&
                         HRRB(22)
      !=(1_y,12|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(12)+&
                         ABx*(2.D0*ABy*HRRA(6)+&
                         ABx*(ABy*HRRA(3)+&
                         HRRA(7))+&
                         2.D0*HRRA(13))+&
                         HRRA(23)
      !=(1,12_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(5)+&
                         ABy*(ABy*HRRB(5)+&
                         2.D0*HRRB(12))+&
                         ABx*(-2.D0*HRR(2)+&
                         ABy*(2.D0*ABy*HRRB(2)+&
                         4.D0*HRRB(6))+&
                         ABx*(-HRR(1)+&
                         ABy*(ABy*HRRB(1)+&
                         2.D0*HRRB(3))+&
                         HRRB(7))+&
                         2.D0*HRRB(13))+&
                         HRRB(23)
      !=(1_z,12|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(15)+&
                         ABx*(2.D0*ABy*HRRA(8)+&
                         ABx*(ABy*HRRA(4)+&
                         HRRA(9))+&
                         2.D0*HRRA(16))+&
                         HRRA(27)
      !=(1,12_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(12)+&
                         ABy*(ABz*HRRB(5)+&
                         HRRB(15))+&
                         ABx*(2.D0*ABz*HRRB(6)+&
                         ABy*(2.D0*ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         ABx*(ABz*HRRB(3)+&
                         ABy*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         HRRB(9))+&
                         2.D0*HRRB(16))+&
                         HRRB(27)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(1_x,13|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*(ABy*HRRA(5)+&
                         2.D0*HRRA(12))+&
                         ABx*(ABy*(ABy*HRRA(2)+&
                         2.D0*HRRA(6))+&
                         HRRA(13))+&
                         HRRA(23)
      !=(1,13_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(7)+&
                         ABy*(-2.D0*HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRB(5))+&
                         2.D0*HRRB(12))+&
                         ABx*(ABy*(2.D0*ABy*HRRB(2)+&
                         4.D0*HRRB(6))+&
                         ABx*(ABy*(ABy*HRRB(1)+&
                         2.D0*HRRB(3))+&
                         HRRB(7))+&
                         2.D0*HRRB(13))+&
                         HRRB(23)
      !=(1_y,13|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*(ABy*HRRA(6)+&
                         2.D0*HRRA(13))+&
                         ABx*(ABy*(ABy*HRRA(3)+&
                         2.D0*HRRA(7))+&
                         HRRA(14))+&
                         HRRA(24)
      !=(1,13_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(6)+&
                         ABy*(-2.D0*HRR(2)+&
                         ABy*(ABy*HRRB(2)+&
                         3.D0*HRRB(6))+&
                         3.D0*HRRB(13))+&
                         ABx*(-2.D0*HRR(3)+&
                         ABy*(-2.D0*HRR(1)+&
                         ABy*(ABy*HRRB(1)+&
                         3.D0*HRRB(3))+&
                         3.D0*HRRB(7))+&
                         HRRB(14))+&
                         HRRB(24)
      !=(1_z,13|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*(ABy*HRRA(8)+&
                         2.D0*HRRA(16))+&
                         ABx*(ABy*(ABy*HRRA(4)+&
                         2.D0*HRRA(9))+&
                         HRRA(17))+&
                         HRRA(28)
      !=(1,13_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(13)+&
                         ABy*(2.D0*ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         2.D0*HRRB(16))+&
                         ABx*(ABz*HRRB(7)+&
                         ABy*(2.D0*ABz*HRRB(3)+&
                         ABy*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         2.D0*HRRB(9))+&
                         HRRB(17))+&
                         HRRB(28)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(1_x,14|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*(ABy*(ABy*HRRA(2)+&
                         3.D0*HRRA(6))+&
                         3.D0*HRRA(13))+&
                         HRRA(24)
      !=(1,14_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*(ABy*HRRB(2)+&
                         3.D0*HRRB(6))+&
                         3.D0*HRRB(13))+&
                         ABx*(ABy*(ABy*(ABy*HRRB(1)+&
                         3.D0*HRRB(3))+&
                         3.D0*HRRB(7))+&
                         HRRB(14))+&
                         HRRB(24)
      !=(1_y,14|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*(ABy*(ABy*HRRA(3)+&
                         3.D0*HRRA(7))+&
                         3.D0*HRRA(14))+&
                         HRRA(25)
      !=(1,14_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-3.D0*HRR(7)+&
                         ABy*(-6.D0*HRR(3)+&
                         ABy*(-3.D0*HRR(1)+&
                         ABy*(ABy*HRRB(1)+&
                         4.D0*HRRB(3))+&
                         6.D0*HRRB(7))+&
                         4.D0*HRRB(14))+&
                         HRRB(25)
      !=(1_z,14|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*(ABy*(ABy*HRRA(4)+&
                         3.D0*HRRA(9))+&
                         3.D0*HRRA(17))+&
                         HRRA(29)
      !=(1,14_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(14)+&
                         ABy*(3.D0*ABz*HRRB(7)+&
                         ABy*(3.D0*ABz*HRRB(3)+&
                         ABy*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         3.D0*HRRB(9))+&
                         3.D0*HRRB(17))+&
                         HRRB(29)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(1_x,15|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(11)+&
                         ABx*(2.D0*ABz*HRRA(5)+&
                         ABx*(ABz*HRRA(2)+&
                         HRRA(8))+&
                         2.D0*HRRA(15))+&
                         HRRA(26)
      !=(1,15_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(8)+&
                         ABz*(-2.D0*HRR(2)+&
                         HRRB(11))+&
                         ABx*(-2.D0*HRR(4)+&
                         ABz*(-2.D0*HRR(1)+&
                         3.D0*HRRB(5))+&
                         ABx*(3.D0*ABz*HRRB(2)+&
                         ABx*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         3.D0*HRRB(8))+&
                         3.D0*HRRB(15))+&
                         HRRB(26)
      !=(1_y,15|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(12)+&
                         ABx*(2.D0*ABz*HRRA(6)+&
                         ABx*(ABz*HRRA(3)+&
                         HRRA(9))+&
                         2.D0*HRRA(16))+&
                         HRRA(27)
      !=(1,15_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(12)+&
                         ABy*(ABz*HRRB(5)+&
                         HRRB(15))+&
                         ABx*(2.D0*ABz*HRRB(6)+&
                         ABy*(2.D0*ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         ABx*(ABz*HRRB(3)+&
                         ABy*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         HRRB(9))+&
                         2.D0*HRRB(16))+&
                         HRRB(27)
      !=(1_z,15|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(15)+&
                         ABx*(2.D0*ABz*HRRA(8)+&
                         ABx*(ABz*HRRA(4)+&
                         HRRA(10))+&
                         2.D0*HRRA(18))+&
                         HRRA(30)
      !=(1,15_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(5)+&
                         ABz*(ABz*HRRB(5)+&
                         2.D0*HRRB(15))+&
                         ABx*(-2.D0*HRR(2)+&
                         ABz*(2.D0*ABz*HRRB(2)+&
                         4.D0*HRRB(8))+&
                         ABx*(-HRR(1)+&
                         ABz*(ABz*HRRB(1)+&
                         2.D0*HRRB(4))+&
                         HRRB(10))+&
                         2.D0*HRRB(18))+&
                         HRRB(30)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(1_x,16|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(12)+&
                         ABy*(ABz*HRRA(5)+&
                         HRRA(15))+&
                         ABx*(ABz*HRRA(6)+&
                         ABy*(ABz*HRRA(2)+&
                         HRRA(8))+&
                         HRRA(16))+&
                         HRRA(27)
      !=(1,16_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRB(12))+&
                         ABy*(-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRB(5))+&
                         HRRB(15))+&
                         ABx*(2.D0*ABz*HRRB(6)+&
                         ABy*(2.D0*ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         ABx*(ABz*HRRB(3)+&
                         ABy*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         HRRB(9))+&
                         2.D0*HRRB(16))+&
                         HRRB(27)
      !=(1_y,16|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(13)+&
                         ABy*(ABz*HRRA(6)+&
                         HRRA(16))+&
                         ABx*(ABz*HRRA(7)+&
                         ABy*(ABz*HRRA(3)+&
                         HRRA(9))+&
                         HRRA(17))+&
                         HRRA(28)
      !=(1,16_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRB(13))+&
                         ABy*(2.D0*ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         2.D0*HRRB(16))+&
                         ABx*(-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRB(7))+&
                         ABy*(2.D0*ABz*HRRB(3)+&
                         ABy*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         2.D0*HRRB(9))+&
                         HRRB(17))+&
                         HRRB(28)
      !=(1_z,16|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(16)+&
                         ABy*(ABz*HRRA(8)+&
                         HRRA(18))+&
                         ABx*(ABz*HRRA(9)+&
                         ABy*(ABz*HRRA(4)+&
                         HRRA(10))+&
                         HRRA(19))+&
                         HRRA(31)
      !=(1,16_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(6)+&
                         ABz*(ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         ABy*(-HRR(2)+&
                         ABz*(ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         HRRB(18))+&
                         ABx*(-HRR(3)+&
                         ABz*(ABz*HRRB(3)+&
                         2.D0*HRRB(9))+&
                         ABy*(-HRR(1)+&
                         ABz*(ABz*HRRB(1)+&
                         2.D0*HRRB(4))+&
                         HRRB(10))+&
                         HRRB(19))+&
                         HRRB(31)
      OffSet=(OA+0)*LDA+(OB+6)*LDB+CDOffSet !=
      !=(1_x,17|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(13)+&
                         ABy*(2.D0*ABz*HRRA(6)+&
                         ABy*(ABz*HRRA(2)+&
                         HRRA(8))+&
                         2.D0*HRRA(16))+&
                         HRRA(28)
      !=(1,17_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(13)+&
                         ABy*(2.D0*ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         2.D0*HRRB(16))+&
                         ABx*(ABz*HRRB(7)+&
                         ABy*(2.D0*ABz*HRRB(3)+&
                         ABy*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         2.D0*HRRB(9))+&
                         HRRB(17))+&
                         HRRB(28)
      !=(1_y,17|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(14)+&
                         ABy*(2.D0*ABz*HRRA(7)+&
                         ABy*(ABz*HRRA(3)+&
                         HRRA(9))+&
                         2.D0*HRRA(17))+&
                         HRRA(29)
      !=(1,17_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(9)+&
                         ABz*(-2.D0*HRR(3)+&
                         HRRB(14))+&
                         ABy*(-2.D0*HRR(4)+&
                         ABz*(-2.D0*HRR(1)+&
                         3.D0*HRRB(7))+&
                         ABy*(3.D0*ABz*HRRB(3)+&
                         ABy*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         3.D0*HRRB(9))+&
                         3.D0*HRRB(17))+&
                         HRRB(29)
      !=(1_z,17|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(17)+&
                         ABy*(2.D0*ABz*HRRA(9)+&
                         ABy*(ABz*HRRA(4)+&
                         HRRA(10))+&
                         2.D0*HRRA(19))+&
                         HRRA(32)
      !=(1,17_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(7)+&
                         ABz*(ABz*HRRB(7)+&
                         2.D0*HRRB(17))+&
                         ABy*(-2.D0*HRR(3)+&
                         ABz*(2.D0*ABz*HRRB(3)+&
                         4.D0*HRRB(9))+&
                         ABy*(-HRR(1)+&
                         ABz*(ABz*HRRB(1)+&
                         2.D0*HRRB(4))+&
                         HRRB(10))+&
                         2.D0*HRRB(19))+&
                         HRRB(32)
      OffSet=(OA+0)*LDA+(OB+7)*LDB+CDOffSet !=
      !=(1_x,18|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*(ABz*HRRA(5)+&
                         2.D0*HRRA(15))+&
                         ABx*(ABz*(ABz*HRRA(2)+&
                         2.D0*HRRA(8))+&
                         HRRA(18))+&
                         HRRA(30)
      !=(1,18_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRB(5))+&
                         2.D0*HRRB(15))+&
                         ABx*(ABz*(2.D0*ABz*HRRB(2)+&
                         4.D0*HRRB(8))+&
                         ABx*(ABz*(ABz*HRRB(1)+&
                         2.D0*HRRB(4))+&
                         HRRB(10))+&
                         2.D0*HRRB(18))+&
                         HRRB(30)
      !=(1_y,18|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*(ABz*HRRA(6)+&
                         2.D0*HRRA(16))+&
                         ABx*(ABz*(ABz*HRRA(3)+&
                         2.D0*HRRA(9))+&
                         HRRA(19))+&
                         HRRA(31)
      !=(1,18_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         ABy*(ABz*(ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         HRRB(18))+&
                         ABx*(ABz*(ABz*HRRB(3)+&
                         2.D0*HRRB(9))+&
                         ABy*(ABz*(ABz*HRRB(1)+&
                         2.D0*HRRB(4))+&
                         HRRB(10))+&
                         HRRB(19))+&
                         HRRB(31)
      !=(1_z,18|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*(ABz*HRRA(8)+&
                         2.D0*HRRA(18))+&
                         ABx*(ABz*(ABz*HRRA(4)+&
                         2.D0*HRRA(10))+&
                         HRRA(20))+&
                         HRRA(33)
      !=(1,18_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(8)+&
                         ABz*(-2.D0*HRR(2)+&
                         ABz*(ABz*HRRB(2)+&
                         3.D0*HRRB(8))+&
                         3.D0*HRRB(18))+&
                         ABx*(-2.D0*HRR(4)+&
                         ABz*(-2.D0*HRR(1)+&
                         ABz*(ABz*HRRB(1)+&
                         3.D0*HRRB(4))+&
                         3.D0*HRRB(10))+&
                         HRRB(20))+&
                         HRRB(33)
      OffSet=(OA+0)*LDA+(OB+8)*LDB+CDOffSet !=
      !=(1_x,19|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*(ABz*HRRA(6)+&
                         2.D0*HRRA(16))+&
                         ABy*(ABz*(ABz*HRRA(2)+&
                         2.D0*HRRA(8))+&
                         HRRA(18))+&
                         HRRA(31)
      !=(1,19_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         ABy*(ABz*(ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         HRRB(18))+&
                         ABx*(ABz*(ABz*HRRB(3)+&
                         2.D0*HRRB(9))+&
                         ABy*(ABz*(ABz*HRRB(1)+&
                         2.D0*HRRB(4))+&
                         HRRB(10))+&
                         HRRB(19))+&
                         HRRB(31)
      !=(1_y,19|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*(ABz*HRRA(7)+&
                         2.D0*HRRA(17))+&
                         ABy*(ABz*(ABz*HRRA(3)+&
                         2.D0*HRRA(9))+&
                         HRRA(19))+&
                         HRRA(32)
      !=(1,19_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRB(7))+&
                         2.D0*HRRB(17))+&
                         ABy*(ABz*(2.D0*ABz*HRRB(3)+&
                         4.D0*HRRB(9))+&
                         ABy*(ABz*(ABz*HRRB(1)+&
                         2.D0*HRRB(4))+&
                         HRRB(10))+&
                         2.D0*HRRB(19))+&
                         HRRB(32)
      !=(1_z,19|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*(ABz*HRRA(9)+&
                         2.D0*HRRA(19))+&
                         ABy*(ABz*(ABz*HRRA(4)+&
                         2.D0*HRRA(10))+&
                         HRRA(20))+&
                         HRRA(34)
      !=(1,19_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(9)+&
                         ABz*(-2.D0*HRR(3)+&
                         ABz*(ABz*HRRB(3)+&
                         3.D0*HRRB(9))+&
                         3.D0*HRRB(19))+&
                         ABy*(-2.D0*HRR(4)+&
                         ABz*(-2.D0*HRR(1)+&
                         ABz*(ABz*HRRB(1)+&
                         3.D0*HRRB(4))+&
                         3.D0*HRRB(10))+&
                         HRRB(20))+&
                         HRRB(34)
      OffSet=(OA+0)*LDA+(OB+9)*LDB+CDOffSet !=
      !=(1_x,20|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*(ABz*(ABz*HRRA(2)+&
                         3.D0*HRRA(8))+&
                         3.D0*HRRA(18))+&
                         HRRA(33)
      !=(1,20_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*(ABz*HRRB(2)+&
                         3.D0*HRRB(8))+&
                         3.D0*HRRB(18))+&
                         ABx*(ABz*(ABz*(ABz*HRRB(1)+&
                         3.D0*HRRB(4))+&
                         3.D0*HRRB(10))+&
                         HRRB(20))+&
                         HRRB(33)
      !=(1_y,20|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*(ABz*(ABz*HRRA(3)+&
                         3.D0*HRRA(9))+&
                         3.D0*HRRA(19))+&
                         HRRA(34)
      !=(1,20_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*(ABz*HRRB(3)+&
                         3.D0*HRRB(9))+&
                         3.D0*HRRB(19))+&
                         ABy*(ABz*(ABz*(ABz*HRRB(1)+&
                         3.D0*HRRB(4))+&
                         3.D0*HRRB(10))+&
                         HRRB(20))+&
                         HRRB(34)
      !=(1_z,20|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*(ABz*(ABz*HRRA(4)+&
                         3.D0*HRRA(10))+&
                         3.D0*HRRA(20))+&
                         HRRA(35)
      !=(1,20_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-3.D0*HRR(10)+&
                         ABz*(-6.D0*HRR(4)+&
                         ABz*(-3.D0*HRR(1)+&
                         ABz*(ABz*HRRB(1)+&
                         4.D0*HRRB(4))+&
                         6.D0*HRRB(10))+&
                         4.D0*HRRB(20))+&
                         HRRB(35)
    END SUBROUTINE BraHRR110ab
    SUBROUTINE BraHRR110cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,Cart,HRR,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,Cart,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(1,11|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*(ABx*HRR(1)+&
                                3.D0*HRR(2))+&
                                3.D0*HRR(5))+&
                                HRR(11)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(1,12|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(5)+&
                                ABx*(2.D0*ABy*HRR(2)+&
                                ABx*(ABy*HRR(1)+&
                                HRR(3))+&
                                2.D0*HRR(6))+&
                                HRR(12)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(1,13|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(2)+&
                                2.D0*HRR(6))+&
                                ABx*(ABy*(ABy*HRR(1)+&
                                2.D0*HRR(3))+&
                                HRR(7))+&
                                HRR(13)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(1,14|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*(ABy*HRR(1)+&
                                3.D0*HRR(3))+&
                                3.D0*HRR(7))+&
                                HRR(14)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(1,15|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(5)+&
                                ABx*(2.D0*ABz*HRR(2)+&
                                ABx*(ABz*HRR(1)+&
                                HRR(4))+&
                                2.D0*HRR(8))+&
                                HRR(15)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(1,16|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(6)+&
                                ABy*(ABz*HRR(2)+&
                                HRR(8))+&
                                ABx*(ABz*HRR(3)+&
                                ABy*(ABz*HRR(1)+&
                                HRR(4))+&
                                HRR(9))+&
                                HRR(16)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+6)*LDB+CDOffSet !=
      !=(1,17|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(7)+&
                                ABy*(2.D0*ABz*HRR(3)+&
                                ABy*(ABz*HRR(1)+&
                                HRR(4))+&
                                2.D0*HRR(9))+&
                                HRR(17)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+7)*LDB+CDOffSet !=
      !=(1,18|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(2)+&
                                2.D0*HRR(8))+&
                                ABx*(ABz*(ABz*HRR(1)+&
                                2.D0*HRR(4))+&
                                HRR(10))+&
                                HRR(18)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+8)*LDB+CDOffSet !=
      !=(1,19|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(3)+&
                                2.D0*HRR(9))+&
                                ABy*(ABz*(ABz*HRR(1)+&
                                2.D0*HRR(4))+&
                                HRR(10))+&
                                HRR(19)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+9)*LDB+CDOffSet !=
      !=(1,20|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*(ABz*HRR(1)+&
                                3.D0*HRR(4))+&
                                3.D0*HRR(10))+&
                                HRR(20)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
    END SUBROUTINE BraHRR110cd
