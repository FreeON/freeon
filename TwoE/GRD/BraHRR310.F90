    SUBROUTINE BraHRR310ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR,HRRA,HRRB,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*),HRRA(*),HRRB(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(2_x,11|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(11)+&
                         ABx*(-3.D0*HRR(5)+&
                         ABx*(-3.D0*HRR(2)+&
                         ABx*(-HRR(1)+&
                         HRRA(5))+&
                         3.D0*HRRA(11))+&
                         3.D0*HRRA(21))+&
                         HRRA(36)
      !=(2,11_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-3.D0*HRR(11)+&
                         ABx*(-6.D0*HRR(5)+&
                         ABx*(-3.D0*HRR(2)+&
                         ABx*(ABx*HRRB(2)+&
                         4.D0*HRRB(5))+&
                         6.D0*HRRB(11))+&
                         4.D0*HRRB(21))+&
                         HRRB(36)
      !=(2_y,11|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*(ABx*(ABx*HRRA(6)+&
                         3.D0*HRRA(12))+&
                         3.D0*HRRA(22))+&
                         HRRA(37)
      !=(2,11_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(21)+&
                         ABx*(3.D0*ABy*HRRB(11)+&
                         ABx*(3.D0*ABy*HRRB(5)+&
                         ABx*(ABy*HRRB(2)+&
                         HRRB(6))+&
                         3.D0*HRRB(12))+&
                         3.D0*HRRB(22))+&
                         HRRB(37)
      !=(2_z,11|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*(ABx*(ABx*HRRA(8)+&
                         3.D0*HRRA(15))+&
                         3.D0*HRRA(26))+&
                         HRRA(42)
      !=(2,11_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(21)+&
                         ABx*(3.D0*ABz*HRRB(11)+&
                         ABx*(3.D0*ABz*HRRB(5)+&
                         ABx*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         3.D0*HRRB(15))+&
                         3.D0*HRRB(26))+&
                         HRRB(42)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(3_x,11|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*(ABx*(ABx*HRRA(6)+&
                         3.D0*HRRA(12))+&
                         3.D0*HRRA(22))+&
                         HRRA(37)
      !=(3,11_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-3.D0*HRR(12)+&
                         ABx*(-6.D0*HRR(6)+&
                         ABx*(-3.D0*HRR(3)+&
                         ABx*(ABx*HRRB(3)+&
                         4.D0*HRRB(6))+&
                         6.D0*HRRB(12))+&
                         4.D0*HRRB(22))+&
                         HRRB(37)
      !=(3_y,11|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(11)+&
                         ABx*(-3.D0*HRR(5)+&
                         ABx*(-3.D0*HRR(2)+&
                         ABx*(-HRR(1)+&
                         HRRA(7))+&
                         3.D0*HRRA(13))+&
                         3.D0*HRRA(23))+&
                         HRRA(38)
      !=(3,11_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(22)+&
                         ABx*(3.D0*ABy*HRRB(12)+&
                         ABx*(3.D0*ABy*HRRB(6)+&
                         ABx*(ABy*HRRB(3)+&
                         HRRB(7))+&
                         3.D0*HRRB(13))+&
                         3.D0*HRRB(23))+&
                         HRRB(38)
      !=(3_z,11|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*(ABx*(ABx*HRRA(9)+&
                         3.D0*HRRA(16))+&
                         3.D0*HRRA(27))+&
                         HRRA(43)
      !=(3,11_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(22)+&
                         ABx*(3.D0*ABz*HRRB(12)+&
                         ABx*(3.D0*ABz*HRRB(6)+&
                         ABx*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         3.D0*HRRB(16))+&
                         3.D0*HRRB(27))+&
                         HRRB(43)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(4_x,11|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*(ABx*(ABx*HRRA(8)+&
                         3.D0*HRRA(15))+&
                         3.D0*HRRA(26))+&
                         HRRA(42)
      !=(4,11_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-3.D0*HRR(15)+&
                         ABx*(-6.D0*HRR(8)+&
                         ABx*(-3.D0*HRR(4)+&
                         ABx*(ABx*HRRB(4)+&
                         4.D0*HRRB(8))+&
                         6.D0*HRRB(15))+&
                         4.D0*HRRB(26))+&
                         HRRB(42)
      !=(4_y,11|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*(ABx*(ABx*HRRA(9)+&
                         3.D0*HRRA(16))+&
                         3.D0*HRRA(27))+&
                         HRRA(43)
      !=(4,11_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(26)+&
                         ABx*(3.D0*ABy*HRRB(15)+&
                         ABx*(3.D0*ABy*HRRB(8)+&
                         ABx*(ABy*HRRB(4)+&
                         HRRB(9))+&
                         3.D0*HRRB(16))+&
                         3.D0*HRRB(27))+&
                         HRRB(43)
      !=(4_z,11|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(11)+&
                         ABx*(-3.D0*HRR(5)+&
                         ABx*(-3.D0*HRR(2)+&
                         ABx*(-HRR(1)+&
                         HRRA(10))+&
                         3.D0*HRRA(18))+&
                         3.D0*HRRA(30))+&
                         HRRA(47)
      !=(4,11_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(26)+&
                         ABx*(3.D0*ABz*HRRB(15)+&
                         ABx*(3.D0*ABz*HRRB(8)+&
                         ABx*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         3.D0*HRRB(18))+&
                         3.D0*HRRB(30))+&
                         HRRB(47)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(2_x,12|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(12)+&
                         ABy*(-HRR(5)+&
                         HRRA(21))+&
                         ABx*(-2.D0*HRR(6)+&
                         ABy*(-2.D0*HRR(2)+&
                         2.D0*HRRA(11))+&
                         ABx*(-HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(5))+&
                         HRRA(12))+&
                         2.D0*HRRA(22))+&
                         HRRA(37)
      !=(2,12_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(12)+&
                         ABy*(-2.D0*HRR(5)+&
                         HRRB(21))+&
                         ABx*(-2.D0*HRR(6)+&
                         ABy*(-2.D0*HRR(2)+&
                         3.D0*HRRB(11))+&
                         ABx*(3.D0*ABy*HRRB(5)+&
                         ABx*(ABy*HRRB(2)+&
                         HRRB(6))+&
                         3.D0*HRRB(12))+&
                         3.D0*HRRB(22))+&
                         HRRB(37)
      !=(2_y,12|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(22)+&
                         ABx*(2.D0*ABy*HRRA(12)+&
                         ABx*(ABy*HRRA(6)+&
                         HRRA(13))+&
                         2.D0*HRRA(23))+&
                         HRRA(38)
      !=(2,12_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(11)+&
                         ABy*(ABy*HRRB(11)+&
                         2.D0*HRRB(22))+&
                         ABx*(-2.D0*HRR(5)+&
                         ABy*(2.D0*ABy*HRRB(5)+&
                         4.D0*HRRB(12))+&
                         ABx*(-HRR(2)+&
                         ABy*(ABy*HRRB(2)+&
                         2.D0*HRRB(6))+&
                         HRRB(13))+&
                         2.D0*HRRB(23))+&
                         HRRB(38)
      !=(2_z,12|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(26)+&
                         ABx*(2.D0*ABy*HRRA(15)+&
                         ABx*(ABy*HRRA(8)+&
                         HRRA(16))+&
                         2.D0*HRRA(27))+&
                         HRRA(43)
      !=(2,12_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(22)+&
                         ABy*(ABz*HRRB(11)+&
                         HRRB(26))+&
                         ABx*(2.D0*ABz*HRRB(12)+&
                         ABy*(2.D0*ABz*HRRB(5)+&
                         2.D0*HRRB(15))+&
                         ABx*(ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         HRRB(16))+&
                         2.D0*HRRB(27))+&
                         HRRB(43)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(3_x,12|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(22)+&
                         ABx*(2.D0*ABy*HRRA(12)+&
                         ABx*(ABy*HRRA(6)+&
                         HRRA(13))+&
                         2.D0*HRRA(23))+&
                         HRRA(38)
      !=(3,12_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(13)+&
                         ABy*(-2.D0*HRR(6)+&
                         HRRB(22))+&
                         ABx*(-2.D0*HRR(7)+&
                         ABy*(-2.D0*HRR(3)+&
                         3.D0*HRRB(12))+&
                         ABx*(3.D0*ABy*HRRB(6)+&
                         ABx*(ABy*HRRB(3)+&
                         HRRB(7))+&
                         3.D0*HRRB(13))+&
                         3.D0*HRRB(23))+&
                         HRRB(38)
      !=(3_y,12|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(12)+&
                         ABy*(-HRR(5)+&
                         HRRA(23))+&
                         ABx*(-2.D0*HRR(6)+&
                         ABy*(-2.D0*HRR(2)+&
                         2.D0*HRRA(13))+&
                         ABx*(-HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(7))+&
                         HRRA(14))+&
                         2.D0*HRRA(24))+&
                         HRRA(39)
      !=(3,12_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(12)+&
                         ABy*(ABy*HRRB(12)+&
                         2.D0*HRRB(23))+&
                         ABx*(-2.D0*HRR(6)+&
                         ABy*(2.D0*ABy*HRRB(6)+&
                         4.D0*HRRB(13))+&
                         ABx*(-HRR(3)+&
                         ABy*(ABy*HRRB(3)+&
                         2.D0*HRRB(7))+&
                         HRRB(14))+&
                         2.D0*HRRB(24))+&
                         HRRB(39)
      !=(3_z,12|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(27)+&
                         ABx*(2.D0*ABy*HRRA(16)+&
                         ABx*(ABy*HRRA(9)+&
                         HRRA(17))+&
                         2.D0*HRRA(28))+&
                         HRRA(44)
      !=(3,12_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(23)+&
                         ABy*(ABz*HRRB(12)+&
                         HRRB(27))+&
                         ABx*(2.D0*ABz*HRRB(13)+&
                         ABy*(2.D0*ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         ABx*(ABz*HRRB(7)+&
                         ABy*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         HRRB(17))+&
                         2.D0*HRRB(28))+&
                         HRRB(44)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(4_x,12|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(26)+&
                         ABx*(2.D0*ABy*HRRA(15)+&
                         ABx*(ABy*HRRA(8)+&
                         HRRA(16))+&
                         2.D0*HRRA(27))+&
                         HRRA(43)
      !=(4,12_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(16)+&
                         ABy*(-2.D0*HRR(8)+&
                         HRRB(26))+&
                         ABx*(-2.D0*HRR(9)+&
                         ABy*(-2.D0*HRR(4)+&
                         3.D0*HRRB(15))+&
                         ABx*(3.D0*ABy*HRRB(8)+&
                         ABx*(ABy*HRRB(4)+&
                         HRRB(9))+&
                         3.D0*HRRB(16))+&
                         3.D0*HRRB(27))+&
                         HRRB(43)
      !=(4_y,12|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(27)+&
                         ABx*(2.D0*ABy*HRRA(16)+&
                         ABx*(ABy*HRRA(9)+&
                         HRRA(17))+&
                         2.D0*HRRA(28))+&
                         HRRA(44)
      !=(4,12_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(15)+&
                         ABy*(ABy*HRRB(15)+&
                         2.D0*HRRB(27))+&
                         ABx*(-2.D0*HRR(8)+&
                         ABy*(2.D0*ABy*HRRB(8)+&
                         4.D0*HRRB(16))+&
                         ABx*(-HRR(4)+&
                         ABy*(ABy*HRRB(4)+&
                         2.D0*HRRB(9))+&
                         HRRB(17))+&
                         2.D0*HRRB(28))+&
                         HRRB(44)
      !=(4_z,12|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(12)+&
                         ABy*(-HRR(5)+&
                         HRRA(30))+&
                         ABx*(-2.D0*HRR(6)+&
                         ABy*(-2.D0*HRR(2)+&
                         2.D0*HRRA(18))+&
                         ABx*(-HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(10))+&
                         HRRA(19))+&
                         2.D0*HRRA(31))+&
                         HRRA(48)
      !=(4,12_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(27)+&
                         ABy*(ABz*HRRB(15)+&
                         HRRB(30))+&
                         ABx*(2.D0*ABz*HRRB(16)+&
                         ABy*(2.D0*ABz*HRRB(8)+&
                         2.D0*HRRB(18))+&
                         ABx*(ABz*HRRB(9)+&
                         ABy*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         HRRB(19))+&
                         2.D0*HRRB(31))+&
                         HRRB(48)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(2_x,13|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(13)+&
                         ABy*(-2.D0*HRR(6)+&
                         ABy*(-HRR(2)+&
                         HRRA(11))+&
                         2.D0*HRRA(22))+&
                         ABx*(-HRR(7)+&
                         ABy*(-2.D0*HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(5))+&
                         2.D0*HRRA(12))+&
                         HRRA(23))+&
                         HRRA(38)
      !=(2,13_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(13)+&
                         ABy*(-2.D0*HRR(6)+&
                         ABy*(-HRR(2)+&
                         HRRB(11))+&
                         2.D0*HRRB(22))+&
                         ABx*(ABy*(2.D0*ABy*HRRB(5)+&
                         4.D0*HRRB(12))+&
                         ABx*(ABy*(ABy*HRRB(2)+&
                         2.D0*HRRB(6))+&
                         HRRB(13))+&
                         2.D0*HRRB(23))+&
                         HRRB(38)
      !=(2_y,13|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*(ABy*HRRA(12)+&
                         2.D0*HRRA(23))+&
                         ABx*(ABy*(ABy*HRRA(6)+&
                         2.D0*HRRA(13))+&
                         HRRA(24))+&
                         HRRA(39)
      !=(2,13_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(12)+&
                         ABy*(-2.D0*HRR(5)+&
                         ABy*(ABy*HRRB(5)+&
                         3.D0*HRRB(12))+&
                         3.D0*HRRB(23))+&
                         ABx*(-2.D0*HRR(6)+&
                         ABy*(-2.D0*HRR(2)+&
                         ABy*(ABy*HRRB(2)+&
                         3.D0*HRRB(6))+&
                         3.D0*HRRB(13))+&
                         HRRB(24))+&
                         HRRB(39)
      !=(2_z,13|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*(ABy*HRRA(15)+&
                         2.D0*HRRA(27))+&
                         ABx*(ABy*(ABy*HRRA(8)+&
                         2.D0*HRRA(16))+&
                         HRRA(28))+&
                         HRRA(44)
      !=(2,13_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(23)+&
                         ABy*(2.D0*ABz*HRRB(12)+&
                         ABy*(ABz*HRRB(5)+&
                         HRRB(15))+&
                         2.D0*HRRB(27))+&
                         ABx*(ABz*HRRB(13)+&
                         ABy*(2.D0*ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         2.D0*HRRB(16))+&
                         HRRB(28))+&
                         HRRB(44)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(3_x,13|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*(ABy*HRRA(12)+&
                         2.D0*HRRA(23))+&
                         ABx*(ABy*(ABy*HRRA(6)+&
                         2.D0*HRRA(13))+&
                         HRRA(24))+&
                         HRRA(39)
      !=(3,13_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(14)+&
                         ABy*(-2.D0*HRR(7)+&
                         ABy*(-HRR(3)+&
                         HRRB(12))+&
                         2.D0*HRRB(23))+&
                         ABx*(ABy*(2.D0*ABy*HRRB(6)+&
                         4.D0*HRRB(13))+&
                         ABx*(ABy*(ABy*HRRB(3)+&
                         2.D0*HRRB(7))+&
                         HRRB(14))+&
                         2.D0*HRRB(24))+&
                         HRRB(39)
      !=(3_y,13|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(13)+&
                         ABy*(-2.D0*HRR(6)+&
                         ABy*(-HRR(2)+&
                         HRRA(13))+&
                         2.D0*HRRA(24))+&
                         ABx*(-HRR(7)+&
                         ABy*(-2.D0*HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(7))+&
                         2.D0*HRRA(14))+&
                         HRRA(25))+&
                         HRRA(40)
      !=(3,13_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(13)+&
                         ABy*(-2.D0*HRR(6)+&
                         ABy*(ABy*HRRB(6)+&
                         3.D0*HRRB(13))+&
                         3.D0*HRRB(24))+&
                         ABx*(-2.D0*HRR(7)+&
                         ABy*(-2.D0*HRR(3)+&
                         ABy*(ABy*HRRB(3)+&
                         3.D0*HRRB(7))+&
                         3.D0*HRRB(14))+&
                         HRRB(25))+&
                         HRRB(40)
      !=(3_z,13|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*(ABy*HRRA(16)+&
                         2.D0*HRRA(28))+&
                         ABx*(ABy*(ABy*HRRA(9)+&
                         2.D0*HRRA(17))+&
                         HRRA(29))+&
                         HRRA(45)
      !=(3,13_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(24)+&
                         ABy*(2.D0*ABz*HRRB(13)+&
                         ABy*(ABz*HRRB(6)+&
                         HRRB(16))+&
                         2.D0*HRRB(28))+&
                         ABx*(ABz*HRRB(14)+&
                         ABy*(2.D0*ABz*HRRB(7)+&
                         ABy*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         2.D0*HRRB(17))+&
                         HRRB(29))+&
                         HRRB(45)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(4_x,13|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*(ABy*HRRA(15)+&
                         2.D0*HRRA(27))+&
                         ABx*(ABy*(ABy*HRRA(8)+&
                         2.D0*HRRA(16))+&
                         HRRA(28))+&
                         HRRA(44)
      !=(4,13_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(17)+&
                         ABy*(-2.D0*HRR(9)+&
                         ABy*(-HRR(4)+&
                         HRRB(15))+&
                         2.D0*HRRB(27))+&
                         ABx*(ABy*(2.D0*ABy*HRRB(8)+&
                         4.D0*HRRB(16))+&
                         ABx*(ABy*(ABy*HRRB(4)+&
                         2.D0*HRRB(9))+&
                         HRRB(17))+&
                         2.D0*HRRB(28))+&
                         HRRB(44)
      !=(4_y,13|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*(ABy*HRRA(16)+&
                         2.D0*HRRA(28))+&
                         ABx*(ABy*(ABy*HRRA(9)+&
                         2.D0*HRRA(17))+&
                         HRRA(29))+&
                         HRRA(45)
      !=(4,13_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(16)+&
                         ABy*(-2.D0*HRR(8)+&
                         ABy*(ABy*HRRB(8)+&
                         3.D0*HRRB(16))+&
                         3.D0*HRRB(28))+&
                         ABx*(-2.D0*HRR(9)+&
                         ABy*(-2.D0*HRR(4)+&
                         ABy*(ABy*HRRB(4)+&
                         3.D0*HRRB(9))+&
                         3.D0*HRRB(17))+&
                         HRRB(29))+&
                         HRRB(45)
      !=(4_z,13|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(13)+&
                         ABy*(-2.D0*HRR(6)+&
                         ABy*(-HRR(2)+&
                         HRRA(18))+&
                         2.D0*HRRA(31))+&
                         ABx*(-HRR(7)+&
                         ABy*(-2.D0*HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(10))+&
                         2.D0*HRRA(19))+&
                         HRRA(32))+&
                         HRRA(49)
      !=(4,13_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(28)+&
                         ABy*(2.D0*ABz*HRRB(16)+&
                         ABy*(ABz*HRRB(8)+&
                         HRRB(18))+&
                         2.D0*HRRB(31))+&
                         ABx*(ABz*HRRB(17)+&
                         ABy*(2.D0*ABz*HRRB(9)+&
                         ABy*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         2.D0*HRRB(19))+&
                         HRRB(32))+&
                         HRRB(49)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(2_x,14|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(14)+&
                         ABy*(-3.D0*HRR(7)+&
                         ABy*(-3.D0*HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(5))+&
                         3.D0*HRRA(12))+&
                         3.D0*HRRA(23))+&
                         HRRA(39)
      !=(2,14_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*(ABy*HRRB(5)+&
                         3.D0*HRRB(12))+&
                         3.D0*HRRB(23))+&
                         ABx*(ABy*(ABy*(ABy*HRRB(2)+&
                         3.D0*HRRB(6))+&
                         3.D0*HRRB(13))+&
                         HRRB(24))+&
                         HRRB(39)
      !=(2_y,14|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*(ABy*(ABy*HRRA(6)+&
                         3.D0*HRRA(13))+&
                         3.D0*HRRA(24))+&
                         HRRA(40)
      !=(2,14_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-3.D0*HRR(13)+&
                         ABy*(-6.D0*HRR(6)+&
                         ABy*(-3.D0*HRR(2)+&
                         ABy*(ABy*HRRB(2)+&
                         4.D0*HRRB(6))+&
                         6.D0*HRRB(13))+&
                         4.D0*HRRB(24))+&
                         HRRB(40)
      !=(2_z,14|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*(ABy*(ABy*HRRA(8)+&
                         3.D0*HRRA(16))+&
                         3.D0*HRRA(28))+&
                         HRRA(45)
      !=(2,14_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(24)+&
                         ABy*(3.D0*ABz*HRRB(13)+&
                         ABy*(3.D0*ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         3.D0*HRRB(16))+&
                         3.D0*HRRB(28))+&
                         HRRB(45)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(3_x,14|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*(ABy*(ABy*HRRA(6)+&
                         3.D0*HRRA(13))+&
                         3.D0*HRRA(24))+&
                         HRRA(40)
      !=(3,14_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*(ABy*HRRB(6)+&
                         3.D0*HRRB(13))+&
                         3.D0*HRRB(24))+&
                         ABx*(ABy*(ABy*(ABy*HRRB(3)+&
                         3.D0*HRRB(7))+&
                         3.D0*HRRB(14))+&
                         HRRB(25))+&
                         HRRB(40)
      !=(3_y,14|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(14)+&
                         ABy*(-3.D0*HRR(7)+&
                         ABy*(-3.D0*HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(7))+&
                         3.D0*HRRA(14))+&
                         3.D0*HRRA(25))+&
                         HRRA(41)
      !=(3,14_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-3.D0*HRR(14)+&
                         ABy*(-6.D0*HRR(7)+&
                         ABy*(-3.D0*HRR(3)+&
                         ABy*(ABy*HRRB(3)+&
                         4.D0*HRRB(7))+&
                         6.D0*HRRB(14))+&
                         4.D0*HRRB(25))+&
                         HRRB(41)
      !=(3_z,14|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*(ABy*(ABy*HRRA(9)+&
                         3.D0*HRRA(17))+&
                         3.D0*HRRA(29))+&
                         HRRA(46)
      !=(3,14_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(25)+&
                         ABy*(3.D0*ABz*HRRB(14)+&
                         ABy*(3.D0*ABz*HRRB(7)+&
                         ABy*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         3.D0*HRRB(17))+&
                         3.D0*HRRB(29))+&
                         HRRB(46)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(4_x,14|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*(ABy*(ABy*HRRA(8)+&
                         3.D0*HRRA(16))+&
                         3.D0*HRRA(28))+&
                         HRRA(45)
      !=(4,14_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*(ABy*HRRB(8)+&
                         3.D0*HRRB(16))+&
                         3.D0*HRRB(28))+&
                         ABx*(ABy*(ABy*(ABy*HRRB(4)+&
                         3.D0*HRRB(9))+&
                         3.D0*HRRB(17))+&
                         HRRB(29))+&
                         HRRB(45)
      !=(4_y,14|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*(ABy*(ABy*HRRA(9)+&
                         3.D0*HRRA(17))+&
                         3.D0*HRRA(29))+&
                         HRRA(46)
      !=(4,14_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-3.D0*HRR(17)+&
                         ABy*(-6.D0*HRR(9)+&
                         ABy*(-3.D0*HRR(4)+&
                         ABy*(ABy*HRRB(4)+&
                         4.D0*HRRB(9))+&
                         6.D0*HRRB(17))+&
                         4.D0*HRRB(29))+&
                         HRRB(46)
      !=(4_z,14|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(14)+&
                         ABy*(-3.D0*HRR(7)+&
                         ABy*(-3.D0*HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(10))+&
                         3.D0*HRRA(19))+&
                         3.D0*HRRA(32))+&
                         HRRA(50)
      !=(4,14_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(29)+&
                         ABy*(3.D0*ABz*HRRB(17)+&
                         ABy*(3.D0*ABz*HRRB(9)+&
                         ABy*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         3.D0*HRRB(19))+&
                         3.D0*HRRB(32))+&
                         HRRB(50)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(2_x,15|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(15)+&
                         ABz*(-HRR(5)+&
                         HRRA(21))+&
                         ABx*(-2.D0*HRR(8)+&
                         ABz*(-2.D0*HRR(2)+&
                         2.D0*HRRA(11))+&
                         ABx*(-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(5))+&
                         HRRA(15))+&
                         2.D0*HRRA(26))+&
                         HRRA(42)
      !=(2,15_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(15)+&
                         ABz*(-2.D0*HRR(5)+&
                         HRRB(21))+&
                         ABx*(-2.D0*HRR(8)+&
                         ABz*(-2.D0*HRR(2)+&
                         3.D0*HRRB(11))+&
                         ABx*(3.D0*ABz*HRRB(5)+&
                         ABx*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         3.D0*HRRB(15))+&
                         3.D0*HRRB(26))+&
                         HRRB(42)
      !=(2_y,15|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(22)+&
                         ABx*(2.D0*ABz*HRRA(12)+&
                         ABx*(ABz*HRRA(6)+&
                         HRRA(16))+&
                         2.D0*HRRA(27))+&
                         HRRA(43)
      !=(2,15_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(22)+&
                         ABy*(ABz*HRRB(11)+&
                         HRRB(26))+&
                         ABx*(2.D0*ABz*HRRB(12)+&
                         ABy*(2.D0*ABz*HRRB(5)+&
                         2.D0*HRRB(15))+&
                         ABx*(ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         HRRB(16))+&
                         2.D0*HRRB(27))+&
                         HRRB(43)
      !=(2_z,15|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(26)+&
                         ABx*(2.D0*ABz*HRRA(15)+&
                         ABx*(ABz*HRRA(8)+&
                         HRRA(18))+&
                         2.D0*HRRA(30))+&
                         HRRA(47)
      !=(2,15_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(11)+&
                         ABz*(ABz*HRRB(11)+&
                         2.D0*HRRB(26))+&
                         ABx*(-2.D0*HRR(5)+&
                         ABz*(2.D0*ABz*HRRB(5)+&
                         4.D0*HRRB(15))+&
                         ABx*(-HRR(2)+&
                         ABz*(ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         HRRB(18))+&
                         2.D0*HRRB(30))+&
                         HRRB(47)
      OffSet=(OA+1)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(3_x,15|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(22)+&
                         ABx*(2.D0*ABz*HRRA(12)+&
                         ABx*(ABz*HRRA(6)+&
                         HRRA(16))+&
                         2.D0*HRRA(27))+&
                         HRRA(43)
      !=(3,15_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(16)+&
                         ABz*(-2.D0*HRR(6)+&
                         HRRB(22))+&
                         ABx*(-2.D0*HRR(9)+&
                         ABz*(-2.D0*HRR(3)+&
                         3.D0*HRRB(12))+&
                         ABx*(3.D0*ABz*HRRB(6)+&
                         ABx*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         3.D0*HRRB(16))+&
                         3.D0*HRRB(27))+&
                         HRRB(43)
      !=(3_y,15|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(15)+&
                         ABz*(-HRR(5)+&
                         HRRA(23))+&
                         ABx*(-2.D0*HRR(8)+&
                         ABz*(-2.D0*HRR(2)+&
                         2.D0*HRRA(13))+&
                         ABx*(-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(7))+&
                         HRRA(17))+&
                         2.D0*HRRA(28))+&
                         HRRA(44)
      !=(3,15_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(23)+&
                         ABy*(ABz*HRRB(12)+&
                         HRRB(27))+&
                         ABx*(2.D0*ABz*HRRB(13)+&
                         ABy*(2.D0*ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         ABx*(ABz*HRRB(7)+&
                         ABy*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         HRRB(17))+&
                         2.D0*HRRB(28))+&
                         HRRB(44)
      !=(3_z,15|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(27)+&
                         ABx*(2.D0*ABz*HRRA(16)+&
                         ABx*(ABz*HRRA(9)+&
                         HRRA(19))+&
                         2.D0*HRRA(31))+&
                         HRRA(48)
      !=(3,15_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(12)+&
                         ABz*(ABz*HRRB(12)+&
                         2.D0*HRRB(27))+&
                         ABx*(-2.D0*HRR(6)+&
                         ABz*(2.D0*ABz*HRRB(6)+&
                         4.D0*HRRB(16))+&
                         ABx*(-HRR(3)+&
                         ABz*(ABz*HRRB(3)+&
                         2.D0*HRRB(9))+&
                         HRRB(19))+&
                         2.D0*HRRB(31))+&
                         HRRB(48)
      OffSet=(OA+2)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(4_x,15|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(26)+&
                         ABx*(2.D0*ABz*HRRA(15)+&
                         ABx*(ABz*HRRA(8)+&
                         HRRA(18))+&
                         2.D0*HRRA(30))+&
                         HRRA(47)
      !=(4,15_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(18)+&
                         ABz*(-2.D0*HRR(8)+&
                         HRRB(26))+&
                         ABx*(-2.D0*HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         3.D0*HRRB(15))+&
                         ABx*(3.D0*ABz*HRRB(8)+&
                         ABx*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         3.D0*HRRB(18))+&
                         3.D0*HRRB(30))+&
                         HRRB(47)
      !=(4_y,15|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(27)+&
                         ABx*(2.D0*ABz*HRRA(16)+&
                         ABx*(ABz*HRRA(9)+&
                         HRRA(19))+&
                         2.D0*HRRA(31))+&
                         HRRA(48)
      !=(4,15_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(27)+&
                         ABy*(ABz*HRRB(15)+&
                         HRRB(30))+&
                         ABx*(2.D0*ABz*HRRB(16)+&
                         ABy*(2.D0*ABz*HRRB(8)+&
                         2.D0*HRRB(18))+&
                         ABx*(ABz*HRRB(9)+&
                         ABy*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         HRRB(19))+&
                         2.D0*HRRB(31))+&
                         HRRB(48)
      !=(4_z,15|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(15)+&
                         ABz*(-HRR(5)+&
                         HRRA(30))+&
                         ABx*(-2.D0*HRR(8)+&
                         ABz*(-2.D0*HRR(2)+&
                         2.D0*HRRA(18))+&
                         ABx*(-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(10))+&
                         HRRA(20))+&
                         2.D0*HRRA(33))+&
                         HRRA(51)
      !=(4,15_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(15)+&
                         ABz*(ABz*HRRB(15)+&
                         2.D0*HRRB(30))+&
                         ABx*(-2.D0*HRR(8)+&
                         ABz*(2.D0*ABz*HRRB(8)+&
                         4.D0*HRRB(18))+&
                         ABx*(-HRR(4)+&
                         ABz*(ABz*HRRB(4)+&
                         2.D0*HRRB(10))+&
                         HRRB(20))+&
                         2.D0*HRRB(33))+&
                         HRRB(51)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(2_x,16|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(16)+&
                         ABz*(-HRR(6)+&
                         HRRA(22))+&
                         ABy*(-HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRA(11))+&
                         HRRA(26))+&
                         ABx*(-HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRA(12))+&
                         ABy*(-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(5))+&
                         HRRA(15))+&
                         HRRA(27))+&
                         HRRA(43)
      !=(2,16_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(16)+&
                         ABz*(-HRR(6)+&
                         HRRB(22))+&
                         ABy*(-HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRB(11))+&
                         HRRB(26))+&
                         ABx*(2.D0*ABz*HRRB(12)+&
                         ABy*(2.D0*ABz*HRRB(5)+&
                         2.D0*HRRB(15))+&
                         ABx*(ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         HRRB(16))+&
                         2.D0*HRRB(27))+&
                         HRRB(43)
      !=(2_y,16|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(23)+&
                         ABy*(ABz*HRRA(12)+&
                         HRRA(27))+&
                         ABx*(ABz*HRRA(13)+&
                         ABy*(ABz*HRRA(6)+&
                         HRRA(16))+&
                         HRRA(28))+&
                         HRRA(44)
      !=(2,16_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(15)+&
                         ABz*(-HRR(5)+&
                         HRRB(23))+&
                         ABy*(2.D0*ABz*HRRB(12)+&
                         ABy*(ABz*HRRB(5)+&
                         HRRB(15))+&
                         2.D0*HRRB(27))+&
                         ABx*(-HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRB(13))+&
                         ABy*(2.D0*ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         2.D0*HRRB(16))+&
                         HRRB(28))+&
                         HRRB(44)
      !=(2_z,16|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(27)+&
                         ABy*(ABz*HRRA(15)+&
                         HRRA(30))+&
                         ABx*(ABz*HRRA(16)+&
                         ABy*(ABz*HRRA(8)+&
                         HRRA(18))+&
                         HRRA(31))+&
                         HRRA(48)
      !=(2,16_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(12)+&
                         ABz*(ABz*HRRB(12)+&
                         2.D0*HRRB(27))+&
                         ABy*(-HRR(5)+&
                         ABz*(ABz*HRRB(5)+&
                         2.D0*HRRB(15))+&
                         HRRB(30))+&
                         ABx*(-HRR(6)+&
                         ABz*(ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         ABy*(-HRR(2)+&
                         ABz*(ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         HRRB(18))+&
                         HRRB(31))+&
                         HRRB(48)
      OffSet=(OA+1)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(3_x,16|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(23)+&
                         ABy*(ABz*HRRA(12)+&
                         HRRA(27))+&
                         ABx*(ABz*HRRA(13)+&
                         ABy*(ABz*HRRA(6)+&
                         HRRA(16))+&
                         HRRA(28))+&
                         HRRA(44)
      !=(3,16_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(17)+&
                         ABz*(-HRR(7)+&
                         HRRB(23))+&
                         ABy*(-HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRB(12))+&
                         HRRB(27))+&
                         ABx*(2.D0*ABz*HRRB(13)+&
                         ABy*(2.D0*ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         ABx*(ABz*HRRB(7)+&
                         ABy*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         HRRB(17))+&
                         2.D0*HRRB(28))+&
                         HRRB(44)
      !=(3_y,16|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(16)+&
                         ABz*(-HRR(6)+&
                         HRRA(24))+&
                         ABy*(-HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRA(13))+&
                         HRRA(28))+&
                         ABx*(-HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRA(14))+&
                         ABy*(-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(7))+&
                         HRRA(17))+&
                         HRRA(29))+&
                         HRRA(45)
      !=(3,16_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(16)+&
                         ABz*(-HRR(6)+&
                         HRRB(24))+&
                         ABy*(2.D0*ABz*HRRB(13)+&
                         ABy*(ABz*HRRB(6)+&
                         HRRB(16))+&
                         2.D0*HRRB(28))+&
                         ABx*(-HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRB(14))+&
                         ABy*(2.D0*ABz*HRRB(7)+&
                         ABy*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         2.D0*HRRB(17))+&
                         HRRB(29))+&
                         HRRB(45)
      !=(3_z,16|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(28)+&
                         ABy*(ABz*HRRA(16)+&
                         HRRA(31))+&
                         ABx*(ABz*HRRA(17)+&
                         ABy*(ABz*HRRA(9)+&
                         HRRA(19))+&
                         HRRA(32))+&
                         HRRA(49)
      !=(3,16_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(13)+&
                         ABz*(ABz*HRRB(13)+&
                         2.D0*HRRB(28))+&
                         ABy*(-HRR(6)+&
                         ABz*(ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         HRRB(31))+&
                         ABx*(-HRR(7)+&
                         ABz*(ABz*HRRB(7)+&
                         2.D0*HRRB(17))+&
                         ABy*(-HRR(3)+&
                         ABz*(ABz*HRRB(3)+&
                         2.D0*HRRB(9))+&
                         HRRB(19))+&
                         HRRB(32))+&
                         HRRB(49)
      OffSet=(OA+2)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(4_x,16|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(27)+&
                         ABy*(ABz*HRRA(15)+&
                         HRRA(30))+&
                         ABx*(ABz*HRRA(16)+&
                         ABy*(ABz*HRRA(8)+&
                         HRRA(18))+&
                         HRRA(31))+&
                         HRRA(48)
      !=(4,16_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(19)+&
                         ABz*(-HRR(9)+&
                         HRRB(27))+&
                         ABy*(-HRR(10)+&
                         ABz*(-HRR(4)+&
                         HRRB(15))+&
                         HRRB(30))+&
                         ABx*(2.D0*ABz*HRRB(16)+&
                         ABy*(2.D0*ABz*HRRB(8)+&
                         2.D0*HRRB(18))+&
                         ABx*(ABz*HRRB(9)+&
                         ABy*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         HRRB(19))+&
                         2.D0*HRRB(31))+&
                         HRRB(48)
      !=(4_y,16|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(28)+&
                         ABy*(ABz*HRRA(16)+&
                         HRRA(31))+&
                         ABx*(ABz*HRRA(17)+&
                         ABy*(ABz*HRRA(9)+&
                         HRRA(19))+&
                         HRRA(32))+&
                         HRRA(49)
      !=(4,16_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(18)+&
                         ABz*(-HRR(8)+&
                         HRRB(28))+&
                         ABy*(2.D0*ABz*HRRB(16)+&
                         ABy*(ABz*HRRB(8)+&
                         HRRB(18))+&
                         2.D0*HRRB(31))+&
                         ABx*(-HRR(10)+&
                         ABz*(-HRR(4)+&
                         HRRB(17))+&
                         ABy*(2.D0*ABz*HRRB(9)+&
                         ABy*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         2.D0*HRRB(19))+&
                         HRRB(32))+&
                         HRRB(49)
      !=(4_z,16|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(16)+&
                         ABz*(-HRR(6)+&
                         HRRA(31))+&
                         ABy*(-HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRA(18))+&
                         HRRA(33))+&
                         ABx*(-HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRA(19))+&
                         ABy*(-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(10))+&
                         HRRA(20))+&
                         HRRA(34))+&
                         HRRA(52)
      !=(4,16_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(16)+&
                         ABz*(ABz*HRRB(16)+&
                         2.D0*HRRB(31))+&
                         ABy*(-HRR(8)+&
                         ABz*(ABz*HRRB(8)+&
                         2.D0*HRRB(18))+&
                         HRRB(33))+&
                         ABx*(-HRR(9)+&
                         ABz*(ABz*HRRB(9)+&
                         2.D0*HRRB(19))+&
                         ABy*(-HRR(4)+&
                         ABz*(ABz*HRRB(4)+&
                         2.D0*HRRB(10))+&
                         HRRB(20))+&
                         HRRB(34))+&
                         HRRB(52)
      OffSet=(OA+0)*LDA+(OB+6)*LDB+CDOffSet !=
      !=(2_x,17|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(17)+&
                         ABz*(-HRR(7)+&
                         HRRA(23))+&
                         ABy*(-2.D0*HRR(9)+&
                         ABz*(-2.D0*HRR(3)+&
                         2.D0*HRRA(12))+&
                         ABy*(-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(5))+&
                         HRRA(15))+&
                         2.D0*HRRA(27))+&
                         HRRA(44)
      !=(2,17_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(23)+&
                         ABy*(2.D0*ABz*HRRB(12)+&
                         ABy*(ABz*HRRB(5)+&
                         HRRB(15))+&
                         2.D0*HRRB(27))+&
                         ABx*(ABz*HRRB(13)+&
                         ABy*(2.D0*ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         2.D0*HRRB(16))+&
                         HRRB(28))+&
                         HRRB(44)
      !=(2_y,17|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(24)+&
                         ABy*(2.D0*ABz*HRRA(13)+&
                         ABy*(ABz*HRRA(6)+&
                         HRRA(16))+&
                         2.D0*HRRA(28))+&
                         HRRA(45)
      !=(2,17_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(16)+&
                         ABz*(-2.D0*HRR(6)+&
                         HRRB(24))+&
                         ABy*(-2.D0*HRR(8)+&
                         ABz*(-2.D0*HRR(2)+&
                         3.D0*HRRB(13))+&
                         ABy*(3.D0*ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         3.D0*HRRB(16))+&
                         3.D0*HRRB(28))+&
                         HRRB(45)
      !=(2_z,17|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(28)+&
                         ABy*(2.D0*ABz*HRRA(16)+&
                         ABy*(ABz*HRRA(8)+&
                         HRRA(18))+&
                         2.D0*HRRA(31))+&
                         HRRA(49)
      !=(2,17_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(13)+&
                         ABz*(ABz*HRRB(13)+&
                         2.D0*HRRB(28))+&
                         ABy*(-2.D0*HRR(6)+&
                         ABz*(2.D0*ABz*HRRB(6)+&
                         4.D0*HRRB(16))+&
                         ABy*(-HRR(2)+&
                         ABz*(ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         HRRB(18))+&
                         2.D0*HRRB(31))+&
                         HRRB(49)
      OffSet=(OA+1)*LDA+(OB+6)*LDB+CDOffSet !=
      !=(3_x,17|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(24)+&
                         ABy*(2.D0*ABz*HRRA(13)+&
                         ABy*(ABz*HRRA(6)+&
                         HRRA(16))+&
                         2.D0*HRRA(28))+&
                         HRRA(45)
      !=(3,17_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(24)+&
                         ABy*(2.D0*ABz*HRRB(13)+&
                         ABy*(ABz*HRRB(6)+&
                         HRRB(16))+&
                         2.D0*HRRB(28))+&
                         ABx*(ABz*HRRB(14)+&
                         ABy*(2.D0*ABz*HRRB(7)+&
                         ABy*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         2.D0*HRRB(17))+&
                         HRRB(29))+&
                         HRRB(45)
      !=(3_y,17|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(17)+&
                         ABz*(-HRR(7)+&
                         HRRA(25))+&
                         ABy*(-2.D0*HRR(9)+&
                         ABz*(-2.D0*HRR(3)+&
                         2.D0*HRRA(14))+&
                         ABy*(-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(7))+&
                         HRRA(17))+&
                         2.D0*HRRA(29))+&
                         HRRA(46)
      !=(3,17_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(17)+&
                         ABz*(-2.D0*HRR(7)+&
                         HRRB(25))+&
                         ABy*(-2.D0*HRR(9)+&
                         ABz*(-2.D0*HRR(3)+&
                         3.D0*HRRB(14))+&
                         ABy*(3.D0*ABz*HRRB(7)+&
                         ABy*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         3.D0*HRRB(17))+&
                         3.D0*HRRB(29))+&
                         HRRB(46)
      !=(3_z,17|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(29)+&
                         ABy*(2.D0*ABz*HRRA(17)+&
                         ABy*(ABz*HRRA(9)+&
                         HRRA(19))+&
                         2.D0*HRRA(32))+&
                         HRRA(50)
      !=(3,17_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(14)+&
                         ABz*(ABz*HRRB(14)+&
                         2.D0*HRRB(29))+&
                         ABy*(-2.D0*HRR(7)+&
                         ABz*(2.D0*ABz*HRRB(7)+&
                         4.D0*HRRB(17))+&
                         ABy*(-HRR(3)+&
                         ABz*(ABz*HRRB(3)+&
                         2.D0*HRRB(9))+&
                         HRRB(19))+&
                         2.D0*HRRB(32))+&
                         HRRB(50)
      OffSet=(OA+2)*LDA+(OB+6)*LDB+CDOffSet !=
      !=(4_x,17|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(28)+&
                         ABy*(2.D0*ABz*HRRA(16)+&
                         ABy*(ABz*HRRA(8)+&
                         HRRA(18))+&
                         2.D0*HRRA(31))+&
                         HRRA(49)
      !=(4,17_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(28)+&
                         ABy*(2.D0*ABz*HRRB(16)+&
                         ABy*(ABz*HRRB(8)+&
                         HRRB(18))+&
                         2.D0*HRRB(31))+&
                         ABx*(ABz*HRRB(17)+&
                         ABy*(2.D0*ABz*HRRB(9)+&
                         ABy*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         2.D0*HRRB(19))+&
                         HRRB(32))+&
                         HRRB(49)
      !=(4_y,17|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(29)+&
                         ABy*(2.D0*ABz*HRRA(17)+&
                         ABy*(ABz*HRRA(9)+&
                         HRRA(19))+&
                         2.D0*HRRA(32))+&
                         HRRA(50)
      !=(4,17_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(19)+&
                         ABz*(-2.D0*HRR(9)+&
                         HRRB(29))+&
                         ABy*(-2.D0*HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         3.D0*HRRB(17))+&
                         ABy*(3.D0*ABz*HRRB(9)+&
                         ABy*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         3.D0*HRRB(19))+&
                         3.D0*HRRB(32))+&
                         HRRB(50)
      !=(4_z,17|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(17)+&
                         ABz*(-HRR(7)+&
                         HRRA(32))+&
                         ABy*(-2.D0*HRR(9)+&
                         ABz*(-2.D0*HRR(3)+&
                         2.D0*HRRA(19))+&
                         ABy*(-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(10))+&
                         HRRA(20))+&
                         2.D0*HRRA(34))+&
                         HRRA(53)
      !=(4,17_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(17)+&
                         ABz*(ABz*HRRB(17)+&
                         2.D0*HRRB(32))+&
                         ABy*(-2.D0*HRR(9)+&
                         ABz*(2.D0*ABz*HRRB(9)+&
                         4.D0*HRRB(19))+&
                         ABy*(-HRR(4)+&
                         ABz*(ABz*HRRB(4)+&
                         2.D0*HRRB(10))+&
                         HRRB(20))+&
                         2.D0*HRRB(34))+&
                         HRRB(53)
      OffSet=(OA+0)*LDA+(OB+7)*LDB+CDOffSet !=
      !=(2_x,18|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(18)+&
                         ABz*(-2.D0*HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRA(11))+&
                         2.D0*HRRA(26))+&
                         ABx*(-HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(5))+&
                         2.D0*HRRA(15))+&
                         HRRA(30))+&
                         HRRA(47)
      !=(2,18_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(18)+&
                         ABz*(-2.D0*HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRB(11))+&
                         2.D0*HRRB(26))+&
                         ABx*(ABz*(2.D0*ABz*HRRB(5)+&
                         4.D0*HRRB(15))+&
                         ABx*(ABz*(ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         HRRB(18))+&
                         2.D0*HRRB(30))+&
                         HRRB(47)
      !=(2_y,18|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*(ABz*HRRA(12)+&
                         2.D0*HRRA(27))+&
                         ABx*(ABz*(ABz*HRRA(6)+&
                         2.D0*HRRA(16))+&
                         HRRA(31))+&
                         HRRA(48)
      !=(2,18_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(12)+&
                         2.D0*HRRB(27))+&
                         ABy*(ABz*(ABz*HRRB(5)+&
                         2.D0*HRRB(15))+&
                         HRRB(30))+&
                         ABx*(ABz*(ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         ABy*(ABz*(ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         HRRB(18))+&
                         HRRB(31))+&
                         HRRB(48)
      !=(2_z,18|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*(ABz*HRRA(15)+&
                         2.D0*HRRA(30))+&
                         ABx*(ABz*(ABz*HRRA(8)+&
                         2.D0*HRRA(18))+&
                         HRRA(33))+&
                         HRRA(51)
      !=(2,18_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(15)+&
                         ABz*(-2.D0*HRR(5)+&
                         ABz*(ABz*HRRB(5)+&
                         3.D0*HRRB(15))+&
                         3.D0*HRRB(30))+&
                         ABx*(-2.D0*HRR(8)+&
                         ABz*(-2.D0*HRR(2)+&
                         ABz*(ABz*HRRB(2)+&
                         3.D0*HRRB(8))+&
                         3.D0*HRRB(18))+&
                         HRRB(33))+&
                         HRRB(51)
      OffSet=(OA+1)*LDA+(OB+7)*LDB+CDOffSet !=
      !=(3_x,18|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*(ABz*HRRA(12)+&
                         2.D0*HRRA(27))+&
                         ABx*(ABz*(ABz*HRRA(6)+&
                         2.D0*HRRA(16))+&
                         HRRA(31))+&
                         HRRA(48)
      !=(3,18_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(19)+&
                         ABz*(-2.D0*HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRB(12))+&
                         2.D0*HRRB(27))+&
                         ABx*(ABz*(2.D0*ABz*HRRB(6)+&
                         4.D0*HRRB(16))+&
                         ABx*(ABz*(ABz*HRRB(3)+&
                         2.D0*HRRB(9))+&
                         HRRB(19))+&
                         2.D0*HRRB(31))+&
                         HRRB(48)
      !=(3_y,18|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(18)+&
                         ABz*(-2.D0*HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRA(13))+&
                         2.D0*HRRA(28))+&
                         ABx*(-HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(7))+&
                         2.D0*HRRA(17))+&
                         HRRA(32))+&
                         HRRA(49)
      !=(3,18_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(13)+&
                         2.D0*HRRB(28))+&
                         ABy*(ABz*(ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         HRRB(31))+&
                         ABx*(ABz*(ABz*HRRB(7)+&
                         2.D0*HRRB(17))+&
                         ABy*(ABz*(ABz*HRRB(3)+&
                         2.D0*HRRB(9))+&
                         HRRB(19))+&
                         HRRB(32))+&
                         HRRB(49)
      !=(3_z,18|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*(ABz*HRRA(16)+&
                         2.D0*HRRA(31))+&
                         ABx*(ABz*(ABz*HRRA(9)+&
                         2.D0*HRRA(19))+&
                         HRRA(34))+&
                         HRRA(52)
      !=(3,18_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(16)+&
                         ABz*(-2.D0*HRR(6)+&
                         ABz*(ABz*HRRB(6)+&
                         3.D0*HRRB(16))+&
                         3.D0*HRRB(31))+&
                         ABx*(-2.D0*HRR(9)+&
                         ABz*(-2.D0*HRR(3)+&
                         ABz*(ABz*HRRB(3)+&
                         3.D0*HRRB(9))+&
                         3.D0*HRRB(19))+&
                         HRRB(34))+&
                         HRRB(52)
      OffSet=(OA+2)*LDA+(OB+7)*LDB+CDOffSet !=
      !=(4_x,18|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*(ABz*HRRA(15)+&
                         2.D0*HRRA(30))+&
                         ABx*(ABz*(ABz*HRRA(8)+&
                         2.D0*HRRA(18))+&
                         HRRA(33))+&
                         HRRA(51)
      !=(4,18_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(20)+&
                         ABz*(-2.D0*HRR(10)+&
                         ABz*(-HRR(4)+&
                         HRRB(15))+&
                         2.D0*HRRB(30))+&
                         ABx*(ABz*(2.D0*ABz*HRRB(8)+&
                         4.D0*HRRB(18))+&
                         ABx*(ABz*(ABz*HRRB(4)+&
                         2.D0*HRRB(10))+&
                         HRRB(20))+&
                         2.D0*HRRB(33))+&
                         HRRB(51)
      !=(4_y,18|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*(ABz*HRRA(16)+&
                         2.D0*HRRA(31))+&
                         ABx*(ABz*(ABz*HRRA(9)+&
                         2.D0*HRRA(19))+&
                         HRRA(34))+&
                         HRRA(52)
      !=(4,18_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(16)+&
                         2.D0*HRRB(31))+&
                         ABy*(ABz*(ABz*HRRB(8)+&
                         2.D0*HRRB(18))+&
                         HRRB(33))+&
                         ABx*(ABz*(ABz*HRRB(9)+&
                         2.D0*HRRB(19))+&
                         ABy*(ABz*(ABz*HRRB(4)+&
                         2.D0*HRRB(10))+&
                         HRRB(20))+&
                         HRRB(34))+&
                         HRRB(52)
      !=(4_z,18|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(18)+&
                         ABz*(-2.D0*HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRA(18))+&
                         2.D0*HRRA(33))+&
                         ABx*(-HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(10))+&
                         2.D0*HRRA(20))+&
                         HRRA(35))+&
                         HRRA(54)
      !=(4,18_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(18)+&
                         ABz*(-2.D0*HRR(8)+&
                         ABz*(ABz*HRRB(8)+&
                         3.D0*HRRB(18))+&
                         3.D0*HRRB(33))+&
                         ABx*(-2.D0*HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         ABz*(ABz*HRRB(4)+&
                         3.D0*HRRB(10))+&
                         3.D0*HRRB(20))+&
                         HRRB(35))+&
                         HRRB(54)
      OffSet=(OA+0)*LDA+(OB+8)*LDB+CDOffSet !=
      !=(2_x,19|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(19)+&
                         ABz*(-2.D0*HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRA(12))+&
                         2.D0*HRRA(27))+&
                         ABy*(-HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(5))+&
                         2.D0*HRRA(15))+&
                         HRRA(30))+&
                         HRRA(48)
      !=(2,19_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(12)+&
                         2.D0*HRRB(27))+&
                         ABy*(ABz*(ABz*HRRB(5)+&
                         2.D0*HRRB(15))+&
                         HRRB(30))+&
                         ABx*(ABz*(ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         ABy*(ABz*(ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         HRRB(18))+&
                         HRRB(31))+&
                         HRRB(48)
      !=(2_y,19|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*(ABz*HRRA(13)+&
                         2.D0*HRRA(28))+&
                         ABy*(ABz*(ABz*HRRA(6)+&
                         2.D0*HRRA(16))+&
                         HRRA(31))+&
                         HRRA(49)
      !=(2,19_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(18)+&
                         ABz*(-2.D0*HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRB(13))+&
                         2.D0*HRRB(28))+&
                         ABy*(ABz*(2.D0*ABz*HRRB(6)+&
                         4.D0*HRRB(16))+&
                         ABy*(ABz*(ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         HRRB(18))+&
                         2.D0*HRRB(31))+&
                         HRRB(49)
      !=(2_z,19|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*(ABz*HRRA(16)+&
                         2.D0*HRRA(31))+&
                         ABy*(ABz*(ABz*HRRA(8)+&
                         2.D0*HRRA(18))+&
                         HRRA(33))+&
                         HRRA(52)
      !=(2,19_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(16)+&
                         ABz*(-2.D0*HRR(6)+&
                         ABz*(ABz*HRRB(6)+&
                         3.D0*HRRB(16))+&
                         3.D0*HRRB(31))+&
                         ABy*(-2.D0*HRR(8)+&
                         ABz*(-2.D0*HRR(2)+&
                         ABz*(ABz*HRRB(2)+&
                         3.D0*HRRB(8))+&
                         3.D0*HRRB(18))+&
                         HRRB(33))+&
                         HRRB(52)
      OffSet=(OA+1)*LDA+(OB+8)*LDB+CDOffSet !=
      !=(3_x,19|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*(ABz*HRRA(13)+&
                         2.D0*HRRA(28))+&
                         ABy*(ABz*(ABz*HRRA(6)+&
                         2.D0*HRRA(16))+&
                         HRRA(31))+&
                         HRRA(49)
      !=(3,19_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(13)+&
                         2.D0*HRRB(28))+&
                         ABy*(ABz*(ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         HRRB(31))+&
                         ABx*(ABz*(ABz*HRRB(7)+&
                         2.D0*HRRB(17))+&
                         ABy*(ABz*(ABz*HRRB(3)+&
                         2.D0*HRRB(9))+&
                         HRRB(19))+&
                         HRRB(32))+&
                         HRRB(49)
      !=(3_y,19|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(19)+&
                         ABz*(-2.D0*HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRA(14))+&
                         2.D0*HRRA(29))+&
                         ABy*(-HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(7))+&
                         2.D0*HRRA(17))+&
                         HRRA(32))+&
                         HRRA(50)
      !=(3,19_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(19)+&
                         ABz*(-2.D0*HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRB(14))+&
                         2.D0*HRRB(29))+&
                         ABy*(ABz*(2.D0*ABz*HRRB(7)+&
                         4.D0*HRRB(17))+&
                         ABy*(ABz*(ABz*HRRB(3)+&
                         2.D0*HRRB(9))+&
                         HRRB(19))+&
                         2.D0*HRRB(32))+&
                         HRRB(50)
      !=(3_z,19|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*(ABz*HRRA(17)+&
                         2.D0*HRRA(32))+&
                         ABy*(ABz*(ABz*HRRA(9)+&
                         2.D0*HRRA(19))+&
                         HRRA(34))+&
                         HRRA(53)
      !=(3,19_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(17)+&
                         ABz*(-2.D0*HRR(7)+&
                         ABz*(ABz*HRRB(7)+&
                         3.D0*HRRB(17))+&
                         3.D0*HRRB(32))+&
                         ABy*(-2.D0*HRR(9)+&
                         ABz*(-2.D0*HRR(3)+&
                         ABz*(ABz*HRRB(3)+&
                         3.D0*HRRB(9))+&
                         3.D0*HRRB(19))+&
                         HRRB(34))+&
                         HRRB(53)
      OffSet=(OA+2)*LDA+(OB+8)*LDB+CDOffSet !=
      !=(4_x,19|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*(ABz*HRRA(16)+&
                         2.D0*HRRA(31))+&
                         ABy*(ABz*(ABz*HRRA(8)+&
                         2.D0*HRRA(18))+&
                         HRRA(33))+&
                         HRRA(52)
      !=(4,19_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(16)+&
                         2.D0*HRRB(31))+&
                         ABy*(ABz*(ABz*HRRB(8)+&
                         2.D0*HRRB(18))+&
                         HRRB(33))+&
                         ABx*(ABz*(ABz*HRRB(9)+&
                         2.D0*HRRB(19))+&
                         ABy*(ABz*(ABz*HRRB(4)+&
                         2.D0*HRRB(10))+&
                         HRRB(20))+&
                         HRRB(34))+&
                         HRRB(52)
      !=(4_y,19|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*(ABz*HRRA(17)+&
                         2.D0*HRRA(32))+&
                         ABy*(ABz*(ABz*HRRA(9)+&
                         2.D0*HRRA(19))+&
                         HRRA(34))+&
                         HRRA(53)
      !=(4,19_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(20)+&
                         ABz*(-2.D0*HRR(10)+&
                         ABz*(-HRR(4)+&
                         HRRB(17))+&
                         2.D0*HRRB(32))+&
                         ABy*(ABz*(2.D0*ABz*HRRB(9)+&
                         4.D0*HRRB(19))+&
                         ABy*(ABz*(ABz*HRRB(4)+&
                         2.D0*HRRB(10))+&
                         HRRB(20))+&
                         2.D0*HRRB(34))+&
                         HRRB(53)
      !=(4_z,19|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(19)+&
                         ABz*(-2.D0*HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRA(19))+&
                         2.D0*HRRA(34))+&
                         ABy*(-HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(10))+&
                         2.D0*HRRA(20))+&
                         HRRA(35))+&
                         HRRA(55)
      !=(4,19_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(19)+&
                         ABz*(-2.D0*HRR(9)+&
                         ABz*(ABz*HRRB(9)+&
                         3.D0*HRRB(19))+&
                         3.D0*HRRB(34))+&
                         ABy*(-2.D0*HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         ABz*(ABz*HRRB(4)+&
                         3.D0*HRRB(10))+&
                         3.D0*HRRB(20))+&
                         HRRB(35))+&
                         HRRB(55)
      OffSet=(OA+0)*LDA+(OB+9)*LDB+CDOffSet !=
      !=(2_x,20|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(20)+&
                         ABz*(-3.D0*HRR(10)+&
                         ABz*(-3.D0*HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(5))+&
                         3.D0*HRRA(15))+&
                         3.D0*HRRA(30))+&
                         HRRA(51)
      !=(2,20_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*(ABz*HRRB(5)+&
                         3.D0*HRRB(15))+&
                         3.D0*HRRB(30))+&
                         ABx*(ABz*(ABz*(ABz*HRRB(2)+&
                         3.D0*HRRB(8))+&
                         3.D0*HRRB(18))+&
                         HRRB(33))+&
                         HRRB(51)
      !=(2_y,20|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*(ABz*(ABz*HRRA(6)+&
                         3.D0*HRRA(16))+&
                         3.D0*HRRA(31))+&
                         HRRA(52)
      !=(2,20_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*(ABz*HRRB(6)+&
                         3.D0*HRRB(16))+&
                         3.D0*HRRB(31))+&
                         ABy*(ABz*(ABz*(ABz*HRRB(2)+&
                         3.D0*HRRB(8))+&
                         3.D0*HRRB(18))+&
                         HRRB(33))+&
                         HRRB(52)
      !=(2_z,20|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*(ABz*(ABz*HRRA(8)+&
                         3.D0*HRRA(18))+&
                         3.D0*HRRA(33))+&
                         HRRA(54)
      !=(2,20_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-3.D0*HRR(18)+&
                         ABz*(-6.D0*HRR(8)+&
                         ABz*(-3.D0*HRR(2)+&
                         ABz*(ABz*HRRB(2)+&
                         4.D0*HRRB(8))+&
                         6.D0*HRRB(18))+&
                         4.D0*HRRB(33))+&
                         HRRB(54)
      OffSet=(OA+1)*LDA+(OB+9)*LDB+CDOffSet !=
      !=(3_x,20|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*(ABz*(ABz*HRRA(6)+&
                         3.D0*HRRA(16))+&
                         3.D0*HRRA(31))+&
                         HRRA(52)
      !=(3,20_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*(ABz*HRRB(6)+&
                         3.D0*HRRB(16))+&
                         3.D0*HRRB(31))+&
                         ABx*(ABz*(ABz*(ABz*HRRB(3)+&
                         3.D0*HRRB(9))+&
                         3.D0*HRRB(19))+&
                         HRRB(34))+&
                         HRRB(52)
      !=(3_y,20|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(20)+&
                         ABz*(-3.D0*HRR(10)+&
                         ABz*(-3.D0*HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(7))+&
                         3.D0*HRRA(17))+&
                         3.D0*HRRA(32))+&
                         HRRA(53)
      !=(3,20_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*(ABz*HRRB(7)+&
                         3.D0*HRRB(17))+&
                         3.D0*HRRB(32))+&
                         ABy*(ABz*(ABz*(ABz*HRRB(3)+&
                         3.D0*HRRB(9))+&
                         3.D0*HRRB(19))+&
                         HRRB(34))+&
                         HRRB(53)
      !=(3_z,20|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*(ABz*(ABz*HRRA(9)+&
                         3.D0*HRRA(19))+&
                         3.D0*HRRA(34))+&
                         HRRA(55)
      !=(3,20_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-3.D0*HRR(19)+&
                         ABz*(-6.D0*HRR(9)+&
                         ABz*(-3.D0*HRR(3)+&
                         ABz*(ABz*HRRB(3)+&
                         4.D0*HRRB(9))+&
                         6.D0*HRRB(19))+&
                         4.D0*HRRB(34))+&
                         HRRB(55)
      OffSet=(OA+2)*LDA+(OB+9)*LDB+CDOffSet !=
      !=(4_x,20|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*(ABz*(ABz*HRRA(8)+&
                         3.D0*HRRA(18))+&
                         3.D0*HRRA(33))+&
                         HRRA(54)
      !=(4,20_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*(ABz*HRRB(8)+&
                         3.D0*HRRB(18))+&
                         3.D0*HRRB(33))+&
                         ABx*(ABz*(ABz*(ABz*HRRB(4)+&
                         3.D0*HRRB(10))+&
                         3.D0*HRRB(20))+&
                         HRRB(35))+&
                         HRRB(54)
      !=(4_y,20|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*(ABz*(ABz*HRRA(9)+&
                         3.D0*HRRA(19))+&
                         3.D0*HRRA(34))+&
                         HRRA(55)
      !=(4,20_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*(ABz*HRRB(9)+&
                         3.D0*HRRB(19))+&
                         3.D0*HRRB(34))+&
                         ABy*(ABz*(ABz*(ABz*HRRB(4)+&
                         3.D0*HRRB(10))+&
                         3.D0*HRRB(20))+&
                         HRRB(35))+&
                         HRRB(55)
      !=(4_z,20|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(20)+&
                         ABz*(-3.D0*HRR(10)+&
                         ABz*(-3.D0*HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(10))+&
                         3.D0*HRRA(20))+&
                         3.D0*HRRA(35))+&
                         HRRA(56)
      !=(4,20_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-3.D0*HRR(20)+&
                         ABz*(-6.D0*HRR(10)+&
                         ABz*(-3.D0*HRR(4)+&
                         ABz*(ABz*HRRB(4)+&
                         4.D0*HRRB(10))+&
                         6.D0*HRRB(20))+&
                         4.D0*HRRB(35))+&
                         HRRB(56)
    END SUBROUTINE BraHRR310ab
    SUBROUTINE BraHRR310cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,Cart,HRR,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,Cart,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(2,11|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*(ABx*HRR(2)+&
                                3.D0*HRR(5))+&
                                3.D0*HRR(11))+&
                                HRR(21)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(3,11|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*(ABx*HRR(3)+&
                                3.D0*HRR(6))+&
                                3.D0*HRR(12))+&
                                HRR(22)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(4,11|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*(ABx*HRR(4)+&
                                3.D0*HRR(8))+&
                                3.D0*HRR(15))+&
                                HRR(26)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(2,12|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(11)+&
                                ABx*(2.D0*ABy*HRR(5)+&
                                ABx*(ABy*HRR(2)+&
                                HRR(6))+&
                                2.D0*HRR(12))+&
                                HRR(22)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(3,12|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(12)+&
                                ABx*(2.D0*ABy*HRR(6)+&
                                ABx*(ABy*HRR(3)+&
                                HRR(7))+&
                                2.D0*HRR(13))+&
                                HRR(23)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(4,12|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(15)+&
                                ABx*(2.D0*ABy*HRR(8)+&
                                ABx*(ABy*HRR(4)+&
                                HRR(9))+&
                                2.D0*HRR(16))+&
                                HRR(27)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(2,13|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(5)+&
                                2.D0*HRR(12))+&
                                ABx*(ABy*(ABy*HRR(2)+&
                                2.D0*HRR(6))+&
                                HRR(13))+&
                                HRR(23)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(3,13|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(6)+&
                                2.D0*HRR(13))+&
                                ABx*(ABy*(ABy*HRR(3)+&
                                2.D0*HRR(7))+&
                                HRR(14))+&
                                HRR(24)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(4,13|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(8)+&
                                2.D0*HRR(16))+&
                                ABx*(ABy*(ABy*HRR(4)+&
                                2.D0*HRR(9))+&
                                HRR(17))+&
                                HRR(28)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(2,14|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*(ABy*HRR(2)+&
                                3.D0*HRR(6))+&
                                3.D0*HRR(13))+&
                                HRR(24)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(3,14|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*(ABy*HRR(3)+&
                                3.D0*HRR(7))+&
                                3.D0*HRR(14))+&
                                HRR(25)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(4,14|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*(ABy*HRR(4)+&
                                3.D0*HRR(9))+&
                                3.D0*HRR(17))+&
                                HRR(29)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(2,15|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(11)+&
                                ABx*(2.D0*ABz*HRR(5)+&
                                ABx*(ABz*HRR(2)+&
                                HRR(8))+&
                                2.D0*HRR(15))+&
                                HRR(26)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(3,15|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(12)+&
                                ABx*(2.D0*ABz*HRR(6)+&
                                ABx*(ABz*HRR(3)+&
                                HRR(9))+&
                                2.D0*HRR(16))+&
                                HRR(27)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(4,15|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(15)+&
                                ABx*(2.D0*ABz*HRR(8)+&
                                ABx*(ABz*HRR(4)+&
                                HRR(10))+&
                                2.D0*HRR(18))+&
                                HRR(30)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(2,16|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(12)+&
                                ABy*(ABz*HRR(5)+&
                                HRR(15))+&
                                ABx*(ABz*HRR(6)+&
                                ABy*(ABz*HRR(2)+&
                                HRR(8))+&
                                HRR(16))+&
                                HRR(27)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(3,16|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(13)+&
                                ABy*(ABz*HRR(6)+&
                                HRR(16))+&
                                ABx*(ABz*HRR(7)+&
                                ABy*(ABz*HRR(3)+&
                                HRR(9))+&
                                HRR(17))+&
                                HRR(28)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(4,16|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(16)+&
                                ABy*(ABz*HRR(8)+&
                                HRR(18))+&
                                ABx*(ABz*HRR(9)+&
                                ABy*(ABz*HRR(4)+&
                                HRR(10))+&
                                HRR(19))+&
                                HRR(31)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+6)*LDB+CDOffSet !=
      !=(2,17|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(13)+&
                                ABy*(2.D0*ABz*HRR(6)+&
                                ABy*(ABz*HRR(2)+&
                                HRR(8))+&
                                2.D0*HRR(16))+&
                                HRR(28)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+6)*LDB+CDOffSet !=
      !=(3,17|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(14)+&
                                ABy*(2.D0*ABz*HRR(7)+&
                                ABy*(ABz*HRR(3)+&
                                HRR(9))+&
                                2.D0*HRR(17))+&
                                HRR(29)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+6)*LDB+CDOffSet !=
      !=(4,17|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(17)+&
                                ABy*(2.D0*ABz*HRR(9)+&
                                ABy*(ABz*HRR(4)+&
                                HRR(10))+&
                                2.D0*HRR(19))+&
                                HRR(32)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+7)*LDB+CDOffSet !=
      !=(2,18|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(5)+&
                                2.D0*HRR(15))+&
                                ABx*(ABz*(ABz*HRR(2)+&
                                2.D0*HRR(8))+&
                                HRR(18))+&
                                HRR(30)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+7)*LDB+CDOffSet !=
      !=(3,18|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(6)+&
                                2.D0*HRR(16))+&
                                ABx*(ABz*(ABz*HRR(3)+&
                                2.D0*HRR(9))+&
                                HRR(19))+&
                                HRR(31)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+7)*LDB+CDOffSet !=
      !=(4,18|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(8)+&
                                2.D0*HRR(18))+&
                                ABx*(ABz*(ABz*HRR(4)+&
                                2.D0*HRR(10))+&
                                HRR(20))+&
                                HRR(33)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+8)*LDB+CDOffSet !=
      !=(2,19|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(6)+&
                                2.D0*HRR(16))+&
                                ABy*(ABz*(ABz*HRR(2)+&
                                2.D0*HRR(8))+&
                                HRR(18))+&
                                HRR(31)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+8)*LDB+CDOffSet !=
      !=(3,19|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(7)+&
                                2.D0*HRR(17))+&
                                ABy*(ABz*(ABz*HRR(3)+&
                                2.D0*HRR(9))+&
                                HRR(19))+&
                                HRR(32)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+8)*LDB+CDOffSet !=
      !=(4,19|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(9)+&
                                2.D0*HRR(19))+&
                                ABy*(ABz*(ABz*HRR(4)+&
                                2.D0*HRR(10))+&
                                HRR(20))+&
                                HRR(34)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+9)*LDB+CDOffSet !=
      !=(2,20|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*(ABz*HRR(2)+&
                                3.D0*HRR(8))+&
                                3.D0*HRR(18))+&
                                HRR(33)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+9)*LDB+CDOffSet !=
      !=(3,20|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*(ABz*HRR(3)+&
                                3.D0*HRR(9))+&
                                3.D0*HRR(19))+&
                                HRR(34)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+9)*LDB+CDOffSet !=
      !=(4,20|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*(ABz*HRR(4)+&
                                3.D0*HRR(10))+&
                                3.D0*HRR(20))+&
                                HRR(35)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
    END SUBROUTINE BraHRR310cd
