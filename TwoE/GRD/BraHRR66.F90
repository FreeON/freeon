    SUBROUTINE BraHRR66ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR,HRRA,HRRB,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*),HRRA(*),HRRB(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(5_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(11)+&
                         ABx*(-4.D0*HRR(5)+&
                         ABx*(-2.D0*HRR(2)+&
                         HRRA(11))+&
                         2.D0*HRRA(21))+&
                         HRRA(36)
      !=(5,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(11)+&
                         ABx*(-2.D0*HRR(5)+&
                         ABx*(ABx*HRRB(5)+&
                         3.D0*HRRB(11))+&
                         3.D0*HRRB(21))+&
                         HRRB(36)
      !=(5_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*(ABx*HRRA(12)+&
                         2.D0*HRRA(22))+&
                         HRRA(37)
      !=(5,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(21)+&
                         ABx*(2.D0*ABy*HRRB(11)+&
                         ABx*(ABy*HRRB(5)+&
                         HRRB(12))+&
                         2.D0*HRRB(22))+&
                         HRRB(37)
      !=(5_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*(ABx*HRRA(15)+&
                         2.D0*HRRA(26))+&
                         HRRA(42)
      !=(5,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(21)+&
                         ABx*(2.D0*ABz*HRRB(11)+&
                         ABx*(ABz*HRRB(5)+&
                         HRRB(15))+&
                         2.D0*HRRB(26))+&
                         HRRB(42)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(6_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(12)+&
                         ABx*(-2.D0*HRR(6)+&
                         ABx*(-HRR(3)+&
                         HRRA(12))+&
                         2.D0*HRRA(22))+&
                         HRRA(37)
      !=(6,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(12)+&
                         ABx*(-2.D0*HRR(6)+&
                         ABx*(ABx*HRRB(6)+&
                         3.D0*HRRB(12))+&
                         3.D0*HRRB(22))+&
                         HRRB(37)
      !=(6_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(11)+&
                         ABx*(-2.D0*HRR(5)+&
                         ABx*(-HRR(2)+&
                         HRRA(13))+&
                         2.D0*HRRA(23))+&
                         HRRA(38)
      !=(6,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(22)+&
                         ABx*(2.D0*ABy*HRRB(12)+&
                         ABx*(ABy*HRRB(6)+&
                         HRRB(13))+&
                         2.D0*HRRB(23))+&
                         HRRB(38)
      !=(6_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*(ABx*HRRA(16)+&
                         2.D0*HRRA(27))+&
                         HRRA(43)
      !=(6,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(22)+&
                         ABx*(2.D0*ABz*HRRB(12)+&
                         ABx*(ABz*HRRB(6)+&
                         HRRB(16))+&
                         2.D0*HRRB(27))+&
                         HRRB(43)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(7_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*(ABx*HRRA(13)+&
                         2.D0*HRRA(23))+&
                         HRRA(38)
      !=(7,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(13)+&
                         ABx*(-2.D0*HRR(7)+&
                         ABx*(ABx*HRRB(7)+&
                         3.D0*HRRB(13))+&
                         3.D0*HRRB(23))+&
                         HRRB(38)
      !=(7_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(12)+&
                         ABx*(-4.D0*HRR(6)+&
                         ABx*(-2.D0*HRR(3)+&
                         HRRA(14))+&
                         2.D0*HRRA(24))+&
                         HRRA(39)
      !=(7,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(23)+&
                         ABx*(2.D0*ABy*HRRB(13)+&
                         ABx*(ABy*HRRB(7)+&
                         HRRB(14))+&
                         2.D0*HRRB(24))+&
                         HRRB(39)
      !=(7_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*(ABx*HRRA(17)+&
                         2.D0*HRRA(28))+&
                         HRRA(44)
      !=(7,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(23)+&
                         ABx*(2.D0*ABz*HRRB(13)+&
                         ABx*(ABz*HRRB(7)+&
                         HRRB(17))+&
                         2.D0*HRRB(28))+&
                         HRRB(44)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(8_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(15)+&
                         ABx*(-2.D0*HRR(8)+&
                         ABx*(-HRR(4)+&
                         HRRA(15))+&
                         2.D0*HRRA(26))+&
                         HRRA(42)
      !=(8,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(15)+&
                         ABx*(-2.D0*HRR(8)+&
                         ABx*(ABx*HRRB(8)+&
                         3.D0*HRRB(15))+&
                         3.D0*HRRB(26))+&
                         HRRB(42)
      !=(8_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*(ABx*HRRA(16)+&
                         2.D0*HRRA(27))+&
                         HRRA(43)
      !=(8,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(26)+&
                         ABx*(2.D0*ABy*HRRB(15)+&
                         ABx*(ABy*HRRB(8)+&
                         HRRB(16))+&
                         2.D0*HRRB(27))+&
                         HRRB(43)
      !=(8_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(11)+&
                         ABx*(-2.D0*HRR(5)+&
                         ABx*(-HRR(2)+&
                         HRRA(18))+&
                         2.D0*HRRA(30))+&
                         HRRA(47)
      !=(8,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(26)+&
                         ABx*(2.D0*ABz*HRRB(15)+&
                         ABx*(ABz*HRRB(8)+&
                         HRRB(18))+&
                         2.D0*HRRB(30))+&
                         HRRB(47)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(9_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*(ABx*HRRA(16)+&
                         2.D0*HRRA(27))+&
                         HRRA(43)
      !=(9,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(16)+&
                         ABx*(-2.D0*HRR(9)+&
                         ABx*(ABx*HRRB(9)+&
                         3.D0*HRRB(16))+&
                         3.D0*HRRB(27))+&
                         HRRB(43)
      !=(9_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(15)+&
                         ABx*(-2.D0*HRR(8)+&
                         ABx*(-HRR(4)+&
                         HRRA(17))+&
                         2.D0*HRRA(28))+&
                         HRRA(44)
      !=(9,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(27)+&
                         ABx*(2.D0*ABy*HRRB(16)+&
                         ABx*(ABy*HRRB(9)+&
                         HRRB(17))+&
                         2.D0*HRRB(28))+&
                         HRRB(44)
      !=(9_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(12)+&
                         ABx*(-2.D0*HRR(6)+&
                         ABx*(-HRR(3)+&
                         HRRA(19))+&
                         2.D0*HRRA(31))+&
                         HRRA(48)
      !=(9,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(27)+&
                         ABx*(2.D0*ABz*HRRB(16)+&
                         ABx*(ABz*HRRB(9)+&
                         HRRB(19))+&
                         2.D0*HRRB(31))+&
                         HRRB(48)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(10_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*(ABx*HRRA(18)+&
                         2.D0*HRRA(30))+&
                         HRRA(47)
      !=(10,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(18)+&
                         ABx*(-2.D0*HRR(10)+&
                         ABx*(ABx*HRRB(10)+&
                         3.D0*HRRB(18))+&
                         3.D0*HRRB(30))+&
                         HRRB(47)
      !=(10_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*(ABx*HRRA(19)+&
                         2.D0*HRRA(31))+&
                         HRRA(48)
      !=(10,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(30)+&
                         ABx*(2.D0*ABy*HRRB(18)+&
                         ABx*(ABy*HRRB(10)+&
                         HRRB(19))+&
                         2.D0*HRRB(31))+&
                         HRRB(48)
      !=(10_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(15)+&
                         ABx*(-4.D0*HRR(8)+&
                         ABx*(-2.D0*HRR(4)+&
                         HRRA(20))+&
                         2.D0*HRRA(33))+&
                         HRRA(51)
      !=(10,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(30)+&
                         ABx*(2.D0*ABz*HRRB(18)+&
                         ABx*(ABz*HRRB(10)+&
                         HRRB(20))+&
                         2.D0*HRRB(33))+&
                         HRRB(51)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(5_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(12)+&
                         ABy*(-2.D0*HRR(5)+&
                         HRRA(21))+&
                         ABx*(-2.D0*HRR(6)+&
                         ABy*(-2.D0*HRR(2)+&
                         HRRA(11))+&
                         HRRA(22))+&
                         HRRA(37)
      !=(5,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(12)+&
                         ABy*(-HRR(5)+&
                         HRRB(21))+&
                         ABx*(2.D0*ABy*HRRB(11)+&
                         ABx*(ABy*HRRB(5)+&
                         HRRB(12))+&
                         2.D0*HRRB(22))+&
                         HRRB(37)
      !=(5_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(22)+&
                         ABx*(ABy*HRRA(12)+&
                         HRRA(23))+&
                         HRRA(38)
      !=(5,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(11)+&
                         ABy*(ABy*HRRB(11)+&
                         2.D0*HRRB(22))+&
                         ABx*(-HRR(5)+&
                         ABy*(ABy*HRRB(5)+&
                         2.D0*HRRB(12))+&
                         HRRB(23))+&
                         HRRB(38)
      !=(5_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(26)+&
                         ABx*(ABy*HRRA(15)+&
                         HRRA(27))+&
                         HRRA(43)
      !=(5,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(22)+&
                         ABy*(ABz*HRRB(11)+&
                         HRRB(26))+&
                         ABx*(ABz*HRRB(12)+&
                         ABy*(ABz*HRRB(5)+&
                         HRRB(15))+&
                         HRRB(27))+&
                         HRRB(43)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(6_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(13)+&
                         ABy*(-HRR(6)+&
                         HRRA(22))+&
                         ABx*(-HRR(7)+&
                         ABy*(-HRR(3)+&
                         HRRA(12))+&
                         HRRA(23))+&
                         HRRA(38)
      !=(6,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(13)+&
                         ABy*(-HRR(6)+&
                         HRRB(22))+&
                         ABx*(2.D0*ABy*HRRB(12)+&
                         ABx*(ABy*HRRB(6)+&
                         HRRB(13))+&
                         2.D0*HRRB(23))+&
                         HRRB(38)
      !=(6_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(12)+&
                         ABy*(-HRR(5)+&
                         HRRA(23))+&
                         ABx*(-HRR(6)+&
                         ABy*(-HRR(2)+&
                         HRRA(13))+&
                         HRRA(24))+&
                         HRRA(39)
      !=(6,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(12)+&
                         ABy*(ABy*HRRB(12)+&
                         2.D0*HRRB(23))+&
                         ABx*(-HRR(6)+&
                         ABy*(ABy*HRRB(6)+&
                         2.D0*HRRB(13))+&
                         HRRB(24))+&
                         HRRB(39)
      !=(6_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(27)+&
                         ABx*(ABy*HRRA(16)+&
                         HRRA(28))+&
                         HRRA(44)
      !=(6,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(23)+&
                         ABy*(ABz*HRRB(12)+&
                         HRRB(27))+&
                         ABx*(ABz*HRRB(13)+&
                         ABy*(ABz*HRRB(6)+&
                         HRRB(16))+&
                         HRRB(28))+&
                         HRRB(44)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(7_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(23)+&
                         ABx*(ABy*HRRA(13)+&
                         HRRA(24))+&
                         HRRA(39)
      !=(7,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(14)+&
                         ABy*(-HRR(7)+&
                         HRRB(23))+&
                         ABx*(2.D0*ABy*HRRB(13)+&
                         ABx*(ABy*HRRB(7)+&
                         HRRB(14))+&
                         2.D0*HRRB(24))+&
                         HRRB(39)
      !=(7_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(13)+&
                         ABy*(-2.D0*HRR(6)+&
                         HRRA(24))+&
                         ABx*(-2.D0*HRR(7)+&
                         ABy*(-2.D0*HRR(3)+&
                         HRRA(14))+&
                         HRRA(25))+&
                         HRRA(40)
      !=(7,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(13)+&
                         ABy*(ABy*HRRB(13)+&
                         2.D0*HRRB(24))+&
                         ABx*(-HRR(7)+&
                         ABy*(ABy*HRRB(7)+&
                         2.D0*HRRB(14))+&
                         HRRB(25))+&
                         HRRB(40)
      !=(7_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(28)+&
                         ABx*(ABy*HRRA(17)+&
                         HRRA(29))+&
                         HRRA(45)
      !=(7,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(24)+&
                         ABy*(ABz*HRRB(13)+&
                         HRRB(28))+&
                         ABx*(ABz*HRRB(14)+&
                         ABy*(ABz*HRRB(7)+&
                         HRRB(17))+&
                         HRRB(29))+&
                         HRRB(45)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(8_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(16)+&
                         ABy*(-HRR(8)+&
                         HRRA(26))+&
                         ABx*(-HRR(9)+&
                         ABy*(-HRR(4)+&
                         HRRA(15))+&
                         HRRA(27))+&
                         HRRA(43)
      !=(8,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(16)+&
                         ABy*(-HRR(8)+&
                         HRRB(26))+&
                         ABx*(2.D0*ABy*HRRB(15)+&
                         ABx*(ABy*HRRB(8)+&
                         HRRB(16))+&
                         2.D0*HRRB(27))+&
                         HRRB(43)
      !=(8_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(27)+&
                         ABx*(ABy*HRRA(16)+&
                         HRRA(28))+&
                         HRRA(44)
      !=(8,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(15)+&
                         ABy*(ABy*HRRB(15)+&
                         2.D0*HRRB(27))+&
                         ABx*(-HRR(8)+&
                         ABy*(ABy*HRRB(8)+&
                         2.D0*HRRB(16))+&
                         HRRB(28))+&
                         HRRB(44)
      !=(8_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(12)+&
                         ABy*(-HRR(5)+&
                         HRRA(30))+&
                         ABx*(-HRR(6)+&
                         ABy*(-HRR(2)+&
                         HRRA(18))+&
                         HRRA(31))+&
                         HRRA(48)
      !=(8,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(27)+&
                         ABy*(ABz*HRRB(15)+&
                         HRRB(30))+&
                         ABx*(ABz*HRRB(16)+&
                         ABy*(ABz*HRRB(8)+&
                         HRRB(18))+&
                         HRRB(31))+&
                         HRRB(48)
      OffSet=(OA+4)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(9_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(27)+&
                         ABx*(ABy*HRRA(16)+&
                         HRRA(28))+&
                         HRRA(44)
      !=(9,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(17)+&
                         ABy*(-HRR(9)+&
                         HRRB(27))+&
                         ABx*(2.D0*ABy*HRRB(16)+&
                         ABx*(ABy*HRRB(9)+&
                         HRRB(17))+&
                         2.D0*HRRB(28))+&
                         HRRB(44)
      !=(9_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(16)+&
                         ABy*(-HRR(8)+&
                         HRRA(28))+&
                         ABx*(-HRR(9)+&
                         ABy*(-HRR(4)+&
                         HRRA(17))+&
                         HRRA(29))+&
                         HRRA(45)
      !=(9,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(16)+&
                         ABy*(ABy*HRRB(16)+&
                         2.D0*HRRB(28))+&
                         ABx*(-HRR(9)+&
                         ABy*(ABy*HRRB(9)+&
                         2.D0*HRRB(17))+&
                         HRRB(29))+&
                         HRRB(45)
      !=(9_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(13)+&
                         ABy*(-HRR(6)+&
                         HRRA(31))+&
                         ABx*(-HRR(7)+&
                         ABy*(-HRR(3)+&
                         HRRA(19))+&
                         HRRA(32))+&
                         HRRA(49)
      !=(9,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(28)+&
                         ABy*(ABz*HRRB(16)+&
                         HRRB(31))+&
                         ABx*(ABz*HRRB(17)+&
                         ABy*(ABz*HRRB(9)+&
                         HRRB(19))+&
                         HRRB(32))+&
                         HRRB(49)
      OffSet=(OA+5)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(10_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(30)+&
                         ABx*(ABy*HRRA(18)+&
                         HRRA(31))+&
                         HRRA(48)
      !=(10,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(19)+&
                         ABy*(-HRR(10)+&
                         HRRB(30))+&
                         ABx*(2.D0*ABy*HRRB(18)+&
                         ABx*(ABy*HRRB(10)+&
                         HRRB(19))+&
                         2.D0*HRRB(31))+&
                         HRRB(48)
      !=(10_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(31)+&
                         ABx*(ABy*HRRA(19)+&
                         HRRA(32))+&
                         HRRA(49)
      !=(10,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(18)+&
                         ABy*(ABy*HRRB(18)+&
                         2.D0*HRRB(31))+&
                         ABx*(-HRR(10)+&
                         ABy*(ABy*HRRB(10)+&
                         2.D0*HRRB(19))+&
                         HRRB(32))+&
                         HRRB(49)
      !=(10_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(16)+&
                         ABy*(-2.D0*HRR(8)+&
                         HRRA(33))+&
                         ABx*(-2.D0*HRR(9)+&
                         ABy*(-2.D0*HRR(4)+&
                         HRRA(20))+&
                         HRRA(34))+&
                         HRRA(52)
      !=(10,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(31)+&
                         ABy*(ABz*HRRB(18)+&
                         HRRB(33))+&
                         ABx*(ABz*HRRB(19)+&
                         ABy*(ABz*HRRB(10)+&
                         HRRB(20))+&
                         HRRB(34))+&
                         HRRB(52)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(5_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(13)+&
                         ABy*(-4.D0*HRR(6)+&
                         ABy*(-2.D0*HRR(2)+&
                         HRRA(11))+&
                         2.D0*HRRA(22))+&
                         HRRA(38)
      !=(5,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(11)+&
                         2.D0*HRRB(22))+&
                         ABx*(ABy*(ABy*HRRB(5)+&
                         2.D0*HRRB(12))+&
                         HRRB(23))+&
                         HRRB(38)
      !=(5_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*(ABy*HRRA(12)+&
                         2.D0*HRRA(23))+&
                         HRRA(39)
      !=(5,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(12)+&
                         ABy*(-2.D0*HRR(5)+&
                         ABy*(ABy*HRRB(5)+&
                         3.D0*HRRB(12))+&
                         3.D0*HRRB(23))+&
                         HRRB(39)
      !=(5_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*(ABy*HRRA(15)+&
                         2.D0*HRRA(27))+&
                         HRRA(44)
      !=(5,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(23)+&
                         ABy*(2.D0*ABz*HRRB(12)+&
                         ABy*(ABz*HRRB(5)+&
                         HRRB(15))+&
                         2.D0*HRRB(27))+&
                         HRRB(44)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(6_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(14)+&
                         ABy*(-2.D0*HRR(7)+&
                         ABy*(-HRR(3)+&
                         HRRA(12))+&
                         2.D0*HRRA(23))+&
                         HRRA(39)
      !=(6,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(12)+&
                         2.D0*HRRB(23))+&
                         ABx*(ABy*(ABy*HRRB(6)+&
                         2.D0*HRRB(13))+&
                         HRRB(24))+&
                         HRRB(39)
      !=(6_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(13)+&
                         ABy*(-2.D0*HRR(6)+&
                         ABy*(-HRR(2)+&
                         HRRA(13))+&
                         2.D0*HRRA(24))+&
                         HRRA(40)
      !=(6,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(13)+&
                         ABy*(-2.D0*HRR(6)+&
                         ABy*(ABy*HRRB(6)+&
                         3.D0*HRRB(13))+&
                         3.D0*HRRB(24))+&
                         HRRB(40)
      !=(6_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*(ABy*HRRA(16)+&
                         2.D0*HRRA(28))+&
                         HRRA(45)
      !=(6,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(24)+&
                         ABy*(2.D0*ABz*HRRB(13)+&
                         ABy*(ABz*HRRB(6)+&
                         HRRB(16))+&
                         2.D0*HRRB(28))+&
                         HRRB(45)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(7_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*(ABy*HRRA(13)+&
                         2.D0*HRRA(24))+&
                         HRRA(40)
      !=(7,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(13)+&
                         2.D0*HRRB(24))+&
                         ABx*(ABy*(ABy*HRRB(7)+&
                         2.D0*HRRB(14))+&
                         HRRB(25))+&
                         HRRB(40)
      !=(7_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(14)+&
                         ABy*(-4.D0*HRR(7)+&
                         ABy*(-2.D0*HRR(3)+&
                         HRRA(14))+&
                         2.D0*HRRA(25))+&
                         HRRA(41)
      !=(7,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(14)+&
                         ABy*(-2.D0*HRR(7)+&
                         ABy*(ABy*HRRB(7)+&
                         3.D0*HRRB(14))+&
                         3.D0*HRRB(25))+&
                         HRRB(41)
      !=(7_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*(ABy*HRRA(17)+&
                         2.D0*HRRA(29))+&
                         HRRA(46)
      !=(7,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(25)+&
                         ABy*(2.D0*ABz*HRRB(14)+&
                         ABy*(ABz*HRRB(7)+&
                         HRRB(17))+&
                         2.D0*HRRB(29))+&
                         HRRB(46)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(8_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(17)+&
                         ABy*(-2.D0*HRR(9)+&
                         ABy*(-HRR(4)+&
                         HRRA(15))+&
                         2.D0*HRRA(27))+&
                         HRRA(44)
      !=(8,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(15)+&
                         2.D0*HRRB(27))+&
                         ABx*(ABy*(ABy*HRRB(8)+&
                         2.D0*HRRB(16))+&
                         HRRB(28))+&
                         HRRB(44)
      !=(8_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*(ABy*HRRA(16)+&
                         2.D0*HRRA(28))+&
                         HRRA(45)
      !=(8,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(16)+&
                         ABy*(-2.D0*HRR(8)+&
                         ABy*(ABy*HRRB(8)+&
                         3.D0*HRRB(16))+&
                         3.D0*HRRB(28))+&
                         HRRB(45)
      !=(8_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(13)+&
                         ABy*(-2.D0*HRR(6)+&
                         ABy*(-HRR(2)+&
                         HRRA(18))+&
                         2.D0*HRRA(31))+&
                         HRRA(49)
      !=(8,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(28)+&
                         ABy*(2.D0*ABz*HRRB(16)+&
                         ABy*(ABz*HRRB(8)+&
                         HRRB(18))+&
                         2.D0*HRRB(31))+&
                         HRRB(49)
      OffSet=(OA+4)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(9_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*(ABy*HRRA(16)+&
                         2.D0*HRRA(28))+&
                         HRRA(45)
      !=(9,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(16)+&
                         2.D0*HRRB(28))+&
                         ABx*(ABy*(ABy*HRRB(9)+&
                         2.D0*HRRB(17))+&
                         HRRB(29))+&
                         HRRB(45)
      !=(9_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(17)+&
                         ABy*(-2.D0*HRR(9)+&
                         ABy*(-HRR(4)+&
                         HRRA(17))+&
                         2.D0*HRRA(29))+&
                         HRRA(46)
      !=(9,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(17)+&
                         ABy*(-2.D0*HRR(9)+&
                         ABy*(ABy*HRRB(9)+&
                         3.D0*HRRB(17))+&
                         3.D0*HRRB(29))+&
                         HRRB(46)
      !=(9_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(14)+&
                         ABy*(-2.D0*HRR(7)+&
                         ABy*(-HRR(3)+&
                         HRRA(19))+&
                         2.D0*HRRA(32))+&
                         HRRA(50)
      !=(9,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(29)+&
                         ABy*(2.D0*ABz*HRRB(17)+&
                         ABy*(ABz*HRRB(9)+&
                         HRRB(19))+&
                         2.D0*HRRB(32))+&
                         HRRB(50)
      OffSet=(OA+5)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(10_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*(ABy*HRRA(18)+&
                         2.D0*HRRA(31))+&
                         HRRA(49)
      !=(10,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(18)+&
                         2.D0*HRRB(31))+&
                         ABx*(ABy*(ABy*HRRB(10)+&
                         2.D0*HRRB(19))+&
                         HRRB(32))+&
                         HRRB(49)
      !=(10_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*(ABy*HRRA(19)+&
                         2.D0*HRRA(32))+&
                         HRRA(50)
      !=(10,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(19)+&
                         ABy*(-2.D0*HRR(10)+&
                         ABy*(ABy*HRRB(10)+&
                         3.D0*HRRB(19))+&
                         3.D0*HRRB(32))+&
                         HRRB(50)
      !=(10_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(17)+&
                         ABy*(-4.D0*HRR(9)+&
                         ABy*(-2.D0*HRR(4)+&
                         HRRA(20))+&
                         2.D0*HRRA(34))+&
                         HRRA(53)
      !=(10,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(32)+&
                         ABy*(2.D0*ABz*HRRB(19)+&
                         ABy*(ABz*HRRB(10)+&
                         HRRB(20))+&
                         2.D0*HRRB(34))+&
                         HRRB(53)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(5_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(15)+&
                         ABz*(-2.D0*HRR(5)+&
                         HRRA(21))+&
                         ABx*(-2.D0*HRR(8)+&
                         ABz*(-2.D0*HRR(2)+&
                         HRRA(11))+&
                         HRRA(26))+&
                         HRRA(42)
      !=(5,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(15)+&
                         ABz*(-HRR(5)+&
                         HRRB(21))+&
                         ABx*(2.D0*ABz*HRRB(11)+&
                         ABx*(ABz*HRRB(5)+&
                         HRRB(15))+&
                         2.D0*HRRB(26))+&
                         HRRB(42)
      !=(5_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(22)+&
                         ABx*(ABz*HRRA(12)+&
                         HRRA(27))+&
                         HRRA(43)
      !=(5,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(22)+&
                         ABy*(ABz*HRRB(11)+&
                         HRRB(26))+&
                         ABx*(ABz*HRRB(12)+&
                         ABy*(ABz*HRRB(5)+&
                         HRRB(15))+&
                         HRRB(27))+&
                         HRRB(43)
      !=(5_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(26)+&
                         ABx*(ABz*HRRA(15)+&
                         HRRA(30))+&
                         HRRA(47)
      !=(5,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(11)+&
                         ABz*(ABz*HRRB(11)+&
                         2.D0*HRRB(26))+&
                         ABx*(-HRR(5)+&
                         ABz*(ABz*HRRB(5)+&
                         2.D0*HRRB(15))+&
                         HRRB(30))+&
                         HRRB(47)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(6_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(16)+&
                         ABz*(-HRR(6)+&
                         HRRA(22))+&
                         ABx*(-HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRA(12))+&
                         HRRA(27))+&
                         HRRA(43)
      !=(6,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(16)+&
                         ABz*(-HRR(6)+&
                         HRRB(22))+&
                         ABx*(2.D0*ABz*HRRB(12)+&
                         ABx*(ABz*HRRB(6)+&
                         HRRB(16))+&
                         2.D0*HRRB(27))+&
                         HRRB(43)
      !=(6_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(15)+&
                         ABz*(-HRR(5)+&
                         HRRA(23))+&
                         ABx*(-HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRA(13))+&
                         HRRA(28))+&
                         HRRA(44)
      !=(6,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(23)+&
                         ABy*(ABz*HRRB(12)+&
                         HRRB(27))+&
                         ABx*(ABz*HRRB(13)+&
                         ABy*(ABz*HRRB(6)+&
                         HRRB(16))+&
                         HRRB(28))+&
                         HRRB(44)
      !=(6_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(27)+&
                         ABx*(ABz*HRRA(16)+&
                         HRRA(31))+&
                         HRRA(48)
      !=(6,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(12)+&
                         ABz*(ABz*HRRB(12)+&
                         2.D0*HRRB(27))+&
                         ABx*(-HRR(6)+&
                         ABz*(ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         HRRB(31))+&
                         HRRB(48)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(7_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(23)+&
                         ABx*(ABz*HRRA(13)+&
                         HRRA(28))+&
                         HRRA(44)
      !=(7,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(17)+&
                         ABz*(-HRR(7)+&
                         HRRB(23))+&
                         ABx*(2.D0*ABz*HRRB(13)+&
                         ABx*(ABz*HRRB(7)+&
                         HRRB(17))+&
                         2.D0*HRRB(28))+&
                         HRRB(44)
      !=(7_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(16)+&
                         ABz*(-2.D0*HRR(6)+&
                         HRRA(24))+&
                         ABx*(-2.D0*HRR(9)+&
                         ABz*(-2.D0*HRR(3)+&
                         HRRA(14))+&
                         HRRA(29))+&
                         HRRA(45)
      !=(7,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(24)+&
                         ABy*(ABz*HRRB(13)+&
                         HRRB(28))+&
                         ABx*(ABz*HRRB(14)+&
                         ABy*(ABz*HRRB(7)+&
                         HRRB(17))+&
                         HRRB(29))+&
                         HRRB(45)
      !=(7_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(28)+&
                         ABx*(ABz*HRRA(17)+&
                         HRRA(32))+&
                         HRRA(49)
      !=(7,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(13)+&
                         ABz*(ABz*HRRB(13)+&
                         2.D0*HRRB(28))+&
                         ABx*(-HRR(7)+&
                         ABz*(ABz*HRRB(7)+&
                         2.D0*HRRB(17))+&
                         HRRB(32))+&
                         HRRB(49)
      OffSet=(OA+3)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(8_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(18)+&
                         ABz*(-HRR(8)+&
                         HRRA(26))+&
                         ABx*(-HRR(10)+&
                         ABz*(-HRR(4)+&
                         HRRA(15))+&
                         HRRA(30))+&
                         HRRA(47)
      !=(8,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(18)+&
                         ABz*(-HRR(8)+&
                         HRRB(26))+&
                         ABx*(2.D0*ABz*HRRB(15)+&
                         ABx*(ABz*HRRB(8)+&
                         HRRB(18))+&
                         2.D0*HRRB(30))+&
                         HRRB(47)
      !=(8_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(27)+&
                         ABx*(ABz*HRRA(16)+&
                         HRRA(31))+&
                         HRRA(48)
      !=(8,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(27)+&
                         ABy*(ABz*HRRB(15)+&
                         HRRB(30))+&
                         ABx*(ABz*HRRB(16)+&
                         ABy*(ABz*HRRB(8)+&
                         HRRB(18))+&
                         HRRB(31))+&
                         HRRB(48)
      !=(8_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(15)+&
                         ABz*(-HRR(5)+&
                         HRRA(30))+&
                         ABx*(-HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRA(18))+&
                         HRRA(33))+&
                         HRRA(51)
      !=(8,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(15)+&
                         ABz*(ABz*HRRB(15)+&
                         2.D0*HRRB(30))+&
                         ABx*(-HRR(8)+&
                         ABz*(ABz*HRRB(8)+&
                         2.D0*HRRB(18))+&
                         HRRB(33))+&
                         HRRB(51)
      OffSet=(OA+4)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(9_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(27)+&
                         ABx*(ABz*HRRA(16)+&
                         HRRA(31))+&
                         HRRA(48)
      !=(9,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(19)+&
                         ABz*(-HRR(9)+&
                         HRRB(27))+&
                         ABx*(2.D0*ABz*HRRB(16)+&
                         ABx*(ABz*HRRB(9)+&
                         HRRB(19))+&
                         2.D0*HRRB(31))+&
                         HRRB(48)
      !=(9_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(18)+&
                         ABz*(-HRR(8)+&
                         HRRA(28))+&
                         ABx*(-HRR(10)+&
                         ABz*(-HRR(4)+&
                         HRRA(17))+&
                         HRRA(32))+&
                         HRRA(49)
      !=(9,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(28)+&
                         ABy*(ABz*HRRB(16)+&
                         HRRB(31))+&
                         ABx*(ABz*HRRB(17)+&
                         ABy*(ABz*HRRB(9)+&
                         HRRB(19))+&
                         HRRB(32))+&
                         HRRB(49)
      !=(9_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(16)+&
                         ABz*(-HRR(6)+&
                         HRRA(31))+&
                         ABx*(-HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRA(19))+&
                         HRRA(34))+&
                         HRRA(52)
      !=(9,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(16)+&
                         ABz*(ABz*HRRB(16)+&
                         2.D0*HRRB(31))+&
                         ABx*(-HRR(9)+&
                         ABz*(ABz*HRRB(9)+&
                         2.D0*HRRB(19))+&
                         HRRB(34))+&
                         HRRB(52)
      OffSet=(OA+5)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(10_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(30)+&
                         ABx*(ABz*HRRA(18)+&
                         HRRA(33))+&
                         HRRA(51)
      !=(10,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(20)+&
                         ABz*(-HRR(10)+&
                         HRRB(30))+&
                         ABx*(2.D0*ABz*HRRB(18)+&
                         ABx*(ABz*HRRB(10)+&
                         HRRB(20))+&
                         2.D0*HRRB(33))+&
                         HRRB(51)
      !=(10_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(31)+&
                         ABx*(ABz*HRRA(19)+&
                         HRRA(34))+&
                         HRRA(52)
      !=(10,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(31)+&
                         ABy*(ABz*HRRB(18)+&
                         HRRB(33))+&
                         ABx*(ABz*HRRB(19)+&
                         ABy*(ABz*HRRB(10)+&
                         HRRB(20))+&
                         HRRB(34))+&
                         HRRB(52)
      !=(10_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(18)+&
                         ABz*(-2.D0*HRR(8)+&
                         HRRA(33))+&
                         ABx*(-2.D0*HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         HRRA(20))+&
                         HRRA(35))+&
                         HRRA(54)
      !=(10,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(18)+&
                         ABz*(ABz*HRRB(18)+&
                         2.D0*HRRB(33))+&
                         ABx*(-HRR(10)+&
                         ABz*(ABz*HRRB(10)+&
                         2.D0*HRRB(20))+&
                         HRRB(35))+&
                         HRRB(54)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(5_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(16)+&
                         ABz*(-2.D0*HRR(6)+&
                         HRRA(22))+&
                         ABy*(-2.D0*HRR(8)+&
                         ABz*(-2.D0*HRR(2)+&
                         HRRA(11))+&
                         HRRA(26))+&
                         HRRA(43)
      !=(5,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(22)+&
                         ABy*(ABz*HRRB(11)+&
                         HRRB(26))+&
                         ABx*(ABz*HRRB(12)+&
                         ABy*(ABz*HRRB(5)+&
                         HRRB(15))+&
                         HRRB(27))+&
                         HRRB(43)
      !=(5_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(23)+&
                         ABy*(ABz*HRRA(12)+&
                         HRRA(27))+&
                         HRRA(44)
      !=(5,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(15)+&
                         ABz*(-HRR(5)+&
                         HRRB(23))+&
                         ABy*(2.D0*ABz*HRRB(12)+&
                         ABy*(ABz*HRRB(5)+&
                         HRRB(15))+&
                         2.D0*HRRB(27))+&
                         HRRB(44)
      !=(5_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(27)+&
                         ABy*(ABz*HRRA(15)+&
                         HRRA(30))+&
                         HRRA(48)
      !=(5,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(12)+&
                         ABz*(ABz*HRRB(12)+&
                         2.D0*HRRB(27))+&
                         ABy*(-HRR(5)+&
                         ABz*(ABz*HRRB(5)+&
                         2.D0*HRRB(15))+&
                         HRRB(30))+&
                         HRRB(48)
      OffSet=(OA+1)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(6_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(17)+&
                         ABz*(-HRR(7)+&
                         HRRA(23))+&
                         ABy*(-HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRA(12))+&
                         HRRA(27))+&
                         HRRA(44)
      !=(6,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(23)+&
                         ABy*(ABz*HRRB(12)+&
                         HRRB(27))+&
                         ABx*(ABz*HRRB(13)+&
                         ABy*(ABz*HRRB(6)+&
                         HRRB(16))+&
                         HRRB(28))+&
                         HRRB(44)
      !=(6_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(16)+&
                         ABz*(-HRR(6)+&
                         HRRA(24))+&
                         ABy*(-HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRA(13))+&
                         HRRA(28))+&
                         HRRA(45)
      !=(6,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(16)+&
                         ABz*(-HRR(6)+&
                         HRRB(24))+&
                         ABy*(2.D0*ABz*HRRB(13)+&
                         ABy*(ABz*HRRB(6)+&
                         HRRB(16))+&
                         2.D0*HRRB(28))+&
                         HRRB(45)
      !=(6_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(28)+&
                         ABy*(ABz*HRRA(16)+&
                         HRRA(31))+&
                         HRRA(49)
      !=(6,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(13)+&
                         ABz*(ABz*HRRB(13)+&
                         2.D0*HRRB(28))+&
                         ABy*(-HRR(6)+&
                         ABz*(ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         HRRB(31))+&
                         HRRB(49)
      OffSet=(OA+2)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(7_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(24)+&
                         ABy*(ABz*HRRA(13)+&
                         HRRA(28))+&
                         HRRA(45)
      !=(7,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(24)+&
                         ABy*(ABz*HRRB(13)+&
                         HRRB(28))+&
                         ABx*(ABz*HRRB(14)+&
                         ABy*(ABz*HRRB(7)+&
                         HRRB(17))+&
                         HRRB(29))+&
                         HRRB(45)
      !=(7_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(17)+&
                         ABz*(-2.D0*HRR(7)+&
                         HRRA(25))+&
                         ABy*(-2.D0*HRR(9)+&
                         ABz*(-2.D0*HRR(3)+&
                         HRRA(14))+&
                         HRRA(29))+&
                         HRRA(46)
      !=(7,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(17)+&
                         ABz*(-HRR(7)+&
                         HRRB(25))+&
                         ABy*(2.D0*ABz*HRRB(14)+&
                         ABy*(ABz*HRRB(7)+&
                         HRRB(17))+&
                         2.D0*HRRB(29))+&
                         HRRB(46)
      !=(7_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(29)+&
                         ABy*(ABz*HRRA(17)+&
                         HRRA(32))+&
                         HRRA(50)
      !=(7,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(14)+&
                         ABz*(ABz*HRRB(14)+&
                         2.D0*HRRB(29))+&
                         ABy*(-HRR(7)+&
                         ABz*(ABz*HRRB(7)+&
                         2.D0*HRRB(17))+&
                         HRRB(32))+&
                         HRRB(50)
      OffSet=(OA+3)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(8_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(19)+&
                         ABz*(-HRR(9)+&
                         HRRA(27))+&
                         ABy*(-HRR(10)+&
                         ABz*(-HRR(4)+&
                         HRRA(15))+&
                         HRRA(30))+&
                         HRRA(48)
      !=(8,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(27)+&
                         ABy*(ABz*HRRB(15)+&
                         HRRB(30))+&
                         ABx*(ABz*HRRB(16)+&
                         ABy*(ABz*HRRB(8)+&
                         HRRB(18))+&
                         HRRB(31))+&
                         HRRB(48)
      !=(8_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(28)+&
                         ABy*(ABz*HRRA(16)+&
                         HRRA(31))+&
                         HRRA(49)
      !=(8,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(18)+&
                         ABz*(-HRR(8)+&
                         HRRB(28))+&
                         ABy*(2.D0*ABz*HRRB(16)+&
                         ABy*(ABz*HRRB(8)+&
                         HRRB(18))+&
                         2.D0*HRRB(31))+&
                         HRRB(49)
      !=(8_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(16)+&
                         ABz*(-HRR(6)+&
                         HRRA(31))+&
                         ABy*(-HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRA(18))+&
                         HRRA(33))+&
                         HRRA(52)
      !=(8,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(16)+&
                         ABz*(ABz*HRRB(16)+&
                         2.D0*HRRB(31))+&
                         ABy*(-HRR(8)+&
                         ABz*(ABz*HRRB(8)+&
                         2.D0*HRRB(18))+&
                         HRRB(33))+&
                         HRRB(52)
      OffSet=(OA+4)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(9_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(28)+&
                         ABy*(ABz*HRRA(16)+&
                         HRRA(31))+&
                         HRRA(49)
      !=(9,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(28)+&
                         ABy*(ABz*HRRB(16)+&
                         HRRB(31))+&
                         ABx*(ABz*HRRB(17)+&
                         ABy*(ABz*HRRB(9)+&
                         HRRB(19))+&
                         HRRB(32))+&
                         HRRB(49)
      !=(9_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(19)+&
                         ABz*(-HRR(9)+&
                         HRRA(29))+&
                         ABy*(-HRR(10)+&
                         ABz*(-HRR(4)+&
                         HRRA(17))+&
                         HRRA(32))+&
                         HRRA(50)
      !=(9,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(19)+&
                         ABz*(-HRR(9)+&
                         HRRB(29))+&
                         ABy*(2.D0*ABz*HRRB(17)+&
                         ABy*(ABz*HRRB(9)+&
                         HRRB(19))+&
                         2.D0*HRRB(32))+&
                         HRRB(50)
      !=(9_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(17)+&
                         ABz*(-HRR(7)+&
                         HRRA(32))+&
                         ABy*(-HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRA(19))+&
                         HRRA(34))+&
                         HRRA(53)
      !=(9,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(17)+&
                         ABz*(ABz*HRRB(17)+&
                         2.D0*HRRB(32))+&
                         ABy*(-HRR(9)+&
                         ABz*(ABz*HRRB(9)+&
                         2.D0*HRRB(19))+&
                         HRRB(34))+&
                         HRRB(53)
      OffSet=(OA+5)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(10_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(31)+&
                         ABy*(ABz*HRRA(18)+&
                         HRRA(33))+&
                         HRRA(52)
      !=(10,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(31)+&
                         ABy*(ABz*HRRB(18)+&
                         HRRB(33))+&
                         ABx*(ABz*HRRB(19)+&
                         ABy*(ABz*HRRB(10)+&
                         HRRB(20))+&
                         HRRB(34))+&
                         HRRB(52)
      !=(10_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(32)+&
                         ABy*(ABz*HRRA(19)+&
                         HRRA(34))+&
                         HRRA(53)
      !=(10,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(20)+&
                         ABz*(-HRR(10)+&
                         HRRB(32))+&
                         ABy*(2.D0*ABz*HRRB(19)+&
                         ABy*(ABz*HRRB(10)+&
                         HRRB(20))+&
                         2.D0*HRRB(34))+&
                         HRRB(53)
      !=(10_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(19)+&
                         ABz*(-2.D0*HRR(9)+&
                         HRRA(34))+&
                         ABy*(-2.D0*HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         HRRA(20))+&
                         HRRA(35))+&
                         HRRA(55)
      !=(10,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(19)+&
                         ABz*(ABz*HRRB(19)+&
                         2.D0*HRRB(34))+&
                         ABy*(-HRR(10)+&
                         ABz*(ABz*HRRB(10)+&
                         2.D0*HRRB(20))+&
                         HRRB(35))+&
                         HRRB(55)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(5_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(18)+&
                         ABz*(-4.D0*HRR(8)+&
                         ABz*(-2.D0*HRR(2)+&
                         HRRA(11))+&
                         2.D0*HRRA(26))+&
                         HRRA(47)
      !=(5,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(11)+&
                         2.D0*HRRB(26))+&
                         ABx*(ABz*(ABz*HRRB(5)+&
                         2.D0*HRRB(15))+&
                         HRRB(30))+&
                         HRRB(47)
      !=(5_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*(ABz*HRRA(12)+&
                         2.D0*HRRA(27))+&
                         HRRA(48)
      !=(5,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(12)+&
                         2.D0*HRRB(27))+&
                         ABy*(ABz*(ABz*HRRB(5)+&
                         2.D0*HRRB(15))+&
                         HRRB(30))+&
                         HRRB(48)
      !=(5_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*(ABz*HRRA(15)+&
                         2.D0*HRRA(30))+&
                         HRRA(51)
      !=(5,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(15)+&
                         ABz*(-2.D0*HRR(5)+&
                         ABz*(ABz*HRRB(5)+&
                         3.D0*HRRB(15))+&
                         3.D0*HRRB(30))+&
                         HRRB(51)
      OffSet=(OA+1)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(6_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(19)+&
                         ABz*(-2.D0*HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRA(12))+&
                         2.D0*HRRA(27))+&
                         HRRA(48)
      !=(6,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(12)+&
                         2.D0*HRRB(27))+&
                         ABx*(ABz*(ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         HRRB(31))+&
                         HRRB(48)
      !=(6_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(18)+&
                         ABz*(-2.D0*HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRA(13))+&
                         2.D0*HRRA(28))+&
                         HRRA(49)
      !=(6,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(13)+&
                         2.D0*HRRB(28))+&
                         ABy*(ABz*(ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         HRRB(31))+&
                         HRRB(49)
      !=(6_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*(ABz*HRRA(16)+&
                         2.D0*HRRA(31))+&
                         HRRA(52)
      !=(6,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(16)+&
                         ABz*(-2.D0*HRR(6)+&
                         ABz*(ABz*HRRB(6)+&
                         3.D0*HRRB(16))+&
                         3.D0*HRRB(31))+&
                         HRRB(52)
      OffSet=(OA+2)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(7_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*(ABz*HRRA(13)+&
                         2.D0*HRRA(28))+&
                         HRRA(49)
      !=(7,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(13)+&
                         2.D0*HRRB(28))+&
                         ABx*(ABz*(ABz*HRRB(7)+&
                         2.D0*HRRB(17))+&
                         HRRB(32))+&
                         HRRB(49)
      !=(7_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(19)+&
                         ABz*(-4.D0*HRR(9)+&
                         ABz*(-2.D0*HRR(3)+&
                         HRRA(14))+&
                         2.D0*HRRA(29))+&
                         HRRA(50)
      !=(7,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(14)+&
                         2.D0*HRRB(29))+&
                         ABy*(ABz*(ABz*HRRB(7)+&
                         2.D0*HRRB(17))+&
                         HRRB(32))+&
                         HRRB(50)
      !=(7_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*(ABz*HRRA(17)+&
                         2.D0*HRRA(32))+&
                         HRRA(53)
      !=(7,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(17)+&
                         ABz*(-2.D0*HRR(7)+&
                         ABz*(ABz*HRRB(7)+&
                         3.D0*HRRB(17))+&
                         3.D0*HRRB(32))+&
                         HRRB(53)
      OffSet=(OA+3)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(8_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(20)+&
                         ABz*(-2.D0*HRR(10)+&
                         ABz*(-HRR(4)+&
                         HRRA(15))+&
                         2.D0*HRRA(30))+&
                         HRRA(51)
      !=(8,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(15)+&
                         2.D0*HRRB(30))+&
                         ABx*(ABz*(ABz*HRRB(8)+&
                         2.D0*HRRB(18))+&
                         HRRB(33))+&
                         HRRB(51)
      !=(8_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*(ABz*HRRA(16)+&
                         2.D0*HRRA(31))+&
                         HRRA(52)
      !=(8,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(16)+&
                         2.D0*HRRB(31))+&
                         ABy*(ABz*(ABz*HRRB(8)+&
                         2.D0*HRRB(18))+&
                         HRRB(33))+&
                         HRRB(52)
      !=(8_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(18)+&
                         ABz*(-2.D0*HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRA(18))+&
                         2.D0*HRRA(33))+&
                         HRRA(54)
      !=(8,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(18)+&
                         ABz*(-2.D0*HRR(8)+&
                         ABz*(ABz*HRRB(8)+&
                         3.D0*HRRB(18))+&
                         3.D0*HRRB(33))+&
                         HRRB(54)
      OffSet=(OA+4)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(9_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*(ABz*HRRA(16)+&
                         2.D0*HRRA(31))+&
                         HRRA(52)
      !=(9,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(16)+&
                         2.D0*HRRB(31))+&
                         ABx*(ABz*(ABz*HRRB(9)+&
                         2.D0*HRRB(19))+&
                         HRRB(34))+&
                         HRRB(52)
      !=(9_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(20)+&
                         ABz*(-2.D0*HRR(10)+&
                         ABz*(-HRR(4)+&
                         HRRA(17))+&
                         2.D0*HRRA(32))+&
                         HRRA(53)
      !=(9,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(17)+&
                         2.D0*HRRB(32))+&
                         ABy*(ABz*(ABz*HRRB(9)+&
                         2.D0*HRRB(19))+&
                         HRRB(34))+&
                         HRRB(53)
      !=(9_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(19)+&
                         ABz*(-2.D0*HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRA(19))+&
                         2.D0*HRRA(34))+&
                         HRRA(55)
      !=(9,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(19)+&
                         ABz*(-2.D0*HRR(9)+&
                         ABz*(ABz*HRRB(9)+&
                         3.D0*HRRB(19))+&
                         3.D0*HRRB(34))+&
                         HRRB(55)
      OffSet=(OA+5)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(10_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*(ABz*HRRA(18)+&
                         2.D0*HRRA(33))+&
                         HRRA(54)
      !=(10,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(18)+&
                         2.D0*HRRB(33))+&
                         ABx*(ABz*(ABz*HRRB(10)+&
                         2.D0*HRRB(20))+&
                         HRRB(35))+&
                         HRRB(54)
      !=(10_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*(ABz*HRRA(19)+&
                         2.D0*HRRA(34))+&
                         HRRA(55)
      !=(10,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(19)+&
                         2.D0*HRRB(34))+&
                         ABy*(ABz*(ABz*HRRB(10)+&
                         2.D0*HRRB(20))+&
                         HRRB(35))+&
                         HRRB(55)
      !=(10_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(20)+&
                         ABz*(-4.D0*HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         HRRA(20))+&
                         2.D0*HRRA(35))+&
                         HRRA(56)
      !=(10,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(20)+&
                         ABz*(-2.D0*HRR(10)+&
                         ABz*(ABz*HRRB(10)+&
                         3.D0*HRRB(20))+&
                         3.D0*HRRB(35))+&
                         HRRB(56)
    END SUBROUTINE BraHRR66ab
    SUBROUTINE BraHRR66cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,Cart,HRR,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,Cart,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(5,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(5)+&
                                2.D0*HRR(11))+&
                                HRR(21)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(6,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(6)+&
                                2.D0*HRR(12))+&
                                HRR(22)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(7,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(7)+&
                                2.D0*HRR(13))+&
                                HRR(23)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(8,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(8)+&
                                2.D0*HRR(15))+&
                                HRR(26)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(9,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(9)+&
                                2.D0*HRR(16))+&
                                HRR(27)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(10,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(10)+&
                                2.D0*HRR(18))+&
                                HRR(30)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(5,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(11)+&
                                ABx*(ABy*HRR(5)+&
                                HRR(12))+&
                                HRR(22)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(6,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(12)+&
                                ABx*(ABy*HRR(6)+&
                                HRR(13))+&
                                HRR(23)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(7,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(13)+&
                                ABx*(ABy*HRR(7)+&
                                HRR(14))+&
                                HRR(24)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(8,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(15)+&
                                ABx*(ABy*HRR(8)+&
                                HRR(16))+&
                                HRR(27)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+4)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(9,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(16)+&
                                ABx*(ABy*HRR(9)+&
                                HRR(17))+&
                                HRR(28)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+5)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(10,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(18)+&
                                ABx*(ABy*HRR(10)+&
                                HRR(19))+&
                                HRR(31)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(5,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(5)+&
                                2.D0*HRR(12))+&
                                HRR(23)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(6,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(6)+&
                                2.D0*HRR(13))+&
                                HRR(24)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(7,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(7)+&
                                2.D0*HRR(14))+&
                                HRR(25)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(8,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(8)+&
                                2.D0*HRR(16))+&
                                HRR(28)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+4)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(9,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(9)+&
                                2.D0*HRR(17))+&
                                HRR(29)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+5)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(10,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(10)+&
                                2.D0*HRR(19))+&
                                HRR(32)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(5,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(11)+&
                                ABx*(ABz*HRR(5)+&
                                HRR(15))+&
                                HRR(26)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(6,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(12)+&
                                ABx*(ABz*HRR(6)+&
                                HRR(16))+&
                                HRR(27)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(7,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(13)+&
                                ABx*(ABz*HRR(7)+&
                                HRR(17))+&
                                HRR(28)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(8,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(15)+&
                                ABx*(ABz*HRR(8)+&
                                HRR(18))+&
                                HRR(30)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+4)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(9,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(16)+&
                                ABx*(ABz*HRR(9)+&
                                HRR(19))+&
                                HRR(31)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+5)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(10,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(18)+&
                                ABx*(ABz*HRR(10)+&
                                HRR(20))+&
                                HRR(33)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(5,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(12)+&
                                ABy*(ABz*HRR(5)+&
                                HRR(15))+&
                                HRR(27)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(6,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(13)+&
                                ABy*(ABz*HRR(6)+&
                                HRR(16))+&
                                HRR(28)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(7,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(14)+&
                                ABy*(ABz*HRR(7)+&
                                HRR(17))+&
                                HRR(29)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(8,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(16)+&
                                ABy*(ABz*HRR(8)+&
                                HRR(18))+&
                                HRR(31)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+4)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(9,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(17)+&
                                ABy*(ABz*HRR(9)+&
                                HRR(19))+&
                                HRR(32)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+5)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(10,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(19)+&
                                ABy*(ABz*HRR(10)+&
                                HRR(20))+&
                                HRR(34)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(5,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(5)+&
                                2.D0*HRR(15))+&
                                HRR(30)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(6,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(6)+&
                                2.D0*HRR(16))+&
                                HRR(31)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(7,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(7)+&
                                2.D0*HRR(17))+&
                                HRR(32)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(8,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(8)+&
                                2.D0*HRR(18))+&
                                HRR(33)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+4)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(9,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(9)+&
                                2.D0*HRR(19))+&
                                HRR(34)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+5)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(10,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(10)+&
                                2.D0*HRR(20))+&
                                HRR(35)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
    END SUBROUTINE BraHRR66cd
