    SUBROUTINE BraHRR103ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR,HRRA,HRRB,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*),HRRA(*),HRRB(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(11_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-3.D0*HRR(11)+&
                         ABx*(-3.D0*HRR(5)+&
                         HRRA(21))+&
                         HRRA(36)
      !=(11,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(11)+&
                         ABx*(ABx*HRRB(11)+&
                         2.D0*HRRB(21))+&
                         HRRB(36)
      !=(11_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*HRRA(22)+&
                         HRRA(37)
      !=(11,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(21)+&
                         ABx*(ABy*HRRB(11)+&
                         HRRB(22))+&
                         HRRB(37)
      !=(11_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*HRRA(26)+&
                         HRRA(42)
      !=(11,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(21)+&
                         ABx*(ABz*HRRB(11)+&
                         HRRB(26))+&
                         HRRB(42)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(12_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(12)+&
                         ABx*(-2.D0*HRR(6)+&
                         HRRA(22))+&
                         HRRA(37)
      !=(12,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(12)+&
                         ABx*(ABx*HRRB(12)+&
                         2.D0*HRRB(22))+&
                         HRRB(37)
      !=(12_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(11)+&
                         ABx*(-HRR(5)+&
                         HRRA(23))+&
                         HRRA(38)
      !=(12,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(22)+&
                         ABx*(ABy*HRRB(12)+&
                         HRRB(23))+&
                         HRRB(38)
      !=(12_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*HRRA(27)+&
                         HRRA(43)
      !=(12,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(22)+&
                         ABx*(ABz*HRRB(12)+&
                         HRRB(27))+&
                         HRRB(43)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(13_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(13)+&
                         ABx*(-HRR(7)+&
                         HRRA(23))+&
                         HRRA(38)
      !=(13,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(13)+&
                         ABx*(ABx*HRRB(13)+&
                         2.D0*HRRB(23))+&
                         HRRB(38)
      !=(13_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(12)+&
                         ABx*(-2.D0*HRR(6)+&
                         HRRA(24))+&
                         HRRA(39)
      !=(13,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(23)+&
                         ABx*(ABy*HRRB(13)+&
                         HRRB(24))+&
                         HRRB(39)
      !=(13_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*HRRA(28)+&
                         HRRA(44)
      !=(13,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(23)+&
                         ABx*(ABz*HRRB(13)+&
                         HRRB(28))+&
                         HRRB(44)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(14_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*HRRA(24)+&
                         HRRA(39)
      !=(14,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(14)+&
                         ABx*(ABx*HRRB(14)+&
                         2.D0*HRRB(24))+&
                         HRRB(39)
      !=(14_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-3.D0*HRR(13)+&
                         ABx*(-3.D0*HRR(7)+&
                         HRRA(25))+&
                         HRRA(40)
      !=(14,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(24)+&
                         ABx*(ABy*HRRB(14)+&
                         HRRB(25))+&
                         HRRB(40)
      !=(14_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*HRRA(29)+&
                         HRRA(45)
      !=(14,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(24)+&
                         ABx*(ABz*HRRB(14)+&
                         HRRB(29))+&
                         HRRB(45)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(15_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(15)+&
                         ABx*(-2.D0*HRR(8)+&
                         HRRA(26))+&
                         HRRA(42)
      !=(15,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(15)+&
                         ABx*(ABx*HRRB(15)+&
                         2.D0*HRRB(26))+&
                         HRRB(42)
      !=(15_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*HRRA(27)+&
                         HRRA(43)
      !=(15,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(26)+&
                         ABx*(ABy*HRRB(15)+&
                         HRRB(27))+&
                         HRRB(43)
      !=(15_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(11)+&
                         ABx*(-HRR(5)+&
                         HRRA(30))+&
                         HRRA(47)
      !=(15,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(26)+&
                         ABx*(ABz*HRRB(15)+&
                         HRRB(30))+&
                         HRRB(47)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(16_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(16)+&
                         ABx*(-HRR(9)+&
                         HRRA(27))+&
                         HRRA(43)
      !=(16,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(16)+&
                         ABx*(ABx*HRRB(16)+&
                         2.D0*HRRB(27))+&
                         HRRB(43)
      !=(16_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(15)+&
                         ABx*(-HRR(8)+&
                         HRRA(28))+&
                         HRRA(44)
      !=(16,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(27)+&
                         ABx*(ABy*HRRB(16)+&
                         HRRB(28))+&
                         HRRB(44)
      !=(16_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(12)+&
                         ABx*(-HRR(6)+&
                         HRRA(31))+&
                         HRRA(48)
      !=(16,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(27)+&
                         ABx*(ABz*HRRB(16)+&
                         HRRB(31))+&
                         HRRB(48)
      OffSet=(OA+6)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(17_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*HRRA(28)+&
                         HRRA(44)
      !=(17,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(17)+&
                         ABx*(ABx*HRRB(17)+&
                         2.D0*HRRB(28))+&
                         HRRB(44)
      !=(17_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(16)+&
                         ABx*(-2.D0*HRR(9)+&
                         HRRA(29))+&
                         HRRA(45)
      !=(17,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(28)+&
                         ABx*(ABy*HRRB(17)+&
                         HRRB(29))+&
                         HRRB(45)
      !=(17_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(13)+&
                         ABx*(-HRR(7)+&
                         HRRA(32))+&
                         HRRA(49)
      !=(17,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(28)+&
                         ABx*(ABz*HRRB(17)+&
                         HRRB(32))+&
                         HRRB(49)
      OffSet=(OA+7)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(18_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(18)+&
                         ABx*(-HRR(10)+&
                         HRRA(30))+&
                         HRRA(47)
      !=(18,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(18)+&
                         ABx*(ABx*HRRB(18)+&
                         2.D0*HRRB(30))+&
                         HRRB(47)
      !=(18_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*HRRA(31)+&
                         HRRA(48)
      !=(18,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(30)+&
                         ABx*(ABy*HRRB(18)+&
                         HRRB(31))+&
                         HRRB(48)
      !=(18_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(15)+&
                         ABx*(-2.D0*HRR(8)+&
                         HRRA(33))+&
                         HRRA(51)
      !=(18,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(30)+&
                         ABx*(ABz*HRRB(18)+&
                         HRRB(33))+&
                         HRRB(51)
      OffSet=(OA+8)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(19_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*HRRA(31)+&
                         HRRA(48)
      !=(19,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(19)+&
                         ABx*(ABx*HRRB(19)+&
                         2.D0*HRRB(31))+&
                         HRRB(48)
      !=(19_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(18)+&
                         ABx*(-HRR(10)+&
                         HRRA(32))+&
                         HRRA(49)
      !=(19,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(31)+&
                         ABx*(ABy*HRRB(19)+&
                         HRRB(32))+&
                         HRRB(49)
      !=(19_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(16)+&
                         ABx*(-2.D0*HRR(9)+&
                         HRRA(34))+&
                         HRRA(52)
      !=(19,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(31)+&
                         ABx*(ABz*HRRB(19)+&
                         HRRB(34))+&
                         HRRB(52)
      OffSet=(OA+9)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(20_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*HRRA(33)+&
                         HRRA(51)
      !=(20,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(20)+&
                         ABx*(ABx*HRRB(20)+&
                         2.D0*HRRB(33))+&
                         HRRB(51)
      !=(20_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*HRRA(34)+&
                         HRRA(52)
      !=(20,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(33)+&
                         ABx*(ABy*HRRB(20)+&
                         HRRB(34))+&
                         HRRB(52)
      !=(20_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-3.D0*HRR(18)+&
                         ABx*(-3.D0*HRR(10)+&
                         HRRA(35))+&
                         HRRA(54)
      !=(20,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(33)+&
                         ABx*(ABz*HRRB(20)+&
                         HRRB(35))+&
                         HRRB(54)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(11_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-3.D0*HRR(12)+&
                         ABy*(-3.D0*HRR(5)+&
                         HRRA(21))+&
                         HRRA(37)
      !=(11,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(21)+&
                         ABx*(ABy*HRRB(11)+&
                         HRRB(22))+&
                         HRRB(37)
      !=(11_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(22)+&
                         HRRA(38)
      !=(11,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(11)+&
                         ABy*(ABy*HRRB(11)+&
                         2.D0*HRRB(22))+&
                         HRRB(38)
      !=(11_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(26)+&
                         HRRA(43)
      !=(11,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(22)+&
                         ABy*(ABz*HRRB(11)+&
                         HRRB(26))+&
                         HRRB(43)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(12_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(13)+&
                         ABy*(-2.D0*HRR(6)+&
                         HRRA(22))+&
                         HRRA(38)
      !=(12,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(22)+&
                         ABx*(ABy*HRRB(12)+&
                         HRRB(23))+&
                         HRRB(38)
      !=(12_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(12)+&
                         ABy*(-HRR(5)+&
                         HRRA(23))+&
                         HRRA(39)
      !=(12,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(12)+&
                         ABy*(ABy*HRRB(12)+&
                         2.D0*HRRB(23))+&
                         HRRB(39)
      !=(12_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(27)+&
                         HRRA(44)
      !=(12,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(23)+&
                         ABy*(ABz*HRRB(12)+&
                         HRRB(27))+&
                         HRRB(44)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(13_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(14)+&
                         ABy*(-HRR(7)+&
                         HRRA(23))+&
                         HRRA(39)
      !=(13,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(23)+&
                         ABx*(ABy*HRRB(13)+&
                         HRRB(24))+&
                         HRRB(39)
      !=(13_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(13)+&
                         ABy*(-2.D0*HRR(6)+&
                         HRRA(24))+&
                         HRRA(40)
      !=(13,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(13)+&
                         ABy*(ABy*HRRB(13)+&
                         2.D0*HRRB(24))+&
                         HRRB(40)
      !=(13_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(28)+&
                         HRRA(45)
      !=(13,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(24)+&
                         ABy*(ABz*HRRB(13)+&
                         HRRB(28))+&
                         HRRB(45)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(14_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(24)+&
                         HRRA(40)
      !=(14,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(24)+&
                         ABx*(ABy*HRRB(14)+&
                         HRRB(25))+&
                         HRRB(40)
      !=(14_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-3.D0*HRR(14)+&
                         ABy*(-3.D0*HRR(7)+&
                         HRRA(25))+&
                         HRRA(41)
      !=(14,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(14)+&
                         ABy*(ABy*HRRB(14)+&
                         2.D0*HRRB(25))+&
                         HRRB(41)
      !=(14_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(29)+&
                         HRRA(46)
      !=(14,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(25)+&
                         ABy*(ABz*HRRB(14)+&
                         HRRB(29))+&
                         HRRB(46)
      OffSet=(OA+4)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(15_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(16)+&
                         ABy*(-2.D0*HRR(8)+&
                         HRRA(26))+&
                         HRRA(43)
      !=(15,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(26)+&
                         ABx*(ABy*HRRB(15)+&
                         HRRB(27))+&
                         HRRB(43)
      !=(15_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(27)+&
                         HRRA(44)
      !=(15,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(15)+&
                         ABy*(ABy*HRRB(15)+&
                         2.D0*HRRB(27))+&
                         HRRB(44)
      !=(15_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(12)+&
                         ABy*(-HRR(5)+&
                         HRRA(30))+&
                         HRRA(48)
      !=(15,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(27)+&
                         ABy*(ABz*HRRB(15)+&
                         HRRB(30))+&
                         HRRB(48)
      OffSet=(OA+5)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(16_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(17)+&
                         ABy*(-HRR(9)+&
                         HRRA(27))+&
                         HRRA(44)
      !=(16,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(27)+&
                         ABx*(ABy*HRRB(16)+&
                         HRRB(28))+&
                         HRRB(44)
      !=(16_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(16)+&
                         ABy*(-HRR(8)+&
                         HRRA(28))+&
                         HRRA(45)
      !=(16,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(16)+&
                         ABy*(ABy*HRRB(16)+&
                         2.D0*HRRB(28))+&
                         HRRB(45)
      !=(16_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(13)+&
                         ABy*(-HRR(6)+&
                         HRRA(31))+&
                         HRRA(49)
      !=(16,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(28)+&
                         ABy*(ABz*HRRB(16)+&
                         HRRB(31))+&
                         HRRB(49)
      OffSet=(OA+6)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(17_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(28)+&
                         HRRA(45)
      !=(17,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(28)+&
                         ABx*(ABy*HRRB(17)+&
                         HRRB(29))+&
                         HRRB(45)
      !=(17_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(17)+&
                         ABy*(-2.D0*HRR(9)+&
                         HRRA(29))+&
                         HRRA(46)
      !=(17,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(17)+&
                         ABy*(ABy*HRRB(17)+&
                         2.D0*HRRB(29))+&
                         HRRB(46)
      !=(17_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(14)+&
                         ABy*(-HRR(7)+&
                         HRRA(32))+&
                         HRRA(50)
      !=(17,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(29)+&
                         ABy*(ABz*HRRB(17)+&
                         HRRB(32))+&
                         HRRB(50)
      OffSet=(OA+7)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(18_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(19)+&
                         ABy*(-HRR(10)+&
                         HRRA(30))+&
                         HRRA(48)
      !=(18,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(30)+&
                         ABx*(ABy*HRRB(18)+&
                         HRRB(31))+&
                         HRRB(48)
      !=(18_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(31)+&
                         HRRA(49)
      !=(18,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(18)+&
                         ABy*(ABy*HRRB(18)+&
                         2.D0*HRRB(31))+&
                         HRRB(49)
      !=(18_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(16)+&
                         ABy*(-2.D0*HRR(8)+&
                         HRRA(33))+&
                         HRRA(52)
      !=(18,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(31)+&
                         ABy*(ABz*HRRB(18)+&
                         HRRB(33))+&
                         HRRB(52)
      OffSet=(OA+8)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(19_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(31)+&
                         HRRA(49)
      !=(19,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(31)+&
                         ABx*(ABy*HRRB(19)+&
                         HRRB(32))+&
                         HRRB(49)
      !=(19_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(19)+&
                         ABy*(-HRR(10)+&
                         HRRA(32))+&
                         HRRA(50)
      !=(19,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(19)+&
                         ABy*(ABy*HRRB(19)+&
                         2.D0*HRRB(32))+&
                         HRRB(50)
      !=(19_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(17)+&
                         ABy*(-2.D0*HRR(9)+&
                         HRRA(34))+&
                         HRRA(53)
      !=(19,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(32)+&
                         ABy*(ABz*HRRB(19)+&
                         HRRB(34))+&
                         HRRB(53)
      OffSet=(OA+9)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(20_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(33)+&
                         HRRA(52)
      !=(20,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(33)+&
                         ABx*(ABy*HRRB(20)+&
                         HRRB(34))+&
                         HRRB(52)
      !=(20_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(34)+&
                         HRRA(53)
      !=(20,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(20)+&
                         ABy*(ABy*HRRB(20)+&
                         2.D0*HRRB(34))+&
                         HRRB(53)
      !=(20_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-3.D0*HRR(19)+&
                         ABy*(-3.D0*HRR(10)+&
                         HRRA(35))+&
                         HRRA(55)
      !=(20,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(34)+&
                         ABy*(ABz*HRRB(20)+&
                         HRRB(35))+&
                         HRRB(55)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(11_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-3.D0*HRR(15)+&
                         ABz*(-3.D0*HRR(5)+&
                         HRRA(21))+&
                         HRRA(42)
      !=(11,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(21)+&
                         ABx*(ABz*HRRB(11)+&
                         HRRB(26))+&
                         HRRB(42)
      !=(11_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(22)+&
                         HRRA(43)
      !=(11,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(22)+&
                         ABy*(ABz*HRRB(11)+&
                         HRRB(26))+&
                         HRRB(43)
      !=(11_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(26)+&
                         HRRA(47)
      !=(11,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(11)+&
                         ABz*(ABz*HRRB(11)+&
                         2.D0*HRRB(26))+&
                         HRRB(47)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(12_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(16)+&
                         ABz*(-2.D0*HRR(6)+&
                         HRRA(22))+&
                         HRRA(43)
      !=(12,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(22)+&
                         ABx*(ABz*HRRB(12)+&
                         HRRB(27))+&
                         HRRB(43)
      !=(12_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(15)+&
                         ABz*(-HRR(5)+&
                         HRRA(23))+&
                         HRRA(44)
      !=(12,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(23)+&
                         ABy*(ABz*HRRB(12)+&
                         HRRB(27))+&
                         HRRB(44)
      !=(12_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(27)+&
                         HRRA(48)
      !=(12,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(12)+&
                         ABz*(ABz*HRRB(12)+&
                         2.D0*HRRB(27))+&
                         HRRB(48)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(13_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(17)+&
                         ABz*(-HRR(7)+&
                         HRRA(23))+&
                         HRRA(44)
      !=(13,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(23)+&
                         ABx*(ABz*HRRB(13)+&
                         HRRB(28))+&
                         HRRB(44)
      !=(13_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(16)+&
                         ABz*(-2.D0*HRR(6)+&
                         HRRA(24))+&
                         HRRA(45)
      !=(13,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(24)+&
                         ABy*(ABz*HRRB(13)+&
                         HRRB(28))+&
                         HRRB(45)
      !=(13_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(28)+&
                         HRRA(49)
      !=(13,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(13)+&
                         ABz*(ABz*HRRB(13)+&
                         2.D0*HRRB(28))+&
                         HRRB(49)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(14_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(24)+&
                         HRRA(45)
      !=(14,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(24)+&
                         ABx*(ABz*HRRB(14)+&
                         HRRB(29))+&
                         HRRB(45)
      !=(14_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-3.D0*HRR(17)+&
                         ABz*(-3.D0*HRR(7)+&
                         HRRA(25))+&
                         HRRA(46)
      !=(14,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(25)+&
                         ABy*(ABz*HRRB(14)+&
                         HRRB(29))+&
                         HRRB(46)
      !=(14_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(29)+&
                         HRRA(50)
      !=(14,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(14)+&
                         ABz*(ABz*HRRB(14)+&
                         2.D0*HRRB(29))+&
                         HRRB(50)
      OffSet=(OA+4)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(15_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(18)+&
                         ABz*(-2.D0*HRR(8)+&
                         HRRA(26))+&
                         HRRA(47)
      !=(15,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(26)+&
                         ABx*(ABz*HRRB(15)+&
                         HRRB(30))+&
                         HRRB(47)
      !=(15_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(27)+&
                         HRRA(48)
      !=(15,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(27)+&
                         ABy*(ABz*HRRB(15)+&
                         HRRB(30))+&
                         HRRB(48)
      !=(15_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(15)+&
                         ABz*(-HRR(5)+&
                         HRRA(30))+&
                         HRRA(51)
      !=(15,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(15)+&
                         ABz*(ABz*HRRB(15)+&
                         2.D0*HRRB(30))+&
                         HRRB(51)
      OffSet=(OA+5)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(16_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(19)+&
                         ABz*(-HRR(9)+&
                         HRRA(27))+&
                         HRRA(48)
      !=(16,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(27)+&
                         ABx*(ABz*HRRB(16)+&
                         HRRB(31))+&
                         HRRB(48)
      !=(16_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(18)+&
                         ABz*(-HRR(8)+&
                         HRRA(28))+&
                         HRRA(49)
      !=(16,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(28)+&
                         ABy*(ABz*HRRB(16)+&
                         HRRB(31))+&
                         HRRB(49)
      !=(16_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(16)+&
                         ABz*(-HRR(6)+&
                         HRRA(31))+&
                         HRRA(52)
      !=(16,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(16)+&
                         ABz*(ABz*HRRB(16)+&
                         2.D0*HRRB(31))+&
                         HRRB(52)
      OffSet=(OA+6)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(17_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(28)+&
                         HRRA(49)
      !=(17,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(28)+&
                         ABx*(ABz*HRRB(17)+&
                         HRRB(32))+&
                         HRRB(49)
      !=(17_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(19)+&
                         ABz*(-2.D0*HRR(9)+&
                         HRRA(29))+&
                         HRRA(50)
      !=(17,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(29)+&
                         ABy*(ABz*HRRB(17)+&
                         HRRB(32))+&
                         HRRB(50)
      !=(17_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(17)+&
                         ABz*(-HRR(7)+&
                         HRRA(32))+&
                         HRRA(53)
      !=(17,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(17)+&
                         ABz*(ABz*HRRB(17)+&
                         2.D0*HRRB(32))+&
                         HRRB(53)
      OffSet=(OA+7)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(18_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(20)+&
                         ABz*(-HRR(10)+&
                         HRRA(30))+&
                         HRRA(51)
      !=(18,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(30)+&
                         ABx*(ABz*HRRB(18)+&
                         HRRB(33))+&
                         HRRB(51)
      !=(18_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(31)+&
                         HRRA(52)
      !=(18,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(31)+&
                         ABy*(ABz*HRRB(18)+&
                         HRRB(33))+&
                         HRRB(52)
      !=(18_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(18)+&
                         ABz*(-2.D0*HRR(8)+&
                         HRRA(33))+&
                         HRRA(54)
      !=(18,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(18)+&
                         ABz*(ABz*HRRB(18)+&
                         2.D0*HRRB(33))+&
                         HRRB(54)
      OffSet=(OA+8)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(19_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(31)+&
                         HRRA(52)
      !=(19,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(31)+&
                         ABx*(ABz*HRRB(19)+&
                         HRRB(34))+&
                         HRRB(52)
      !=(19_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(20)+&
                         ABz*(-HRR(10)+&
                         HRRA(32))+&
                         HRRA(53)
      !=(19,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(32)+&
                         ABy*(ABz*HRRB(19)+&
                         HRRB(34))+&
                         HRRB(53)
      !=(19_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(19)+&
                         ABz*(-2.D0*HRR(9)+&
                         HRRA(34))+&
                         HRRA(55)
      !=(19,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(19)+&
                         ABz*(ABz*HRRB(19)+&
                         2.D0*HRRB(34))+&
                         HRRB(55)
      OffSet=(OA+9)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(20_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(33)+&
                         HRRA(54)
      !=(20,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(33)+&
                         ABx*(ABz*HRRB(20)+&
                         HRRB(35))+&
                         HRRB(54)
      !=(20_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(34)+&
                         HRRA(55)
      !=(20,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(34)+&
                         ABy*(ABz*HRRB(20)+&
                         HRRB(35))+&
                         HRRB(55)
      !=(20_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-3.D0*HRR(20)+&
                         ABz*(-3.D0*HRR(10)+&
                         HRRA(35))+&
                         HRRA(56)
      !=(20,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(20)+&
                         ABz*(ABz*HRRB(20)+&
                         2.D0*HRRB(35))+&
                         HRRB(56)
    END SUBROUTINE BraHRR103ab
    SUBROUTINE BraHRR103cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,Cart,HRR,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,Cart,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(11,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(11)+&
                                HRR(21)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(12,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(12)+&
                                HRR(22)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(13,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(13)+&
                                HRR(23)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(14,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(14)+&
                                HRR(24)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(15,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(15)+&
                                HRR(26)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(16,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(16)+&
                                HRR(27)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+6)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(17,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(17)+&
                                HRR(28)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+7)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(18,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(18)+&
                                HRR(30)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+8)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(19,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(19)+&
                                HRR(31)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+9)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(20,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(20)+&
                                HRR(33)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(11,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(11)+&
                                HRR(22)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(12,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(12)+&
                                HRR(23)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(13,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(13)+&
                                HRR(24)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(14,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(14)+&
                                HRR(25)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+4)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(15,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(15)+&
                                HRR(27)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+5)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(16,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(16)+&
                                HRR(28)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+6)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(17,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(17)+&
                                HRR(29)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+7)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(18,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(18)+&
                                HRR(31)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+8)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(19,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(19)+&
                                HRR(32)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+9)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(20,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(20)+&
                                HRR(34)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(11,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(11)+&
                                HRR(26)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(12,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(12)+&
                                HRR(27)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(13,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(13)+&
                                HRR(28)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(14,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(14)+&
                                HRR(29)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+4)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(15,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(15)+&
                                HRR(30)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+5)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(16,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(16)+&
                                HRR(31)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+6)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(17,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(17)+&
                                HRR(32)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+7)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(18,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(18)+&
                                HRR(33)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+8)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(19,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(19)+&
                                HRR(34)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+9)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(20,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(20)+&
                                HRR(35)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
    END SUBROUTINE BraHRR103cd
