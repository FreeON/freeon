    SUBROUTINE BraHRR106ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR,HRRA,HRRB,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*),HRRA(*),HRRB(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(11_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-3.D0*HRR(21)+&
                         ABx*(-6.D0*HRR(11)+&
                         ABx*(-3.D0*HRR(5)+&
                         HRRA(21))+&
                         2.D0*HRRA(36))+&
                         HRRA(57)
      !=(11,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(21)+&
                         ABx*(-2.D0*HRR(11)+&
                         ABx*(ABx*HRRB(11)+&
                         3.D0*HRRB(21))+&
                         3.D0*HRRB(36))+&
                         HRRB(57)
      !=(11_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*(ABx*HRRA(22)+&
                         2.D0*HRRA(37))+&
                         HRRA(58)
      !=(11,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(36)+&
                         ABx*(2.D0*ABy*HRRB(21)+&
                         ABx*(ABy*HRRB(11)+&
                         HRRB(22))+&
                         2.D0*HRRB(37))+&
                         HRRB(58)
      !=(11_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*(ABx*HRRA(26)+&
                         2.D0*HRRA(42))+&
                         HRRA(64)
      !=(11,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(36)+&
                         ABx*(2.D0*ABz*HRRB(21)+&
                         ABx*(ABz*HRRB(11)+&
                         HRRB(26))+&
                         2.D0*HRRB(42))+&
                         HRRB(64)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(12_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(22)+&
                         ABx*(-4.D0*HRR(12)+&
                         ABx*(-2.D0*HRR(6)+&
                         HRRA(22))+&
                         2.D0*HRRA(37))+&
                         HRRA(58)
      !=(12,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(22)+&
                         ABx*(-2.D0*HRR(12)+&
                         ABx*(ABx*HRRB(12)+&
                         3.D0*HRRB(22))+&
                         3.D0*HRRB(37))+&
                         HRRB(58)
      !=(12_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(21)+&
                         ABx*(-2.D0*HRR(11)+&
                         ABx*(-HRR(5)+&
                         HRRA(23))+&
                         2.D0*HRRA(38))+&
                         HRRA(59)
      !=(12,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(37)+&
                         ABx*(2.D0*ABy*HRRB(22)+&
                         ABx*(ABy*HRRB(12)+&
                         HRRB(23))+&
                         2.D0*HRRB(38))+&
                         HRRB(59)
      !=(12_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*(ABx*HRRA(27)+&
                         2.D0*HRRA(43))+&
                         HRRA(65)
      !=(12,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(37)+&
                         ABx*(2.D0*ABz*HRRB(22)+&
                         ABx*(ABz*HRRB(12)+&
                         HRRB(27))+&
                         2.D0*HRRB(43))+&
                         HRRB(65)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(13_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(23)+&
                         ABx*(-2.D0*HRR(13)+&
                         ABx*(-HRR(7)+&
                         HRRA(23))+&
                         2.D0*HRRA(38))+&
                         HRRA(59)
      !=(13,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(23)+&
                         ABx*(-2.D0*HRR(13)+&
                         ABx*(ABx*HRRB(13)+&
                         3.D0*HRRB(23))+&
                         3.D0*HRRB(38))+&
                         HRRB(59)
      !=(13_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(22)+&
                         ABx*(-4.D0*HRR(12)+&
                         ABx*(-2.D0*HRR(6)+&
                         HRRA(24))+&
                         2.D0*HRRA(39))+&
                         HRRA(60)
      !=(13,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(38)+&
                         ABx*(2.D0*ABy*HRRB(23)+&
                         ABx*(ABy*HRRB(13)+&
                         HRRB(24))+&
                         2.D0*HRRB(39))+&
                         HRRB(60)
      !=(13_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*(ABx*HRRA(28)+&
                         2.D0*HRRA(44))+&
                         HRRA(66)
      !=(13,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(38)+&
                         ABx*(2.D0*ABz*HRRB(23)+&
                         ABx*(ABz*HRRB(13)+&
                         HRRB(28))+&
                         2.D0*HRRB(44))+&
                         HRRB(66)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(14_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*(ABx*HRRA(24)+&
                         2.D0*HRRA(39))+&
                         HRRA(60)
      !=(14,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(24)+&
                         ABx*(-2.D0*HRR(14)+&
                         ABx*(ABx*HRRB(14)+&
                         3.D0*HRRB(24))+&
                         3.D0*HRRB(39))+&
                         HRRB(60)
      !=(14_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-3.D0*HRR(23)+&
                         ABx*(-6.D0*HRR(13)+&
                         ABx*(-3.D0*HRR(7)+&
                         HRRA(25))+&
                         2.D0*HRRA(40))+&
                         HRRA(61)
      !=(14,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(39)+&
                         ABx*(2.D0*ABy*HRRB(24)+&
                         ABx*(ABy*HRRB(14)+&
                         HRRB(25))+&
                         2.D0*HRRB(40))+&
                         HRRB(61)
      !=(14_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*(ABx*HRRA(29)+&
                         2.D0*HRRA(45))+&
                         HRRA(67)
      !=(14,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(39)+&
                         ABx*(2.D0*ABz*HRRB(24)+&
                         ABx*(ABz*HRRB(14)+&
                         HRRB(29))+&
                         2.D0*HRRB(45))+&
                         HRRB(67)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(15_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(26)+&
                         ABx*(-4.D0*HRR(15)+&
                         ABx*(-2.D0*HRR(8)+&
                         HRRA(26))+&
                         2.D0*HRRA(42))+&
                         HRRA(64)
      !=(15,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(26)+&
                         ABx*(-2.D0*HRR(15)+&
                         ABx*(ABx*HRRB(15)+&
                         3.D0*HRRB(26))+&
                         3.D0*HRRB(42))+&
                         HRRB(64)
      !=(15_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*(ABx*HRRA(27)+&
                         2.D0*HRRA(43))+&
                         HRRA(65)
      !=(15,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(42)+&
                         ABx*(2.D0*ABy*HRRB(26)+&
                         ABx*(ABy*HRRB(15)+&
                         HRRB(27))+&
                         2.D0*HRRB(43))+&
                         HRRB(65)
      !=(15_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(21)+&
                         ABx*(-2.D0*HRR(11)+&
                         ABx*(-HRR(5)+&
                         HRRA(30))+&
                         2.D0*HRRA(47))+&
                         HRRA(70)
      !=(15,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(42)+&
                         ABx*(2.D0*ABz*HRRB(26)+&
                         ABx*(ABz*HRRB(15)+&
                         HRRB(30))+&
                         2.D0*HRRB(47))+&
                         HRRB(70)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(16_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(27)+&
                         ABx*(-2.D0*HRR(16)+&
                         ABx*(-HRR(9)+&
                         HRRA(27))+&
                         2.D0*HRRA(43))+&
                         HRRA(65)
      !=(16,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(27)+&
                         ABx*(-2.D0*HRR(16)+&
                         ABx*(ABx*HRRB(16)+&
                         3.D0*HRRB(27))+&
                         3.D0*HRRB(43))+&
                         HRRB(65)
      !=(16_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(26)+&
                         ABx*(-2.D0*HRR(15)+&
                         ABx*(-HRR(8)+&
                         HRRA(28))+&
                         2.D0*HRRA(44))+&
                         HRRA(66)
      !=(16,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(43)+&
                         ABx*(2.D0*ABy*HRRB(27)+&
                         ABx*(ABy*HRRB(16)+&
                         HRRB(28))+&
                         2.D0*HRRB(44))+&
                         HRRB(66)
      !=(16_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(22)+&
                         ABx*(-2.D0*HRR(12)+&
                         ABx*(-HRR(6)+&
                         HRRA(31))+&
                         2.D0*HRRA(48))+&
                         HRRA(71)
      !=(16,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(43)+&
                         ABx*(2.D0*ABz*HRRB(27)+&
                         ABx*(ABz*HRRB(16)+&
                         HRRB(31))+&
                         2.D0*HRRB(48))+&
                         HRRB(71)
      OffSet=(OA+6)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(17_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*(ABx*HRRA(28)+&
                         2.D0*HRRA(44))+&
                         HRRA(66)
      !=(17,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(28)+&
                         ABx*(-2.D0*HRR(17)+&
                         ABx*(ABx*HRRB(17)+&
                         3.D0*HRRB(28))+&
                         3.D0*HRRB(44))+&
                         HRRB(66)
      !=(17_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(27)+&
                         ABx*(-4.D0*HRR(16)+&
                         ABx*(-2.D0*HRR(9)+&
                         HRRA(29))+&
                         2.D0*HRRA(45))+&
                         HRRA(67)
      !=(17,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(44)+&
                         ABx*(2.D0*ABy*HRRB(28)+&
                         ABx*(ABy*HRRB(17)+&
                         HRRB(29))+&
                         2.D0*HRRB(45))+&
                         HRRB(67)
      !=(17_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(23)+&
                         ABx*(-2.D0*HRR(13)+&
                         ABx*(-HRR(7)+&
                         HRRA(32))+&
                         2.D0*HRRA(49))+&
                         HRRA(72)
      !=(17,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(44)+&
                         ABx*(2.D0*ABz*HRRB(28)+&
                         ABx*(ABz*HRRB(17)+&
                         HRRB(32))+&
                         2.D0*HRRB(49))+&
                         HRRB(72)
      OffSet=(OA+7)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(18_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(30)+&
                         ABx*(-2.D0*HRR(18)+&
                         ABx*(-HRR(10)+&
                         HRRA(30))+&
                         2.D0*HRRA(47))+&
                         HRRA(70)
      !=(18,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(30)+&
                         ABx*(-2.D0*HRR(18)+&
                         ABx*(ABx*HRRB(18)+&
                         3.D0*HRRB(30))+&
                         3.D0*HRRB(47))+&
                         HRRB(70)
      !=(18_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*(ABx*HRRA(31)+&
                         2.D0*HRRA(48))+&
                         HRRA(71)
      !=(18,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(47)+&
                         ABx*(2.D0*ABy*HRRB(30)+&
                         ABx*(ABy*HRRB(18)+&
                         HRRB(31))+&
                         2.D0*HRRB(48))+&
                         HRRB(71)
      !=(18_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(26)+&
                         ABx*(-4.D0*HRR(15)+&
                         ABx*(-2.D0*HRR(8)+&
                         HRRA(33))+&
                         2.D0*HRRA(51))+&
                         HRRA(75)
      !=(18,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(47)+&
                         ABx*(2.D0*ABz*HRRB(30)+&
                         ABx*(ABz*HRRB(18)+&
                         HRRB(33))+&
                         2.D0*HRRB(51))+&
                         HRRB(75)
      OffSet=(OA+8)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(19_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*(ABx*HRRA(31)+&
                         2.D0*HRRA(48))+&
                         HRRA(71)
      !=(19,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(31)+&
                         ABx*(-2.D0*HRR(19)+&
                         ABx*(ABx*HRRB(19)+&
                         3.D0*HRRB(31))+&
                         3.D0*HRRB(48))+&
                         HRRB(71)
      !=(19_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(30)+&
                         ABx*(-2.D0*HRR(18)+&
                         ABx*(-HRR(10)+&
                         HRRA(32))+&
                         2.D0*HRRA(49))+&
                         HRRA(72)
      !=(19,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(48)+&
                         ABx*(2.D0*ABy*HRRB(31)+&
                         ABx*(ABy*HRRB(19)+&
                         HRRB(32))+&
                         2.D0*HRRB(49))+&
                         HRRB(72)
      !=(19_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(27)+&
                         ABx*(-4.D0*HRR(16)+&
                         ABx*(-2.D0*HRR(9)+&
                         HRRA(34))+&
                         2.D0*HRRA(52))+&
                         HRRA(76)
      !=(19,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(48)+&
                         ABx*(2.D0*ABz*HRRB(31)+&
                         ABx*(ABz*HRRB(19)+&
                         HRRB(34))+&
                         2.D0*HRRB(52))+&
                         HRRB(76)
      OffSet=(OA+9)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(20_x,5|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*(ABx*HRRA(33)+&
                         2.D0*HRRA(51))+&
                         HRRA(75)
      !=(20,5_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-2.D0*HRR(33)+&
                         ABx*(-2.D0*HRR(20)+&
                         ABx*(ABx*HRRB(20)+&
                         3.D0*HRRB(33))+&
                         3.D0*HRRB(51))+&
                         HRRB(75)
      !=(20_y,5|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*(ABx*HRRA(34)+&
                         2.D0*HRRA(52))+&
                         HRRA(76)
      !=(20,5_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(51)+&
                         ABx*(2.D0*ABy*HRRB(33)+&
                         ABx*(ABy*HRRB(20)+&
                         HRRB(34))+&
                         2.D0*HRRB(52))+&
                         HRRB(76)
      !=(20_z,5|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-3.D0*HRR(30)+&
                         ABx*(-6.D0*HRR(18)+&
                         ABx*(-3.D0*HRR(10)+&
                         HRRA(35))+&
                         2.D0*HRRA(54))+&
                         HRRA(79)
      !=(20,5_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(51)+&
                         ABx*(2.D0*ABz*HRRB(33)+&
                         ABx*(ABz*HRRB(20)+&
                         HRRB(35))+&
                         2.D0*HRRB(54))+&
                         HRRB(79)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(11_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-3.D0*HRR(22)+&
                         ABy*(-3.D0*HRR(11)+&
                         HRRA(36))+&
                         ABx*(-3.D0*HRR(12)+&
                         ABy*(-3.D0*HRR(5)+&
                         HRRA(21))+&
                         HRRA(37))+&
                         HRRA(58)
      !=(11,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(22)+&
                         ABy*(-HRR(11)+&
                         HRRB(36))+&
                         ABx*(2.D0*ABy*HRRB(21)+&
                         ABx*(ABy*HRRB(11)+&
                         HRRB(22))+&
                         2.D0*HRRB(37))+&
                         HRRB(58)
      !=(11_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(37)+&
                         ABx*(ABy*HRRA(22)+&
                         HRRA(38))+&
                         HRRA(59)
      !=(11,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(21)+&
                         ABy*(ABy*HRRB(21)+&
                         2.D0*HRRB(37))+&
                         ABx*(-HRR(11)+&
                         ABy*(ABy*HRRB(11)+&
                         2.D0*HRRB(22))+&
                         HRRB(38))+&
                         HRRB(59)
      !=(11_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(42)+&
                         ABx*(ABy*HRRA(26)+&
                         HRRA(43))+&
                         HRRA(65)
      !=(11,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(37)+&
                         ABy*(ABz*HRRB(21)+&
                         HRRB(42))+&
                         ABx*(ABz*HRRB(22)+&
                         ABy*(ABz*HRRB(11)+&
                         HRRB(26))+&
                         HRRB(43))+&
                         HRRB(65)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(12_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(23)+&
                         ABy*(-2.D0*HRR(12)+&
                         HRRA(37))+&
                         ABx*(-2.D0*HRR(13)+&
                         ABy*(-2.D0*HRR(6)+&
                         HRRA(22))+&
                         HRRA(38))+&
                         HRRA(59)
      !=(12,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(23)+&
                         ABy*(-HRR(12)+&
                         HRRB(37))+&
                         ABx*(2.D0*ABy*HRRB(22)+&
                         ABx*(ABy*HRRB(12)+&
                         HRRB(23))+&
                         2.D0*HRRB(38))+&
                         HRRB(59)
      !=(12_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(22)+&
                         ABy*(-HRR(11)+&
                         HRRA(38))+&
                         ABx*(-HRR(12)+&
                         ABy*(-HRR(5)+&
                         HRRA(23))+&
                         HRRA(39))+&
                         HRRA(60)
      !=(12,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(22)+&
                         ABy*(ABy*HRRB(22)+&
                         2.D0*HRRB(38))+&
                         ABx*(-HRR(12)+&
                         ABy*(ABy*HRRB(12)+&
                         2.D0*HRRB(23))+&
                         HRRB(39))+&
                         HRRB(60)
      !=(12_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(43)+&
                         ABx*(ABy*HRRA(27)+&
                         HRRA(44))+&
                         HRRA(66)
      !=(12,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(38)+&
                         ABy*(ABz*HRRB(22)+&
                         HRRB(43))+&
                         ABx*(ABz*HRRB(23)+&
                         ABy*(ABz*HRRB(12)+&
                         HRRB(27))+&
                         HRRB(44))+&
                         HRRB(66)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(13_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(24)+&
                         ABy*(-HRR(13)+&
                         HRRA(38))+&
                         ABx*(-HRR(14)+&
                         ABy*(-HRR(7)+&
                         HRRA(23))+&
                         HRRA(39))+&
                         HRRA(60)
      !=(13,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(24)+&
                         ABy*(-HRR(13)+&
                         HRRB(38))+&
                         ABx*(2.D0*ABy*HRRB(23)+&
                         ABx*(ABy*HRRB(13)+&
                         HRRB(24))+&
                         2.D0*HRRB(39))+&
                         HRRB(60)
      !=(13_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(23)+&
                         ABy*(-2.D0*HRR(12)+&
                         HRRA(39))+&
                         ABx*(-2.D0*HRR(13)+&
                         ABy*(-2.D0*HRR(6)+&
                         HRRA(24))+&
                         HRRA(40))+&
                         HRRA(61)
      !=(13,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(23)+&
                         ABy*(ABy*HRRB(23)+&
                         2.D0*HRRB(39))+&
                         ABx*(-HRR(13)+&
                         ABy*(ABy*HRRB(13)+&
                         2.D0*HRRB(24))+&
                         HRRB(40))+&
                         HRRB(61)
      !=(13_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(44)+&
                         ABx*(ABy*HRRA(28)+&
                         HRRA(45))+&
                         HRRA(67)
      !=(13,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(39)+&
                         ABy*(ABz*HRRB(23)+&
                         HRRB(44))+&
                         ABx*(ABz*HRRB(24)+&
                         ABy*(ABz*HRRB(13)+&
                         HRRB(28))+&
                         HRRB(45))+&
                         HRRB(67)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(14_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(39)+&
                         ABx*(ABy*HRRA(24)+&
                         HRRA(40))+&
                         HRRA(61)
      !=(14,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(25)+&
                         ABy*(-HRR(14)+&
                         HRRB(39))+&
                         ABx*(2.D0*ABy*HRRB(24)+&
                         ABx*(ABy*HRRB(14)+&
                         HRRB(25))+&
                         2.D0*HRRB(40))+&
                         HRRB(61)
      !=(14_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-3.D0*HRR(24)+&
                         ABy*(-3.D0*HRR(13)+&
                         HRRA(40))+&
                         ABx*(-3.D0*HRR(14)+&
                         ABy*(-3.D0*HRR(7)+&
                         HRRA(25))+&
                         HRRA(41))+&
                         HRRA(62)
      !=(14,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(24)+&
                         ABy*(ABy*HRRB(24)+&
                         2.D0*HRRB(40))+&
                         ABx*(-HRR(14)+&
                         ABy*(ABy*HRRB(14)+&
                         2.D0*HRRB(25))+&
                         HRRB(41))+&
                         HRRB(62)
      !=(14_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(45)+&
                         ABx*(ABy*HRRA(29)+&
                         HRRA(46))+&
                         HRRA(68)
      !=(14,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(40)+&
                         ABy*(ABz*HRRB(24)+&
                         HRRB(45))+&
                         ABx*(ABz*HRRB(25)+&
                         ABy*(ABz*HRRB(14)+&
                         HRRB(29))+&
                         HRRB(46))+&
                         HRRB(68)
      OffSet=(OA+4)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(15_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(27)+&
                         ABy*(-2.D0*HRR(15)+&
                         HRRA(42))+&
                         ABx*(-2.D0*HRR(16)+&
                         ABy*(-2.D0*HRR(8)+&
                         HRRA(26))+&
                         HRRA(43))+&
                         HRRA(65)
      !=(15,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(27)+&
                         ABy*(-HRR(15)+&
                         HRRB(42))+&
                         ABx*(2.D0*ABy*HRRB(26)+&
                         ABx*(ABy*HRRB(15)+&
                         HRRB(27))+&
                         2.D0*HRRB(43))+&
                         HRRB(65)
      !=(15_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(43)+&
                         ABx*(ABy*HRRA(27)+&
                         HRRA(44))+&
                         HRRA(66)
      !=(15,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(26)+&
                         ABy*(ABy*HRRB(26)+&
                         2.D0*HRRB(43))+&
                         ABx*(-HRR(15)+&
                         ABy*(ABy*HRRB(15)+&
                         2.D0*HRRB(27))+&
                         HRRB(44))+&
                         HRRB(66)
      !=(15_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(22)+&
                         ABy*(-HRR(11)+&
                         HRRA(47))+&
                         ABx*(-HRR(12)+&
                         ABy*(-HRR(5)+&
                         HRRA(30))+&
                         HRRA(48))+&
                         HRRA(71)
      !=(15,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(43)+&
                         ABy*(ABz*HRRB(26)+&
                         HRRB(47))+&
                         ABx*(ABz*HRRB(27)+&
                         ABy*(ABz*HRRB(15)+&
                         HRRB(30))+&
                         HRRB(48))+&
                         HRRB(71)
      OffSet=(OA+5)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(16_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(28)+&
                         ABy*(-HRR(16)+&
                         HRRA(43))+&
                         ABx*(-HRR(17)+&
                         ABy*(-HRR(9)+&
                         HRRA(27))+&
                         HRRA(44))+&
                         HRRA(66)
      !=(16,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(28)+&
                         ABy*(-HRR(16)+&
                         HRRB(43))+&
                         ABx*(2.D0*ABy*HRRB(27)+&
                         ABx*(ABy*HRRB(16)+&
                         HRRB(28))+&
                         2.D0*HRRB(44))+&
                         HRRB(66)
      !=(16_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(27)+&
                         ABy*(-HRR(15)+&
                         HRRA(44))+&
                         ABx*(-HRR(16)+&
                         ABy*(-HRR(8)+&
                         HRRA(28))+&
                         HRRA(45))+&
                         HRRA(67)
      !=(16,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(27)+&
                         ABy*(ABy*HRRB(27)+&
                         2.D0*HRRB(44))+&
                         ABx*(-HRR(16)+&
                         ABy*(ABy*HRRB(16)+&
                         2.D0*HRRB(28))+&
                         HRRB(45))+&
                         HRRB(67)
      !=(16_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(23)+&
                         ABy*(-HRR(12)+&
                         HRRA(48))+&
                         ABx*(-HRR(13)+&
                         ABy*(-HRR(6)+&
                         HRRA(31))+&
                         HRRA(49))+&
                         HRRA(72)
      !=(16,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(44)+&
                         ABy*(ABz*HRRB(27)+&
                         HRRB(48))+&
                         ABx*(ABz*HRRB(28)+&
                         ABy*(ABz*HRRB(16)+&
                         HRRB(31))+&
                         HRRB(49))+&
                         HRRB(72)
      OffSet=(OA+6)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(17_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(44)+&
                         ABx*(ABy*HRRA(28)+&
                         HRRA(45))+&
                         HRRA(67)
      !=(17,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(29)+&
                         ABy*(-HRR(17)+&
                         HRRB(44))+&
                         ABx*(2.D0*ABy*HRRB(28)+&
                         ABx*(ABy*HRRB(17)+&
                         HRRB(29))+&
                         2.D0*HRRB(45))+&
                         HRRB(67)
      !=(17_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(28)+&
                         ABy*(-2.D0*HRR(16)+&
                         HRRA(45))+&
                         ABx*(-2.D0*HRR(17)+&
                         ABy*(-2.D0*HRR(9)+&
                         HRRA(29))+&
                         HRRA(46))+&
                         HRRA(68)
      !=(17,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(28)+&
                         ABy*(ABy*HRRB(28)+&
                         2.D0*HRRB(45))+&
                         ABx*(-HRR(17)+&
                         ABy*(ABy*HRRB(17)+&
                         2.D0*HRRB(29))+&
                         HRRB(46))+&
                         HRRB(68)
      !=(17_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(24)+&
                         ABy*(-HRR(13)+&
                         HRRA(49))+&
                         ABx*(-HRR(14)+&
                         ABy*(-HRR(7)+&
                         HRRA(32))+&
                         HRRA(50))+&
                         HRRA(73)
      !=(17,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(45)+&
                         ABy*(ABz*HRRB(28)+&
                         HRRB(49))+&
                         ABx*(ABz*HRRB(29)+&
                         ABy*(ABz*HRRB(17)+&
                         HRRB(32))+&
                         HRRB(50))+&
                         HRRB(73)
      OffSet=(OA+7)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(18_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(31)+&
                         ABy*(-HRR(18)+&
                         HRRA(47))+&
                         ABx*(-HRR(19)+&
                         ABy*(-HRR(10)+&
                         HRRA(30))+&
                         HRRA(48))+&
                         HRRA(71)
      !=(18,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(31)+&
                         ABy*(-HRR(18)+&
                         HRRB(47))+&
                         ABx*(2.D0*ABy*HRRB(30)+&
                         ABx*(ABy*HRRB(18)+&
                         HRRB(31))+&
                         2.D0*HRRB(48))+&
                         HRRB(71)
      !=(18_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(48)+&
                         ABx*(ABy*HRRA(31)+&
                         HRRA(49))+&
                         HRRA(72)
      !=(18,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(30)+&
                         ABy*(ABy*HRRB(30)+&
                         2.D0*HRRB(48))+&
                         ABx*(-HRR(18)+&
                         ABy*(ABy*HRRB(18)+&
                         2.D0*HRRB(31))+&
                         HRRB(49))+&
                         HRRB(72)
      !=(18_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(27)+&
                         ABy*(-2.D0*HRR(15)+&
                         HRRA(51))+&
                         ABx*(-2.D0*HRR(16)+&
                         ABy*(-2.D0*HRR(8)+&
                         HRRA(33))+&
                         HRRA(52))+&
                         HRRA(76)
      !=(18,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(48)+&
                         ABy*(ABz*HRRB(30)+&
                         HRRB(51))+&
                         ABx*(ABz*HRRB(31)+&
                         ABy*(ABz*HRRB(18)+&
                         HRRB(33))+&
                         HRRB(52))+&
                         HRRB(76)
      OffSet=(OA+8)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(19_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(48)+&
                         ABx*(ABy*HRRA(31)+&
                         HRRA(49))+&
                         HRRA(72)
      !=(19,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(32)+&
                         ABy*(-HRR(19)+&
                         HRRB(48))+&
                         ABx*(2.D0*ABy*HRRB(31)+&
                         ABx*(ABy*HRRB(19)+&
                         HRRB(32))+&
                         2.D0*HRRB(49))+&
                         HRRB(72)
      !=(19_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(31)+&
                         ABy*(-HRR(18)+&
                         HRRA(49))+&
                         ABx*(-HRR(19)+&
                         ABy*(-HRR(10)+&
                         HRRA(32))+&
                         HRRA(50))+&
                         HRRA(73)
      !=(19,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(31)+&
                         ABy*(ABy*HRRB(31)+&
                         2.D0*HRRB(49))+&
                         ABx*(-HRR(19)+&
                         ABy*(ABy*HRRB(19)+&
                         2.D0*HRRB(32))+&
                         HRRB(50))+&
                         HRRB(73)
      !=(19_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(28)+&
                         ABy*(-2.D0*HRR(16)+&
                         HRRA(52))+&
                         ABx*(-2.D0*HRR(17)+&
                         ABy*(-2.D0*HRR(9)+&
                         HRRA(34))+&
                         HRRA(53))+&
                         HRRA(77)
      !=(19,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(49)+&
                         ABy*(ABz*HRRB(31)+&
                         HRRB(52))+&
                         ABx*(ABz*HRRB(32)+&
                         ABy*(ABz*HRRB(19)+&
                         HRRB(34))+&
                         HRRB(53))+&
                         HRRB(77)
      OffSet=(OA+9)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(20_x,6|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(51)+&
                         ABx*(ABy*HRRA(33)+&
                         HRRA(52))+&
                         HRRA(76)
      !=(20,6_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(34)+&
                         ABy*(-HRR(20)+&
                         HRRB(51))+&
                         ABx*(2.D0*ABy*HRRB(33)+&
                         ABx*(ABy*HRRB(20)+&
                         HRRB(34))+&
                         2.D0*HRRB(52))+&
                         HRRB(76)
      !=(20_y,6|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(52)+&
                         ABx*(ABy*HRRA(34)+&
                         HRRA(53))+&
                         HRRA(77)
      !=(20,6_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(33)+&
                         ABy*(ABy*HRRB(33)+&
                         2.D0*HRRB(52))+&
                         ABx*(-HRR(20)+&
                         ABy*(ABy*HRRB(20)+&
                         2.D0*HRRB(34))+&
                         HRRB(53))+&
                         HRRB(77)
      !=(20_z,6|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-3.D0*HRR(31)+&
                         ABy*(-3.D0*HRR(18)+&
                         HRRA(54))+&
                         ABx*(-3.D0*HRR(19)+&
                         ABy*(-3.D0*HRR(10)+&
                         HRRA(35))+&
                         HRRA(55))+&
                         HRRA(80)
      !=(20,6_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(52)+&
                         ABy*(ABz*HRRB(33)+&
                         HRRB(54))+&
                         ABx*(ABz*HRRB(34)+&
                         ABy*(ABz*HRRB(20)+&
                         HRRB(35))+&
                         HRRB(55))+&
                         HRRB(80)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(11_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-3.D0*HRR(23)+&
                         ABy*(-6.D0*HRR(12)+&
                         ABy*(-3.D0*HRR(5)+&
                         HRRA(21))+&
                         2.D0*HRRA(37))+&
                         HRRA(59)
      !=(11,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(21)+&
                         2.D0*HRRB(37))+&
                         ABx*(ABy*(ABy*HRRB(11)+&
                         2.D0*HRRB(22))+&
                         HRRB(38))+&
                         HRRB(59)
      !=(11_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*(ABy*HRRA(22)+&
                         2.D0*HRRA(38))+&
                         HRRA(60)
      !=(11,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(22)+&
                         ABy*(-2.D0*HRR(11)+&
                         ABy*(ABy*HRRB(11)+&
                         3.D0*HRRB(22))+&
                         3.D0*HRRB(38))+&
                         HRRB(60)
      !=(11_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*(ABy*HRRA(26)+&
                         2.D0*HRRA(43))+&
                         HRRA(66)
      !=(11,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(38)+&
                         ABy*(2.D0*ABz*HRRB(22)+&
                         ABy*(ABz*HRRB(11)+&
                         HRRB(26))+&
                         2.D0*HRRB(43))+&
                         HRRB(66)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(12_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(24)+&
                         ABy*(-4.D0*HRR(13)+&
                         ABy*(-2.D0*HRR(6)+&
                         HRRA(22))+&
                         2.D0*HRRA(38))+&
                         HRRA(60)
      !=(12,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(22)+&
                         2.D0*HRRB(38))+&
                         ABx*(ABy*(ABy*HRRB(12)+&
                         2.D0*HRRB(23))+&
                         HRRB(39))+&
                         HRRB(60)
      !=(12_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(23)+&
                         ABy*(-2.D0*HRR(12)+&
                         ABy*(-HRR(5)+&
                         HRRA(23))+&
                         2.D0*HRRA(39))+&
                         HRRA(61)
      !=(12,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(23)+&
                         ABy*(-2.D0*HRR(12)+&
                         ABy*(ABy*HRRB(12)+&
                         3.D0*HRRB(23))+&
                         3.D0*HRRB(39))+&
                         HRRB(61)
      !=(12_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*(ABy*HRRA(27)+&
                         2.D0*HRRA(44))+&
                         HRRA(67)
      !=(12,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(39)+&
                         ABy*(2.D0*ABz*HRRB(23)+&
                         ABy*(ABz*HRRB(12)+&
                         HRRB(27))+&
                         2.D0*HRRB(44))+&
                         HRRB(67)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(13_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(25)+&
                         ABy*(-2.D0*HRR(14)+&
                         ABy*(-HRR(7)+&
                         HRRA(23))+&
                         2.D0*HRRA(39))+&
                         HRRA(61)
      !=(13,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(23)+&
                         2.D0*HRRB(39))+&
                         ABx*(ABy*(ABy*HRRB(13)+&
                         2.D0*HRRB(24))+&
                         HRRB(40))+&
                         HRRB(61)
      !=(13_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(24)+&
                         ABy*(-4.D0*HRR(13)+&
                         ABy*(-2.D0*HRR(6)+&
                         HRRA(24))+&
                         2.D0*HRRA(40))+&
                         HRRA(62)
      !=(13,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(24)+&
                         ABy*(-2.D0*HRR(13)+&
                         ABy*(ABy*HRRB(13)+&
                         3.D0*HRRB(24))+&
                         3.D0*HRRB(40))+&
                         HRRB(62)
      !=(13_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*(ABy*HRRA(28)+&
                         2.D0*HRRA(45))+&
                         HRRA(68)
      !=(13,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(40)+&
                         ABy*(2.D0*ABz*HRRB(24)+&
                         ABy*(ABz*HRRB(13)+&
                         HRRB(28))+&
                         2.D0*HRRB(45))+&
                         HRRB(68)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(14_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*(ABy*HRRA(24)+&
                         2.D0*HRRA(40))+&
                         HRRA(62)
      !=(14,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(24)+&
                         2.D0*HRRB(40))+&
                         ABx*(ABy*(ABy*HRRB(14)+&
                         2.D0*HRRB(25))+&
                         HRRB(41))+&
                         HRRB(62)
      !=(14_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-3.D0*HRR(25)+&
                         ABy*(-6.D0*HRR(14)+&
                         ABy*(-3.D0*HRR(7)+&
                         HRRA(25))+&
                         2.D0*HRRA(41))+&
                         HRRA(63)
      !=(14,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(25)+&
                         ABy*(-2.D0*HRR(14)+&
                         ABy*(ABy*HRRB(14)+&
                         3.D0*HRRB(25))+&
                         3.D0*HRRB(41))+&
                         HRRB(63)
      !=(14_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*(ABy*HRRA(29)+&
                         2.D0*HRRA(46))+&
                         HRRA(69)
      !=(14,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(41)+&
                         ABy*(2.D0*ABz*HRRB(25)+&
                         ABy*(ABz*HRRB(14)+&
                         HRRB(29))+&
                         2.D0*HRRB(46))+&
                         HRRB(69)
      OffSet=(OA+4)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(15_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(28)+&
                         ABy*(-4.D0*HRR(16)+&
                         ABy*(-2.D0*HRR(8)+&
                         HRRA(26))+&
                         2.D0*HRRA(43))+&
                         HRRA(66)
      !=(15,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(26)+&
                         2.D0*HRRB(43))+&
                         ABx*(ABy*(ABy*HRRB(15)+&
                         2.D0*HRRB(27))+&
                         HRRB(44))+&
                         HRRB(66)
      !=(15_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*(ABy*HRRA(27)+&
                         2.D0*HRRA(44))+&
                         HRRA(67)
      !=(15,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(27)+&
                         ABy*(-2.D0*HRR(15)+&
                         ABy*(ABy*HRRB(15)+&
                         3.D0*HRRB(27))+&
                         3.D0*HRRB(44))+&
                         HRRB(67)
      !=(15_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(23)+&
                         ABy*(-2.D0*HRR(12)+&
                         ABy*(-HRR(5)+&
                         HRRA(30))+&
                         2.D0*HRRA(48))+&
                         HRRA(72)
      !=(15,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(44)+&
                         ABy*(2.D0*ABz*HRRB(27)+&
                         ABy*(ABz*HRRB(15)+&
                         HRRB(30))+&
                         2.D0*HRRB(48))+&
                         HRRB(72)
      OffSet=(OA+5)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(16_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(29)+&
                         ABy*(-2.D0*HRR(17)+&
                         ABy*(-HRR(9)+&
                         HRRA(27))+&
                         2.D0*HRRA(44))+&
                         HRRA(67)
      !=(16,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(27)+&
                         2.D0*HRRB(44))+&
                         ABx*(ABy*(ABy*HRRB(16)+&
                         2.D0*HRRB(28))+&
                         HRRB(45))+&
                         HRRB(67)
      !=(16_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(28)+&
                         ABy*(-2.D0*HRR(16)+&
                         ABy*(-HRR(8)+&
                         HRRA(28))+&
                         2.D0*HRRA(45))+&
                         HRRA(68)
      !=(16,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(28)+&
                         ABy*(-2.D0*HRR(16)+&
                         ABy*(ABy*HRRB(16)+&
                         3.D0*HRRB(28))+&
                         3.D0*HRRB(45))+&
                         HRRB(68)
      !=(16_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(24)+&
                         ABy*(-2.D0*HRR(13)+&
                         ABy*(-HRR(6)+&
                         HRRA(31))+&
                         2.D0*HRRA(49))+&
                         HRRA(73)
      !=(16,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(45)+&
                         ABy*(2.D0*ABz*HRRB(28)+&
                         ABy*(ABz*HRRB(16)+&
                         HRRB(31))+&
                         2.D0*HRRB(49))+&
                         HRRB(73)
      OffSet=(OA+6)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(17_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*(ABy*HRRA(28)+&
                         2.D0*HRRA(45))+&
                         HRRA(68)
      !=(17,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(28)+&
                         2.D0*HRRB(45))+&
                         ABx*(ABy*(ABy*HRRB(17)+&
                         2.D0*HRRB(29))+&
                         HRRB(46))+&
                         HRRB(68)
      !=(17_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(29)+&
                         ABy*(-4.D0*HRR(17)+&
                         ABy*(-2.D0*HRR(9)+&
                         HRRA(29))+&
                         2.D0*HRRA(46))+&
                         HRRA(69)
      !=(17,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(29)+&
                         ABy*(-2.D0*HRR(17)+&
                         ABy*(ABy*HRRB(17)+&
                         3.D0*HRRB(29))+&
                         3.D0*HRRB(46))+&
                         HRRB(69)
      !=(17_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(25)+&
                         ABy*(-2.D0*HRR(14)+&
                         ABy*(-HRR(7)+&
                         HRRA(32))+&
                         2.D0*HRRA(50))+&
                         HRRA(74)
      !=(17,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(46)+&
                         ABy*(2.D0*ABz*HRRB(29)+&
                         ABy*(ABz*HRRB(17)+&
                         HRRB(32))+&
                         2.D0*HRRB(50))+&
                         HRRB(74)
      OffSet=(OA+7)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(18_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(32)+&
                         ABy*(-2.D0*HRR(19)+&
                         ABy*(-HRR(10)+&
                         HRRA(30))+&
                         2.D0*HRRA(48))+&
                         HRRA(72)
      !=(18,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(30)+&
                         2.D0*HRRB(48))+&
                         ABx*(ABy*(ABy*HRRB(18)+&
                         2.D0*HRRB(31))+&
                         HRRB(49))+&
                         HRRB(72)
      !=(18_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*(ABy*HRRA(31)+&
                         2.D0*HRRA(49))+&
                         HRRA(73)
      !=(18,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(31)+&
                         ABy*(-2.D0*HRR(18)+&
                         ABy*(ABy*HRRB(18)+&
                         3.D0*HRRB(31))+&
                         3.D0*HRRB(49))+&
                         HRRB(73)
      !=(18_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(28)+&
                         ABy*(-4.D0*HRR(16)+&
                         ABy*(-2.D0*HRR(8)+&
                         HRRA(33))+&
                         2.D0*HRRA(52))+&
                         HRRA(77)
      !=(18,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(49)+&
                         ABy*(2.D0*ABz*HRRB(31)+&
                         ABy*(ABz*HRRB(18)+&
                         HRRB(33))+&
                         2.D0*HRRB(52))+&
                         HRRB(77)
      OffSet=(OA+8)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(19_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*(ABy*HRRA(31)+&
                         2.D0*HRRA(49))+&
                         HRRA(73)
      !=(19,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(31)+&
                         2.D0*HRRB(49))+&
                         ABx*(ABy*(ABy*HRRB(19)+&
                         2.D0*HRRB(32))+&
                         HRRB(50))+&
                         HRRB(73)
      !=(19_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(32)+&
                         ABy*(-2.D0*HRR(19)+&
                         ABy*(-HRR(10)+&
                         HRRA(32))+&
                         2.D0*HRRA(50))+&
                         HRRA(74)
      !=(19,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(32)+&
                         ABy*(-2.D0*HRR(19)+&
                         ABy*(ABy*HRRB(19)+&
                         3.D0*HRRB(32))+&
                         3.D0*HRRB(50))+&
                         HRRB(74)
      !=(19_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(29)+&
                         ABy*(-4.D0*HRR(17)+&
                         ABy*(-2.D0*HRR(9)+&
                         HRRA(34))+&
                         2.D0*HRRA(53))+&
                         HRRA(78)
      !=(19,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(50)+&
                         ABy*(2.D0*ABz*HRRB(32)+&
                         ABy*(ABz*HRRB(19)+&
                         HRRB(34))+&
                         2.D0*HRRB(53))+&
                         HRRB(78)
      OffSet=(OA+9)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(20_x,7|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*(ABy*HRRA(33)+&
                         2.D0*HRRA(52))+&
                         HRRA(77)
      !=(20,7_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*(ABy*HRRB(33)+&
                         2.D0*HRRB(52))+&
                         ABx*(ABy*(ABy*HRRB(20)+&
                         2.D0*HRRB(34))+&
                         HRRB(53))+&
                         HRRB(77)
      !=(20_y,7|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*(ABy*HRRA(34)+&
                         2.D0*HRRA(53))+&
                         HRRA(78)
      !=(20,7_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-2.D0*HRR(34)+&
                         ABy*(-2.D0*HRR(20)+&
                         ABy*(ABy*HRRB(20)+&
                         3.D0*HRRB(34))+&
                         3.D0*HRRB(53))+&
                         HRRB(78)
      !=(20_z,7|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-3.D0*HRR(32)+&
                         ABy*(-6.D0*HRR(19)+&
                         ABy*(-3.D0*HRR(10)+&
                         HRRA(35))+&
                         2.D0*HRRA(55))+&
                         HRRA(81)
      !=(20,7_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(53)+&
                         ABy*(2.D0*ABz*HRRB(34)+&
                         ABy*(ABz*HRRB(20)+&
                         HRRB(35))+&
                         2.D0*HRRB(55))+&
                         HRRB(81)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(11_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-3.D0*HRR(26)+&
                         ABz*(-3.D0*HRR(11)+&
                         HRRA(36))+&
                         ABx*(-3.D0*HRR(15)+&
                         ABz*(-3.D0*HRR(5)+&
                         HRRA(21))+&
                         HRRA(42))+&
                         HRRA(64)
      !=(11,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(26)+&
                         ABz*(-HRR(11)+&
                         HRRB(36))+&
                         ABx*(2.D0*ABz*HRRB(21)+&
                         ABx*(ABz*HRRB(11)+&
                         HRRB(26))+&
                         2.D0*HRRB(42))+&
                         HRRB(64)
      !=(11_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(37)+&
                         ABx*(ABz*HRRA(22)+&
                         HRRA(43))+&
                         HRRA(65)
      !=(11,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(37)+&
                         ABy*(ABz*HRRB(21)+&
                         HRRB(42))+&
                         ABx*(ABz*HRRB(22)+&
                         ABy*(ABz*HRRB(11)+&
                         HRRB(26))+&
                         HRRB(43))+&
                         HRRB(65)
      !=(11_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(42)+&
                         ABx*(ABz*HRRA(26)+&
                         HRRA(47))+&
                         HRRA(70)
      !=(11,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(21)+&
                         ABz*(ABz*HRRB(21)+&
                         2.D0*HRRB(42))+&
                         ABx*(-HRR(11)+&
                         ABz*(ABz*HRRB(11)+&
                         2.D0*HRRB(26))+&
                         HRRB(47))+&
                         HRRB(70)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(12_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(27)+&
                         ABz*(-2.D0*HRR(12)+&
                         HRRA(37))+&
                         ABx*(-2.D0*HRR(16)+&
                         ABz*(-2.D0*HRR(6)+&
                         HRRA(22))+&
                         HRRA(43))+&
                         HRRA(65)
      !=(12,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(27)+&
                         ABz*(-HRR(12)+&
                         HRRB(37))+&
                         ABx*(2.D0*ABz*HRRB(22)+&
                         ABx*(ABz*HRRB(12)+&
                         HRRB(27))+&
                         2.D0*HRRB(43))+&
                         HRRB(65)
      !=(12_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(26)+&
                         ABz*(-HRR(11)+&
                         HRRA(38))+&
                         ABx*(-HRR(15)+&
                         ABz*(-HRR(5)+&
                         HRRA(23))+&
                         HRRA(44))+&
                         HRRA(66)
      !=(12,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(38)+&
                         ABy*(ABz*HRRB(22)+&
                         HRRB(43))+&
                         ABx*(ABz*HRRB(23)+&
                         ABy*(ABz*HRRB(12)+&
                         HRRB(27))+&
                         HRRB(44))+&
                         HRRB(66)
      !=(12_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(43)+&
                         ABx*(ABz*HRRA(27)+&
                         HRRA(48))+&
                         HRRA(71)
      !=(12,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(22)+&
                         ABz*(ABz*HRRB(22)+&
                         2.D0*HRRB(43))+&
                         ABx*(-HRR(12)+&
                         ABz*(ABz*HRRB(12)+&
                         2.D0*HRRB(27))+&
                         HRRB(48))+&
                         HRRB(71)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(13_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(28)+&
                         ABz*(-HRR(13)+&
                         HRRA(38))+&
                         ABx*(-HRR(17)+&
                         ABz*(-HRR(7)+&
                         HRRA(23))+&
                         HRRA(44))+&
                         HRRA(66)
      !=(13,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(28)+&
                         ABz*(-HRR(13)+&
                         HRRB(38))+&
                         ABx*(2.D0*ABz*HRRB(23)+&
                         ABx*(ABz*HRRB(13)+&
                         HRRB(28))+&
                         2.D0*HRRB(44))+&
                         HRRB(66)
      !=(13_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(27)+&
                         ABz*(-2.D0*HRR(12)+&
                         HRRA(39))+&
                         ABx*(-2.D0*HRR(16)+&
                         ABz*(-2.D0*HRR(6)+&
                         HRRA(24))+&
                         HRRA(45))+&
                         HRRA(67)
      !=(13,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(39)+&
                         ABy*(ABz*HRRB(23)+&
                         HRRB(44))+&
                         ABx*(ABz*HRRB(24)+&
                         ABy*(ABz*HRRB(13)+&
                         HRRB(28))+&
                         HRRB(45))+&
                         HRRB(67)
      !=(13_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(44)+&
                         ABx*(ABz*HRRA(28)+&
                         HRRA(49))+&
                         HRRA(72)
      !=(13,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(23)+&
                         ABz*(ABz*HRRB(23)+&
                         2.D0*HRRB(44))+&
                         ABx*(-HRR(13)+&
                         ABz*(ABz*HRRB(13)+&
                         2.D0*HRRB(28))+&
                         HRRB(49))+&
                         HRRB(72)
      OffSet=(OA+3)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(14_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(39)+&
                         ABx*(ABz*HRRA(24)+&
                         HRRA(45))+&
                         HRRA(67)
      !=(14,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(29)+&
                         ABz*(-HRR(14)+&
                         HRRB(39))+&
                         ABx*(2.D0*ABz*HRRB(24)+&
                         ABx*(ABz*HRRB(14)+&
                         HRRB(29))+&
                         2.D0*HRRB(45))+&
                         HRRB(67)
      !=(14_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-3.D0*HRR(28)+&
                         ABz*(-3.D0*HRR(13)+&
                         HRRA(40))+&
                         ABx*(-3.D0*HRR(17)+&
                         ABz*(-3.D0*HRR(7)+&
                         HRRA(25))+&
                         HRRA(46))+&
                         HRRA(68)
      !=(14,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(40)+&
                         ABy*(ABz*HRRB(24)+&
                         HRRB(45))+&
                         ABx*(ABz*HRRB(25)+&
                         ABy*(ABz*HRRB(14)+&
                         HRRB(29))+&
                         HRRB(46))+&
                         HRRB(68)
      !=(14_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(45)+&
                         ABx*(ABz*HRRA(29)+&
                         HRRA(50))+&
                         HRRA(73)
      !=(14,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(24)+&
                         ABz*(ABz*HRRB(24)+&
                         2.D0*HRRB(45))+&
                         ABx*(-HRR(14)+&
                         ABz*(ABz*HRRB(14)+&
                         2.D0*HRRB(29))+&
                         HRRB(50))+&
                         HRRB(73)
      OffSet=(OA+4)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(15_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(30)+&
                         ABz*(-2.D0*HRR(15)+&
                         HRRA(42))+&
                         ABx*(-2.D0*HRR(18)+&
                         ABz*(-2.D0*HRR(8)+&
                         HRRA(26))+&
                         HRRA(47))+&
                         HRRA(70)
      !=(15,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(30)+&
                         ABz*(-HRR(15)+&
                         HRRB(42))+&
                         ABx*(2.D0*ABz*HRRB(26)+&
                         ABx*(ABz*HRRB(15)+&
                         HRRB(30))+&
                         2.D0*HRRB(47))+&
                         HRRB(70)
      !=(15_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(43)+&
                         ABx*(ABz*HRRA(27)+&
                         HRRA(48))+&
                         HRRA(71)
      !=(15,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(43)+&
                         ABy*(ABz*HRRB(26)+&
                         HRRB(47))+&
                         ABx*(ABz*HRRB(27)+&
                         ABy*(ABz*HRRB(15)+&
                         HRRB(30))+&
                         HRRB(48))+&
                         HRRB(71)
      !=(15_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(26)+&
                         ABz*(-HRR(11)+&
                         HRRA(47))+&
                         ABx*(-HRR(15)+&
                         ABz*(-HRR(5)+&
                         HRRA(30))+&
                         HRRA(51))+&
                         HRRA(75)
      !=(15,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(26)+&
                         ABz*(ABz*HRRB(26)+&
                         2.D0*HRRB(47))+&
                         ABx*(-HRR(15)+&
                         ABz*(ABz*HRRB(15)+&
                         2.D0*HRRB(30))+&
                         HRRB(51))+&
                         HRRB(75)
      OffSet=(OA+5)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(16_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(31)+&
                         ABz*(-HRR(16)+&
                         HRRA(43))+&
                         ABx*(-HRR(19)+&
                         ABz*(-HRR(9)+&
                         HRRA(27))+&
                         HRRA(48))+&
                         HRRA(71)
      !=(16,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(31)+&
                         ABz*(-HRR(16)+&
                         HRRB(43))+&
                         ABx*(2.D0*ABz*HRRB(27)+&
                         ABx*(ABz*HRRB(16)+&
                         HRRB(31))+&
                         2.D0*HRRB(48))+&
                         HRRB(71)
      !=(16_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(30)+&
                         ABz*(-HRR(15)+&
                         HRRA(44))+&
                         ABx*(-HRR(18)+&
                         ABz*(-HRR(8)+&
                         HRRA(28))+&
                         HRRA(49))+&
                         HRRA(72)
      !=(16,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(44)+&
                         ABy*(ABz*HRRB(27)+&
                         HRRB(48))+&
                         ABx*(ABz*HRRB(28)+&
                         ABy*(ABz*HRRB(16)+&
                         HRRB(31))+&
                         HRRB(49))+&
                         HRRB(72)
      !=(16_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(27)+&
                         ABz*(-HRR(12)+&
                         HRRA(48))+&
                         ABx*(-HRR(16)+&
                         ABz*(-HRR(6)+&
                         HRRA(31))+&
                         HRRA(52))+&
                         HRRA(76)
      !=(16,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(27)+&
                         ABz*(ABz*HRRB(27)+&
                         2.D0*HRRB(48))+&
                         ABx*(-HRR(16)+&
                         ABz*(ABz*HRRB(16)+&
                         2.D0*HRRB(31))+&
                         HRRB(52))+&
                         HRRB(76)
      OffSet=(OA+6)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(17_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(44)+&
                         ABx*(ABz*HRRA(28)+&
                         HRRA(49))+&
                         HRRA(72)
      !=(17,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(32)+&
                         ABz*(-HRR(17)+&
                         HRRB(44))+&
                         ABx*(2.D0*ABz*HRRB(28)+&
                         ABx*(ABz*HRRB(17)+&
                         HRRB(32))+&
                         2.D0*HRRB(49))+&
                         HRRB(72)
      !=(17_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(31)+&
                         ABz*(-2.D0*HRR(16)+&
                         HRRA(45))+&
                         ABx*(-2.D0*HRR(19)+&
                         ABz*(-2.D0*HRR(9)+&
                         HRRA(29))+&
                         HRRA(50))+&
                         HRRA(73)
      !=(17,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(45)+&
                         ABy*(ABz*HRRB(28)+&
                         HRRB(49))+&
                         ABx*(ABz*HRRB(29)+&
                         ABy*(ABz*HRRB(17)+&
                         HRRB(32))+&
                         HRRB(50))+&
                         HRRB(73)
      !=(17_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(28)+&
                         ABz*(-HRR(13)+&
                         HRRA(49))+&
                         ABx*(-HRR(17)+&
                         ABz*(-HRR(7)+&
                         HRRA(32))+&
                         HRRA(53))+&
                         HRRA(77)
      !=(17,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(28)+&
                         ABz*(ABz*HRRB(28)+&
                         2.D0*HRRB(49))+&
                         ABx*(-HRR(17)+&
                         ABz*(ABz*HRRB(17)+&
                         2.D0*HRRB(32))+&
                         HRRB(53))+&
                         HRRB(77)
      OffSet=(OA+7)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(18_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(33)+&
                         ABz*(-HRR(18)+&
                         HRRA(47))+&
                         ABx*(-HRR(20)+&
                         ABz*(-HRR(10)+&
                         HRRA(30))+&
                         HRRA(51))+&
                         HRRA(75)
      !=(18,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(33)+&
                         ABz*(-HRR(18)+&
                         HRRB(47))+&
                         ABx*(2.D0*ABz*HRRB(30)+&
                         ABx*(ABz*HRRB(18)+&
                         HRRB(33))+&
                         2.D0*HRRB(51))+&
                         HRRB(75)
      !=(18_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(48)+&
                         ABx*(ABz*HRRA(31)+&
                         HRRA(52))+&
                         HRRA(76)
      !=(18,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(48)+&
                         ABy*(ABz*HRRB(30)+&
                         HRRB(51))+&
                         ABx*(ABz*HRRB(31)+&
                         ABy*(ABz*HRRB(18)+&
                         HRRB(33))+&
                         HRRB(52))+&
                         HRRB(76)
      !=(18_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(30)+&
                         ABz*(-2.D0*HRR(15)+&
                         HRRA(51))+&
                         ABx*(-2.D0*HRR(18)+&
                         ABz*(-2.D0*HRR(8)+&
                         HRRA(33))+&
                         HRRA(54))+&
                         HRRA(79)
      !=(18,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(30)+&
                         ABz*(ABz*HRRB(30)+&
                         2.D0*HRRB(51))+&
                         ABx*(-HRR(18)+&
                         ABz*(ABz*HRRB(18)+&
                         2.D0*HRRB(33))+&
                         HRRB(54))+&
                         HRRB(79)
      OffSet=(OA+8)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(19_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(48)+&
                         ABx*(ABz*HRRA(31)+&
                         HRRA(52))+&
                         HRRA(76)
      !=(19,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(34)+&
                         ABz*(-HRR(19)+&
                         HRRB(48))+&
                         ABx*(2.D0*ABz*HRRB(31)+&
                         ABx*(ABz*HRRB(19)+&
                         HRRB(34))+&
                         2.D0*HRRB(52))+&
                         HRRB(76)
      !=(19_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(33)+&
                         ABz*(-HRR(18)+&
                         HRRA(49))+&
                         ABx*(-HRR(20)+&
                         ABz*(-HRR(10)+&
                         HRRA(32))+&
                         HRRA(53))+&
                         HRRA(77)
      !=(19,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(49)+&
                         ABy*(ABz*HRRB(31)+&
                         HRRB(52))+&
                         ABx*(ABz*HRRB(32)+&
                         ABy*(ABz*HRRB(19)+&
                         HRRB(34))+&
                         HRRB(53))+&
                         HRRB(77)
      !=(19_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(31)+&
                         ABz*(-2.D0*HRR(16)+&
                         HRRA(52))+&
                         ABx*(-2.D0*HRR(19)+&
                         ABz*(-2.D0*HRR(9)+&
                         HRRA(34))+&
                         HRRA(55))+&
                         HRRA(80)
      !=(19,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(31)+&
                         ABz*(ABz*HRRB(31)+&
                         2.D0*HRRB(52))+&
                         ABx*(-HRR(19)+&
                         ABz*(ABz*HRRB(19)+&
                         2.D0*HRRB(34))+&
                         HRRB(55))+&
                         HRRB(80)
      OffSet=(OA+9)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(20_x,8|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(51)+&
                         ABx*(ABz*HRRA(33)+&
                         HRRA(54))+&
                         HRRA(79)
      !=(20,8_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(35)+&
                         ABz*(-HRR(20)+&
                         HRRB(51))+&
                         ABx*(2.D0*ABz*HRRB(33)+&
                         ABx*(ABz*HRRB(20)+&
                         HRRB(35))+&
                         2.D0*HRRB(54))+&
                         HRRB(79)
      !=(20_y,8|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(52)+&
                         ABx*(ABz*HRRA(34)+&
                         HRRA(55))+&
                         HRRA(80)
      !=(20,8_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(52)+&
                         ABy*(ABz*HRRB(33)+&
                         HRRB(54))+&
                         ABx*(ABz*HRRB(34)+&
                         ABy*(ABz*HRRB(20)+&
                         HRRB(35))+&
                         HRRB(55))+&
                         HRRB(80)
      !=(20_z,8|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-3.D0*HRR(33)+&
                         ABz*(-3.D0*HRR(18)+&
                         HRRA(54))+&
                         ABx*(-3.D0*HRR(20)+&
                         ABz*(-3.D0*HRR(10)+&
                         HRRA(35))+&
                         HRRA(56))+&
                         HRRA(82)
      !=(20,8_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(33)+&
                         ABz*(ABz*HRRB(33)+&
                         2.D0*HRRB(54))+&
                         ABx*(-HRR(20)+&
                         ABz*(ABz*HRRB(20)+&
                         2.D0*HRRB(35))+&
                         HRRB(56))+&
                         HRRB(82)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(11_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-3.D0*HRR(27)+&
                         ABz*(-3.D0*HRR(12)+&
                         HRRA(37))+&
                         ABy*(-3.D0*HRR(15)+&
                         ABz*(-3.D0*HRR(5)+&
                         HRRA(21))+&
                         HRRA(42))+&
                         HRRA(65)
      !=(11,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(37)+&
                         ABy*(ABz*HRRB(21)+&
                         HRRB(42))+&
                         ABx*(ABz*HRRB(22)+&
                         ABy*(ABz*HRRB(11)+&
                         HRRB(26))+&
                         HRRB(43))+&
                         HRRB(65)
      !=(11_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(38)+&
                         ABy*(ABz*HRRA(22)+&
                         HRRA(43))+&
                         HRRA(66)
      !=(11,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(26)+&
                         ABz*(-HRR(11)+&
                         HRRB(38))+&
                         ABy*(2.D0*ABz*HRRB(22)+&
                         ABy*(ABz*HRRB(11)+&
                         HRRB(26))+&
                         2.D0*HRRB(43))+&
                         HRRB(66)
      !=(11_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(43)+&
                         ABy*(ABz*HRRA(26)+&
                         HRRA(47))+&
                         HRRA(71)
      !=(11,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(22)+&
                         ABz*(ABz*HRRB(22)+&
                         2.D0*HRRB(43))+&
                         ABy*(-HRR(11)+&
                         ABz*(ABz*HRRB(11)+&
                         2.D0*HRRB(26))+&
                         HRRB(47))+&
                         HRRB(71)
      OffSet=(OA+1)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(12_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(28)+&
                         ABz*(-2.D0*HRR(13)+&
                         HRRA(38))+&
                         ABy*(-2.D0*HRR(16)+&
                         ABz*(-2.D0*HRR(6)+&
                         HRRA(22))+&
                         HRRA(43))+&
                         HRRA(66)
      !=(12,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(38)+&
                         ABy*(ABz*HRRB(22)+&
                         HRRB(43))+&
                         ABx*(ABz*HRRB(23)+&
                         ABy*(ABz*HRRB(12)+&
                         HRRB(27))+&
                         HRRB(44))+&
                         HRRB(66)
      !=(12_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(27)+&
                         ABz*(-HRR(12)+&
                         HRRA(39))+&
                         ABy*(-HRR(15)+&
                         ABz*(-HRR(5)+&
                         HRRA(23))+&
                         HRRA(44))+&
                         HRRA(67)
      !=(12,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(27)+&
                         ABz*(-HRR(12)+&
                         HRRB(39))+&
                         ABy*(2.D0*ABz*HRRB(23)+&
                         ABy*(ABz*HRRB(12)+&
                         HRRB(27))+&
                         2.D0*HRRB(44))+&
                         HRRB(67)
      !=(12_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(44)+&
                         ABy*(ABz*HRRA(27)+&
                         HRRA(48))+&
                         HRRA(72)
      !=(12,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(23)+&
                         ABz*(ABz*HRRB(23)+&
                         2.D0*HRRB(44))+&
                         ABy*(-HRR(12)+&
                         ABz*(ABz*HRRB(12)+&
                         2.D0*HRRB(27))+&
                         HRRB(48))+&
                         HRRB(72)
      OffSet=(OA+2)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(13_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(29)+&
                         ABz*(-HRR(14)+&
                         HRRA(39))+&
                         ABy*(-HRR(17)+&
                         ABz*(-HRR(7)+&
                         HRRA(23))+&
                         HRRA(44))+&
                         HRRA(67)
      !=(13,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(39)+&
                         ABy*(ABz*HRRB(23)+&
                         HRRB(44))+&
                         ABx*(ABz*HRRB(24)+&
                         ABy*(ABz*HRRB(13)+&
                         HRRB(28))+&
                         HRRB(45))+&
                         HRRB(67)
      !=(13_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(28)+&
                         ABz*(-2.D0*HRR(13)+&
                         HRRA(40))+&
                         ABy*(-2.D0*HRR(16)+&
                         ABz*(-2.D0*HRR(6)+&
                         HRRA(24))+&
                         HRRA(45))+&
                         HRRA(68)
      !=(13,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(28)+&
                         ABz*(-HRR(13)+&
                         HRRB(40))+&
                         ABy*(2.D0*ABz*HRRB(24)+&
                         ABy*(ABz*HRRB(13)+&
                         HRRB(28))+&
                         2.D0*HRRB(45))+&
                         HRRB(68)
      !=(13_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(45)+&
                         ABy*(ABz*HRRA(28)+&
                         HRRA(49))+&
                         HRRA(73)
      !=(13,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(24)+&
                         ABz*(ABz*HRRB(24)+&
                         2.D0*HRRB(45))+&
                         ABy*(-HRR(13)+&
                         ABz*(ABz*HRRB(13)+&
                         2.D0*HRRB(28))+&
                         HRRB(49))+&
                         HRRB(73)
      OffSet=(OA+3)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(14_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(40)+&
                         ABy*(ABz*HRRA(24)+&
                         HRRA(45))+&
                         HRRA(68)
      !=(14,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(40)+&
                         ABy*(ABz*HRRB(24)+&
                         HRRB(45))+&
                         ABx*(ABz*HRRB(25)+&
                         ABy*(ABz*HRRB(14)+&
                         HRRB(29))+&
                         HRRB(46))+&
                         HRRB(68)
      !=(14_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-3.D0*HRR(29)+&
                         ABz*(-3.D0*HRR(14)+&
                         HRRA(41))+&
                         ABy*(-3.D0*HRR(17)+&
                         ABz*(-3.D0*HRR(7)+&
                         HRRA(25))+&
                         HRRA(46))+&
                         HRRA(69)
      !=(14,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(29)+&
                         ABz*(-HRR(14)+&
                         HRRB(41))+&
                         ABy*(2.D0*ABz*HRRB(25)+&
                         ABy*(ABz*HRRB(14)+&
                         HRRB(29))+&
                         2.D0*HRRB(46))+&
                         HRRB(69)
      !=(14_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(46)+&
                         ABy*(ABz*HRRA(29)+&
                         HRRA(50))+&
                         HRRA(74)
      !=(14,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(25)+&
                         ABz*(ABz*HRRB(25)+&
                         2.D0*HRRB(46))+&
                         ABy*(-HRR(14)+&
                         ABz*(ABz*HRRB(14)+&
                         2.D0*HRRB(29))+&
                         HRRB(50))+&
                         HRRB(74)
      OffSet=(OA+4)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(15_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(31)+&
                         ABz*(-2.D0*HRR(16)+&
                         HRRA(43))+&
                         ABy*(-2.D0*HRR(18)+&
                         ABz*(-2.D0*HRR(8)+&
                         HRRA(26))+&
                         HRRA(47))+&
                         HRRA(71)
      !=(15,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(43)+&
                         ABy*(ABz*HRRB(26)+&
                         HRRB(47))+&
                         ABx*(ABz*HRRB(27)+&
                         ABy*(ABz*HRRB(15)+&
                         HRRB(30))+&
                         HRRB(48))+&
                         HRRB(71)
      !=(15_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(44)+&
                         ABy*(ABz*HRRA(27)+&
                         HRRA(48))+&
                         HRRA(72)
      !=(15,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(30)+&
                         ABz*(-HRR(15)+&
                         HRRB(44))+&
                         ABy*(2.D0*ABz*HRRB(27)+&
                         ABy*(ABz*HRRB(15)+&
                         HRRB(30))+&
                         2.D0*HRRB(48))+&
                         HRRB(72)
      !=(15_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(27)+&
                         ABz*(-HRR(12)+&
                         HRRA(48))+&
                         ABy*(-HRR(15)+&
                         ABz*(-HRR(5)+&
                         HRRA(30))+&
                         HRRA(51))+&
                         HRRA(76)
      !=(15,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(27)+&
                         ABz*(ABz*HRRB(27)+&
                         2.D0*HRRB(48))+&
                         ABy*(-HRR(15)+&
                         ABz*(ABz*HRRB(15)+&
                         2.D0*HRRB(30))+&
                         HRRB(51))+&
                         HRRB(76)
      OffSet=(OA+5)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(16_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(32)+&
                         ABz*(-HRR(17)+&
                         HRRA(44))+&
                         ABy*(-HRR(19)+&
                         ABz*(-HRR(9)+&
                         HRRA(27))+&
                         HRRA(48))+&
                         HRRA(72)
      !=(16,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(44)+&
                         ABy*(ABz*HRRB(27)+&
                         HRRB(48))+&
                         ABx*(ABz*HRRB(28)+&
                         ABy*(ABz*HRRB(16)+&
                         HRRB(31))+&
                         HRRB(49))+&
                         HRRB(72)
      !=(16_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(31)+&
                         ABz*(-HRR(16)+&
                         HRRA(45))+&
                         ABy*(-HRR(18)+&
                         ABz*(-HRR(8)+&
                         HRRA(28))+&
                         HRRA(49))+&
                         HRRA(73)
      !=(16,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(31)+&
                         ABz*(-HRR(16)+&
                         HRRB(45))+&
                         ABy*(2.D0*ABz*HRRB(28)+&
                         ABy*(ABz*HRRB(16)+&
                         HRRB(31))+&
                         2.D0*HRRB(49))+&
                         HRRB(73)
      !=(16_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(28)+&
                         ABz*(-HRR(13)+&
                         HRRA(49))+&
                         ABy*(-HRR(16)+&
                         ABz*(-HRR(6)+&
                         HRRA(31))+&
                         HRRA(52))+&
                         HRRA(77)
      !=(16,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(28)+&
                         ABz*(ABz*HRRB(28)+&
                         2.D0*HRRB(49))+&
                         ABy*(-HRR(16)+&
                         ABz*(ABz*HRRB(16)+&
                         2.D0*HRRB(31))+&
                         HRRB(52))+&
                         HRRB(77)
      OffSet=(OA+6)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(17_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(45)+&
                         ABy*(ABz*HRRA(28)+&
                         HRRA(49))+&
                         HRRA(73)
      !=(17,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(45)+&
                         ABy*(ABz*HRRB(28)+&
                         HRRB(49))+&
                         ABx*(ABz*HRRB(29)+&
                         ABy*(ABz*HRRB(17)+&
                         HRRB(32))+&
                         HRRB(50))+&
                         HRRB(73)
      !=(17_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(32)+&
                         ABz*(-2.D0*HRR(17)+&
                         HRRA(46))+&
                         ABy*(-2.D0*HRR(19)+&
                         ABz*(-2.D0*HRR(9)+&
                         HRRA(29))+&
                         HRRA(50))+&
                         HRRA(74)
      !=(17,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(32)+&
                         ABz*(-HRR(17)+&
                         HRRB(46))+&
                         ABy*(2.D0*ABz*HRRB(29)+&
                         ABy*(ABz*HRRB(17)+&
                         HRRB(32))+&
                         2.D0*HRRB(50))+&
                         HRRB(74)
      !=(17_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(29)+&
                         ABz*(-HRR(14)+&
                         HRRA(50))+&
                         ABy*(-HRR(17)+&
                         ABz*(-HRR(7)+&
                         HRRA(32))+&
                         HRRA(53))+&
                         HRRA(78)
      !=(17,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(29)+&
                         ABz*(ABz*HRRB(29)+&
                         2.D0*HRRB(50))+&
                         ABy*(-HRR(17)+&
                         ABz*(ABz*HRRB(17)+&
                         2.D0*HRRB(32))+&
                         HRRB(53))+&
                         HRRB(78)
      OffSet=(OA+7)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(18_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(34)+&
                         ABz*(-HRR(19)+&
                         HRRA(48))+&
                         ABy*(-HRR(20)+&
                         ABz*(-HRR(10)+&
                         HRRA(30))+&
                         HRRA(51))+&
                         HRRA(76)
      !=(18,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(48)+&
                         ABy*(ABz*HRRB(30)+&
                         HRRB(51))+&
                         ABx*(ABz*HRRB(31)+&
                         ABy*(ABz*HRRB(18)+&
                         HRRB(33))+&
                         HRRB(52))+&
                         HRRB(76)
      !=(18_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(49)+&
                         ABy*(ABz*HRRA(31)+&
                         HRRA(52))+&
                         HRRA(77)
      !=(18,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(33)+&
                         ABz*(-HRR(18)+&
                         HRRB(49))+&
                         ABy*(2.D0*ABz*HRRB(31)+&
                         ABy*(ABz*HRRB(18)+&
                         HRRB(33))+&
                         2.D0*HRRB(52))+&
                         HRRB(77)
      !=(18_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(31)+&
                         ABz*(-2.D0*HRR(16)+&
                         HRRA(52))+&
                         ABy*(-2.D0*HRR(18)+&
                         ABz*(-2.D0*HRR(8)+&
                         HRRA(33))+&
                         HRRA(54))+&
                         HRRA(80)
      !=(18,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(31)+&
                         ABz*(ABz*HRRB(31)+&
                         2.D0*HRRB(52))+&
                         ABy*(-HRR(18)+&
                         ABz*(ABz*HRRB(18)+&
                         2.D0*HRRB(33))+&
                         HRRB(54))+&
                         HRRB(80)
      OffSet=(OA+8)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(19_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(49)+&
                         ABy*(ABz*HRRA(31)+&
                         HRRA(52))+&
                         HRRA(77)
      !=(19,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(49)+&
                         ABy*(ABz*HRRB(31)+&
                         HRRB(52))+&
                         ABx*(ABz*HRRB(32)+&
                         ABy*(ABz*HRRB(19)+&
                         HRRB(34))+&
                         HRRB(53))+&
                         HRRB(77)
      !=(19_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(34)+&
                         ABz*(-HRR(19)+&
                         HRRA(50))+&
                         ABy*(-HRR(20)+&
                         ABz*(-HRR(10)+&
                         HRRA(32))+&
                         HRRA(53))+&
                         HRRA(78)
      !=(19,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(34)+&
                         ABz*(-HRR(19)+&
                         HRRB(50))+&
                         ABy*(2.D0*ABz*HRRB(32)+&
                         ABy*(ABz*HRRB(19)+&
                         HRRB(34))+&
                         2.D0*HRRB(53))+&
                         HRRB(78)
      !=(19_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(32)+&
                         ABz*(-2.D0*HRR(17)+&
                         HRRA(53))+&
                         ABy*(-2.D0*HRR(19)+&
                         ABz*(-2.D0*HRR(9)+&
                         HRRA(34))+&
                         HRRA(55))+&
                         HRRA(81)
      !=(19,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(32)+&
                         ABz*(ABz*HRRB(32)+&
                         2.D0*HRRB(53))+&
                         ABy*(-HRR(19)+&
                         ABz*(ABz*HRRB(19)+&
                         2.D0*HRRB(34))+&
                         HRRB(55))+&
                         HRRB(81)
      OffSet=(OA+9)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(20_x,9|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(52)+&
                         ABy*(ABz*HRRA(33)+&
                         HRRA(54))+&
                         HRRA(80)
      !=(20,9_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(52)+&
                         ABy*(ABz*HRRB(33)+&
                         HRRB(54))+&
                         ABx*(ABz*HRRB(34)+&
                         ABy*(ABz*HRRB(20)+&
                         HRRB(35))+&
                         HRRB(55))+&
                         HRRB(80)
      !=(20_y,9|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(53)+&
                         ABy*(ABz*HRRA(34)+&
                         HRRA(55))+&
                         HRRA(81)
      !=(20,9_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(35)+&
                         ABz*(-HRR(20)+&
                         HRRB(53))+&
                         ABy*(2.D0*ABz*HRRB(34)+&
                         ABy*(ABz*HRRB(20)+&
                         HRRB(35))+&
                         2.D0*HRRB(55))+&
                         HRRB(81)
      !=(20_z,9|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-3.D0*HRR(34)+&
                         ABz*(-3.D0*HRR(19)+&
                         HRRA(55))+&
                         ABy*(-3.D0*HRR(20)+&
                         ABz*(-3.D0*HRR(10)+&
                         HRRA(35))+&
                         HRRA(56))+&
                         HRRA(83)
      !=(20,9_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(34)+&
                         ABz*(ABz*HRRB(34)+&
                         2.D0*HRRB(55))+&
                         ABy*(-HRR(20)+&
                         ABz*(ABz*HRRB(20)+&
                         2.D0*HRRB(35))+&
                         HRRB(56))+&
                         HRRB(83)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(11_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-3.D0*HRR(30)+&
                         ABz*(-6.D0*HRR(15)+&
                         ABz*(-3.D0*HRR(5)+&
                         HRRA(21))+&
                         2.D0*HRRA(42))+&
                         HRRA(70)
      !=(11,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(21)+&
                         2.D0*HRRB(42))+&
                         ABx*(ABz*(ABz*HRRB(11)+&
                         2.D0*HRRB(26))+&
                         HRRB(47))+&
                         HRRB(70)
      !=(11_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*(ABz*HRRA(22)+&
                         2.D0*HRRA(43))+&
                         HRRA(71)
      !=(11,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(22)+&
                         2.D0*HRRB(43))+&
                         ABy*(ABz*(ABz*HRRB(11)+&
                         2.D0*HRRB(26))+&
                         HRRB(47))+&
                         HRRB(71)
      !=(11_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*(ABz*HRRA(26)+&
                         2.D0*HRRA(47))+&
                         HRRA(75)
      !=(11,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(26)+&
                         ABz*(-2.D0*HRR(11)+&
                         ABz*(ABz*HRRB(11)+&
                         3.D0*HRRB(26))+&
                         3.D0*HRRB(47))+&
                         HRRB(75)
      OffSet=(OA+1)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(12_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(31)+&
                         ABz*(-4.D0*HRR(16)+&
                         ABz*(-2.D0*HRR(6)+&
                         HRRA(22))+&
                         2.D0*HRRA(43))+&
                         HRRA(71)
      !=(12,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(22)+&
                         2.D0*HRRB(43))+&
                         ABx*(ABz*(ABz*HRRB(12)+&
                         2.D0*HRRB(27))+&
                         HRRB(48))+&
                         HRRB(71)
      !=(12_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(30)+&
                         ABz*(-2.D0*HRR(15)+&
                         ABz*(-HRR(5)+&
                         HRRA(23))+&
                         2.D0*HRRA(44))+&
                         HRRA(72)
      !=(12,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(23)+&
                         2.D0*HRRB(44))+&
                         ABy*(ABz*(ABz*HRRB(12)+&
                         2.D0*HRRB(27))+&
                         HRRB(48))+&
                         HRRB(72)
      !=(12_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*(ABz*HRRA(27)+&
                         2.D0*HRRA(48))+&
                         HRRA(76)
      !=(12,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(27)+&
                         ABz*(-2.D0*HRR(12)+&
                         ABz*(ABz*HRRB(12)+&
                         3.D0*HRRB(27))+&
                         3.D0*HRRB(48))+&
                         HRRB(76)
      OffSet=(OA+2)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(13_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(32)+&
                         ABz*(-2.D0*HRR(17)+&
                         ABz*(-HRR(7)+&
                         HRRA(23))+&
                         2.D0*HRRA(44))+&
                         HRRA(72)
      !=(13,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(23)+&
                         2.D0*HRRB(44))+&
                         ABx*(ABz*(ABz*HRRB(13)+&
                         2.D0*HRRB(28))+&
                         HRRB(49))+&
                         HRRB(72)
      !=(13_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(31)+&
                         ABz*(-4.D0*HRR(16)+&
                         ABz*(-2.D0*HRR(6)+&
                         HRRA(24))+&
                         2.D0*HRRA(45))+&
                         HRRA(73)
      !=(13,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(24)+&
                         2.D0*HRRB(45))+&
                         ABy*(ABz*(ABz*HRRB(13)+&
                         2.D0*HRRB(28))+&
                         HRRB(49))+&
                         HRRB(73)
      !=(13_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*(ABz*HRRA(28)+&
                         2.D0*HRRA(49))+&
                         HRRA(77)
      !=(13,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(28)+&
                         ABz*(-2.D0*HRR(13)+&
                         ABz*(ABz*HRRB(13)+&
                         3.D0*HRRB(28))+&
                         3.D0*HRRB(49))+&
                         HRRB(77)
      OffSet=(OA+3)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(14_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*(ABz*HRRA(24)+&
                         2.D0*HRRA(45))+&
                         HRRA(73)
      !=(14,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(24)+&
                         2.D0*HRRB(45))+&
                         ABx*(ABz*(ABz*HRRB(14)+&
                         2.D0*HRRB(29))+&
                         HRRB(50))+&
                         HRRB(73)
      !=(14_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-3.D0*HRR(32)+&
                         ABz*(-6.D0*HRR(17)+&
                         ABz*(-3.D0*HRR(7)+&
                         HRRA(25))+&
                         2.D0*HRRA(46))+&
                         HRRA(74)
      !=(14,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(25)+&
                         2.D0*HRRB(46))+&
                         ABy*(ABz*(ABz*HRRB(14)+&
                         2.D0*HRRB(29))+&
                         HRRB(50))+&
                         HRRB(74)
      !=(14_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*(ABz*HRRA(29)+&
                         2.D0*HRRA(50))+&
                         HRRA(78)
      !=(14,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(29)+&
                         ABz*(-2.D0*HRR(14)+&
                         ABz*(ABz*HRRB(14)+&
                         3.D0*HRRB(29))+&
                         3.D0*HRRB(50))+&
                         HRRB(78)
      OffSet=(OA+4)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(15_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(33)+&
                         ABz*(-4.D0*HRR(18)+&
                         ABz*(-2.D0*HRR(8)+&
                         HRRA(26))+&
                         2.D0*HRRA(47))+&
                         HRRA(75)
      !=(15,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(26)+&
                         2.D0*HRRB(47))+&
                         ABx*(ABz*(ABz*HRRB(15)+&
                         2.D0*HRRB(30))+&
                         HRRB(51))+&
                         HRRB(75)
      !=(15_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*(ABz*HRRA(27)+&
                         2.D0*HRRA(48))+&
                         HRRA(76)
      !=(15,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(27)+&
                         2.D0*HRRB(48))+&
                         ABy*(ABz*(ABz*HRRB(15)+&
                         2.D0*HRRB(30))+&
                         HRRB(51))+&
                         HRRB(76)
      !=(15_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(30)+&
                         ABz*(-2.D0*HRR(15)+&
                         ABz*(-HRR(5)+&
                         HRRA(30))+&
                         2.D0*HRRA(51))+&
                         HRRA(79)
      !=(15,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(30)+&
                         ABz*(-2.D0*HRR(15)+&
                         ABz*(ABz*HRRB(15)+&
                         3.D0*HRRB(30))+&
                         3.D0*HRRB(51))+&
                         HRRB(79)
      OffSet=(OA+5)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(16_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(34)+&
                         ABz*(-2.D0*HRR(19)+&
                         ABz*(-HRR(9)+&
                         HRRA(27))+&
                         2.D0*HRRA(48))+&
                         HRRA(76)
      !=(16,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(27)+&
                         2.D0*HRRB(48))+&
                         ABx*(ABz*(ABz*HRRB(16)+&
                         2.D0*HRRB(31))+&
                         HRRB(52))+&
                         HRRB(76)
      !=(16_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(33)+&
                         ABz*(-2.D0*HRR(18)+&
                         ABz*(-HRR(8)+&
                         HRRA(28))+&
                         2.D0*HRRA(49))+&
                         HRRA(77)
      !=(16,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(28)+&
                         2.D0*HRRB(49))+&
                         ABy*(ABz*(ABz*HRRB(16)+&
                         2.D0*HRRB(31))+&
                         HRRB(52))+&
                         HRRB(77)
      !=(16_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(31)+&
                         ABz*(-2.D0*HRR(16)+&
                         ABz*(-HRR(6)+&
                         HRRA(31))+&
                         2.D0*HRRA(52))+&
                         HRRA(80)
      !=(16,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(31)+&
                         ABz*(-2.D0*HRR(16)+&
                         ABz*(ABz*HRRB(16)+&
                         3.D0*HRRB(31))+&
                         3.D0*HRRB(52))+&
                         HRRB(80)
      OffSet=(OA+6)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(17_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*(ABz*HRRA(28)+&
                         2.D0*HRRA(49))+&
                         HRRA(77)
      !=(17,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(28)+&
                         2.D0*HRRB(49))+&
                         ABx*(ABz*(ABz*HRRB(17)+&
                         2.D0*HRRB(32))+&
                         HRRB(53))+&
                         HRRB(77)
      !=(17_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(34)+&
                         ABz*(-4.D0*HRR(19)+&
                         ABz*(-2.D0*HRR(9)+&
                         HRRA(29))+&
                         2.D0*HRRA(50))+&
                         HRRA(78)
      !=(17,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(29)+&
                         2.D0*HRRB(50))+&
                         ABy*(ABz*(ABz*HRRB(17)+&
                         2.D0*HRRB(32))+&
                         HRRB(53))+&
                         HRRB(78)
      !=(17_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(32)+&
                         ABz*(-2.D0*HRR(17)+&
                         ABz*(-HRR(7)+&
                         HRRA(32))+&
                         2.D0*HRRA(53))+&
                         HRRA(81)
      !=(17,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(32)+&
                         ABz*(-2.D0*HRR(17)+&
                         ABz*(ABz*HRRB(17)+&
                         3.D0*HRRB(32))+&
                         3.D0*HRRB(53))+&
                         HRRB(81)
      OffSet=(OA+7)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(18_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(35)+&
                         ABz*(-2.D0*HRR(20)+&
                         ABz*(-HRR(10)+&
                         HRRA(30))+&
                         2.D0*HRRA(51))+&
                         HRRA(79)
      !=(18,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(30)+&
                         2.D0*HRRB(51))+&
                         ABx*(ABz*(ABz*HRRB(18)+&
                         2.D0*HRRB(33))+&
                         HRRB(54))+&
                         HRRB(79)
      !=(18_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*(ABz*HRRA(31)+&
                         2.D0*HRRA(52))+&
                         HRRA(80)
      !=(18,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(31)+&
                         2.D0*HRRB(52))+&
                         ABy*(ABz*(ABz*HRRB(18)+&
                         2.D0*HRRB(33))+&
                         HRRB(54))+&
                         HRRB(80)
      !=(18_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(33)+&
                         ABz*(-4.D0*HRR(18)+&
                         ABz*(-2.D0*HRR(8)+&
                         HRRA(33))+&
                         2.D0*HRRA(54))+&
                         HRRA(82)
      !=(18,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(33)+&
                         ABz*(-2.D0*HRR(18)+&
                         ABz*(ABz*HRRB(18)+&
                         3.D0*HRRB(33))+&
                         3.D0*HRRB(54))+&
                         HRRB(82)
      OffSet=(OA+8)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(19_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*(ABz*HRRA(31)+&
                         2.D0*HRRA(52))+&
                         HRRA(80)
      !=(19,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(31)+&
                         2.D0*HRRB(52))+&
                         ABx*(ABz*(ABz*HRRB(19)+&
                         2.D0*HRRB(34))+&
                         HRRB(55))+&
                         HRRB(80)
      !=(19_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(35)+&
                         ABz*(-2.D0*HRR(20)+&
                         ABz*(-HRR(10)+&
                         HRRA(32))+&
                         2.D0*HRRA(53))+&
                         HRRA(81)
      !=(19,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(32)+&
                         2.D0*HRRB(53))+&
                         ABy*(ABz*(ABz*HRRB(19)+&
                         2.D0*HRRB(34))+&
                         HRRB(55))+&
                         HRRB(81)
      !=(19_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(34)+&
                         ABz*(-4.D0*HRR(19)+&
                         ABz*(-2.D0*HRR(9)+&
                         HRRA(34))+&
                         2.D0*HRRA(55))+&
                         HRRA(83)
      !=(19,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(34)+&
                         ABz*(-2.D0*HRR(19)+&
                         ABz*(ABz*HRRB(19)+&
                         3.D0*HRRB(34))+&
                         3.D0*HRRB(55))+&
                         HRRB(83)
      OffSet=(OA+9)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(20_x,10|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*(ABz*HRRA(33)+&
                         2.D0*HRRA(54))+&
                         HRRA(82)
      !=(20,10_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*(ABz*HRRB(33)+&
                         2.D0*HRRB(54))+&
                         ABx*(ABz*(ABz*HRRB(20)+&
                         2.D0*HRRB(35))+&
                         HRRB(56))+&
                         HRRB(82)
      !=(20_y,10|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*(ABz*HRRA(34)+&
                         2.D0*HRRA(55))+&
                         HRRA(83)
      !=(20,10_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*(ABz*HRRB(34)+&
                         2.D0*HRRB(55))+&
                         ABy*(ABz*(ABz*HRRB(20)+&
                         2.D0*HRRB(35))+&
                         HRRB(56))+&
                         HRRB(83)
      !=(20_z,10|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-3.D0*HRR(35)+&
                         ABz*(-6.D0*HRR(20)+&
                         ABz*(-3.D0*HRR(10)+&
                         HRRA(35))+&
                         2.D0*HRRA(56))+&
                         HRRA(84)
      !=(20,10_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-2.D0*HRR(35)+&
                         ABz*(-2.D0*HRR(20)+&
                         ABz*(ABz*HRRB(20)+&
                         3.D0*HRRB(35))+&
                         3.D0*HRRB(56))+&
                         HRRB(84)
    END SUBROUTINE BraHRR106ab
    SUBROUTINE BraHRR106cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,Cart,HRR,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,Cart,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(11,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(11)+&
                                2.D0*HRR(21))+&
                                HRR(36)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(12,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(12)+&
                                2.D0*HRR(22))+&
                                HRR(37)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(13,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(13)+&
                                2.D0*HRR(23))+&
                                HRR(38)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(14,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(14)+&
                                2.D0*HRR(24))+&
                                HRR(39)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(15,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(15)+&
                                2.D0*HRR(26))+&
                                HRR(42)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(16,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(16)+&
                                2.D0*HRR(27))+&
                                HRR(43)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+6)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(17,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(17)+&
                                2.D0*HRR(28))+&
                                HRR(44)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+7)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(18,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(18)+&
                                2.D0*HRR(30))+&
                                HRR(47)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+8)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(19,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(19)+&
                                2.D0*HRR(31))+&
                                HRR(48)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+9)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(20,5|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*(ABx*HRR(20)+&
                                2.D0*HRR(33))+&
                                HRR(51)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(11,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(21)+&
                                ABx*(ABy*HRR(11)+&
                                HRR(22))+&
                                HRR(37)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(12,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(22)+&
                                ABx*(ABy*HRR(12)+&
                                HRR(23))+&
                                HRR(38)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(13,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(23)+&
                                ABx*(ABy*HRR(13)+&
                                HRR(24))+&
                                HRR(39)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(14,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(24)+&
                                ABx*(ABy*HRR(14)+&
                                HRR(25))+&
                                HRR(40)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+4)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(15,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(26)+&
                                ABx*(ABy*HRR(15)+&
                                HRR(27))+&
                                HRR(43)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+5)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(16,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(27)+&
                                ABx*(ABy*HRR(16)+&
                                HRR(28))+&
                                HRR(44)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+6)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(17,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(28)+&
                                ABx*(ABy*HRR(17)+&
                                HRR(29))+&
                                HRR(45)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+7)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(18,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(30)+&
                                ABx*(ABy*HRR(18)+&
                                HRR(31))+&
                                HRR(48)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+8)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(19,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(31)+&
                                ABx*(ABy*HRR(19)+&
                                HRR(32))+&
                                HRR(49)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+9)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(20,6|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(33)+&
                                ABx*(ABy*HRR(20)+&
                                HRR(34))+&
                                HRR(52)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(11,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(11)+&
                                2.D0*HRR(22))+&
                                HRR(38)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(12,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(12)+&
                                2.D0*HRR(23))+&
                                HRR(39)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(13,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(13)+&
                                2.D0*HRR(24))+&
                                HRR(40)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(14,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(14)+&
                                2.D0*HRR(25))+&
                                HRR(41)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+4)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(15,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(15)+&
                                2.D0*HRR(27))+&
                                HRR(44)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+5)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(16,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(16)+&
                                2.D0*HRR(28))+&
                                HRR(45)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+6)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(17,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(17)+&
                                2.D0*HRR(29))+&
                                HRR(46)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+7)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(18,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(18)+&
                                2.D0*HRR(31))+&
                                HRR(49)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+8)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(19,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(19)+&
                                2.D0*HRR(32))+&
                                HRR(50)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+9)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(20,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(20)+&
                                2.D0*HRR(34))+&
                                HRR(53)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(11,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(21)+&
                                ABx*(ABz*HRR(11)+&
                                HRR(26))+&
                                HRR(42)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(12,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(22)+&
                                ABx*(ABz*HRR(12)+&
                                HRR(27))+&
                                HRR(43)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(13,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(23)+&
                                ABx*(ABz*HRR(13)+&
                                HRR(28))+&
                                HRR(44)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(14,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(24)+&
                                ABx*(ABz*HRR(14)+&
                                HRR(29))+&
                                HRR(45)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+4)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(15,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(26)+&
                                ABx*(ABz*HRR(15)+&
                                HRR(30))+&
                                HRR(47)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+5)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(16,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(27)+&
                                ABx*(ABz*HRR(16)+&
                                HRR(31))+&
                                HRR(48)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+6)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(17,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(28)+&
                                ABx*(ABz*HRR(17)+&
                                HRR(32))+&
                                HRR(49)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+7)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(18,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(30)+&
                                ABx*(ABz*HRR(18)+&
                                HRR(33))+&
                                HRR(51)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+8)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(19,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(31)+&
                                ABx*(ABz*HRR(19)+&
                                HRR(34))+&
                                HRR(52)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+9)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(20,8|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(33)+&
                                ABx*(ABz*HRR(20)+&
                                HRR(35))+&
                                HRR(54)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(11,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(22)+&
                                ABy*(ABz*HRR(11)+&
                                HRR(26))+&
                                HRR(43)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(12,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(23)+&
                                ABy*(ABz*HRR(12)+&
                                HRR(27))+&
                                HRR(44)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(13,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(24)+&
                                ABy*(ABz*HRR(13)+&
                                HRR(28))+&
                                HRR(45)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(14,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(25)+&
                                ABy*(ABz*HRR(14)+&
                                HRR(29))+&
                                HRR(46)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+4)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(15,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(27)+&
                                ABy*(ABz*HRR(15)+&
                                HRR(30))+&
                                HRR(48)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+5)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(16,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(28)+&
                                ABy*(ABz*HRR(16)+&
                                HRR(31))+&
                                HRR(49)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+6)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(17,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(29)+&
                                ABy*(ABz*HRR(17)+&
                                HRR(32))+&
                                HRR(50)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+7)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(18,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(31)+&
                                ABy*(ABz*HRR(18)+&
                                HRR(33))+&
                                HRR(52)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+8)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(19,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(32)+&
                                ABy*(ABz*HRR(19)+&
                                HRR(34))+&
                                HRR(53)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+9)*LDA+(OB+4)*LDB+CDOffSet !=
      !=(20,9|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(34)+&
                                ABy*(ABz*HRR(20)+&
                                HRR(35))+&
                                HRR(55)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(11,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(11)+&
                                2.D0*HRR(26))+&
                                HRR(47)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(12,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(12)+&
                                2.D0*HRR(27))+&
                                HRR(48)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(13,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(13)+&
                                2.D0*HRR(28))+&
                                HRR(49)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(14,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(14)+&
                                2.D0*HRR(29))+&
                                HRR(50)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+4)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(15,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(15)+&
                                2.D0*HRR(30))+&
                                HRR(51)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+5)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(16,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(16)+&
                                2.D0*HRR(31))+&
                                HRR(52)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+6)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(17,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(17)+&
                                2.D0*HRR(32))+&
                                HRR(53)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+7)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(18,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(18)+&
                                2.D0*HRR(33))+&
                                HRR(54)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+8)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(19,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(19)+&
                                2.D0*HRR(34))+&
                                HRR(55)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+9)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(20,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(20)+&
                                2.D0*HRR(35))+&
                                HRR(56)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
    END SUBROUTINE BraHRR106cd
