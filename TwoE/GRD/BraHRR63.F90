    SUBROUTINE BraHRR63ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR,HRRA,HRRB,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*),HRRA(*),HRRB(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(5_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(5)+&
                         ABx*(-2.D0*HRR(2)+&
                         HRRA(11))+&
                         HRRA(21)
      !=(5,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(5)+&
                         ABx*(ABx*HRRB(5)+&
                         2.D0*HRRB(11))+&
                         HRRB(21)
      !=(5_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*HRRA(12)+&
                         HRRA(22)
      !=(5,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(11)+&
                         ABx*(ABy*HRRB(5)+&
                         HRRB(12))+&
                         HRRB(22)
      !=(5_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*HRRA(15)+&
                         HRRA(26)
      !=(5,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(11)+&
                         ABx*(ABz*HRRB(5)+&
                         HRRB(15))+&
                         HRRB(26)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(6_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(6)+&
                         ABx*(-HRR(3)+&
                         HRRA(12))+&
                         HRRA(22)
      !=(6,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(6)+&
                         ABx*(ABx*HRRB(6)+&
                         2.D0*HRRB(12))+&
                         HRRB(22)
      !=(6_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(5)+&
                         ABx*(-HRR(2)+&
                         HRRA(13))+&
                         HRRA(23)
      !=(6,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(12)+&
                         ABx*(ABy*HRRB(6)+&
                         HRRB(13))+&
                         HRRB(23)
      !=(6_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*HRRA(16)+&
                         HRRA(27)
      !=(6,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(12)+&
                         ABx*(ABz*HRRB(6)+&
                         HRRB(16))+&
                         HRRB(27)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(7_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*HRRA(13)+&
                         HRRA(23)
      !=(7,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(7)+&
                         ABx*(ABx*HRRB(7)+&
                         2.D0*HRRB(13))+&
                         HRRB(23)
      !=(7_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(6)+&
                         ABx*(-2.D0*HRR(3)+&
                         HRRA(14))+&
                         HRRA(24)
      !=(7,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(13)+&
                         ABx*(ABy*HRRB(7)+&
                         HRRB(14))+&
                         HRRB(24)
      !=(7_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*HRRA(17)+&
                         HRRA(28)
      !=(7,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(13)+&
                         ABx*(ABz*HRRB(7)+&
                         HRRB(17))+&
                         HRRB(28)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(8_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(8)+&
                         ABx*(-HRR(4)+&
                         HRRA(15))+&
                         HRRA(26)
      !=(8,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(8)+&
                         ABx*(ABx*HRRB(8)+&
                         2.D0*HRRB(15))+&
                         HRRB(26)
      !=(8_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*HRRA(16)+&
                         HRRA(27)
      !=(8,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(15)+&
                         ABx*(ABy*HRRB(8)+&
                         HRRB(16))+&
                         HRRB(27)
      !=(8_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(5)+&
                         ABx*(-HRR(2)+&
                         HRRA(18))+&
                         HRRA(30)
      !=(8,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(15)+&
                         ABx*(ABz*HRRB(8)+&
                         HRRB(18))+&
                         HRRB(30)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(9_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*HRRA(16)+&
                         HRRA(27)
      !=(9,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(9)+&
                         ABx*(ABx*HRRB(9)+&
                         2.D0*HRRB(16))+&
                         HRRB(27)
      !=(9_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(8)+&
                         ABx*(-HRR(4)+&
                         HRRA(17))+&
                         HRRA(28)
      !=(9,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(16)+&
                         ABx*(ABy*HRRB(9)+&
                         HRRB(17))+&
                         HRRB(28)
      !=(9_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(6)+&
                         ABx*(-HRR(3)+&
                         HRRA(19))+&
                         HRRA(31)
      !=(9,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(16)+&
                         ABx*(ABz*HRRB(9)+&
                         HRRB(19))+&
                         HRRB(31)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(10_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*HRRA(18)+&
                         HRRA(30)
      !=(10,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(10)+&
                         ABx*(ABx*HRRB(10)+&
                         2.D0*HRRB(18))+&
                         HRRB(30)
      !=(10_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*HRRA(19)+&
                         HRRA(31)
      !=(10,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(18)+&
                         ABx*(ABy*HRRB(10)+&
                         HRRB(19))+&
                         HRRB(31)
      !=(10_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(8)+&
                         ABx*(-2.D0*HRR(4)+&
                         HRRA(20))+&
                         HRRA(33)
      !=(10,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(18)+&
                         ABx*(ABz*HRRB(10)+&
                         HRRB(20))+&
                         HRRB(33)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(5_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(6)+&
                         ABy*(-2.D0*HRR(2)+&
                         HRRA(11))+&
                         HRRA(22)
      !=(5,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(11)+&
                         ABx*(ABy*HRRB(5)+&
                         HRRB(12))+&
                         HRRB(22)
      !=(5_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(12)+&
                         HRRA(23)
      !=(5,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(5)+&
                         ABy*(ABy*HRRB(5)+&
                         2.D0*HRRB(12))+&
                         HRRB(23)
      !=(5_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(15)+&
                         HRRA(27)
      !=(5,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(12)+&
                         ABy*(ABz*HRRB(5)+&
                         HRRB(15))+&
                         HRRB(27)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(6_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(7)+&
                         ABy*(-HRR(3)+&
                         HRRA(12))+&
                         HRRA(23)
      !=(6,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(12)+&
                         ABx*(ABy*HRRB(6)+&
                         HRRB(13))+&
                         HRRB(23)
      !=(6_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(6)+&
                         ABy*(-HRR(2)+&
                         HRRA(13))+&
                         HRRA(24)
      !=(6,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(6)+&
                         ABy*(ABy*HRRB(6)+&
                         2.D0*HRRB(13))+&
                         HRRB(24)
      !=(6_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(16)+&
                         HRRA(28)
      !=(6,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(13)+&
                         ABy*(ABz*HRRB(6)+&
                         HRRB(16))+&
                         HRRB(28)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(7_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(13)+&
                         HRRA(24)
      !=(7,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(13)+&
                         ABx*(ABy*HRRB(7)+&
                         HRRB(14))+&
                         HRRB(24)
      !=(7_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(7)+&
                         ABy*(-2.D0*HRR(3)+&
                         HRRA(14))+&
                         HRRA(25)
      !=(7,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(7)+&
                         ABy*(ABy*HRRB(7)+&
                         2.D0*HRRB(14))+&
                         HRRB(25)
      !=(7_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(17)+&
                         HRRA(29)
      !=(7,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(14)+&
                         ABy*(ABz*HRRB(7)+&
                         HRRB(17))+&
                         HRRB(29)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(8_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(9)+&
                         ABy*(-HRR(4)+&
                         HRRA(15))+&
                         HRRA(27)
      !=(8,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(15)+&
                         ABx*(ABy*HRRB(8)+&
                         HRRB(16))+&
                         HRRB(27)
      !=(8_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(16)+&
                         HRRA(28)
      !=(8,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(8)+&
                         ABy*(ABy*HRRB(8)+&
                         2.D0*HRRB(16))+&
                         HRRB(28)
      !=(8_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(6)+&
                         ABy*(-HRR(2)+&
                         HRRA(18))+&
                         HRRA(31)
      !=(8,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(16)+&
                         ABy*(ABz*HRRB(8)+&
                         HRRB(18))+&
                         HRRB(31)
      OffSet=(OA+4)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(9_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(16)+&
                         HRRA(28)
      !=(9,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(16)+&
                         ABx*(ABy*HRRB(9)+&
                         HRRB(17))+&
                         HRRB(28)
      !=(9_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(9)+&
                         ABy*(-HRR(4)+&
                         HRRA(17))+&
                         HRRA(29)
      !=(9,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(9)+&
                         ABy*(ABy*HRRB(9)+&
                         2.D0*HRRB(17))+&
                         HRRB(29)
      !=(9_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(7)+&
                         ABy*(-HRR(3)+&
                         HRRA(19))+&
                         HRRA(32)
      !=(9,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(17)+&
                         ABy*(ABz*HRRB(9)+&
                         HRRB(19))+&
                         HRRB(32)
      OffSet=(OA+5)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(10_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(18)+&
                         HRRA(31)
      !=(10,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(18)+&
                         ABx*(ABy*HRRB(10)+&
                         HRRB(19))+&
                         HRRB(31)
      !=(10_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(19)+&
                         HRRA(32)
      !=(10,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(10)+&
                         ABy*(ABy*HRRB(10)+&
                         2.D0*HRRB(19))+&
                         HRRB(32)
      !=(10_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(9)+&
                         ABy*(-2.D0*HRR(4)+&
                         HRRA(20))+&
                         HRRA(34)
      !=(10,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(19)+&
                         ABy*(ABz*HRRB(10)+&
                         HRRB(20))+&
                         HRRB(34)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(5_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(8)+&
                         ABz*(-2.D0*HRR(2)+&
                         HRRA(11))+&
                         HRRA(26)
      !=(5,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(11)+&
                         ABx*(ABz*HRRB(5)+&
                         HRRB(15))+&
                         HRRB(26)
      !=(5_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(12)+&
                         HRRA(27)
      !=(5,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(12)+&
                         ABy*(ABz*HRRB(5)+&
                         HRRB(15))+&
                         HRRB(27)
      !=(5_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(15)+&
                         HRRA(30)
      !=(5,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(5)+&
                         ABz*(ABz*HRRB(5)+&
                         2.D0*HRRB(15))+&
                         HRRB(30)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(6_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRA(12))+&
                         HRRA(27)
      !=(6,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(12)+&
                         ABx*(ABz*HRRB(6)+&
                         HRRB(16))+&
                         HRRB(27)
      !=(6_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRA(13))+&
                         HRRA(28)
      !=(6,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(13)+&
                         ABy*(ABz*HRRB(6)+&
                         HRRB(16))+&
                         HRRB(28)
      !=(6_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(16)+&
                         HRRA(31)
      !=(6,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(6)+&
                         ABz*(ABz*HRRB(6)+&
                         2.D0*HRRB(16))+&
                         HRRB(31)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(7_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(13)+&
                         HRRA(28)
      !=(7,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(13)+&
                         ABx*(ABz*HRRB(7)+&
                         HRRB(17))+&
                         HRRB(28)
      !=(7_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(9)+&
                         ABz*(-2.D0*HRR(3)+&
                         HRRA(14))+&
                         HRRA(29)
      !=(7,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(14)+&
                         ABy*(ABz*HRRB(7)+&
                         HRRB(17))+&
                         HRRB(29)
      !=(7_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(17)+&
                         HRRA(32)
      !=(7,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(7)+&
                         ABz*(ABz*HRRB(7)+&
                         2.D0*HRRB(17))+&
                         HRRB(32)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(8_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(10)+&
                         ABz*(-HRR(4)+&
                         HRRA(15))+&
                         HRRA(30)
      !=(8,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(15)+&
                         ABx*(ABz*HRRB(8)+&
                         HRRB(18))+&
                         HRRB(30)
      !=(8_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(16)+&
                         HRRA(31)
      !=(8,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(16)+&
                         ABy*(ABz*HRRB(8)+&
                         HRRB(18))+&
                         HRRB(31)
      !=(8_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(8)+&
                         ABz*(-HRR(2)+&
                         HRRA(18))+&
                         HRRA(33)
      !=(8,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(8)+&
                         ABz*(ABz*HRRB(8)+&
                         2.D0*HRRB(18))+&
                         HRRB(33)
      OffSet=(OA+4)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(9_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(16)+&
                         HRRA(31)
      !=(9,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(16)+&
                         ABx*(ABz*HRRB(9)+&
                         HRRB(19))+&
                         HRRB(31)
      !=(9_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(10)+&
                         ABz*(-HRR(4)+&
                         HRRA(17))+&
                         HRRA(32)
      !=(9,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(17)+&
                         ABy*(ABz*HRRB(9)+&
                         HRRB(19))+&
                         HRRB(32)
      !=(9_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(9)+&
                         ABz*(-HRR(3)+&
                         HRRA(19))+&
                         HRRA(34)
      !=(9,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(9)+&
                         ABz*(ABz*HRRB(9)+&
                         2.D0*HRRB(19))+&
                         HRRB(34)
      OffSet=(OA+5)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(10_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(18)+&
                         HRRA(33)
      !=(10,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(18)+&
                         ABx*(ABz*HRRB(10)+&
                         HRRB(20))+&
                         HRRB(33)
      !=(10_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(19)+&
                         HRRA(34)
      !=(10,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(19)+&
                         ABy*(ABz*HRRB(10)+&
                         HRRB(20))+&
                         HRRB(34)
      !=(10_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(10)+&
                         ABz*(-2.D0*HRR(4)+&
                         HRRA(20))+&
                         HRRA(35)
      !=(10,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(10)+&
                         ABz*(ABz*HRRB(10)+&
                         2.D0*HRRB(20))+&
                         HRRB(35)
    END SUBROUTINE BraHRR63ab
    SUBROUTINE BraHRR63cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,Cart,HRR,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,Cart,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(5,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(5)+&
                                HRR(11)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(6,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(6)+&
                                HRR(12)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(7,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(7)+&
                                HRR(13)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(8,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(8)+&
                                HRR(15)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(9,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(9)+&
                                HRR(16)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(10,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(10)+&
                                HRR(18)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(5,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(5)+&
                                HRR(12)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(6,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(6)+&
                                HRR(13)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(7,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(7)+&
                                HRR(14)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(8,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(8)+&
                                HRR(16)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+4)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(9,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(9)+&
                                HRR(17)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+5)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(10,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(10)+&
                                HRR(19)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(5,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(5)+&
                                HRR(15)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(6,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(6)+&
                                HRR(16)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(7,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(7)+&
                                HRR(17)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(8,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(8)+&
                                HRR(18)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+4)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(9,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(9)+&
                                HRR(19)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+5)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(10,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(10)+&
                                HRR(20)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
    END SUBROUTINE BraHRR63cd
