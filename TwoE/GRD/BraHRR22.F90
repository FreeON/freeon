    SUBROUTINE BraHRR22ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR,HRRA,HRRB,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*),HRRA(*),HRRB(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(1_x,1|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         HRRA(2)
      !=(1,1_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABx*HRRB(1)+&
                         HRRB(2)
      !=(1_y,1|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         HRRA(3)
      !=(1,1_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(1)+&
                         HRRB(3)
      !=(1_z,1|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         HRRA(4)
      !=(1,1_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(1)+&
                         HRRB(4)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(2_x,1|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(1)+&
                         HRRA(5)
      !=(2,1_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABx*HRRB(2)+&
                         HRRB(5)
      !=(2_y,1|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         HRRA(6)
      !=(2,1_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(2)+&
                         HRRB(6)
      !=(2_z,1|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         HRRA(8)
      !=(2,1_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(2)+&
                         HRRB(8)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(3_x,1|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         HRRA(6)
      !=(3,1_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABx*HRRB(3)+&
                         HRRB(6)
      !=(3_y,1|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(1)+&
                         HRRA(7)
      !=(3,1_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(3)+&
                         HRRB(7)
      !=(3_z,1|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         HRRA(9)
      !=(3,1_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(3)+&
                         HRRB(9)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(4_x,1|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         HRRA(8)
      !=(4,1_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABx*HRRB(4)+&
                         HRRB(8)
      !=(4_y,1|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         HRRA(9)
      !=(4,1_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(4)+&
                         HRRB(9)
      !=(4_z,1|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(1)+&
                         HRRA(10)
      !=(4,1_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(4)+&
                         HRRB(10)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(1_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*HRRA(2)+&
                         HRRA(5)
      !=(1,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(1)+&
                         ABx*(ABx*HRRB(1)+&
                         2.D0*HRRB(2))+&
                         HRRB(5)
      !=(1_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*HRRA(3)+&
                         HRRA(6)
      !=(1,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(2)+&
                         ABx*(ABy*HRRB(1)+&
                         HRRB(3))+&
                         HRRB(6)
      !=(1_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*HRRA(4)+&
                         HRRA(8)
      !=(1,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(2)+&
                         ABx*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         HRRB(8)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(1_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(2)+&
                         HRRA(6)
      !=(1,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(2)+&
                         ABx*(ABy*HRRB(1)+&
                         HRRB(3))+&
                         HRRB(6)
      !=(1_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(3)+&
                         HRRA(7)
      !=(1,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(1)+&
                         ABy*(ABy*HRRB(1)+&
                         2.D0*HRRB(3))+&
                         HRRB(7)
      !=(1_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(4)+&
                         HRRA(9)
      !=(1,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(3)+&
                         ABy*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         HRRB(9)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(1_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(2)+&
                         HRRA(8)
      !=(1,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(2)+&
                         ABx*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         HRRB(8)
      !=(1_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(3)+&
                         HRRA(9)
      !=(1,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(3)+&
                         ABy*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         HRRB(9)
      !=(1_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(4)+&
                         HRRA(10)
      !=(1,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(1)+&
                         ABz*(ABz*HRRB(1)+&
                         2.D0*HRRB(4))+&
                         HRRB(10)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(2_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(2)+&
                         ABx*(-HRR(1)+&
                         HRRA(5))+&
                         HRRA(11)
      !=(2,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(2)+&
                         ABx*(ABx*HRRB(2)+&
                         2.D0*HRRB(5))+&
                         HRRB(11)
      !=(2_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*HRRA(6)+&
                         HRRA(12)
      !=(2,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(5)+&
                         ABx*(ABy*HRRB(2)+&
                         HRRB(6))+&
                         HRRB(12)
      !=(2_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*HRRA(8)+&
                         HRRA(15)
      !=(2,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(5)+&
                         ABx*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         HRRB(15)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(3_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*HRRA(6)+&
                         HRRA(12)
      !=(3,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(3)+&
                         ABx*(ABx*HRRB(3)+&
                         2.D0*HRRB(6))+&
                         HRRB(12)
      !=(3_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(2)+&
                         ABx*(-HRR(1)+&
                         HRRA(7))+&
                         HRRA(13)
      !=(3,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(6)+&
                         ABx*(ABy*HRRB(3)+&
                         HRRB(7))+&
                         HRRB(13)
      !=(3_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABx*HRRA(9)+&
                         HRRA(16)
      !=(3,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(6)+&
                         ABx*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         HRRB(16)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(4_x,2|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABx*HRRA(8)+&
                         HRRA(15)
      !=(4,2_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)-HRR(4)+&
                         ABx*(ABx*HRRB(4)+&
                         2.D0*HRRB(8))+&
                         HRRB(15)
      !=(4_y,2|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABx*HRRA(9)+&
                         HRRA(16)
      !=(4,2_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(8)+&
                         ABx*(ABy*HRRB(4)+&
                         HRRB(9))+&
                         HRRB(16)
      !=(4_z,2|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(2)+&
                         ABx*(-HRR(1)+&
                         HRRA(10))+&
                         HRRA(18)
      !=(4,2_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(8)+&
                         ABx*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         HRRB(18)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(2_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(5))+&
                         HRRA(12)
      !=(2,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(5)+&
                         ABx*(ABy*HRRB(2)+&
                         HRRB(6))+&
                         HRRB(12)
      !=(2_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(6)+&
                         HRRA(13)
      !=(2,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(2)+&
                         ABy*(ABy*HRRB(2)+&
                         2.D0*HRRB(6))+&
                         HRRB(13)
      !=(2_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(8)+&
                         HRRA(16)
      !=(2,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         HRRB(16)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(3_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(6)+&
                         HRRA(13)
      !=(3,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(6)+&
                         ABx*(ABy*HRRB(3)+&
                         HRRB(7))+&
                         HRRB(13)
      !=(3_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(7))+&
                         HRRA(14)
      !=(3,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(3)+&
                         ABy*(ABy*HRRB(3)+&
                         2.D0*HRRB(7))+&
                         HRRB(14)
      !=(3_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABy*HRRA(9)+&
                         HRRA(17)
      !=(3,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(7)+&
                         ABy*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         HRRB(17)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(4_x,3|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABy*HRRA(8)+&
                         HRRA(16)
      !=(4,3_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABy*HRRB(8)+&
                         ABx*(ABy*HRRB(4)+&
                         HRRB(9))+&
                         HRRB(16)
      !=(4_y,3|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABy*HRRA(9)+&
                         HRRA(17)
      !=(4,3_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)-HRR(4)+&
                         ABy*(ABy*HRRB(4)+&
                         2.D0*HRRB(9))+&
                         HRRB(17)
      !=(4_z,3|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(10))+&
                         HRRA(19)
      !=(4,3_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(9)+&
                         ABy*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         HRRB(19)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(2_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(5))+&
                         HRRA(15)
      !=(2,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(5)+&
                         ABx*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         HRRB(15)
      !=(2_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(6)+&
                         HRRA(16)
      !=(2,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         HRRB(16)
      !=(2_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(8)+&
                         HRRA(18)
      !=(2,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(2)+&
                         ABz*(ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         HRRB(18)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(3_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(6)+&
                         HRRA(16)
      !=(3,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(6)+&
                         ABx*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         HRRB(16)
      !=(3_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(7))+&
                         HRRA(17)
      !=(3,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(7)+&
                         ABy*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         HRRB(17)
      !=(3_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         ABz*HRRA(9)+&
                         HRRA(19)
      !=(3,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(3)+&
                         ABz*(ABz*HRRB(3)+&
                         2.D0*HRRB(9))+&
                         HRRB(19)
      OffSet=(OA+3)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(4_x,4|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         ABz*HRRA(8)+&
                         HRRA(18)
      !=(4,4_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABz*HRRB(8)+&
                         ABx*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         HRRB(18)
      !=(4_y,4|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         ABz*HRRA(9)+&
                         HRRA(19)
      !=(4,4_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABz*HRRB(9)+&
                         ABy*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         HRRB(19)
      !=(4_z,4|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(10))+&
                         HRRA(20)
      !=(4,4_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)-HRR(4)+&
                         ABz*(ABz*HRRB(4)+&
                         2.D0*HRRB(10))+&
                         HRRB(20)
    END SUBROUTINE BraHRR22ab
    SUBROUTINE BraHRR22cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,Cart,HRR,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,Cart,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(1,1|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                HRR(1)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(2,1|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                HRR(2)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(3,1|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                HRR(3)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(4,1|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                HRR(4)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(1,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(1)+&
                                HRR(2)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(1,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(1)+&
                                HRR(3)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(1,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(1)+&
                                HRR(4)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(2,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(2)+&
                                HRR(5)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(3,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(3)+&
                                HRR(6)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(4,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(4)+&
                                HRR(8)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(2,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(2)+&
                                HRR(6)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(3,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(3)+&
                                HRR(7)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(4,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(4)+&
                                HRR(9)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(2,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(2)+&
                                HRR(8)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(3,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(3)+&
                                HRR(9)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(4,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(4)+&
                                HRR(10)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
    END SUBROUTINE BraHRR22cd
