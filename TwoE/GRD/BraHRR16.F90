    SUBROUTINE BraHRR16ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR,HRRA,HRRB,GRADIENT)
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
    END SUBROUTINE BraHRR16ab
    SUBROUTINE BraHRR16cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,Cart,HRR,GRADIENT)
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
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
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
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(1,7|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*(ABy*HRR(1)+&
                                2.D0*HRR(3))+&
                                HRR(7)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
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
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
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
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=
      !=(1,10|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*(ABz*HRR(1)+&
                                2.D0*HRR(4))+&
                                HRR(10)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
    END SUBROUTINE BraHRR16cd
