    SUBROUTINE BraHRR101ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR,HRRA,HRRB,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*),HRRA(*),HRRB(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(11_x,1|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-3.D0*HRR(5)+&
                         HRRA(21)
      !=(11,1_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABx*HRRB(11)+&
                         HRRB(21)
      !=(11_y,1|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         HRRA(22)
      !=(11,1_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(11)+&
                         HRRB(22)
      !=(11_z,1|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         HRRA(26)
      !=(11,1_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(11)+&
                         HRRB(26)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(12_x,1|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(6)+&
                         HRRA(22)
      !=(12,1_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABx*HRRB(12)+&
                         HRRB(22)
      !=(12_y,1|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(5)+&
                         HRRA(23)
      !=(12,1_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(12)+&
                         HRRB(23)
      !=(12_z,1|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         HRRA(27)
      !=(12,1_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(12)+&
                         HRRB(27)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(13_x,1|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(7)+&
                         HRRA(23)
      !=(13,1_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABx*HRRB(13)+&
                         HRRB(23)
      !=(13_y,1|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(6)+&
                         HRRA(24)
      !=(13,1_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(13)+&
                         HRRB(24)
      !=(13_z,1|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         HRRA(28)
      !=(13,1_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(13)+&
                         HRRB(28)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(14_x,1|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         HRRA(24)
      !=(14,1_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABx*HRRB(14)+&
                         HRRB(24)
      !=(14_y,1|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-3.D0*HRR(7)+&
                         HRRA(25)
      !=(14,1_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(14)+&
                         HRRB(25)
      !=(14_z,1|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)+&
                         HRRA(29)
      !=(14,1_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(14)+&
                         HRRB(29)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(15_x,1|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-2.D0*HRR(8)+&
                         HRRA(26)
      !=(15,1_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABx*HRRB(15)+&
                         HRRB(26)
      !=(15_y,1|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         HRRA(27)
      !=(15,1_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(15)+&
                         HRRB(27)
      !=(15_z,1|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(5)+&
                         HRRA(30)
      !=(15,1_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(15)+&
                         HRRB(30)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(16_x,1|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(9)+&
                         HRRA(27)
      !=(16,1_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABx*HRRB(16)+&
                         HRRB(27)
      !=(16_y,1|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(8)+&
                         HRRA(28)
      !=(16,1_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(16)+&
                         HRRB(28)
      !=(16_z,1|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(6)+&
                         HRRA(31)
      !=(16,1_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(16)+&
                         HRRB(31)
      OffSet=(OA+6)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(17_x,1|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         HRRA(28)
      !=(17,1_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABx*HRRB(17)+&
                         HRRB(28)
      !=(17_y,1|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-2.D0*HRR(9)+&
                         HRRA(29)
      !=(17,1_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(17)+&
                         HRRB(29)
      !=(17_z,1|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-HRR(7)+&
                         HRRA(32)
      !=(17,1_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(17)+&
                         HRRB(32)
      OffSet=(OA+7)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(18_x,1|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)-HRR(10)+&
                         HRRA(30)
      !=(18,1_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABx*HRRB(18)+&
                         HRRB(30)
      !=(18_y,1|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         HRRA(31)
      !=(18,1_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(18)+&
                         HRRB(31)
      !=(18_z,1|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(8)+&
                         HRRA(33)
      !=(18,1_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(18)+&
                         HRRB(33)
      OffSet=(OA+8)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(19_x,1|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         HRRA(31)
      !=(19,1_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABx*HRRB(19)+&
                         HRRB(31)
      !=(19_y,1|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)-HRR(10)+&
                         HRRA(32)
      !=(19,1_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(19)+&
                         HRRB(32)
      !=(19_z,1|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-2.D0*HRR(9)+&
                         HRRA(34)
      !=(19,1_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(19)+&
                         HRRB(34)
      OffSet=(OA+9)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(20_x,1|
      GRADIENT(OffSet,GOA)=GRADIENT(OffSet,GOA)+&
                         HRRA(33)
      !=(20,1_x|
      GRADIENT(OffSet,GOB)=GRADIENT(OffSet,GOB)+&
                         ABx*HRRB(20)+&
                         HRRB(33)
      !=(20_y,1|
      GRADIENT(OffSet,1 + GOA)=GRADIENT(OffSet,1+&
                         GOA)+&
                         HRRA(34)
      !=(20,1_y|
      GRADIENT(OffSet,1 + GOB)=GRADIENT(OffSet,1+&
                         GOB)+&
                         ABy*HRRB(20)+&
                         HRRB(34)
      !=(20_z,1|
      GRADIENT(OffSet,2 + GOA)=GRADIENT(OffSet,2+&
                         GOA)-3.D0*HRR(10)+&
                         HRRA(35)
      !=(20,1_z|
      GRADIENT(OffSet,2 + GOB)=GRADIENT(OffSet,2+&
                         GOB)+&
                         ABz*HRRB(20)+&
                         HRRB(35)
    END SUBROUTINE BraHRR101ab
    SUBROUTINE BraHRR101cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,Cart,HRR,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,Cart,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(11,1|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                HRR(11)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(12,1|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                HRR(12)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(13,1|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                HRR(13)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(14,1|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                HRR(14)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(15,1|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                HRR(15)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(16,1|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                HRR(16)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+6)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(17,1|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                HRR(17)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+7)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(18,1|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                HRR(18)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+8)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(19,1|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                HRR(19)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+9)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(20,1|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                HRR(20)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
    END SUBROUTINE BraHRR101cd
