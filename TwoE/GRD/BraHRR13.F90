    SUBROUTINE BraHRR13ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR,HRRA,HRRB,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*),HRRA(*),HRRB(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
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
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
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
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
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
    END SUBROUTINE BraHRR13ab
    SUBROUTINE BraHRR13cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,Cart,HRR,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,Cart,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(1,2|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABx*HRR(1)+&
                                HRR(2)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(1,3|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABy*HRR(1)+&
                                HRR(3)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(1,4|
      GRADIENT(OffSet,Cart + GOC)=GRADIENT(OffSet,cart+&
                                GOC)+&
                                ABz*HRR(1)+&
                                HRR(4)
      GRADIENT(OffSet,Cart + GOD)= GRADIENT(OffSet,Cart+GOD)&
                                -GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
    END SUBROUTINE BraHRR13cd
