    SUBROUTINE BraHRR31ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR,HRRA,HRRB,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*),HRRA(*),HRRB(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(2_x,1|
      GRADIENT(OffSet,1 + GOA)=-HRR(1)+&
                         HRRA(5)
      !=(2,1_x|
      GRADIENT(OffSet,1 + GOB)=ABx*HRRB(2)+&
                         HRRB(5)
      !=(2_y,1|
      GRADIENT(OffSet,2 + GOA)=HRRA(6)
      !=(2,1_y|
      GRADIENT(OffSet,2 + GOB)=ABy*HRRB(2)+&
                         HRRB(6)
      !=(2_z,1|
      GRADIENT(OffSet,3 + GOA)=HRRA(8)
      !=(2,1_z|
      GRADIENT(OffSet,3 + GOB)=ABz*HRRB(2)+&
                         HRRB(8)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(3_x,1|
      GRADIENT(OffSet,1 + GOA)=HRRA(6)
      !=(3,1_x|
      GRADIENT(OffSet,1 + GOB)=ABx*HRRB(3)+&
                         HRRB(6)
      !=(3_y,1|
      GRADIENT(OffSet,2 + GOA)=-HRR(1)+&
                         HRRA(7)
      !=(3,1_y|
      GRADIENT(OffSet,2 + GOB)=ABy*HRRB(3)+&
                         HRRB(7)
      !=(3_z,1|
      GRADIENT(OffSet,3 + GOA)=HRRA(9)
      !=(3,1_z|
      GRADIENT(OffSet,3 + GOB)=ABz*HRRB(3)+&
                         HRRB(9)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(4_x,1|
      GRADIENT(OffSet,1 + GOA)=HRRA(8)
      !=(4,1_x|
      GRADIENT(OffSet,1 + GOB)=ABx*HRRB(4)+&
                         HRRB(8)
      !=(4_y,1|
      GRADIENT(OffSet,2 + GOA)=HRRA(9)
      !=(4,1_y|
      GRADIENT(OffSet,2 + GOB)=ABy*HRRB(4)+&
                         HRRB(9)
      !=(4_z,1|
      GRADIENT(OffSet,3 + GOA)=-HRR(1)+&
                         HRRA(10)
      !=(4,1_z|
      GRADIENT(OffSet,3 + GOB)=ABz*HRRB(4)+&
                         HRRB(10)
    END SUBROUTINE BraHRR31ab
    SUBROUTINE BraHRR31cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,Cart,HRR,GRADIENT)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,Cart,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(2,1|
      GRADIENT(OffSet,Cart + GOC)=HRR(2)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(3,1|
      GRADIENT(OffSet,Cart + GOC)=HRR(3)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(4,1|
      GRADIENT(OffSet,Cart + GOC)=HRR(4)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
    END SUBROUTINE BraHRR31cd
