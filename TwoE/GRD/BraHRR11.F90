    SUBROUTINE BraHRR11ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR,HRRA,HRRB,GRADIENT)
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
    END SUBROUTINE BraHRR11ab
    SUBROUTINE BraHRR11cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,Cart,HRR,GRADIENT)
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
    END SUBROUTINE BraHRR11cd
