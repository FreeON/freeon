    SUBROUTINE BraHRR32ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR,HRRA,HRRB,GRADIENT,FP,STRESS)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT NONE
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*),HRRA(*),HRRB(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      REAL(DOUBLE)  :: STRESS(NINT,9),FP(9),DUM
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(2_x,1|
      DUM=-HRR(1)+&
                         HRRA(5)
      GRADIENT(OffSet,GOA)=DUM+&
                         GRADIENT(OffSet,GOA)
      STRESS(OffSet,1)=DUM*FP(1)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(2)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(3)+&
                         STRESS(OffSet,7)
      !=(2,1_x|
      DUM=ABx*HRRB(2)+&
                         HRRB(5)
      GRADIENT(OffSet,GOB)=DUM+&
                         GRADIENT(OffSet,GOB)
      STRESS(OffSet,1)=DUM*FP(7)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(8)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(9)+&
                         STRESS(OffSet,7)
      !=(2_y,1|
      DUM=HRRA(6)
      GRADIENT(OffSet,1 + GOA)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOA)
      STRESS(OffSet,2)=DUM*FP(1)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(2)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(3)+&
                         STRESS(OffSet,8)
      !=(2,1_y|
      DUM=ABy*HRRB(2)+&
                         HRRB(6)
      GRADIENT(OffSet,1 + GOB)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOB)
      STRESS(OffSet,2)=DUM*FP(7)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(8)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(9)+&
                         STRESS(OffSet,8)
      !=(2_z,1|
      DUM=HRRA(8)
      GRADIENT(OffSet,2 + GOA)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOA)
      STRESS(OffSet,3)=DUM*FP(1)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(2)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(3)+&
                         STRESS(OffSet,9)
      !=(2,1_z|
      DUM=ABz*HRRB(2)+&
                         HRRB(8)
      GRADIENT(OffSet,2 + GOB)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOB)
      STRESS(OffSet,3)=DUM*FP(7)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(8)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(9)+&
                         STRESS(OffSet,9)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(3_x,1|
      DUM=HRRA(6)
      GRADIENT(OffSet,GOA)=DUM+&
                         GRADIENT(OffSet,GOA)
      STRESS(OffSet,1)=DUM*FP(1)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(2)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(3)+&
                         STRESS(OffSet,7)
      !=(3,1_x|
      DUM=ABx*HRRB(3)+&
                         HRRB(6)
      GRADIENT(OffSet,GOB)=DUM+&
                         GRADIENT(OffSet,GOB)
      STRESS(OffSet,1)=DUM*FP(7)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(8)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(9)+&
                         STRESS(OffSet,7)
      !=(3_y,1|
      DUM=-HRR(1)+&
                         HRRA(7)
      GRADIENT(OffSet,1 + GOA)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOA)
      STRESS(OffSet,2)=DUM*FP(1)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(2)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(3)+&
                         STRESS(OffSet,8)
      !=(3,1_y|
      DUM=ABy*HRRB(3)+&
                         HRRB(7)
      GRADIENT(OffSet,1 + GOB)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOB)
      STRESS(OffSet,2)=DUM*FP(7)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(8)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(9)+&
                         STRESS(OffSet,8)
      !=(3_z,1|
      DUM=HRRA(9)
      GRADIENT(OffSet,2 + GOA)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOA)
      STRESS(OffSet,3)=DUM*FP(1)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(2)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(3)+&
                         STRESS(OffSet,9)
      !=(3,1_z|
      DUM=ABz*HRRB(3)+&
                         HRRB(9)
      GRADIENT(OffSet,2 + GOB)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOB)
      STRESS(OffSet,3)=DUM*FP(7)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(8)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(9)+&
                         STRESS(OffSet,9)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(4_x,1|
      DUM=HRRA(8)
      GRADIENT(OffSet,GOA)=DUM+&
                         GRADIENT(OffSet,GOA)
      STRESS(OffSet,1)=DUM*FP(1)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(2)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(3)+&
                         STRESS(OffSet,7)
      !=(4,1_x|
      DUM=ABx*HRRB(4)+&
                         HRRB(8)
      GRADIENT(OffSet,GOB)=DUM+&
                         GRADIENT(OffSet,GOB)
      STRESS(OffSet,1)=DUM*FP(7)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(8)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(9)+&
                         STRESS(OffSet,7)
      !=(4_y,1|
      DUM=HRRA(9)
      GRADIENT(OffSet,1 + GOA)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOA)
      STRESS(OffSet,2)=DUM*FP(1)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(2)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(3)+&
                         STRESS(OffSet,8)
      !=(4,1_y|
      DUM=ABy*HRRB(4)+&
                         HRRB(9)
      GRADIENT(OffSet,1 + GOB)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOB)
      STRESS(OffSet,2)=DUM*FP(7)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(8)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(9)+&
                         STRESS(OffSet,8)
      !=(4_z,1|
      DUM=-HRR(1)+&
                         HRRA(10)
      GRADIENT(OffSet,2 + GOA)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOA)
      STRESS(OffSet,3)=DUM*FP(1)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(2)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(3)+&
                         STRESS(OffSet,9)
      !=(4,1_z|
      DUM=ABz*HRRB(4)+&
                         HRRB(10)
      GRADIENT(OffSet,2 + GOB)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOB)
      STRESS(OffSet,3)=DUM*FP(7)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(8)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(9)+&
                         STRESS(OffSet,9)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(2_x,2|
      DUM=-HRR(2)+&
                         ABx*(-HRR(1)+&
                         HRRA(5))+&
                         HRRA(11)
      GRADIENT(OffSet,GOA)=DUM+&
                         GRADIENT(OffSet,GOA)
      STRESS(OffSet,1)=DUM*FP(1)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(2)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(3)+&
                         STRESS(OffSet,7)
      !=(2,2_x|
      DUM=-HRR(2)+&
                         ABx*(ABx*HRRB(2)+&
                         2.D0*HRRB(5))+&
                         HRRB(11)
      GRADIENT(OffSet,GOB)=DUM+&
                         GRADIENT(OffSet,GOB)
      STRESS(OffSet,1)=DUM*FP(7)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(8)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(9)+&
                         STRESS(OffSet,7)
      !=(2_y,2|
      DUM=ABx*HRRA(6)+&
                         HRRA(12)
      GRADIENT(OffSet,1 + GOA)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOA)
      STRESS(OffSet,2)=DUM*FP(1)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(2)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(3)+&
                         STRESS(OffSet,8)
      !=(2,2_y|
      DUM=ABy*HRRB(5)+&
                         ABx*(ABy*HRRB(2)+&
                         HRRB(6))+&
                         HRRB(12)
      GRADIENT(OffSet,1 + GOB)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOB)
      STRESS(OffSet,2)=DUM*FP(7)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(8)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(9)+&
                         STRESS(OffSet,8)
      !=(2_z,2|
      DUM=ABx*HRRA(8)+&
                         HRRA(15)
      GRADIENT(OffSet,2 + GOA)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOA)
      STRESS(OffSet,3)=DUM*FP(1)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(2)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(3)+&
                         STRESS(OffSet,9)
      !=(2,2_z|
      DUM=ABz*HRRB(5)+&
                         ABx*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         HRRB(15)
      GRADIENT(OffSet,2 + GOB)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOB)
      STRESS(OffSet,3)=DUM*FP(7)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(8)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(9)+&
                         STRESS(OffSet,9)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(3_x,2|
      DUM=ABx*HRRA(6)+&
                         HRRA(12)
      GRADIENT(OffSet,GOA)=DUM+&
                         GRADIENT(OffSet,GOA)
      STRESS(OffSet,1)=DUM*FP(1)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(2)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(3)+&
                         STRESS(OffSet,7)
      !=(3,2_x|
      DUM=-HRR(3)+&
                         ABx*(ABx*HRRB(3)+&
                         2.D0*HRRB(6))+&
                         HRRB(12)
      GRADIENT(OffSet,GOB)=DUM+&
                         GRADIENT(OffSet,GOB)
      STRESS(OffSet,1)=DUM*FP(7)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(8)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(9)+&
                         STRESS(OffSet,7)
      !=(3_y,2|
      DUM=-HRR(2)+&
                         ABx*(-HRR(1)+&
                         HRRA(7))+&
                         HRRA(13)
      GRADIENT(OffSet,1 + GOA)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOA)
      STRESS(OffSet,2)=DUM*FP(1)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(2)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(3)+&
                         STRESS(OffSet,8)
      !=(3,2_y|
      DUM=ABy*HRRB(6)+&
                         ABx*(ABy*HRRB(3)+&
                         HRRB(7))+&
                         HRRB(13)
      GRADIENT(OffSet,1 + GOB)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOB)
      STRESS(OffSet,2)=DUM*FP(7)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(8)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(9)+&
                         STRESS(OffSet,8)
      !=(3_z,2|
      DUM=ABx*HRRA(9)+&
                         HRRA(16)
      GRADIENT(OffSet,2 + GOA)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOA)
      STRESS(OffSet,3)=DUM*FP(1)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(2)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(3)+&
                         STRESS(OffSet,9)
      !=(3,2_z|
      DUM=ABz*HRRB(6)+&
                         ABx*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         HRRB(16)
      GRADIENT(OffSet,2 + GOB)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOB)
      STRESS(OffSet,3)=DUM*FP(7)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(8)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(9)+&
                         STRESS(OffSet,9)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(4_x,2|
      DUM=ABx*HRRA(8)+&
                         HRRA(15)
      GRADIENT(OffSet,GOA)=DUM+&
                         GRADIENT(OffSet,GOA)
      STRESS(OffSet,1)=DUM*FP(1)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(2)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(3)+&
                         STRESS(OffSet,7)
      !=(4,2_x|
      DUM=-HRR(4)+&
                         ABx*(ABx*HRRB(4)+&
                         2.D0*HRRB(8))+&
                         HRRB(15)
      GRADIENT(OffSet,GOB)=DUM+&
                         GRADIENT(OffSet,GOB)
      STRESS(OffSet,1)=DUM*FP(7)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(8)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(9)+&
                         STRESS(OffSet,7)
      !=(4_y,2|
      DUM=ABx*HRRA(9)+&
                         HRRA(16)
      GRADIENT(OffSet,1 + GOA)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOA)
      STRESS(OffSet,2)=DUM*FP(1)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(2)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(3)+&
                         STRESS(OffSet,8)
      !=(4,2_y|
      DUM=ABy*HRRB(8)+&
                         ABx*(ABy*HRRB(4)+&
                         HRRB(9))+&
                         HRRB(16)
      GRADIENT(OffSet,1 + GOB)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOB)
      STRESS(OffSet,2)=DUM*FP(7)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(8)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(9)+&
                         STRESS(OffSet,8)
      !=(4_z,2|
      DUM=-HRR(2)+&
                         ABx*(-HRR(1)+&
                         HRRA(10))+&
                         HRRA(18)
      GRADIENT(OffSet,2 + GOA)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOA)
      STRESS(OffSet,3)=DUM*FP(1)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(2)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(3)+&
                         STRESS(OffSet,9)
      !=(4,2_z|
      DUM=ABz*HRRB(8)+&
                         ABx*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         HRRB(18)
      GRADIENT(OffSet,2 + GOB)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOB)
      STRESS(OffSet,3)=DUM*FP(7)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(8)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(9)+&
                         STRESS(OffSet,9)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(2_x,3|
      DUM=-HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(5))+&
                         HRRA(12)
      GRADIENT(OffSet,GOA)=DUM+&
                         GRADIENT(OffSet,GOA)
      STRESS(OffSet,1)=DUM*FP(1)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(2)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(3)+&
                         STRESS(OffSet,7)
      !=(2,3_x|
      DUM=ABy*HRRB(5)+&
                         ABx*(ABy*HRRB(2)+&
                         HRRB(6))+&
                         HRRB(12)
      GRADIENT(OffSet,GOB)=DUM+&
                         GRADIENT(OffSet,GOB)
      STRESS(OffSet,1)=DUM*FP(7)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(8)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(9)+&
                         STRESS(OffSet,7)
      !=(2_y,3|
      DUM=ABy*HRRA(6)+&
                         HRRA(13)
      GRADIENT(OffSet,1 + GOA)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOA)
      STRESS(OffSet,2)=DUM*FP(1)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(2)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(3)+&
                         STRESS(OffSet,8)
      !=(2,3_y|
      DUM=-HRR(2)+&
                         ABy*(ABy*HRRB(2)+&
                         2.D0*HRRB(6))+&
                         HRRB(13)
      GRADIENT(OffSet,1 + GOB)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOB)
      STRESS(OffSet,2)=DUM*FP(7)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(8)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(9)+&
                         STRESS(OffSet,8)
      !=(2_z,3|
      DUM=ABy*HRRA(8)+&
                         HRRA(16)
      GRADIENT(OffSet,2 + GOA)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOA)
      STRESS(OffSet,3)=DUM*FP(1)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(2)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(3)+&
                         STRESS(OffSet,9)
      !=(2,3_z|
      DUM=ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         HRRB(16)
      GRADIENT(OffSet,2 + GOB)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOB)
      STRESS(OffSet,3)=DUM*FP(7)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(8)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(9)+&
                         STRESS(OffSet,9)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(3_x,3|
      DUM=ABy*HRRA(6)+&
                         HRRA(13)
      GRADIENT(OffSet,GOA)=DUM+&
                         GRADIENT(OffSet,GOA)
      STRESS(OffSet,1)=DUM*FP(1)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(2)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(3)+&
                         STRESS(OffSet,7)
      !=(3,3_x|
      DUM=ABy*HRRB(6)+&
                         ABx*(ABy*HRRB(3)+&
                         HRRB(7))+&
                         HRRB(13)
      GRADIENT(OffSet,GOB)=DUM+&
                         GRADIENT(OffSet,GOB)
      STRESS(OffSet,1)=DUM*FP(7)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(8)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(9)+&
                         STRESS(OffSet,7)
      !=(3_y,3|
      DUM=-HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(7))+&
                         HRRA(14)
      GRADIENT(OffSet,1 + GOA)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOA)
      STRESS(OffSet,2)=DUM*FP(1)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(2)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(3)+&
                         STRESS(OffSet,8)
      !=(3,3_y|
      DUM=-HRR(3)+&
                         ABy*(ABy*HRRB(3)+&
                         2.D0*HRRB(7))+&
                         HRRB(14)
      GRADIENT(OffSet,1 + GOB)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOB)
      STRESS(OffSet,2)=DUM*FP(7)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(8)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(9)+&
                         STRESS(OffSet,8)
      !=(3_z,3|
      DUM=ABy*HRRA(9)+&
                         HRRA(17)
      GRADIENT(OffSet,2 + GOA)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOA)
      STRESS(OffSet,3)=DUM*FP(1)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(2)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(3)+&
                         STRESS(OffSet,9)
      !=(3,3_z|
      DUM=ABz*HRRB(7)+&
                         ABy*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         HRRB(17)
      GRADIENT(OffSet,2 + GOB)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOB)
      STRESS(OffSet,3)=DUM*FP(7)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(8)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(9)+&
                         STRESS(OffSet,9)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(4_x,3|
      DUM=ABy*HRRA(8)+&
                         HRRA(16)
      GRADIENT(OffSet,GOA)=DUM+&
                         GRADIENT(OffSet,GOA)
      STRESS(OffSet,1)=DUM*FP(1)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(2)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(3)+&
                         STRESS(OffSet,7)
      !=(4,3_x|
      DUM=ABy*HRRB(8)+&
                         ABx*(ABy*HRRB(4)+&
                         HRRB(9))+&
                         HRRB(16)
      GRADIENT(OffSet,GOB)=DUM+&
                         GRADIENT(OffSet,GOB)
      STRESS(OffSet,1)=DUM*FP(7)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(8)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(9)+&
                         STRESS(OffSet,7)
      !=(4_y,3|
      DUM=ABy*HRRA(9)+&
                         HRRA(17)
      GRADIENT(OffSet,1 + GOA)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOA)
      STRESS(OffSet,2)=DUM*FP(1)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(2)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(3)+&
                         STRESS(OffSet,8)
      !=(4,3_y|
      DUM=-HRR(4)+&
                         ABy*(ABy*HRRB(4)+&
                         2.D0*HRRB(9))+&
                         HRRB(17)
      GRADIENT(OffSet,1 + GOB)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOB)
      STRESS(OffSet,2)=DUM*FP(7)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(8)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(9)+&
                         STRESS(OffSet,8)
      !=(4_z,3|
      DUM=-HRR(3)+&
                         ABy*(-HRR(1)+&
                         HRRA(10))+&
                         HRRA(19)
      GRADIENT(OffSet,2 + GOA)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOA)
      STRESS(OffSet,3)=DUM*FP(1)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(2)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(3)+&
                         STRESS(OffSet,9)
      !=(4,3_z|
      DUM=ABz*HRRB(9)+&
                         ABy*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         HRRB(19)
      GRADIENT(OffSet,2 + GOB)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOB)
      STRESS(OffSet,3)=DUM*FP(7)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(8)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(9)+&
                         STRESS(OffSet,9)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(2_x,4|
      DUM=-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(5))+&
                         HRRA(15)
      GRADIENT(OffSet,GOA)=DUM+&
                         GRADIENT(OffSet,GOA)
      STRESS(OffSet,1)=DUM*FP(1)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(2)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(3)+&
                         STRESS(OffSet,7)
      !=(2,4_x|
      DUM=ABz*HRRB(5)+&
                         ABx*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         HRRB(15)
      GRADIENT(OffSet,GOB)=DUM+&
                         GRADIENT(OffSet,GOB)
      STRESS(OffSet,1)=DUM*FP(7)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(8)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(9)+&
                         STRESS(OffSet,7)
      !=(2_y,4|
      DUM=ABz*HRRA(6)+&
                         HRRA(16)
      GRADIENT(OffSet,1 + GOA)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOA)
      STRESS(OffSet,2)=DUM*FP(1)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(2)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(3)+&
                         STRESS(OffSet,8)
      !=(2,4_y|
      DUM=ABz*HRRB(6)+&
                         ABy*(ABz*HRRB(2)+&
                         HRRB(8))+&
                         HRRB(16)
      GRADIENT(OffSet,1 + GOB)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOB)
      STRESS(OffSet,2)=DUM*FP(7)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(8)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(9)+&
                         STRESS(OffSet,8)
      !=(2_z,4|
      DUM=ABz*HRRA(8)+&
                         HRRA(18)
      GRADIENT(OffSet,2 + GOA)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOA)
      STRESS(OffSet,3)=DUM*FP(1)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(2)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(3)+&
                         STRESS(OffSet,9)
      !=(2,4_z|
      DUM=-HRR(2)+&
                         ABz*(ABz*HRRB(2)+&
                         2.D0*HRRB(8))+&
                         HRRB(18)
      GRADIENT(OffSet,2 + GOB)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOB)
      STRESS(OffSet,3)=DUM*FP(7)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(8)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(9)+&
                         STRESS(OffSet,9)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(3_x,4|
      DUM=ABz*HRRA(6)+&
                         HRRA(16)
      GRADIENT(OffSet,GOA)=DUM+&
                         GRADIENT(OffSet,GOA)
      STRESS(OffSet,1)=DUM*FP(1)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(2)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(3)+&
                         STRESS(OffSet,7)
      !=(3,4_x|
      DUM=ABz*HRRB(6)+&
                         ABx*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         HRRB(16)
      GRADIENT(OffSet,GOB)=DUM+&
                         GRADIENT(OffSet,GOB)
      STRESS(OffSet,1)=DUM*FP(7)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(8)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(9)+&
                         STRESS(OffSet,7)
      !=(3_y,4|
      DUM=-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(7))+&
                         HRRA(17)
      GRADIENT(OffSet,1 + GOA)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOA)
      STRESS(OffSet,2)=DUM*FP(1)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(2)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(3)+&
                         STRESS(OffSet,8)
      !=(3,4_y|
      DUM=ABz*HRRB(7)+&
                         ABy*(ABz*HRRB(3)+&
                         HRRB(9))+&
                         HRRB(17)
      GRADIENT(OffSet,1 + GOB)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOB)
      STRESS(OffSet,2)=DUM*FP(7)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(8)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(9)+&
                         STRESS(OffSet,8)
      !=(3_z,4|
      DUM=ABz*HRRA(9)+&
                         HRRA(19)
      GRADIENT(OffSet,2 + GOA)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOA)
      STRESS(OffSet,3)=DUM*FP(1)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(2)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(3)+&
                         STRESS(OffSet,9)
      !=(3,4_z|
      DUM=-HRR(3)+&
                         ABz*(ABz*HRRB(3)+&
                         2.D0*HRRB(9))+&
                         HRRB(19)
      GRADIENT(OffSet,2 + GOB)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOB)
      STRESS(OffSet,3)=DUM*FP(7)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(8)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(9)+&
                         STRESS(OffSet,9)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(4_x,4|
      DUM=ABz*HRRA(8)+&
                         HRRA(18)
      GRADIENT(OffSet,GOA)=DUM+&
                         GRADIENT(OffSet,GOA)
      STRESS(OffSet,1)=DUM*FP(1)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(2)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(3)+&
                         STRESS(OffSet,7)
      !=(4,4_x|
      DUM=ABz*HRRB(8)+&
                         ABx*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         HRRB(18)
      GRADIENT(OffSet,GOB)=DUM+&
                         GRADIENT(OffSet,GOB)
      STRESS(OffSet,1)=DUM*FP(7)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(8)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(9)+&
                         STRESS(OffSet,7)
      !=(4_y,4|
      DUM=ABz*HRRA(9)+&
                         HRRA(19)
      GRADIENT(OffSet,1 + GOA)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOA)
      STRESS(OffSet,2)=DUM*FP(1)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(2)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(3)+&
                         STRESS(OffSet,8)
      !=(4,4_y|
      DUM=ABz*HRRB(9)+&
                         ABy*(ABz*HRRB(4)+&
                         HRRB(10))+&
                         HRRB(19)
      GRADIENT(OffSet,1 + GOB)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOB)
      STRESS(OffSet,2)=DUM*FP(7)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(8)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(9)+&
                         STRESS(OffSet,8)
      !=(4_z,4|
      DUM=-HRR(4)+&
                         ABz*(-HRR(1)+&
                         HRRA(10))+&
                         HRRA(20)
      GRADIENT(OffSet,2 + GOA)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOA)
      STRESS(OffSet,3)=DUM*FP(1)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(2)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(3)+&
                         STRESS(OffSet,9)
      !=(4,4_z|
      DUM=-HRR(4)+&
                         ABz*(ABz*HRRB(4)+&
                         2.D0*HRRB(10))+&
                         HRRB(20)
      GRADIENT(OffSet,2 + GOB)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOB)
      STRESS(OffSet,3)=DUM*FP(7)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(8)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(9)+&
                         STRESS(OffSet,9)
    END SUBROUTINE BraHRR32ab
    SUBROUTINE BraHRR32cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,Cart,HRR,GRADIENT,FP,STRESS)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT NONE
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,Cart,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      REAL(DOUBLE)  :: STRESS(NINT,9),FP(9),DUM
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(2,1|
      DUM=HRR(2)
      GRADIENT(OffSet,Cart + GOC)=DUM+&
                                GRADIENT(OffSet,Cart+&
                                GOC)
      STRESS(OffSet,1+Cart)=DUM*FP(4)+&
                                STRESS(OffSet,1+&
                                Cart)
      STRESS(OffSet,4+Cart)=DUM*FP(5)+&
                                STRESS(OffSet,4+&
                                Cart)
      STRESS(OffSet,7+Cart)=DUM*FP(6)+&
                                STRESS(OffSet,7+&
                                Cart)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(3,1|
      DUM=HRR(3)
      GRADIENT(OffSet,Cart + GOC)=DUM+&
                                GRADIENT(OffSet,Cart+&
                                GOC)
      STRESS(OffSet,1+Cart)=DUM*FP(4)+&
                                STRESS(OffSet,1+&
                                Cart)
      STRESS(OffSet,4+Cart)=DUM*FP(5)+&
                                STRESS(OffSet,4+&
                                Cart)
      STRESS(OffSet,7+Cart)=DUM*FP(6)+&
                                STRESS(OffSet,7+&
                                Cart)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(4,1|
      DUM=HRR(4)
      GRADIENT(OffSet,Cart + GOC)=DUM+&
                                GRADIENT(OffSet,Cart+&
                                GOC)
      STRESS(OffSet,1+Cart)=DUM*FP(4)+&
                                STRESS(OffSet,1+&
                                Cart)
      STRESS(OffSet,4+Cart)=DUM*FP(5)+&
                                STRESS(OffSet,4+&
                                Cart)
      STRESS(OffSet,7+Cart)=DUM*FP(6)+&
                                STRESS(OffSet,7+&
                                Cart)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(2,2|
      DUM=ABx*HRR(2)+&
                                HRR(5)
      GRADIENT(OffSet,Cart + GOC)=DUM+&
                                GRADIENT(OffSet,Cart+&
                                GOC)
      STRESS(OffSet,1+Cart)=DUM*FP(4)+&
                                STRESS(OffSet,1+&
                                Cart)
      STRESS(OffSet,4+Cart)=DUM*FP(5)+&
                                STRESS(OffSet,4+&
                                Cart)
      STRESS(OffSet,7+Cart)=DUM*FP(6)+&
                                STRESS(OffSet,7+&
                                Cart)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(3,2|
      DUM=ABx*HRR(3)+&
                                HRR(6)
      GRADIENT(OffSet,Cart + GOC)=DUM+&
                                GRADIENT(OffSet,Cart+&
                                GOC)
      STRESS(OffSet,1+Cart)=DUM*FP(4)+&
                                STRESS(OffSet,1+&
                                Cart)
      STRESS(OffSet,4+Cart)=DUM*FP(5)+&
                                STRESS(OffSet,4+&
                                Cart)
      STRESS(OffSet,7+Cart)=DUM*FP(6)+&
                                STRESS(OffSet,7+&
                                Cart)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(4,2|
      DUM=ABx*HRR(4)+&
                                HRR(8)
      GRADIENT(OffSet,Cart + GOC)=DUM+&
                                GRADIENT(OffSet,Cart+&
                                GOC)
      STRESS(OffSet,1+Cart)=DUM*FP(4)+&
                                STRESS(OffSet,1+&
                                Cart)
      STRESS(OffSet,4+Cart)=DUM*FP(5)+&
                                STRESS(OffSet,4+&
                                Cart)
      STRESS(OffSet,7+Cart)=DUM*FP(6)+&
                                STRESS(OffSet,7+&
                                Cart)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(2,3|
      DUM=ABy*HRR(2)+&
                                HRR(6)
      GRADIENT(OffSet,Cart + GOC)=DUM+&
                                GRADIENT(OffSet,Cart+&
                                GOC)
      STRESS(OffSet,1+Cart)=DUM*FP(4)+&
                                STRESS(OffSet,1+&
                                Cart)
      STRESS(OffSet,4+Cart)=DUM*FP(5)+&
                                STRESS(OffSet,4+&
                                Cart)
      STRESS(OffSet,7+Cart)=DUM*FP(6)+&
                                STRESS(OffSet,7+&
                                Cart)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(3,3|
      DUM=ABy*HRR(3)+&
                                HRR(7)
      GRADIENT(OffSet,Cart + GOC)=DUM+&
                                GRADIENT(OffSet,Cart+&
                                GOC)
      STRESS(OffSet,1+Cart)=DUM*FP(4)+&
                                STRESS(OffSet,1+&
                                Cart)
      STRESS(OffSet,4+Cart)=DUM*FP(5)+&
                                STRESS(OffSet,4+&
                                Cart)
      STRESS(OffSet,7+Cart)=DUM*FP(6)+&
                                STRESS(OffSet,7+&
                                Cart)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(4,3|
      DUM=ABy*HRR(4)+&
                                HRR(9)
      GRADIENT(OffSet,Cart + GOC)=DUM+&
                                GRADIENT(OffSet,Cart+&
                                GOC)
      STRESS(OffSet,1+Cart)=DUM*FP(4)+&
                                STRESS(OffSet,1+&
                                Cart)
      STRESS(OffSet,4+Cart)=DUM*FP(5)+&
                                STRESS(OffSet,4+&
                                Cart)
      STRESS(OffSet,7+Cart)=DUM*FP(6)+&
                                STRESS(OffSet,7+&
                                Cart)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(2,4|
      DUM=ABz*HRR(2)+&
                                HRR(8)
      GRADIENT(OffSet,Cart + GOC)=DUM+&
                                GRADIENT(OffSet,Cart+&
                                GOC)
      STRESS(OffSet,1+Cart)=DUM*FP(4)+&
                                STRESS(OffSet,1+&
                                Cart)
      STRESS(OffSet,4+Cart)=DUM*FP(5)+&
                                STRESS(OffSet,4+&
                                Cart)
      STRESS(OffSet,7+Cart)=DUM*FP(6)+&
                                STRESS(OffSet,7+&
                                Cart)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(3,4|
      DUM=ABz*HRR(3)+&
                                HRR(9)
      GRADIENT(OffSet,Cart + GOC)=DUM+&
                                GRADIENT(OffSet,Cart+&
                                GOC)
      STRESS(OffSet,1+Cart)=DUM*FP(4)+&
                                STRESS(OffSet,1+&
                                Cart)
      STRESS(OffSet,4+Cart)=DUM*FP(5)+&
                                STRESS(OffSet,4+&
                                Cart)
      STRESS(OffSet,7+Cart)=DUM*FP(6)+&
                                STRESS(OffSet,7+&
                                Cart)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=
      !=(4,4|
      DUM=ABz*HRR(4)+&
                                HRR(10)
      GRADIENT(OffSet,Cart + GOC)=DUM+&
                                GRADIENT(OffSet,Cart+&
                                GOC)
      STRESS(OffSet,1+Cart)=DUM*FP(4)+&
                                STRESS(OffSet,1+&
                                Cart)
      STRESS(OffSet,4+Cart)=DUM*FP(5)+&
                                STRESS(OffSet,4+&
                                Cart)
      STRESS(OffSet,7+Cart)=DUM*FP(6)+&
                                STRESS(OffSet,7+&
                                Cart)
      GRADIENT(OffSet,Cart + GOD)=-GRADIENT(OffSet,Cart+GOA)&
                                -GRADIENT(OffSet,Cart+GOB)&
                                -GRADIENT(OffSet,Cart+GOC)
    END SUBROUTINE BraHRR32cd
