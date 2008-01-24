!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
    SUBROUTINE BraHRR13ab(NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,HRR,HRRA,HRRB,GRADIENT,FP,STRESS)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT NONE
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*),HRRA(*),HRRB(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      REAL(DOUBLE)  :: STRESS(NINT,9),FP(9),DUM
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(1_x,2|
      DUM=ABx*HRRA(2)+&
                         HRRA(5)
      GRADIENT(OffSet,GOA)=DUM+&
                         GRADIENT(OffSet,GOA)
      STRESS(OffSet,1)=DUM*FP(1)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(2)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(3)+&
                         STRESS(OffSet,7)
      !=(1,2_x|
      DUM=-HRR(1)+&
                         ABx*(ABx*HRRB(1)+&
                         2.D0*HRRB(2))+&
                         HRRB(5)
      GRADIENT(OffSet,GOB)=DUM+&
                         GRADIENT(OffSet,GOB)
      STRESS(OffSet,1)=DUM*FP(7)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(8)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(9)+&
                         STRESS(OffSet,7)
      !=(1_y,2|
      DUM=ABx*HRRA(3)+&
                         HRRA(6)
      GRADIENT(OffSet,1 + GOA)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOA)
      STRESS(OffSet,2)=DUM*FP(1)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(2)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(3)+&
                         STRESS(OffSet,8)
      !=(1,2_y|
      DUM=ABy*HRRB(2)+&
                         ABx*(ABy*HRRB(1)+&
                         HRRB(3))+&
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
      !=(1_z,2|
      DUM=ABx*HRRA(4)+&
                         HRRA(8)
      GRADIENT(OffSet,2 + GOA)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOA)
      STRESS(OffSet,3)=DUM*FP(1)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(2)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(3)+&
                         STRESS(OffSet,9)
      !=(1,2_z|
      DUM=ABz*HRRB(2)+&
                         ABx*(ABz*HRRB(1)+&
                         HRRB(4))+&
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
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=
      !=(1_x,3|
      DUM=ABy*HRRA(2)+&
                         HRRA(6)
      GRADIENT(OffSet,GOA)=DUM+&
                         GRADIENT(OffSet,GOA)
      STRESS(OffSet,1)=DUM*FP(1)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(2)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(3)+&
                         STRESS(OffSet,7)
      !=(1,3_x|
      DUM=ABy*HRRB(2)+&
                         ABx*(ABy*HRRB(1)+&
                         HRRB(3))+&
                         HRRB(6)
      GRADIENT(OffSet,GOB)=DUM+&
                         GRADIENT(OffSet,GOB)
      STRESS(OffSet,1)=DUM*FP(7)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(8)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(9)+&
                         STRESS(OffSet,7)
      !=(1_y,3|
      DUM=ABy*HRRA(3)+&
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
      !=(1,3_y|
      DUM=-HRR(1)+&
                         ABy*(ABy*HRRB(1)+&
                         2.D0*HRRB(3))+&
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
      !=(1_z,3|
      DUM=ABy*HRRA(4)+&
                         HRRA(9)
      GRADIENT(OffSet,2 + GOA)=DUM+&
                         GRADIENT(OffSet,2+&
                         GOA)
      STRESS(OffSet,3)=DUM*FP(1)+&
                         STRESS(OffSet,3)
      STRESS(OffSet,6)=DUM*FP(2)+&
                         STRESS(OffSet,6)
      STRESS(OffSet,9)=DUM*FP(3)+&
                         STRESS(OffSet,9)
      !=(1,3_z|
      DUM=ABz*HRRB(3)+&
                         ABy*(ABz*HRRB(1)+&
                         HRRB(4))+&
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
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=
      !=(1_x,4|
      DUM=ABz*HRRA(2)+&
                         HRRA(8)
      GRADIENT(OffSet,GOA)=DUM+&
                         GRADIENT(OffSet,GOA)
      STRESS(OffSet,1)=DUM*FP(1)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(2)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(3)+&
                         STRESS(OffSet,7)
      !=(1,4_x|
      DUM=ABz*HRRB(2)+&
                         ABx*(ABz*HRRB(1)+&
                         HRRB(4))+&
                         HRRB(8)
      GRADIENT(OffSet,GOB)=DUM+&
                         GRADIENT(OffSet,GOB)
      STRESS(OffSet,1)=DUM*FP(7)+&
                         STRESS(OffSet,1)
      STRESS(OffSet,4)=DUM*FP(8)+&
                         STRESS(OffSet,4)
      STRESS(OffSet,7)=DUM*FP(9)+&
                         STRESS(OffSet,7)
      !=(1_y,4|
      DUM=ABz*HRRA(3)+&
                         HRRA(9)
      GRADIENT(OffSet,1 + GOA)=DUM+&
                         GRADIENT(OffSet,1+&
                         GOA)
      STRESS(OffSet,2)=DUM*FP(1)+&
                         STRESS(OffSet,2)
      STRESS(OffSet,5)=DUM*FP(2)+&
                         STRESS(OffSet,5)
      STRESS(OffSet,8)=DUM*FP(3)+&
                         STRESS(OffSet,8)
      !=(1,4_y|
      DUM=ABz*HRRB(3)+&
                         ABy*(ABz*HRRB(1)+&
                         HRRB(4))+&
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
      !=(1_z,4|
      DUM=ABz*HRRA(4)+&
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
      !=(1,4_z|
      DUM=-HRR(1)+&
                         ABz*(ABz*HRRB(1)+&
                         2.D0*HRRB(4))+&
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
    END SUBROUTINE BraHRR13ab
    SUBROUTINE BraHRR13cd(NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,CDOffSet,Cart,HRR,GRADIENT,FP,STRESS)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT NONE
      INTEGER       :: NINT,LDA,LDB,OA,OB,GOA,GOB,GOC,GOD,Cart,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: GRADIENT(NINT,12)
      REAL(DOUBLE)  :: STRESS(NINT,9),FP(9),DUM
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=
      !=(1,2|
      DUM=ABx*HRR(1)+&
                                HRR(2)
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
      !=(1,3|
      DUM=ABy*HRR(1)+&
                                HRR(3)
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
      !=(1,4|
      DUM=ABz*HRR(1)+&
                                HRR(4)
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
    END SUBROUTINE BraHRR13cd
