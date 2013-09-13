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
   SUBROUTINE BraHRR26(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(1,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(21)+  &
        2.D0*HRR(22))+  &
        HRR(25)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(1,6|
      INTGRL(OffSet)=ABy*HRR(22)+  &
        ABx*(ABy*HRR(21)+  &
        HRR(23))+  &
        HRR(26)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(1,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(21)+  &
        2.D0*HRR(23))+  &
        HRR(27)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=(1,8|
      INTGRL(OffSet)=ABz*HRR(22)+  &
        ABx*(ABz*HRR(21)+  &
        HRR(24))+  &
        HRR(28)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=(1,9|
      INTGRL(OffSet)=ABz*HRR(23)+  &
        ABy*(ABz*HRR(21)+  &
        HRR(24))+  &
        HRR(29)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=(1,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(21)+  &
        2.D0*HRR(24))+  &
        HRR(30)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(2,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(2)+  &
        2.D0*HRR(5))+  &
        HRR(11)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(3,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(3)+  &
        2.D0*HRR(6))+  &
        HRR(12)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=(4,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(4)+  &
        2.D0*HRR(8))+  &
        HRR(15)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=(2,6|
      INTGRL(OffSet)=ABy*HRR(5)+  &
        ABx*(ABy*HRR(2)+  &
        HRR(6))+  &
        HRR(12)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=(3,6|
      INTGRL(OffSet)=ABy*HRR(6)+  &
        ABx*(ABy*HRR(3)+  &
        HRR(7))+  &
        HRR(13)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=(4,6|
      INTGRL(OffSet)=ABy*HRR(8)+  &
        ABx*(ABy*HRR(4)+  &
        HRR(9))+  &
        HRR(16)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=(2,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(2)+  &
        2.D0*HRR(6))+  &
        HRR(13)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=(3,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(3)+  &
        2.D0*HRR(7))+  &
        HRR(14)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=(4,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(4)+  &
        2.D0*HRR(9))+  &
        HRR(17)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=(2,8|
      INTGRL(OffSet)=ABz*HRR(5)+  &
        ABx*(ABz*HRR(2)+  &
        HRR(8))+  &
        HRR(15)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=(3,8|
      INTGRL(OffSet)=ABz*HRR(6)+  &
        ABx*(ABz*HRR(3)+  &
        HRR(9))+  &
        HRR(16)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+3)*LDB+CDOffSet !=(4,8|
      INTGRL(OffSet)=ABz*HRR(8)+  &
        ABx*(ABz*HRR(4)+  &
        HRR(10))+  &
        HRR(18)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+4)*LDB+CDOffSet !=(2,9|
      INTGRL(OffSet)=ABz*HRR(6)+  &
        ABy*(ABz*HRR(2)+  &
        HRR(8))+  &
        HRR(16)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+4)*LDB+CDOffSet !=(3,9|
      INTGRL(OffSet)=ABz*HRR(7)+  &
        ABy*(ABz*HRR(3)+  &
        HRR(9))+  &
        HRR(17)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+4)*LDB+CDOffSet !=(4,9|
      INTGRL(OffSet)=ABz*HRR(9)+  &
        ABy*(ABz*HRR(4)+  &
        HRR(10))+  &
        HRR(19)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+5)*LDB+CDOffSet !=(2,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(2)+  &
        2.D0*HRR(8))+  &
        HRR(18)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+5)*LDB+CDOffSet !=(3,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(3)+  &
        2.D0*HRR(9))+  &
        HRR(19)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+5)*LDB+CDOffSet !=(4,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(4)+  &
        2.D0*HRR(10))+  &
        HRR(20)+  &
        INTGRL(OffSet)
END SUBROUTINE BraHRR26
