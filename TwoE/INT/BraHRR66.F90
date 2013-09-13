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
   SUBROUTINE BraHRR66(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(5,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(5)+  &
        2.D0*HRR(11))+  &
        HRR(21)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(6,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(6)+  &
        2.D0*HRR(12))+  &
        HRR(22)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(7,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(7)+  &
        2.D0*HRR(13))+  &
        HRR(23)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=(8,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(8)+  &
        2.D0*HRR(15))+  &
        HRR(26)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=(9,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(9)+  &
        2.D0*HRR(16))+  &
        HRR(27)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=(10,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(10)+  &
        2.D0*HRR(18))+  &
        HRR(30)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(5,6|
      INTGRL(OffSet)=ABy*HRR(11)+  &
        ABx*(ABy*HRR(5)+  &
        HRR(12))+  &
        HRR(22)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=(6,6|
      INTGRL(OffSet)=ABy*HRR(12)+  &
        ABx*(ABy*HRR(6)+  &
        HRR(13))+  &
        HRR(23)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=(7,6|
      INTGRL(OffSet)=ABy*HRR(13)+  &
        ABx*(ABy*HRR(7)+  &
        HRR(14))+  &
        HRR(24)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=(8,6|
      INTGRL(OffSet)=ABy*HRR(15)+  &
        ABx*(ABy*HRR(8)+  &
        HRR(16))+  &
        HRR(27)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+1)*LDB+CDOffSet !=(9,6|
      INTGRL(OffSet)=ABy*HRR(16)+  &
        ABx*(ABy*HRR(9)+  &
        HRR(17))+  &
        HRR(28)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+1)*LDB+CDOffSet !=(10,6|
      INTGRL(OffSet)=ABy*HRR(18)+  &
        ABx*(ABy*HRR(10)+  &
        HRR(19))+  &
        HRR(31)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(5,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(5)+  &
        2.D0*HRR(12))+  &
        HRR(23)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=(6,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(6)+  &
        2.D0*HRR(13))+  &
        HRR(24)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=(7,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(7)+  &
        2.D0*HRR(14))+  &
        HRR(25)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=(8,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(8)+  &
        2.D0*HRR(16))+  &
        HRR(28)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+2)*LDB+CDOffSet !=(9,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(9)+  &
        2.D0*HRR(17))+  &
        HRR(29)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+2)*LDB+CDOffSet !=(10,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(10)+  &
        2.D0*HRR(19))+  &
        HRR(32)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=(5,8|
      INTGRL(OffSet)=ABz*HRR(11)+  &
        ABx*(ABz*HRR(5)+  &
        HRR(15))+  &
        HRR(26)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=(6,8|
      INTGRL(OffSet)=ABz*HRR(12)+  &
        ABx*(ABz*HRR(6)+  &
        HRR(16))+  &
        HRR(27)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=(7,8|
      INTGRL(OffSet)=ABz*HRR(13)+  &
        ABx*(ABz*HRR(7)+  &
        HRR(17))+  &
        HRR(28)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+3)*LDB+CDOffSet !=(8,8|
      INTGRL(OffSet)=ABz*HRR(15)+  &
        ABx*(ABz*HRR(8)+  &
        HRR(18))+  &
        HRR(30)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+3)*LDB+CDOffSet !=(9,8|
      INTGRL(OffSet)=ABz*HRR(16)+  &
        ABx*(ABz*HRR(9)+  &
        HRR(19))+  &
        HRR(31)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+3)*LDB+CDOffSet !=(10,8|
      INTGRL(OffSet)=ABz*HRR(18)+  &
        ABx*(ABz*HRR(10)+  &
        HRR(20))+  &
        HRR(33)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=(5,9|
      INTGRL(OffSet)=ABz*HRR(12)+  &
        ABy*(ABz*HRR(5)+  &
        HRR(15))+  &
        HRR(27)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+4)*LDB+CDOffSet !=(6,9|
      INTGRL(OffSet)=ABz*HRR(13)+  &
        ABy*(ABz*HRR(6)+  &
        HRR(16))+  &
        HRR(28)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+4)*LDB+CDOffSet !=(7,9|
      INTGRL(OffSet)=ABz*HRR(14)+  &
        ABy*(ABz*HRR(7)+  &
        HRR(17))+  &
        HRR(29)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+4)*LDB+CDOffSet !=(8,9|
      INTGRL(OffSet)=ABz*HRR(16)+  &
        ABy*(ABz*HRR(8)+  &
        HRR(18))+  &
        HRR(31)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+4)*LDB+CDOffSet !=(9,9|
      INTGRL(OffSet)=ABz*HRR(17)+  &
        ABy*(ABz*HRR(9)+  &
        HRR(19))+  &
        HRR(32)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+4)*LDB+CDOffSet !=(10,9|
      INTGRL(OffSet)=ABz*HRR(19)+  &
        ABy*(ABz*HRR(10)+  &
        HRR(20))+  &
        HRR(34)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=(5,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(5)+  &
        2.D0*HRR(15))+  &
        HRR(30)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+5)*LDB+CDOffSet !=(6,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(6)+  &
        2.D0*HRR(16))+  &
        HRR(31)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+5)*LDB+CDOffSet !=(7,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(7)+  &
        2.D0*HRR(17))+  &
        HRR(32)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+5)*LDB+CDOffSet !=(8,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(8)+  &
        2.D0*HRR(18))+  &
        HRR(33)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+5)*LDB+CDOffSet !=(9,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(9)+  &
        2.D0*HRR(19))+  &
        HRR(34)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+5)*LDB+CDOffSet !=(10,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(10)+  &
        2.D0*HRR(20))+  &
        HRR(35)+  &
        INTGRL(OffSet)
END SUBROUTINE BraHRR66
