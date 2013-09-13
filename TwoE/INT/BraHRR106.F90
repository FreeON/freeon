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
   SUBROUTINE BraHRR106(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(11,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(11)+  &
        2.D0*HRR(21))+  &
        HRR(36)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(12,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(12)+  &
        2.D0*HRR(22))+  &
        HRR(37)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(13,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(13)+  &
        2.D0*HRR(23))+  &
        HRR(38)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=(14,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(14)+  &
        2.D0*HRR(24))+  &
        HRR(39)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=(15,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(15)+  &
        2.D0*HRR(26))+  &
        HRR(42)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=(16,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(16)+  &
        2.D0*HRR(27))+  &
        HRR(43)+  &
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+0)*LDB+CDOffSet !=(17,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(17)+  &
        2.D0*HRR(28))+  &
        HRR(44)+  &
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+0)*LDB+CDOffSet !=(18,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(18)+  &
        2.D0*HRR(30))+  &
        HRR(47)+  &
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+0)*LDB+CDOffSet !=(19,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(19)+  &
        2.D0*HRR(31))+  &
        HRR(48)+  &
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+0)*LDB+CDOffSet !=(20,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(20)+  &
        2.D0*HRR(33))+  &
        HRR(51)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(11,6|
      INTGRL(OffSet)=ABy*HRR(21)+  &
        ABx*(ABy*HRR(11)+  &
        HRR(22))+  &
        HRR(37)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=(12,6|
      INTGRL(OffSet)=ABy*HRR(22)+  &
        ABx*(ABy*HRR(12)+  &
        HRR(23))+  &
        HRR(38)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=(13,6|
      INTGRL(OffSet)=ABy*HRR(23)+  &
        ABx*(ABy*HRR(13)+  &
        HRR(24))+  &
        HRR(39)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=(14,6|
      INTGRL(OffSet)=ABy*HRR(24)+  &
        ABx*(ABy*HRR(14)+  &
        HRR(25))+  &
        HRR(40)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+1)*LDB+CDOffSet !=(15,6|
      INTGRL(OffSet)=ABy*HRR(26)+  &
        ABx*(ABy*HRR(15)+  &
        HRR(27))+  &
        HRR(43)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+1)*LDB+CDOffSet !=(16,6|
      INTGRL(OffSet)=ABy*HRR(27)+  &
        ABx*(ABy*HRR(16)+  &
        HRR(28))+  &
        HRR(44)+  &
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+1)*LDB+CDOffSet !=(17,6|
      INTGRL(OffSet)=ABy*HRR(28)+  &
        ABx*(ABy*HRR(17)+  &
        HRR(29))+  &
        HRR(45)+  &
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+1)*LDB+CDOffSet !=(18,6|
      INTGRL(OffSet)=ABy*HRR(30)+  &
        ABx*(ABy*HRR(18)+  &
        HRR(31))+  &
        HRR(48)+  &
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+1)*LDB+CDOffSet !=(19,6|
      INTGRL(OffSet)=ABy*HRR(31)+  &
        ABx*(ABy*HRR(19)+  &
        HRR(32))+  &
        HRR(49)+  &
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+1)*LDB+CDOffSet !=(20,6|
      INTGRL(OffSet)=ABy*HRR(33)+  &
        ABx*(ABy*HRR(20)+  &
        HRR(34))+  &
        HRR(52)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(11,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(11)+  &
        2.D0*HRR(22))+  &
        HRR(38)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=(12,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(12)+  &
        2.D0*HRR(23))+  &
        HRR(39)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=(13,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(13)+  &
        2.D0*HRR(24))+  &
        HRR(40)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=(14,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(14)+  &
        2.D0*HRR(25))+  &
        HRR(41)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+2)*LDB+CDOffSet !=(15,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(15)+  &
        2.D0*HRR(27))+  &
        HRR(44)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+2)*LDB+CDOffSet !=(16,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(16)+  &
        2.D0*HRR(28))+  &
        HRR(45)+  &
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+2)*LDB+CDOffSet !=(17,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(17)+  &
        2.D0*HRR(29))+  &
        HRR(46)+  &
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+2)*LDB+CDOffSet !=(18,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(18)+  &
        2.D0*HRR(31))+  &
        HRR(49)+  &
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+2)*LDB+CDOffSet !=(19,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(19)+  &
        2.D0*HRR(32))+  &
        HRR(50)+  &
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+2)*LDB+CDOffSet !=(20,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(20)+  &
        2.D0*HRR(34))+  &
        HRR(53)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=(11,8|
      INTGRL(OffSet)=ABz*HRR(21)+  &
        ABx*(ABz*HRR(11)+  &
        HRR(26))+  &
        HRR(42)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=(12,8|
      INTGRL(OffSet)=ABz*HRR(22)+  &
        ABx*(ABz*HRR(12)+  &
        HRR(27))+  &
        HRR(43)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=(13,8|
      INTGRL(OffSet)=ABz*HRR(23)+  &
        ABx*(ABz*HRR(13)+  &
        HRR(28))+  &
        HRR(44)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+3)*LDB+CDOffSet !=(14,8|
      INTGRL(OffSet)=ABz*HRR(24)+  &
        ABx*(ABz*HRR(14)+  &
        HRR(29))+  &
        HRR(45)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+3)*LDB+CDOffSet !=(15,8|
      INTGRL(OffSet)=ABz*HRR(26)+  &
        ABx*(ABz*HRR(15)+  &
        HRR(30))+  &
        HRR(47)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+3)*LDB+CDOffSet !=(16,8|
      INTGRL(OffSet)=ABz*HRR(27)+  &
        ABx*(ABz*HRR(16)+  &
        HRR(31))+  &
        HRR(48)+  &
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+3)*LDB+CDOffSet !=(17,8|
      INTGRL(OffSet)=ABz*HRR(28)+  &
        ABx*(ABz*HRR(17)+  &
        HRR(32))+  &
        HRR(49)+  &
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+3)*LDB+CDOffSet !=(18,8|
      INTGRL(OffSet)=ABz*HRR(30)+  &
        ABx*(ABz*HRR(18)+  &
        HRR(33))+  &
        HRR(51)+  &
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+3)*LDB+CDOffSet !=(19,8|
      INTGRL(OffSet)=ABz*HRR(31)+  &
        ABx*(ABz*HRR(19)+  &
        HRR(34))+  &
        HRR(52)+  &
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+3)*LDB+CDOffSet !=(20,8|
      INTGRL(OffSet)=ABz*HRR(33)+  &
        ABx*(ABz*HRR(20)+  &
        HRR(35))+  &
        HRR(54)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=(11,9|
      INTGRL(OffSet)=ABz*HRR(22)+  &
        ABy*(ABz*HRR(11)+  &
        HRR(26))+  &
        HRR(43)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+4)*LDB+CDOffSet !=(12,9|
      INTGRL(OffSet)=ABz*HRR(23)+  &
        ABy*(ABz*HRR(12)+  &
        HRR(27))+  &
        HRR(44)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+4)*LDB+CDOffSet !=(13,9|
      INTGRL(OffSet)=ABz*HRR(24)+  &
        ABy*(ABz*HRR(13)+  &
        HRR(28))+  &
        HRR(45)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+4)*LDB+CDOffSet !=(14,9|
      INTGRL(OffSet)=ABz*HRR(25)+  &
        ABy*(ABz*HRR(14)+  &
        HRR(29))+  &
        HRR(46)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+4)*LDB+CDOffSet !=(15,9|
      INTGRL(OffSet)=ABz*HRR(27)+  &
        ABy*(ABz*HRR(15)+  &
        HRR(30))+  &
        HRR(48)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+4)*LDB+CDOffSet !=(16,9|
      INTGRL(OffSet)=ABz*HRR(28)+  &
        ABy*(ABz*HRR(16)+  &
        HRR(31))+  &
        HRR(49)+  &
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+4)*LDB+CDOffSet !=(17,9|
      INTGRL(OffSet)=ABz*HRR(29)+  &
        ABy*(ABz*HRR(17)+  &
        HRR(32))+  &
        HRR(50)+  &
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+4)*LDB+CDOffSet !=(18,9|
      INTGRL(OffSet)=ABz*HRR(31)+  &
        ABy*(ABz*HRR(18)+  &
        HRR(33))+  &
        HRR(52)+  &
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+4)*LDB+CDOffSet !=(19,9|
      INTGRL(OffSet)=ABz*HRR(32)+  &
        ABy*(ABz*HRR(19)+  &
        HRR(34))+  &
        HRR(53)+  &
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+4)*LDB+CDOffSet !=(20,9|
      INTGRL(OffSet)=ABz*HRR(34)+  &
        ABy*(ABz*HRR(20)+  &
        HRR(35))+  &
        HRR(55)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=(11,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(11)+  &
        2.D0*HRR(26))+  &
        HRR(47)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+5)*LDB+CDOffSet !=(12,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(12)+  &
        2.D0*HRR(27))+  &
        HRR(48)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+5)*LDB+CDOffSet !=(13,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(13)+  &
        2.D0*HRR(28))+  &
        HRR(49)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+5)*LDB+CDOffSet !=(14,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(14)+  &
        2.D0*HRR(29))+  &
        HRR(50)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+5)*LDB+CDOffSet !=(15,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(15)+  &
        2.D0*HRR(30))+  &
        HRR(51)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+5)*LDB+CDOffSet !=(16,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(16)+  &
        2.D0*HRR(31))+  &
        HRR(52)+  &
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+5)*LDB+CDOffSet !=(17,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(17)+  &
        2.D0*HRR(32))+  &
        HRR(53)+  &
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+5)*LDB+CDOffSet !=(18,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(18)+  &
        2.D0*HRR(33))+  &
        HRR(54)+  &
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+5)*LDB+CDOffSet !=(19,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(19)+  &
        2.D0*HRR(34))+  &
        HRR(55)+  &
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+5)*LDB+CDOffSet !=(20,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(20)+  &
        2.D0*HRR(35))+  &
        HRR(56)+  &
        INTGRL(OffSet)
END SUBROUTINE BraHRR106
