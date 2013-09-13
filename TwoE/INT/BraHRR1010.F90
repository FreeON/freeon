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
   SUBROUTINE BraHRR1010(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(11,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(11)+  &
        3.D0*HRR(21))+  &
        3.D0*HRR(36))+  &
        HRR(57)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(12,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(12)+  &
        3.D0*HRR(22))+  &
        3.D0*HRR(37))+  &
        HRR(58)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(13,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(13)+  &
        3.D0*HRR(23))+  &
        3.D0*HRR(38))+  &
        HRR(59)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=(14,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(14)+  &
        3.D0*HRR(24))+  &
        3.D0*HRR(39))+  &
        HRR(60)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=(15,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(15)+  &
        3.D0*HRR(26))+  &
        3.D0*HRR(42))+  &
        HRR(64)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=(16,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(16)+  &
        3.D0*HRR(27))+  &
        3.D0*HRR(43))+  &
        HRR(65)+  &
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+0)*LDB+CDOffSet !=(17,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(17)+  &
        3.D0*HRR(28))+  &
        3.D0*HRR(44))+  &
        HRR(66)+  &
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+0)*LDB+CDOffSet !=(18,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(18)+  &
        3.D0*HRR(30))+  &
        3.D0*HRR(47))+  &
        HRR(70)+  &
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+0)*LDB+CDOffSet !=(19,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(19)+  &
        3.D0*HRR(31))+  &
        3.D0*HRR(48))+  &
        HRR(71)+  &
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+0)*LDB+CDOffSet !=(20,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(20)+  &
        3.D0*HRR(33))+  &
        3.D0*HRR(51))+  &
        HRR(75)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(11,12|
      INTGRL(OffSet)=ABy*HRR(36)+  &
        ABx*(2.D0*ABy*HRR(21)+  &
        ABx*(ABy*HRR(11)+  &
        HRR(22))+  &
        2.D0*HRR(37))+  &
        HRR(58)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=(12,12|
      INTGRL(OffSet)=ABy*HRR(37)+  &
        ABx*(2.D0*ABy*HRR(22)+  &
        ABx*(ABy*HRR(12)+  &
        HRR(23))+  &
        2.D0*HRR(38))+  &
        HRR(59)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=(13,12|
      INTGRL(OffSet)=ABy*HRR(38)+  &
        ABx*(2.D0*ABy*HRR(23)+  &
        ABx*(ABy*HRR(13)+  &
        HRR(24))+  &
        2.D0*HRR(39))+  &
        HRR(60)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=(14,12|
      INTGRL(OffSet)=ABy*HRR(39)+  &
        ABx*(2.D0*ABy*HRR(24)+  &
        ABx*(ABy*HRR(14)+  &
        HRR(25))+  &
        2.D0*HRR(40))+  &
        HRR(61)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+1)*LDB+CDOffSet !=(15,12|
      INTGRL(OffSet)=ABy*HRR(42)+  &
        ABx*(2.D0*ABy*HRR(26)+  &
        ABx*(ABy*HRR(15)+  &
        HRR(27))+  &
        2.D0*HRR(43))+  &
        HRR(65)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+1)*LDB+CDOffSet !=(16,12|
      INTGRL(OffSet)=ABy*HRR(43)+  &
        ABx*(2.D0*ABy*HRR(27)+  &
        ABx*(ABy*HRR(16)+  &
        HRR(28))+  &
        2.D0*HRR(44))+  &
        HRR(66)+  &
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+1)*LDB+CDOffSet !=(17,12|
      INTGRL(OffSet)=ABy*HRR(44)+  &
        ABx*(2.D0*ABy*HRR(28)+  &
        ABx*(ABy*HRR(17)+  &
        HRR(29))+  &
        2.D0*HRR(45))+  &
        HRR(67)+  &
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+1)*LDB+CDOffSet !=(18,12|
      INTGRL(OffSet)=ABy*HRR(47)+  &
        ABx*(2.D0*ABy*HRR(30)+  &
        ABx*(ABy*HRR(18)+  &
        HRR(31))+  &
        2.D0*HRR(48))+  &
        HRR(71)+  &
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+1)*LDB+CDOffSet !=(19,12|
      INTGRL(OffSet)=ABy*HRR(48)+  &
        ABx*(2.D0*ABy*HRR(31)+  &
        ABx*(ABy*HRR(19)+  &
        HRR(32))+  &
        2.D0*HRR(49))+  &
        HRR(72)+  &
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+1)*LDB+CDOffSet !=(20,12|
      INTGRL(OffSet)=ABy*HRR(51)+  &
        ABx*(2.D0*ABy*HRR(33)+  &
        ABx*(ABy*HRR(20)+  &
        HRR(34))+  &
        2.D0*HRR(52))+  &
        HRR(76)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(11,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(21)+  &
        2.D0*HRR(37))+  &
        ABx*(ABy*(ABy*HRR(11)+  &
        2.D0*HRR(22))+  &
        HRR(38))+  &
        HRR(59)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=(12,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(22)+  &
        2.D0*HRR(38))+  &
        ABx*(ABy*(ABy*HRR(12)+  &
        2.D0*HRR(23))+  &
        HRR(39))+  &
        HRR(60)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=(13,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(23)+  &
        2.D0*HRR(39))+  &
        ABx*(ABy*(ABy*HRR(13)+  &
        2.D0*HRR(24))+  &
        HRR(40))+  &
        HRR(61)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=(14,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(24)+  &
        2.D0*HRR(40))+  &
        ABx*(ABy*(ABy*HRR(14)+  &
        2.D0*HRR(25))+  &
        HRR(41))+  &
        HRR(62)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+2)*LDB+CDOffSet !=(15,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(26)+  &
        2.D0*HRR(43))+  &
        ABx*(ABy*(ABy*HRR(15)+  &
        2.D0*HRR(27))+  &
        HRR(44))+  &
        HRR(66)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+2)*LDB+CDOffSet !=(16,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(27)+  &
        2.D0*HRR(44))+  &
        ABx*(ABy*(ABy*HRR(16)+  &
        2.D0*HRR(28))+  &
        HRR(45))+  &
        HRR(67)+  &
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+2)*LDB+CDOffSet !=(17,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(28)+  &
        2.D0*HRR(45))+  &
        ABx*(ABy*(ABy*HRR(17)+  &
        2.D0*HRR(29))+  &
        HRR(46))+  &
        HRR(68)+  &
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+2)*LDB+CDOffSet !=(18,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(30)+  &
        2.D0*HRR(48))+  &
        ABx*(ABy*(ABy*HRR(18)+  &
        2.D0*HRR(31))+  &
        HRR(49))+  &
        HRR(72)+  &
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+2)*LDB+CDOffSet !=(19,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(31)+  &
        2.D0*HRR(49))+  &
        ABx*(ABy*(ABy*HRR(19)+  &
        2.D0*HRR(32))+  &
        HRR(50))+  &
        HRR(73)+  &
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+2)*LDB+CDOffSet !=(20,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(33)+  &
        2.D0*HRR(52))+  &
        ABx*(ABy*(ABy*HRR(20)+  &
        2.D0*HRR(34))+  &
        HRR(53))+  &
        HRR(77)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=(11,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(11)+  &
        3.D0*HRR(22))+  &
        3.D0*HRR(38))+  &
        HRR(60)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=(12,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(12)+  &
        3.D0*HRR(23))+  &
        3.D0*HRR(39))+  &
        HRR(61)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=(13,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(13)+  &
        3.D0*HRR(24))+  &
        3.D0*HRR(40))+  &
        HRR(62)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+3)*LDB+CDOffSet !=(14,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(14)+  &
        3.D0*HRR(25))+  &
        3.D0*HRR(41))+  &
        HRR(63)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+3)*LDB+CDOffSet !=(15,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(15)+  &
        3.D0*HRR(27))+  &
        3.D0*HRR(44))+  &
        HRR(67)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+3)*LDB+CDOffSet !=(16,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(16)+  &
        3.D0*HRR(28))+  &
        3.D0*HRR(45))+  &
        HRR(68)+  &
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+3)*LDB+CDOffSet !=(17,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(17)+  &
        3.D0*HRR(29))+  &
        3.D0*HRR(46))+  &
        HRR(69)+  &
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+3)*LDB+CDOffSet !=(18,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(18)+  &
        3.D0*HRR(31))+  &
        3.D0*HRR(49))+  &
        HRR(73)+  &
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+3)*LDB+CDOffSet !=(19,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(19)+  &
        3.D0*HRR(32))+  &
        3.D0*HRR(50))+  &
        HRR(74)+  &
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+3)*LDB+CDOffSet !=(20,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(20)+  &
        3.D0*HRR(34))+  &
        3.D0*HRR(53))+  &
        HRR(78)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=(11,15|
      INTGRL(OffSet)=ABz*HRR(36)+  &
        ABx*(2.D0*ABz*HRR(21)+  &
        ABx*(ABz*HRR(11)+  &
        HRR(26))+  &
        2.D0*HRR(42))+  &
        HRR(64)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+4)*LDB+CDOffSet !=(12,15|
      INTGRL(OffSet)=ABz*HRR(37)+  &
        ABx*(2.D0*ABz*HRR(22)+  &
        ABx*(ABz*HRR(12)+  &
        HRR(27))+  &
        2.D0*HRR(43))+  &
        HRR(65)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+4)*LDB+CDOffSet !=(13,15|
      INTGRL(OffSet)=ABz*HRR(38)+  &
        ABx*(2.D0*ABz*HRR(23)+  &
        ABx*(ABz*HRR(13)+  &
        HRR(28))+  &
        2.D0*HRR(44))+  &
        HRR(66)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+4)*LDB+CDOffSet !=(14,15|
      INTGRL(OffSet)=ABz*HRR(39)+  &
        ABx*(2.D0*ABz*HRR(24)+  &
        ABx*(ABz*HRR(14)+  &
        HRR(29))+  &
        2.D0*HRR(45))+  &
        HRR(67)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+4)*LDB+CDOffSet !=(15,15|
      INTGRL(OffSet)=ABz*HRR(42)+  &
        ABx*(2.D0*ABz*HRR(26)+  &
        ABx*(ABz*HRR(15)+  &
        HRR(30))+  &
        2.D0*HRR(47))+  &
        HRR(70)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+4)*LDB+CDOffSet !=(16,15|
      INTGRL(OffSet)=ABz*HRR(43)+  &
        ABx*(2.D0*ABz*HRR(27)+  &
        ABx*(ABz*HRR(16)+  &
        HRR(31))+  &
        2.D0*HRR(48))+  &
        HRR(71)+  &
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+4)*LDB+CDOffSet !=(17,15|
      INTGRL(OffSet)=ABz*HRR(44)+  &
        ABx*(2.D0*ABz*HRR(28)+  &
        ABx*(ABz*HRR(17)+  &
        HRR(32))+  &
        2.D0*HRR(49))+  &
        HRR(72)+  &
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+4)*LDB+CDOffSet !=(18,15|
      INTGRL(OffSet)=ABz*HRR(47)+  &
        ABx*(2.D0*ABz*HRR(30)+  &
        ABx*(ABz*HRR(18)+  &
        HRR(33))+  &
        2.D0*HRR(51))+  &
        HRR(75)+  &
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+4)*LDB+CDOffSet !=(19,15|
      INTGRL(OffSet)=ABz*HRR(48)+  &
        ABx*(2.D0*ABz*HRR(31)+  &
        ABx*(ABz*HRR(19)+  &
        HRR(34))+  &
        2.D0*HRR(52))+  &
        HRR(76)+  &
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+4)*LDB+CDOffSet !=(20,15|
      INTGRL(OffSet)=ABz*HRR(51)+  &
        ABx*(2.D0*ABz*HRR(33)+  &
        ABx*(ABz*HRR(20)+  &
        HRR(35))+  &
        2.D0*HRR(54))+  &
        HRR(79)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=(11,16|
      INTGRL(OffSet)=ABz*HRR(37)+  &
        ABy*(ABz*HRR(21)+  &
        HRR(42))+  &
        ABx*(ABz*HRR(22)+  &
        ABy*(ABz*HRR(11)+  &
        HRR(26))+  &
        HRR(43))+  &
        HRR(65)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+5)*LDB+CDOffSet !=(12,16|
      INTGRL(OffSet)=ABz*HRR(38)+  &
        ABy*(ABz*HRR(22)+  &
        HRR(43))+  &
        ABx*(ABz*HRR(23)+  &
        ABy*(ABz*HRR(12)+  &
        HRR(27))+  &
        HRR(44))+  &
        HRR(66)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+5)*LDB+CDOffSet !=(13,16|
      INTGRL(OffSet)=ABz*HRR(39)+  &
        ABy*(ABz*HRR(23)+  &
        HRR(44))+  &
        ABx*(ABz*HRR(24)+  &
        ABy*(ABz*HRR(13)+  &
        HRR(28))+  &
        HRR(45))+  &
        HRR(67)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+5)*LDB+CDOffSet !=(14,16|
      INTGRL(OffSet)=ABz*HRR(40)+  &
        ABy*(ABz*HRR(24)+  &
        HRR(45))+  &
        ABx*(ABz*HRR(25)+  &
        ABy*(ABz*HRR(14)+  &
        HRR(29))+  &
        HRR(46))+  &
        HRR(68)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+5)*LDB+CDOffSet !=(15,16|
      INTGRL(OffSet)=ABz*HRR(43)+  &
        ABy*(ABz*HRR(26)+  &
        HRR(47))+  &
        ABx*(ABz*HRR(27)+  &
        ABy*(ABz*HRR(15)+  &
        HRR(30))+  &
        HRR(48))+  &
        HRR(71)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+5)*LDB+CDOffSet !=(16,16|
      INTGRL(OffSet)=ABz*HRR(44)+  &
        ABy*(ABz*HRR(27)+  &
        HRR(48))+  &
        ABx*(ABz*HRR(28)+  &
        ABy*(ABz*HRR(16)+  &
        HRR(31))+  &
        HRR(49))+  &
        HRR(72)+  &
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+5)*LDB+CDOffSet !=(17,16|
      INTGRL(OffSet)=ABz*HRR(45)+  &
        ABy*(ABz*HRR(28)+  &
        HRR(49))+  &
        ABx*(ABz*HRR(29)+  &
        ABy*(ABz*HRR(17)+  &
        HRR(32))+  &
        HRR(50))+  &
        HRR(73)+  &
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+5)*LDB+CDOffSet !=(18,16|
      INTGRL(OffSet)=ABz*HRR(48)+  &
        ABy*(ABz*HRR(30)+  &
        HRR(51))+  &
        ABx*(ABz*HRR(31)+  &
        ABy*(ABz*HRR(18)+  &
        HRR(33))+  &
        HRR(52))+  &
        HRR(76)+  &
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+5)*LDB+CDOffSet !=(19,16|
      INTGRL(OffSet)=ABz*HRR(49)+  &
        ABy*(ABz*HRR(31)+  &
        HRR(52))+  &
        ABx*(ABz*HRR(32)+  &
        ABy*(ABz*HRR(19)+  &
        HRR(34))+  &
        HRR(53))+  &
        HRR(77)+  &
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+5)*LDB+CDOffSet !=(20,16|
      INTGRL(OffSet)=ABz*HRR(52)+  &
        ABy*(ABz*HRR(33)+  &
        HRR(54))+  &
        ABx*(ABz*HRR(34)+  &
        ABy*(ABz*HRR(20)+  &
        HRR(35))+  &
        HRR(55))+  &
        HRR(80)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+6)*LDB+CDOffSet !=(11,17|
      INTGRL(OffSet)=ABz*HRR(38)+  &
        ABy*(2.D0*ABz*HRR(22)+  &
        ABy*(ABz*HRR(11)+  &
        HRR(26))+  &
        2.D0*HRR(43))+  &
        HRR(66)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+6)*LDB+CDOffSet !=(12,17|
      INTGRL(OffSet)=ABz*HRR(39)+  &
        ABy*(2.D0*ABz*HRR(23)+  &
        ABy*(ABz*HRR(12)+  &
        HRR(27))+  &
        2.D0*HRR(44))+  &
        HRR(67)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+6)*LDB+CDOffSet !=(13,17|
      INTGRL(OffSet)=ABz*HRR(40)+  &
        ABy*(2.D0*ABz*HRR(24)+  &
        ABy*(ABz*HRR(13)+  &
        HRR(28))+  &
        2.D0*HRR(45))+  &
        HRR(68)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+6)*LDB+CDOffSet !=(14,17|
      INTGRL(OffSet)=ABz*HRR(41)+  &
        ABy*(2.D0*ABz*HRR(25)+  &
        ABy*(ABz*HRR(14)+  &
        HRR(29))+  &
        2.D0*HRR(46))+  &
        HRR(69)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+6)*LDB+CDOffSet !=(15,17|
      INTGRL(OffSet)=ABz*HRR(44)+  &
        ABy*(2.D0*ABz*HRR(27)+  &
        ABy*(ABz*HRR(15)+  &
        HRR(30))+  &
        2.D0*HRR(48))+  &
        HRR(72)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+6)*LDB+CDOffSet !=(16,17|
      INTGRL(OffSet)=ABz*HRR(45)+  &
        ABy*(2.D0*ABz*HRR(28)+  &
        ABy*(ABz*HRR(16)+  &
        HRR(31))+  &
        2.D0*HRR(49))+  &
        HRR(73)+  &
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+6)*LDB+CDOffSet !=(17,17|
      INTGRL(OffSet)=ABz*HRR(46)+  &
        ABy*(2.D0*ABz*HRR(29)+  &
        ABy*(ABz*HRR(17)+  &
        HRR(32))+  &
        2.D0*HRR(50))+  &
        HRR(74)+  &
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+6)*LDB+CDOffSet !=(18,17|
      INTGRL(OffSet)=ABz*HRR(49)+  &
        ABy*(2.D0*ABz*HRR(31)+  &
        ABy*(ABz*HRR(18)+  &
        HRR(33))+  &
        2.D0*HRR(52))+  &
        HRR(77)+  &
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+6)*LDB+CDOffSet !=(19,17|
      INTGRL(OffSet)=ABz*HRR(50)+  &
        ABy*(2.D0*ABz*HRR(32)+  &
        ABy*(ABz*HRR(19)+  &
        HRR(34))+  &
        2.D0*HRR(53))+  &
        HRR(78)+  &
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+6)*LDB+CDOffSet !=(20,17|
      INTGRL(OffSet)=ABz*HRR(53)+  &
        ABy*(2.D0*ABz*HRR(34)+  &
        ABy*(ABz*HRR(20)+  &
        HRR(35))+  &
        2.D0*HRR(55))+  &
        HRR(81)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+7)*LDB+CDOffSet !=(11,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(21)+  &
        2.D0*HRR(42))+  &
        ABx*(ABz*(ABz*HRR(11)+  &
        2.D0*HRR(26))+  &
        HRR(47))+  &
        HRR(70)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+7)*LDB+CDOffSet !=(12,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(22)+  &
        2.D0*HRR(43))+  &
        ABx*(ABz*(ABz*HRR(12)+  &
        2.D0*HRR(27))+  &
        HRR(48))+  &
        HRR(71)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+7)*LDB+CDOffSet !=(13,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(23)+  &
        2.D0*HRR(44))+  &
        ABx*(ABz*(ABz*HRR(13)+  &
        2.D0*HRR(28))+  &
        HRR(49))+  &
        HRR(72)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+7)*LDB+CDOffSet !=(14,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(24)+  &
        2.D0*HRR(45))+  &
        ABx*(ABz*(ABz*HRR(14)+  &
        2.D0*HRR(29))+  &
        HRR(50))+  &
        HRR(73)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+7)*LDB+CDOffSet !=(15,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(26)+  &
        2.D0*HRR(47))+  &
        ABx*(ABz*(ABz*HRR(15)+  &
        2.D0*HRR(30))+  &
        HRR(51))+  &
        HRR(75)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+7)*LDB+CDOffSet !=(16,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(27)+  &
        2.D0*HRR(48))+  &
        ABx*(ABz*(ABz*HRR(16)+  &
        2.D0*HRR(31))+  &
        HRR(52))+  &
        HRR(76)+  &
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+7)*LDB+CDOffSet !=(17,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(28)+  &
        2.D0*HRR(49))+  &
        ABx*(ABz*(ABz*HRR(17)+  &
        2.D0*HRR(32))+  &
        HRR(53))+  &
        HRR(77)+  &
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+7)*LDB+CDOffSet !=(18,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(30)+  &
        2.D0*HRR(51))+  &
        ABx*(ABz*(ABz*HRR(18)+  &
        2.D0*HRR(33))+  &
        HRR(54))+  &
        HRR(79)+  &
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+7)*LDB+CDOffSet !=(19,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(31)+  &
        2.D0*HRR(52))+  &
        ABx*(ABz*(ABz*HRR(19)+  &
        2.D0*HRR(34))+  &
        HRR(55))+  &
        HRR(80)+  &
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+7)*LDB+CDOffSet !=(20,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(33)+  &
        2.D0*HRR(54))+  &
        ABx*(ABz*(ABz*HRR(20)+  &
        2.D0*HRR(35))+  &
        HRR(56))+  &
        HRR(82)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+8)*LDB+CDOffSet !=(11,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(22)+  &
        2.D0*HRR(43))+  &
        ABy*(ABz*(ABz*HRR(11)+  &
        2.D0*HRR(26))+  &
        HRR(47))+  &
        HRR(71)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+8)*LDB+CDOffSet !=(12,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(23)+  &
        2.D0*HRR(44))+  &
        ABy*(ABz*(ABz*HRR(12)+  &
        2.D0*HRR(27))+  &
        HRR(48))+  &
        HRR(72)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+8)*LDB+CDOffSet !=(13,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(24)+  &
        2.D0*HRR(45))+  &
        ABy*(ABz*(ABz*HRR(13)+  &
        2.D0*HRR(28))+  &
        HRR(49))+  &
        HRR(73)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+8)*LDB+CDOffSet !=(14,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(25)+  &
        2.D0*HRR(46))+  &
        ABy*(ABz*(ABz*HRR(14)+  &
        2.D0*HRR(29))+  &
        HRR(50))+  &
        HRR(74)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+8)*LDB+CDOffSet !=(15,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(27)+  &
        2.D0*HRR(48))+  &
        ABy*(ABz*(ABz*HRR(15)+  &
        2.D0*HRR(30))+  &
        HRR(51))+  &
        HRR(76)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+8)*LDB+CDOffSet !=(16,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(28)+  &
        2.D0*HRR(49))+  &
        ABy*(ABz*(ABz*HRR(16)+  &
        2.D0*HRR(31))+  &
        HRR(52))+  &
        HRR(77)+  &
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+8)*LDB+CDOffSet !=(17,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(29)+  &
        2.D0*HRR(50))+  &
        ABy*(ABz*(ABz*HRR(17)+  &
        2.D0*HRR(32))+  &
        HRR(53))+  &
        HRR(78)+  &
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+8)*LDB+CDOffSet !=(18,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(31)+  &
        2.D0*HRR(52))+  &
        ABy*(ABz*(ABz*HRR(18)+  &
        2.D0*HRR(33))+  &
        HRR(54))+  &
        HRR(80)+  &
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+8)*LDB+CDOffSet !=(19,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(32)+  &
        2.D0*HRR(53))+  &
        ABy*(ABz*(ABz*HRR(19)+  &
        2.D0*HRR(34))+  &
        HRR(55))+  &
        HRR(81)+  &
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+8)*LDB+CDOffSet !=(20,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(34)+  &
        2.D0*HRR(55))+  &
        ABy*(ABz*(ABz*HRR(20)+  &
        2.D0*HRR(35))+  &
        HRR(56))+  &
        HRR(83)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+9)*LDB+CDOffSet !=(11,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(11)+  &
        3.D0*HRR(26))+  &
        3.D0*HRR(47))+  &
        HRR(75)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+9)*LDB+CDOffSet !=(12,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(12)+  &
        3.D0*HRR(27))+  &
        3.D0*HRR(48))+  &
        HRR(76)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+9)*LDB+CDOffSet !=(13,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(13)+  &
        3.D0*HRR(28))+  &
        3.D0*HRR(49))+  &
        HRR(77)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+9)*LDB+CDOffSet !=(14,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(14)+  &
        3.D0*HRR(29))+  &
        3.D0*HRR(50))+  &
        HRR(78)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+9)*LDB+CDOffSet !=(15,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(15)+  &
        3.D0*HRR(30))+  &
        3.D0*HRR(51))+  &
        HRR(79)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+9)*LDB+CDOffSet !=(16,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(16)+  &
        3.D0*HRR(31))+  &
        3.D0*HRR(52))+  &
        HRR(80)+  &
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+9)*LDB+CDOffSet !=(17,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(17)+  &
        3.D0*HRR(32))+  &
        3.D0*HRR(53))+  &
        HRR(81)+  &
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+9)*LDB+CDOffSet !=(18,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(18)+  &
        3.D0*HRR(33))+  &
        3.D0*HRR(54))+  &
        HRR(82)+  &
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+9)*LDB+CDOffSet !=(19,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(19)+  &
        3.D0*HRR(34))+  &
        3.D0*HRR(55))+  &
        HRR(83)+  &
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+9)*LDB+CDOffSet !=(20,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(20)+  &
        3.D0*HRR(35))+  &
        3.D0*HRR(56))+  &
        HRR(84)+  &
        INTGRL(OffSet)
END SUBROUTINE BraHRR1010
