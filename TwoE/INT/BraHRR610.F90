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
   SUBROUTINE BraHRR610(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(5,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(5)+  &
        3.D0*HRR(11))+  &
        3.D0*HRR(21))+  &
        HRR(36)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(6,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(6)+  &
        3.D0*HRR(12))+  &
        3.D0*HRR(22))+  &
        HRR(37)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(7,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(7)+  &
        3.D0*HRR(13))+  &
        3.D0*HRR(23))+  &
        HRR(38)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=(8,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(8)+  &
        3.D0*HRR(15))+  &
        3.D0*HRR(26))+  &
        HRR(42)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=(9,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(9)+  &
        3.D0*HRR(16))+  &
        3.D0*HRR(27))+  &
        HRR(43)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=(10,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(10)+  &
        3.D0*HRR(18))+  &
        3.D0*HRR(30))+  &
        HRR(47)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(5,12|
      INTGRL(OffSet)=ABy*HRR(21)+  &
        ABx*(2.D0*ABy*HRR(11)+  &
        ABx*(ABy*HRR(5)+  &
        HRR(12))+  &
        2.D0*HRR(22))+  &
        HRR(37)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=(6,12|
      INTGRL(OffSet)=ABy*HRR(22)+  &
        ABx*(2.D0*ABy*HRR(12)+  &
        ABx*(ABy*HRR(6)+  &
        HRR(13))+  &
        2.D0*HRR(23))+  &
        HRR(38)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=(7,12|
      INTGRL(OffSet)=ABy*HRR(23)+  &
        ABx*(2.D0*ABy*HRR(13)+  &
        ABx*(ABy*HRR(7)+  &
        HRR(14))+  &
        2.D0*HRR(24))+  &
        HRR(39)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=(8,12|
      INTGRL(OffSet)=ABy*HRR(26)+  &
        ABx*(2.D0*ABy*HRR(15)+  &
        ABx*(ABy*HRR(8)+  &
        HRR(16))+  &
        2.D0*HRR(27))+  &
        HRR(43)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+1)*LDB+CDOffSet !=(9,12|
      INTGRL(OffSet)=ABy*HRR(27)+  &
        ABx*(2.D0*ABy*HRR(16)+  &
        ABx*(ABy*HRR(9)+  &
        HRR(17))+  &
        2.D0*HRR(28))+  &
        HRR(44)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+1)*LDB+CDOffSet !=(10,12|
      INTGRL(OffSet)=ABy*HRR(30)+  &
        ABx*(2.D0*ABy*HRR(18)+  &
        ABx*(ABy*HRR(10)+  &
        HRR(19))+  &
        2.D0*HRR(31))+  &
        HRR(48)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(5,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(11)+  &
        2.D0*HRR(22))+  &
        ABx*(ABy*(ABy*HRR(5)+  &
        2.D0*HRR(12))+  &
        HRR(23))+  &
        HRR(38)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=(6,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(12)+  &
        2.D0*HRR(23))+  &
        ABx*(ABy*(ABy*HRR(6)+  &
        2.D0*HRR(13))+  &
        HRR(24))+  &
        HRR(39)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=(7,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(13)+  &
        2.D0*HRR(24))+  &
        ABx*(ABy*(ABy*HRR(7)+  &
        2.D0*HRR(14))+  &
        HRR(25))+  &
        HRR(40)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=(8,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(15)+  &
        2.D0*HRR(27))+  &
        ABx*(ABy*(ABy*HRR(8)+  &
        2.D0*HRR(16))+  &
        HRR(28))+  &
        HRR(44)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+2)*LDB+CDOffSet !=(9,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(16)+  &
        2.D0*HRR(28))+  &
        ABx*(ABy*(ABy*HRR(9)+  &
        2.D0*HRR(17))+  &
        HRR(29))+  &
        HRR(45)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+2)*LDB+CDOffSet !=(10,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(18)+  &
        2.D0*HRR(31))+  &
        ABx*(ABy*(ABy*HRR(10)+  &
        2.D0*HRR(19))+  &
        HRR(32))+  &
        HRR(49)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=(5,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(5)+  &
        3.D0*HRR(12))+  &
        3.D0*HRR(23))+  &
        HRR(39)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=(6,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(6)+  &
        3.D0*HRR(13))+  &
        3.D0*HRR(24))+  &
        HRR(40)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=(7,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(7)+  &
        3.D0*HRR(14))+  &
        3.D0*HRR(25))+  &
        HRR(41)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+3)*LDB+CDOffSet !=(8,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(8)+  &
        3.D0*HRR(16))+  &
        3.D0*HRR(28))+  &
        HRR(45)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+3)*LDB+CDOffSet !=(9,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(9)+  &
        3.D0*HRR(17))+  &
        3.D0*HRR(29))+  &
        HRR(46)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+3)*LDB+CDOffSet !=(10,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(10)+  &
        3.D0*HRR(19))+  &
        3.D0*HRR(32))+  &
        HRR(50)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=(5,15|
      INTGRL(OffSet)=ABz*HRR(21)+  &
        ABx*(2.D0*ABz*HRR(11)+  &
        ABx*(ABz*HRR(5)+  &
        HRR(15))+  &
        2.D0*HRR(26))+  &
        HRR(42)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+4)*LDB+CDOffSet !=(6,15|
      INTGRL(OffSet)=ABz*HRR(22)+  &
        ABx*(2.D0*ABz*HRR(12)+  &
        ABx*(ABz*HRR(6)+  &
        HRR(16))+  &
        2.D0*HRR(27))+  &
        HRR(43)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+4)*LDB+CDOffSet !=(7,15|
      INTGRL(OffSet)=ABz*HRR(23)+  &
        ABx*(2.D0*ABz*HRR(13)+  &
        ABx*(ABz*HRR(7)+  &
        HRR(17))+  &
        2.D0*HRR(28))+  &
        HRR(44)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+4)*LDB+CDOffSet !=(8,15|
      INTGRL(OffSet)=ABz*HRR(26)+  &
        ABx*(2.D0*ABz*HRR(15)+  &
        ABx*(ABz*HRR(8)+  &
        HRR(18))+  &
        2.D0*HRR(30))+  &
        HRR(47)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+4)*LDB+CDOffSet !=(9,15|
      INTGRL(OffSet)=ABz*HRR(27)+  &
        ABx*(2.D0*ABz*HRR(16)+  &
        ABx*(ABz*HRR(9)+  &
        HRR(19))+  &
        2.D0*HRR(31))+  &
        HRR(48)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+4)*LDB+CDOffSet !=(10,15|
      INTGRL(OffSet)=ABz*HRR(30)+  &
        ABx*(2.D0*ABz*HRR(18)+  &
        ABx*(ABz*HRR(10)+  &
        HRR(20))+  &
        2.D0*HRR(33))+  &
        HRR(51)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=(5,16|
      INTGRL(OffSet)=ABz*HRR(22)+  &
        ABy*(ABz*HRR(11)+  &
        HRR(26))+  &
        ABx*(ABz*HRR(12)+  &
        ABy*(ABz*HRR(5)+  &
        HRR(15))+  &
        HRR(27))+  &
        HRR(43)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+5)*LDB+CDOffSet !=(6,16|
      INTGRL(OffSet)=ABz*HRR(23)+  &
        ABy*(ABz*HRR(12)+  &
        HRR(27))+  &
        ABx*(ABz*HRR(13)+  &
        ABy*(ABz*HRR(6)+  &
        HRR(16))+  &
        HRR(28))+  &
        HRR(44)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+5)*LDB+CDOffSet !=(7,16|
      INTGRL(OffSet)=ABz*HRR(24)+  &
        ABy*(ABz*HRR(13)+  &
        HRR(28))+  &
        ABx*(ABz*HRR(14)+  &
        ABy*(ABz*HRR(7)+  &
        HRR(17))+  &
        HRR(29))+  &
        HRR(45)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+5)*LDB+CDOffSet !=(8,16|
      INTGRL(OffSet)=ABz*HRR(27)+  &
        ABy*(ABz*HRR(15)+  &
        HRR(30))+  &
        ABx*(ABz*HRR(16)+  &
        ABy*(ABz*HRR(8)+  &
        HRR(18))+  &
        HRR(31))+  &
        HRR(48)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+5)*LDB+CDOffSet !=(9,16|
      INTGRL(OffSet)=ABz*HRR(28)+  &
        ABy*(ABz*HRR(16)+  &
        HRR(31))+  &
        ABx*(ABz*HRR(17)+  &
        ABy*(ABz*HRR(9)+  &
        HRR(19))+  &
        HRR(32))+  &
        HRR(49)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+5)*LDB+CDOffSet !=(10,16|
      INTGRL(OffSet)=ABz*HRR(31)+  &
        ABy*(ABz*HRR(18)+  &
        HRR(33))+  &
        ABx*(ABz*HRR(19)+  &
        ABy*(ABz*HRR(10)+  &
        HRR(20))+  &
        HRR(34))+  &
        HRR(52)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+6)*LDB+CDOffSet !=(5,17|
      INTGRL(OffSet)=ABz*HRR(23)+  &
        ABy*(2.D0*ABz*HRR(12)+  &
        ABy*(ABz*HRR(5)+  &
        HRR(15))+  &
        2.D0*HRR(27))+  &
        HRR(44)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+6)*LDB+CDOffSet !=(6,17|
      INTGRL(OffSet)=ABz*HRR(24)+  &
        ABy*(2.D0*ABz*HRR(13)+  &
        ABy*(ABz*HRR(6)+  &
        HRR(16))+  &
        2.D0*HRR(28))+  &
        HRR(45)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+6)*LDB+CDOffSet !=(7,17|
      INTGRL(OffSet)=ABz*HRR(25)+  &
        ABy*(2.D0*ABz*HRR(14)+  &
        ABy*(ABz*HRR(7)+  &
        HRR(17))+  &
        2.D0*HRR(29))+  &
        HRR(46)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+6)*LDB+CDOffSet !=(8,17|
      INTGRL(OffSet)=ABz*HRR(28)+  &
        ABy*(2.D0*ABz*HRR(16)+  &
        ABy*(ABz*HRR(8)+  &
        HRR(18))+  &
        2.D0*HRR(31))+  &
        HRR(49)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+6)*LDB+CDOffSet !=(9,17|
      INTGRL(OffSet)=ABz*HRR(29)+  &
        ABy*(2.D0*ABz*HRR(17)+  &
        ABy*(ABz*HRR(9)+  &
        HRR(19))+  &
        2.D0*HRR(32))+  &
        HRR(50)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+6)*LDB+CDOffSet !=(10,17|
      INTGRL(OffSet)=ABz*HRR(32)+  &
        ABy*(2.D0*ABz*HRR(19)+  &
        ABy*(ABz*HRR(10)+  &
        HRR(20))+  &
        2.D0*HRR(34))+  &
        HRR(53)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+7)*LDB+CDOffSet !=(5,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(11)+  &
        2.D0*HRR(26))+  &
        ABx*(ABz*(ABz*HRR(5)+  &
        2.D0*HRR(15))+  &
        HRR(30))+  &
        HRR(47)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+7)*LDB+CDOffSet !=(6,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(12)+  &
        2.D0*HRR(27))+  &
        ABx*(ABz*(ABz*HRR(6)+  &
        2.D0*HRR(16))+  &
        HRR(31))+  &
        HRR(48)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+7)*LDB+CDOffSet !=(7,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(13)+  &
        2.D0*HRR(28))+  &
        ABx*(ABz*(ABz*HRR(7)+  &
        2.D0*HRR(17))+  &
        HRR(32))+  &
        HRR(49)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+7)*LDB+CDOffSet !=(8,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(15)+  &
        2.D0*HRR(30))+  &
        ABx*(ABz*(ABz*HRR(8)+  &
        2.D0*HRR(18))+  &
        HRR(33))+  &
        HRR(51)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+7)*LDB+CDOffSet !=(9,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(16)+  &
        2.D0*HRR(31))+  &
        ABx*(ABz*(ABz*HRR(9)+  &
        2.D0*HRR(19))+  &
        HRR(34))+  &
        HRR(52)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+7)*LDB+CDOffSet !=(10,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(18)+  &
        2.D0*HRR(33))+  &
        ABx*(ABz*(ABz*HRR(10)+  &
        2.D0*HRR(20))+  &
        HRR(35))+  &
        HRR(54)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+8)*LDB+CDOffSet !=(5,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(12)+  &
        2.D0*HRR(27))+  &
        ABy*(ABz*(ABz*HRR(5)+  &
        2.D0*HRR(15))+  &
        HRR(30))+  &
        HRR(48)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+8)*LDB+CDOffSet !=(6,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(13)+  &
        2.D0*HRR(28))+  &
        ABy*(ABz*(ABz*HRR(6)+  &
        2.D0*HRR(16))+  &
        HRR(31))+  &
        HRR(49)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+8)*LDB+CDOffSet !=(7,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(14)+  &
        2.D0*HRR(29))+  &
        ABy*(ABz*(ABz*HRR(7)+  &
        2.D0*HRR(17))+  &
        HRR(32))+  &
        HRR(50)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+8)*LDB+CDOffSet !=(8,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(16)+  &
        2.D0*HRR(31))+  &
        ABy*(ABz*(ABz*HRR(8)+  &
        2.D0*HRR(18))+  &
        HRR(33))+  &
        HRR(52)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+8)*LDB+CDOffSet !=(9,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(17)+  &
        2.D0*HRR(32))+  &
        ABy*(ABz*(ABz*HRR(9)+  &
        2.D0*HRR(19))+  &
        HRR(34))+  &
        HRR(53)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+8)*LDB+CDOffSet !=(10,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(19)+  &
        2.D0*HRR(34))+  &
        ABy*(ABz*(ABz*HRR(10)+  &
        2.D0*HRR(20))+  &
        HRR(35))+  &
        HRR(55)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+9)*LDB+CDOffSet !=(5,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(5)+  &
        3.D0*HRR(15))+  &
        3.D0*HRR(30))+  &
        HRR(51)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+9)*LDB+CDOffSet !=(6,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(6)+  &
        3.D0*HRR(16))+  &
        3.D0*HRR(31))+  &
        HRR(52)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+9)*LDB+CDOffSet !=(7,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(7)+  &
        3.D0*HRR(17))+  &
        3.D0*HRR(32))+  &
        HRR(53)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+9)*LDB+CDOffSet !=(8,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(8)+  &
        3.D0*HRR(18))+  &
        3.D0*HRR(33))+  &
        HRR(54)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+9)*LDB+CDOffSet !=(9,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(9)+  &
        3.D0*HRR(19))+  &
        3.D0*HRR(34))+  &
        HRR(55)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+9)*LDB+CDOffSet !=(10,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(10)+  &
        3.D0*HRR(20))+  &
        3.D0*HRR(35))+  &
        HRR(56)+  &
        INTGRL(OffSet)
END SUBROUTINE BraHRR610
