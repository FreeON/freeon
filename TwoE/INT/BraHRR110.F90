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
   SUBROUTINE BraHRR110(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(1,11|
      INTGRL(OffSet)=ABx*(ABx*(ABx*HRR(1)+  &
        3.D0*HRR(2))+  &
        3.D0*HRR(5))+  &
        HRR(11)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(1,12|
      INTGRL(OffSet)=ABy*HRR(5)+  &
        ABx*(2.D0*ABy*HRR(2)+  &
        ABx*(ABy*HRR(1)+  &
        HRR(3))+  &
        2.D0*HRR(6))+  &
        HRR(12)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(1,13|
      INTGRL(OffSet)=ABy*(ABy*HRR(2)+  &
        2.D0*HRR(6))+  &
        ABx*(ABy*(ABy*HRR(1)+  &
        2.D0*HRR(3))+  &
        HRR(7))+  &
        HRR(13)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=(1,14|
      INTGRL(OffSet)=ABy*(ABy*(ABy*HRR(1)+  &
        3.D0*HRR(3))+  &
        3.D0*HRR(7))+  &
        HRR(14)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=(1,15|
      INTGRL(OffSet)=ABz*HRR(5)+  &
        ABx*(2.D0*ABz*HRR(2)+  &
        ABx*(ABz*HRR(1)+  &
        HRR(4))+  &
        2.D0*HRR(8))+  &
        HRR(15)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=(1,16|
      INTGRL(OffSet)=ABz*HRR(6)+  &
        ABy*(ABz*HRR(2)+  &
        HRR(8))+  &
        ABx*(ABz*HRR(3)+  &
        ABy*(ABz*HRR(1)+  &
        HRR(4))+  &
        HRR(9))+  &
        HRR(16)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+6)*LDB+CDOffSet !=(1,17|
      INTGRL(OffSet)=ABz*HRR(7)+  &
        ABy*(2.D0*ABz*HRR(3)+  &
        ABy*(ABz*HRR(1)+  &
        HRR(4))+  &
        2.D0*HRR(9))+  &
        HRR(17)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+7)*LDB+CDOffSet !=(1,18|
      INTGRL(OffSet)=ABz*(ABz*HRR(2)+  &
        2.D0*HRR(8))+  &
        ABx*(ABz*(ABz*HRR(1)+  &
        2.D0*HRR(4))+  &
        HRR(10))+  &
        HRR(18)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+8)*LDB+CDOffSet !=(1,19|
      INTGRL(OffSet)=ABz*(ABz*HRR(3)+  &
        2.D0*HRR(9))+  &
        ABy*(ABz*(ABz*HRR(1)+  &
        2.D0*HRR(4))+  &
        HRR(10))+  &
        HRR(19)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+9)*LDB+CDOffSet !=(1,20|
      INTGRL(OffSet)=ABz*(ABz*(ABz*HRR(1)+  &
        3.D0*HRR(4))+  &
        3.D0*HRR(10))+  &
        HRR(20)+  &
        INTGRL(OffSet)
END SUBROUTINE BraHRR110
