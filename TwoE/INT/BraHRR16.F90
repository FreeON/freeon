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
   SUBROUTINE BraHRR16(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(1,5|
      INTGRL(OffSet)=ABx*(ABx*HRR(1)+  &
        2.D0*HRR(2))+  &
        HRR(5)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(1,6|
      INTGRL(OffSet)=ABy*HRR(2)+  &
        ABx*(ABy*HRR(1)+  &
        HRR(3))+  &
        HRR(6)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(1,7|
      INTGRL(OffSet)=ABy*(ABy*HRR(1)+  &
        2.D0*HRR(3))+  &
        HRR(7)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=(1,8|
      INTGRL(OffSet)=ABz*HRR(2)+  &
        ABx*(ABz*HRR(1)+  &
        HRR(4))+  &
        HRR(8)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+4)*LDB+CDOffSet !=(1,9|
      INTGRL(OffSet)=ABz*HRR(3)+  &
        ABy*(ABz*HRR(1)+  &
        HRR(4))+  &
        HRR(9)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+5)*LDB+CDOffSet !=(1,10|
      INTGRL(OffSet)=ABz*(ABz*HRR(1)+  &
        2.D0*HRR(4))+  &
        HRR(10)+  &
        INTGRL(OffSet)
END SUBROUTINE BraHRR16
