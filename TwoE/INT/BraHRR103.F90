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
   SUBROUTINE BraHRR103(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(11,2|
      INTGRL(OffSet)=ABx*HRR(11)+  &
        HRR(21)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(12,2|
      INTGRL(OffSet)=ABx*HRR(12)+  &
        HRR(22)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(13,2|
      INTGRL(OffSet)=ABx*HRR(13)+  &
        HRR(23)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=(14,2|
      INTGRL(OffSet)=ABx*HRR(14)+  &
        HRR(24)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=(15,2|
      INTGRL(OffSet)=ABx*HRR(15)+  &
        HRR(26)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=(16,2|
      INTGRL(OffSet)=ABx*HRR(16)+  &
        HRR(27)+  &
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+0)*LDB+CDOffSet !=(17,2|
      INTGRL(OffSet)=ABx*HRR(17)+  &
        HRR(28)+  &
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+0)*LDB+CDOffSet !=(18,2|
      INTGRL(OffSet)=ABx*HRR(18)+  &
        HRR(30)+  &
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+0)*LDB+CDOffSet !=(19,2|
      INTGRL(OffSet)=ABx*HRR(19)+  &
        HRR(31)+  &
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+0)*LDB+CDOffSet !=(20,2|
      INTGRL(OffSet)=ABx*HRR(20)+  &
        HRR(33)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(11,3|
      INTGRL(OffSet)=ABy*HRR(11)+  &
        HRR(22)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=(12,3|
      INTGRL(OffSet)=ABy*HRR(12)+  &
        HRR(23)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=(13,3|
      INTGRL(OffSet)=ABy*HRR(13)+  &
        HRR(24)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=(14,3|
      INTGRL(OffSet)=ABy*HRR(14)+  &
        HRR(25)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+1)*LDB+CDOffSet !=(15,3|
      INTGRL(OffSet)=ABy*HRR(15)+  &
        HRR(27)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+1)*LDB+CDOffSet !=(16,3|
      INTGRL(OffSet)=ABy*HRR(16)+  &
        HRR(28)+  &
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+1)*LDB+CDOffSet !=(17,3|
      INTGRL(OffSet)=ABy*HRR(17)+  &
        HRR(29)+  &
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+1)*LDB+CDOffSet !=(18,3|
      INTGRL(OffSet)=ABy*HRR(18)+  &
        HRR(31)+  &
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+1)*LDB+CDOffSet !=(19,3|
      INTGRL(OffSet)=ABy*HRR(19)+  &
        HRR(32)+  &
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+1)*LDB+CDOffSet !=(20,3|
      INTGRL(OffSet)=ABy*HRR(20)+  &
        HRR(34)+  &
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(11,4|
      INTGRL(OffSet)=ABz*HRR(11)+  &
        HRR(26)+  &
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=(12,4|
      INTGRL(OffSet)=ABz*HRR(12)+  &
        HRR(27)+  &
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=(13,4|
      INTGRL(OffSet)=ABz*HRR(13)+  &
        HRR(28)+  &
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=(14,4|
      INTGRL(OffSet)=ABz*HRR(14)+  &
        HRR(29)+  &
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+2)*LDB+CDOffSet !=(15,4|
      INTGRL(OffSet)=ABz*HRR(15)+  &
        HRR(30)+  &
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+2)*LDB+CDOffSet !=(16,4|
      INTGRL(OffSet)=ABz*HRR(16)+  &
        HRR(31)+  &
        INTGRL(OffSet)
      OffSet=(OA+6)*LDA+(OB+2)*LDB+CDOffSet !=(17,4|
      INTGRL(OffSet)=ABz*HRR(17)+  &
        HRR(32)+  &
        INTGRL(OffSet)
      OffSet=(OA+7)*LDA+(OB+2)*LDB+CDOffSet !=(18,4|
      INTGRL(OffSet)=ABz*HRR(18)+  &
        HRR(33)+  &
        INTGRL(OffSet)
      OffSet=(OA+8)*LDA+(OB+2)*LDB+CDOffSet !=(19,4|
      INTGRL(OffSet)=ABz*HRR(19)+  &
        HRR(34)+  &
        INTGRL(OffSet)
      OffSet=(OA+9)*LDA+(OB+2)*LDB+CDOffSet !=(20,4|
      INTGRL(OffSet)=ABz*HRR(20)+  &
        HRR(35)+  &
        INTGRL(OffSet)
END SUBROUTINE BraHRR103
