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
   SUBROUTINE BraHRR62(OA,OB,LDA,LDB,CDOffSet,HRR,INTGRL) 
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER       :: OA,OB,LDA,LDB,CDOffSet,OffSet
      REAL(DOUBLE)  :: HRR(*)
      REAL(DOUBLE)  :: INTGRL(*)
      OffSet=(OA+0)*LDA+(OB+0)*LDB+CDOffSet !=(5,1|
      INTGRL(OffSet)=HRR(21)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+0)*LDB+CDOffSet !=(6,1|
      INTGRL(OffSet)=HRR(22)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+0)*LDB+CDOffSet !=(7,1|
      INTGRL(OffSet)=HRR(23)+  & 
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+0)*LDB+CDOffSet !=(8,1|
      INTGRL(OffSet)=HRR(24)+  & 
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+0)*LDB+CDOffSet !=(9,1|
      INTGRL(OffSet)=HRR(25)+  & 
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+0)*LDB+CDOffSet !=(10,1|
      INTGRL(OffSet)=HRR(26)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+1)*LDB+CDOffSet !=(5,2|
      INTGRL(OffSet)=ABx*HRR(5)+  & 
        HRR(11)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+1)*LDB+CDOffSet !=(6,2|
      INTGRL(OffSet)=ABx*HRR(6)+  & 
        HRR(12)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+1)*LDB+CDOffSet !=(7,2|
      INTGRL(OffSet)=ABx*HRR(7)+  & 
        HRR(13)+  & 
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+1)*LDB+CDOffSet !=(8,2|
      INTGRL(OffSet)=ABx*HRR(8)+  & 
        HRR(15)+  & 
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+1)*LDB+CDOffSet !=(9,2|
      INTGRL(OffSet)=ABx*HRR(9)+  & 
        HRR(16)+  & 
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+1)*LDB+CDOffSet !=(10,2|
      INTGRL(OffSet)=ABx*HRR(10)+  & 
        HRR(18)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+2)*LDB+CDOffSet !=(5,3|
      INTGRL(OffSet)=ABy*HRR(5)+  & 
        HRR(12)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+2)*LDB+CDOffSet !=(6,3|
      INTGRL(OffSet)=ABy*HRR(6)+  & 
        HRR(13)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+2)*LDB+CDOffSet !=(7,3|
      INTGRL(OffSet)=ABy*HRR(7)+  & 
        HRR(14)+  & 
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+2)*LDB+CDOffSet !=(8,3|
      INTGRL(OffSet)=ABy*HRR(8)+  & 
        HRR(16)+  & 
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+2)*LDB+CDOffSet !=(9,3|
      INTGRL(OffSet)=ABy*HRR(9)+  & 
        HRR(17)+  & 
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+2)*LDB+CDOffSet !=(10,3|
      INTGRL(OffSet)=ABy*HRR(10)+  & 
        HRR(19)+  & 
        INTGRL(OffSet)
      OffSet=(OA+0)*LDA+(OB+3)*LDB+CDOffSet !=(5,4|
      INTGRL(OffSet)=ABz*HRR(5)+  & 
        HRR(15)+  & 
        INTGRL(OffSet)
      OffSet=(OA+1)*LDA+(OB+3)*LDB+CDOffSet !=(6,4|
      INTGRL(OffSet)=ABz*HRR(6)+  & 
        HRR(16)+  & 
        INTGRL(OffSet)
      OffSet=(OA+2)*LDA+(OB+3)*LDB+CDOffSet !=(7,4|
      INTGRL(OffSet)=ABz*HRR(7)+  & 
        HRR(17)+  & 
        INTGRL(OffSet)
      OffSet=(OA+3)*LDA+(OB+3)*LDB+CDOffSet !=(8,4|
      INTGRL(OffSet)=ABz*HRR(8)+  & 
        HRR(18)+  & 
        INTGRL(OffSet)
      OffSet=(OA+4)*LDA+(OB+3)*LDB+CDOffSet !=(9,4|
      INTGRL(OffSet)=ABz*HRR(9)+  & 
        HRR(19)+  & 
        INTGRL(OffSet)
      OffSet=(OA+5)*LDA+(OB+3)*LDB+CDOffSet !=(10,4|
      INTGRL(OffSet)=ABz*HRR(10)+  & 
        HRR(20)+  & 
        INTGRL(OffSet)
END SUBROUTINE BraHRR62
