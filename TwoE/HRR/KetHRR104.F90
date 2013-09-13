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
   SUBROUTINE KetHRR104(LB,HRR)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB
      REAL(DOUBLE) :: HRR(1:LB,56,10)
      !=|11,5)
      HRR(1:LB,11,5)=CDx*(CDx*HRR(1:LB,11,1)+  &
                        2.D0*HRR(1:LB,21,1))+  &
                        HRR(1:LB,36,1)
      !=|12,5)
      HRR(1:LB,12,5)=CDx*(CDx*HRR(1:LB,12,1)+  &
                        2.D0*HRR(1:LB,22,1))+  &
                        HRR(1:LB,37,1)
      !=|13,5)
      HRR(1:LB,13,5)=CDx*(CDx*HRR(1:LB,13,1)+  &
                        2.D0*HRR(1:LB,23,1))+  &
                        HRR(1:LB,38,1)
      !=|14,5)
      HRR(1:LB,14,5)=CDx*(CDx*HRR(1:LB,14,1)+  &
                        2.D0*HRR(1:LB,24,1))+  &
                        HRR(1:LB,39,1)
      !=|15,5)
      HRR(1:LB,15,5)=CDx*(CDx*HRR(1:LB,15,1)+  &
                        2.D0*HRR(1:LB,26,1))+  &
                        HRR(1:LB,42,1)
      !=|16,5)
      HRR(1:LB,16,5)=CDx*(CDx*HRR(1:LB,16,1)+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        HRR(1:LB,43,1)
      !=|17,5)
      HRR(1:LB,17,5)=CDx*(CDx*HRR(1:LB,17,1)+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,44,1)
      !=|18,5)
      HRR(1:LB,18,5)=CDx*(CDx*HRR(1:LB,18,1)+  &
                        2.D0*HRR(1:LB,30,1))+  &
                        HRR(1:LB,47,1)
      !=|19,5)
      HRR(1:LB,19,5)=CDx*(CDx*HRR(1:LB,19,1)+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,48,1)
      !=|20,5)
      HRR(1:LB,20,5)=CDx*(CDx*HRR(1:LB,20,1)+  &
                        2.D0*HRR(1:LB,33,1))+  &
                        HRR(1:LB,51,1)
      !=|11,6)
      HRR(1:LB,11,6)=CDy*HRR(1:LB,21,1)+  &
                        CDx*(CDy*HRR(1:LB,11,1)+  &
                        HRR(1:LB,22,1))+  &
                        HRR(1:LB,37,1)
      !=|12,6)
      HRR(1:LB,12,6)=CDy*HRR(1:LB,22,1)+  &
                        CDx*(CDy*HRR(1:LB,12,1)+  &
                        HRR(1:LB,23,1))+  &
                        HRR(1:LB,38,1)
      !=|13,6)
      HRR(1:LB,13,6)=CDy*HRR(1:LB,23,1)+  &
                        CDx*(CDy*HRR(1:LB,13,1)+  &
                        HRR(1:LB,24,1))+  &
                        HRR(1:LB,39,1)
      !=|14,6)
      HRR(1:LB,14,6)=CDy*HRR(1:LB,24,1)+  &
                        CDx*(CDy*HRR(1:LB,14,1)+  &
                        HRR(1:LB,25,1))+  &
                        HRR(1:LB,40,1)
      !=|15,6)
      HRR(1:LB,15,6)=CDy*HRR(1:LB,26,1)+  &
                        CDx*(CDy*HRR(1:LB,15,1)+  &
                        HRR(1:LB,27,1))+  &
                        HRR(1:LB,43,1)
      !=|16,6)
      HRR(1:LB,16,6)=CDy*HRR(1:LB,27,1)+  &
                        CDx*(CDy*HRR(1:LB,16,1)+  &
                        HRR(1:LB,28,1))+  &
                        HRR(1:LB,44,1)
      !=|17,6)
      HRR(1:LB,17,6)=CDy*HRR(1:LB,28,1)+  &
                        CDx*(CDy*HRR(1:LB,17,1)+  &
                        HRR(1:LB,29,1))+  &
                        HRR(1:LB,45,1)
      !=|18,6)
      HRR(1:LB,18,6)=CDy*HRR(1:LB,30,1)+  &
                        CDx*(CDy*HRR(1:LB,18,1)+  &
                        HRR(1:LB,31,1))+  &
                        HRR(1:LB,48,1)
      !=|19,6)
      HRR(1:LB,19,6)=CDy*HRR(1:LB,31,1)+  &
                        CDx*(CDy*HRR(1:LB,19,1)+  &
                        HRR(1:LB,32,1))+  &
                        HRR(1:LB,49,1)
      !=|20,6)
      HRR(1:LB,20,6)=CDy*HRR(1:LB,33,1)+  &
                        CDx*(CDy*HRR(1:LB,20,1)+  &
                        HRR(1:LB,34,1))+  &
                        HRR(1:LB,52,1)
      !=|11,7)
      HRR(1:LB,11,7)=CDy*(CDy*HRR(1:LB,11,1)+  &
                        2.D0*HRR(1:LB,22,1))+  &
                        HRR(1:LB,38,1)
      !=|12,7)
      HRR(1:LB,12,7)=CDy*(CDy*HRR(1:LB,12,1)+  &
                        2.D0*HRR(1:LB,23,1))+  &
                        HRR(1:LB,39,1)
      !=|13,7)
      HRR(1:LB,13,7)=CDy*(CDy*HRR(1:LB,13,1)+  &
                        2.D0*HRR(1:LB,24,1))+  &
                        HRR(1:LB,40,1)
      !=|14,7)
      HRR(1:LB,14,7)=CDy*(CDy*HRR(1:LB,14,1)+  &
                        2.D0*HRR(1:LB,25,1))+  &
                        HRR(1:LB,41,1)
      !=|15,7)
      HRR(1:LB,15,7)=CDy*(CDy*HRR(1:LB,15,1)+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        HRR(1:LB,44,1)
      !=|16,7)
      HRR(1:LB,16,7)=CDy*(CDy*HRR(1:LB,16,1)+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,45,1)
      !=|17,7)
      HRR(1:LB,17,7)=CDy*(CDy*HRR(1:LB,17,1)+  &
                        2.D0*HRR(1:LB,29,1))+  &
                        HRR(1:LB,46,1)
      !=|18,7)
      HRR(1:LB,18,7)=CDy*(CDy*HRR(1:LB,18,1)+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,49,1)
      !=|19,7)
      HRR(1:LB,19,7)=CDy*(CDy*HRR(1:LB,19,1)+  &
                        2.D0*HRR(1:LB,32,1))+  &
                        HRR(1:LB,50,1)
      !=|20,7)
      HRR(1:LB,20,7)=CDy*(CDy*HRR(1:LB,20,1)+  &
                        2.D0*HRR(1:LB,34,1))+  &
                        HRR(1:LB,53,1)
      !=|11,8)
      HRR(1:LB,11,8)=CDz*HRR(1:LB,21,1)+  &
                        CDx*(CDz*HRR(1:LB,11,1)+  &
                        HRR(1:LB,26,1))+  &
                        HRR(1:LB,42,1)
      !=|12,8)
      HRR(1:LB,12,8)=CDz*HRR(1:LB,22,1)+  &
                        CDx*(CDz*HRR(1:LB,12,1)+  &
                        HRR(1:LB,27,1))+  &
                        HRR(1:LB,43,1)
      !=|13,8)
      HRR(1:LB,13,8)=CDz*HRR(1:LB,23,1)+  &
                        CDx*(CDz*HRR(1:LB,13,1)+  &
                        HRR(1:LB,28,1))+  &
                        HRR(1:LB,44,1)
      !=|14,8)
      HRR(1:LB,14,8)=CDz*HRR(1:LB,24,1)+  &
                        CDx*(CDz*HRR(1:LB,14,1)+  &
                        HRR(1:LB,29,1))+  &
                        HRR(1:LB,45,1)
      !=|15,8)
      HRR(1:LB,15,8)=CDz*HRR(1:LB,26,1)+  &
                        CDx*(CDz*HRR(1:LB,15,1)+  &
                        HRR(1:LB,30,1))+  &
                        HRR(1:LB,47,1)
      !=|16,8)
      HRR(1:LB,16,8)=CDz*HRR(1:LB,27,1)+  &
                        CDx*(CDz*HRR(1:LB,16,1)+  &
                        HRR(1:LB,31,1))+  &
                        HRR(1:LB,48,1)
      !=|17,8)
      HRR(1:LB,17,8)=CDz*HRR(1:LB,28,1)+  &
                        CDx*(CDz*HRR(1:LB,17,1)+  &
                        HRR(1:LB,32,1))+  &
                        HRR(1:LB,49,1)
      !=|18,8)
      HRR(1:LB,18,8)=CDz*HRR(1:LB,30,1)+  &
                        CDx*(CDz*HRR(1:LB,18,1)+  &
                        HRR(1:LB,33,1))+  &
                        HRR(1:LB,51,1)
      !=|19,8)
      HRR(1:LB,19,8)=CDz*HRR(1:LB,31,1)+  &
                        CDx*(CDz*HRR(1:LB,19,1)+  &
                        HRR(1:LB,34,1))+  &
                        HRR(1:LB,52,1)
      !=|20,8)
      HRR(1:LB,20,8)=CDz*HRR(1:LB,33,1)+  &
                        CDx*(CDz*HRR(1:LB,20,1)+  &
                        HRR(1:LB,35,1))+  &
                        HRR(1:LB,54,1)
      !=|11,9)
      HRR(1:LB,11,9)=CDz*HRR(1:LB,22,1)+  &
                        CDy*(CDz*HRR(1:LB,11,1)+  &
                        HRR(1:LB,26,1))+  &
                        HRR(1:LB,43,1)
      !=|12,9)
      HRR(1:LB,12,9)=CDz*HRR(1:LB,23,1)+  &
                        CDy*(CDz*HRR(1:LB,12,1)+  &
                        HRR(1:LB,27,1))+  &
                        HRR(1:LB,44,1)
      !=|13,9)
      HRR(1:LB,13,9)=CDz*HRR(1:LB,24,1)+  &
                        CDy*(CDz*HRR(1:LB,13,1)+  &
                        HRR(1:LB,28,1))+  &
                        HRR(1:LB,45,1)
      !=|14,9)
      HRR(1:LB,14,9)=CDz*HRR(1:LB,25,1)+  &
                        CDy*(CDz*HRR(1:LB,14,1)+  &
                        HRR(1:LB,29,1))+  &
                        HRR(1:LB,46,1)
      !=|15,9)
      HRR(1:LB,15,9)=CDz*HRR(1:LB,27,1)+  &
                        CDy*(CDz*HRR(1:LB,15,1)+  &
                        HRR(1:LB,30,1))+  &
                        HRR(1:LB,48,1)
      !=|16,9)
      HRR(1:LB,16,9)=CDz*HRR(1:LB,28,1)+  &
                        CDy*(CDz*HRR(1:LB,16,1)+  &
                        HRR(1:LB,31,1))+  &
                        HRR(1:LB,49,1)
      !=|17,9)
      HRR(1:LB,17,9)=CDz*HRR(1:LB,29,1)+  &
                        CDy*(CDz*HRR(1:LB,17,1)+  &
                        HRR(1:LB,32,1))+  &
                        HRR(1:LB,50,1)
      !=|18,9)
      HRR(1:LB,18,9)=CDz*HRR(1:LB,31,1)+  &
                        CDy*(CDz*HRR(1:LB,18,1)+  &
                        HRR(1:LB,33,1))+  &
                        HRR(1:LB,52,1)
      !=|19,9)
      HRR(1:LB,19,9)=CDz*HRR(1:LB,32,1)+  &
                        CDy*(CDz*HRR(1:LB,19,1)+  &
                        HRR(1:LB,34,1))+  &
                        HRR(1:LB,53,1)
      !=|20,9)
      HRR(1:LB,20,9)=CDz*HRR(1:LB,34,1)+  &
                        CDy*(CDz*HRR(1:LB,20,1)+  &
                        HRR(1:LB,35,1))+  &
                        HRR(1:LB,55,1)
      !=|11,10)
      HRR(1:LB,11,10)=CDz*(CDz*HRR(1:LB,11,1)+  &
                        2.D0*HRR(1:LB,26,1))+  &
                        HRR(1:LB,47,1)
      !=|12,10)
      HRR(1:LB,12,10)=CDz*(CDz*HRR(1:LB,12,1)+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        HRR(1:LB,48,1)
      !=|13,10)
      HRR(1:LB,13,10)=CDz*(CDz*HRR(1:LB,13,1)+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,49,1)
      !=|14,10)
      HRR(1:LB,14,10)=CDz*(CDz*HRR(1:LB,14,1)+  &
                        2.D0*HRR(1:LB,29,1))+  &
                        HRR(1:LB,50,1)
      !=|15,10)
      HRR(1:LB,15,10)=CDz*(CDz*HRR(1:LB,15,1)+  &
                        2.D0*HRR(1:LB,30,1))+  &
                        HRR(1:LB,51,1)
      !=|16,10)
      HRR(1:LB,16,10)=CDz*(CDz*HRR(1:LB,16,1)+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,52,1)
      !=|17,10)
      HRR(1:LB,17,10)=CDz*(CDz*HRR(1:LB,17,1)+  &
                        2.D0*HRR(1:LB,32,1))+  &
                        HRR(1:LB,53,1)
      !=|18,10)
      HRR(1:LB,18,10)=CDz*(CDz*HRR(1:LB,18,1)+  &
                        2.D0*HRR(1:LB,33,1))+  &
                        HRR(1:LB,54,1)
      !=|19,10)
      HRR(1:LB,19,10)=CDz*(CDz*HRR(1:LB,19,1)+  &
                        2.D0*HRR(1:LB,34,1))+  &
                        HRR(1:LB,55,1)
      !=|20,10)
      HRR(1:LB,20,10)=CDz*(CDz*HRR(1:LB,20,1)+  &
                        2.D0*HRR(1:LB,35,1))+  &
                        HRR(1:LB,56,1)
      !=|5,5)
      HRR(1:LB,5,5)=CDx*(CDx*HRR(1:LB,5,1)+  &
                        2.D0*HRR(1:LB,11,1))+  &
                        HRR(1:LB,21,1)
      !=|6,5)
      HRR(1:LB,6,5)=CDx*(CDx*HRR(1:LB,6,1)+  &
                        2.D0*HRR(1:LB,12,1))+  &
                        HRR(1:LB,22,1)
      !=|7,5)
      HRR(1:LB,7,5)=CDx*(CDx*HRR(1:LB,7,1)+  &
                        2.D0*HRR(1:LB,13,1))+  &
                        HRR(1:LB,23,1)
      !=|8,5)
      HRR(1:LB,8,5)=CDx*(CDx*HRR(1:LB,8,1)+  &
                        2.D0*HRR(1:LB,15,1))+  &
                        HRR(1:LB,26,1)
      !=|9,5)
      HRR(1:LB,9,5)=CDx*(CDx*HRR(1:LB,9,1)+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        HRR(1:LB,27,1)
      !=|10,5)
      HRR(1:LB,10,5)=CDx*(CDx*HRR(1:LB,10,1)+  &
                        2.D0*HRR(1:LB,18,1))+  &
                        HRR(1:LB,30,1)
      !=|5,6)
      HRR(1:LB,5,6)=CDy*HRR(1:LB,11,1)+  &
                        CDx*(CDy*HRR(1:LB,5,1)+  &
                        HRR(1:LB,12,1))+  &
                        HRR(1:LB,22,1)
      !=|6,6)
      HRR(1:LB,6,6)=CDy*HRR(1:LB,12,1)+  &
                        CDx*(CDy*HRR(1:LB,6,1)+  &
                        HRR(1:LB,13,1))+  &
                        HRR(1:LB,23,1)
      !=|7,6)
      HRR(1:LB,7,6)=CDy*HRR(1:LB,13,1)+  &
                        CDx*(CDy*HRR(1:LB,7,1)+  &
                        HRR(1:LB,14,1))+  &
                        HRR(1:LB,24,1)
      !=|8,6)
      HRR(1:LB,8,6)=CDy*HRR(1:LB,15,1)+  &
                        CDx*(CDy*HRR(1:LB,8,1)+  &
                        HRR(1:LB,16,1))+  &
                        HRR(1:LB,27,1)
      !=|9,6)
      HRR(1:LB,9,6)=CDy*HRR(1:LB,16,1)+  &
                        CDx*(CDy*HRR(1:LB,9,1)+  &
                        HRR(1:LB,17,1))+  &
                        HRR(1:LB,28,1)
      !=|10,6)
      HRR(1:LB,10,6)=CDy*HRR(1:LB,18,1)+  &
                        CDx*(CDy*HRR(1:LB,10,1)+  &
                        HRR(1:LB,19,1))+  &
                        HRR(1:LB,31,1)
      !=|5,7)
      HRR(1:LB,5,7)=CDy*(CDy*HRR(1:LB,5,1)+  &
                        2.D0*HRR(1:LB,12,1))+  &
                        HRR(1:LB,23,1)
      !=|6,7)
      HRR(1:LB,6,7)=CDy*(CDy*HRR(1:LB,6,1)+  &
                        2.D0*HRR(1:LB,13,1))+  &
                        HRR(1:LB,24,1)
      !=|7,7)
      HRR(1:LB,7,7)=CDy*(CDy*HRR(1:LB,7,1)+  &
                        2.D0*HRR(1:LB,14,1))+  &
                        HRR(1:LB,25,1)
      !=|8,7)
      HRR(1:LB,8,7)=CDy*(CDy*HRR(1:LB,8,1)+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        HRR(1:LB,28,1)
      !=|9,7)
      HRR(1:LB,9,7)=CDy*(CDy*HRR(1:LB,9,1)+  &
                        2.D0*HRR(1:LB,17,1))+  &
                        HRR(1:LB,29,1)
      !=|10,7)
      HRR(1:LB,10,7)=CDy*(CDy*HRR(1:LB,10,1)+  &
                        2.D0*HRR(1:LB,19,1))+  &
                        HRR(1:LB,32,1)
      !=|5,8)
      HRR(1:LB,5,8)=CDz*HRR(1:LB,11,1)+  &
                        CDx*(CDz*HRR(1:LB,5,1)+  &
                        HRR(1:LB,15,1))+  &
                        HRR(1:LB,26,1)
      !=|6,8)
      HRR(1:LB,6,8)=CDz*HRR(1:LB,12,1)+  &
                        CDx*(CDz*HRR(1:LB,6,1)+  &
                        HRR(1:LB,16,1))+  &
                        HRR(1:LB,27,1)
      !=|7,8)
      HRR(1:LB,7,8)=CDz*HRR(1:LB,13,1)+  &
                        CDx*(CDz*HRR(1:LB,7,1)+  &
                        HRR(1:LB,17,1))+  &
                        HRR(1:LB,28,1)
      !=|8,8)
      HRR(1:LB,8,8)=CDz*HRR(1:LB,15,1)+  &
                        CDx*(CDz*HRR(1:LB,8,1)+  &
                        HRR(1:LB,18,1))+  &
                        HRR(1:LB,30,1)
      !=|9,8)
      HRR(1:LB,9,8)=CDz*HRR(1:LB,16,1)+  &
                        CDx*(CDz*HRR(1:LB,9,1)+  &
                        HRR(1:LB,19,1))+  &
                        HRR(1:LB,31,1)
      !=|10,8)
      HRR(1:LB,10,8)=CDz*HRR(1:LB,18,1)+  &
                        CDx*(CDz*HRR(1:LB,10,1)+  &
                        HRR(1:LB,20,1))+  &
                        HRR(1:LB,33,1)
      !=|5,9)
      HRR(1:LB,5,9)=CDz*HRR(1:LB,12,1)+  &
                        CDy*(CDz*HRR(1:LB,5,1)+  &
                        HRR(1:LB,15,1))+  &
                        HRR(1:LB,27,1)
      !=|6,9)
      HRR(1:LB,6,9)=CDz*HRR(1:LB,13,1)+  &
                        CDy*(CDz*HRR(1:LB,6,1)+  &
                        HRR(1:LB,16,1))+  &
                        HRR(1:LB,28,1)
      !=|7,9)
      HRR(1:LB,7,9)=CDz*HRR(1:LB,14,1)+  &
                        CDy*(CDz*HRR(1:LB,7,1)+  &
                        HRR(1:LB,17,1))+  &
                        HRR(1:LB,29,1)
      !=|8,9)
      HRR(1:LB,8,9)=CDz*HRR(1:LB,16,1)+  &
                        CDy*(CDz*HRR(1:LB,8,1)+  &
                        HRR(1:LB,18,1))+  &
                        HRR(1:LB,31,1)
      !=|9,9)
      HRR(1:LB,9,9)=CDz*HRR(1:LB,17,1)+  &
                        CDy*(CDz*HRR(1:LB,9,1)+  &
                        HRR(1:LB,19,1))+  &
                        HRR(1:LB,32,1)
      !=|10,9)
      HRR(1:LB,10,9)=CDz*HRR(1:LB,19,1)+  &
                        CDy*(CDz*HRR(1:LB,10,1)+  &
                        HRR(1:LB,20,1))+  &
                        HRR(1:LB,34,1)
      !=|5,10)
      HRR(1:LB,5,10)=CDz*(CDz*HRR(1:LB,5,1)+  &
                        2.D0*HRR(1:LB,15,1))+  &
                        HRR(1:LB,30,1)
      !=|6,10)
      HRR(1:LB,6,10)=CDz*(CDz*HRR(1:LB,6,1)+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        HRR(1:LB,31,1)
      !=|7,10)
      HRR(1:LB,7,10)=CDz*(CDz*HRR(1:LB,7,1)+  &
                        2.D0*HRR(1:LB,17,1))+  &
                        HRR(1:LB,32,1)
      !=|8,10)
      HRR(1:LB,8,10)=CDz*(CDz*HRR(1:LB,8,1)+  &
                        2.D0*HRR(1:LB,18,1))+  &
                        HRR(1:LB,33,1)
      !=|9,10)
      HRR(1:LB,9,10)=CDz*(CDz*HRR(1:LB,9,1)+  &
                        2.D0*HRR(1:LB,19,1))+  &
                        HRR(1:LB,34,1)
      !=|10,10)
      HRR(1:LB,10,10)=CDz*(CDz*HRR(1:LB,10,1)+  &
                        2.D0*HRR(1:LB,20,1))+  &
                        HRR(1:LB,35,1)
      !=|2,5)
      HRR(1:LB,2,5)=CDx*(CDx*HRR(1:LB,2,1)+  &
                        2.D0*HRR(1:LB,5,1))+  &
                        HRR(1:LB,11,1)
      !=|3,5)
      HRR(1:LB,3,5)=CDx*(CDx*HRR(1:LB,3,1)+  &
                        2.D0*HRR(1:LB,6,1))+  &
                        HRR(1:LB,12,1)
      !=|4,5)
      HRR(1:LB,4,5)=CDx*(CDx*HRR(1:LB,4,1)+  &
                        2.D0*HRR(1:LB,8,1))+  &
                        HRR(1:LB,15,1)
      !=|2,6)
      HRR(1:LB,2,6)=CDy*HRR(1:LB,5,1)+  &
                        CDx*(CDy*HRR(1:LB,2,1)+  &
                        HRR(1:LB,6,1))+  &
                        HRR(1:LB,12,1)
      !=|3,6)
      HRR(1:LB,3,6)=CDy*HRR(1:LB,6,1)+  &
                        CDx*(CDy*HRR(1:LB,3,1)+  &
                        HRR(1:LB,7,1))+  &
                        HRR(1:LB,13,1)
      !=|4,6)
      HRR(1:LB,4,6)=CDy*HRR(1:LB,8,1)+  &
                        CDx*(CDy*HRR(1:LB,4,1)+  &
                        HRR(1:LB,9,1))+  &
                        HRR(1:LB,16,1)
      !=|2,7)
      HRR(1:LB,2,7)=CDy*(CDy*HRR(1:LB,2,1)+  &
                        2.D0*HRR(1:LB,6,1))+  &
                        HRR(1:LB,13,1)
      !=|3,7)
      HRR(1:LB,3,7)=CDy*(CDy*HRR(1:LB,3,1)+  &
                        2.D0*HRR(1:LB,7,1))+  &
                        HRR(1:LB,14,1)
      !=|4,7)
      HRR(1:LB,4,7)=CDy*(CDy*HRR(1:LB,4,1)+  &
                        2.D0*HRR(1:LB,9,1))+  &
                        HRR(1:LB,17,1)
      !=|2,8)
      HRR(1:LB,2,8)=CDz*HRR(1:LB,5,1)+  &
                        CDx*(CDz*HRR(1:LB,2,1)+  &
                        HRR(1:LB,8,1))+  &
                        HRR(1:LB,15,1)
      !=|3,8)
      HRR(1:LB,3,8)=CDz*HRR(1:LB,6,1)+  &
                        CDx*(CDz*HRR(1:LB,3,1)+  &
                        HRR(1:LB,9,1))+  &
                        HRR(1:LB,16,1)
      !=|4,8)
      HRR(1:LB,4,8)=CDz*HRR(1:LB,8,1)+  &
                        CDx*(CDz*HRR(1:LB,4,1)+  &
                        HRR(1:LB,10,1))+  &
                        HRR(1:LB,18,1)
      !=|2,9)
      HRR(1:LB,2,9)=CDz*HRR(1:LB,6,1)+  &
                        CDy*(CDz*HRR(1:LB,2,1)+  &
                        HRR(1:LB,8,1))+  &
                        HRR(1:LB,16,1)
      !=|3,9)
      HRR(1:LB,3,9)=CDz*HRR(1:LB,7,1)+  &
                        CDy*(CDz*HRR(1:LB,3,1)+  &
                        HRR(1:LB,9,1))+  &
                        HRR(1:LB,17,1)
      !=|4,9)
      HRR(1:LB,4,9)=CDz*HRR(1:LB,9,1)+  &
                        CDy*(CDz*HRR(1:LB,4,1)+  &
                        HRR(1:LB,10,1))+  &
                        HRR(1:LB,19,1)
      !=|2,10)
      HRR(1:LB,2,10)=CDz*(CDz*HRR(1:LB,2,1)+  &
                        2.D0*HRR(1:LB,8,1))+  &
                        HRR(1:LB,18,1)
      !=|3,10)
      HRR(1:LB,3,10)=CDz*(CDz*HRR(1:LB,3,1)+  &
                        2.D0*HRR(1:LB,9,1))+  &
                        HRR(1:LB,19,1)
      !=|4,10)
      HRR(1:LB,4,10)=CDz*(CDz*HRR(1:LB,4,1)+  &
                        2.D0*HRR(1:LB,10,1))+  &
                        HRR(1:LB,20,1)
      !=|1,5)
      HRR(1:LB,1,5)=CDx*(CDx*HRR(1:LB,1,1)+  &
                        2.D0*HRR(1:LB,2,1))+  &
                        HRR(1:LB,5,1)
      !=|1,6)
      HRR(1:LB,1,6)=CDy*HRR(1:LB,2,1)+  &
                        CDx*(CDy*HRR(1:LB,1,1)+  &
                        HRR(1:LB,3,1))+  &
                        HRR(1:LB,6,1)
      !=|1,7)
      HRR(1:LB,1,7)=CDy*(CDy*HRR(1:LB,1,1)+  &
                        2.D0*HRR(1:LB,3,1))+  &
                        HRR(1:LB,7,1)
      !=|1,8)
      HRR(1:LB,1,8)=CDz*HRR(1:LB,2,1)+  &
                        CDx*(CDz*HRR(1:LB,1,1)+  &
                        HRR(1:LB,4,1))+  &
                        HRR(1:LB,8,1)
      !=|1,9)
      HRR(1:LB,1,9)=CDz*HRR(1:LB,3,1)+  &
                        CDy*(CDz*HRR(1:LB,1,1)+  &
                        HRR(1:LB,4,1))+  &
                        HRR(1:LB,9,1)
      !=|1,10)
      HRR(1:LB,1,10)=CDz*(CDz*HRR(1:LB,1,1)+  &
                        2.D0*HRR(1:LB,4,1))+  &
                        HRR(1:LB,10,1)
      !=|11,2)
      HRR(1:LB,11,2)=CDx*HRR(1:LB,11,1)+  &
                        HRR(1:LB,21,1)
      !=|12,2)
      HRR(1:LB,12,2)=CDx*HRR(1:LB,12,1)+  &
                        HRR(1:LB,22,1)
      !=|13,2)
      HRR(1:LB,13,2)=CDx*HRR(1:LB,13,1)+  &
                        HRR(1:LB,23,1)
      !=|14,2)
      HRR(1:LB,14,2)=CDx*HRR(1:LB,14,1)+  &
                        HRR(1:LB,24,1)
      !=|15,2)
      HRR(1:LB,15,2)=CDx*HRR(1:LB,15,1)+  &
                        HRR(1:LB,26,1)
      !=|16,2)
      HRR(1:LB,16,2)=CDx*HRR(1:LB,16,1)+  &
                        HRR(1:LB,27,1)
      !=|17,2)
      HRR(1:LB,17,2)=CDx*HRR(1:LB,17,1)+  &
                        HRR(1:LB,28,1)
      !=|18,2)
      HRR(1:LB,18,2)=CDx*HRR(1:LB,18,1)+  &
                        HRR(1:LB,30,1)
      !=|19,2)
      HRR(1:LB,19,2)=CDx*HRR(1:LB,19,1)+  &
                        HRR(1:LB,31,1)
      !=|20,2)
      HRR(1:LB,20,2)=CDx*HRR(1:LB,20,1)+  &
                        HRR(1:LB,33,1)
      !=|11,3)
      HRR(1:LB,11,3)=CDy*HRR(1:LB,11,1)+  &
                        HRR(1:LB,22,1)
      !=|12,3)
      HRR(1:LB,12,3)=CDy*HRR(1:LB,12,1)+  &
                        HRR(1:LB,23,1)
      !=|13,3)
      HRR(1:LB,13,3)=CDy*HRR(1:LB,13,1)+  &
                        HRR(1:LB,24,1)
      !=|14,3)
      HRR(1:LB,14,3)=CDy*HRR(1:LB,14,1)+  &
                        HRR(1:LB,25,1)
      !=|15,3)
      HRR(1:LB,15,3)=CDy*HRR(1:LB,15,1)+  &
                        HRR(1:LB,27,1)
      !=|16,3)
      HRR(1:LB,16,3)=CDy*HRR(1:LB,16,1)+  &
                        HRR(1:LB,28,1)
      !=|17,3)
      HRR(1:LB,17,3)=CDy*HRR(1:LB,17,1)+  &
                        HRR(1:LB,29,1)
      !=|18,3)
      HRR(1:LB,18,3)=CDy*HRR(1:LB,18,1)+  &
                        HRR(1:LB,31,1)
      !=|19,3)
      HRR(1:LB,19,3)=CDy*HRR(1:LB,19,1)+  &
                        HRR(1:LB,32,1)
      !=|20,3)
      HRR(1:LB,20,3)=CDy*HRR(1:LB,20,1)+  &
                        HRR(1:LB,34,1)
      !=|11,4)
      HRR(1:LB,11,4)=CDz*HRR(1:LB,11,1)+  &
                        HRR(1:LB,26,1)
      !=|12,4)
      HRR(1:LB,12,4)=CDz*HRR(1:LB,12,1)+  &
                        HRR(1:LB,27,1)
      !=|13,4)
      HRR(1:LB,13,4)=CDz*HRR(1:LB,13,1)+  &
                        HRR(1:LB,28,1)
      !=|14,4)
      HRR(1:LB,14,4)=CDz*HRR(1:LB,14,1)+  &
                        HRR(1:LB,29,1)
      !=|15,4)
      HRR(1:LB,15,4)=CDz*HRR(1:LB,15,1)+  &
                        HRR(1:LB,30,1)
      !=|16,4)
      HRR(1:LB,16,4)=CDz*HRR(1:LB,16,1)+  &
                        HRR(1:LB,31,1)
      !=|17,4)
      HRR(1:LB,17,4)=CDz*HRR(1:LB,17,1)+  &
                        HRR(1:LB,32,1)
      !=|18,4)
      HRR(1:LB,18,4)=CDz*HRR(1:LB,18,1)+  &
                        HRR(1:LB,33,1)
      !=|19,4)
      HRR(1:LB,19,4)=CDz*HRR(1:LB,19,1)+  &
                        HRR(1:LB,34,1)
      !=|20,4)
      HRR(1:LB,20,4)=CDz*HRR(1:LB,20,1)+  &
                        HRR(1:LB,35,1)
      !=|5,2)
      HRR(1:LB,5,2)=CDx*HRR(1:LB,5,1)+  &
                        HRR(1:LB,11,1)
      !=|6,2)
      HRR(1:LB,6,2)=CDx*HRR(1:LB,6,1)+  &
                        HRR(1:LB,12,1)
      !=|7,2)
      HRR(1:LB,7,2)=CDx*HRR(1:LB,7,1)+  &
                        HRR(1:LB,13,1)
      !=|8,2)
      HRR(1:LB,8,2)=CDx*HRR(1:LB,8,1)+  &
                        HRR(1:LB,15,1)
      !=|9,2)
      HRR(1:LB,9,2)=CDx*HRR(1:LB,9,1)+  &
                        HRR(1:LB,16,1)
      !=|10,2)
      HRR(1:LB,10,2)=CDx*HRR(1:LB,10,1)+  &
                        HRR(1:LB,18,1)
      !=|5,3)
      HRR(1:LB,5,3)=CDy*HRR(1:LB,5,1)+  &
                        HRR(1:LB,12,1)
      !=|6,3)
      HRR(1:LB,6,3)=CDy*HRR(1:LB,6,1)+  &
                        HRR(1:LB,13,1)
      !=|7,3)
      HRR(1:LB,7,3)=CDy*HRR(1:LB,7,1)+  &
                        HRR(1:LB,14,1)
      !=|8,3)
      HRR(1:LB,8,3)=CDy*HRR(1:LB,8,1)+  &
                        HRR(1:LB,16,1)
      !=|9,3)
      HRR(1:LB,9,3)=CDy*HRR(1:LB,9,1)+  &
                        HRR(1:LB,17,1)
      !=|10,3)
      HRR(1:LB,10,3)=CDy*HRR(1:LB,10,1)+  &
                        HRR(1:LB,19,1)
      !=|5,4)
      HRR(1:LB,5,4)=CDz*HRR(1:LB,5,1)+  &
                        HRR(1:LB,15,1)
      !=|6,4)
      HRR(1:LB,6,4)=CDz*HRR(1:LB,6,1)+  &
                        HRR(1:LB,16,1)
      !=|7,4)
      HRR(1:LB,7,4)=CDz*HRR(1:LB,7,1)+  &
                        HRR(1:LB,17,1)
      !=|8,4)
      HRR(1:LB,8,4)=CDz*HRR(1:LB,8,1)+  &
                        HRR(1:LB,18,1)
      !=|9,4)
      HRR(1:LB,9,4)=CDz*HRR(1:LB,9,1)+  &
                        HRR(1:LB,19,1)
      !=|10,4)
      HRR(1:LB,10,4)=CDz*HRR(1:LB,10,1)+  &
                        HRR(1:LB,20,1)
      !=|2,2)
      HRR(1:LB,2,2)=CDx*HRR(1:LB,2,1)+  &
                        HRR(1:LB,5,1)
      !=|3,2)
      HRR(1:LB,3,2)=CDx*HRR(1:LB,3,1)+  &
                        HRR(1:LB,6,1)
      !=|4,2)
      HRR(1:LB,4,2)=CDx*HRR(1:LB,4,1)+  &
                        HRR(1:LB,8,1)
      !=|2,3)
      HRR(1:LB,2,3)=CDy*HRR(1:LB,2,1)+  &
                        HRR(1:LB,6,1)
      !=|3,3)
      HRR(1:LB,3,3)=CDy*HRR(1:LB,3,1)+  &
                        HRR(1:LB,7,1)
      !=|4,3)
      HRR(1:LB,4,3)=CDy*HRR(1:LB,4,1)+  &
                        HRR(1:LB,9,1)
      !=|2,4)
      HRR(1:LB,2,4)=CDz*HRR(1:LB,2,1)+  &
                        HRR(1:LB,8,1)
      !=|3,4)
      HRR(1:LB,3,4)=CDz*HRR(1:LB,3,1)+  &
                        HRR(1:LB,9,1)
      !=|4,4)
      HRR(1:LB,4,4)=CDz*HRR(1:LB,4,1)+  &
                        HRR(1:LB,10,1)
      !=|1,2)
      HRR(1:LB,1,2)=CDx*HRR(1:LB,1,1)+  &
                        HRR(1:LB,2,1)
      !=|1,3)
      HRR(1:LB,1,3)=CDy*HRR(1:LB,1,1)+  &
                        HRR(1:LB,3,1)
      !=|1,4)
      HRR(1:LB,1,4)=CDz*HRR(1:LB,1,1)+  &
                        HRR(1:LB,4,1)
END SUBROUTINE KetHRR104
