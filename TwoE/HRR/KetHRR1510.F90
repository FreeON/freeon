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
   SUBROUTINE KetHRR1510(LB,HRR)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB
      REAL(DOUBLE) :: HRR(1:LB,120,20)
      !=|21,11)
      HRR(1:LB,21,11)=CDx*(CDx*(CDx*HRR(1:LB,21,1)+  &
                        3.D0*HRR(1:LB,36,1))+  &
                        3.D0*HRR(1:LB,57,1))+  &
                        HRR(1:LB,85,1)
      !=|22,11)
      HRR(1:LB,22,11)=CDx*(CDx*(CDx*HRR(1:LB,22,1)+  &
                        3.D0*HRR(1:LB,37,1))+  &
                        3.D0*HRR(1:LB,58,1))+  &
                        HRR(1:LB,86,1)
      !=|23,11)
      HRR(1:LB,23,11)=CDx*(CDx*(CDx*HRR(1:LB,23,1)+  &
                        3.D0*HRR(1:LB,38,1))+  &
                        3.D0*HRR(1:LB,59,1))+  &
                        HRR(1:LB,87,1)
      !=|24,11)
      HRR(1:LB,24,11)=CDx*(CDx*(CDx*HRR(1:LB,24,1)+  &
                        3.D0*HRR(1:LB,39,1))+  &
                        3.D0*HRR(1:LB,60,1))+  &
                        HRR(1:LB,88,1)
      !=|25,11)
      HRR(1:LB,25,11)=CDx*(CDx*(CDx*HRR(1:LB,25,1)+  &
                        3.D0*HRR(1:LB,40,1))+  &
                        3.D0*HRR(1:LB,61,1))+  &
                        HRR(1:LB,89,1)
      !=|26,11)
      HRR(1:LB,26,11)=CDx*(CDx*(CDx*HRR(1:LB,26,1)+  &
                        3.D0*HRR(1:LB,42,1))+  &
                        3.D0*HRR(1:LB,64,1))+  &
                        HRR(1:LB,93,1)
      !=|27,11)
      HRR(1:LB,27,11)=CDx*(CDx*(CDx*HRR(1:LB,27,1)+  &
                        3.D0*HRR(1:LB,43,1))+  &
                        3.D0*HRR(1:LB,65,1))+  &
                        HRR(1:LB,94,1)
      !=|28,11)
      HRR(1:LB,28,11)=CDx*(CDx*(CDx*HRR(1:LB,28,1)+  &
                        3.D0*HRR(1:LB,44,1))+  &
                        3.D0*HRR(1:LB,66,1))+  &
                        HRR(1:LB,95,1)
      !=|29,11)
      HRR(1:LB,29,11)=CDx*(CDx*(CDx*HRR(1:LB,29,1)+  &
                        3.D0*HRR(1:LB,45,1))+  &
                        3.D0*HRR(1:LB,67,1))+  &
                        HRR(1:LB,96,1)
      !=|30,11)
      HRR(1:LB,30,11)=CDx*(CDx*(CDx*HRR(1:LB,30,1)+  &
                        3.D0*HRR(1:LB,47,1))+  &
                        3.D0*HRR(1:LB,70,1))+  &
                        HRR(1:LB,100,1)
      !=|31,11)
      HRR(1:LB,31,11)=CDx*(CDx*(CDx*HRR(1:LB,31,1)+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        3.D0*HRR(1:LB,71,1))+  &
                        HRR(1:LB,101,1)
      !=|32,11)
      HRR(1:LB,32,11)=CDx*(CDx*(CDx*HRR(1:LB,32,1)+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        3.D0*HRR(1:LB,72,1))+  &
                        HRR(1:LB,102,1)
      !=|33,11)
      HRR(1:LB,33,11)=CDx*(CDx*(CDx*HRR(1:LB,33,1)+  &
                        3.D0*HRR(1:LB,51,1))+  &
                        3.D0*HRR(1:LB,75,1))+  &
                        HRR(1:LB,106,1)
      !=|34,11)
      HRR(1:LB,34,11)=CDx*(CDx*(CDx*HRR(1:LB,34,1)+  &
                        3.D0*HRR(1:LB,52,1))+  &
                        3.D0*HRR(1:LB,76,1))+  &
                        HRR(1:LB,107,1)
      !=|35,11)
      HRR(1:LB,35,11)=CDx*(CDx*(CDx*HRR(1:LB,35,1)+  &
                        3.D0*HRR(1:LB,54,1))+  &
                        3.D0*HRR(1:LB,79,1))+  &
                        HRR(1:LB,111,1)
      !=|21,12)
      HRR(1:LB,21,12)=CDy*HRR(1:LB,57,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,36,1)+  &
                        CDx*(CDy*HRR(1:LB,21,1)+  &
                        HRR(1:LB,37,1))+  &
                        2.D0*HRR(1:LB,58,1))+  &
                        HRR(1:LB,86,1)
      !=|22,12)
      HRR(1:LB,22,12)=CDy*HRR(1:LB,58,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,37,1)+  &
                        CDx*(CDy*HRR(1:LB,22,1)+  &
                        HRR(1:LB,38,1))+  &
                        2.D0*HRR(1:LB,59,1))+  &
                        HRR(1:LB,87,1)
      !=|23,12)
      HRR(1:LB,23,12)=CDy*HRR(1:LB,59,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,38,1)+  &
                        CDx*(CDy*HRR(1:LB,23,1)+  &
                        HRR(1:LB,39,1))+  &
                        2.D0*HRR(1:LB,60,1))+  &
                        HRR(1:LB,88,1)
      !=|24,12)
      HRR(1:LB,24,12)=CDy*HRR(1:LB,60,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,39,1)+  &
                        CDx*(CDy*HRR(1:LB,24,1)+  &
                        HRR(1:LB,40,1))+  &
                        2.D0*HRR(1:LB,61,1))+  &
                        HRR(1:LB,89,1)
      !=|25,12)
      HRR(1:LB,25,12)=CDy*HRR(1:LB,61,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,40,1)+  &
                        CDx*(CDy*HRR(1:LB,25,1)+  &
                        HRR(1:LB,41,1))+  &
                        2.D0*HRR(1:LB,62,1))+  &
                        HRR(1:LB,90,1)
      !=|26,12)
      HRR(1:LB,26,12)=CDy*HRR(1:LB,64,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,42,1)+  &
                        CDx*(CDy*HRR(1:LB,26,1)+  &
                        HRR(1:LB,43,1))+  &
                        2.D0*HRR(1:LB,65,1))+  &
                        HRR(1:LB,94,1)
      !=|27,12)
      HRR(1:LB,27,12)=CDy*HRR(1:LB,65,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,43,1)+  &
                        CDx*(CDy*HRR(1:LB,27,1)+  &
                        HRR(1:LB,44,1))+  &
                        2.D0*HRR(1:LB,66,1))+  &
                        HRR(1:LB,95,1)
      !=|28,12)
      HRR(1:LB,28,12)=CDy*HRR(1:LB,66,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,44,1)+  &
                        CDx*(CDy*HRR(1:LB,28,1)+  &
                        HRR(1:LB,45,1))+  &
                        2.D0*HRR(1:LB,67,1))+  &
                        HRR(1:LB,96,1)
      !=|29,12)
      HRR(1:LB,29,12)=CDy*HRR(1:LB,67,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,45,1)+  &
                        CDx*(CDy*HRR(1:LB,29,1)+  &
                        HRR(1:LB,46,1))+  &
                        2.D0*HRR(1:LB,68,1))+  &
                        HRR(1:LB,97,1)
      !=|30,12)
      HRR(1:LB,30,12)=CDy*HRR(1:LB,70,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,47,1)+  &
                        CDx*(CDy*HRR(1:LB,30,1)+  &
                        HRR(1:LB,48,1))+  &
                        2.D0*HRR(1:LB,71,1))+  &
                        HRR(1:LB,101,1)
      !=|31,12)
      HRR(1:LB,31,12)=CDy*HRR(1:LB,71,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,48,1)+  &
                        CDx*(CDy*HRR(1:LB,31,1)+  &
                        HRR(1:LB,49,1))+  &
                        2.D0*HRR(1:LB,72,1))+  &
                        HRR(1:LB,102,1)
      !=|32,12)
      HRR(1:LB,32,12)=CDy*HRR(1:LB,72,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,49,1)+  &
                        CDx*(CDy*HRR(1:LB,32,1)+  &
                        HRR(1:LB,50,1))+  &
                        2.D0*HRR(1:LB,73,1))+  &
                        HRR(1:LB,103,1)
      !=|33,12)
      HRR(1:LB,33,12)=CDy*HRR(1:LB,75,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,51,1)+  &
                        CDx*(CDy*HRR(1:LB,33,1)+  &
                        HRR(1:LB,52,1))+  &
                        2.D0*HRR(1:LB,76,1))+  &
                        HRR(1:LB,107,1)
      !=|34,12)
      HRR(1:LB,34,12)=CDy*HRR(1:LB,76,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,52,1)+  &
                        CDx*(CDy*HRR(1:LB,34,1)+  &
                        HRR(1:LB,53,1))+  &
                        2.D0*HRR(1:LB,77,1))+  &
                        HRR(1:LB,108,1)
      !=|35,12)
      HRR(1:LB,35,12)=CDy*HRR(1:LB,79,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,54,1)+  &
                        CDx*(CDy*HRR(1:LB,35,1)+  &
                        HRR(1:LB,55,1))+  &
                        2.D0*HRR(1:LB,80,1))+  &
                        HRR(1:LB,112,1)
      !=|21,13)
      HRR(1:LB,21,13)=CDy*(CDy*HRR(1:LB,36,1)+  &
                        2.D0*HRR(1:LB,58,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,21,1)+  &
                        2.D0*HRR(1:LB,37,1))+  &
                        HRR(1:LB,59,1))+  &
                        HRR(1:LB,87,1)
      !=|22,13)
      HRR(1:LB,22,13)=CDy*(CDy*HRR(1:LB,37,1)+  &
                        2.D0*HRR(1:LB,59,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,22,1)+  &
                        2.D0*HRR(1:LB,38,1))+  &
                        HRR(1:LB,60,1))+  &
                        HRR(1:LB,88,1)
      !=|23,13)
      HRR(1:LB,23,13)=CDy*(CDy*HRR(1:LB,38,1)+  &
                        2.D0*HRR(1:LB,60,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,39,1))+  &
                        HRR(1:LB,61,1))+  &
                        HRR(1:LB,89,1)
      !=|24,13)
      HRR(1:LB,24,13)=CDy*(CDy*HRR(1:LB,39,1)+  &
                        2.D0*HRR(1:LB,61,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,24,1)+  &
                        2.D0*HRR(1:LB,40,1))+  &
                        HRR(1:LB,62,1))+  &
                        HRR(1:LB,90,1)
      !=|25,13)
      HRR(1:LB,25,13)=CDy*(CDy*HRR(1:LB,40,1)+  &
                        2.D0*HRR(1:LB,62,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,25,1)+  &
                        2.D0*HRR(1:LB,41,1))+  &
                        HRR(1:LB,63,1))+  &
                        HRR(1:LB,91,1)
      !=|26,13)
      HRR(1:LB,26,13)=CDy*(CDy*HRR(1:LB,42,1)+  &
                        2.D0*HRR(1:LB,65,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,26,1)+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,66,1))+  &
                        HRR(1:LB,95,1)
      !=|27,13)
      HRR(1:LB,27,13)=CDy*(CDy*HRR(1:LB,43,1)+  &
                        2.D0*HRR(1:LB,66,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,27,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,67,1))+  &
                        HRR(1:LB,96,1)
      !=|28,13)
      HRR(1:LB,28,13)=CDy*(CDy*HRR(1:LB,44,1)+  &
                        2.D0*HRR(1:LB,67,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,28,1)+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,68,1))+  &
                        HRR(1:LB,97,1)
      !=|29,13)
      HRR(1:LB,29,13)=CDy*(CDy*HRR(1:LB,45,1)+  &
                        2.D0*HRR(1:LB,68,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,29,1)+  &
                        2.D0*HRR(1:LB,46,1))+  &
                        HRR(1:LB,69,1))+  &
                        HRR(1:LB,98,1)
      !=|30,13)
      HRR(1:LB,30,13)=CDy*(CDy*HRR(1:LB,47,1)+  &
                        2.D0*HRR(1:LB,71,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,30,1)+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,72,1))+  &
                        HRR(1:LB,102,1)
      !=|31,13)
      HRR(1:LB,31,13)=CDy*(CDy*HRR(1:LB,48,1)+  &
                        2.D0*HRR(1:LB,72,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,31,1)+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,73,1))+  &
                        HRR(1:LB,103,1)
      !=|32,13)
      HRR(1:LB,32,13)=CDy*(CDy*HRR(1:LB,49,1)+  &
                        2.D0*HRR(1:LB,73,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,32,1)+  &
                        2.D0*HRR(1:LB,50,1))+  &
                        HRR(1:LB,74,1))+  &
                        HRR(1:LB,104,1)
      !=|33,13)
      HRR(1:LB,33,13)=CDy*(CDy*HRR(1:LB,51,1)+  &
                        2.D0*HRR(1:LB,76,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,33,1)+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,77,1))+  &
                        HRR(1:LB,108,1)
      !=|34,13)
      HRR(1:LB,34,13)=CDy*(CDy*HRR(1:LB,52,1)+  &
                        2.D0*HRR(1:LB,77,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,34,1)+  &
                        2.D0*HRR(1:LB,53,1))+  &
                        HRR(1:LB,78,1))+  &
                        HRR(1:LB,109,1)
      !=|35,13)
      HRR(1:LB,35,13)=CDy*(CDy*HRR(1:LB,54,1)+  &
                        2.D0*HRR(1:LB,80,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,35,1)+  &
                        2.D0*HRR(1:LB,55,1))+  &
                        HRR(1:LB,81,1))+  &
                        HRR(1:LB,113,1)
      !=|21,14)
      HRR(1:LB,21,14)=CDy*(CDy*(CDy*HRR(1:LB,21,1)+  &
                        3.D0*HRR(1:LB,37,1))+  &
                        3.D0*HRR(1:LB,59,1))+  &
                        HRR(1:LB,88,1)
      !=|22,14)
      HRR(1:LB,22,14)=CDy*(CDy*(CDy*HRR(1:LB,22,1)+  &
                        3.D0*HRR(1:LB,38,1))+  &
                        3.D0*HRR(1:LB,60,1))+  &
                        HRR(1:LB,89,1)
      !=|23,14)
      HRR(1:LB,23,14)=CDy*(CDy*(CDy*HRR(1:LB,23,1)+  &
                        3.D0*HRR(1:LB,39,1))+  &
                        3.D0*HRR(1:LB,61,1))+  &
                        HRR(1:LB,90,1)
      !=|24,14)
      HRR(1:LB,24,14)=CDy*(CDy*(CDy*HRR(1:LB,24,1)+  &
                        3.D0*HRR(1:LB,40,1))+  &
                        3.D0*HRR(1:LB,62,1))+  &
                        HRR(1:LB,91,1)
      !=|25,14)
      HRR(1:LB,25,14)=CDy*(CDy*(CDy*HRR(1:LB,25,1)+  &
                        3.D0*HRR(1:LB,41,1))+  &
                        3.D0*HRR(1:LB,63,1))+  &
                        HRR(1:LB,92,1)
      !=|26,14)
      HRR(1:LB,26,14)=CDy*(CDy*(CDy*HRR(1:LB,26,1)+  &
                        3.D0*HRR(1:LB,43,1))+  &
                        3.D0*HRR(1:LB,66,1))+  &
                        HRR(1:LB,96,1)
      !=|27,14)
      HRR(1:LB,27,14)=CDy*(CDy*(CDy*HRR(1:LB,27,1)+  &
                        3.D0*HRR(1:LB,44,1))+  &
                        3.D0*HRR(1:LB,67,1))+  &
                        HRR(1:LB,97,1)
      !=|28,14)
      HRR(1:LB,28,14)=CDy*(CDy*(CDy*HRR(1:LB,28,1)+  &
                        3.D0*HRR(1:LB,45,1))+  &
                        3.D0*HRR(1:LB,68,1))+  &
                        HRR(1:LB,98,1)
      !=|29,14)
      HRR(1:LB,29,14)=CDy*(CDy*(CDy*HRR(1:LB,29,1)+  &
                        3.D0*HRR(1:LB,46,1))+  &
                        3.D0*HRR(1:LB,69,1))+  &
                        HRR(1:LB,99,1)
      !=|30,14)
      HRR(1:LB,30,14)=CDy*(CDy*(CDy*HRR(1:LB,30,1)+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        3.D0*HRR(1:LB,72,1))+  &
                        HRR(1:LB,103,1)
      !=|31,14)
      HRR(1:LB,31,14)=CDy*(CDy*(CDy*HRR(1:LB,31,1)+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        3.D0*HRR(1:LB,73,1))+  &
                        HRR(1:LB,104,1)
      !=|32,14)
      HRR(1:LB,32,14)=CDy*(CDy*(CDy*HRR(1:LB,32,1)+  &
                        3.D0*HRR(1:LB,50,1))+  &
                        3.D0*HRR(1:LB,74,1))+  &
                        HRR(1:LB,105,1)
      !=|33,14)
      HRR(1:LB,33,14)=CDy*(CDy*(CDy*HRR(1:LB,33,1)+  &
                        3.D0*HRR(1:LB,52,1))+  &
                        3.D0*HRR(1:LB,77,1))+  &
                        HRR(1:LB,109,1)
      !=|34,14)
      HRR(1:LB,34,14)=CDy*(CDy*(CDy*HRR(1:LB,34,1)+  &
                        3.D0*HRR(1:LB,53,1))+  &
                        3.D0*HRR(1:LB,78,1))+  &
                        HRR(1:LB,110,1)
      !=|35,14)
      HRR(1:LB,35,14)=CDy*(CDy*(CDy*HRR(1:LB,35,1)+  &
                        3.D0*HRR(1:LB,55,1))+  &
                        3.D0*HRR(1:LB,81,1))+  &
                        HRR(1:LB,114,1)
      !=|21,15)
      HRR(1:LB,21,15)=CDz*HRR(1:LB,57,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,36,1)+  &
                        CDx*(CDz*HRR(1:LB,21,1)+  &
                        HRR(1:LB,42,1))+  &
                        2.D0*HRR(1:LB,64,1))+  &
                        HRR(1:LB,93,1)
      !=|22,15)
      HRR(1:LB,22,15)=CDz*HRR(1:LB,58,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,37,1)+  &
                        CDx*(CDz*HRR(1:LB,22,1)+  &
                        HRR(1:LB,43,1))+  &
                        2.D0*HRR(1:LB,65,1))+  &
                        HRR(1:LB,94,1)
      !=|23,15)
      HRR(1:LB,23,15)=CDz*HRR(1:LB,59,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,38,1)+  &
                        CDx*(CDz*HRR(1:LB,23,1)+  &
                        HRR(1:LB,44,1))+  &
                        2.D0*HRR(1:LB,66,1))+  &
                        HRR(1:LB,95,1)
      !=|24,15)
      HRR(1:LB,24,15)=CDz*HRR(1:LB,60,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,39,1)+  &
                        CDx*(CDz*HRR(1:LB,24,1)+  &
                        HRR(1:LB,45,1))+  &
                        2.D0*HRR(1:LB,67,1))+  &
                        HRR(1:LB,96,1)
      !=|25,15)
      HRR(1:LB,25,15)=CDz*HRR(1:LB,61,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,40,1)+  &
                        CDx*(CDz*HRR(1:LB,25,1)+  &
                        HRR(1:LB,46,1))+  &
                        2.D0*HRR(1:LB,68,1))+  &
                        HRR(1:LB,97,1)
      !=|26,15)
      HRR(1:LB,26,15)=CDz*HRR(1:LB,64,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,42,1)+  &
                        CDx*(CDz*HRR(1:LB,26,1)+  &
                        HRR(1:LB,47,1))+  &
                        2.D0*HRR(1:LB,70,1))+  &
                        HRR(1:LB,100,1)
      !=|27,15)
      HRR(1:LB,27,15)=CDz*HRR(1:LB,65,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,43,1)+  &
                        CDx*(CDz*HRR(1:LB,27,1)+  &
                        HRR(1:LB,48,1))+  &
                        2.D0*HRR(1:LB,71,1))+  &
                        HRR(1:LB,101,1)
      !=|28,15)
      HRR(1:LB,28,15)=CDz*HRR(1:LB,66,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,44,1)+  &
                        CDx*(CDz*HRR(1:LB,28,1)+  &
                        HRR(1:LB,49,1))+  &
                        2.D0*HRR(1:LB,72,1))+  &
                        HRR(1:LB,102,1)
      !=|29,15)
      HRR(1:LB,29,15)=CDz*HRR(1:LB,67,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,45,1)+  &
                        CDx*(CDz*HRR(1:LB,29,1)+  &
                        HRR(1:LB,50,1))+  &
                        2.D0*HRR(1:LB,73,1))+  &
                        HRR(1:LB,103,1)
      !=|30,15)
      HRR(1:LB,30,15)=CDz*HRR(1:LB,70,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,47,1)+  &
                        CDx*(CDz*HRR(1:LB,30,1)+  &
                        HRR(1:LB,51,1))+  &
                        2.D0*HRR(1:LB,75,1))+  &
                        HRR(1:LB,106,1)
      !=|31,15)
      HRR(1:LB,31,15)=CDz*HRR(1:LB,71,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,48,1)+  &
                        CDx*(CDz*HRR(1:LB,31,1)+  &
                        HRR(1:LB,52,1))+  &
                        2.D0*HRR(1:LB,76,1))+  &
                        HRR(1:LB,107,1)
      !=|32,15)
      HRR(1:LB,32,15)=CDz*HRR(1:LB,72,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,49,1)+  &
                        CDx*(CDz*HRR(1:LB,32,1)+  &
                        HRR(1:LB,53,1))+  &
                        2.D0*HRR(1:LB,77,1))+  &
                        HRR(1:LB,108,1)
      !=|33,15)
      HRR(1:LB,33,15)=CDz*HRR(1:LB,75,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,51,1)+  &
                        CDx*(CDz*HRR(1:LB,33,1)+  &
                        HRR(1:LB,54,1))+  &
                        2.D0*HRR(1:LB,79,1))+  &
                        HRR(1:LB,111,1)
      !=|34,15)
      HRR(1:LB,34,15)=CDz*HRR(1:LB,76,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,52,1)+  &
                        CDx*(CDz*HRR(1:LB,34,1)+  &
                        HRR(1:LB,55,1))+  &
                        2.D0*HRR(1:LB,80,1))+  &
                        HRR(1:LB,112,1)
      !=|35,15)
      HRR(1:LB,35,15)=CDz*HRR(1:LB,79,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,54,1)+  &
                        CDx*(CDz*HRR(1:LB,35,1)+  &
                        HRR(1:LB,56,1))+  &
                        2.D0*HRR(1:LB,82,1))+  &
                        HRR(1:LB,115,1)
      !=|21,16)
      HRR(1:LB,21,16)=CDz*HRR(1:LB,58,1)+  &
                        CDy*(CDz*HRR(1:LB,36,1)+  &
                        HRR(1:LB,64,1))+  &
                        CDx*(CDz*HRR(1:LB,37,1)+  &
                        CDy*(CDz*HRR(1:LB,21,1)+  &
                        HRR(1:LB,42,1))+  &
                        HRR(1:LB,65,1))+  &
                        HRR(1:LB,94,1)
      !=|22,16)
      HRR(1:LB,22,16)=CDz*HRR(1:LB,59,1)+  &
                        CDy*(CDz*HRR(1:LB,37,1)+  &
                        HRR(1:LB,65,1))+  &
                        CDx*(CDz*HRR(1:LB,38,1)+  &
                        CDy*(CDz*HRR(1:LB,22,1)+  &
                        HRR(1:LB,43,1))+  &
                        HRR(1:LB,66,1))+  &
                        HRR(1:LB,95,1)
      !=|23,16)
      HRR(1:LB,23,16)=CDz*HRR(1:LB,60,1)+  &
                        CDy*(CDz*HRR(1:LB,38,1)+  &
                        HRR(1:LB,66,1))+  &
                        CDx*(CDz*HRR(1:LB,39,1)+  &
                        CDy*(CDz*HRR(1:LB,23,1)+  &
                        HRR(1:LB,44,1))+  &
                        HRR(1:LB,67,1))+  &
                        HRR(1:LB,96,1)
      !=|24,16)
      HRR(1:LB,24,16)=CDz*HRR(1:LB,61,1)+  &
                        CDy*(CDz*HRR(1:LB,39,1)+  &
                        HRR(1:LB,67,1))+  &
                        CDx*(CDz*HRR(1:LB,40,1)+  &
                        CDy*(CDz*HRR(1:LB,24,1)+  &
                        HRR(1:LB,45,1))+  &
                        HRR(1:LB,68,1))+  &
                        HRR(1:LB,97,1)
      !=|25,16)
      HRR(1:LB,25,16)=CDz*HRR(1:LB,62,1)+  &
                        CDy*(CDz*HRR(1:LB,40,1)+  &
                        HRR(1:LB,68,1))+  &
                        CDx*(CDz*HRR(1:LB,41,1)+  &
                        CDy*(CDz*HRR(1:LB,25,1)+  &
                        HRR(1:LB,46,1))+  &
                        HRR(1:LB,69,1))+  &
                        HRR(1:LB,98,1)
      !=|26,16)
      HRR(1:LB,26,16)=CDz*HRR(1:LB,65,1)+  &
                        CDy*(CDz*HRR(1:LB,42,1)+  &
                        HRR(1:LB,70,1))+  &
                        CDx*(CDz*HRR(1:LB,43,1)+  &
                        CDy*(CDz*HRR(1:LB,26,1)+  &
                        HRR(1:LB,47,1))+  &
                        HRR(1:LB,71,1))+  &
                        HRR(1:LB,101,1)
      !=|27,16)
      HRR(1:LB,27,16)=CDz*HRR(1:LB,66,1)+  &
                        CDy*(CDz*HRR(1:LB,43,1)+  &
                        HRR(1:LB,71,1))+  &
                        CDx*(CDz*HRR(1:LB,44,1)+  &
                        CDy*(CDz*HRR(1:LB,27,1)+  &
                        HRR(1:LB,48,1))+  &
                        HRR(1:LB,72,1))+  &
                        HRR(1:LB,102,1)
      !=|28,16)
      HRR(1:LB,28,16)=CDz*HRR(1:LB,67,1)+  &
                        CDy*(CDz*HRR(1:LB,44,1)+  &
                        HRR(1:LB,72,1))+  &
                        CDx*(CDz*HRR(1:LB,45,1)+  &
                        CDy*(CDz*HRR(1:LB,28,1)+  &
                        HRR(1:LB,49,1))+  &
                        HRR(1:LB,73,1))+  &
                        HRR(1:LB,103,1)
      !=|29,16)
      HRR(1:LB,29,16)=CDz*HRR(1:LB,68,1)+  &
                        CDy*(CDz*HRR(1:LB,45,1)+  &
                        HRR(1:LB,73,1))+  &
                        CDx*(CDz*HRR(1:LB,46,1)+  &
                        CDy*(CDz*HRR(1:LB,29,1)+  &
                        HRR(1:LB,50,1))+  &
                        HRR(1:LB,74,1))+  &
                        HRR(1:LB,104,1)
      !=|30,16)
      HRR(1:LB,30,16)=CDz*HRR(1:LB,71,1)+  &
                        CDy*(CDz*HRR(1:LB,47,1)+  &
                        HRR(1:LB,75,1))+  &
                        CDx*(CDz*HRR(1:LB,48,1)+  &
                        CDy*(CDz*HRR(1:LB,30,1)+  &
                        HRR(1:LB,51,1))+  &
                        HRR(1:LB,76,1))+  &
                        HRR(1:LB,107,1)
      !=|31,16)
      HRR(1:LB,31,16)=CDz*HRR(1:LB,72,1)+  &
                        CDy*(CDz*HRR(1:LB,48,1)+  &
                        HRR(1:LB,76,1))+  &
                        CDx*(CDz*HRR(1:LB,49,1)+  &
                        CDy*(CDz*HRR(1:LB,31,1)+  &
                        HRR(1:LB,52,1))+  &
                        HRR(1:LB,77,1))+  &
                        HRR(1:LB,108,1)
      !=|32,16)
      HRR(1:LB,32,16)=CDz*HRR(1:LB,73,1)+  &
                        CDy*(CDz*HRR(1:LB,49,1)+  &
                        HRR(1:LB,77,1))+  &
                        CDx*(CDz*HRR(1:LB,50,1)+  &
                        CDy*(CDz*HRR(1:LB,32,1)+  &
                        HRR(1:LB,53,1))+  &
                        HRR(1:LB,78,1))+  &
                        HRR(1:LB,109,1)
      !=|33,16)
      HRR(1:LB,33,16)=CDz*HRR(1:LB,76,1)+  &
                        CDy*(CDz*HRR(1:LB,51,1)+  &
                        HRR(1:LB,79,1))+  &
                        CDx*(CDz*HRR(1:LB,52,1)+  &
                        CDy*(CDz*HRR(1:LB,33,1)+  &
                        HRR(1:LB,54,1))+  &
                        HRR(1:LB,80,1))+  &
                        HRR(1:LB,112,1)
      !=|34,16)
      HRR(1:LB,34,16)=CDz*HRR(1:LB,77,1)+  &
                        CDy*(CDz*HRR(1:LB,52,1)+  &
                        HRR(1:LB,80,1))+  &
                        CDx*(CDz*HRR(1:LB,53,1)+  &
                        CDy*(CDz*HRR(1:LB,34,1)+  &
                        HRR(1:LB,55,1))+  &
                        HRR(1:LB,81,1))+  &
                        HRR(1:LB,113,1)
      !=|35,16)
      HRR(1:LB,35,16)=CDz*HRR(1:LB,80,1)+  &
                        CDy*(CDz*HRR(1:LB,54,1)+  &
                        HRR(1:LB,82,1))+  &
                        CDx*(CDz*HRR(1:LB,55,1)+  &
                        CDy*(CDz*HRR(1:LB,35,1)+  &
                        HRR(1:LB,56,1))+  &
                        HRR(1:LB,83,1))+  &
                        HRR(1:LB,116,1)
      !=|21,17)
      HRR(1:LB,21,17)=CDz*HRR(1:LB,59,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,37,1)+  &
                        CDy*(CDz*HRR(1:LB,21,1)+  &
                        HRR(1:LB,42,1))+  &
                        2.D0*HRR(1:LB,65,1))+  &
                        HRR(1:LB,95,1)
      !=|22,17)
      HRR(1:LB,22,17)=CDz*HRR(1:LB,60,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,38,1)+  &
                        CDy*(CDz*HRR(1:LB,22,1)+  &
                        HRR(1:LB,43,1))+  &
                        2.D0*HRR(1:LB,66,1))+  &
                        HRR(1:LB,96,1)
      !=|23,17)
      HRR(1:LB,23,17)=CDz*HRR(1:LB,61,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,39,1)+  &
                        CDy*(CDz*HRR(1:LB,23,1)+  &
                        HRR(1:LB,44,1))+  &
                        2.D0*HRR(1:LB,67,1))+  &
                        HRR(1:LB,97,1)
      !=|24,17)
      HRR(1:LB,24,17)=CDz*HRR(1:LB,62,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,40,1)+  &
                        CDy*(CDz*HRR(1:LB,24,1)+  &
                        HRR(1:LB,45,1))+  &
                        2.D0*HRR(1:LB,68,1))+  &
                        HRR(1:LB,98,1)
      !=|25,17)
      HRR(1:LB,25,17)=CDz*HRR(1:LB,63,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,41,1)+  &
                        CDy*(CDz*HRR(1:LB,25,1)+  &
                        HRR(1:LB,46,1))+  &
                        2.D0*HRR(1:LB,69,1))+  &
                        HRR(1:LB,99,1)
      !=|26,17)
      HRR(1:LB,26,17)=CDz*HRR(1:LB,66,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,43,1)+  &
                        CDy*(CDz*HRR(1:LB,26,1)+  &
                        HRR(1:LB,47,1))+  &
                        2.D0*HRR(1:LB,71,1))+  &
                        HRR(1:LB,102,1)
      !=|27,17)
      HRR(1:LB,27,17)=CDz*HRR(1:LB,67,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,44,1)+  &
                        CDy*(CDz*HRR(1:LB,27,1)+  &
                        HRR(1:LB,48,1))+  &
                        2.D0*HRR(1:LB,72,1))+  &
                        HRR(1:LB,103,1)
      !=|28,17)
      HRR(1:LB,28,17)=CDz*HRR(1:LB,68,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,45,1)+  &
                        CDy*(CDz*HRR(1:LB,28,1)+  &
                        HRR(1:LB,49,1))+  &
                        2.D0*HRR(1:LB,73,1))+  &
                        HRR(1:LB,104,1)
      !=|29,17)
      HRR(1:LB,29,17)=CDz*HRR(1:LB,69,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,46,1)+  &
                        CDy*(CDz*HRR(1:LB,29,1)+  &
                        HRR(1:LB,50,1))+  &
                        2.D0*HRR(1:LB,74,1))+  &
                        HRR(1:LB,105,1)
      !=|30,17)
      HRR(1:LB,30,17)=CDz*HRR(1:LB,72,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,48,1)+  &
                        CDy*(CDz*HRR(1:LB,30,1)+  &
                        HRR(1:LB,51,1))+  &
                        2.D0*HRR(1:LB,76,1))+  &
                        HRR(1:LB,108,1)
      !=|31,17)
      HRR(1:LB,31,17)=CDz*HRR(1:LB,73,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,49,1)+  &
                        CDy*(CDz*HRR(1:LB,31,1)+  &
                        HRR(1:LB,52,1))+  &
                        2.D0*HRR(1:LB,77,1))+  &
                        HRR(1:LB,109,1)
      !=|32,17)
      HRR(1:LB,32,17)=CDz*HRR(1:LB,74,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,50,1)+  &
                        CDy*(CDz*HRR(1:LB,32,1)+  &
                        HRR(1:LB,53,1))+  &
                        2.D0*HRR(1:LB,78,1))+  &
                        HRR(1:LB,110,1)
      !=|33,17)
      HRR(1:LB,33,17)=CDz*HRR(1:LB,77,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,52,1)+  &
                        CDy*(CDz*HRR(1:LB,33,1)+  &
                        HRR(1:LB,54,1))+  &
                        2.D0*HRR(1:LB,80,1))+  &
                        HRR(1:LB,113,1)
      !=|34,17)
      HRR(1:LB,34,17)=CDz*HRR(1:LB,78,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,53,1)+  &
                        CDy*(CDz*HRR(1:LB,34,1)+  &
                        HRR(1:LB,55,1))+  &
                        2.D0*HRR(1:LB,81,1))+  &
                        HRR(1:LB,114,1)
      !=|35,17)
      HRR(1:LB,35,17)=CDz*HRR(1:LB,81,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,55,1)+  &
                        CDy*(CDz*HRR(1:LB,35,1)+  &
                        HRR(1:LB,56,1))+  &
                        2.D0*HRR(1:LB,83,1))+  &
                        HRR(1:LB,117,1)
      !=|21,18)
      HRR(1:LB,21,18)=CDz*(CDz*HRR(1:LB,36,1)+  &
                        2.D0*HRR(1:LB,64,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,21,1)+  &
                        2.D0*HRR(1:LB,42,1))+  &
                        HRR(1:LB,70,1))+  &
                        HRR(1:LB,100,1)
      !=|22,18)
      HRR(1:LB,22,18)=CDz*(CDz*HRR(1:LB,37,1)+  &
                        2.D0*HRR(1:LB,65,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,22,1)+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,71,1))+  &
                        HRR(1:LB,101,1)
      !=|23,18)
      HRR(1:LB,23,18)=CDz*(CDz*HRR(1:LB,38,1)+  &
                        2.D0*HRR(1:LB,66,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,72,1))+  &
                        HRR(1:LB,102,1)
      !=|24,18)
      HRR(1:LB,24,18)=CDz*(CDz*HRR(1:LB,39,1)+  &
                        2.D0*HRR(1:LB,67,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,24,1)+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,73,1))+  &
                        HRR(1:LB,103,1)
      !=|25,18)
      HRR(1:LB,25,18)=CDz*(CDz*HRR(1:LB,40,1)+  &
                        2.D0*HRR(1:LB,68,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,25,1)+  &
                        2.D0*HRR(1:LB,46,1))+  &
                        HRR(1:LB,74,1))+  &
                        HRR(1:LB,104,1)
      !=|26,18)
      HRR(1:LB,26,18)=CDz*(CDz*HRR(1:LB,42,1)+  &
                        2.D0*HRR(1:LB,70,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,26,1)+  &
                        2.D0*HRR(1:LB,47,1))+  &
                        HRR(1:LB,75,1))+  &
                        HRR(1:LB,106,1)
      !=|27,18)
      HRR(1:LB,27,18)=CDz*(CDz*HRR(1:LB,43,1)+  &
                        2.D0*HRR(1:LB,71,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,27,1)+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,76,1))+  &
                        HRR(1:LB,107,1)
      !=|28,18)
      HRR(1:LB,28,18)=CDz*(CDz*HRR(1:LB,44,1)+  &
                        2.D0*HRR(1:LB,72,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,28,1)+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,77,1))+  &
                        HRR(1:LB,108,1)
      !=|29,18)
      HRR(1:LB,29,18)=CDz*(CDz*HRR(1:LB,45,1)+  &
                        2.D0*HRR(1:LB,73,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,29,1)+  &
                        2.D0*HRR(1:LB,50,1))+  &
                        HRR(1:LB,78,1))+  &
                        HRR(1:LB,109,1)
      !=|30,18)
      HRR(1:LB,30,18)=CDz*(CDz*HRR(1:LB,47,1)+  &
                        2.D0*HRR(1:LB,75,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,30,1)+  &
                        2.D0*HRR(1:LB,51,1))+  &
                        HRR(1:LB,79,1))+  &
                        HRR(1:LB,111,1)
      !=|31,18)
      HRR(1:LB,31,18)=CDz*(CDz*HRR(1:LB,48,1)+  &
                        2.D0*HRR(1:LB,76,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,31,1)+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,80,1))+  &
                        HRR(1:LB,112,1)
      !=|32,18)
      HRR(1:LB,32,18)=CDz*(CDz*HRR(1:LB,49,1)+  &
                        2.D0*HRR(1:LB,77,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,32,1)+  &
                        2.D0*HRR(1:LB,53,1))+  &
                        HRR(1:LB,81,1))+  &
                        HRR(1:LB,113,1)
      !=|33,18)
      HRR(1:LB,33,18)=CDz*(CDz*HRR(1:LB,51,1)+  &
                        2.D0*HRR(1:LB,79,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,33,1)+  &
                        2.D0*HRR(1:LB,54,1))+  &
                        HRR(1:LB,82,1))+  &
                        HRR(1:LB,115,1)
      !=|34,18)
      HRR(1:LB,34,18)=CDz*(CDz*HRR(1:LB,52,1)+  &
                        2.D0*HRR(1:LB,80,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,34,1)+  &
                        2.D0*HRR(1:LB,55,1))+  &
                        HRR(1:LB,83,1))+  &
                        HRR(1:LB,116,1)
      !=|35,18)
      HRR(1:LB,35,18)=CDz*(CDz*HRR(1:LB,54,1)+  &
                        2.D0*HRR(1:LB,82,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,35,1)+  &
                        2.D0*HRR(1:LB,56,1))+  &
                        HRR(1:LB,84,1))+  &
                        HRR(1:LB,118,1)
      !=|21,19)
      HRR(1:LB,21,19)=CDz*(CDz*HRR(1:LB,37,1)+  &
                        2.D0*HRR(1:LB,65,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,21,1)+  &
                        2.D0*HRR(1:LB,42,1))+  &
                        HRR(1:LB,70,1))+  &
                        HRR(1:LB,101,1)
      !=|22,19)
      HRR(1:LB,22,19)=CDz*(CDz*HRR(1:LB,38,1)+  &
                        2.D0*HRR(1:LB,66,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,22,1)+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,71,1))+  &
                        HRR(1:LB,102,1)
      !=|23,19)
      HRR(1:LB,23,19)=CDz*(CDz*HRR(1:LB,39,1)+  &
                        2.D0*HRR(1:LB,67,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,72,1))+  &
                        HRR(1:LB,103,1)
      !=|24,19)
      HRR(1:LB,24,19)=CDz*(CDz*HRR(1:LB,40,1)+  &
                        2.D0*HRR(1:LB,68,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,24,1)+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,73,1))+  &
                        HRR(1:LB,104,1)
      !=|25,19)
      HRR(1:LB,25,19)=CDz*(CDz*HRR(1:LB,41,1)+  &
                        2.D0*HRR(1:LB,69,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,25,1)+  &
                        2.D0*HRR(1:LB,46,1))+  &
                        HRR(1:LB,74,1))+  &
                        HRR(1:LB,105,1)
      !=|26,19)
      HRR(1:LB,26,19)=CDz*(CDz*HRR(1:LB,43,1)+  &
                        2.D0*HRR(1:LB,71,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,26,1)+  &
                        2.D0*HRR(1:LB,47,1))+  &
                        HRR(1:LB,75,1))+  &
                        HRR(1:LB,107,1)
      !=|27,19)
      HRR(1:LB,27,19)=CDz*(CDz*HRR(1:LB,44,1)+  &
                        2.D0*HRR(1:LB,72,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,27,1)+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,76,1))+  &
                        HRR(1:LB,108,1)
      !=|28,19)
      HRR(1:LB,28,19)=CDz*(CDz*HRR(1:LB,45,1)+  &
                        2.D0*HRR(1:LB,73,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,28,1)+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,77,1))+  &
                        HRR(1:LB,109,1)
      !=|29,19)
      HRR(1:LB,29,19)=CDz*(CDz*HRR(1:LB,46,1)+  &
                        2.D0*HRR(1:LB,74,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,29,1)+  &
                        2.D0*HRR(1:LB,50,1))+  &
                        HRR(1:LB,78,1))+  &
                        HRR(1:LB,110,1)
      !=|30,19)
      HRR(1:LB,30,19)=CDz*(CDz*HRR(1:LB,48,1)+  &
                        2.D0*HRR(1:LB,76,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,30,1)+  &
                        2.D0*HRR(1:LB,51,1))+  &
                        HRR(1:LB,79,1))+  &
                        HRR(1:LB,112,1)
      !=|31,19)
      HRR(1:LB,31,19)=CDz*(CDz*HRR(1:LB,49,1)+  &
                        2.D0*HRR(1:LB,77,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,31,1)+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,80,1))+  &
                        HRR(1:LB,113,1)
      !=|32,19)
      HRR(1:LB,32,19)=CDz*(CDz*HRR(1:LB,50,1)+  &
                        2.D0*HRR(1:LB,78,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,32,1)+  &
                        2.D0*HRR(1:LB,53,1))+  &
                        HRR(1:LB,81,1))+  &
                        HRR(1:LB,114,1)
      !=|33,19)
      HRR(1:LB,33,19)=CDz*(CDz*HRR(1:LB,52,1)+  &
                        2.D0*HRR(1:LB,80,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,33,1)+  &
                        2.D0*HRR(1:LB,54,1))+  &
                        HRR(1:LB,82,1))+  &
                        HRR(1:LB,116,1)
      !=|34,19)
      HRR(1:LB,34,19)=CDz*(CDz*HRR(1:LB,53,1)+  &
                        2.D0*HRR(1:LB,81,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,34,1)+  &
                        2.D0*HRR(1:LB,55,1))+  &
                        HRR(1:LB,83,1))+  &
                        HRR(1:LB,117,1)
      !=|35,19)
      HRR(1:LB,35,19)=CDz*(CDz*HRR(1:LB,55,1)+  &
                        2.D0*HRR(1:LB,83,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,35,1)+  &
                        2.D0*HRR(1:LB,56,1))+  &
                        HRR(1:LB,84,1))+  &
                        HRR(1:LB,119,1)
      !=|21,20)
      HRR(1:LB,21,20)=CDz*(CDz*(CDz*HRR(1:LB,21,1)+  &
                        3.D0*HRR(1:LB,42,1))+  &
                        3.D0*HRR(1:LB,70,1))+  &
                        HRR(1:LB,106,1)
      !=|22,20)
      HRR(1:LB,22,20)=CDz*(CDz*(CDz*HRR(1:LB,22,1)+  &
                        3.D0*HRR(1:LB,43,1))+  &
                        3.D0*HRR(1:LB,71,1))+  &
                        HRR(1:LB,107,1)
      !=|23,20)
      HRR(1:LB,23,20)=CDz*(CDz*(CDz*HRR(1:LB,23,1)+  &
                        3.D0*HRR(1:LB,44,1))+  &
                        3.D0*HRR(1:LB,72,1))+  &
                        HRR(1:LB,108,1)
      !=|24,20)
      HRR(1:LB,24,20)=CDz*(CDz*(CDz*HRR(1:LB,24,1)+  &
                        3.D0*HRR(1:LB,45,1))+  &
                        3.D0*HRR(1:LB,73,1))+  &
                        HRR(1:LB,109,1)
      !=|25,20)
      HRR(1:LB,25,20)=CDz*(CDz*(CDz*HRR(1:LB,25,1)+  &
                        3.D0*HRR(1:LB,46,1))+  &
                        3.D0*HRR(1:LB,74,1))+  &
                        HRR(1:LB,110,1)
      !=|26,20)
      HRR(1:LB,26,20)=CDz*(CDz*(CDz*HRR(1:LB,26,1)+  &
                        3.D0*HRR(1:LB,47,1))+  &
                        3.D0*HRR(1:LB,75,1))+  &
                        HRR(1:LB,111,1)
      !=|27,20)
      HRR(1:LB,27,20)=CDz*(CDz*(CDz*HRR(1:LB,27,1)+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        3.D0*HRR(1:LB,76,1))+  &
                        HRR(1:LB,112,1)
      !=|28,20)
      HRR(1:LB,28,20)=CDz*(CDz*(CDz*HRR(1:LB,28,1)+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        3.D0*HRR(1:LB,77,1))+  &
                        HRR(1:LB,113,1)
      !=|29,20)
      HRR(1:LB,29,20)=CDz*(CDz*(CDz*HRR(1:LB,29,1)+  &
                        3.D0*HRR(1:LB,50,1))+  &
                        3.D0*HRR(1:LB,78,1))+  &
                        HRR(1:LB,114,1)
      !=|30,20)
      HRR(1:LB,30,20)=CDz*(CDz*(CDz*HRR(1:LB,30,1)+  &
                        3.D0*HRR(1:LB,51,1))+  &
                        3.D0*HRR(1:LB,79,1))+  &
                        HRR(1:LB,115,1)
      !=|31,20)
      HRR(1:LB,31,20)=CDz*(CDz*(CDz*HRR(1:LB,31,1)+  &
                        3.D0*HRR(1:LB,52,1))+  &
                        3.D0*HRR(1:LB,80,1))+  &
                        HRR(1:LB,116,1)
      !=|32,20)
      HRR(1:LB,32,20)=CDz*(CDz*(CDz*HRR(1:LB,32,1)+  &
                        3.D0*HRR(1:LB,53,1))+  &
                        3.D0*HRR(1:LB,81,1))+  &
                        HRR(1:LB,117,1)
      !=|33,20)
      HRR(1:LB,33,20)=CDz*(CDz*(CDz*HRR(1:LB,33,1)+  &
                        3.D0*HRR(1:LB,54,1))+  &
                        3.D0*HRR(1:LB,82,1))+  &
                        HRR(1:LB,118,1)
      !=|34,20)
      HRR(1:LB,34,20)=CDz*(CDz*(CDz*HRR(1:LB,34,1)+  &
                        3.D0*HRR(1:LB,55,1))+  &
                        3.D0*HRR(1:LB,83,1))+  &
                        HRR(1:LB,119,1)
      !=|35,20)
      HRR(1:LB,35,20)=CDz*(CDz*(CDz*HRR(1:LB,35,1)+  &
                        3.D0*HRR(1:LB,56,1))+  &
                        3.D0*HRR(1:LB,84,1))+  &
                        HRR(1:LB,120,1)
      !=|11,11)
      HRR(1:LB,11,11)=CDx*(CDx*(CDx*HRR(1:LB,11,1)+  &
                        3.D0*HRR(1:LB,21,1))+  &
                        3.D0*HRR(1:LB,36,1))+  &
                        HRR(1:LB,57,1)
      !=|12,11)
      HRR(1:LB,12,11)=CDx*(CDx*(CDx*HRR(1:LB,12,1)+  &
                        3.D0*HRR(1:LB,22,1))+  &
                        3.D0*HRR(1:LB,37,1))+  &
                        HRR(1:LB,58,1)
      !=|13,11)
      HRR(1:LB,13,11)=CDx*(CDx*(CDx*HRR(1:LB,13,1)+  &
                        3.D0*HRR(1:LB,23,1))+  &
                        3.D0*HRR(1:LB,38,1))+  &
                        HRR(1:LB,59,1)
      !=|14,11)
      HRR(1:LB,14,11)=CDx*(CDx*(CDx*HRR(1:LB,14,1)+  &
                        3.D0*HRR(1:LB,24,1))+  &
                        3.D0*HRR(1:LB,39,1))+  &
                        HRR(1:LB,60,1)
      !=|15,11)
      HRR(1:LB,15,11)=CDx*(CDx*(CDx*HRR(1:LB,15,1)+  &
                        3.D0*HRR(1:LB,26,1))+  &
                        3.D0*HRR(1:LB,42,1))+  &
                        HRR(1:LB,64,1)
      !=|16,11)
      HRR(1:LB,16,11)=CDx*(CDx*(CDx*HRR(1:LB,16,1)+  &
                        3.D0*HRR(1:LB,27,1))+  &
                        3.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,65,1)
      !=|17,11)
      HRR(1:LB,17,11)=CDx*(CDx*(CDx*HRR(1:LB,17,1)+  &
                        3.D0*HRR(1:LB,28,1))+  &
                        3.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,66,1)
      !=|18,11)
      HRR(1:LB,18,11)=CDx*(CDx*(CDx*HRR(1:LB,18,1)+  &
                        3.D0*HRR(1:LB,30,1))+  &
                        3.D0*HRR(1:LB,47,1))+  &
                        HRR(1:LB,70,1)
      !=|19,11)
      HRR(1:LB,19,11)=CDx*(CDx*(CDx*HRR(1:LB,19,1)+  &
                        3.D0*HRR(1:LB,31,1))+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,71,1)
      !=|20,11)
      HRR(1:LB,20,11)=CDx*(CDx*(CDx*HRR(1:LB,20,1)+  &
                        3.D0*HRR(1:LB,33,1))+  &
                        3.D0*HRR(1:LB,51,1))+  &
                        HRR(1:LB,75,1)
      !=|11,12)
      HRR(1:LB,11,12)=CDy*HRR(1:LB,36,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,21,1)+  &
                        CDx*(CDy*HRR(1:LB,11,1)+  &
                        HRR(1:LB,22,1))+  &
                        2.D0*HRR(1:LB,37,1))+  &
                        HRR(1:LB,58,1)
      !=|12,12)
      HRR(1:LB,12,12)=CDy*HRR(1:LB,37,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,22,1)+  &
                        CDx*(CDy*HRR(1:LB,12,1)+  &
                        HRR(1:LB,23,1))+  &
                        2.D0*HRR(1:LB,38,1))+  &
                        HRR(1:LB,59,1)
      !=|13,12)
      HRR(1:LB,13,12)=CDy*HRR(1:LB,38,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,23,1)+  &
                        CDx*(CDy*HRR(1:LB,13,1)+  &
                        HRR(1:LB,24,1))+  &
                        2.D0*HRR(1:LB,39,1))+  &
                        HRR(1:LB,60,1)
      !=|14,12)
      HRR(1:LB,14,12)=CDy*HRR(1:LB,39,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,24,1)+  &
                        CDx*(CDy*HRR(1:LB,14,1)+  &
                        HRR(1:LB,25,1))+  &
                        2.D0*HRR(1:LB,40,1))+  &
                        HRR(1:LB,61,1)
      !=|15,12)
      HRR(1:LB,15,12)=CDy*HRR(1:LB,42,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,26,1)+  &
                        CDx*(CDy*HRR(1:LB,15,1)+  &
                        HRR(1:LB,27,1))+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,65,1)
      !=|16,12)
      HRR(1:LB,16,12)=CDy*HRR(1:LB,43,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,27,1)+  &
                        CDx*(CDy*HRR(1:LB,16,1)+  &
                        HRR(1:LB,28,1))+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,66,1)
      !=|17,12)
      HRR(1:LB,17,12)=CDy*HRR(1:LB,44,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,28,1)+  &
                        CDx*(CDy*HRR(1:LB,17,1)+  &
                        HRR(1:LB,29,1))+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,67,1)
      !=|18,12)
      HRR(1:LB,18,12)=CDy*HRR(1:LB,47,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,30,1)+  &
                        CDx*(CDy*HRR(1:LB,18,1)+  &
                        HRR(1:LB,31,1))+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,71,1)
      !=|19,12)
      HRR(1:LB,19,12)=CDy*HRR(1:LB,48,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,31,1)+  &
                        CDx*(CDy*HRR(1:LB,19,1)+  &
                        HRR(1:LB,32,1))+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,72,1)
      !=|20,12)
      HRR(1:LB,20,12)=CDy*HRR(1:LB,51,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,33,1)+  &
                        CDx*(CDy*HRR(1:LB,20,1)+  &
                        HRR(1:LB,34,1))+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,76,1)
      !=|11,13)
      HRR(1:LB,11,13)=CDy*(CDy*HRR(1:LB,21,1)+  &
                        2.D0*HRR(1:LB,37,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,11,1)+  &
                        2.D0*HRR(1:LB,22,1))+  &
                        HRR(1:LB,38,1))+  &
                        HRR(1:LB,59,1)
      !=|12,13)
      HRR(1:LB,12,13)=CDy*(CDy*HRR(1:LB,22,1)+  &
                        2.D0*HRR(1:LB,38,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,12,1)+  &
                        2.D0*HRR(1:LB,23,1))+  &
                        HRR(1:LB,39,1))+  &
                        HRR(1:LB,60,1)
      !=|13,13)
      HRR(1:LB,13,13)=CDy*(CDy*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,39,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,13,1)+  &
                        2.D0*HRR(1:LB,24,1))+  &
                        HRR(1:LB,40,1))+  &
                        HRR(1:LB,61,1)
      !=|14,13)
      HRR(1:LB,14,13)=CDy*(CDy*HRR(1:LB,24,1)+  &
                        2.D0*HRR(1:LB,40,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,14,1)+  &
                        2.D0*HRR(1:LB,25,1))+  &
                        HRR(1:LB,41,1))+  &
                        HRR(1:LB,62,1)
      !=|15,13)
      HRR(1:LB,15,13)=CDy*(CDy*HRR(1:LB,26,1)+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,15,1)+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        HRR(1:LB,44,1))+  &
                        HRR(1:LB,66,1)
      !=|16,13)
      HRR(1:LB,16,13)=CDy*(CDy*HRR(1:LB,27,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,16,1)+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,45,1))+  &
                        HRR(1:LB,67,1)
      !=|17,13)
      HRR(1:LB,17,13)=CDy*(CDy*HRR(1:LB,28,1)+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,17,1)+  &
                        2.D0*HRR(1:LB,29,1))+  &
                        HRR(1:LB,46,1))+  &
                        HRR(1:LB,68,1)
      !=|18,13)
      HRR(1:LB,18,13)=CDy*(CDy*HRR(1:LB,30,1)+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,18,1)+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,49,1))+  &
                        HRR(1:LB,72,1)
      !=|19,13)
      HRR(1:LB,19,13)=CDy*(CDy*HRR(1:LB,31,1)+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,19,1)+  &
                        2.D0*HRR(1:LB,32,1))+  &
                        HRR(1:LB,50,1))+  &
                        HRR(1:LB,73,1)
      !=|20,13)
      HRR(1:LB,20,13)=CDy*(CDy*HRR(1:LB,33,1)+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,20,1)+  &
                        2.D0*HRR(1:LB,34,1))+  &
                        HRR(1:LB,53,1))+  &
                        HRR(1:LB,77,1)
      !=|11,14)
      HRR(1:LB,11,14)=CDy*(CDy*(CDy*HRR(1:LB,11,1)+  &
                        3.D0*HRR(1:LB,22,1))+  &
                        3.D0*HRR(1:LB,38,1))+  &
                        HRR(1:LB,60,1)
      !=|12,14)
      HRR(1:LB,12,14)=CDy*(CDy*(CDy*HRR(1:LB,12,1)+  &
                        3.D0*HRR(1:LB,23,1))+  &
                        3.D0*HRR(1:LB,39,1))+  &
                        HRR(1:LB,61,1)
      !=|13,14)
      HRR(1:LB,13,14)=CDy*(CDy*(CDy*HRR(1:LB,13,1)+  &
                        3.D0*HRR(1:LB,24,1))+  &
                        3.D0*HRR(1:LB,40,1))+  &
                        HRR(1:LB,62,1)
      !=|14,14)
      HRR(1:LB,14,14)=CDy*(CDy*(CDy*HRR(1:LB,14,1)+  &
                        3.D0*HRR(1:LB,25,1))+  &
                        3.D0*HRR(1:LB,41,1))+  &
                        HRR(1:LB,63,1)
      !=|15,14)
      HRR(1:LB,15,14)=CDy*(CDy*(CDy*HRR(1:LB,15,1)+  &
                        3.D0*HRR(1:LB,27,1))+  &
                        3.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,67,1)
      !=|16,14)
      HRR(1:LB,16,14)=CDy*(CDy*(CDy*HRR(1:LB,16,1)+  &
                        3.D0*HRR(1:LB,28,1))+  &
                        3.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,68,1)
      !=|17,14)
      HRR(1:LB,17,14)=CDy*(CDy*(CDy*HRR(1:LB,17,1)+  &
                        3.D0*HRR(1:LB,29,1))+  &
                        3.D0*HRR(1:LB,46,1))+  &
                        HRR(1:LB,69,1)
      !=|18,14)
      HRR(1:LB,18,14)=CDy*(CDy*(CDy*HRR(1:LB,18,1)+  &
                        3.D0*HRR(1:LB,31,1))+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,73,1)
      !=|19,14)
      HRR(1:LB,19,14)=CDy*(CDy*(CDy*HRR(1:LB,19,1)+  &
                        3.D0*HRR(1:LB,32,1))+  &
                        3.D0*HRR(1:LB,50,1))+  &
                        HRR(1:LB,74,1)
      !=|20,14)
      HRR(1:LB,20,14)=CDy*(CDy*(CDy*HRR(1:LB,20,1)+  &
                        3.D0*HRR(1:LB,34,1))+  &
                        3.D0*HRR(1:LB,53,1))+  &
                        HRR(1:LB,78,1)
      !=|11,15)
      HRR(1:LB,11,15)=CDz*HRR(1:LB,36,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,21,1)+  &
                        CDx*(CDz*HRR(1:LB,11,1)+  &
                        HRR(1:LB,26,1))+  &
                        2.D0*HRR(1:LB,42,1))+  &
                        HRR(1:LB,64,1)
      !=|12,15)
      HRR(1:LB,12,15)=CDz*HRR(1:LB,37,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,22,1)+  &
                        CDx*(CDz*HRR(1:LB,12,1)+  &
                        HRR(1:LB,27,1))+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,65,1)
      !=|13,15)
      HRR(1:LB,13,15)=CDz*HRR(1:LB,38,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,23,1)+  &
                        CDx*(CDz*HRR(1:LB,13,1)+  &
                        HRR(1:LB,28,1))+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,66,1)
      !=|14,15)
      HRR(1:LB,14,15)=CDz*HRR(1:LB,39,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,24,1)+  &
                        CDx*(CDz*HRR(1:LB,14,1)+  &
                        HRR(1:LB,29,1))+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,67,1)
      !=|15,15)
      HRR(1:LB,15,15)=CDz*HRR(1:LB,42,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,26,1)+  &
                        CDx*(CDz*HRR(1:LB,15,1)+  &
                        HRR(1:LB,30,1))+  &
                        2.D0*HRR(1:LB,47,1))+  &
                        HRR(1:LB,70,1)
      !=|16,15)
      HRR(1:LB,16,15)=CDz*HRR(1:LB,43,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,27,1)+  &
                        CDx*(CDz*HRR(1:LB,16,1)+  &
                        HRR(1:LB,31,1))+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,71,1)
      !=|17,15)
      HRR(1:LB,17,15)=CDz*HRR(1:LB,44,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,28,1)+  &
                        CDx*(CDz*HRR(1:LB,17,1)+  &
                        HRR(1:LB,32,1))+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,72,1)
      !=|18,15)
      HRR(1:LB,18,15)=CDz*HRR(1:LB,47,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,30,1)+  &
                        CDx*(CDz*HRR(1:LB,18,1)+  &
                        HRR(1:LB,33,1))+  &
                        2.D0*HRR(1:LB,51,1))+  &
                        HRR(1:LB,75,1)
      !=|19,15)
      HRR(1:LB,19,15)=CDz*HRR(1:LB,48,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,31,1)+  &
                        CDx*(CDz*HRR(1:LB,19,1)+  &
                        HRR(1:LB,34,1))+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,76,1)
      !=|20,15)
      HRR(1:LB,20,15)=CDz*HRR(1:LB,51,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,33,1)+  &
                        CDx*(CDz*HRR(1:LB,20,1)+  &
                        HRR(1:LB,35,1))+  &
                        2.D0*HRR(1:LB,54,1))+  &
                        HRR(1:LB,79,1)
      !=|11,16)
      HRR(1:LB,11,16)=CDz*HRR(1:LB,37,1)+  &
                        CDy*(CDz*HRR(1:LB,21,1)+  &
                        HRR(1:LB,42,1))+  &
                        CDx*(CDz*HRR(1:LB,22,1)+  &
                        CDy*(CDz*HRR(1:LB,11,1)+  &
                        HRR(1:LB,26,1))+  &
                        HRR(1:LB,43,1))+  &
                        HRR(1:LB,65,1)
      !=|12,16)
      HRR(1:LB,12,16)=CDz*HRR(1:LB,38,1)+  &
                        CDy*(CDz*HRR(1:LB,22,1)+  &
                        HRR(1:LB,43,1))+  &
                        CDx*(CDz*HRR(1:LB,23,1)+  &
                        CDy*(CDz*HRR(1:LB,12,1)+  &
                        HRR(1:LB,27,1))+  &
                        HRR(1:LB,44,1))+  &
                        HRR(1:LB,66,1)
      !=|13,16)
      HRR(1:LB,13,16)=CDz*HRR(1:LB,39,1)+  &
                        CDy*(CDz*HRR(1:LB,23,1)+  &
                        HRR(1:LB,44,1))+  &
                        CDx*(CDz*HRR(1:LB,24,1)+  &
                        CDy*(CDz*HRR(1:LB,13,1)+  &
                        HRR(1:LB,28,1))+  &
                        HRR(1:LB,45,1))+  &
                        HRR(1:LB,67,1)
      !=|14,16)
      HRR(1:LB,14,16)=CDz*HRR(1:LB,40,1)+  &
                        CDy*(CDz*HRR(1:LB,24,1)+  &
                        HRR(1:LB,45,1))+  &
                        CDx*(CDz*HRR(1:LB,25,1)+  &
                        CDy*(CDz*HRR(1:LB,14,1)+  &
                        HRR(1:LB,29,1))+  &
                        HRR(1:LB,46,1))+  &
                        HRR(1:LB,68,1)
      !=|15,16)
      HRR(1:LB,15,16)=CDz*HRR(1:LB,43,1)+  &
                        CDy*(CDz*HRR(1:LB,26,1)+  &
                        HRR(1:LB,47,1))+  &
                        CDx*(CDz*HRR(1:LB,27,1)+  &
                        CDy*(CDz*HRR(1:LB,15,1)+  &
                        HRR(1:LB,30,1))+  &
                        HRR(1:LB,48,1))+  &
                        HRR(1:LB,71,1)
      !=|16,16)
      HRR(1:LB,16,16)=CDz*HRR(1:LB,44,1)+  &
                        CDy*(CDz*HRR(1:LB,27,1)+  &
                        HRR(1:LB,48,1))+  &
                        CDx*(CDz*HRR(1:LB,28,1)+  &
                        CDy*(CDz*HRR(1:LB,16,1)+  &
                        HRR(1:LB,31,1))+  &
                        HRR(1:LB,49,1))+  &
                        HRR(1:LB,72,1)
      !=|17,16)
      HRR(1:LB,17,16)=CDz*HRR(1:LB,45,1)+  &
                        CDy*(CDz*HRR(1:LB,28,1)+  &
                        HRR(1:LB,49,1))+  &
                        CDx*(CDz*HRR(1:LB,29,1)+  &
                        CDy*(CDz*HRR(1:LB,17,1)+  &
                        HRR(1:LB,32,1))+  &
                        HRR(1:LB,50,1))+  &
                        HRR(1:LB,73,1)
      !=|18,16)
      HRR(1:LB,18,16)=CDz*HRR(1:LB,48,1)+  &
                        CDy*(CDz*HRR(1:LB,30,1)+  &
                        HRR(1:LB,51,1))+  &
                        CDx*(CDz*HRR(1:LB,31,1)+  &
                        CDy*(CDz*HRR(1:LB,18,1)+  &
                        HRR(1:LB,33,1))+  &
                        HRR(1:LB,52,1))+  &
                        HRR(1:LB,76,1)
      !=|19,16)
      HRR(1:LB,19,16)=CDz*HRR(1:LB,49,1)+  &
                        CDy*(CDz*HRR(1:LB,31,1)+  &
                        HRR(1:LB,52,1))+  &
                        CDx*(CDz*HRR(1:LB,32,1)+  &
                        CDy*(CDz*HRR(1:LB,19,1)+  &
                        HRR(1:LB,34,1))+  &
                        HRR(1:LB,53,1))+  &
                        HRR(1:LB,77,1)
      !=|20,16)
      HRR(1:LB,20,16)=CDz*HRR(1:LB,52,1)+  &
                        CDy*(CDz*HRR(1:LB,33,1)+  &
                        HRR(1:LB,54,1))+  &
                        CDx*(CDz*HRR(1:LB,34,1)+  &
                        CDy*(CDz*HRR(1:LB,20,1)+  &
                        HRR(1:LB,35,1))+  &
                        HRR(1:LB,55,1))+  &
                        HRR(1:LB,80,1)
      !=|11,17)
      HRR(1:LB,11,17)=CDz*HRR(1:LB,38,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,22,1)+  &
                        CDy*(CDz*HRR(1:LB,11,1)+  &
                        HRR(1:LB,26,1))+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,66,1)
      !=|12,17)
      HRR(1:LB,12,17)=CDz*HRR(1:LB,39,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,23,1)+  &
                        CDy*(CDz*HRR(1:LB,12,1)+  &
                        HRR(1:LB,27,1))+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,67,1)
      !=|13,17)
      HRR(1:LB,13,17)=CDz*HRR(1:LB,40,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,24,1)+  &
                        CDy*(CDz*HRR(1:LB,13,1)+  &
                        HRR(1:LB,28,1))+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,68,1)
      !=|14,17)
      HRR(1:LB,14,17)=CDz*HRR(1:LB,41,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,25,1)+  &
                        CDy*(CDz*HRR(1:LB,14,1)+  &
                        HRR(1:LB,29,1))+  &
                        2.D0*HRR(1:LB,46,1))+  &
                        HRR(1:LB,69,1)
      !=|15,17)
      HRR(1:LB,15,17)=CDz*HRR(1:LB,44,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,27,1)+  &
                        CDy*(CDz*HRR(1:LB,15,1)+  &
                        HRR(1:LB,30,1))+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,72,1)
      !=|16,17)
      HRR(1:LB,16,17)=CDz*HRR(1:LB,45,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,28,1)+  &
                        CDy*(CDz*HRR(1:LB,16,1)+  &
                        HRR(1:LB,31,1))+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,73,1)
      !=|17,17)
      HRR(1:LB,17,17)=CDz*HRR(1:LB,46,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,29,1)+  &
                        CDy*(CDz*HRR(1:LB,17,1)+  &
                        HRR(1:LB,32,1))+  &
                        2.D0*HRR(1:LB,50,1))+  &
                        HRR(1:LB,74,1)
      !=|18,17)
      HRR(1:LB,18,17)=CDz*HRR(1:LB,49,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,31,1)+  &
                        CDy*(CDz*HRR(1:LB,18,1)+  &
                        HRR(1:LB,33,1))+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,77,1)
      !=|19,17)
      HRR(1:LB,19,17)=CDz*HRR(1:LB,50,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,32,1)+  &
                        CDy*(CDz*HRR(1:LB,19,1)+  &
                        HRR(1:LB,34,1))+  &
                        2.D0*HRR(1:LB,53,1))+  &
                        HRR(1:LB,78,1)
      !=|20,17)
      HRR(1:LB,20,17)=CDz*HRR(1:LB,53,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,34,1)+  &
                        CDy*(CDz*HRR(1:LB,20,1)+  &
                        HRR(1:LB,35,1))+  &
                        2.D0*HRR(1:LB,55,1))+  &
                        HRR(1:LB,81,1)
      !=|11,18)
      HRR(1:LB,11,18)=CDz*(CDz*HRR(1:LB,21,1)+  &
                        2.D0*HRR(1:LB,42,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,11,1)+  &
                        2.D0*HRR(1:LB,26,1))+  &
                        HRR(1:LB,47,1))+  &
                        HRR(1:LB,70,1)
      !=|12,18)
      HRR(1:LB,12,18)=CDz*(CDz*HRR(1:LB,22,1)+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,12,1)+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        HRR(1:LB,48,1))+  &
                        HRR(1:LB,71,1)
      !=|13,18)
      HRR(1:LB,13,18)=CDz*(CDz*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,13,1)+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,49,1))+  &
                        HRR(1:LB,72,1)
      !=|14,18)
      HRR(1:LB,14,18)=CDz*(CDz*HRR(1:LB,24,1)+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,14,1)+  &
                        2.D0*HRR(1:LB,29,1))+  &
                        HRR(1:LB,50,1))+  &
                        HRR(1:LB,73,1)
      !=|15,18)
      HRR(1:LB,15,18)=CDz*(CDz*HRR(1:LB,26,1)+  &
                        2.D0*HRR(1:LB,47,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,15,1)+  &
                        2.D0*HRR(1:LB,30,1))+  &
                        HRR(1:LB,51,1))+  &
                        HRR(1:LB,75,1)
      !=|16,18)
      HRR(1:LB,16,18)=CDz*(CDz*HRR(1:LB,27,1)+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,16,1)+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,52,1))+  &
                        HRR(1:LB,76,1)
      !=|17,18)
      HRR(1:LB,17,18)=CDz*(CDz*HRR(1:LB,28,1)+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,17,1)+  &
                        2.D0*HRR(1:LB,32,1))+  &
                        HRR(1:LB,53,1))+  &
                        HRR(1:LB,77,1)
      !=|18,18)
      HRR(1:LB,18,18)=CDz*(CDz*HRR(1:LB,30,1)+  &
                        2.D0*HRR(1:LB,51,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,18,1)+  &
                        2.D0*HRR(1:LB,33,1))+  &
                        HRR(1:LB,54,1))+  &
                        HRR(1:LB,79,1)
      !=|19,18)
      HRR(1:LB,19,18)=CDz*(CDz*HRR(1:LB,31,1)+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,19,1)+  &
                        2.D0*HRR(1:LB,34,1))+  &
                        HRR(1:LB,55,1))+  &
                        HRR(1:LB,80,1)
      !=|20,18)
      HRR(1:LB,20,18)=CDz*(CDz*HRR(1:LB,33,1)+  &
                        2.D0*HRR(1:LB,54,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,20,1)+  &
                        2.D0*HRR(1:LB,35,1))+  &
                        HRR(1:LB,56,1))+  &
                        HRR(1:LB,82,1)
      !=|11,19)
      HRR(1:LB,11,19)=CDz*(CDz*HRR(1:LB,22,1)+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,11,1)+  &
                        2.D0*HRR(1:LB,26,1))+  &
                        HRR(1:LB,47,1))+  &
                        HRR(1:LB,71,1)
      !=|12,19)
      HRR(1:LB,12,19)=CDz*(CDz*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,12,1)+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        HRR(1:LB,48,1))+  &
                        HRR(1:LB,72,1)
      !=|13,19)
      HRR(1:LB,13,19)=CDz*(CDz*HRR(1:LB,24,1)+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,13,1)+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,49,1))+  &
                        HRR(1:LB,73,1)
      !=|14,19)
      HRR(1:LB,14,19)=CDz*(CDz*HRR(1:LB,25,1)+  &
                        2.D0*HRR(1:LB,46,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,14,1)+  &
                        2.D0*HRR(1:LB,29,1))+  &
                        HRR(1:LB,50,1))+  &
                        HRR(1:LB,74,1)
      !=|15,19)
      HRR(1:LB,15,19)=CDz*(CDz*HRR(1:LB,27,1)+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,15,1)+  &
                        2.D0*HRR(1:LB,30,1))+  &
                        HRR(1:LB,51,1))+  &
                        HRR(1:LB,76,1)
      !=|16,19)
      HRR(1:LB,16,19)=CDz*(CDz*HRR(1:LB,28,1)+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,16,1)+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,52,1))+  &
                        HRR(1:LB,77,1)
      !=|17,19)
      HRR(1:LB,17,19)=CDz*(CDz*HRR(1:LB,29,1)+  &
                        2.D0*HRR(1:LB,50,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,17,1)+  &
                        2.D0*HRR(1:LB,32,1))+  &
                        HRR(1:LB,53,1))+  &
                        HRR(1:LB,78,1)
      !=|18,19)
      HRR(1:LB,18,19)=CDz*(CDz*HRR(1:LB,31,1)+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,18,1)+  &
                        2.D0*HRR(1:LB,33,1))+  &
                        HRR(1:LB,54,1))+  &
                        HRR(1:LB,80,1)
      !=|19,19)
      HRR(1:LB,19,19)=CDz*(CDz*HRR(1:LB,32,1)+  &
                        2.D0*HRR(1:LB,53,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,19,1)+  &
                        2.D0*HRR(1:LB,34,1))+  &
                        HRR(1:LB,55,1))+  &
                        HRR(1:LB,81,1)
      !=|20,19)
      HRR(1:LB,20,19)=CDz*(CDz*HRR(1:LB,34,1)+  &
                        2.D0*HRR(1:LB,55,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,20,1)+  &
                        2.D0*HRR(1:LB,35,1))+  &
                        HRR(1:LB,56,1))+  &
                        HRR(1:LB,83,1)
      !=|11,20)
      HRR(1:LB,11,20)=CDz*(CDz*(CDz*HRR(1:LB,11,1)+  &
                        3.D0*HRR(1:LB,26,1))+  &
                        3.D0*HRR(1:LB,47,1))+  &
                        HRR(1:LB,75,1)
      !=|12,20)
      HRR(1:LB,12,20)=CDz*(CDz*(CDz*HRR(1:LB,12,1)+  &
                        3.D0*HRR(1:LB,27,1))+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,76,1)
      !=|13,20)
      HRR(1:LB,13,20)=CDz*(CDz*(CDz*HRR(1:LB,13,1)+  &
                        3.D0*HRR(1:LB,28,1))+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,77,1)
      !=|14,20)
      HRR(1:LB,14,20)=CDz*(CDz*(CDz*HRR(1:LB,14,1)+  &
                        3.D0*HRR(1:LB,29,1))+  &
                        3.D0*HRR(1:LB,50,1))+  &
                        HRR(1:LB,78,1)
      !=|15,20)
      HRR(1:LB,15,20)=CDz*(CDz*(CDz*HRR(1:LB,15,1)+  &
                        3.D0*HRR(1:LB,30,1))+  &
                        3.D0*HRR(1:LB,51,1))+  &
                        HRR(1:LB,79,1)
      !=|16,20)
      HRR(1:LB,16,20)=CDz*(CDz*(CDz*HRR(1:LB,16,1)+  &
                        3.D0*HRR(1:LB,31,1))+  &
                        3.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,80,1)
      !=|17,20)
      HRR(1:LB,17,20)=CDz*(CDz*(CDz*HRR(1:LB,17,1)+  &
                        3.D0*HRR(1:LB,32,1))+  &
                        3.D0*HRR(1:LB,53,1))+  &
                        HRR(1:LB,81,1)
      !=|18,20)
      HRR(1:LB,18,20)=CDz*(CDz*(CDz*HRR(1:LB,18,1)+  &
                        3.D0*HRR(1:LB,33,1))+  &
                        3.D0*HRR(1:LB,54,1))+  &
                        HRR(1:LB,82,1)
      !=|19,20)
      HRR(1:LB,19,20)=CDz*(CDz*(CDz*HRR(1:LB,19,1)+  &
                        3.D0*HRR(1:LB,34,1))+  &
                        3.D0*HRR(1:LB,55,1))+  &
                        HRR(1:LB,83,1)
      !=|20,20)
      HRR(1:LB,20,20)=CDz*(CDz*(CDz*HRR(1:LB,20,1)+  &
                        3.D0*HRR(1:LB,35,1))+  &
                        3.D0*HRR(1:LB,56,1))+  &
                        HRR(1:LB,84,1)
      !=|5,11)
      HRR(1:LB,5,11)=CDx*(CDx*(CDx*HRR(1:LB,5,1)+  &
                        3.D0*HRR(1:LB,11,1))+  &
                        3.D0*HRR(1:LB,21,1))+  &
                        HRR(1:LB,36,1)
      !=|6,11)
      HRR(1:LB,6,11)=CDx*(CDx*(CDx*HRR(1:LB,6,1)+  &
                        3.D0*HRR(1:LB,12,1))+  &
                        3.D0*HRR(1:LB,22,1))+  &
                        HRR(1:LB,37,1)
      !=|7,11)
      HRR(1:LB,7,11)=CDx*(CDx*(CDx*HRR(1:LB,7,1)+  &
                        3.D0*HRR(1:LB,13,1))+  &
                        3.D0*HRR(1:LB,23,1))+  &
                        HRR(1:LB,38,1)
      !=|8,11)
      HRR(1:LB,8,11)=CDx*(CDx*(CDx*HRR(1:LB,8,1)+  &
                        3.D0*HRR(1:LB,15,1))+  &
                        3.D0*HRR(1:LB,26,1))+  &
                        HRR(1:LB,42,1)
      !=|9,11)
      HRR(1:LB,9,11)=CDx*(CDx*(CDx*HRR(1:LB,9,1)+  &
                        3.D0*HRR(1:LB,16,1))+  &
                        3.D0*HRR(1:LB,27,1))+  &
                        HRR(1:LB,43,1)
      !=|10,11)
      HRR(1:LB,10,11)=CDx*(CDx*(CDx*HRR(1:LB,10,1)+  &
                        3.D0*HRR(1:LB,18,1))+  &
                        3.D0*HRR(1:LB,30,1))+  &
                        HRR(1:LB,47,1)
      !=|5,12)
      HRR(1:LB,5,12)=CDy*HRR(1:LB,21,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,11,1)+  &
                        CDx*(CDy*HRR(1:LB,5,1)+  &
                        HRR(1:LB,12,1))+  &
                        2.D0*HRR(1:LB,22,1))+  &
                        HRR(1:LB,37,1)
      !=|6,12)
      HRR(1:LB,6,12)=CDy*HRR(1:LB,22,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,12,1)+  &
                        CDx*(CDy*HRR(1:LB,6,1)+  &
                        HRR(1:LB,13,1))+  &
                        2.D0*HRR(1:LB,23,1))+  &
                        HRR(1:LB,38,1)
      !=|7,12)
      HRR(1:LB,7,12)=CDy*HRR(1:LB,23,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,13,1)+  &
                        CDx*(CDy*HRR(1:LB,7,1)+  &
                        HRR(1:LB,14,1))+  &
                        2.D0*HRR(1:LB,24,1))+  &
                        HRR(1:LB,39,1)
      !=|8,12)
      HRR(1:LB,8,12)=CDy*HRR(1:LB,26,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,15,1)+  &
                        CDx*(CDy*HRR(1:LB,8,1)+  &
                        HRR(1:LB,16,1))+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        HRR(1:LB,43,1)
      !=|9,12)
      HRR(1:LB,9,12)=CDy*HRR(1:LB,27,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,16,1)+  &
                        CDx*(CDy*HRR(1:LB,9,1)+  &
                        HRR(1:LB,17,1))+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,44,1)
      !=|10,12)
      HRR(1:LB,10,12)=CDy*HRR(1:LB,30,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,18,1)+  &
                        CDx*(CDy*HRR(1:LB,10,1)+  &
                        HRR(1:LB,19,1))+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,48,1)
      !=|5,13)
      HRR(1:LB,5,13)=CDy*(CDy*HRR(1:LB,11,1)+  &
                        2.D0*HRR(1:LB,22,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,5,1)+  &
                        2.D0*HRR(1:LB,12,1))+  &
                        HRR(1:LB,23,1))+  &
                        HRR(1:LB,38,1)
      !=|6,13)
      HRR(1:LB,6,13)=CDy*(CDy*HRR(1:LB,12,1)+  &
                        2.D0*HRR(1:LB,23,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,6,1)+  &
                        2.D0*HRR(1:LB,13,1))+  &
                        HRR(1:LB,24,1))+  &
                        HRR(1:LB,39,1)
      !=|7,13)
      HRR(1:LB,7,13)=CDy*(CDy*HRR(1:LB,13,1)+  &
                        2.D0*HRR(1:LB,24,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,7,1)+  &
                        2.D0*HRR(1:LB,14,1))+  &
                        HRR(1:LB,25,1))+  &
                        HRR(1:LB,40,1)
      !=|8,13)
      HRR(1:LB,8,13)=CDy*(CDy*HRR(1:LB,15,1)+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,8,1)+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        HRR(1:LB,28,1))+  &
                        HRR(1:LB,44,1)
      !=|9,13)
      HRR(1:LB,9,13)=CDy*(CDy*HRR(1:LB,16,1)+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,9,1)+  &
                        2.D0*HRR(1:LB,17,1))+  &
                        HRR(1:LB,29,1))+  &
                        HRR(1:LB,45,1)
      !=|10,13)
      HRR(1:LB,10,13)=CDy*(CDy*HRR(1:LB,18,1)+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,10,1)+  &
                        2.D0*HRR(1:LB,19,1))+  &
                        HRR(1:LB,32,1))+  &
                        HRR(1:LB,49,1)
      !=|5,14)
      HRR(1:LB,5,14)=CDy*(CDy*(CDy*HRR(1:LB,5,1)+  &
                        3.D0*HRR(1:LB,12,1))+  &
                        3.D0*HRR(1:LB,23,1))+  &
                        HRR(1:LB,39,1)
      !=|6,14)
      HRR(1:LB,6,14)=CDy*(CDy*(CDy*HRR(1:LB,6,1)+  &
                        3.D0*HRR(1:LB,13,1))+  &
                        3.D0*HRR(1:LB,24,1))+  &
                        HRR(1:LB,40,1)
      !=|7,14)
      HRR(1:LB,7,14)=CDy*(CDy*(CDy*HRR(1:LB,7,1)+  &
                        3.D0*HRR(1:LB,14,1))+  &
                        3.D0*HRR(1:LB,25,1))+  &
                        HRR(1:LB,41,1)
      !=|8,14)
      HRR(1:LB,8,14)=CDy*(CDy*(CDy*HRR(1:LB,8,1)+  &
                        3.D0*HRR(1:LB,16,1))+  &
                        3.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,45,1)
      !=|9,14)
      HRR(1:LB,9,14)=CDy*(CDy*(CDy*HRR(1:LB,9,1)+  &
                        3.D0*HRR(1:LB,17,1))+  &
                        3.D0*HRR(1:LB,29,1))+  &
                        HRR(1:LB,46,1)
      !=|10,14)
      HRR(1:LB,10,14)=CDy*(CDy*(CDy*HRR(1:LB,10,1)+  &
                        3.D0*HRR(1:LB,19,1))+  &
                        3.D0*HRR(1:LB,32,1))+  &
                        HRR(1:LB,50,1)
      !=|5,15)
      HRR(1:LB,5,15)=CDz*HRR(1:LB,21,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,11,1)+  &
                        CDx*(CDz*HRR(1:LB,5,1)+  &
                        HRR(1:LB,15,1))+  &
                        2.D0*HRR(1:LB,26,1))+  &
                        HRR(1:LB,42,1)
      !=|6,15)
      HRR(1:LB,6,15)=CDz*HRR(1:LB,22,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,12,1)+  &
                        CDx*(CDz*HRR(1:LB,6,1)+  &
                        HRR(1:LB,16,1))+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        HRR(1:LB,43,1)
      !=|7,15)
      HRR(1:LB,7,15)=CDz*HRR(1:LB,23,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,13,1)+  &
                        CDx*(CDz*HRR(1:LB,7,1)+  &
                        HRR(1:LB,17,1))+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,44,1)
      !=|8,15)
      HRR(1:LB,8,15)=CDz*HRR(1:LB,26,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,15,1)+  &
                        CDx*(CDz*HRR(1:LB,8,1)+  &
                        HRR(1:LB,18,1))+  &
                        2.D0*HRR(1:LB,30,1))+  &
                        HRR(1:LB,47,1)
      !=|9,15)
      HRR(1:LB,9,15)=CDz*HRR(1:LB,27,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,16,1)+  &
                        CDx*(CDz*HRR(1:LB,9,1)+  &
                        HRR(1:LB,19,1))+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,48,1)
      !=|10,15)
      HRR(1:LB,10,15)=CDz*HRR(1:LB,30,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,18,1)+  &
                        CDx*(CDz*HRR(1:LB,10,1)+  &
                        HRR(1:LB,20,1))+  &
                        2.D0*HRR(1:LB,33,1))+  &
                        HRR(1:LB,51,1)
      !=|5,16)
      HRR(1:LB,5,16)=CDz*HRR(1:LB,22,1)+  &
                        CDy*(CDz*HRR(1:LB,11,1)+  &
                        HRR(1:LB,26,1))+  &
                        CDx*(CDz*HRR(1:LB,12,1)+  &
                        CDy*(CDz*HRR(1:LB,5,1)+  &
                        HRR(1:LB,15,1))+  &
                        HRR(1:LB,27,1))+  &
                        HRR(1:LB,43,1)
      !=|6,16)
      HRR(1:LB,6,16)=CDz*HRR(1:LB,23,1)+  &
                        CDy*(CDz*HRR(1:LB,12,1)+  &
                        HRR(1:LB,27,1))+  &
                        CDx*(CDz*HRR(1:LB,13,1)+  &
                        CDy*(CDz*HRR(1:LB,6,1)+  &
                        HRR(1:LB,16,1))+  &
                        HRR(1:LB,28,1))+  &
                        HRR(1:LB,44,1)
      !=|7,16)
      HRR(1:LB,7,16)=CDz*HRR(1:LB,24,1)+  &
                        CDy*(CDz*HRR(1:LB,13,1)+  &
                        HRR(1:LB,28,1))+  &
                        CDx*(CDz*HRR(1:LB,14,1)+  &
                        CDy*(CDz*HRR(1:LB,7,1)+  &
                        HRR(1:LB,17,1))+  &
                        HRR(1:LB,29,1))+  &
                        HRR(1:LB,45,1)
      !=|8,16)
      HRR(1:LB,8,16)=CDz*HRR(1:LB,27,1)+  &
                        CDy*(CDz*HRR(1:LB,15,1)+  &
                        HRR(1:LB,30,1))+  &
                        CDx*(CDz*HRR(1:LB,16,1)+  &
                        CDy*(CDz*HRR(1:LB,8,1)+  &
                        HRR(1:LB,18,1))+  &
                        HRR(1:LB,31,1))+  &
                        HRR(1:LB,48,1)
      !=|9,16)
      HRR(1:LB,9,16)=CDz*HRR(1:LB,28,1)+  &
                        CDy*(CDz*HRR(1:LB,16,1)+  &
                        HRR(1:LB,31,1))+  &
                        CDx*(CDz*HRR(1:LB,17,1)+  &
                        CDy*(CDz*HRR(1:LB,9,1)+  &
                        HRR(1:LB,19,1))+  &
                        HRR(1:LB,32,1))+  &
                        HRR(1:LB,49,1)
      !=|10,16)
      HRR(1:LB,10,16)=CDz*HRR(1:LB,31,1)+  &
                        CDy*(CDz*HRR(1:LB,18,1)+  &
                        HRR(1:LB,33,1))+  &
                        CDx*(CDz*HRR(1:LB,19,1)+  &
                        CDy*(CDz*HRR(1:LB,10,1)+  &
                        HRR(1:LB,20,1))+  &
                        HRR(1:LB,34,1))+  &
                        HRR(1:LB,52,1)
      !=|5,17)
      HRR(1:LB,5,17)=CDz*HRR(1:LB,23,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,12,1)+  &
                        CDy*(CDz*HRR(1:LB,5,1)+  &
                        HRR(1:LB,15,1))+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        HRR(1:LB,44,1)
      !=|6,17)
      HRR(1:LB,6,17)=CDz*HRR(1:LB,24,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,13,1)+  &
                        CDy*(CDz*HRR(1:LB,6,1)+  &
                        HRR(1:LB,16,1))+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,45,1)
      !=|7,17)
      HRR(1:LB,7,17)=CDz*HRR(1:LB,25,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,14,1)+  &
                        CDy*(CDz*HRR(1:LB,7,1)+  &
                        HRR(1:LB,17,1))+  &
                        2.D0*HRR(1:LB,29,1))+  &
                        HRR(1:LB,46,1)
      !=|8,17)
      HRR(1:LB,8,17)=CDz*HRR(1:LB,28,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,16,1)+  &
                        CDy*(CDz*HRR(1:LB,8,1)+  &
                        HRR(1:LB,18,1))+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,49,1)
      !=|9,17)
      HRR(1:LB,9,17)=CDz*HRR(1:LB,29,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,17,1)+  &
                        CDy*(CDz*HRR(1:LB,9,1)+  &
                        HRR(1:LB,19,1))+  &
                        2.D0*HRR(1:LB,32,1))+  &
                        HRR(1:LB,50,1)
      !=|10,17)
      HRR(1:LB,10,17)=CDz*HRR(1:LB,32,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,19,1)+  &
                        CDy*(CDz*HRR(1:LB,10,1)+  &
                        HRR(1:LB,20,1))+  &
                        2.D0*HRR(1:LB,34,1))+  &
                        HRR(1:LB,53,1)
      !=|5,18)
      HRR(1:LB,5,18)=CDz*(CDz*HRR(1:LB,11,1)+  &
                        2.D0*HRR(1:LB,26,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,5,1)+  &
                        2.D0*HRR(1:LB,15,1))+  &
                        HRR(1:LB,30,1))+  &
                        HRR(1:LB,47,1)
      !=|6,18)
      HRR(1:LB,6,18)=CDz*(CDz*HRR(1:LB,12,1)+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,6,1)+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        HRR(1:LB,31,1))+  &
                        HRR(1:LB,48,1)
      !=|7,18)
      HRR(1:LB,7,18)=CDz*(CDz*HRR(1:LB,13,1)+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,7,1)+  &
                        2.D0*HRR(1:LB,17,1))+  &
                        HRR(1:LB,32,1))+  &
                        HRR(1:LB,49,1)
      !=|8,18)
      HRR(1:LB,8,18)=CDz*(CDz*HRR(1:LB,15,1)+  &
                        2.D0*HRR(1:LB,30,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,8,1)+  &
                        2.D0*HRR(1:LB,18,1))+  &
                        HRR(1:LB,33,1))+  &
                        HRR(1:LB,51,1)
      !=|9,18)
      HRR(1:LB,9,18)=CDz*(CDz*HRR(1:LB,16,1)+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,9,1)+  &
                        2.D0*HRR(1:LB,19,1))+  &
                        HRR(1:LB,34,1))+  &
                        HRR(1:LB,52,1)
      !=|10,18)
      HRR(1:LB,10,18)=CDz*(CDz*HRR(1:LB,18,1)+  &
                        2.D0*HRR(1:LB,33,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,10,1)+  &
                        2.D0*HRR(1:LB,20,1))+  &
                        HRR(1:LB,35,1))+  &
                        HRR(1:LB,54,1)
      !=|5,19)
      HRR(1:LB,5,19)=CDz*(CDz*HRR(1:LB,12,1)+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,5,1)+  &
                        2.D0*HRR(1:LB,15,1))+  &
                        HRR(1:LB,30,1))+  &
                        HRR(1:LB,48,1)
      !=|6,19)
      HRR(1:LB,6,19)=CDz*(CDz*HRR(1:LB,13,1)+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,6,1)+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        HRR(1:LB,31,1))+  &
                        HRR(1:LB,49,1)
      !=|7,19)
      HRR(1:LB,7,19)=CDz*(CDz*HRR(1:LB,14,1)+  &
                        2.D0*HRR(1:LB,29,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,7,1)+  &
                        2.D0*HRR(1:LB,17,1))+  &
                        HRR(1:LB,32,1))+  &
                        HRR(1:LB,50,1)
      !=|8,19)
      HRR(1:LB,8,19)=CDz*(CDz*HRR(1:LB,16,1)+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,8,1)+  &
                        2.D0*HRR(1:LB,18,1))+  &
                        HRR(1:LB,33,1))+  &
                        HRR(1:LB,52,1)
      !=|9,19)
      HRR(1:LB,9,19)=CDz*(CDz*HRR(1:LB,17,1)+  &
                        2.D0*HRR(1:LB,32,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,9,1)+  &
                        2.D0*HRR(1:LB,19,1))+  &
                        HRR(1:LB,34,1))+  &
                        HRR(1:LB,53,1)
      !=|10,19)
      HRR(1:LB,10,19)=CDz*(CDz*HRR(1:LB,19,1)+  &
                        2.D0*HRR(1:LB,34,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,10,1)+  &
                        2.D0*HRR(1:LB,20,1))+  &
                        HRR(1:LB,35,1))+  &
                        HRR(1:LB,55,1)
      !=|5,20)
      HRR(1:LB,5,20)=CDz*(CDz*(CDz*HRR(1:LB,5,1)+  &
                        3.D0*HRR(1:LB,15,1))+  &
                        3.D0*HRR(1:LB,30,1))+  &
                        HRR(1:LB,51,1)
      !=|6,20)
      HRR(1:LB,6,20)=CDz*(CDz*(CDz*HRR(1:LB,6,1)+  &
                        3.D0*HRR(1:LB,16,1))+  &
                        3.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,52,1)
      !=|7,20)
      HRR(1:LB,7,20)=CDz*(CDz*(CDz*HRR(1:LB,7,1)+  &
                        3.D0*HRR(1:LB,17,1))+  &
                        3.D0*HRR(1:LB,32,1))+  &
                        HRR(1:LB,53,1)
      !=|8,20)
      HRR(1:LB,8,20)=CDz*(CDz*(CDz*HRR(1:LB,8,1)+  &
                        3.D0*HRR(1:LB,18,1))+  &
                        3.D0*HRR(1:LB,33,1))+  &
                        HRR(1:LB,54,1)
      !=|9,20)
      HRR(1:LB,9,20)=CDz*(CDz*(CDz*HRR(1:LB,9,1)+  &
                        3.D0*HRR(1:LB,19,1))+  &
                        3.D0*HRR(1:LB,34,1))+  &
                        HRR(1:LB,55,1)
      !=|10,20)
      HRR(1:LB,10,20)=CDz*(CDz*(CDz*HRR(1:LB,10,1)+  &
                        3.D0*HRR(1:LB,20,1))+  &
                        3.D0*HRR(1:LB,35,1))+  &
                        HRR(1:LB,56,1)
      !=|2,11)
      HRR(1:LB,2,11)=CDx*(CDx*(CDx*HRR(1:LB,2,1)+  &
                        3.D0*HRR(1:LB,5,1))+  &
                        3.D0*HRR(1:LB,11,1))+  &
                        HRR(1:LB,21,1)
      !=|3,11)
      HRR(1:LB,3,11)=CDx*(CDx*(CDx*HRR(1:LB,3,1)+  &
                        3.D0*HRR(1:LB,6,1))+  &
                        3.D0*HRR(1:LB,12,1))+  &
                        HRR(1:LB,22,1)
      !=|4,11)
      HRR(1:LB,4,11)=CDx*(CDx*(CDx*HRR(1:LB,4,1)+  &
                        3.D0*HRR(1:LB,8,1))+  &
                        3.D0*HRR(1:LB,15,1))+  &
                        HRR(1:LB,26,1)
      !=|2,12)
      HRR(1:LB,2,12)=CDy*HRR(1:LB,11,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,5,1)+  &
                        CDx*(CDy*HRR(1:LB,2,1)+  &
                        HRR(1:LB,6,1))+  &
                        2.D0*HRR(1:LB,12,1))+  &
                        HRR(1:LB,22,1)
      !=|3,12)
      HRR(1:LB,3,12)=CDy*HRR(1:LB,12,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,6,1)+  &
                        CDx*(CDy*HRR(1:LB,3,1)+  &
                        HRR(1:LB,7,1))+  &
                        2.D0*HRR(1:LB,13,1))+  &
                        HRR(1:LB,23,1)
      !=|4,12)
      HRR(1:LB,4,12)=CDy*HRR(1:LB,15,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,8,1)+  &
                        CDx*(CDy*HRR(1:LB,4,1)+  &
                        HRR(1:LB,9,1))+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        HRR(1:LB,27,1)
      !=|2,13)
      HRR(1:LB,2,13)=CDy*(CDy*HRR(1:LB,5,1)+  &
                        2.D0*HRR(1:LB,12,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,2,1)+  &
                        2.D0*HRR(1:LB,6,1))+  &
                        HRR(1:LB,13,1))+  &
                        HRR(1:LB,23,1)
      !=|3,13)
      HRR(1:LB,3,13)=CDy*(CDy*HRR(1:LB,6,1)+  &
                        2.D0*HRR(1:LB,13,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,3,1)+  &
                        2.D0*HRR(1:LB,7,1))+  &
                        HRR(1:LB,14,1))+  &
                        HRR(1:LB,24,1)
      !=|4,13)
      HRR(1:LB,4,13)=CDy*(CDy*HRR(1:LB,8,1)+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,4,1)+  &
                        2.D0*HRR(1:LB,9,1))+  &
                        HRR(1:LB,17,1))+  &
                        HRR(1:LB,28,1)
      !=|2,14)
      HRR(1:LB,2,14)=CDy*(CDy*(CDy*HRR(1:LB,2,1)+  &
                        3.D0*HRR(1:LB,6,1))+  &
                        3.D0*HRR(1:LB,13,1))+  &
                        HRR(1:LB,24,1)
      !=|3,14)
      HRR(1:LB,3,14)=CDy*(CDy*(CDy*HRR(1:LB,3,1)+  &
                        3.D0*HRR(1:LB,7,1))+  &
                        3.D0*HRR(1:LB,14,1))+  &
                        HRR(1:LB,25,1)
      !=|4,14)
      HRR(1:LB,4,14)=CDy*(CDy*(CDy*HRR(1:LB,4,1)+  &
                        3.D0*HRR(1:LB,9,1))+  &
                        3.D0*HRR(1:LB,17,1))+  &
                        HRR(1:LB,29,1)
      !=|2,15)
      HRR(1:LB,2,15)=CDz*HRR(1:LB,11,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,5,1)+  &
                        CDx*(CDz*HRR(1:LB,2,1)+  &
                        HRR(1:LB,8,1))+  &
                        2.D0*HRR(1:LB,15,1))+  &
                        HRR(1:LB,26,1)
      !=|3,15)
      HRR(1:LB,3,15)=CDz*HRR(1:LB,12,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,6,1)+  &
                        CDx*(CDz*HRR(1:LB,3,1)+  &
                        HRR(1:LB,9,1))+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        HRR(1:LB,27,1)
      !=|4,15)
      HRR(1:LB,4,15)=CDz*HRR(1:LB,15,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,8,1)+  &
                        CDx*(CDz*HRR(1:LB,4,1)+  &
                        HRR(1:LB,10,1))+  &
                        2.D0*HRR(1:LB,18,1))+  &
                        HRR(1:LB,30,1)
      !=|2,16)
      HRR(1:LB,2,16)=CDz*HRR(1:LB,12,1)+  &
                        CDy*(CDz*HRR(1:LB,5,1)+  &
                        HRR(1:LB,15,1))+  &
                        CDx*(CDz*HRR(1:LB,6,1)+  &
                        CDy*(CDz*HRR(1:LB,2,1)+  &
                        HRR(1:LB,8,1))+  &
                        HRR(1:LB,16,1))+  &
                        HRR(1:LB,27,1)
      !=|3,16)
      HRR(1:LB,3,16)=CDz*HRR(1:LB,13,1)+  &
                        CDy*(CDz*HRR(1:LB,6,1)+  &
                        HRR(1:LB,16,1))+  &
                        CDx*(CDz*HRR(1:LB,7,1)+  &
                        CDy*(CDz*HRR(1:LB,3,1)+  &
                        HRR(1:LB,9,1))+  &
                        HRR(1:LB,17,1))+  &
                        HRR(1:LB,28,1)
      !=|4,16)
      HRR(1:LB,4,16)=CDz*HRR(1:LB,16,1)+  &
                        CDy*(CDz*HRR(1:LB,8,1)+  &
                        HRR(1:LB,18,1))+  &
                        CDx*(CDz*HRR(1:LB,9,1)+  &
                        CDy*(CDz*HRR(1:LB,4,1)+  &
                        HRR(1:LB,10,1))+  &
                        HRR(1:LB,19,1))+  &
                        HRR(1:LB,31,1)
      !=|2,17)
      HRR(1:LB,2,17)=CDz*HRR(1:LB,13,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,6,1)+  &
                        CDy*(CDz*HRR(1:LB,2,1)+  &
                        HRR(1:LB,8,1))+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        HRR(1:LB,28,1)
      !=|3,17)
      HRR(1:LB,3,17)=CDz*HRR(1:LB,14,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,7,1)+  &
                        CDy*(CDz*HRR(1:LB,3,1)+  &
                        HRR(1:LB,9,1))+  &
                        2.D0*HRR(1:LB,17,1))+  &
                        HRR(1:LB,29,1)
      !=|4,17)
      HRR(1:LB,4,17)=CDz*HRR(1:LB,17,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,9,1)+  &
                        CDy*(CDz*HRR(1:LB,4,1)+  &
                        HRR(1:LB,10,1))+  &
                        2.D0*HRR(1:LB,19,1))+  &
                        HRR(1:LB,32,1)
      !=|2,18)
      HRR(1:LB,2,18)=CDz*(CDz*HRR(1:LB,5,1)+  &
                        2.D0*HRR(1:LB,15,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,2,1)+  &
                        2.D0*HRR(1:LB,8,1))+  &
                        HRR(1:LB,18,1))+  &
                        HRR(1:LB,30,1)
      !=|3,18)
      HRR(1:LB,3,18)=CDz*(CDz*HRR(1:LB,6,1)+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,3,1)+  &
                        2.D0*HRR(1:LB,9,1))+  &
                        HRR(1:LB,19,1))+  &
                        HRR(1:LB,31,1)
      !=|4,18)
      HRR(1:LB,4,18)=CDz*(CDz*HRR(1:LB,8,1)+  &
                        2.D0*HRR(1:LB,18,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,4,1)+  &
                        2.D0*HRR(1:LB,10,1))+  &
                        HRR(1:LB,20,1))+  &
                        HRR(1:LB,33,1)
      !=|2,19)
      HRR(1:LB,2,19)=CDz*(CDz*HRR(1:LB,6,1)+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,2,1)+  &
                        2.D0*HRR(1:LB,8,1))+  &
                        HRR(1:LB,18,1))+  &
                        HRR(1:LB,31,1)
      !=|3,19)
      HRR(1:LB,3,19)=CDz*(CDz*HRR(1:LB,7,1)+  &
                        2.D0*HRR(1:LB,17,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,3,1)+  &
                        2.D0*HRR(1:LB,9,1))+  &
                        HRR(1:LB,19,1))+  &
                        HRR(1:LB,32,1)
      !=|4,19)
      HRR(1:LB,4,19)=CDz*(CDz*HRR(1:LB,9,1)+  &
                        2.D0*HRR(1:LB,19,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,4,1)+  &
                        2.D0*HRR(1:LB,10,1))+  &
                        HRR(1:LB,20,1))+  &
                        HRR(1:LB,34,1)
      !=|2,20)
      HRR(1:LB,2,20)=CDz*(CDz*(CDz*HRR(1:LB,2,1)+  &
                        3.D0*HRR(1:LB,8,1))+  &
                        3.D0*HRR(1:LB,18,1))+  &
                        HRR(1:LB,33,1)
      !=|3,20)
      HRR(1:LB,3,20)=CDz*(CDz*(CDz*HRR(1:LB,3,1)+  &
                        3.D0*HRR(1:LB,9,1))+  &
                        3.D0*HRR(1:LB,19,1))+  &
                        HRR(1:LB,34,1)
      !=|4,20)
      HRR(1:LB,4,20)=CDz*(CDz*(CDz*HRR(1:LB,4,1)+  &
                        3.D0*HRR(1:LB,10,1))+  &
                        3.D0*HRR(1:LB,20,1))+  &
                        HRR(1:LB,35,1)
      !=|1,11)
      HRR(1:LB,1,11)=CDx*(CDx*(CDx*HRR(1:LB,1,1)+  &
                        3.D0*HRR(1:LB,2,1))+  &
                        3.D0*HRR(1:LB,5,1))+  &
                        HRR(1:LB,11,1)
      !=|1,12)
      HRR(1:LB,1,12)=CDy*HRR(1:LB,5,1)+  &
                        CDx*(2.D0*CDy*HRR(1:LB,2,1)+  &
                        CDx*(CDy*HRR(1:LB,1,1)+  &
                        HRR(1:LB,3,1))+  &
                        2.D0*HRR(1:LB,6,1))+  &
                        HRR(1:LB,12,1)
      !=|1,13)
      HRR(1:LB,1,13)=CDy*(CDy*HRR(1:LB,2,1)+  &
                        2.D0*HRR(1:LB,6,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,1,1)+  &
                        2.D0*HRR(1:LB,3,1))+  &
                        HRR(1:LB,7,1))+  &
                        HRR(1:LB,13,1)
      !=|1,14)
      HRR(1:LB,1,14)=CDy*(CDy*(CDy*HRR(1:LB,1,1)+  &
                        3.D0*HRR(1:LB,3,1))+  &
                        3.D0*HRR(1:LB,7,1))+  &
                        HRR(1:LB,14,1)
      !=|1,15)
      HRR(1:LB,1,15)=CDz*HRR(1:LB,5,1)+  &
                        CDx*(2.D0*CDz*HRR(1:LB,2,1)+  &
                        CDx*(CDz*HRR(1:LB,1,1)+  &
                        HRR(1:LB,4,1))+  &
                        2.D0*HRR(1:LB,8,1))+  &
                        HRR(1:LB,15,1)
      !=|1,16)
      HRR(1:LB,1,16)=CDz*HRR(1:LB,6,1)+  &
                        CDy*(CDz*HRR(1:LB,2,1)+  &
                        HRR(1:LB,8,1))+  &
                        CDx*(CDz*HRR(1:LB,3,1)+  &
                        CDy*(CDz*HRR(1:LB,1,1)+  &
                        HRR(1:LB,4,1))+  &
                        HRR(1:LB,9,1))+  &
                        HRR(1:LB,16,1)
      !=|1,17)
      HRR(1:LB,1,17)=CDz*HRR(1:LB,7,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,3,1)+  &
                        CDy*(CDz*HRR(1:LB,1,1)+  &
                        HRR(1:LB,4,1))+  &
                        2.D0*HRR(1:LB,9,1))+  &
                        HRR(1:LB,17,1)
      !=|1,18)
      HRR(1:LB,1,18)=CDz*(CDz*HRR(1:LB,2,1)+  &
                        2.D0*HRR(1:LB,8,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,1,1)+  &
                        2.D0*HRR(1:LB,4,1))+  &
                        HRR(1:LB,10,1))+  &
                        HRR(1:LB,18,1)
      !=|1,19)
      HRR(1:LB,1,19)=CDz*(CDz*HRR(1:LB,3,1)+  &
                        2.D0*HRR(1:LB,9,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,1,1)+  &
                        2.D0*HRR(1:LB,4,1))+  &
                        HRR(1:LB,10,1))+  &
                        HRR(1:LB,19,1)
      !=|1,20)
      HRR(1:LB,1,20)=CDz*(CDz*(CDz*HRR(1:LB,1,1)+  &
                        3.D0*HRR(1:LB,4,1))+  &
                        3.D0*HRR(1:LB,10,1))+  &
                        HRR(1:LB,20,1)
END SUBROUTINE KetHRR1510
