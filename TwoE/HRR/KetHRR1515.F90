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
   SUBROUTINE KetHRR1515(LB,HRR)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB
      REAL(DOUBLE) :: HRR(1:LB,165,35)
      !=|21,21)
      HRR(1:LB,21,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,21,1)+  &
                        4.D0*HRR(1:LB,36,1))+  &
                        6.D0*HRR(1:LB,57,1))+  &
                        4.D0*HRR(1:LB,85,1))+  &
                        HRR(1:LB,121,1)
      !=|22,21)
      HRR(1:LB,22,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,22,1)+  &
                        4.D0*HRR(1:LB,37,1))+  &
                        6.D0*HRR(1:LB,58,1))+  &
                        4.D0*HRR(1:LB,86,1))+  &
                        HRR(1:LB,122,1)
      !=|23,21)
      HRR(1:LB,23,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,23,1)+  &
                        4.D0*HRR(1:LB,38,1))+  &
                        6.D0*HRR(1:LB,59,1))+  &
                        4.D0*HRR(1:LB,87,1))+  &
                        HRR(1:LB,123,1)
      !=|24,21)
      HRR(1:LB,24,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,24,1)+  &
                        4.D0*HRR(1:LB,39,1))+  &
                        6.D0*HRR(1:LB,60,1))+  &
                        4.D0*HRR(1:LB,88,1))+  &
                        HRR(1:LB,124,1)
      !=|25,21)
      HRR(1:LB,25,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,25,1)+  &
                        4.D0*HRR(1:LB,40,1))+  &
                        6.D0*HRR(1:LB,61,1))+  &
                        4.D0*HRR(1:LB,89,1))+  &
                        HRR(1:LB,125,1)
      !=|26,21)
      HRR(1:LB,26,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,26,1)+  &
                        4.D0*HRR(1:LB,42,1))+  &
                        6.D0*HRR(1:LB,64,1))+  &
                        4.D0*HRR(1:LB,93,1))+  &
                        HRR(1:LB,130,1)
      !=|27,21)
      HRR(1:LB,27,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,27,1)+  &
                        4.D0*HRR(1:LB,43,1))+  &
                        6.D0*HRR(1:LB,65,1))+  &
                        4.D0*HRR(1:LB,94,1))+  &
                        HRR(1:LB,131,1)
      !=|28,21)
      HRR(1:LB,28,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,28,1)+  &
                        4.D0*HRR(1:LB,44,1))+  &
                        6.D0*HRR(1:LB,66,1))+  &
                        4.D0*HRR(1:LB,95,1))+  &
                        HRR(1:LB,132,1)
      !=|29,21)
      HRR(1:LB,29,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,29,1)+  &
                        4.D0*HRR(1:LB,45,1))+  &
                        6.D0*HRR(1:LB,67,1))+  &
                        4.D0*HRR(1:LB,96,1))+  &
                        HRR(1:LB,133,1)
      !=|30,21)
      HRR(1:LB,30,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,30,1)+  &
                        4.D0*HRR(1:LB,47,1))+  &
                        6.D0*HRR(1:LB,70,1))+  &
                        4.D0*HRR(1:LB,100,1))+  &
                        HRR(1:LB,138,1)
      !=|31,21)
      HRR(1:LB,31,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,31,1)+  &
                        4.D0*HRR(1:LB,48,1))+  &
                        6.D0*HRR(1:LB,71,1))+  &
                        4.D0*HRR(1:LB,101,1))+  &
                        HRR(1:LB,139,1)
      !=|32,21)
      HRR(1:LB,32,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,32,1)+  &
                        4.D0*HRR(1:LB,49,1))+  &
                        6.D0*HRR(1:LB,72,1))+  &
                        4.D0*HRR(1:LB,102,1))+  &
                        HRR(1:LB,140,1)
      !=|33,21)
      HRR(1:LB,33,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,33,1)+  &
                        4.D0*HRR(1:LB,51,1))+  &
                        6.D0*HRR(1:LB,75,1))+  &
                        4.D0*HRR(1:LB,106,1))+  &
                        HRR(1:LB,145,1)
      !=|34,21)
      HRR(1:LB,34,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,34,1)+  &
                        4.D0*HRR(1:LB,52,1))+  &
                        6.D0*HRR(1:LB,76,1))+  &
                        4.D0*HRR(1:LB,107,1))+  &
                        HRR(1:LB,146,1)
      !=|35,21)
      HRR(1:LB,35,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,35,1)+  &
                        4.D0*HRR(1:LB,54,1))+  &
                        6.D0*HRR(1:LB,79,1))+  &
                        4.D0*HRR(1:LB,111,1))+  &
                        HRR(1:LB,151,1)
      !=|21,22)
      HRR(1:LB,21,22)=CDy*HRR(1:LB,85,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,57,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,36,1)+  &
                        CDx*(CDy*HRR(1:LB,21,1)+  &
                        HRR(1:LB,37,1))+  &
                        3.D0*HRR(1:LB,58,1))+  &
                        3.D0*HRR(1:LB,86,1))+  &
                        HRR(1:LB,122,1)
      !=|22,22)
      HRR(1:LB,22,22)=CDy*HRR(1:LB,86,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,58,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,37,1)+  &
                        CDx*(CDy*HRR(1:LB,22,1)+  &
                        HRR(1:LB,38,1))+  &
                        3.D0*HRR(1:LB,59,1))+  &
                        3.D0*HRR(1:LB,87,1))+  &
                        HRR(1:LB,123,1)
      !=|23,22)
      HRR(1:LB,23,22)=CDy*HRR(1:LB,87,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,59,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,38,1)+  &
                        CDx*(CDy*HRR(1:LB,23,1)+  &
                        HRR(1:LB,39,1))+  &
                        3.D0*HRR(1:LB,60,1))+  &
                        3.D0*HRR(1:LB,88,1))+  &
                        HRR(1:LB,124,1)
      !=|24,22)
      HRR(1:LB,24,22)=CDy*HRR(1:LB,88,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,60,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,39,1)+  &
                        CDx*(CDy*HRR(1:LB,24,1)+  &
                        HRR(1:LB,40,1))+  &
                        3.D0*HRR(1:LB,61,1))+  &
                        3.D0*HRR(1:LB,89,1))+  &
                        HRR(1:LB,125,1)
      !=|25,22)
      HRR(1:LB,25,22)=CDy*HRR(1:LB,89,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,61,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,40,1)+  &
                        CDx*(CDy*HRR(1:LB,25,1)+  &
                        HRR(1:LB,41,1))+  &
                        3.D0*HRR(1:LB,62,1))+  &
                        3.D0*HRR(1:LB,90,1))+  &
                        HRR(1:LB,126,1)
      !=|26,22)
      HRR(1:LB,26,22)=CDy*HRR(1:LB,93,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,64,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,42,1)+  &
                        CDx*(CDy*HRR(1:LB,26,1)+  &
                        HRR(1:LB,43,1))+  &
                        3.D0*HRR(1:LB,65,1))+  &
                        3.D0*HRR(1:LB,94,1))+  &
                        HRR(1:LB,131,1)
      !=|27,22)
      HRR(1:LB,27,22)=CDy*HRR(1:LB,94,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,65,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,43,1)+  &
                        CDx*(CDy*HRR(1:LB,27,1)+  &
                        HRR(1:LB,44,1))+  &
                        3.D0*HRR(1:LB,66,1))+  &
                        3.D0*HRR(1:LB,95,1))+  &
                        HRR(1:LB,132,1)
      !=|28,22)
      HRR(1:LB,28,22)=CDy*HRR(1:LB,95,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,66,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,44,1)+  &
                        CDx*(CDy*HRR(1:LB,28,1)+  &
                        HRR(1:LB,45,1))+  &
                        3.D0*HRR(1:LB,67,1))+  &
                        3.D0*HRR(1:LB,96,1))+  &
                        HRR(1:LB,133,1)
      !=|29,22)
      HRR(1:LB,29,22)=CDy*HRR(1:LB,96,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,67,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,45,1)+  &
                        CDx*(CDy*HRR(1:LB,29,1)+  &
                        HRR(1:LB,46,1))+  &
                        3.D0*HRR(1:LB,68,1))+  &
                        3.D0*HRR(1:LB,97,1))+  &
                        HRR(1:LB,134,1)
      !=|30,22)
      HRR(1:LB,30,22)=CDy*HRR(1:LB,100,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,70,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,47,1)+  &
                        CDx*(CDy*HRR(1:LB,30,1)+  &
                        HRR(1:LB,48,1))+  &
                        3.D0*HRR(1:LB,71,1))+  &
                        3.D0*HRR(1:LB,101,1))+  &
                        HRR(1:LB,139,1)
      !=|31,22)
      HRR(1:LB,31,22)=CDy*HRR(1:LB,101,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,71,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,48,1)+  &
                        CDx*(CDy*HRR(1:LB,31,1)+  &
                        HRR(1:LB,49,1))+  &
                        3.D0*HRR(1:LB,72,1))+  &
                        3.D0*HRR(1:LB,102,1))+  &
                        HRR(1:LB,140,1)
      !=|32,22)
      HRR(1:LB,32,22)=CDy*HRR(1:LB,102,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,72,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,49,1)+  &
                        CDx*(CDy*HRR(1:LB,32,1)+  &
                        HRR(1:LB,50,1))+  &
                        3.D0*HRR(1:LB,73,1))+  &
                        3.D0*HRR(1:LB,103,1))+  &
                        HRR(1:LB,141,1)
      !=|33,22)
      HRR(1:LB,33,22)=CDy*HRR(1:LB,106,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,75,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,51,1)+  &
                        CDx*(CDy*HRR(1:LB,33,1)+  &
                        HRR(1:LB,52,1))+  &
                        3.D0*HRR(1:LB,76,1))+  &
                        3.D0*HRR(1:LB,107,1))+  &
                        HRR(1:LB,146,1)
      !=|34,22)
      HRR(1:LB,34,22)=CDy*HRR(1:LB,107,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,76,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,52,1)+  &
                        CDx*(CDy*HRR(1:LB,34,1)+  &
                        HRR(1:LB,53,1))+  &
                        3.D0*HRR(1:LB,77,1))+  &
                        3.D0*HRR(1:LB,108,1))+  &
                        HRR(1:LB,147,1)
      !=|35,22)
      HRR(1:LB,35,22)=CDy*HRR(1:LB,111,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,79,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,54,1)+  &
                        CDx*(CDy*HRR(1:LB,35,1)+  &
                        HRR(1:LB,55,1))+  &
                        3.D0*HRR(1:LB,80,1))+  &
                        3.D0*HRR(1:LB,112,1))+  &
                        HRR(1:LB,152,1)
      !=|21,23)
      HRR(1:LB,21,23)=CDy*(CDy*HRR(1:LB,57,1)+  &
                        2.D0*HRR(1:LB,86,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,36,1)+  &
                        4.D0*HRR(1:LB,58,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,21,1)+  &
                        2.D0*HRR(1:LB,37,1))+  &
                        HRR(1:LB,59,1))+  &
                        2.D0*HRR(1:LB,87,1))+  &
                        HRR(1:LB,123,1)
      !=|22,23)
      HRR(1:LB,22,23)=CDy*(CDy*HRR(1:LB,58,1)+  &
                        2.D0*HRR(1:LB,87,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,37,1)+  &
                        4.D0*HRR(1:LB,59,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,22,1)+  &
                        2.D0*HRR(1:LB,38,1))+  &
                        HRR(1:LB,60,1))+  &
                        2.D0*HRR(1:LB,88,1))+  &
                        HRR(1:LB,124,1)
      !=|23,23)
      HRR(1:LB,23,23)=CDy*(CDy*HRR(1:LB,59,1)+  &
                        2.D0*HRR(1:LB,88,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,38,1)+  &
                        4.D0*HRR(1:LB,60,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,39,1))+  &
                        HRR(1:LB,61,1))+  &
                        2.D0*HRR(1:LB,89,1))+  &
                        HRR(1:LB,125,1)
      !=|24,23)
      HRR(1:LB,24,23)=CDy*(CDy*HRR(1:LB,60,1)+  &
                        2.D0*HRR(1:LB,89,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,39,1)+  &
                        4.D0*HRR(1:LB,61,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,24,1)+  &
                        2.D0*HRR(1:LB,40,1))+  &
                        HRR(1:LB,62,1))+  &
                        2.D0*HRR(1:LB,90,1))+  &
                        HRR(1:LB,126,1)
      !=|25,23)
      HRR(1:LB,25,23)=CDy*(CDy*HRR(1:LB,61,1)+  &
                        2.D0*HRR(1:LB,90,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,40,1)+  &
                        4.D0*HRR(1:LB,62,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,25,1)+  &
                        2.D0*HRR(1:LB,41,1))+  &
                        HRR(1:LB,63,1))+  &
                        2.D0*HRR(1:LB,91,1))+  &
                        HRR(1:LB,127,1)
      !=|26,23)
      HRR(1:LB,26,23)=CDy*(CDy*HRR(1:LB,64,1)+  &
                        2.D0*HRR(1:LB,94,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,42,1)+  &
                        4.D0*HRR(1:LB,65,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,26,1)+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,66,1))+  &
                        2.D0*HRR(1:LB,95,1))+  &
                        HRR(1:LB,132,1)
      !=|27,23)
      HRR(1:LB,27,23)=CDy*(CDy*HRR(1:LB,65,1)+  &
                        2.D0*HRR(1:LB,95,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,43,1)+  &
                        4.D0*HRR(1:LB,66,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,27,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,67,1))+  &
                        2.D0*HRR(1:LB,96,1))+  &
                        HRR(1:LB,133,1)
      !=|28,23)
      HRR(1:LB,28,23)=CDy*(CDy*HRR(1:LB,66,1)+  &
                        2.D0*HRR(1:LB,96,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,44,1)+  &
                        4.D0*HRR(1:LB,67,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,28,1)+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,68,1))+  &
                        2.D0*HRR(1:LB,97,1))+  &
                        HRR(1:LB,134,1)
      !=|29,23)
      HRR(1:LB,29,23)=CDy*(CDy*HRR(1:LB,67,1)+  &
                        2.D0*HRR(1:LB,97,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,45,1)+  &
                        4.D0*HRR(1:LB,68,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,29,1)+  &
                        2.D0*HRR(1:LB,46,1))+  &
                        HRR(1:LB,69,1))+  &
                        2.D0*HRR(1:LB,98,1))+  &
                        HRR(1:LB,135,1)
      !=|30,23)
      HRR(1:LB,30,23)=CDy*(CDy*HRR(1:LB,70,1)+  &
                        2.D0*HRR(1:LB,101,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,47,1)+  &
                        4.D0*HRR(1:LB,71,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,30,1)+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,72,1))+  &
                        2.D0*HRR(1:LB,102,1))+  &
                        HRR(1:LB,140,1)
      !=|31,23)
      HRR(1:LB,31,23)=CDy*(CDy*HRR(1:LB,71,1)+  &
                        2.D0*HRR(1:LB,102,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,48,1)+  &
                        4.D0*HRR(1:LB,72,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,31,1)+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,73,1))+  &
                        2.D0*HRR(1:LB,103,1))+  &
                        HRR(1:LB,141,1)
      !=|32,23)
      HRR(1:LB,32,23)=CDy*(CDy*HRR(1:LB,72,1)+  &
                        2.D0*HRR(1:LB,103,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,49,1)+  &
                        4.D0*HRR(1:LB,73,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,32,1)+  &
                        2.D0*HRR(1:LB,50,1))+  &
                        HRR(1:LB,74,1))+  &
                        2.D0*HRR(1:LB,104,1))+  &
                        HRR(1:LB,142,1)
      !=|33,23)
      HRR(1:LB,33,23)=CDy*(CDy*HRR(1:LB,75,1)+  &
                        2.D0*HRR(1:LB,107,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,51,1)+  &
                        4.D0*HRR(1:LB,76,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,33,1)+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,77,1))+  &
                        2.D0*HRR(1:LB,108,1))+  &
                        HRR(1:LB,147,1)
      !=|34,23)
      HRR(1:LB,34,23)=CDy*(CDy*HRR(1:LB,76,1)+  &
                        2.D0*HRR(1:LB,108,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,52,1)+  &
                        4.D0*HRR(1:LB,77,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,34,1)+  &
                        2.D0*HRR(1:LB,53,1))+  &
                        HRR(1:LB,78,1))+  &
                        2.D0*HRR(1:LB,109,1))+  &
                        HRR(1:LB,148,1)
      !=|35,23)
      HRR(1:LB,35,23)=CDy*(CDy*HRR(1:LB,79,1)+  &
                        2.D0*HRR(1:LB,112,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,54,1)+  &
                        4.D0*HRR(1:LB,80,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,35,1)+  &
                        2.D0*HRR(1:LB,55,1))+  &
                        HRR(1:LB,81,1))+  &
                        2.D0*HRR(1:LB,113,1))+  &
                        HRR(1:LB,153,1)
      !=|21,24)
      HRR(1:LB,21,24)=CDy*(CDy*(CDy*HRR(1:LB,36,1)+  &
                        3.D0*HRR(1:LB,58,1))+  &
                        3.D0*HRR(1:LB,87,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,21,1)+  &
                        3.D0*HRR(1:LB,37,1))+  &
                        3.D0*HRR(1:LB,59,1))+  &
                        HRR(1:LB,88,1))+  &
                        HRR(1:LB,124,1)
      !=|22,24)
      HRR(1:LB,22,24)=CDy*(CDy*(CDy*HRR(1:LB,37,1)+  &
                        3.D0*HRR(1:LB,59,1))+  &
                        3.D0*HRR(1:LB,88,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,22,1)+  &
                        3.D0*HRR(1:LB,38,1))+  &
                        3.D0*HRR(1:LB,60,1))+  &
                        HRR(1:LB,89,1))+  &
                        HRR(1:LB,125,1)
      !=|23,24)
      HRR(1:LB,23,24)=CDy*(CDy*(CDy*HRR(1:LB,38,1)+  &
                        3.D0*HRR(1:LB,60,1))+  &
                        3.D0*HRR(1:LB,89,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,23,1)+  &
                        3.D0*HRR(1:LB,39,1))+  &
                        3.D0*HRR(1:LB,61,1))+  &
                        HRR(1:LB,90,1))+  &
                        HRR(1:LB,126,1)
      !=|24,24)
      HRR(1:LB,24,24)=CDy*(CDy*(CDy*HRR(1:LB,39,1)+  &
                        3.D0*HRR(1:LB,61,1))+  &
                        3.D0*HRR(1:LB,90,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,24,1)+  &
                        3.D0*HRR(1:LB,40,1))+  &
                        3.D0*HRR(1:LB,62,1))+  &
                        HRR(1:LB,91,1))+  &
                        HRR(1:LB,127,1)
      !=|25,24)
      HRR(1:LB,25,24)=CDy*(CDy*(CDy*HRR(1:LB,40,1)+  &
                        3.D0*HRR(1:LB,62,1))+  &
                        3.D0*HRR(1:LB,91,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,25,1)+  &
                        3.D0*HRR(1:LB,41,1))+  &
                        3.D0*HRR(1:LB,63,1))+  &
                        HRR(1:LB,92,1))+  &
                        HRR(1:LB,128,1)
      !=|26,24)
      HRR(1:LB,26,24)=CDy*(CDy*(CDy*HRR(1:LB,42,1)+  &
                        3.D0*HRR(1:LB,65,1))+  &
                        3.D0*HRR(1:LB,95,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,26,1)+  &
                        3.D0*HRR(1:LB,43,1))+  &
                        3.D0*HRR(1:LB,66,1))+  &
                        HRR(1:LB,96,1))+  &
                        HRR(1:LB,133,1)
      !=|27,24)
      HRR(1:LB,27,24)=CDy*(CDy*(CDy*HRR(1:LB,43,1)+  &
                        3.D0*HRR(1:LB,66,1))+  &
                        3.D0*HRR(1:LB,96,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,27,1)+  &
                        3.D0*HRR(1:LB,44,1))+  &
                        3.D0*HRR(1:LB,67,1))+  &
                        HRR(1:LB,97,1))+  &
                        HRR(1:LB,134,1)
      !=|28,24)
      HRR(1:LB,28,24)=CDy*(CDy*(CDy*HRR(1:LB,44,1)+  &
                        3.D0*HRR(1:LB,67,1))+  &
                        3.D0*HRR(1:LB,97,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,28,1)+  &
                        3.D0*HRR(1:LB,45,1))+  &
                        3.D0*HRR(1:LB,68,1))+  &
                        HRR(1:LB,98,1))+  &
                        HRR(1:LB,135,1)
      !=|29,24)
      HRR(1:LB,29,24)=CDy*(CDy*(CDy*HRR(1:LB,45,1)+  &
                        3.D0*HRR(1:LB,68,1))+  &
                        3.D0*HRR(1:LB,98,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,29,1)+  &
                        3.D0*HRR(1:LB,46,1))+  &
                        3.D0*HRR(1:LB,69,1))+  &
                        HRR(1:LB,99,1))+  &
                        HRR(1:LB,136,1)
      !=|30,24)
      HRR(1:LB,30,24)=CDy*(CDy*(CDy*HRR(1:LB,47,1)+  &
                        3.D0*HRR(1:LB,71,1))+  &
                        3.D0*HRR(1:LB,102,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,30,1)+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        3.D0*HRR(1:LB,72,1))+  &
                        HRR(1:LB,103,1))+  &
                        HRR(1:LB,141,1)
      !=|31,24)
      HRR(1:LB,31,24)=CDy*(CDy*(CDy*HRR(1:LB,48,1)+  &
                        3.D0*HRR(1:LB,72,1))+  &
                        3.D0*HRR(1:LB,103,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,31,1)+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        3.D0*HRR(1:LB,73,1))+  &
                        HRR(1:LB,104,1))+  &
                        HRR(1:LB,142,1)
      !=|32,24)
      HRR(1:LB,32,24)=CDy*(CDy*(CDy*HRR(1:LB,49,1)+  &
                        3.D0*HRR(1:LB,73,1))+  &
                        3.D0*HRR(1:LB,104,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,32,1)+  &
                        3.D0*HRR(1:LB,50,1))+  &
                        3.D0*HRR(1:LB,74,1))+  &
                        HRR(1:LB,105,1))+  &
                        HRR(1:LB,143,1)
      !=|33,24)
      HRR(1:LB,33,24)=CDy*(CDy*(CDy*HRR(1:LB,51,1)+  &
                        3.D0*HRR(1:LB,76,1))+  &
                        3.D0*HRR(1:LB,108,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,33,1)+  &
                        3.D0*HRR(1:LB,52,1))+  &
                        3.D0*HRR(1:LB,77,1))+  &
                        HRR(1:LB,109,1))+  &
                        HRR(1:LB,148,1)
      !=|34,24)
      HRR(1:LB,34,24)=CDy*(CDy*(CDy*HRR(1:LB,52,1)+  &
                        3.D0*HRR(1:LB,77,1))+  &
                        3.D0*HRR(1:LB,109,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,34,1)+  &
                        3.D0*HRR(1:LB,53,1))+  &
                        3.D0*HRR(1:LB,78,1))+  &
                        HRR(1:LB,110,1))+  &
                        HRR(1:LB,149,1)
      !=|35,24)
      HRR(1:LB,35,24)=CDy*(CDy*(CDy*HRR(1:LB,54,1)+  &
                        3.D0*HRR(1:LB,80,1))+  &
                        3.D0*HRR(1:LB,113,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,35,1)+  &
                        3.D0*HRR(1:LB,55,1))+  &
                        3.D0*HRR(1:LB,81,1))+  &
                        HRR(1:LB,114,1))+  &
                        HRR(1:LB,154,1)
      !=|21,25)
      HRR(1:LB,21,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,21,1)+  &
                        4.D0*HRR(1:LB,37,1))+  &
                        6.D0*HRR(1:LB,59,1))+  &
                        4.D0*HRR(1:LB,88,1))+  &
                        HRR(1:LB,125,1)
      !=|22,25)
      HRR(1:LB,22,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,22,1)+  &
                        4.D0*HRR(1:LB,38,1))+  &
                        6.D0*HRR(1:LB,60,1))+  &
                        4.D0*HRR(1:LB,89,1))+  &
                        HRR(1:LB,126,1)
      !=|23,25)
      HRR(1:LB,23,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,23,1)+  &
                        4.D0*HRR(1:LB,39,1))+  &
                        6.D0*HRR(1:LB,61,1))+  &
                        4.D0*HRR(1:LB,90,1))+  &
                        HRR(1:LB,127,1)
      !=|24,25)
      HRR(1:LB,24,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,24,1)+  &
                        4.D0*HRR(1:LB,40,1))+  &
                        6.D0*HRR(1:LB,62,1))+  &
                        4.D0*HRR(1:LB,91,1))+  &
                        HRR(1:LB,128,1)
      !=|25,25)
      HRR(1:LB,25,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,25,1)+  &
                        4.D0*HRR(1:LB,41,1))+  &
                        6.D0*HRR(1:LB,63,1))+  &
                        4.D0*HRR(1:LB,92,1))+  &
                        HRR(1:LB,129,1)
      !=|26,25)
      HRR(1:LB,26,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,26,1)+  &
                        4.D0*HRR(1:LB,43,1))+  &
                        6.D0*HRR(1:LB,66,1))+  &
                        4.D0*HRR(1:LB,96,1))+  &
                        HRR(1:LB,134,1)
      !=|27,25)
      HRR(1:LB,27,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,27,1)+  &
                        4.D0*HRR(1:LB,44,1))+  &
                        6.D0*HRR(1:LB,67,1))+  &
                        4.D0*HRR(1:LB,97,1))+  &
                        HRR(1:LB,135,1)
      !=|28,25)
      HRR(1:LB,28,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,28,1)+  &
                        4.D0*HRR(1:LB,45,1))+  &
                        6.D0*HRR(1:LB,68,1))+  &
                        4.D0*HRR(1:LB,98,1))+  &
                        HRR(1:LB,136,1)
      !=|29,25)
      HRR(1:LB,29,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,29,1)+  &
                        4.D0*HRR(1:LB,46,1))+  &
                        6.D0*HRR(1:LB,69,1))+  &
                        4.D0*HRR(1:LB,99,1))+  &
                        HRR(1:LB,137,1)
      !=|30,25)
      HRR(1:LB,30,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,30,1)+  &
                        4.D0*HRR(1:LB,48,1))+  &
                        6.D0*HRR(1:LB,72,1))+  &
                        4.D0*HRR(1:LB,103,1))+  &
                        HRR(1:LB,142,1)
      !=|31,25)
      HRR(1:LB,31,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,31,1)+  &
                        4.D0*HRR(1:LB,49,1))+  &
                        6.D0*HRR(1:LB,73,1))+  &
                        4.D0*HRR(1:LB,104,1))+  &
                        HRR(1:LB,143,1)
      !=|32,25)
      HRR(1:LB,32,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,32,1)+  &
                        4.D0*HRR(1:LB,50,1))+  &
                        6.D0*HRR(1:LB,74,1))+  &
                        4.D0*HRR(1:LB,105,1))+  &
                        HRR(1:LB,144,1)
      !=|33,25)
      HRR(1:LB,33,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,33,1)+  &
                        4.D0*HRR(1:LB,52,1))+  &
                        6.D0*HRR(1:LB,77,1))+  &
                        4.D0*HRR(1:LB,109,1))+  &
                        HRR(1:LB,149,1)
      !=|34,25)
      HRR(1:LB,34,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,34,1)+  &
                        4.D0*HRR(1:LB,53,1))+  &
                        6.D0*HRR(1:LB,78,1))+  &
                        4.D0*HRR(1:LB,110,1))+  &
                        HRR(1:LB,150,1)
      !=|35,25)
      HRR(1:LB,35,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,35,1)+  &
                        4.D0*HRR(1:LB,55,1))+  &
                        6.D0*HRR(1:LB,81,1))+  &
                        4.D0*HRR(1:LB,114,1))+  &
                        HRR(1:LB,155,1)
      !=|21,26)
      HRR(1:LB,21,26)=CDz*HRR(1:LB,85,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,57,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,36,1)+  &
                        CDx*(CDz*HRR(1:LB,21,1)+  &
                        HRR(1:LB,42,1))+  &
                        3.D0*HRR(1:LB,64,1))+  &
                        3.D0*HRR(1:LB,93,1))+  &
                        HRR(1:LB,130,1)
      !=|22,26)
      HRR(1:LB,22,26)=CDz*HRR(1:LB,86,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,58,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,37,1)+  &
                        CDx*(CDz*HRR(1:LB,22,1)+  &
                        HRR(1:LB,43,1))+  &
                        3.D0*HRR(1:LB,65,1))+  &
                        3.D0*HRR(1:LB,94,1))+  &
                        HRR(1:LB,131,1)
      !=|23,26)
      HRR(1:LB,23,26)=CDz*HRR(1:LB,87,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,59,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,38,1)+  &
                        CDx*(CDz*HRR(1:LB,23,1)+  &
                        HRR(1:LB,44,1))+  &
                        3.D0*HRR(1:LB,66,1))+  &
                        3.D0*HRR(1:LB,95,1))+  &
                        HRR(1:LB,132,1)
      !=|24,26)
      HRR(1:LB,24,26)=CDz*HRR(1:LB,88,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,60,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,39,1)+  &
                        CDx*(CDz*HRR(1:LB,24,1)+  &
                        HRR(1:LB,45,1))+  &
                        3.D0*HRR(1:LB,67,1))+  &
                        3.D0*HRR(1:LB,96,1))+  &
                        HRR(1:LB,133,1)
      !=|25,26)
      HRR(1:LB,25,26)=CDz*HRR(1:LB,89,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,61,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,40,1)+  &
                        CDx*(CDz*HRR(1:LB,25,1)+  &
                        HRR(1:LB,46,1))+  &
                        3.D0*HRR(1:LB,68,1))+  &
                        3.D0*HRR(1:LB,97,1))+  &
                        HRR(1:LB,134,1)
      !=|26,26)
      HRR(1:LB,26,26)=CDz*HRR(1:LB,93,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,64,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,42,1)+  &
                        CDx*(CDz*HRR(1:LB,26,1)+  &
                        HRR(1:LB,47,1))+  &
                        3.D0*HRR(1:LB,70,1))+  &
                        3.D0*HRR(1:LB,100,1))+  &
                        HRR(1:LB,138,1)
      !=|27,26)
      HRR(1:LB,27,26)=CDz*HRR(1:LB,94,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,65,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,43,1)+  &
                        CDx*(CDz*HRR(1:LB,27,1)+  &
                        HRR(1:LB,48,1))+  &
                        3.D0*HRR(1:LB,71,1))+  &
                        3.D0*HRR(1:LB,101,1))+  &
                        HRR(1:LB,139,1)
      !=|28,26)
      HRR(1:LB,28,26)=CDz*HRR(1:LB,95,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,66,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,44,1)+  &
                        CDx*(CDz*HRR(1:LB,28,1)+  &
                        HRR(1:LB,49,1))+  &
                        3.D0*HRR(1:LB,72,1))+  &
                        3.D0*HRR(1:LB,102,1))+  &
                        HRR(1:LB,140,1)
      !=|29,26)
      HRR(1:LB,29,26)=CDz*HRR(1:LB,96,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,67,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,45,1)+  &
                        CDx*(CDz*HRR(1:LB,29,1)+  &
                        HRR(1:LB,50,1))+  &
                        3.D0*HRR(1:LB,73,1))+  &
                        3.D0*HRR(1:LB,103,1))+  &
                        HRR(1:LB,141,1)
      !=|30,26)
      HRR(1:LB,30,26)=CDz*HRR(1:LB,100,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,70,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,47,1)+  &
                        CDx*(CDz*HRR(1:LB,30,1)+  &
                        HRR(1:LB,51,1))+  &
                        3.D0*HRR(1:LB,75,1))+  &
                        3.D0*HRR(1:LB,106,1))+  &
                        HRR(1:LB,145,1)
      !=|31,26)
      HRR(1:LB,31,26)=CDz*HRR(1:LB,101,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,71,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,48,1)+  &
                        CDx*(CDz*HRR(1:LB,31,1)+  &
                        HRR(1:LB,52,1))+  &
                        3.D0*HRR(1:LB,76,1))+  &
                        3.D0*HRR(1:LB,107,1))+  &
                        HRR(1:LB,146,1)
      !=|32,26)
      HRR(1:LB,32,26)=CDz*HRR(1:LB,102,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,72,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,49,1)+  &
                        CDx*(CDz*HRR(1:LB,32,1)+  &
                        HRR(1:LB,53,1))+  &
                        3.D0*HRR(1:LB,77,1))+  &
                        3.D0*HRR(1:LB,108,1))+  &
                        HRR(1:LB,147,1)
      !=|33,26)
      HRR(1:LB,33,26)=CDz*HRR(1:LB,106,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,75,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,51,1)+  &
                        CDx*(CDz*HRR(1:LB,33,1)+  &
                        HRR(1:LB,54,1))+  &
                        3.D0*HRR(1:LB,79,1))+  &
                        3.D0*HRR(1:LB,111,1))+  &
                        HRR(1:LB,151,1)
      !=|34,26)
      HRR(1:LB,34,26)=CDz*HRR(1:LB,107,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,76,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,52,1)+  &
                        CDx*(CDz*HRR(1:LB,34,1)+  &
                        HRR(1:LB,55,1))+  &
                        3.D0*HRR(1:LB,80,1))+  &
                        3.D0*HRR(1:LB,112,1))+  &
                        HRR(1:LB,152,1)
      !=|35,26)
      HRR(1:LB,35,26)=CDz*HRR(1:LB,111,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,79,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,54,1)+  &
                        CDx*(CDz*HRR(1:LB,35,1)+  &
                        HRR(1:LB,56,1))+  &
                        3.D0*HRR(1:LB,82,1))+  &
                        3.D0*HRR(1:LB,115,1))+  &
                        HRR(1:LB,156,1)
      !=|21,27)
      HRR(1:LB,21,27)=CDz*HRR(1:LB,86,1)+  &
                        CDy*(CDz*HRR(1:LB,57,1)+  &
                        HRR(1:LB,93,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,58,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,36,1)+  &
                        2.D0*HRR(1:LB,64,1))+  &
                        CDx*(CDz*HRR(1:LB,37,1)+  &
                        CDy*(CDz*HRR(1:LB,21,1)+  &
                        HRR(1:LB,42,1))+  &
                        HRR(1:LB,65,1))+  &
                        2.D0*HRR(1:LB,94,1))+  &
                        HRR(1:LB,131,1)
      !=|22,27)
      HRR(1:LB,22,27)=CDz*HRR(1:LB,87,1)+  &
                        CDy*(CDz*HRR(1:LB,58,1)+  &
                        HRR(1:LB,94,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,59,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,37,1)+  &
                        2.D0*HRR(1:LB,65,1))+  &
                        CDx*(CDz*HRR(1:LB,38,1)+  &
                        CDy*(CDz*HRR(1:LB,22,1)+  &
                        HRR(1:LB,43,1))+  &
                        HRR(1:LB,66,1))+  &
                        2.D0*HRR(1:LB,95,1))+  &
                        HRR(1:LB,132,1)
      !=|23,27)
      HRR(1:LB,23,27)=CDz*HRR(1:LB,88,1)+  &
                        CDy*(CDz*HRR(1:LB,59,1)+  &
                        HRR(1:LB,95,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,60,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,38,1)+  &
                        2.D0*HRR(1:LB,66,1))+  &
                        CDx*(CDz*HRR(1:LB,39,1)+  &
                        CDy*(CDz*HRR(1:LB,23,1)+  &
                        HRR(1:LB,44,1))+  &
                        HRR(1:LB,67,1))+  &
                        2.D0*HRR(1:LB,96,1))+  &
                        HRR(1:LB,133,1)
      !=|24,27)
      HRR(1:LB,24,27)=CDz*HRR(1:LB,89,1)+  &
                        CDy*(CDz*HRR(1:LB,60,1)+  &
                        HRR(1:LB,96,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,61,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,39,1)+  &
                        2.D0*HRR(1:LB,67,1))+  &
                        CDx*(CDz*HRR(1:LB,40,1)+  &
                        CDy*(CDz*HRR(1:LB,24,1)+  &
                        HRR(1:LB,45,1))+  &
                        HRR(1:LB,68,1))+  &
                        2.D0*HRR(1:LB,97,1))+  &
                        HRR(1:LB,134,1)
      !=|25,27)
      HRR(1:LB,25,27)=CDz*HRR(1:LB,90,1)+  &
                        CDy*(CDz*HRR(1:LB,61,1)+  &
                        HRR(1:LB,97,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,62,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,40,1)+  &
                        2.D0*HRR(1:LB,68,1))+  &
                        CDx*(CDz*HRR(1:LB,41,1)+  &
                        CDy*(CDz*HRR(1:LB,25,1)+  &
                        HRR(1:LB,46,1))+  &
                        HRR(1:LB,69,1))+  &
                        2.D0*HRR(1:LB,98,1))+  &
                        HRR(1:LB,135,1)
      !=|26,27)
      HRR(1:LB,26,27)=CDz*HRR(1:LB,94,1)+  &
                        CDy*(CDz*HRR(1:LB,64,1)+  &
                        HRR(1:LB,100,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,65,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,42,1)+  &
                        2.D0*HRR(1:LB,70,1))+  &
                        CDx*(CDz*HRR(1:LB,43,1)+  &
                        CDy*(CDz*HRR(1:LB,26,1)+  &
                        HRR(1:LB,47,1))+  &
                        HRR(1:LB,71,1))+  &
                        2.D0*HRR(1:LB,101,1))+  &
                        HRR(1:LB,139,1)
      !=|27,27)
      HRR(1:LB,27,27)=CDz*HRR(1:LB,95,1)+  &
                        CDy*(CDz*HRR(1:LB,65,1)+  &
                        HRR(1:LB,101,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,66,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,43,1)+  &
                        2.D0*HRR(1:LB,71,1))+  &
                        CDx*(CDz*HRR(1:LB,44,1)+  &
                        CDy*(CDz*HRR(1:LB,27,1)+  &
                        HRR(1:LB,48,1))+  &
                        HRR(1:LB,72,1))+  &
                        2.D0*HRR(1:LB,102,1))+  &
                        HRR(1:LB,140,1)
      !=|28,27)
      HRR(1:LB,28,27)=CDz*HRR(1:LB,96,1)+  &
                        CDy*(CDz*HRR(1:LB,66,1)+  &
                        HRR(1:LB,102,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,67,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,44,1)+  &
                        2.D0*HRR(1:LB,72,1))+  &
                        CDx*(CDz*HRR(1:LB,45,1)+  &
                        CDy*(CDz*HRR(1:LB,28,1)+  &
                        HRR(1:LB,49,1))+  &
                        HRR(1:LB,73,1))+  &
                        2.D0*HRR(1:LB,103,1))+  &
                        HRR(1:LB,141,1)
      !=|29,27)
      HRR(1:LB,29,27)=CDz*HRR(1:LB,97,1)+  &
                        CDy*(CDz*HRR(1:LB,67,1)+  &
                        HRR(1:LB,103,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,68,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,45,1)+  &
                        2.D0*HRR(1:LB,73,1))+  &
                        CDx*(CDz*HRR(1:LB,46,1)+  &
                        CDy*(CDz*HRR(1:LB,29,1)+  &
                        HRR(1:LB,50,1))+  &
                        HRR(1:LB,74,1))+  &
                        2.D0*HRR(1:LB,104,1))+  &
                        HRR(1:LB,142,1)
      !=|30,27)
      HRR(1:LB,30,27)=CDz*HRR(1:LB,101,1)+  &
                        CDy*(CDz*HRR(1:LB,70,1)+  &
                        HRR(1:LB,106,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,71,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,47,1)+  &
                        2.D0*HRR(1:LB,75,1))+  &
                        CDx*(CDz*HRR(1:LB,48,1)+  &
                        CDy*(CDz*HRR(1:LB,30,1)+  &
                        HRR(1:LB,51,1))+  &
                        HRR(1:LB,76,1))+  &
                        2.D0*HRR(1:LB,107,1))+  &
                        HRR(1:LB,146,1)
      !=|31,27)
      HRR(1:LB,31,27)=CDz*HRR(1:LB,102,1)+  &
                        CDy*(CDz*HRR(1:LB,71,1)+  &
                        HRR(1:LB,107,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,72,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,48,1)+  &
                        2.D0*HRR(1:LB,76,1))+  &
                        CDx*(CDz*HRR(1:LB,49,1)+  &
                        CDy*(CDz*HRR(1:LB,31,1)+  &
                        HRR(1:LB,52,1))+  &
                        HRR(1:LB,77,1))+  &
                        2.D0*HRR(1:LB,108,1))+  &
                        HRR(1:LB,147,1)
      !=|32,27)
      HRR(1:LB,32,27)=CDz*HRR(1:LB,103,1)+  &
                        CDy*(CDz*HRR(1:LB,72,1)+  &
                        HRR(1:LB,108,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,73,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,49,1)+  &
                        2.D0*HRR(1:LB,77,1))+  &
                        CDx*(CDz*HRR(1:LB,50,1)+  &
                        CDy*(CDz*HRR(1:LB,32,1)+  &
                        HRR(1:LB,53,1))+  &
                        HRR(1:LB,78,1))+  &
                        2.D0*HRR(1:LB,109,1))+  &
                        HRR(1:LB,148,1)
      !=|33,27)
      HRR(1:LB,33,27)=CDz*HRR(1:LB,107,1)+  &
                        CDy*(CDz*HRR(1:LB,75,1)+  &
                        HRR(1:LB,111,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,76,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,51,1)+  &
                        2.D0*HRR(1:LB,79,1))+  &
                        CDx*(CDz*HRR(1:LB,52,1)+  &
                        CDy*(CDz*HRR(1:LB,33,1)+  &
                        HRR(1:LB,54,1))+  &
                        HRR(1:LB,80,1))+  &
                        2.D0*HRR(1:LB,112,1))+  &
                        HRR(1:LB,152,1)
      !=|34,27)
      HRR(1:LB,34,27)=CDz*HRR(1:LB,108,1)+  &
                        CDy*(CDz*HRR(1:LB,76,1)+  &
                        HRR(1:LB,112,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,77,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,52,1)+  &
                        2.D0*HRR(1:LB,80,1))+  &
                        CDx*(CDz*HRR(1:LB,53,1)+  &
                        CDy*(CDz*HRR(1:LB,34,1)+  &
                        HRR(1:LB,55,1))+  &
                        HRR(1:LB,81,1))+  &
                        2.D0*HRR(1:LB,113,1))+  &
                        HRR(1:LB,153,1)
      !=|35,27)
      HRR(1:LB,35,27)=CDz*HRR(1:LB,112,1)+  &
                        CDy*(CDz*HRR(1:LB,79,1)+  &
                        HRR(1:LB,115,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,80,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,54,1)+  &
                        2.D0*HRR(1:LB,82,1))+  &
                        CDx*(CDz*HRR(1:LB,55,1)+  &
                        CDy*(CDz*HRR(1:LB,35,1)+  &
                        HRR(1:LB,56,1))+  &
                        HRR(1:LB,83,1))+  &
                        2.D0*HRR(1:LB,116,1))+  &
                        HRR(1:LB,157,1)
      !=|21,28)
      HRR(1:LB,21,28)=CDz*HRR(1:LB,87,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,58,1)+  &
                        CDy*(CDz*HRR(1:LB,36,1)+  &
                        HRR(1:LB,64,1))+  &
                        2.D0*HRR(1:LB,94,1))+  &
                        CDx*(CDz*HRR(1:LB,59,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,37,1)+  &
                        CDy*(CDz*HRR(1:LB,21,1)+  &
                        HRR(1:LB,42,1))+  &
                        2.D0*HRR(1:LB,65,1))+  &
                        HRR(1:LB,95,1))+  &
                        HRR(1:LB,132,1)
      !=|22,28)
      HRR(1:LB,22,28)=CDz*HRR(1:LB,88,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,59,1)+  &
                        CDy*(CDz*HRR(1:LB,37,1)+  &
                        HRR(1:LB,65,1))+  &
                        2.D0*HRR(1:LB,95,1))+  &
                        CDx*(CDz*HRR(1:LB,60,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,38,1)+  &
                        CDy*(CDz*HRR(1:LB,22,1)+  &
                        HRR(1:LB,43,1))+  &
                        2.D0*HRR(1:LB,66,1))+  &
                        HRR(1:LB,96,1))+  &
                        HRR(1:LB,133,1)
      !=|23,28)
      HRR(1:LB,23,28)=CDz*HRR(1:LB,89,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,60,1)+  &
                        CDy*(CDz*HRR(1:LB,38,1)+  &
                        HRR(1:LB,66,1))+  &
                        2.D0*HRR(1:LB,96,1))+  &
                        CDx*(CDz*HRR(1:LB,61,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,39,1)+  &
                        CDy*(CDz*HRR(1:LB,23,1)+  &
                        HRR(1:LB,44,1))+  &
                        2.D0*HRR(1:LB,67,1))+  &
                        HRR(1:LB,97,1))+  &
                        HRR(1:LB,134,1)
      !=|24,28)
      HRR(1:LB,24,28)=CDz*HRR(1:LB,90,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,61,1)+  &
                        CDy*(CDz*HRR(1:LB,39,1)+  &
                        HRR(1:LB,67,1))+  &
                        2.D0*HRR(1:LB,97,1))+  &
                        CDx*(CDz*HRR(1:LB,62,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,40,1)+  &
                        CDy*(CDz*HRR(1:LB,24,1)+  &
                        HRR(1:LB,45,1))+  &
                        2.D0*HRR(1:LB,68,1))+  &
                        HRR(1:LB,98,1))+  &
                        HRR(1:LB,135,1)
      !=|25,28)
      HRR(1:LB,25,28)=CDz*HRR(1:LB,91,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,62,1)+  &
                        CDy*(CDz*HRR(1:LB,40,1)+  &
                        HRR(1:LB,68,1))+  &
                        2.D0*HRR(1:LB,98,1))+  &
                        CDx*(CDz*HRR(1:LB,63,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,41,1)+  &
                        CDy*(CDz*HRR(1:LB,25,1)+  &
                        HRR(1:LB,46,1))+  &
                        2.D0*HRR(1:LB,69,1))+  &
                        HRR(1:LB,99,1))+  &
                        HRR(1:LB,136,1)
      !=|26,28)
      HRR(1:LB,26,28)=CDz*HRR(1:LB,95,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,65,1)+  &
                        CDy*(CDz*HRR(1:LB,42,1)+  &
                        HRR(1:LB,70,1))+  &
                        2.D0*HRR(1:LB,101,1))+  &
                        CDx*(CDz*HRR(1:LB,66,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,43,1)+  &
                        CDy*(CDz*HRR(1:LB,26,1)+  &
                        HRR(1:LB,47,1))+  &
                        2.D0*HRR(1:LB,71,1))+  &
                        HRR(1:LB,102,1))+  &
                        HRR(1:LB,140,1)
      !=|27,28)
      HRR(1:LB,27,28)=CDz*HRR(1:LB,96,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,66,1)+  &
                        CDy*(CDz*HRR(1:LB,43,1)+  &
                        HRR(1:LB,71,1))+  &
                        2.D0*HRR(1:LB,102,1))+  &
                        CDx*(CDz*HRR(1:LB,67,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,44,1)+  &
                        CDy*(CDz*HRR(1:LB,27,1)+  &
                        HRR(1:LB,48,1))+  &
                        2.D0*HRR(1:LB,72,1))+  &
                        HRR(1:LB,103,1))+  &
                        HRR(1:LB,141,1)
      !=|28,28)
      HRR(1:LB,28,28)=CDz*HRR(1:LB,97,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,67,1)+  &
                        CDy*(CDz*HRR(1:LB,44,1)+  &
                        HRR(1:LB,72,1))+  &
                        2.D0*HRR(1:LB,103,1))+  &
                        CDx*(CDz*HRR(1:LB,68,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,45,1)+  &
                        CDy*(CDz*HRR(1:LB,28,1)+  &
                        HRR(1:LB,49,1))+  &
                        2.D0*HRR(1:LB,73,1))+  &
                        HRR(1:LB,104,1))+  &
                        HRR(1:LB,142,1)
      !=|29,28)
      HRR(1:LB,29,28)=CDz*HRR(1:LB,98,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,68,1)+  &
                        CDy*(CDz*HRR(1:LB,45,1)+  &
                        HRR(1:LB,73,1))+  &
                        2.D0*HRR(1:LB,104,1))+  &
                        CDx*(CDz*HRR(1:LB,69,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,46,1)+  &
                        CDy*(CDz*HRR(1:LB,29,1)+  &
                        HRR(1:LB,50,1))+  &
                        2.D0*HRR(1:LB,74,1))+  &
                        HRR(1:LB,105,1))+  &
                        HRR(1:LB,143,1)
      !=|30,28)
      HRR(1:LB,30,28)=CDz*HRR(1:LB,102,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,71,1)+  &
                        CDy*(CDz*HRR(1:LB,47,1)+  &
                        HRR(1:LB,75,1))+  &
                        2.D0*HRR(1:LB,107,1))+  &
                        CDx*(CDz*HRR(1:LB,72,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,48,1)+  &
                        CDy*(CDz*HRR(1:LB,30,1)+  &
                        HRR(1:LB,51,1))+  &
                        2.D0*HRR(1:LB,76,1))+  &
                        HRR(1:LB,108,1))+  &
                        HRR(1:LB,147,1)
      !=|31,28)
      HRR(1:LB,31,28)=CDz*HRR(1:LB,103,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,72,1)+  &
                        CDy*(CDz*HRR(1:LB,48,1)+  &
                        HRR(1:LB,76,1))+  &
                        2.D0*HRR(1:LB,108,1))+  &
                        CDx*(CDz*HRR(1:LB,73,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,49,1)+  &
                        CDy*(CDz*HRR(1:LB,31,1)+  &
                        HRR(1:LB,52,1))+  &
                        2.D0*HRR(1:LB,77,1))+  &
                        HRR(1:LB,109,1))+  &
                        HRR(1:LB,148,1)
      !=|32,28)
      HRR(1:LB,32,28)=CDz*HRR(1:LB,104,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,73,1)+  &
                        CDy*(CDz*HRR(1:LB,49,1)+  &
                        HRR(1:LB,77,1))+  &
                        2.D0*HRR(1:LB,109,1))+  &
                        CDx*(CDz*HRR(1:LB,74,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,50,1)+  &
                        CDy*(CDz*HRR(1:LB,32,1)+  &
                        HRR(1:LB,53,1))+  &
                        2.D0*HRR(1:LB,78,1))+  &
                        HRR(1:LB,110,1))+  &
                        HRR(1:LB,149,1)
      !=|33,28)
      HRR(1:LB,33,28)=CDz*HRR(1:LB,108,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,76,1)+  &
                        CDy*(CDz*HRR(1:LB,51,1)+  &
                        HRR(1:LB,79,1))+  &
                        2.D0*HRR(1:LB,112,1))+  &
                        CDx*(CDz*HRR(1:LB,77,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,52,1)+  &
                        CDy*(CDz*HRR(1:LB,33,1)+  &
                        HRR(1:LB,54,1))+  &
                        2.D0*HRR(1:LB,80,1))+  &
                        HRR(1:LB,113,1))+  &
                        HRR(1:LB,153,1)
      !=|34,28)
      HRR(1:LB,34,28)=CDz*HRR(1:LB,109,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,77,1)+  &
                        CDy*(CDz*HRR(1:LB,52,1)+  &
                        HRR(1:LB,80,1))+  &
                        2.D0*HRR(1:LB,113,1))+  &
                        CDx*(CDz*HRR(1:LB,78,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,53,1)+  &
                        CDy*(CDz*HRR(1:LB,34,1)+  &
                        HRR(1:LB,55,1))+  &
                        2.D0*HRR(1:LB,81,1))+  &
                        HRR(1:LB,114,1))+  &
                        HRR(1:LB,154,1)
      !=|35,28)
      HRR(1:LB,35,28)=CDz*HRR(1:LB,113,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,80,1)+  &
                        CDy*(CDz*HRR(1:LB,54,1)+  &
                        HRR(1:LB,82,1))+  &
                        2.D0*HRR(1:LB,116,1))+  &
                        CDx*(CDz*HRR(1:LB,81,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,55,1)+  &
                        CDy*(CDz*HRR(1:LB,35,1)+  &
                        HRR(1:LB,56,1))+  &
                        2.D0*HRR(1:LB,83,1))+  &
                        HRR(1:LB,117,1))+  &
                        HRR(1:LB,158,1)
      !=|21,29)
      HRR(1:LB,21,29)=CDz*HRR(1:LB,88,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,59,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,37,1)+  &
                        CDy*(CDz*HRR(1:LB,21,1)+  &
                        HRR(1:LB,42,1))+  &
                        3.D0*HRR(1:LB,65,1))+  &
                        3.D0*HRR(1:LB,95,1))+  &
                        HRR(1:LB,133,1)
      !=|22,29)
      HRR(1:LB,22,29)=CDz*HRR(1:LB,89,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,60,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,38,1)+  &
                        CDy*(CDz*HRR(1:LB,22,1)+  &
                        HRR(1:LB,43,1))+  &
                        3.D0*HRR(1:LB,66,1))+  &
                        3.D0*HRR(1:LB,96,1))+  &
                        HRR(1:LB,134,1)
      !=|23,29)
      HRR(1:LB,23,29)=CDz*HRR(1:LB,90,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,61,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,39,1)+  &
                        CDy*(CDz*HRR(1:LB,23,1)+  &
                        HRR(1:LB,44,1))+  &
                        3.D0*HRR(1:LB,67,1))+  &
                        3.D0*HRR(1:LB,97,1))+  &
                        HRR(1:LB,135,1)
      !=|24,29)
      HRR(1:LB,24,29)=CDz*HRR(1:LB,91,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,62,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,40,1)+  &
                        CDy*(CDz*HRR(1:LB,24,1)+  &
                        HRR(1:LB,45,1))+  &
                        3.D0*HRR(1:LB,68,1))+  &
                        3.D0*HRR(1:LB,98,1))+  &
                        HRR(1:LB,136,1)
      !=|25,29)
      HRR(1:LB,25,29)=CDz*HRR(1:LB,92,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,63,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,41,1)+  &
                        CDy*(CDz*HRR(1:LB,25,1)+  &
                        HRR(1:LB,46,1))+  &
                        3.D0*HRR(1:LB,69,1))+  &
                        3.D0*HRR(1:LB,99,1))+  &
                        HRR(1:LB,137,1)
      !=|26,29)
      HRR(1:LB,26,29)=CDz*HRR(1:LB,96,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,66,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,43,1)+  &
                        CDy*(CDz*HRR(1:LB,26,1)+  &
                        HRR(1:LB,47,1))+  &
                        3.D0*HRR(1:LB,71,1))+  &
                        3.D0*HRR(1:LB,102,1))+  &
                        HRR(1:LB,141,1)
      !=|27,29)
      HRR(1:LB,27,29)=CDz*HRR(1:LB,97,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,67,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,44,1)+  &
                        CDy*(CDz*HRR(1:LB,27,1)+  &
                        HRR(1:LB,48,1))+  &
                        3.D0*HRR(1:LB,72,1))+  &
                        3.D0*HRR(1:LB,103,1))+  &
                        HRR(1:LB,142,1)
      !=|28,29)
      HRR(1:LB,28,29)=CDz*HRR(1:LB,98,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,68,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,45,1)+  &
                        CDy*(CDz*HRR(1:LB,28,1)+  &
                        HRR(1:LB,49,1))+  &
                        3.D0*HRR(1:LB,73,1))+  &
                        3.D0*HRR(1:LB,104,1))+  &
                        HRR(1:LB,143,1)
      !=|29,29)
      HRR(1:LB,29,29)=CDz*HRR(1:LB,99,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,69,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,46,1)+  &
                        CDy*(CDz*HRR(1:LB,29,1)+  &
                        HRR(1:LB,50,1))+  &
                        3.D0*HRR(1:LB,74,1))+  &
                        3.D0*HRR(1:LB,105,1))+  &
                        HRR(1:LB,144,1)
      !=|30,29)
      HRR(1:LB,30,29)=CDz*HRR(1:LB,103,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,72,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,48,1)+  &
                        CDy*(CDz*HRR(1:LB,30,1)+  &
                        HRR(1:LB,51,1))+  &
                        3.D0*HRR(1:LB,76,1))+  &
                        3.D0*HRR(1:LB,108,1))+  &
                        HRR(1:LB,148,1)
      !=|31,29)
      HRR(1:LB,31,29)=CDz*HRR(1:LB,104,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,73,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,49,1)+  &
                        CDy*(CDz*HRR(1:LB,31,1)+  &
                        HRR(1:LB,52,1))+  &
                        3.D0*HRR(1:LB,77,1))+  &
                        3.D0*HRR(1:LB,109,1))+  &
                        HRR(1:LB,149,1)
      !=|32,29)
      HRR(1:LB,32,29)=CDz*HRR(1:LB,105,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,74,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,50,1)+  &
                        CDy*(CDz*HRR(1:LB,32,1)+  &
                        HRR(1:LB,53,1))+  &
                        3.D0*HRR(1:LB,78,1))+  &
                        3.D0*HRR(1:LB,110,1))+  &
                        HRR(1:LB,150,1)
      !=|33,29)
      HRR(1:LB,33,29)=CDz*HRR(1:LB,109,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,77,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,52,1)+  &
                        CDy*(CDz*HRR(1:LB,33,1)+  &
                        HRR(1:LB,54,1))+  &
                        3.D0*HRR(1:LB,80,1))+  &
                        3.D0*HRR(1:LB,113,1))+  &
                        HRR(1:LB,154,1)
      !=|34,29)
      HRR(1:LB,34,29)=CDz*HRR(1:LB,110,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,78,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,53,1)+  &
                        CDy*(CDz*HRR(1:LB,34,1)+  &
                        HRR(1:LB,55,1))+  &
                        3.D0*HRR(1:LB,81,1))+  &
                        3.D0*HRR(1:LB,114,1))+  &
                        HRR(1:LB,155,1)
      !=|35,29)
      HRR(1:LB,35,29)=CDz*HRR(1:LB,114,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,81,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,55,1)+  &
                        CDy*(CDz*HRR(1:LB,35,1)+  &
                        HRR(1:LB,56,1))+  &
                        3.D0*HRR(1:LB,83,1))+  &
                        3.D0*HRR(1:LB,117,1))+  &
                        HRR(1:LB,159,1)
      !=|21,30)
      HRR(1:LB,21,30)=CDz*(CDz*HRR(1:LB,57,1)+  &
                        2.D0*HRR(1:LB,93,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,36,1)+  &
                        4.D0*HRR(1:LB,64,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,21,1)+  &
                        2.D0*HRR(1:LB,42,1))+  &
                        HRR(1:LB,70,1))+  &
                        2.D0*HRR(1:LB,100,1))+  &
                        HRR(1:LB,138,1)
      !=|22,30)
      HRR(1:LB,22,30)=CDz*(CDz*HRR(1:LB,58,1)+  &
                        2.D0*HRR(1:LB,94,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,37,1)+  &
                        4.D0*HRR(1:LB,65,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,22,1)+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,71,1))+  &
                        2.D0*HRR(1:LB,101,1))+  &
                        HRR(1:LB,139,1)
      !=|23,30)
      HRR(1:LB,23,30)=CDz*(CDz*HRR(1:LB,59,1)+  &
                        2.D0*HRR(1:LB,95,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,38,1)+  &
                        4.D0*HRR(1:LB,66,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,72,1))+  &
                        2.D0*HRR(1:LB,102,1))+  &
                        HRR(1:LB,140,1)
      !=|24,30)
      HRR(1:LB,24,30)=CDz*(CDz*HRR(1:LB,60,1)+  &
                        2.D0*HRR(1:LB,96,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,39,1)+  &
                        4.D0*HRR(1:LB,67,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,24,1)+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,73,1))+  &
                        2.D0*HRR(1:LB,103,1))+  &
                        HRR(1:LB,141,1)
      !=|25,30)
      HRR(1:LB,25,30)=CDz*(CDz*HRR(1:LB,61,1)+  &
                        2.D0*HRR(1:LB,97,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,40,1)+  &
                        4.D0*HRR(1:LB,68,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,25,1)+  &
                        2.D0*HRR(1:LB,46,1))+  &
                        HRR(1:LB,74,1))+  &
                        2.D0*HRR(1:LB,104,1))+  &
                        HRR(1:LB,142,1)
      !=|26,30)
      HRR(1:LB,26,30)=CDz*(CDz*HRR(1:LB,64,1)+  &
                        2.D0*HRR(1:LB,100,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,42,1)+  &
                        4.D0*HRR(1:LB,70,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,26,1)+  &
                        2.D0*HRR(1:LB,47,1))+  &
                        HRR(1:LB,75,1))+  &
                        2.D0*HRR(1:LB,106,1))+  &
                        HRR(1:LB,145,1)
      !=|27,30)
      HRR(1:LB,27,30)=CDz*(CDz*HRR(1:LB,65,1)+  &
                        2.D0*HRR(1:LB,101,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,43,1)+  &
                        4.D0*HRR(1:LB,71,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,27,1)+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,76,1))+  &
                        2.D0*HRR(1:LB,107,1))+  &
                        HRR(1:LB,146,1)
      !=|28,30)
      HRR(1:LB,28,30)=CDz*(CDz*HRR(1:LB,66,1)+  &
                        2.D0*HRR(1:LB,102,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,44,1)+  &
                        4.D0*HRR(1:LB,72,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,28,1)+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,77,1))+  &
                        2.D0*HRR(1:LB,108,1))+  &
                        HRR(1:LB,147,1)
      !=|29,30)
      HRR(1:LB,29,30)=CDz*(CDz*HRR(1:LB,67,1)+  &
                        2.D0*HRR(1:LB,103,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,45,1)+  &
                        4.D0*HRR(1:LB,73,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,29,1)+  &
                        2.D0*HRR(1:LB,50,1))+  &
                        HRR(1:LB,78,1))+  &
                        2.D0*HRR(1:LB,109,1))+  &
                        HRR(1:LB,148,1)
      !=|30,30)
      HRR(1:LB,30,30)=CDz*(CDz*HRR(1:LB,70,1)+  &
                        2.D0*HRR(1:LB,106,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,47,1)+  &
                        4.D0*HRR(1:LB,75,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,30,1)+  &
                        2.D0*HRR(1:LB,51,1))+  &
                        HRR(1:LB,79,1))+  &
                        2.D0*HRR(1:LB,111,1))+  &
                        HRR(1:LB,151,1)
      !=|31,30)
      HRR(1:LB,31,30)=CDz*(CDz*HRR(1:LB,71,1)+  &
                        2.D0*HRR(1:LB,107,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,48,1)+  &
                        4.D0*HRR(1:LB,76,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,31,1)+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,80,1))+  &
                        2.D0*HRR(1:LB,112,1))+  &
                        HRR(1:LB,152,1)
      !=|32,30)
      HRR(1:LB,32,30)=CDz*(CDz*HRR(1:LB,72,1)+  &
                        2.D0*HRR(1:LB,108,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,49,1)+  &
                        4.D0*HRR(1:LB,77,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,32,1)+  &
                        2.D0*HRR(1:LB,53,1))+  &
                        HRR(1:LB,81,1))+  &
                        2.D0*HRR(1:LB,113,1))+  &
                        HRR(1:LB,153,1)
      !=|33,30)
      HRR(1:LB,33,30)=CDz*(CDz*HRR(1:LB,75,1)+  &
                        2.D0*HRR(1:LB,111,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,51,1)+  &
                        4.D0*HRR(1:LB,79,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,33,1)+  &
                        2.D0*HRR(1:LB,54,1))+  &
                        HRR(1:LB,82,1))+  &
                        2.D0*HRR(1:LB,115,1))+  &
                        HRR(1:LB,156,1)
      !=|34,30)
      HRR(1:LB,34,30)=CDz*(CDz*HRR(1:LB,76,1)+  &
                        2.D0*HRR(1:LB,112,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,52,1)+  &
                        4.D0*HRR(1:LB,80,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,34,1)+  &
                        2.D0*HRR(1:LB,55,1))+  &
                        HRR(1:LB,83,1))+  &
                        2.D0*HRR(1:LB,116,1))+  &
                        HRR(1:LB,157,1)
      !=|35,30)
      HRR(1:LB,35,30)=CDz*(CDz*HRR(1:LB,79,1)+  &
                        2.D0*HRR(1:LB,115,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,54,1)+  &
                        4.D0*HRR(1:LB,82,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,35,1)+  &
                        2.D0*HRR(1:LB,56,1))+  &
                        HRR(1:LB,84,1))+  &
                        2.D0*HRR(1:LB,118,1))+  &
                        HRR(1:LB,160,1)
      !=|21,31)
      HRR(1:LB,21,31)=CDz*(CDz*HRR(1:LB,58,1)+  &
                        2.D0*HRR(1:LB,94,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,36,1)+  &
                        2.D0*HRR(1:LB,64,1))+  &
                        HRR(1:LB,100,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,37,1)+  &
                        2.D0*HRR(1:LB,65,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,21,1)+  &
                        2.D0*HRR(1:LB,42,1))+  &
                        HRR(1:LB,70,1))+  &
                        HRR(1:LB,101,1))+  &
                        HRR(1:LB,139,1)
      !=|22,31)
      HRR(1:LB,22,31)=CDz*(CDz*HRR(1:LB,59,1)+  &
                        2.D0*HRR(1:LB,95,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,37,1)+  &
                        2.D0*HRR(1:LB,65,1))+  &
                        HRR(1:LB,101,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,38,1)+  &
                        2.D0*HRR(1:LB,66,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,22,1)+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,71,1))+  &
                        HRR(1:LB,102,1))+  &
                        HRR(1:LB,140,1)
      !=|23,31)
      HRR(1:LB,23,31)=CDz*(CDz*HRR(1:LB,60,1)+  &
                        2.D0*HRR(1:LB,96,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,38,1)+  &
                        2.D0*HRR(1:LB,66,1))+  &
                        HRR(1:LB,102,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,39,1)+  &
                        2.D0*HRR(1:LB,67,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,72,1))+  &
                        HRR(1:LB,103,1))+  &
                        HRR(1:LB,141,1)
      !=|24,31)
      HRR(1:LB,24,31)=CDz*(CDz*HRR(1:LB,61,1)+  &
                        2.D0*HRR(1:LB,97,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,39,1)+  &
                        2.D0*HRR(1:LB,67,1))+  &
                        HRR(1:LB,103,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,40,1)+  &
                        2.D0*HRR(1:LB,68,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,24,1)+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,73,1))+  &
                        HRR(1:LB,104,1))+  &
                        HRR(1:LB,142,1)
      !=|25,31)
      HRR(1:LB,25,31)=CDz*(CDz*HRR(1:LB,62,1)+  &
                        2.D0*HRR(1:LB,98,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,40,1)+  &
                        2.D0*HRR(1:LB,68,1))+  &
                        HRR(1:LB,104,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,41,1)+  &
                        2.D0*HRR(1:LB,69,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,25,1)+  &
                        2.D0*HRR(1:LB,46,1))+  &
                        HRR(1:LB,74,1))+  &
                        HRR(1:LB,105,1))+  &
                        HRR(1:LB,143,1)
      !=|26,31)
      HRR(1:LB,26,31)=CDz*(CDz*HRR(1:LB,65,1)+  &
                        2.D0*HRR(1:LB,101,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,42,1)+  &
                        2.D0*HRR(1:LB,70,1))+  &
                        HRR(1:LB,106,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,43,1)+  &
                        2.D0*HRR(1:LB,71,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,26,1)+  &
                        2.D0*HRR(1:LB,47,1))+  &
                        HRR(1:LB,75,1))+  &
                        HRR(1:LB,107,1))+  &
                        HRR(1:LB,146,1)
      !=|27,31)
      HRR(1:LB,27,31)=CDz*(CDz*HRR(1:LB,66,1)+  &
                        2.D0*HRR(1:LB,102,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,43,1)+  &
                        2.D0*HRR(1:LB,71,1))+  &
                        HRR(1:LB,107,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,44,1)+  &
                        2.D0*HRR(1:LB,72,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,27,1)+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,76,1))+  &
                        HRR(1:LB,108,1))+  &
                        HRR(1:LB,147,1)
      !=|28,31)
      HRR(1:LB,28,31)=CDz*(CDz*HRR(1:LB,67,1)+  &
                        2.D0*HRR(1:LB,103,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,44,1)+  &
                        2.D0*HRR(1:LB,72,1))+  &
                        HRR(1:LB,108,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,45,1)+  &
                        2.D0*HRR(1:LB,73,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,28,1)+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,77,1))+  &
                        HRR(1:LB,109,1))+  &
                        HRR(1:LB,148,1)
      !=|29,31)
      HRR(1:LB,29,31)=CDz*(CDz*HRR(1:LB,68,1)+  &
                        2.D0*HRR(1:LB,104,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,45,1)+  &
                        2.D0*HRR(1:LB,73,1))+  &
                        HRR(1:LB,109,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,46,1)+  &
                        2.D0*HRR(1:LB,74,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,29,1)+  &
                        2.D0*HRR(1:LB,50,1))+  &
                        HRR(1:LB,78,1))+  &
                        HRR(1:LB,110,1))+  &
                        HRR(1:LB,149,1)
      !=|30,31)
      HRR(1:LB,30,31)=CDz*(CDz*HRR(1:LB,71,1)+  &
                        2.D0*HRR(1:LB,107,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,47,1)+  &
                        2.D0*HRR(1:LB,75,1))+  &
                        HRR(1:LB,111,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,48,1)+  &
                        2.D0*HRR(1:LB,76,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,30,1)+  &
                        2.D0*HRR(1:LB,51,1))+  &
                        HRR(1:LB,79,1))+  &
                        HRR(1:LB,112,1))+  &
                        HRR(1:LB,152,1)
      !=|31,31)
      HRR(1:LB,31,31)=CDz*(CDz*HRR(1:LB,72,1)+  &
                        2.D0*HRR(1:LB,108,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,48,1)+  &
                        2.D0*HRR(1:LB,76,1))+  &
                        HRR(1:LB,112,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,49,1)+  &
                        2.D0*HRR(1:LB,77,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,31,1)+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,80,1))+  &
                        HRR(1:LB,113,1))+  &
                        HRR(1:LB,153,1)
      !=|32,31)
      HRR(1:LB,32,31)=CDz*(CDz*HRR(1:LB,73,1)+  &
                        2.D0*HRR(1:LB,109,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,49,1)+  &
                        2.D0*HRR(1:LB,77,1))+  &
                        HRR(1:LB,113,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,50,1)+  &
                        2.D0*HRR(1:LB,78,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,32,1)+  &
                        2.D0*HRR(1:LB,53,1))+  &
                        HRR(1:LB,81,1))+  &
                        HRR(1:LB,114,1))+  &
                        HRR(1:LB,154,1)
      !=|33,31)
      HRR(1:LB,33,31)=CDz*(CDz*HRR(1:LB,76,1)+  &
                        2.D0*HRR(1:LB,112,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,51,1)+  &
                        2.D0*HRR(1:LB,79,1))+  &
                        HRR(1:LB,115,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,52,1)+  &
                        2.D0*HRR(1:LB,80,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,33,1)+  &
                        2.D0*HRR(1:LB,54,1))+  &
                        HRR(1:LB,82,1))+  &
                        HRR(1:LB,116,1))+  &
                        HRR(1:LB,157,1)
      !=|34,31)
      HRR(1:LB,34,31)=CDz*(CDz*HRR(1:LB,77,1)+  &
                        2.D0*HRR(1:LB,113,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,52,1)+  &
                        2.D0*HRR(1:LB,80,1))+  &
                        HRR(1:LB,116,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,53,1)+  &
                        2.D0*HRR(1:LB,81,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,34,1)+  &
                        2.D0*HRR(1:LB,55,1))+  &
                        HRR(1:LB,83,1))+  &
                        HRR(1:LB,117,1))+  &
                        HRR(1:LB,158,1)
      !=|35,31)
      HRR(1:LB,35,31)=CDz*(CDz*HRR(1:LB,80,1)+  &
                        2.D0*HRR(1:LB,116,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,54,1)+  &
                        2.D0*HRR(1:LB,82,1))+  &
                        HRR(1:LB,118,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,55,1)+  &
                        2.D0*HRR(1:LB,83,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,35,1)+  &
                        2.D0*HRR(1:LB,56,1))+  &
                        HRR(1:LB,84,1))+  &
                        HRR(1:LB,119,1))+  &
                        HRR(1:LB,161,1)
      !=|21,32)
      HRR(1:LB,21,32)=CDz*(CDz*HRR(1:LB,59,1)+  &
                        2.D0*HRR(1:LB,95,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,37,1)+  &
                        4.D0*HRR(1:LB,65,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,21,1)+  &
                        2.D0*HRR(1:LB,42,1))+  &
                        HRR(1:LB,70,1))+  &
                        2.D0*HRR(1:LB,101,1))+  &
                        HRR(1:LB,140,1)
      !=|22,32)
      HRR(1:LB,22,32)=CDz*(CDz*HRR(1:LB,60,1)+  &
                        2.D0*HRR(1:LB,96,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,38,1)+  &
                        4.D0*HRR(1:LB,66,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,22,1)+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,71,1))+  &
                        2.D0*HRR(1:LB,102,1))+  &
                        HRR(1:LB,141,1)
      !=|23,32)
      HRR(1:LB,23,32)=CDz*(CDz*HRR(1:LB,61,1)+  &
                        2.D0*HRR(1:LB,97,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,39,1)+  &
                        4.D0*HRR(1:LB,67,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,72,1))+  &
                        2.D0*HRR(1:LB,103,1))+  &
                        HRR(1:LB,142,1)
      !=|24,32)
      HRR(1:LB,24,32)=CDz*(CDz*HRR(1:LB,62,1)+  &
                        2.D0*HRR(1:LB,98,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,40,1)+  &
                        4.D0*HRR(1:LB,68,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,24,1)+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,73,1))+  &
                        2.D0*HRR(1:LB,104,1))+  &
                        HRR(1:LB,143,1)
      !=|25,32)
      HRR(1:LB,25,32)=CDz*(CDz*HRR(1:LB,63,1)+  &
                        2.D0*HRR(1:LB,99,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,41,1)+  &
                        4.D0*HRR(1:LB,69,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,25,1)+  &
                        2.D0*HRR(1:LB,46,1))+  &
                        HRR(1:LB,74,1))+  &
                        2.D0*HRR(1:LB,105,1))+  &
                        HRR(1:LB,144,1)
      !=|26,32)
      HRR(1:LB,26,32)=CDz*(CDz*HRR(1:LB,66,1)+  &
                        2.D0*HRR(1:LB,102,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,43,1)+  &
                        4.D0*HRR(1:LB,71,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,26,1)+  &
                        2.D0*HRR(1:LB,47,1))+  &
                        HRR(1:LB,75,1))+  &
                        2.D0*HRR(1:LB,107,1))+  &
                        HRR(1:LB,147,1)
      !=|27,32)
      HRR(1:LB,27,32)=CDz*(CDz*HRR(1:LB,67,1)+  &
                        2.D0*HRR(1:LB,103,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,44,1)+  &
                        4.D0*HRR(1:LB,72,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,27,1)+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,76,1))+  &
                        2.D0*HRR(1:LB,108,1))+  &
                        HRR(1:LB,148,1)
      !=|28,32)
      HRR(1:LB,28,32)=CDz*(CDz*HRR(1:LB,68,1)+  &
                        2.D0*HRR(1:LB,104,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,45,1)+  &
                        4.D0*HRR(1:LB,73,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,28,1)+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,77,1))+  &
                        2.D0*HRR(1:LB,109,1))+  &
                        HRR(1:LB,149,1)
      !=|29,32)
      HRR(1:LB,29,32)=CDz*(CDz*HRR(1:LB,69,1)+  &
                        2.D0*HRR(1:LB,105,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,46,1)+  &
                        4.D0*HRR(1:LB,74,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,29,1)+  &
                        2.D0*HRR(1:LB,50,1))+  &
                        HRR(1:LB,78,1))+  &
                        2.D0*HRR(1:LB,110,1))+  &
                        HRR(1:LB,150,1)
      !=|30,32)
      HRR(1:LB,30,32)=CDz*(CDz*HRR(1:LB,72,1)+  &
                        2.D0*HRR(1:LB,108,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,48,1)+  &
                        4.D0*HRR(1:LB,76,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,30,1)+  &
                        2.D0*HRR(1:LB,51,1))+  &
                        HRR(1:LB,79,1))+  &
                        2.D0*HRR(1:LB,112,1))+  &
                        HRR(1:LB,153,1)
      !=|31,32)
      HRR(1:LB,31,32)=CDz*(CDz*HRR(1:LB,73,1)+  &
                        2.D0*HRR(1:LB,109,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,49,1)+  &
                        4.D0*HRR(1:LB,77,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,31,1)+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,80,1))+  &
                        2.D0*HRR(1:LB,113,1))+  &
                        HRR(1:LB,154,1)
      !=|32,32)
      HRR(1:LB,32,32)=CDz*(CDz*HRR(1:LB,74,1)+  &
                        2.D0*HRR(1:LB,110,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,50,1)+  &
                        4.D0*HRR(1:LB,78,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,32,1)+  &
                        2.D0*HRR(1:LB,53,1))+  &
                        HRR(1:LB,81,1))+  &
                        2.D0*HRR(1:LB,114,1))+  &
                        HRR(1:LB,155,1)
      !=|33,32)
      HRR(1:LB,33,32)=CDz*(CDz*HRR(1:LB,77,1)+  &
                        2.D0*HRR(1:LB,113,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,52,1)+  &
                        4.D0*HRR(1:LB,80,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,33,1)+  &
                        2.D0*HRR(1:LB,54,1))+  &
                        HRR(1:LB,82,1))+  &
                        2.D0*HRR(1:LB,116,1))+  &
                        HRR(1:LB,158,1)
      !=|34,32)
      HRR(1:LB,34,32)=CDz*(CDz*HRR(1:LB,78,1)+  &
                        2.D0*HRR(1:LB,114,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,53,1)+  &
                        4.D0*HRR(1:LB,81,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,34,1)+  &
                        2.D0*HRR(1:LB,55,1))+  &
                        HRR(1:LB,83,1))+  &
                        2.D0*HRR(1:LB,117,1))+  &
                        HRR(1:LB,159,1)
      !=|35,32)
      HRR(1:LB,35,32)=CDz*(CDz*HRR(1:LB,81,1)+  &
                        2.D0*HRR(1:LB,117,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,55,1)+  &
                        4.D0*HRR(1:LB,83,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,35,1)+  &
                        2.D0*HRR(1:LB,56,1))+  &
                        HRR(1:LB,84,1))+  &
                        2.D0*HRR(1:LB,119,1))+  &
                        HRR(1:LB,162,1)
      !=|21,33)
      HRR(1:LB,21,33)=CDz*(CDz*(CDz*HRR(1:LB,36,1)+  &
                        3.D0*HRR(1:LB,64,1))+  &
                        3.D0*HRR(1:LB,100,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,21,1)+  &
                        3.D0*HRR(1:LB,42,1))+  &
                        3.D0*HRR(1:LB,70,1))+  &
                        HRR(1:LB,106,1))+  &
                        HRR(1:LB,145,1)
      !=|22,33)
      HRR(1:LB,22,33)=CDz*(CDz*(CDz*HRR(1:LB,37,1)+  &
                        3.D0*HRR(1:LB,65,1))+  &
                        3.D0*HRR(1:LB,101,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,22,1)+  &
                        3.D0*HRR(1:LB,43,1))+  &
                        3.D0*HRR(1:LB,71,1))+  &
                        HRR(1:LB,107,1))+  &
                        HRR(1:LB,146,1)
      !=|23,33)
      HRR(1:LB,23,33)=CDz*(CDz*(CDz*HRR(1:LB,38,1)+  &
                        3.D0*HRR(1:LB,66,1))+  &
                        3.D0*HRR(1:LB,102,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,23,1)+  &
                        3.D0*HRR(1:LB,44,1))+  &
                        3.D0*HRR(1:LB,72,1))+  &
                        HRR(1:LB,108,1))+  &
                        HRR(1:LB,147,1)
      !=|24,33)
      HRR(1:LB,24,33)=CDz*(CDz*(CDz*HRR(1:LB,39,1)+  &
                        3.D0*HRR(1:LB,67,1))+  &
                        3.D0*HRR(1:LB,103,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,24,1)+  &
                        3.D0*HRR(1:LB,45,1))+  &
                        3.D0*HRR(1:LB,73,1))+  &
                        HRR(1:LB,109,1))+  &
                        HRR(1:LB,148,1)
      !=|25,33)
      HRR(1:LB,25,33)=CDz*(CDz*(CDz*HRR(1:LB,40,1)+  &
                        3.D0*HRR(1:LB,68,1))+  &
                        3.D0*HRR(1:LB,104,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,25,1)+  &
                        3.D0*HRR(1:LB,46,1))+  &
                        3.D0*HRR(1:LB,74,1))+  &
                        HRR(1:LB,110,1))+  &
                        HRR(1:LB,149,1)
      !=|26,33)
      HRR(1:LB,26,33)=CDz*(CDz*(CDz*HRR(1:LB,42,1)+  &
                        3.D0*HRR(1:LB,70,1))+  &
                        3.D0*HRR(1:LB,106,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,26,1)+  &
                        3.D0*HRR(1:LB,47,1))+  &
                        3.D0*HRR(1:LB,75,1))+  &
                        HRR(1:LB,111,1))+  &
                        HRR(1:LB,151,1)
      !=|27,33)
      HRR(1:LB,27,33)=CDz*(CDz*(CDz*HRR(1:LB,43,1)+  &
                        3.D0*HRR(1:LB,71,1))+  &
                        3.D0*HRR(1:LB,107,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,27,1)+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        3.D0*HRR(1:LB,76,1))+  &
                        HRR(1:LB,112,1))+  &
                        HRR(1:LB,152,1)
      !=|28,33)
      HRR(1:LB,28,33)=CDz*(CDz*(CDz*HRR(1:LB,44,1)+  &
                        3.D0*HRR(1:LB,72,1))+  &
                        3.D0*HRR(1:LB,108,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,28,1)+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        3.D0*HRR(1:LB,77,1))+  &
                        HRR(1:LB,113,1))+  &
                        HRR(1:LB,153,1)
      !=|29,33)
      HRR(1:LB,29,33)=CDz*(CDz*(CDz*HRR(1:LB,45,1)+  &
                        3.D0*HRR(1:LB,73,1))+  &
                        3.D0*HRR(1:LB,109,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,29,1)+  &
                        3.D0*HRR(1:LB,50,1))+  &
                        3.D0*HRR(1:LB,78,1))+  &
                        HRR(1:LB,114,1))+  &
                        HRR(1:LB,154,1)
      !=|30,33)
      HRR(1:LB,30,33)=CDz*(CDz*(CDz*HRR(1:LB,47,1)+  &
                        3.D0*HRR(1:LB,75,1))+  &
                        3.D0*HRR(1:LB,111,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,30,1)+  &
                        3.D0*HRR(1:LB,51,1))+  &
                        3.D0*HRR(1:LB,79,1))+  &
                        HRR(1:LB,115,1))+  &
                        HRR(1:LB,156,1)
      !=|31,33)
      HRR(1:LB,31,33)=CDz*(CDz*(CDz*HRR(1:LB,48,1)+  &
                        3.D0*HRR(1:LB,76,1))+  &
                        3.D0*HRR(1:LB,112,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,31,1)+  &
                        3.D0*HRR(1:LB,52,1))+  &
                        3.D0*HRR(1:LB,80,1))+  &
                        HRR(1:LB,116,1))+  &
                        HRR(1:LB,157,1)
      !=|32,33)
      HRR(1:LB,32,33)=CDz*(CDz*(CDz*HRR(1:LB,49,1)+  &
                        3.D0*HRR(1:LB,77,1))+  &
                        3.D0*HRR(1:LB,113,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,32,1)+  &
                        3.D0*HRR(1:LB,53,1))+  &
                        3.D0*HRR(1:LB,81,1))+  &
                        HRR(1:LB,117,1))+  &
                        HRR(1:LB,158,1)
      !=|33,33)
      HRR(1:LB,33,33)=CDz*(CDz*(CDz*HRR(1:LB,51,1)+  &
                        3.D0*HRR(1:LB,79,1))+  &
                        3.D0*HRR(1:LB,115,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,33,1)+  &
                        3.D0*HRR(1:LB,54,1))+  &
                        3.D0*HRR(1:LB,82,1))+  &
                        HRR(1:LB,118,1))+  &
                        HRR(1:LB,160,1)
      !=|34,33)
      HRR(1:LB,34,33)=CDz*(CDz*(CDz*HRR(1:LB,52,1)+  &
                        3.D0*HRR(1:LB,80,1))+  &
                        3.D0*HRR(1:LB,116,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,34,1)+  &
                        3.D0*HRR(1:LB,55,1))+  &
                        3.D0*HRR(1:LB,83,1))+  &
                        HRR(1:LB,119,1))+  &
                        HRR(1:LB,161,1)
      !=|35,33)
      HRR(1:LB,35,33)=CDz*(CDz*(CDz*HRR(1:LB,54,1)+  &
                        3.D0*HRR(1:LB,82,1))+  &
                        3.D0*HRR(1:LB,118,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,35,1)+  &
                        3.D0*HRR(1:LB,56,1))+  &
                        3.D0*HRR(1:LB,84,1))+  &
                        HRR(1:LB,120,1))+  &
                        HRR(1:LB,163,1)
      !=|21,34)
      HRR(1:LB,21,34)=CDz*(CDz*(CDz*HRR(1:LB,37,1)+  &
                        3.D0*HRR(1:LB,65,1))+  &
                        3.D0*HRR(1:LB,101,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,21,1)+  &
                        3.D0*HRR(1:LB,42,1))+  &
                        3.D0*HRR(1:LB,70,1))+  &
                        HRR(1:LB,106,1))+  &
                        HRR(1:LB,146,1)
      !=|22,34)
      HRR(1:LB,22,34)=CDz*(CDz*(CDz*HRR(1:LB,38,1)+  &
                        3.D0*HRR(1:LB,66,1))+  &
                        3.D0*HRR(1:LB,102,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,22,1)+  &
                        3.D0*HRR(1:LB,43,1))+  &
                        3.D0*HRR(1:LB,71,1))+  &
                        HRR(1:LB,107,1))+  &
                        HRR(1:LB,147,1)
      !=|23,34)
      HRR(1:LB,23,34)=CDz*(CDz*(CDz*HRR(1:LB,39,1)+  &
                        3.D0*HRR(1:LB,67,1))+  &
                        3.D0*HRR(1:LB,103,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,23,1)+  &
                        3.D0*HRR(1:LB,44,1))+  &
                        3.D0*HRR(1:LB,72,1))+  &
                        HRR(1:LB,108,1))+  &
                        HRR(1:LB,148,1)
      !=|24,34)
      HRR(1:LB,24,34)=CDz*(CDz*(CDz*HRR(1:LB,40,1)+  &
                        3.D0*HRR(1:LB,68,1))+  &
                        3.D0*HRR(1:LB,104,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,24,1)+  &
                        3.D0*HRR(1:LB,45,1))+  &
                        3.D0*HRR(1:LB,73,1))+  &
                        HRR(1:LB,109,1))+  &
                        HRR(1:LB,149,1)
      !=|25,34)
      HRR(1:LB,25,34)=CDz*(CDz*(CDz*HRR(1:LB,41,1)+  &
                        3.D0*HRR(1:LB,69,1))+  &
                        3.D0*HRR(1:LB,105,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,25,1)+  &
                        3.D0*HRR(1:LB,46,1))+  &
                        3.D0*HRR(1:LB,74,1))+  &
                        HRR(1:LB,110,1))+  &
                        HRR(1:LB,150,1)
      !=|26,34)
      HRR(1:LB,26,34)=CDz*(CDz*(CDz*HRR(1:LB,43,1)+  &
                        3.D0*HRR(1:LB,71,1))+  &
                        3.D0*HRR(1:LB,107,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,26,1)+  &
                        3.D0*HRR(1:LB,47,1))+  &
                        3.D0*HRR(1:LB,75,1))+  &
                        HRR(1:LB,111,1))+  &
                        HRR(1:LB,152,1)
      !=|27,34)
      HRR(1:LB,27,34)=CDz*(CDz*(CDz*HRR(1:LB,44,1)+  &
                        3.D0*HRR(1:LB,72,1))+  &
                        3.D0*HRR(1:LB,108,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,27,1)+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        3.D0*HRR(1:LB,76,1))+  &
                        HRR(1:LB,112,1))+  &
                        HRR(1:LB,153,1)
      !=|28,34)
      HRR(1:LB,28,34)=CDz*(CDz*(CDz*HRR(1:LB,45,1)+  &
                        3.D0*HRR(1:LB,73,1))+  &
                        3.D0*HRR(1:LB,109,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,28,1)+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        3.D0*HRR(1:LB,77,1))+  &
                        HRR(1:LB,113,1))+  &
                        HRR(1:LB,154,1)
      !=|29,34)
      HRR(1:LB,29,34)=CDz*(CDz*(CDz*HRR(1:LB,46,1)+  &
                        3.D0*HRR(1:LB,74,1))+  &
                        3.D0*HRR(1:LB,110,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,29,1)+  &
                        3.D0*HRR(1:LB,50,1))+  &
                        3.D0*HRR(1:LB,78,1))+  &
                        HRR(1:LB,114,1))+  &
                        HRR(1:LB,155,1)
      !=|30,34)
      HRR(1:LB,30,34)=CDz*(CDz*(CDz*HRR(1:LB,48,1)+  &
                        3.D0*HRR(1:LB,76,1))+  &
                        3.D0*HRR(1:LB,112,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,30,1)+  &
                        3.D0*HRR(1:LB,51,1))+  &
                        3.D0*HRR(1:LB,79,1))+  &
                        HRR(1:LB,115,1))+  &
                        HRR(1:LB,157,1)
      !=|31,34)
      HRR(1:LB,31,34)=CDz*(CDz*(CDz*HRR(1:LB,49,1)+  &
                        3.D0*HRR(1:LB,77,1))+  &
                        3.D0*HRR(1:LB,113,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,31,1)+  &
                        3.D0*HRR(1:LB,52,1))+  &
                        3.D0*HRR(1:LB,80,1))+  &
                        HRR(1:LB,116,1))+  &
                        HRR(1:LB,158,1)
      !=|32,34)
      HRR(1:LB,32,34)=CDz*(CDz*(CDz*HRR(1:LB,50,1)+  &
                        3.D0*HRR(1:LB,78,1))+  &
                        3.D0*HRR(1:LB,114,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,32,1)+  &
                        3.D0*HRR(1:LB,53,1))+  &
                        3.D0*HRR(1:LB,81,1))+  &
                        HRR(1:LB,117,1))+  &
                        HRR(1:LB,159,1)
      !=|33,34)
      HRR(1:LB,33,34)=CDz*(CDz*(CDz*HRR(1:LB,52,1)+  &
                        3.D0*HRR(1:LB,80,1))+  &
                        3.D0*HRR(1:LB,116,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,33,1)+  &
                        3.D0*HRR(1:LB,54,1))+  &
                        3.D0*HRR(1:LB,82,1))+  &
                        HRR(1:LB,118,1))+  &
                        HRR(1:LB,161,1)
      !=|34,34)
      HRR(1:LB,34,34)=CDz*(CDz*(CDz*HRR(1:LB,53,1)+  &
                        3.D0*HRR(1:LB,81,1))+  &
                        3.D0*HRR(1:LB,117,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,34,1)+  &
                        3.D0*HRR(1:LB,55,1))+  &
                        3.D0*HRR(1:LB,83,1))+  &
                        HRR(1:LB,119,1))+  &
                        HRR(1:LB,162,1)
      !=|35,34)
      HRR(1:LB,35,34)=CDz*(CDz*(CDz*HRR(1:LB,55,1)+  &
                        3.D0*HRR(1:LB,83,1))+  &
                        3.D0*HRR(1:LB,119,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,35,1)+  &
                        3.D0*HRR(1:LB,56,1))+  &
                        3.D0*HRR(1:LB,84,1))+  &
                        HRR(1:LB,120,1))+  &
                        HRR(1:LB,164,1)
      !=|21,35)
      HRR(1:LB,21,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,21,1)+  &
                        4.D0*HRR(1:LB,42,1))+  &
                        6.D0*HRR(1:LB,70,1))+  &
                        4.D0*HRR(1:LB,106,1))+  &
                        HRR(1:LB,151,1)
      !=|22,35)
      HRR(1:LB,22,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,22,1)+  &
                        4.D0*HRR(1:LB,43,1))+  &
                        6.D0*HRR(1:LB,71,1))+  &
                        4.D0*HRR(1:LB,107,1))+  &
                        HRR(1:LB,152,1)
      !=|23,35)
      HRR(1:LB,23,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,23,1)+  &
                        4.D0*HRR(1:LB,44,1))+  &
                        6.D0*HRR(1:LB,72,1))+  &
                        4.D0*HRR(1:LB,108,1))+  &
                        HRR(1:LB,153,1)
      !=|24,35)
      HRR(1:LB,24,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,24,1)+  &
                        4.D0*HRR(1:LB,45,1))+  &
                        6.D0*HRR(1:LB,73,1))+  &
                        4.D0*HRR(1:LB,109,1))+  &
                        HRR(1:LB,154,1)
      !=|25,35)
      HRR(1:LB,25,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,25,1)+  &
                        4.D0*HRR(1:LB,46,1))+  &
                        6.D0*HRR(1:LB,74,1))+  &
                        4.D0*HRR(1:LB,110,1))+  &
                        HRR(1:LB,155,1)
      !=|26,35)
      HRR(1:LB,26,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,26,1)+  &
                        4.D0*HRR(1:LB,47,1))+  &
                        6.D0*HRR(1:LB,75,1))+  &
                        4.D0*HRR(1:LB,111,1))+  &
                        HRR(1:LB,156,1)
      !=|27,35)
      HRR(1:LB,27,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,27,1)+  &
                        4.D0*HRR(1:LB,48,1))+  &
                        6.D0*HRR(1:LB,76,1))+  &
                        4.D0*HRR(1:LB,112,1))+  &
                        HRR(1:LB,157,1)
      !=|28,35)
      HRR(1:LB,28,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,28,1)+  &
                        4.D0*HRR(1:LB,49,1))+  &
                        6.D0*HRR(1:LB,77,1))+  &
                        4.D0*HRR(1:LB,113,1))+  &
                        HRR(1:LB,158,1)
      !=|29,35)
      HRR(1:LB,29,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,29,1)+  &
                        4.D0*HRR(1:LB,50,1))+  &
                        6.D0*HRR(1:LB,78,1))+  &
                        4.D0*HRR(1:LB,114,1))+  &
                        HRR(1:LB,159,1)
      !=|30,35)
      HRR(1:LB,30,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,30,1)+  &
                        4.D0*HRR(1:LB,51,1))+  &
                        6.D0*HRR(1:LB,79,1))+  &
                        4.D0*HRR(1:LB,115,1))+  &
                        HRR(1:LB,160,1)
      !=|31,35)
      HRR(1:LB,31,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,31,1)+  &
                        4.D0*HRR(1:LB,52,1))+  &
                        6.D0*HRR(1:LB,80,1))+  &
                        4.D0*HRR(1:LB,116,1))+  &
                        HRR(1:LB,161,1)
      !=|32,35)
      HRR(1:LB,32,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,32,1)+  &
                        4.D0*HRR(1:LB,53,1))+  &
                        6.D0*HRR(1:LB,81,1))+  &
                        4.D0*HRR(1:LB,117,1))+  &
                        HRR(1:LB,162,1)
      !=|33,35)
      HRR(1:LB,33,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,33,1)+  &
                        4.D0*HRR(1:LB,54,1))+  &
                        6.D0*HRR(1:LB,82,1))+  &
                        4.D0*HRR(1:LB,118,1))+  &
                        HRR(1:LB,163,1)
      !=|34,35)
      HRR(1:LB,34,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,34,1)+  &
                        4.D0*HRR(1:LB,55,1))+  &
                        6.D0*HRR(1:LB,83,1))+  &
                        4.D0*HRR(1:LB,119,1))+  &
                        HRR(1:LB,164,1)
      !=|35,35)
      HRR(1:LB,35,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,35,1)+  &
                        4.D0*HRR(1:LB,56,1))+  &
                        6.D0*HRR(1:LB,84,1))+  &
                        4.D0*HRR(1:LB,120,1))+  &
                        HRR(1:LB,165,1)
      !=|11,21)
      HRR(1:LB,11,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,11,1)+  &
                        4.D0*HRR(1:LB,21,1))+  &
                        6.D0*HRR(1:LB,36,1))+  &
                        4.D0*HRR(1:LB,57,1))+  &
                        HRR(1:LB,85,1)
      !=|12,21)
      HRR(1:LB,12,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,12,1)+  &
                        4.D0*HRR(1:LB,22,1))+  &
                        6.D0*HRR(1:LB,37,1))+  &
                        4.D0*HRR(1:LB,58,1))+  &
                        HRR(1:LB,86,1)
      !=|13,21)
      HRR(1:LB,13,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,13,1)+  &
                        4.D0*HRR(1:LB,23,1))+  &
                        6.D0*HRR(1:LB,38,1))+  &
                        4.D0*HRR(1:LB,59,1))+  &
                        HRR(1:LB,87,1)
      !=|14,21)
      HRR(1:LB,14,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,14,1)+  &
                        4.D0*HRR(1:LB,24,1))+  &
                        6.D0*HRR(1:LB,39,1))+  &
                        4.D0*HRR(1:LB,60,1))+  &
                        HRR(1:LB,88,1)
      !=|15,21)
      HRR(1:LB,15,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,15,1)+  &
                        4.D0*HRR(1:LB,26,1))+  &
                        6.D0*HRR(1:LB,42,1))+  &
                        4.D0*HRR(1:LB,64,1))+  &
                        HRR(1:LB,93,1)
      !=|16,21)
      HRR(1:LB,16,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,16,1)+  &
                        4.D0*HRR(1:LB,27,1))+  &
                        6.D0*HRR(1:LB,43,1))+  &
                        4.D0*HRR(1:LB,65,1))+  &
                        HRR(1:LB,94,1)
      !=|17,21)
      HRR(1:LB,17,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,17,1)+  &
                        4.D0*HRR(1:LB,28,1))+  &
                        6.D0*HRR(1:LB,44,1))+  &
                        4.D0*HRR(1:LB,66,1))+  &
                        HRR(1:LB,95,1)
      !=|18,21)
      HRR(1:LB,18,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,18,1)+  &
                        4.D0*HRR(1:LB,30,1))+  &
                        6.D0*HRR(1:LB,47,1))+  &
                        4.D0*HRR(1:LB,70,1))+  &
                        HRR(1:LB,100,1)
      !=|19,21)
      HRR(1:LB,19,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,19,1)+  &
                        4.D0*HRR(1:LB,31,1))+  &
                        6.D0*HRR(1:LB,48,1))+  &
                        4.D0*HRR(1:LB,71,1))+  &
                        HRR(1:LB,101,1)
      !=|20,21)
      HRR(1:LB,20,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,20,1)+  &
                        4.D0*HRR(1:LB,33,1))+  &
                        6.D0*HRR(1:LB,51,1))+  &
                        4.D0*HRR(1:LB,75,1))+  &
                        HRR(1:LB,106,1)
      !=|11,22)
      HRR(1:LB,11,22)=CDy*HRR(1:LB,57,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,36,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,21,1)+  &
                        CDx*(CDy*HRR(1:LB,11,1)+  &
                        HRR(1:LB,22,1))+  &
                        3.D0*HRR(1:LB,37,1))+  &
                        3.D0*HRR(1:LB,58,1))+  &
                        HRR(1:LB,86,1)
      !=|12,22)
      HRR(1:LB,12,22)=CDy*HRR(1:LB,58,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,37,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,22,1)+  &
                        CDx*(CDy*HRR(1:LB,12,1)+  &
                        HRR(1:LB,23,1))+  &
                        3.D0*HRR(1:LB,38,1))+  &
                        3.D0*HRR(1:LB,59,1))+  &
                        HRR(1:LB,87,1)
      !=|13,22)
      HRR(1:LB,13,22)=CDy*HRR(1:LB,59,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,38,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,23,1)+  &
                        CDx*(CDy*HRR(1:LB,13,1)+  &
                        HRR(1:LB,24,1))+  &
                        3.D0*HRR(1:LB,39,1))+  &
                        3.D0*HRR(1:LB,60,1))+  &
                        HRR(1:LB,88,1)
      !=|14,22)
      HRR(1:LB,14,22)=CDy*HRR(1:LB,60,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,39,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,24,1)+  &
                        CDx*(CDy*HRR(1:LB,14,1)+  &
                        HRR(1:LB,25,1))+  &
                        3.D0*HRR(1:LB,40,1))+  &
                        3.D0*HRR(1:LB,61,1))+  &
                        HRR(1:LB,89,1)
      !=|15,22)
      HRR(1:LB,15,22)=CDy*HRR(1:LB,64,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,42,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,26,1)+  &
                        CDx*(CDy*HRR(1:LB,15,1)+  &
                        HRR(1:LB,27,1))+  &
                        3.D0*HRR(1:LB,43,1))+  &
                        3.D0*HRR(1:LB,65,1))+  &
                        HRR(1:LB,94,1)
      !=|16,22)
      HRR(1:LB,16,22)=CDy*HRR(1:LB,65,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,43,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,27,1)+  &
                        CDx*(CDy*HRR(1:LB,16,1)+  &
                        HRR(1:LB,28,1))+  &
                        3.D0*HRR(1:LB,44,1))+  &
                        3.D0*HRR(1:LB,66,1))+  &
                        HRR(1:LB,95,1)
      !=|17,22)
      HRR(1:LB,17,22)=CDy*HRR(1:LB,66,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,44,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,28,1)+  &
                        CDx*(CDy*HRR(1:LB,17,1)+  &
                        HRR(1:LB,29,1))+  &
                        3.D0*HRR(1:LB,45,1))+  &
                        3.D0*HRR(1:LB,67,1))+  &
                        HRR(1:LB,96,1)
      !=|18,22)
      HRR(1:LB,18,22)=CDy*HRR(1:LB,70,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,47,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,30,1)+  &
                        CDx*(CDy*HRR(1:LB,18,1)+  &
                        HRR(1:LB,31,1))+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        3.D0*HRR(1:LB,71,1))+  &
                        HRR(1:LB,101,1)
      !=|19,22)
      HRR(1:LB,19,22)=CDy*HRR(1:LB,71,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,48,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,31,1)+  &
                        CDx*(CDy*HRR(1:LB,19,1)+  &
                        HRR(1:LB,32,1))+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        3.D0*HRR(1:LB,72,1))+  &
                        HRR(1:LB,102,1)
      !=|20,22)
      HRR(1:LB,20,22)=CDy*HRR(1:LB,75,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,51,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,33,1)+  &
                        CDx*(CDy*HRR(1:LB,20,1)+  &
                        HRR(1:LB,34,1))+  &
                        3.D0*HRR(1:LB,52,1))+  &
                        3.D0*HRR(1:LB,76,1))+  &
                        HRR(1:LB,107,1)
      !=|11,23)
      HRR(1:LB,11,23)=CDy*(CDy*HRR(1:LB,36,1)+  &
                        2.D0*HRR(1:LB,58,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,21,1)+  &
                        4.D0*HRR(1:LB,37,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,11,1)+  &
                        2.D0*HRR(1:LB,22,1))+  &
                        HRR(1:LB,38,1))+  &
                        2.D0*HRR(1:LB,59,1))+  &
                        HRR(1:LB,87,1)
      !=|12,23)
      HRR(1:LB,12,23)=CDy*(CDy*HRR(1:LB,37,1)+  &
                        2.D0*HRR(1:LB,59,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,22,1)+  &
                        4.D0*HRR(1:LB,38,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,12,1)+  &
                        2.D0*HRR(1:LB,23,1))+  &
                        HRR(1:LB,39,1))+  &
                        2.D0*HRR(1:LB,60,1))+  &
                        HRR(1:LB,88,1)
      !=|13,23)
      HRR(1:LB,13,23)=CDy*(CDy*HRR(1:LB,38,1)+  &
                        2.D0*HRR(1:LB,60,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,23,1)+  &
                        4.D0*HRR(1:LB,39,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,13,1)+  &
                        2.D0*HRR(1:LB,24,1))+  &
                        HRR(1:LB,40,1))+  &
                        2.D0*HRR(1:LB,61,1))+  &
                        HRR(1:LB,89,1)
      !=|14,23)
      HRR(1:LB,14,23)=CDy*(CDy*HRR(1:LB,39,1)+  &
                        2.D0*HRR(1:LB,61,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,24,1)+  &
                        4.D0*HRR(1:LB,40,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,14,1)+  &
                        2.D0*HRR(1:LB,25,1))+  &
                        HRR(1:LB,41,1))+  &
                        2.D0*HRR(1:LB,62,1))+  &
                        HRR(1:LB,90,1)
      !=|15,23)
      HRR(1:LB,15,23)=CDy*(CDy*HRR(1:LB,42,1)+  &
                        2.D0*HRR(1:LB,65,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,26,1)+  &
                        4.D0*HRR(1:LB,43,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,15,1)+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        HRR(1:LB,44,1))+  &
                        2.D0*HRR(1:LB,66,1))+  &
                        HRR(1:LB,95,1)
      !=|16,23)
      HRR(1:LB,16,23)=CDy*(CDy*HRR(1:LB,43,1)+  &
                        2.D0*HRR(1:LB,66,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,27,1)+  &
                        4.D0*HRR(1:LB,44,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,16,1)+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,45,1))+  &
                        2.D0*HRR(1:LB,67,1))+  &
                        HRR(1:LB,96,1)
      !=|17,23)
      HRR(1:LB,17,23)=CDy*(CDy*HRR(1:LB,44,1)+  &
                        2.D0*HRR(1:LB,67,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,28,1)+  &
                        4.D0*HRR(1:LB,45,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,17,1)+  &
                        2.D0*HRR(1:LB,29,1))+  &
                        HRR(1:LB,46,1))+  &
                        2.D0*HRR(1:LB,68,1))+  &
                        HRR(1:LB,97,1)
      !=|18,23)
      HRR(1:LB,18,23)=CDy*(CDy*HRR(1:LB,47,1)+  &
                        2.D0*HRR(1:LB,71,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,30,1)+  &
                        4.D0*HRR(1:LB,48,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,18,1)+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,49,1))+  &
                        2.D0*HRR(1:LB,72,1))+  &
                        HRR(1:LB,102,1)
      !=|19,23)
      HRR(1:LB,19,23)=CDy*(CDy*HRR(1:LB,48,1)+  &
                        2.D0*HRR(1:LB,72,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,31,1)+  &
                        4.D0*HRR(1:LB,49,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,19,1)+  &
                        2.D0*HRR(1:LB,32,1))+  &
                        HRR(1:LB,50,1))+  &
                        2.D0*HRR(1:LB,73,1))+  &
                        HRR(1:LB,103,1)
      !=|20,23)
      HRR(1:LB,20,23)=CDy*(CDy*HRR(1:LB,51,1)+  &
                        2.D0*HRR(1:LB,76,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,33,1)+  &
                        4.D0*HRR(1:LB,52,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,20,1)+  &
                        2.D0*HRR(1:LB,34,1))+  &
                        HRR(1:LB,53,1))+  &
                        2.D0*HRR(1:LB,77,1))+  &
                        HRR(1:LB,108,1)
      !=|11,24)
      HRR(1:LB,11,24)=CDy*(CDy*(CDy*HRR(1:LB,21,1)+  &
                        3.D0*HRR(1:LB,37,1))+  &
                        3.D0*HRR(1:LB,59,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,11,1)+  &
                        3.D0*HRR(1:LB,22,1))+  &
                        3.D0*HRR(1:LB,38,1))+  &
                        HRR(1:LB,60,1))+  &
                        HRR(1:LB,88,1)
      !=|12,24)
      HRR(1:LB,12,24)=CDy*(CDy*(CDy*HRR(1:LB,22,1)+  &
                        3.D0*HRR(1:LB,38,1))+  &
                        3.D0*HRR(1:LB,60,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,12,1)+  &
                        3.D0*HRR(1:LB,23,1))+  &
                        3.D0*HRR(1:LB,39,1))+  &
                        HRR(1:LB,61,1))+  &
                        HRR(1:LB,89,1)
      !=|13,24)
      HRR(1:LB,13,24)=CDy*(CDy*(CDy*HRR(1:LB,23,1)+  &
                        3.D0*HRR(1:LB,39,1))+  &
                        3.D0*HRR(1:LB,61,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,13,1)+  &
                        3.D0*HRR(1:LB,24,1))+  &
                        3.D0*HRR(1:LB,40,1))+  &
                        HRR(1:LB,62,1))+  &
                        HRR(1:LB,90,1)
      !=|14,24)
      HRR(1:LB,14,24)=CDy*(CDy*(CDy*HRR(1:LB,24,1)+  &
                        3.D0*HRR(1:LB,40,1))+  &
                        3.D0*HRR(1:LB,62,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,14,1)+  &
                        3.D0*HRR(1:LB,25,1))+  &
                        3.D0*HRR(1:LB,41,1))+  &
                        HRR(1:LB,63,1))+  &
                        HRR(1:LB,91,1)
      !=|15,24)
      HRR(1:LB,15,24)=CDy*(CDy*(CDy*HRR(1:LB,26,1)+  &
                        3.D0*HRR(1:LB,43,1))+  &
                        3.D0*HRR(1:LB,66,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,15,1)+  &
                        3.D0*HRR(1:LB,27,1))+  &
                        3.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,67,1))+  &
                        HRR(1:LB,96,1)
      !=|16,24)
      HRR(1:LB,16,24)=CDy*(CDy*(CDy*HRR(1:LB,27,1)+  &
                        3.D0*HRR(1:LB,44,1))+  &
                        3.D0*HRR(1:LB,67,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,16,1)+  &
                        3.D0*HRR(1:LB,28,1))+  &
                        3.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,68,1))+  &
                        HRR(1:LB,97,1)
      !=|17,24)
      HRR(1:LB,17,24)=CDy*(CDy*(CDy*HRR(1:LB,28,1)+  &
                        3.D0*HRR(1:LB,45,1))+  &
                        3.D0*HRR(1:LB,68,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,17,1)+  &
                        3.D0*HRR(1:LB,29,1))+  &
                        3.D0*HRR(1:LB,46,1))+  &
                        HRR(1:LB,69,1))+  &
                        HRR(1:LB,98,1)
      !=|18,24)
      HRR(1:LB,18,24)=CDy*(CDy*(CDy*HRR(1:LB,30,1)+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        3.D0*HRR(1:LB,72,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,18,1)+  &
                        3.D0*HRR(1:LB,31,1))+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,73,1))+  &
                        HRR(1:LB,103,1)
      !=|19,24)
      HRR(1:LB,19,24)=CDy*(CDy*(CDy*HRR(1:LB,31,1)+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        3.D0*HRR(1:LB,73,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,19,1)+  &
                        3.D0*HRR(1:LB,32,1))+  &
                        3.D0*HRR(1:LB,50,1))+  &
                        HRR(1:LB,74,1))+  &
                        HRR(1:LB,104,1)
      !=|20,24)
      HRR(1:LB,20,24)=CDy*(CDy*(CDy*HRR(1:LB,33,1)+  &
                        3.D0*HRR(1:LB,52,1))+  &
                        3.D0*HRR(1:LB,77,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,20,1)+  &
                        3.D0*HRR(1:LB,34,1))+  &
                        3.D0*HRR(1:LB,53,1))+  &
                        HRR(1:LB,78,1))+  &
                        HRR(1:LB,109,1)
      !=|11,25)
      HRR(1:LB,11,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,11,1)+  &
                        4.D0*HRR(1:LB,22,1))+  &
                        6.D0*HRR(1:LB,38,1))+  &
                        4.D0*HRR(1:LB,60,1))+  &
                        HRR(1:LB,89,1)
      !=|12,25)
      HRR(1:LB,12,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,12,1)+  &
                        4.D0*HRR(1:LB,23,1))+  &
                        6.D0*HRR(1:LB,39,1))+  &
                        4.D0*HRR(1:LB,61,1))+  &
                        HRR(1:LB,90,1)
      !=|13,25)
      HRR(1:LB,13,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,13,1)+  &
                        4.D0*HRR(1:LB,24,1))+  &
                        6.D0*HRR(1:LB,40,1))+  &
                        4.D0*HRR(1:LB,62,1))+  &
                        HRR(1:LB,91,1)
      !=|14,25)
      HRR(1:LB,14,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,14,1)+  &
                        4.D0*HRR(1:LB,25,1))+  &
                        6.D0*HRR(1:LB,41,1))+  &
                        4.D0*HRR(1:LB,63,1))+  &
                        HRR(1:LB,92,1)
      !=|15,25)
      HRR(1:LB,15,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,15,1)+  &
                        4.D0*HRR(1:LB,27,1))+  &
                        6.D0*HRR(1:LB,44,1))+  &
                        4.D0*HRR(1:LB,67,1))+  &
                        HRR(1:LB,97,1)
      !=|16,25)
      HRR(1:LB,16,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,16,1)+  &
                        4.D0*HRR(1:LB,28,1))+  &
                        6.D0*HRR(1:LB,45,1))+  &
                        4.D0*HRR(1:LB,68,1))+  &
                        HRR(1:LB,98,1)
      !=|17,25)
      HRR(1:LB,17,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,17,1)+  &
                        4.D0*HRR(1:LB,29,1))+  &
                        6.D0*HRR(1:LB,46,1))+  &
                        4.D0*HRR(1:LB,69,1))+  &
                        HRR(1:LB,99,1)
      !=|18,25)
      HRR(1:LB,18,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,18,1)+  &
                        4.D0*HRR(1:LB,31,1))+  &
                        6.D0*HRR(1:LB,49,1))+  &
                        4.D0*HRR(1:LB,73,1))+  &
                        HRR(1:LB,104,1)
      !=|19,25)
      HRR(1:LB,19,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,19,1)+  &
                        4.D0*HRR(1:LB,32,1))+  &
                        6.D0*HRR(1:LB,50,1))+  &
                        4.D0*HRR(1:LB,74,1))+  &
                        HRR(1:LB,105,1)
      !=|20,25)
      HRR(1:LB,20,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,20,1)+  &
                        4.D0*HRR(1:LB,34,1))+  &
                        6.D0*HRR(1:LB,53,1))+  &
                        4.D0*HRR(1:LB,78,1))+  &
                        HRR(1:LB,110,1)
      !=|11,26)
      HRR(1:LB,11,26)=CDz*HRR(1:LB,57,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,36,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,21,1)+  &
                        CDx*(CDz*HRR(1:LB,11,1)+  &
                        HRR(1:LB,26,1))+  &
                        3.D0*HRR(1:LB,42,1))+  &
                        3.D0*HRR(1:LB,64,1))+  &
                        HRR(1:LB,93,1)
      !=|12,26)
      HRR(1:LB,12,26)=CDz*HRR(1:LB,58,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,37,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,22,1)+  &
                        CDx*(CDz*HRR(1:LB,12,1)+  &
                        HRR(1:LB,27,1))+  &
                        3.D0*HRR(1:LB,43,1))+  &
                        3.D0*HRR(1:LB,65,1))+  &
                        HRR(1:LB,94,1)
      !=|13,26)
      HRR(1:LB,13,26)=CDz*HRR(1:LB,59,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,38,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,23,1)+  &
                        CDx*(CDz*HRR(1:LB,13,1)+  &
                        HRR(1:LB,28,1))+  &
                        3.D0*HRR(1:LB,44,1))+  &
                        3.D0*HRR(1:LB,66,1))+  &
                        HRR(1:LB,95,1)
      !=|14,26)
      HRR(1:LB,14,26)=CDz*HRR(1:LB,60,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,39,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,24,1)+  &
                        CDx*(CDz*HRR(1:LB,14,1)+  &
                        HRR(1:LB,29,1))+  &
                        3.D0*HRR(1:LB,45,1))+  &
                        3.D0*HRR(1:LB,67,1))+  &
                        HRR(1:LB,96,1)
      !=|15,26)
      HRR(1:LB,15,26)=CDz*HRR(1:LB,64,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,42,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,26,1)+  &
                        CDx*(CDz*HRR(1:LB,15,1)+  &
                        HRR(1:LB,30,1))+  &
                        3.D0*HRR(1:LB,47,1))+  &
                        3.D0*HRR(1:LB,70,1))+  &
                        HRR(1:LB,100,1)
      !=|16,26)
      HRR(1:LB,16,26)=CDz*HRR(1:LB,65,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,43,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,27,1)+  &
                        CDx*(CDz*HRR(1:LB,16,1)+  &
                        HRR(1:LB,31,1))+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        3.D0*HRR(1:LB,71,1))+  &
                        HRR(1:LB,101,1)
      !=|17,26)
      HRR(1:LB,17,26)=CDz*HRR(1:LB,66,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,44,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,28,1)+  &
                        CDx*(CDz*HRR(1:LB,17,1)+  &
                        HRR(1:LB,32,1))+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        3.D0*HRR(1:LB,72,1))+  &
                        HRR(1:LB,102,1)
      !=|18,26)
      HRR(1:LB,18,26)=CDz*HRR(1:LB,70,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,47,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,30,1)+  &
                        CDx*(CDz*HRR(1:LB,18,1)+  &
                        HRR(1:LB,33,1))+  &
                        3.D0*HRR(1:LB,51,1))+  &
                        3.D0*HRR(1:LB,75,1))+  &
                        HRR(1:LB,106,1)
      !=|19,26)
      HRR(1:LB,19,26)=CDz*HRR(1:LB,71,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,48,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,31,1)+  &
                        CDx*(CDz*HRR(1:LB,19,1)+  &
                        HRR(1:LB,34,1))+  &
                        3.D0*HRR(1:LB,52,1))+  &
                        3.D0*HRR(1:LB,76,1))+  &
                        HRR(1:LB,107,1)
      !=|20,26)
      HRR(1:LB,20,26)=CDz*HRR(1:LB,75,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,51,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,33,1)+  &
                        CDx*(CDz*HRR(1:LB,20,1)+  &
                        HRR(1:LB,35,1))+  &
                        3.D0*HRR(1:LB,54,1))+  &
                        3.D0*HRR(1:LB,79,1))+  &
                        HRR(1:LB,111,1)
      !=|11,27)
      HRR(1:LB,11,27)=CDz*HRR(1:LB,58,1)+  &
                        CDy*(CDz*HRR(1:LB,36,1)+  &
                        HRR(1:LB,64,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,37,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,21,1)+  &
                        2.D0*HRR(1:LB,42,1))+  &
                        CDx*(CDz*HRR(1:LB,22,1)+  &
                        CDy*(CDz*HRR(1:LB,11,1)+  &
                        HRR(1:LB,26,1))+  &
                        HRR(1:LB,43,1))+  &
                        2.D0*HRR(1:LB,65,1))+  &
                        HRR(1:LB,94,1)
      !=|12,27)
      HRR(1:LB,12,27)=CDz*HRR(1:LB,59,1)+  &
                        CDy*(CDz*HRR(1:LB,37,1)+  &
                        HRR(1:LB,65,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,38,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,22,1)+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        CDx*(CDz*HRR(1:LB,23,1)+  &
                        CDy*(CDz*HRR(1:LB,12,1)+  &
                        HRR(1:LB,27,1))+  &
                        HRR(1:LB,44,1))+  &
                        2.D0*HRR(1:LB,66,1))+  &
                        HRR(1:LB,95,1)
      !=|13,27)
      HRR(1:LB,13,27)=CDz*HRR(1:LB,60,1)+  &
                        CDy*(CDz*HRR(1:LB,38,1)+  &
                        HRR(1:LB,66,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,39,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        CDx*(CDz*HRR(1:LB,24,1)+  &
                        CDy*(CDz*HRR(1:LB,13,1)+  &
                        HRR(1:LB,28,1))+  &
                        HRR(1:LB,45,1))+  &
                        2.D0*HRR(1:LB,67,1))+  &
                        HRR(1:LB,96,1)
      !=|14,27)
      HRR(1:LB,14,27)=CDz*HRR(1:LB,61,1)+  &
                        CDy*(CDz*HRR(1:LB,39,1)+  &
                        HRR(1:LB,67,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,40,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,24,1)+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        CDx*(CDz*HRR(1:LB,25,1)+  &
                        CDy*(CDz*HRR(1:LB,14,1)+  &
                        HRR(1:LB,29,1))+  &
                        HRR(1:LB,46,1))+  &
                        2.D0*HRR(1:LB,68,1))+  &
                        HRR(1:LB,97,1)
      !=|15,27)
      HRR(1:LB,15,27)=CDz*HRR(1:LB,65,1)+  &
                        CDy*(CDz*HRR(1:LB,42,1)+  &
                        HRR(1:LB,70,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,43,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,26,1)+  &
                        2.D0*HRR(1:LB,47,1))+  &
                        CDx*(CDz*HRR(1:LB,27,1)+  &
                        CDy*(CDz*HRR(1:LB,15,1)+  &
                        HRR(1:LB,30,1))+  &
                        HRR(1:LB,48,1))+  &
                        2.D0*HRR(1:LB,71,1))+  &
                        HRR(1:LB,101,1)
      !=|16,27)
      HRR(1:LB,16,27)=CDz*HRR(1:LB,66,1)+  &
                        CDy*(CDz*HRR(1:LB,43,1)+  &
                        HRR(1:LB,71,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,44,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,27,1)+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        CDx*(CDz*HRR(1:LB,28,1)+  &
                        CDy*(CDz*HRR(1:LB,16,1)+  &
                        HRR(1:LB,31,1))+  &
                        HRR(1:LB,49,1))+  &
                        2.D0*HRR(1:LB,72,1))+  &
                        HRR(1:LB,102,1)
      !=|17,27)
      HRR(1:LB,17,27)=CDz*HRR(1:LB,67,1)+  &
                        CDy*(CDz*HRR(1:LB,44,1)+  &
                        HRR(1:LB,72,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,45,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,28,1)+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        CDx*(CDz*HRR(1:LB,29,1)+  &
                        CDy*(CDz*HRR(1:LB,17,1)+  &
                        HRR(1:LB,32,1))+  &
                        HRR(1:LB,50,1))+  &
                        2.D0*HRR(1:LB,73,1))+  &
                        HRR(1:LB,103,1)
      !=|18,27)
      HRR(1:LB,18,27)=CDz*HRR(1:LB,71,1)+  &
                        CDy*(CDz*HRR(1:LB,47,1)+  &
                        HRR(1:LB,75,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,48,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,30,1)+  &
                        2.D0*HRR(1:LB,51,1))+  &
                        CDx*(CDz*HRR(1:LB,31,1)+  &
                        CDy*(CDz*HRR(1:LB,18,1)+  &
                        HRR(1:LB,33,1))+  &
                        HRR(1:LB,52,1))+  &
                        2.D0*HRR(1:LB,76,1))+  &
                        HRR(1:LB,107,1)
      !=|19,27)
      HRR(1:LB,19,27)=CDz*HRR(1:LB,72,1)+  &
                        CDy*(CDz*HRR(1:LB,48,1)+  &
                        HRR(1:LB,76,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,49,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,31,1)+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        CDx*(CDz*HRR(1:LB,32,1)+  &
                        CDy*(CDz*HRR(1:LB,19,1)+  &
                        HRR(1:LB,34,1))+  &
                        HRR(1:LB,53,1))+  &
                        2.D0*HRR(1:LB,77,1))+  &
                        HRR(1:LB,108,1)
      !=|20,27)
      HRR(1:LB,20,27)=CDz*HRR(1:LB,76,1)+  &
                        CDy*(CDz*HRR(1:LB,51,1)+  &
                        HRR(1:LB,79,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,52,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,33,1)+  &
                        2.D0*HRR(1:LB,54,1))+  &
                        CDx*(CDz*HRR(1:LB,34,1)+  &
                        CDy*(CDz*HRR(1:LB,20,1)+  &
                        HRR(1:LB,35,1))+  &
                        HRR(1:LB,55,1))+  &
                        2.D0*HRR(1:LB,80,1))+  &
                        HRR(1:LB,112,1)
      !=|11,28)
      HRR(1:LB,11,28)=CDz*HRR(1:LB,59,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,37,1)+  &
                        CDy*(CDz*HRR(1:LB,21,1)+  &
                        HRR(1:LB,42,1))+  &
                        2.D0*HRR(1:LB,65,1))+  &
                        CDx*(CDz*HRR(1:LB,38,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,22,1)+  &
                        CDy*(CDz*HRR(1:LB,11,1)+  &
                        HRR(1:LB,26,1))+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,66,1))+  &
                        HRR(1:LB,95,1)
      !=|12,28)
      HRR(1:LB,12,28)=CDz*HRR(1:LB,60,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,38,1)+  &
                        CDy*(CDz*HRR(1:LB,22,1)+  &
                        HRR(1:LB,43,1))+  &
                        2.D0*HRR(1:LB,66,1))+  &
                        CDx*(CDz*HRR(1:LB,39,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,23,1)+  &
                        CDy*(CDz*HRR(1:LB,12,1)+  &
                        HRR(1:LB,27,1))+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,67,1))+  &
                        HRR(1:LB,96,1)
      !=|13,28)
      HRR(1:LB,13,28)=CDz*HRR(1:LB,61,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,39,1)+  &
                        CDy*(CDz*HRR(1:LB,23,1)+  &
                        HRR(1:LB,44,1))+  &
                        2.D0*HRR(1:LB,67,1))+  &
                        CDx*(CDz*HRR(1:LB,40,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,24,1)+  &
                        CDy*(CDz*HRR(1:LB,13,1)+  &
                        HRR(1:LB,28,1))+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,68,1))+  &
                        HRR(1:LB,97,1)
      !=|14,28)
      HRR(1:LB,14,28)=CDz*HRR(1:LB,62,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,40,1)+  &
                        CDy*(CDz*HRR(1:LB,24,1)+  &
                        HRR(1:LB,45,1))+  &
                        2.D0*HRR(1:LB,68,1))+  &
                        CDx*(CDz*HRR(1:LB,41,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,25,1)+  &
                        CDy*(CDz*HRR(1:LB,14,1)+  &
                        HRR(1:LB,29,1))+  &
                        2.D0*HRR(1:LB,46,1))+  &
                        HRR(1:LB,69,1))+  &
                        HRR(1:LB,98,1)
      !=|15,28)
      HRR(1:LB,15,28)=CDz*HRR(1:LB,66,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,43,1)+  &
                        CDy*(CDz*HRR(1:LB,26,1)+  &
                        HRR(1:LB,47,1))+  &
                        2.D0*HRR(1:LB,71,1))+  &
                        CDx*(CDz*HRR(1:LB,44,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,27,1)+  &
                        CDy*(CDz*HRR(1:LB,15,1)+  &
                        HRR(1:LB,30,1))+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,72,1))+  &
                        HRR(1:LB,102,1)
      !=|16,28)
      HRR(1:LB,16,28)=CDz*HRR(1:LB,67,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,44,1)+  &
                        CDy*(CDz*HRR(1:LB,27,1)+  &
                        HRR(1:LB,48,1))+  &
                        2.D0*HRR(1:LB,72,1))+  &
                        CDx*(CDz*HRR(1:LB,45,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,28,1)+  &
                        CDy*(CDz*HRR(1:LB,16,1)+  &
                        HRR(1:LB,31,1))+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,73,1))+  &
                        HRR(1:LB,103,1)
      !=|17,28)
      HRR(1:LB,17,28)=CDz*HRR(1:LB,68,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,45,1)+  &
                        CDy*(CDz*HRR(1:LB,28,1)+  &
                        HRR(1:LB,49,1))+  &
                        2.D0*HRR(1:LB,73,1))+  &
                        CDx*(CDz*HRR(1:LB,46,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,29,1)+  &
                        CDy*(CDz*HRR(1:LB,17,1)+  &
                        HRR(1:LB,32,1))+  &
                        2.D0*HRR(1:LB,50,1))+  &
                        HRR(1:LB,74,1))+  &
                        HRR(1:LB,104,1)
      !=|18,28)
      HRR(1:LB,18,28)=CDz*HRR(1:LB,72,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,48,1)+  &
                        CDy*(CDz*HRR(1:LB,30,1)+  &
                        HRR(1:LB,51,1))+  &
                        2.D0*HRR(1:LB,76,1))+  &
                        CDx*(CDz*HRR(1:LB,49,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,31,1)+  &
                        CDy*(CDz*HRR(1:LB,18,1)+  &
                        HRR(1:LB,33,1))+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,77,1))+  &
                        HRR(1:LB,108,1)
      !=|19,28)
      HRR(1:LB,19,28)=CDz*HRR(1:LB,73,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,49,1)+  &
                        CDy*(CDz*HRR(1:LB,31,1)+  &
                        HRR(1:LB,52,1))+  &
                        2.D0*HRR(1:LB,77,1))+  &
                        CDx*(CDz*HRR(1:LB,50,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,32,1)+  &
                        CDy*(CDz*HRR(1:LB,19,1)+  &
                        HRR(1:LB,34,1))+  &
                        2.D0*HRR(1:LB,53,1))+  &
                        HRR(1:LB,78,1))+  &
                        HRR(1:LB,109,1)
      !=|20,28)
      HRR(1:LB,20,28)=CDz*HRR(1:LB,77,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,52,1)+  &
                        CDy*(CDz*HRR(1:LB,33,1)+  &
                        HRR(1:LB,54,1))+  &
                        2.D0*HRR(1:LB,80,1))+  &
                        CDx*(CDz*HRR(1:LB,53,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,34,1)+  &
                        CDy*(CDz*HRR(1:LB,20,1)+  &
                        HRR(1:LB,35,1))+  &
                        2.D0*HRR(1:LB,55,1))+  &
                        HRR(1:LB,81,1))+  &
                        HRR(1:LB,113,1)
      !=|11,29)
      HRR(1:LB,11,29)=CDz*HRR(1:LB,60,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,38,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,22,1)+  &
                        CDy*(CDz*HRR(1:LB,11,1)+  &
                        HRR(1:LB,26,1))+  &
                        3.D0*HRR(1:LB,43,1))+  &
                        3.D0*HRR(1:LB,66,1))+  &
                        HRR(1:LB,96,1)
      !=|12,29)
      HRR(1:LB,12,29)=CDz*HRR(1:LB,61,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,39,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,23,1)+  &
                        CDy*(CDz*HRR(1:LB,12,1)+  &
                        HRR(1:LB,27,1))+  &
                        3.D0*HRR(1:LB,44,1))+  &
                        3.D0*HRR(1:LB,67,1))+  &
                        HRR(1:LB,97,1)
      !=|13,29)
      HRR(1:LB,13,29)=CDz*HRR(1:LB,62,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,40,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,24,1)+  &
                        CDy*(CDz*HRR(1:LB,13,1)+  &
                        HRR(1:LB,28,1))+  &
                        3.D0*HRR(1:LB,45,1))+  &
                        3.D0*HRR(1:LB,68,1))+  &
                        HRR(1:LB,98,1)
      !=|14,29)
      HRR(1:LB,14,29)=CDz*HRR(1:LB,63,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,41,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,25,1)+  &
                        CDy*(CDz*HRR(1:LB,14,1)+  &
                        HRR(1:LB,29,1))+  &
                        3.D0*HRR(1:LB,46,1))+  &
                        3.D0*HRR(1:LB,69,1))+  &
                        HRR(1:LB,99,1)
      !=|15,29)
      HRR(1:LB,15,29)=CDz*HRR(1:LB,67,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,44,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,27,1)+  &
                        CDy*(CDz*HRR(1:LB,15,1)+  &
                        HRR(1:LB,30,1))+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        3.D0*HRR(1:LB,72,1))+  &
                        HRR(1:LB,103,1)
      !=|16,29)
      HRR(1:LB,16,29)=CDz*HRR(1:LB,68,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,45,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,28,1)+  &
                        CDy*(CDz*HRR(1:LB,16,1)+  &
                        HRR(1:LB,31,1))+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        3.D0*HRR(1:LB,73,1))+  &
                        HRR(1:LB,104,1)
      !=|17,29)
      HRR(1:LB,17,29)=CDz*HRR(1:LB,69,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,46,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,29,1)+  &
                        CDy*(CDz*HRR(1:LB,17,1)+  &
                        HRR(1:LB,32,1))+  &
                        3.D0*HRR(1:LB,50,1))+  &
                        3.D0*HRR(1:LB,74,1))+  &
                        HRR(1:LB,105,1)
      !=|18,29)
      HRR(1:LB,18,29)=CDz*HRR(1:LB,73,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,49,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,31,1)+  &
                        CDy*(CDz*HRR(1:LB,18,1)+  &
                        HRR(1:LB,33,1))+  &
                        3.D0*HRR(1:LB,52,1))+  &
                        3.D0*HRR(1:LB,77,1))+  &
                        HRR(1:LB,109,1)
      !=|19,29)
      HRR(1:LB,19,29)=CDz*HRR(1:LB,74,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,50,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,32,1)+  &
                        CDy*(CDz*HRR(1:LB,19,1)+  &
                        HRR(1:LB,34,1))+  &
                        3.D0*HRR(1:LB,53,1))+  &
                        3.D0*HRR(1:LB,78,1))+  &
                        HRR(1:LB,110,1)
      !=|20,29)
      HRR(1:LB,20,29)=CDz*HRR(1:LB,78,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,53,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,34,1)+  &
                        CDy*(CDz*HRR(1:LB,20,1)+  &
                        HRR(1:LB,35,1))+  &
                        3.D0*HRR(1:LB,55,1))+  &
                        3.D0*HRR(1:LB,81,1))+  &
                        HRR(1:LB,114,1)
      !=|11,30)
      HRR(1:LB,11,30)=CDz*(CDz*HRR(1:LB,36,1)+  &
                        2.D0*HRR(1:LB,64,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,21,1)+  &
                        4.D0*HRR(1:LB,42,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,11,1)+  &
                        2.D0*HRR(1:LB,26,1))+  &
                        HRR(1:LB,47,1))+  &
                        2.D0*HRR(1:LB,70,1))+  &
                        HRR(1:LB,100,1)
      !=|12,30)
      HRR(1:LB,12,30)=CDz*(CDz*HRR(1:LB,37,1)+  &
                        2.D0*HRR(1:LB,65,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,22,1)+  &
                        4.D0*HRR(1:LB,43,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,12,1)+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        HRR(1:LB,48,1))+  &
                        2.D0*HRR(1:LB,71,1))+  &
                        HRR(1:LB,101,1)
      !=|13,30)
      HRR(1:LB,13,30)=CDz*(CDz*HRR(1:LB,38,1)+  &
                        2.D0*HRR(1:LB,66,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,23,1)+  &
                        4.D0*HRR(1:LB,44,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,13,1)+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,49,1))+  &
                        2.D0*HRR(1:LB,72,1))+  &
                        HRR(1:LB,102,1)
      !=|14,30)
      HRR(1:LB,14,30)=CDz*(CDz*HRR(1:LB,39,1)+  &
                        2.D0*HRR(1:LB,67,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,24,1)+  &
                        4.D0*HRR(1:LB,45,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,14,1)+  &
                        2.D0*HRR(1:LB,29,1))+  &
                        HRR(1:LB,50,1))+  &
                        2.D0*HRR(1:LB,73,1))+  &
                        HRR(1:LB,103,1)
      !=|15,30)
      HRR(1:LB,15,30)=CDz*(CDz*HRR(1:LB,42,1)+  &
                        2.D0*HRR(1:LB,70,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,26,1)+  &
                        4.D0*HRR(1:LB,47,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,15,1)+  &
                        2.D0*HRR(1:LB,30,1))+  &
                        HRR(1:LB,51,1))+  &
                        2.D0*HRR(1:LB,75,1))+  &
                        HRR(1:LB,106,1)
      !=|16,30)
      HRR(1:LB,16,30)=CDz*(CDz*HRR(1:LB,43,1)+  &
                        2.D0*HRR(1:LB,71,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,27,1)+  &
                        4.D0*HRR(1:LB,48,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,16,1)+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,52,1))+  &
                        2.D0*HRR(1:LB,76,1))+  &
                        HRR(1:LB,107,1)
      !=|17,30)
      HRR(1:LB,17,30)=CDz*(CDz*HRR(1:LB,44,1)+  &
                        2.D0*HRR(1:LB,72,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,28,1)+  &
                        4.D0*HRR(1:LB,49,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,17,1)+  &
                        2.D0*HRR(1:LB,32,1))+  &
                        HRR(1:LB,53,1))+  &
                        2.D0*HRR(1:LB,77,1))+  &
                        HRR(1:LB,108,1)
      !=|18,30)
      HRR(1:LB,18,30)=CDz*(CDz*HRR(1:LB,47,1)+  &
                        2.D0*HRR(1:LB,75,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,30,1)+  &
                        4.D0*HRR(1:LB,51,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,18,1)+  &
                        2.D0*HRR(1:LB,33,1))+  &
                        HRR(1:LB,54,1))+  &
                        2.D0*HRR(1:LB,79,1))+  &
                        HRR(1:LB,111,1)
      !=|19,30)
      HRR(1:LB,19,30)=CDz*(CDz*HRR(1:LB,48,1)+  &
                        2.D0*HRR(1:LB,76,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,31,1)+  &
                        4.D0*HRR(1:LB,52,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,19,1)+  &
                        2.D0*HRR(1:LB,34,1))+  &
                        HRR(1:LB,55,1))+  &
                        2.D0*HRR(1:LB,80,1))+  &
                        HRR(1:LB,112,1)
      !=|20,30)
      HRR(1:LB,20,30)=CDz*(CDz*HRR(1:LB,51,1)+  &
                        2.D0*HRR(1:LB,79,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,33,1)+  &
                        4.D0*HRR(1:LB,54,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,20,1)+  &
                        2.D0*HRR(1:LB,35,1))+  &
                        HRR(1:LB,56,1))+  &
                        2.D0*HRR(1:LB,82,1))+  &
                        HRR(1:LB,115,1)
      !=|11,31)
      HRR(1:LB,11,31)=CDz*(CDz*HRR(1:LB,37,1)+  &
                        2.D0*HRR(1:LB,65,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,21,1)+  &
                        2.D0*HRR(1:LB,42,1))+  &
                        HRR(1:LB,70,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,22,1)+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,11,1)+  &
                        2.D0*HRR(1:LB,26,1))+  &
                        HRR(1:LB,47,1))+  &
                        HRR(1:LB,71,1))+  &
                        HRR(1:LB,101,1)
      !=|12,31)
      HRR(1:LB,12,31)=CDz*(CDz*HRR(1:LB,38,1)+  &
                        2.D0*HRR(1:LB,66,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,22,1)+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,71,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,12,1)+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        HRR(1:LB,48,1))+  &
                        HRR(1:LB,72,1))+  &
                        HRR(1:LB,102,1)
      !=|13,31)
      HRR(1:LB,13,31)=CDz*(CDz*HRR(1:LB,39,1)+  &
                        2.D0*HRR(1:LB,67,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,72,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,24,1)+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,13,1)+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,49,1))+  &
                        HRR(1:LB,73,1))+  &
                        HRR(1:LB,103,1)
      !=|14,31)
      HRR(1:LB,14,31)=CDz*(CDz*HRR(1:LB,40,1)+  &
                        2.D0*HRR(1:LB,68,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,24,1)+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,73,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,25,1)+  &
                        2.D0*HRR(1:LB,46,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,14,1)+  &
                        2.D0*HRR(1:LB,29,1))+  &
                        HRR(1:LB,50,1))+  &
                        HRR(1:LB,74,1))+  &
                        HRR(1:LB,104,1)
      !=|15,31)
      HRR(1:LB,15,31)=CDz*(CDz*HRR(1:LB,43,1)+  &
                        2.D0*HRR(1:LB,71,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,26,1)+  &
                        2.D0*HRR(1:LB,47,1))+  &
                        HRR(1:LB,75,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,27,1)+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,15,1)+  &
                        2.D0*HRR(1:LB,30,1))+  &
                        HRR(1:LB,51,1))+  &
                        HRR(1:LB,76,1))+  &
                        HRR(1:LB,107,1)
      !=|16,31)
      HRR(1:LB,16,31)=CDz*(CDz*HRR(1:LB,44,1)+  &
                        2.D0*HRR(1:LB,72,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,27,1)+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,76,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,28,1)+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,16,1)+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,52,1))+  &
                        HRR(1:LB,77,1))+  &
                        HRR(1:LB,108,1)
      !=|17,31)
      HRR(1:LB,17,31)=CDz*(CDz*HRR(1:LB,45,1)+  &
                        2.D0*HRR(1:LB,73,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,28,1)+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,77,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,29,1)+  &
                        2.D0*HRR(1:LB,50,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,17,1)+  &
                        2.D0*HRR(1:LB,32,1))+  &
                        HRR(1:LB,53,1))+  &
                        HRR(1:LB,78,1))+  &
                        HRR(1:LB,109,1)
      !=|18,31)
      HRR(1:LB,18,31)=CDz*(CDz*HRR(1:LB,48,1)+  &
                        2.D0*HRR(1:LB,76,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,30,1)+  &
                        2.D0*HRR(1:LB,51,1))+  &
                        HRR(1:LB,79,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,31,1)+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,18,1)+  &
                        2.D0*HRR(1:LB,33,1))+  &
                        HRR(1:LB,54,1))+  &
                        HRR(1:LB,80,1))+  &
                        HRR(1:LB,112,1)
      !=|19,31)
      HRR(1:LB,19,31)=CDz*(CDz*HRR(1:LB,49,1)+  &
                        2.D0*HRR(1:LB,77,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,31,1)+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,80,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,32,1)+  &
                        2.D0*HRR(1:LB,53,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,19,1)+  &
                        2.D0*HRR(1:LB,34,1))+  &
                        HRR(1:LB,55,1))+  &
                        HRR(1:LB,81,1))+  &
                        HRR(1:LB,113,1)
      !=|20,31)
      HRR(1:LB,20,31)=CDz*(CDz*HRR(1:LB,52,1)+  &
                        2.D0*HRR(1:LB,80,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,33,1)+  &
                        2.D0*HRR(1:LB,54,1))+  &
                        HRR(1:LB,82,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,34,1)+  &
                        2.D0*HRR(1:LB,55,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,20,1)+  &
                        2.D0*HRR(1:LB,35,1))+  &
                        HRR(1:LB,56,1))+  &
                        HRR(1:LB,83,1))+  &
                        HRR(1:LB,116,1)
      !=|11,32)
      HRR(1:LB,11,32)=CDz*(CDz*HRR(1:LB,38,1)+  &
                        2.D0*HRR(1:LB,66,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,22,1)+  &
                        4.D0*HRR(1:LB,43,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,11,1)+  &
                        2.D0*HRR(1:LB,26,1))+  &
                        HRR(1:LB,47,1))+  &
                        2.D0*HRR(1:LB,71,1))+  &
                        HRR(1:LB,102,1)
      !=|12,32)
      HRR(1:LB,12,32)=CDz*(CDz*HRR(1:LB,39,1)+  &
                        2.D0*HRR(1:LB,67,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,23,1)+  &
                        4.D0*HRR(1:LB,44,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,12,1)+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        HRR(1:LB,48,1))+  &
                        2.D0*HRR(1:LB,72,1))+  &
                        HRR(1:LB,103,1)
      !=|13,32)
      HRR(1:LB,13,32)=CDz*(CDz*HRR(1:LB,40,1)+  &
                        2.D0*HRR(1:LB,68,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,24,1)+  &
                        4.D0*HRR(1:LB,45,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,13,1)+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,49,1))+  &
                        2.D0*HRR(1:LB,73,1))+  &
                        HRR(1:LB,104,1)
      !=|14,32)
      HRR(1:LB,14,32)=CDz*(CDz*HRR(1:LB,41,1)+  &
                        2.D0*HRR(1:LB,69,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,25,1)+  &
                        4.D0*HRR(1:LB,46,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,14,1)+  &
                        2.D0*HRR(1:LB,29,1))+  &
                        HRR(1:LB,50,1))+  &
                        2.D0*HRR(1:LB,74,1))+  &
                        HRR(1:LB,105,1)
      !=|15,32)
      HRR(1:LB,15,32)=CDz*(CDz*HRR(1:LB,44,1)+  &
                        2.D0*HRR(1:LB,72,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,27,1)+  &
                        4.D0*HRR(1:LB,48,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,15,1)+  &
                        2.D0*HRR(1:LB,30,1))+  &
                        HRR(1:LB,51,1))+  &
                        2.D0*HRR(1:LB,76,1))+  &
                        HRR(1:LB,108,1)
      !=|16,32)
      HRR(1:LB,16,32)=CDz*(CDz*HRR(1:LB,45,1)+  &
                        2.D0*HRR(1:LB,73,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,28,1)+  &
                        4.D0*HRR(1:LB,49,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,16,1)+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,52,1))+  &
                        2.D0*HRR(1:LB,77,1))+  &
                        HRR(1:LB,109,1)
      !=|17,32)
      HRR(1:LB,17,32)=CDz*(CDz*HRR(1:LB,46,1)+  &
                        2.D0*HRR(1:LB,74,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,29,1)+  &
                        4.D0*HRR(1:LB,50,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,17,1)+  &
                        2.D0*HRR(1:LB,32,1))+  &
                        HRR(1:LB,53,1))+  &
                        2.D0*HRR(1:LB,78,1))+  &
                        HRR(1:LB,110,1)
      !=|18,32)
      HRR(1:LB,18,32)=CDz*(CDz*HRR(1:LB,49,1)+  &
                        2.D0*HRR(1:LB,77,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,31,1)+  &
                        4.D0*HRR(1:LB,52,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,18,1)+  &
                        2.D0*HRR(1:LB,33,1))+  &
                        HRR(1:LB,54,1))+  &
                        2.D0*HRR(1:LB,80,1))+  &
                        HRR(1:LB,113,1)
      !=|19,32)
      HRR(1:LB,19,32)=CDz*(CDz*HRR(1:LB,50,1)+  &
                        2.D0*HRR(1:LB,78,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,32,1)+  &
                        4.D0*HRR(1:LB,53,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,19,1)+  &
                        2.D0*HRR(1:LB,34,1))+  &
                        HRR(1:LB,55,1))+  &
                        2.D0*HRR(1:LB,81,1))+  &
                        HRR(1:LB,114,1)
      !=|20,32)
      HRR(1:LB,20,32)=CDz*(CDz*HRR(1:LB,53,1)+  &
                        2.D0*HRR(1:LB,81,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,34,1)+  &
                        4.D0*HRR(1:LB,55,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,20,1)+  &
                        2.D0*HRR(1:LB,35,1))+  &
                        HRR(1:LB,56,1))+  &
                        2.D0*HRR(1:LB,83,1))+  &
                        HRR(1:LB,117,1)
      !=|11,33)
      HRR(1:LB,11,33)=CDz*(CDz*(CDz*HRR(1:LB,21,1)+  &
                        3.D0*HRR(1:LB,42,1))+  &
                        3.D0*HRR(1:LB,70,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,11,1)+  &
                        3.D0*HRR(1:LB,26,1))+  &
                        3.D0*HRR(1:LB,47,1))+  &
                        HRR(1:LB,75,1))+  &
                        HRR(1:LB,106,1)
      !=|12,33)
      HRR(1:LB,12,33)=CDz*(CDz*(CDz*HRR(1:LB,22,1)+  &
                        3.D0*HRR(1:LB,43,1))+  &
                        3.D0*HRR(1:LB,71,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,12,1)+  &
                        3.D0*HRR(1:LB,27,1))+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,76,1))+  &
                        HRR(1:LB,107,1)
      !=|13,33)
      HRR(1:LB,13,33)=CDz*(CDz*(CDz*HRR(1:LB,23,1)+  &
                        3.D0*HRR(1:LB,44,1))+  &
                        3.D0*HRR(1:LB,72,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,13,1)+  &
                        3.D0*HRR(1:LB,28,1))+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,77,1))+  &
                        HRR(1:LB,108,1)
      !=|14,33)
      HRR(1:LB,14,33)=CDz*(CDz*(CDz*HRR(1:LB,24,1)+  &
                        3.D0*HRR(1:LB,45,1))+  &
                        3.D0*HRR(1:LB,73,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,14,1)+  &
                        3.D0*HRR(1:LB,29,1))+  &
                        3.D0*HRR(1:LB,50,1))+  &
                        HRR(1:LB,78,1))+  &
                        HRR(1:LB,109,1)
      !=|15,33)
      HRR(1:LB,15,33)=CDz*(CDz*(CDz*HRR(1:LB,26,1)+  &
                        3.D0*HRR(1:LB,47,1))+  &
                        3.D0*HRR(1:LB,75,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,15,1)+  &
                        3.D0*HRR(1:LB,30,1))+  &
                        3.D0*HRR(1:LB,51,1))+  &
                        HRR(1:LB,79,1))+  &
                        HRR(1:LB,111,1)
      !=|16,33)
      HRR(1:LB,16,33)=CDz*(CDz*(CDz*HRR(1:LB,27,1)+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        3.D0*HRR(1:LB,76,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,16,1)+  &
                        3.D0*HRR(1:LB,31,1))+  &
                        3.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,80,1))+  &
                        HRR(1:LB,112,1)
      !=|17,33)
      HRR(1:LB,17,33)=CDz*(CDz*(CDz*HRR(1:LB,28,1)+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        3.D0*HRR(1:LB,77,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,17,1)+  &
                        3.D0*HRR(1:LB,32,1))+  &
                        3.D0*HRR(1:LB,53,1))+  &
                        HRR(1:LB,81,1))+  &
                        HRR(1:LB,113,1)
      !=|18,33)
      HRR(1:LB,18,33)=CDz*(CDz*(CDz*HRR(1:LB,30,1)+  &
                        3.D0*HRR(1:LB,51,1))+  &
                        3.D0*HRR(1:LB,79,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,18,1)+  &
                        3.D0*HRR(1:LB,33,1))+  &
                        3.D0*HRR(1:LB,54,1))+  &
                        HRR(1:LB,82,1))+  &
                        HRR(1:LB,115,1)
      !=|19,33)
      HRR(1:LB,19,33)=CDz*(CDz*(CDz*HRR(1:LB,31,1)+  &
                        3.D0*HRR(1:LB,52,1))+  &
                        3.D0*HRR(1:LB,80,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,19,1)+  &
                        3.D0*HRR(1:LB,34,1))+  &
                        3.D0*HRR(1:LB,55,1))+  &
                        HRR(1:LB,83,1))+  &
                        HRR(1:LB,116,1)
      !=|20,33)
      HRR(1:LB,20,33)=CDz*(CDz*(CDz*HRR(1:LB,33,1)+  &
                        3.D0*HRR(1:LB,54,1))+  &
                        3.D0*HRR(1:LB,82,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,20,1)+  &
                        3.D0*HRR(1:LB,35,1))+  &
                        3.D0*HRR(1:LB,56,1))+  &
                        HRR(1:LB,84,1))+  &
                        HRR(1:LB,118,1)
      !=|11,34)
      HRR(1:LB,11,34)=CDz*(CDz*(CDz*HRR(1:LB,22,1)+  &
                        3.D0*HRR(1:LB,43,1))+  &
                        3.D0*HRR(1:LB,71,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,11,1)+  &
                        3.D0*HRR(1:LB,26,1))+  &
                        3.D0*HRR(1:LB,47,1))+  &
                        HRR(1:LB,75,1))+  &
                        HRR(1:LB,107,1)
      !=|12,34)
      HRR(1:LB,12,34)=CDz*(CDz*(CDz*HRR(1:LB,23,1)+  &
                        3.D0*HRR(1:LB,44,1))+  &
                        3.D0*HRR(1:LB,72,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,12,1)+  &
                        3.D0*HRR(1:LB,27,1))+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,76,1))+  &
                        HRR(1:LB,108,1)
      !=|13,34)
      HRR(1:LB,13,34)=CDz*(CDz*(CDz*HRR(1:LB,24,1)+  &
                        3.D0*HRR(1:LB,45,1))+  &
                        3.D0*HRR(1:LB,73,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,13,1)+  &
                        3.D0*HRR(1:LB,28,1))+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,77,1))+  &
                        HRR(1:LB,109,1)
      !=|14,34)
      HRR(1:LB,14,34)=CDz*(CDz*(CDz*HRR(1:LB,25,1)+  &
                        3.D0*HRR(1:LB,46,1))+  &
                        3.D0*HRR(1:LB,74,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,14,1)+  &
                        3.D0*HRR(1:LB,29,1))+  &
                        3.D0*HRR(1:LB,50,1))+  &
                        HRR(1:LB,78,1))+  &
                        HRR(1:LB,110,1)
      !=|15,34)
      HRR(1:LB,15,34)=CDz*(CDz*(CDz*HRR(1:LB,27,1)+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        3.D0*HRR(1:LB,76,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,15,1)+  &
                        3.D0*HRR(1:LB,30,1))+  &
                        3.D0*HRR(1:LB,51,1))+  &
                        HRR(1:LB,79,1))+  &
                        HRR(1:LB,112,1)
      !=|16,34)
      HRR(1:LB,16,34)=CDz*(CDz*(CDz*HRR(1:LB,28,1)+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        3.D0*HRR(1:LB,77,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,16,1)+  &
                        3.D0*HRR(1:LB,31,1))+  &
                        3.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,80,1))+  &
                        HRR(1:LB,113,1)
      !=|17,34)
      HRR(1:LB,17,34)=CDz*(CDz*(CDz*HRR(1:LB,29,1)+  &
                        3.D0*HRR(1:LB,50,1))+  &
                        3.D0*HRR(1:LB,78,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,17,1)+  &
                        3.D0*HRR(1:LB,32,1))+  &
                        3.D0*HRR(1:LB,53,1))+  &
                        HRR(1:LB,81,1))+  &
                        HRR(1:LB,114,1)
      !=|18,34)
      HRR(1:LB,18,34)=CDz*(CDz*(CDz*HRR(1:LB,31,1)+  &
                        3.D0*HRR(1:LB,52,1))+  &
                        3.D0*HRR(1:LB,80,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,18,1)+  &
                        3.D0*HRR(1:LB,33,1))+  &
                        3.D0*HRR(1:LB,54,1))+  &
                        HRR(1:LB,82,1))+  &
                        HRR(1:LB,116,1)
      !=|19,34)
      HRR(1:LB,19,34)=CDz*(CDz*(CDz*HRR(1:LB,32,1)+  &
                        3.D0*HRR(1:LB,53,1))+  &
                        3.D0*HRR(1:LB,81,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,19,1)+  &
                        3.D0*HRR(1:LB,34,1))+  &
                        3.D0*HRR(1:LB,55,1))+  &
                        HRR(1:LB,83,1))+  &
                        HRR(1:LB,117,1)
      !=|20,34)
      HRR(1:LB,20,34)=CDz*(CDz*(CDz*HRR(1:LB,34,1)+  &
                        3.D0*HRR(1:LB,55,1))+  &
                        3.D0*HRR(1:LB,83,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,20,1)+  &
                        3.D0*HRR(1:LB,35,1))+  &
                        3.D0*HRR(1:LB,56,1))+  &
                        HRR(1:LB,84,1))+  &
                        HRR(1:LB,119,1)
      !=|11,35)
      HRR(1:LB,11,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,11,1)+  &
                        4.D0*HRR(1:LB,26,1))+  &
                        6.D0*HRR(1:LB,47,1))+  &
                        4.D0*HRR(1:LB,75,1))+  &
                        HRR(1:LB,111,1)
      !=|12,35)
      HRR(1:LB,12,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,12,1)+  &
                        4.D0*HRR(1:LB,27,1))+  &
                        6.D0*HRR(1:LB,48,1))+  &
                        4.D0*HRR(1:LB,76,1))+  &
                        HRR(1:LB,112,1)
      !=|13,35)
      HRR(1:LB,13,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,13,1)+  &
                        4.D0*HRR(1:LB,28,1))+  &
                        6.D0*HRR(1:LB,49,1))+  &
                        4.D0*HRR(1:LB,77,1))+  &
                        HRR(1:LB,113,1)
      !=|14,35)
      HRR(1:LB,14,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,14,1)+  &
                        4.D0*HRR(1:LB,29,1))+  &
                        6.D0*HRR(1:LB,50,1))+  &
                        4.D0*HRR(1:LB,78,1))+  &
                        HRR(1:LB,114,1)
      !=|15,35)
      HRR(1:LB,15,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,15,1)+  &
                        4.D0*HRR(1:LB,30,1))+  &
                        6.D0*HRR(1:LB,51,1))+  &
                        4.D0*HRR(1:LB,79,1))+  &
                        HRR(1:LB,115,1)
      !=|16,35)
      HRR(1:LB,16,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,16,1)+  &
                        4.D0*HRR(1:LB,31,1))+  &
                        6.D0*HRR(1:LB,52,1))+  &
                        4.D0*HRR(1:LB,80,1))+  &
                        HRR(1:LB,116,1)
      !=|17,35)
      HRR(1:LB,17,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,17,1)+  &
                        4.D0*HRR(1:LB,32,1))+  &
                        6.D0*HRR(1:LB,53,1))+  &
                        4.D0*HRR(1:LB,81,1))+  &
                        HRR(1:LB,117,1)
      !=|18,35)
      HRR(1:LB,18,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,18,1)+  &
                        4.D0*HRR(1:LB,33,1))+  &
                        6.D0*HRR(1:LB,54,1))+  &
                        4.D0*HRR(1:LB,82,1))+  &
                        HRR(1:LB,118,1)
      !=|19,35)
      HRR(1:LB,19,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,19,1)+  &
                        4.D0*HRR(1:LB,34,1))+  &
                        6.D0*HRR(1:LB,55,1))+  &
                        4.D0*HRR(1:LB,83,1))+  &
                        HRR(1:LB,119,1)
      !=|20,35)
      HRR(1:LB,20,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,20,1)+  &
                        4.D0*HRR(1:LB,35,1))+  &
                        6.D0*HRR(1:LB,56,1))+  &
                        4.D0*HRR(1:LB,84,1))+  &
                        HRR(1:LB,120,1)
      !=|5,21)
      HRR(1:LB,5,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,5,1)+  &
                        4.D0*HRR(1:LB,11,1))+  &
                        6.D0*HRR(1:LB,21,1))+  &
                        4.D0*HRR(1:LB,36,1))+  &
                        HRR(1:LB,57,1)
      !=|6,21)
      HRR(1:LB,6,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,6,1)+  &
                        4.D0*HRR(1:LB,12,1))+  &
                        6.D0*HRR(1:LB,22,1))+  &
                        4.D0*HRR(1:LB,37,1))+  &
                        HRR(1:LB,58,1)
      !=|7,21)
      HRR(1:LB,7,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,7,1)+  &
                        4.D0*HRR(1:LB,13,1))+  &
                        6.D0*HRR(1:LB,23,1))+  &
                        4.D0*HRR(1:LB,38,1))+  &
                        HRR(1:LB,59,1)
      !=|8,21)
      HRR(1:LB,8,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,8,1)+  &
                        4.D0*HRR(1:LB,15,1))+  &
                        6.D0*HRR(1:LB,26,1))+  &
                        4.D0*HRR(1:LB,42,1))+  &
                        HRR(1:LB,64,1)
      !=|9,21)
      HRR(1:LB,9,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,9,1)+  &
                        4.D0*HRR(1:LB,16,1))+  &
                        6.D0*HRR(1:LB,27,1))+  &
                        4.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,65,1)
      !=|10,21)
      HRR(1:LB,10,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,10,1)+  &
                        4.D0*HRR(1:LB,18,1))+  &
                        6.D0*HRR(1:LB,30,1))+  &
                        4.D0*HRR(1:LB,47,1))+  &
                        HRR(1:LB,70,1)
      !=|5,22)
      HRR(1:LB,5,22)=CDy*HRR(1:LB,36,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,21,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,11,1)+  &
                        CDx*(CDy*HRR(1:LB,5,1)+  &
                        HRR(1:LB,12,1))+  &
                        3.D0*HRR(1:LB,22,1))+  &
                        3.D0*HRR(1:LB,37,1))+  &
                        HRR(1:LB,58,1)
      !=|6,22)
      HRR(1:LB,6,22)=CDy*HRR(1:LB,37,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,22,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,12,1)+  &
                        CDx*(CDy*HRR(1:LB,6,1)+  &
                        HRR(1:LB,13,1))+  &
                        3.D0*HRR(1:LB,23,1))+  &
                        3.D0*HRR(1:LB,38,1))+  &
                        HRR(1:LB,59,1)
      !=|7,22)
      HRR(1:LB,7,22)=CDy*HRR(1:LB,38,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,23,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,13,1)+  &
                        CDx*(CDy*HRR(1:LB,7,1)+  &
                        HRR(1:LB,14,1))+  &
                        3.D0*HRR(1:LB,24,1))+  &
                        3.D0*HRR(1:LB,39,1))+  &
                        HRR(1:LB,60,1)
      !=|8,22)
      HRR(1:LB,8,22)=CDy*HRR(1:LB,42,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,26,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,15,1)+  &
                        CDx*(CDy*HRR(1:LB,8,1)+  &
                        HRR(1:LB,16,1))+  &
                        3.D0*HRR(1:LB,27,1))+  &
                        3.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,65,1)
      !=|9,22)
      HRR(1:LB,9,22)=CDy*HRR(1:LB,43,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,27,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,16,1)+  &
                        CDx*(CDy*HRR(1:LB,9,1)+  &
                        HRR(1:LB,17,1))+  &
                        3.D0*HRR(1:LB,28,1))+  &
                        3.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,66,1)
      !=|10,22)
      HRR(1:LB,10,22)=CDy*HRR(1:LB,47,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,30,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,18,1)+  &
                        CDx*(CDy*HRR(1:LB,10,1)+  &
                        HRR(1:LB,19,1))+  &
                        3.D0*HRR(1:LB,31,1))+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,71,1)
      !=|5,23)
      HRR(1:LB,5,23)=CDy*(CDy*HRR(1:LB,21,1)+  &
                        2.D0*HRR(1:LB,37,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,11,1)+  &
                        4.D0*HRR(1:LB,22,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,5,1)+  &
                        2.D0*HRR(1:LB,12,1))+  &
                        HRR(1:LB,23,1))+  &
                        2.D0*HRR(1:LB,38,1))+  &
                        HRR(1:LB,59,1)
      !=|6,23)
      HRR(1:LB,6,23)=CDy*(CDy*HRR(1:LB,22,1)+  &
                        2.D0*HRR(1:LB,38,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,12,1)+  &
                        4.D0*HRR(1:LB,23,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,6,1)+  &
                        2.D0*HRR(1:LB,13,1))+  &
                        HRR(1:LB,24,1))+  &
                        2.D0*HRR(1:LB,39,1))+  &
                        HRR(1:LB,60,1)
      !=|7,23)
      HRR(1:LB,7,23)=CDy*(CDy*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,39,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,13,1)+  &
                        4.D0*HRR(1:LB,24,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,7,1)+  &
                        2.D0*HRR(1:LB,14,1))+  &
                        HRR(1:LB,25,1))+  &
                        2.D0*HRR(1:LB,40,1))+  &
                        HRR(1:LB,61,1)
      !=|8,23)
      HRR(1:LB,8,23)=CDy*(CDy*HRR(1:LB,26,1)+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,15,1)+  &
                        4.D0*HRR(1:LB,27,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,8,1)+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        HRR(1:LB,28,1))+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,66,1)
      !=|9,23)
      HRR(1:LB,9,23)=CDy*(CDy*HRR(1:LB,27,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,16,1)+  &
                        4.D0*HRR(1:LB,28,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,9,1)+  &
                        2.D0*HRR(1:LB,17,1))+  &
                        HRR(1:LB,29,1))+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,67,1)
      !=|10,23)
      HRR(1:LB,10,23)=CDy*(CDy*HRR(1:LB,30,1)+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,18,1)+  &
                        4.D0*HRR(1:LB,31,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,10,1)+  &
                        2.D0*HRR(1:LB,19,1))+  &
                        HRR(1:LB,32,1))+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,72,1)
      !=|5,24)
      HRR(1:LB,5,24)=CDy*(CDy*(CDy*HRR(1:LB,11,1)+  &
                        3.D0*HRR(1:LB,22,1))+  &
                        3.D0*HRR(1:LB,38,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,5,1)+  &
                        3.D0*HRR(1:LB,12,1))+  &
                        3.D0*HRR(1:LB,23,1))+  &
                        HRR(1:LB,39,1))+  &
                        HRR(1:LB,60,1)
      !=|6,24)
      HRR(1:LB,6,24)=CDy*(CDy*(CDy*HRR(1:LB,12,1)+  &
                        3.D0*HRR(1:LB,23,1))+  &
                        3.D0*HRR(1:LB,39,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,6,1)+  &
                        3.D0*HRR(1:LB,13,1))+  &
                        3.D0*HRR(1:LB,24,1))+  &
                        HRR(1:LB,40,1))+  &
                        HRR(1:LB,61,1)
      !=|7,24)
      HRR(1:LB,7,24)=CDy*(CDy*(CDy*HRR(1:LB,13,1)+  &
                        3.D0*HRR(1:LB,24,1))+  &
                        3.D0*HRR(1:LB,40,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,7,1)+  &
                        3.D0*HRR(1:LB,14,1))+  &
                        3.D0*HRR(1:LB,25,1))+  &
                        HRR(1:LB,41,1))+  &
                        HRR(1:LB,62,1)
      !=|8,24)
      HRR(1:LB,8,24)=CDy*(CDy*(CDy*HRR(1:LB,15,1)+  &
                        3.D0*HRR(1:LB,27,1))+  &
                        3.D0*HRR(1:LB,44,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,8,1)+  &
                        3.D0*HRR(1:LB,16,1))+  &
                        3.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,45,1))+  &
                        HRR(1:LB,67,1)
      !=|9,24)
      HRR(1:LB,9,24)=CDy*(CDy*(CDy*HRR(1:LB,16,1)+  &
                        3.D0*HRR(1:LB,28,1))+  &
                        3.D0*HRR(1:LB,45,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,9,1)+  &
                        3.D0*HRR(1:LB,17,1))+  &
                        3.D0*HRR(1:LB,29,1))+  &
                        HRR(1:LB,46,1))+  &
                        HRR(1:LB,68,1)
      !=|10,24)
      HRR(1:LB,10,24)=CDy*(CDy*(CDy*HRR(1:LB,18,1)+  &
                        3.D0*HRR(1:LB,31,1))+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,10,1)+  &
                        3.D0*HRR(1:LB,19,1))+  &
                        3.D0*HRR(1:LB,32,1))+  &
                        HRR(1:LB,50,1))+  &
                        HRR(1:LB,73,1)
      !=|5,25)
      HRR(1:LB,5,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,5,1)+  &
                        4.D0*HRR(1:LB,12,1))+  &
                        6.D0*HRR(1:LB,23,1))+  &
                        4.D0*HRR(1:LB,39,1))+  &
                        HRR(1:LB,61,1)
      !=|6,25)
      HRR(1:LB,6,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,6,1)+  &
                        4.D0*HRR(1:LB,13,1))+  &
                        6.D0*HRR(1:LB,24,1))+  &
                        4.D0*HRR(1:LB,40,1))+  &
                        HRR(1:LB,62,1)
      !=|7,25)
      HRR(1:LB,7,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,7,1)+  &
                        4.D0*HRR(1:LB,14,1))+  &
                        6.D0*HRR(1:LB,25,1))+  &
                        4.D0*HRR(1:LB,41,1))+  &
                        HRR(1:LB,63,1)
      !=|8,25)
      HRR(1:LB,8,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,8,1)+  &
                        4.D0*HRR(1:LB,16,1))+  &
                        6.D0*HRR(1:LB,28,1))+  &
                        4.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,68,1)
      !=|9,25)
      HRR(1:LB,9,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,9,1)+  &
                        4.D0*HRR(1:LB,17,1))+  &
                        6.D0*HRR(1:LB,29,1))+  &
                        4.D0*HRR(1:LB,46,1))+  &
                        HRR(1:LB,69,1)
      !=|10,25)
      HRR(1:LB,10,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,10,1)+  &
                        4.D0*HRR(1:LB,19,1))+  &
                        6.D0*HRR(1:LB,32,1))+  &
                        4.D0*HRR(1:LB,50,1))+  &
                        HRR(1:LB,74,1)
      !=|5,26)
      HRR(1:LB,5,26)=CDz*HRR(1:LB,36,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,21,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,11,1)+  &
                        CDx*(CDz*HRR(1:LB,5,1)+  &
                        HRR(1:LB,15,1))+  &
                        3.D0*HRR(1:LB,26,1))+  &
                        3.D0*HRR(1:LB,42,1))+  &
                        HRR(1:LB,64,1)
      !=|6,26)
      HRR(1:LB,6,26)=CDz*HRR(1:LB,37,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,22,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,12,1)+  &
                        CDx*(CDz*HRR(1:LB,6,1)+  &
                        HRR(1:LB,16,1))+  &
                        3.D0*HRR(1:LB,27,1))+  &
                        3.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,65,1)
      !=|7,26)
      HRR(1:LB,7,26)=CDz*HRR(1:LB,38,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,23,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,13,1)+  &
                        CDx*(CDz*HRR(1:LB,7,1)+  &
                        HRR(1:LB,17,1))+  &
                        3.D0*HRR(1:LB,28,1))+  &
                        3.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,66,1)
      !=|8,26)
      HRR(1:LB,8,26)=CDz*HRR(1:LB,42,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,26,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,15,1)+  &
                        CDx*(CDz*HRR(1:LB,8,1)+  &
                        HRR(1:LB,18,1))+  &
                        3.D0*HRR(1:LB,30,1))+  &
                        3.D0*HRR(1:LB,47,1))+  &
                        HRR(1:LB,70,1)
      !=|9,26)
      HRR(1:LB,9,26)=CDz*HRR(1:LB,43,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,27,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,16,1)+  &
                        CDx*(CDz*HRR(1:LB,9,1)+  &
                        HRR(1:LB,19,1))+  &
                        3.D0*HRR(1:LB,31,1))+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,71,1)
      !=|10,26)
      HRR(1:LB,10,26)=CDz*HRR(1:LB,47,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,30,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,18,1)+  &
                        CDx*(CDz*HRR(1:LB,10,1)+  &
                        HRR(1:LB,20,1))+  &
                        3.D0*HRR(1:LB,33,1))+  &
                        3.D0*HRR(1:LB,51,1))+  &
                        HRR(1:LB,75,1)
      !=|5,27)
      HRR(1:LB,5,27)=CDz*HRR(1:LB,37,1)+  &
                        CDy*(CDz*HRR(1:LB,21,1)+  &
                        HRR(1:LB,42,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,22,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,11,1)+  &
                        2.D0*HRR(1:LB,26,1))+  &
                        CDx*(CDz*HRR(1:LB,12,1)+  &
                        CDy*(CDz*HRR(1:LB,5,1)+  &
                        HRR(1:LB,15,1))+  &
                        HRR(1:LB,27,1))+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,65,1)
      !=|6,27)
      HRR(1:LB,6,27)=CDz*HRR(1:LB,38,1)+  &
                        CDy*(CDz*HRR(1:LB,22,1)+  &
                        HRR(1:LB,43,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,23,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,12,1)+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        CDx*(CDz*HRR(1:LB,13,1)+  &
                        CDy*(CDz*HRR(1:LB,6,1)+  &
                        HRR(1:LB,16,1))+  &
                        HRR(1:LB,28,1))+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,66,1)
      !=|7,27)
      HRR(1:LB,7,27)=CDz*HRR(1:LB,39,1)+  &
                        CDy*(CDz*HRR(1:LB,23,1)+  &
                        HRR(1:LB,44,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,24,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,13,1)+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        CDx*(CDz*HRR(1:LB,14,1)+  &
                        CDy*(CDz*HRR(1:LB,7,1)+  &
                        HRR(1:LB,17,1))+  &
                        HRR(1:LB,29,1))+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,67,1)
      !=|8,27)
      HRR(1:LB,8,27)=CDz*HRR(1:LB,43,1)+  &
                        CDy*(CDz*HRR(1:LB,26,1)+  &
                        HRR(1:LB,47,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,27,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,15,1)+  &
                        2.D0*HRR(1:LB,30,1))+  &
                        CDx*(CDz*HRR(1:LB,16,1)+  &
                        CDy*(CDz*HRR(1:LB,8,1)+  &
                        HRR(1:LB,18,1))+  &
                        HRR(1:LB,31,1))+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,71,1)
      !=|9,27)
      HRR(1:LB,9,27)=CDz*HRR(1:LB,44,1)+  &
                        CDy*(CDz*HRR(1:LB,27,1)+  &
                        HRR(1:LB,48,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,28,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,16,1)+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        CDx*(CDz*HRR(1:LB,17,1)+  &
                        CDy*(CDz*HRR(1:LB,9,1)+  &
                        HRR(1:LB,19,1))+  &
                        HRR(1:LB,32,1))+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,72,1)
      !=|10,27)
      HRR(1:LB,10,27)=CDz*HRR(1:LB,48,1)+  &
                        CDy*(CDz*HRR(1:LB,30,1)+  &
                        HRR(1:LB,51,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,31,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,18,1)+  &
                        2.D0*HRR(1:LB,33,1))+  &
                        CDx*(CDz*HRR(1:LB,19,1)+  &
                        CDy*(CDz*HRR(1:LB,10,1)+  &
                        HRR(1:LB,20,1))+  &
                        HRR(1:LB,34,1))+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,76,1)
      !=|5,28)
      HRR(1:LB,5,28)=CDz*HRR(1:LB,38,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,22,1)+  &
                        CDy*(CDz*HRR(1:LB,11,1)+  &
                        HRR(1:LB,26,1))+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        CDx*(CDz*HRR(1:LB,23,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,12,1)+  &
                        CDy*(CDz*HRR(1:LB,5,1)+  &
                        HRR(1:LB,15,1))+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        HRR(1:LB,44,1))+  &
                        HRR(1:LB,66,1)
      !=|6,28)
      HRR(1:LB,6,28)=CDz*HRR(1:LB,39,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,23,1)+  &
                        CDy*(CDz*HRR(1:LB,12,1)+  &
                        HRR(1:LB,27,1))+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        CDx*(CDz*HRR(1:LB,24,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,13,1)+  &
                        CDy*(CDz*HRR(1:LB,6,1)+  &
                        HRR(1:LB,16,1))+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,45,1))+  &
                        HRR(1:LB,67,1)
      !=|7,28)
      HRR(1:LB,7,28)=CDz*HRR(1:LB,40,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,24,1)+  &
                        CDy*(CDz*HRR(1:LB,13,1)+  &
                        HRR(1:LB,28,1))+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        CDx*(CDz*HRR(1:LB,25,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,14,1)+  &
                        CDy*(CDz*HRR(1:LB,7,1)+  &
                        HRR(1:LB,17,1))+  &
                        2.D0*HRR(1:LB,29,1))+  &
                        HRR(1:LB,46,1))+  &
                        HRR(1:LB,68,1)
      !=|8,28)
      HRR(1:LB,8,28)=CDz*HRR(1:LB,44,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,27,1)+  &
                        CDy*(CDz*HRR(1:LB,15,1)+  &
                        HRR(1:LB,30,1))+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        CDx*(CDz*HRR(1:LB,28,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,16,1)+  &
                        CDy*(CDz*HRR(1:LB,8,1)+  &
                        HRR(1:LB,18,1))+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,49,1))+  &
                        HRR(1:LB,72,1)
      !=|9,28)
      HRR(1:LB,9,28)=CDz*HRR(1:LB,45,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,28,1)+  &
                        CDy*(CDz*HRR(1:LB,16,1)+  &
                        HRR(1:LB,31,1))+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        CDx*(CDz*HRR(1:LB,29,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,17,1)+  &
                        CDy*(CDz*HRR(1:LB,9,1)+  &
                        HRR(1:LB,19,1))+  &
                        2.D0*HRR(1:LB,32,1))+  &
                        HRR(1:LB,50,1))+  &
                        HRR(1:LB,73,1)
      !=|10,28)
      HRR(1:LB,10,28)=CDz*HRR(1:LB,49,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,31,1)+  &
                        CDy*(CDz*HRR(1:LB,18,1)+  &
                        HRR(1:LB,33,1))+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        CDx*(CDz*HRR(1:LB,32,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,19,1)+  &
                        CDy*(CDz*HRR(1:LB,10,1)+  &
                        HRR(1:LB,20,1))+  &
                        2.D0*HRR(1:LB,34,1))+  &
                        HRR(1:LB,53,1))+  &
                        HRR(1:LB,77,1)
      !=|5,29)
      HRR(1:LB,5,29)=CDz*HRR(1:LB,39,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,23,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,12,1)+  &
                        CDy*(CDz*HRR(1:LB,5,1)+  &
                        HRR(1:LB,15,1))+  &
                        3.D0*HRR(1:LB,27,1))+  &
                        3.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,67,1)
      !=|6,29)
      HRR(1:LB,6,29)=CDz*HRR(1:LB,40,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,24,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,13,1)+  &
                        CDy*(CDz*HRR(1:LB,6,1)+  &
                        HRR(1:LB,16,1))+  &
                        3.D0*HRR(1:LB,28,1))+  &
                        3.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,68,1)
      !=|7,29)
      HRR(1:LB,7,29)=CDz*HRR(1:LB,41,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,25,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,14,1)+  &
                        CDy*(CDz*HRR(1:LB,7,1)+  &
                        HRR(1:LB,17,1))+  &
                        3.D0*HRR(1:LB,29,1))+  &
                        3.D0*HRR(1:LB,46,1))+  &
                        HRR(1:LB,69,1)
      !=|8,29)
      HRR(1:LB,8,29)=CDz*HRR(1:LB,45,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,28,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,16,1)+  &
                        CDy*(CDz*HRR(1:LB,8,1)+  &
                        HRR(1:LB,18,1))+  &
                        3.D0*HRR(1:LB,31,1))+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,73,1)
      !=|9,29)
      HRR(1:LB,9,29)=CDz*HRR(1:LB,46,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,29,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,17,1)+  &
                        CDy*(CDz*HRR(1:LB,9,1)+  &
                        HRR(1:LB,19,1))+  &
                        3.D0*HRR(1:LB,32,1))+  &
                        3.D0*HRR(1:LB,50,1))+  &
                        HRR(1:LB,74,1)
      !=|10,29)
      HRR(1:LB,10,29)=CDz*HRR(1:LB,50,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,32,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,19,1)+  &
                        CDy*(CDz*HRR(1:LB,10,1)+  &
                        HRR(1:LB,20,1))+  &
                        3.D0*HRR(1:LB,34,1))+  &
                        3.D0*HRR(1:LB,53,1))+  &
                        HRR(1:LB,78,1)
      !=|5,30)
      HRR(1:LB,5,30)=CDz*(CDz*HRR(1:LB,21,1)+  &
                        2.D0*HRR(1:LB,42,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,11,1)+  &
                        4.D0*HRR(1:LB,26,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,5,1)+  &
                        2.D0*HRR(1:LB,15,1))+  &
                        HRR(1:LB,30,1))+  &
                        2.D0*HRR(1:LB,47,1))+  &
                        HRR(1:LB,70,1)
      !=|6,30)
      HRR(1:LB,6,30)=CDz*(CDz*HRR(1:LB,22,1)+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,12,1)+  &
                        4.D0*HRR(1:LB,27,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,6,1)+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        HRR(1:LB,31,1))+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,71,1)
      !=|7,30)
      HRR(1:LB,7,30)=CDz*(CDz*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,13,1)+  &
                        4.D0*HRR(1:LB,28,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,7,1)+  &
                        2.D0*HRR(1:LB,17,1))+  &
                        HRR(1:LB,32,1))+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,72,1)
      !=|8,30)
      HRR(1:LB,8,30)=CDz*(CDz*HRR(1:LB,26,1)+  &
                        2.D0*HRR(1:LB,47,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,15,1)+  &
                        4.D0*HRR(1:LB,30,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,8,1)+  &
                        2.D0*HRR(1:LB,18,1))+  &
                        HRR(1:LB,33,1))+  &
                        2.D0*HRR(1:LB,51,1))+  &
                        HRR(1:LB,75,1)
      !=|9,30)
      HRR(1:LB,9,30)=CDz*(CDz*HRR(1:LB,27,1)+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,16,1)+  &
                        4.D0*HRR(1:LB,31,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,9,1)+  &
                        2.D0*HRR(1:LB,19,1))+  &
                        HRR(1:LB,34,1))+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,76,1)
      !=|10,30)
      HRR(1:LB,10,30)=CDz*(CDz*HRR(1:LB,30,1)+  &
                        2.D0*HRR(1:LB,51,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,18,1)+  &
                        4.D0*HRR(1:LB,33,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,10,1)+  &
                        2.D0*HRR(1:LB,20,1))+  &
                        HRR(1:LB,35,1))+  &
                        2.D0*HRR(1:LB,54,1))+  &
                        HRR(1:LB,79,1)
      !=|5,31)
      HRR(1:LB,5,31)=CDz*(CDz*HRR(1:LB,22,1)+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,11,1)+  &
                        2.D0*HRR(1:LB,26,1))+  &
                        HRR(1:LB,47,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,12,1)+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,5,1)+  &
                        2.D0*HRR(1:LB,15,1))+  &
                        HRR(1:LB,30,1))+  &
                        HRR(1:LB,48,1))+  &
                        HRR(1:LB,71,1)
      !=|6,31)
      HRR(1:LB,6,31)=CDz*(CDz*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,12,1)+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        HRR(1:LB,48,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,13,1)+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,6,1)+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        HRR(1:LB,31,1))+  &
                        HRR(1:LB,49,1))+  &
                        HRR(1:LB,72,1)
      !=|7,31)
      HRR(1:LB,7,31)=CDz*(CDz*HRR(1:LB,24,1)+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,13,1)+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,49,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,14,1)+  &
                        2.D0*HRR(1:LB,29,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,7,1)+  &
                        2.D0*HRR(1:LB,17,1))+  &
                        HRR(1:LB,32,1))+  &
                        HRR(1:LB,50,1))+  &
                        HRR(1:LB,73,1)
      !=|8,31)
      HRR(1:LB,8,31)=CDz*(CDz*HRR(1:LB,27,1)+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,15,1)+  &
                        2.D0*HRR(1:LB,30,1))+  &
                        HRR(1:LB,51,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,16,1)+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,8,1)+  &
                        2.D0*HRR(1:LB,18,1))+  &
                        HRR(1:LB,33,1))+  &
                        HRR(1:LB,52,1))+  &
                        HRR(1:LB,76,1)
      !=|9,31)
      HRR(1:LB,9,31)=CDz*(CDz*HRR(1:LB,28,1)+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,16,1)+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,52,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,17,1)+  &
                        2.D0*HRR(1:LB,32,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,9,1)+  &
                        2.D0*HRR(1:LB,19,1))+  &
                        HRR(1:LB,34,1))+  &
                        HRR(1:LB,53,1))+  &
                        HRR(1:LB,77,1)
      !=|10,31)
      HRR(1:LB,10,31)=CDz*(CDz*HRR(1:LB,31,1)+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,18,1)+  &
                        2.D0*HRR(1:LB,33,1))+  &
                        HRR(1:LB,54,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,19,1)+  &
                        2.D0*HRR(1:LB,34,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,10,1)+  &
                        2.D0*HRR(1:LB,20,1))+  &
                        HRR(1:LB,35,1))+  &
                        HRR(1:LB,55,1))+  &
                        HRR(1:LB,80,1)
      !=|5,32)
      HRR(1:LB,5,32)=CDz*(CDz*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,12,1)+  &
                        4.D0*HRR(1:LB,27,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,5,1)+  &
                        2.D0*HRR(1:LB,15,1))+  &
                        HRR(1:LB,30,1))+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,72,1)
      !=|6,32)
      HRR(1:LB,6,32)=CDz*(CDz*HRR(1:LB,24,1)+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,13,1)+  &
                        4.D0*HRR(1:LB,28,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,6,1)+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        HRR(1:LB,31,1))+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,73,1)
      !=|7,32)
      HRR(1:LB,7,32)=CDz*(CDz*HRR(1:LB,25,1)+  &
                        2.D0*HRR(1:LB,46,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,14,1)+  &
                        4.D0*HRR(1:LB,29,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,7,1)+  &
                        2.D0*HRR(1:LB,17,1))+  &
                        HRR(1:LB,32,1))+  &
                        2.D0*HRR(1:LB,50,1))+  &
                        HRR(1:LB,74,1)
      !=|8,32)
      HRR(1:LB,8,32)=CDz*(CDz*HRR(1:LB,28,1)+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,16,1)+  &
                        4.D0*HRR(1:LB,31,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,8,1)+  &
                        2.D0*HRR(1:LB,18,1))+  &
                        HRR(1:LB,33,1))+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,77,1)
      !=|9,32)
      HRR(1:LB,9,32)=CDz*(CDz*HRR(1:LB,29,1)+  &
                        2.D0*HRR(1:LB,50,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,17,1)+  &
                        4.D0*HRR(1:LB,32,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,9,1)+  &
                        2.D0*HRR(1:LB,19,1))+  &
                        HRR(1:LB,34,1))+  &
                        2.D0*HRR(1:LB,53,1))+  &
                        HRR(1:LB,78,1)
      !=|10,32)
      HRR(1:LB,10,32)=CDz*(CDz*HRR(1:LB,32,1)+  &
                        2.D0*HRR(1:LB,53,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,19,1)+  &
                        4.D0*HRR(1:LB,34,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,10,1)+  &
                        2.D0*HRR(1:LB,20,1))+  &
                        HRR(1:LB,35,1))+  &
                        2.D0*HRR(1:LB,55,1))+  &
                        HRR(1:LB,81,1)
      !=|5,33)
      HRR(1:LB,5,33)=CDz*(CDz*(CDz*HRR(1:LB,11,1)+  &
                        3.D0*HRR(1:LB,26,1))+  &
                        3.D0*HRR(1:LB,47,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,5,1)+  &
                        3.D0*HRR(1:LB,15,1))+  &
                        3.D0*HRR(1:LB,30,1))+  &
                        HRR(1:LB,51,1))+  &
                        HRR(1:LB,75,1)
      !=|6,33)
      HRR(1:LB,6,33)=CDz*(CDz*(CDz*HRR(1:LB,12,1)+  &
                        3.D0*HRR(1:LB,27,1))+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,6,1)+  &
                        3.D0*HRR(1:LB,16,1))+  &
                        3.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,52,1))+  &
                        HRR(1:LB,76,1)
      !=|7,33)
      HRR(1:LB,7,33)=CDz*(CDz*(CDz*HRR(1:LB,13,1)+  &
                        3.D0*HRR(1:LB,28,1))+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,7,1)+  &
                        3.D0*HRR(1:LB,17,1))+  &
                        3.D0*HRR(1:LB,32,1))+  &
                        HRR(1:LB,53,1))+  &
                        HRR(1:LB,77,1)
      !=|8,33)
      HRR(1:LB,8,33)=CDz*(CDz*(CDz*HRR(1:LB,15,1)+  &
                        3.D0*HRR(1:LB,30,1))+  &
                        3.D0*HRR(1:LB,51,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,8,1)+  &
                        3.D0*HRR(1:LB,18,1))+  &
                        3.D0*HRR(1:LB,33,1))+  &
                        HRR(1:LB,54,1))+  &
                        HRR(1:LB,79,1)
      !=|9,33)
      HRR(1:LB,9,33)=CDz*(CDz*(CDz*HRR(1:LB,16,1)+  &
                        3.D0*HRR(1:LB,31,1))+  &
                        3.D0*HRR(1:LB,52,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,9,1)+  &
                        3.D0*HRR(1:LB,19,1))+  &
                        3.D0*HRR(1:LB,34,1))+  &
                        HRR(1:LB,55,1))+  &
                        HRR(1:LB,80,1)
      !=|10,33)
      HRR(1:LB,10,33)=CDz*(CDz*(CDz*HRR(1:LB,18,1)+  &
                        3.D0*HRR(1:LB,33,1))+  &
                        3.D0*HRR(1:LB,54,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,10,1)+  &
                        3.D0*HRR(1:LB,20,1))+  &
                        3.D0*HRR(1:LB,35,1))+  &
                        HRR(1:LB,56,1))+  &
                        HRR(1:LB,82,1)
      !=|5,34)
      HRR(1:LB,5,34)=CDz*(CDz*(CDz*HRR(1:LB,12,1)+  &
                        3.D0*HRR(1:LB,27,1))+  &
                        3.D0*HRR(1:LB,48,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,5,1)+  &
                        3.D0*HRR(1:LB,15,1))+  &
                        3.D0*HRR(1:LB,30,1))+  &
                        HRR(1:LB,51,1))+  &
                        HRR(1:LB,76,1)
      !=|6,34)
      HRR(1:LB,6,34)=CDz*(CDz*(CDz*HRR(1:LB,13,1)+  &
                        3.D0*HRR(1:LB,28,1))+  &
                        3.D0*HRR(1:LB,49,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,6,1)+  &
                        3.D0*HRR(1:LB,16,1))+  &
                        3.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,52,1))+  &
                        HRR(1:LB,77,1)
      !=|7,34)
      HRR(1:LB,7,34)=CDz*(CDz*(CDz*HRR(1:LB,14,1)+  &
                        3.D0*HRR(1:LB,29,1))+  &
                        3.D0*HRR(1:LB,50,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,7,1)+  &
                        3.D0*HRR(1:LB,17,1))+  &
                        3.D0*HRR(1:LB,32,1))+  &
                        HRR(1:LB,53,1))+  &
                        HRR(1:LB,78,1)
      !=|8,34)
      HRR(1:LB,8,34)=CDz*(CDz*(CDz*HRR(1:LB,16,1)+  &
                        3.D0*HRR(1:LB,31,1))+  &
                        3.D0*HRR(1:LB,52,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,8,1)+  &
                        3.D0*HRR(1:LB,18,1))+  &
                        3.D0*HRR(1:LB,33,1))+  &
                        HRR(1:LB,54,1))+  &
                        HRR(1:LB,80,1)
      !=|9,34)
      HRR(1:LB,9,34)=CDz*(CDz*(CDz*HRR(1:LB,17,1)+  &
                        3.D0*HRR(1:LB,32,1))+  &
                        3.D0*HRR(1:LB,53,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,9,1)+  &
                        3.D0*HRR(1:LB,19,1))+  &
                        3.D0*HRR(1:LB,34,1))+  &
                        HRR(1:LB,55,1))+  &
                        HRR(1:LB,81,1)
      !=|10,34)
      HRR(1:LB,10,34)=CDz*(CDz*(CDz*HRR(1:LB,19,1)+  &
                        3.D0*HRR(1:LB,34,1))+  &
                        3.D0*HRR(1:LB,55,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,10,1)+  &
                        3.D0*HRR(1:LB,20,1))+  &
                        3.D0*HRR(1:LB,35,1))+  &
                        HRR(1:LB,56,1))+  &
                        HRR(1:LB,83,1)
      !=|5,35)
      HRR(1:LB,5,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,5,1)+  &
                        4.D0*HRR(1:LB,15,1))+  &
                        6.D0*HRR(1:LB,30,1))+  &
                        4.D0*HRR(1:LB,51,1))+  &
                        HRR(1:LB,79,1)
      !=|6,35)
      HRR(1:LB,6,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,6,1)+  &
                        4.D0*HRR(1:LB,16,1))+  &
                        6.D0*HRR(1:LB,31,1))+  &
                        4.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,80,1)
      !=|7,35)
      HRR(1:LB,7,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,7,1)+  &
                        4.D0*HRR(1:LB,17,1))+  &
                        6.D0*HRR(1:LB,32,1))+  &
                        4.D0*HRR(1:LB,53,1))+  &
                        HRR(1:LB,81,1)
      !=|8,35)
      HRR(1:LB,8,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,8,1)+  &
                        4.D0*HRR(1:LB,18,1))+  &
                        6.D0*HRR(1:LB,33,1))+  &
                        4.D0*HRR(1:LB,54,1))+  &
                        HRR(1:LB,82,1)
      !=|9,35)
      HRR(1:LB,9,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,9,1)+  &
                        4.D0*HRR(1:LB,19,1))+  &
                        6.D0*HRR(1:LB,34,1))+  &
                        4.D0*HRR(1:LB,55,1))+  &
                        HRR(1:LB,83,1)
      !=|10,35)
      HRR(1:LB,10,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,10,1)+  &
                        4.D0*HRR(1:LB,20,1))+  &
                        6.D0*HRR(1:LB,35,1))+  &
                        4.D0*HRR(1:LB,56,1))+  &
                        HRR(1:LB,84,1)
      !=|2,21)
      HRR(1:LB,2,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,2,1)+  &
                        4.D0*HRR(1:LB,5,1))+  &
                        6.D0*HRR(1:LB,11,1))+  &
                        4.D0*HRR(1:LB,21,1))+  &
                        HRR(1:LB,36,1)
      !=|3,21)
      HRR(1:LB,3,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,3,1)+  &
                        4.D0*HRR(1:LB,6,1))+  &
                        6.D0*HRR(1:LB,12,1))+  &
                        4.D0*HRR(1:LB,22,1))+  &
                        HRR(1:LB,37,1)
      !=|4,21)
      HRR(1:LB,4,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,4,1)+  &
                        4.D0*HRR(1:LB,8,1))+  &
                        6.D0*HRR(1:LB,15,1))+  &
                        4.D0*HRR(1:LB,26,1))+  &
                        HRR(1:LB,42,1)
      !=|2,22)
      HRR(1:LB,2,22)=CDy*HRR(1:LB,21,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,11,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,5,1)+  &
                        CDx*(CDy*HRR(1:LB,2,1)+  &
                        HRR(1:LB,6,1))+  &
                        3.D0*HRR(1:LB,12,1))+  &
                        3.D0*HRR(1:LB,22,1))+  &
                        HRR(1:LB,37,1)
      !=|3,22)
      HRR(1:LB,3,22)=CDy*HRR(1:LB,22,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,12,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,6,1)+  &
                        CDx*(CDy*HRR(1:LB,3,1)+  &
                        HRR(1:LB,7,1))+  &
                        3.D0*HRR(1:LB,13,1))+  &
                        3.D0*HRR(1:LB,23,1))+  &
                        HRR(1:LB,38,1)
      !=|4,22)
      HRR(1:LB,4,22)=CDy*HRR(1:LB,26,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,15,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,8,1)+  &
                        CDx*(CDy*HRR(1:LB,4,1)+  &
                        HRR(1:LB,9,1))+  &
                        3.D0*HRR(1:LB,16,1))+  &
                        3.D0*HRR(1:LB,27,1))+  &
                        HRR(1:LB,43,1)
      !=|2,23)
      HRR(1:LB,2,23)=CDy*(CDy*HRR(1:LB,11,1)+  &
                        2.D0*HRR(1:LB,22,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,5,1)+  &
                        4.D0*HRR(1:LB,12,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,2,1)+  &
                        2.D0*HRR(1:LB,6,1))+  &
                        HRR(1:LB,13,1))+  &
                        2.D0*HRR(1:LB,23,1))+  &
                        HRR(1:LB,38,1)
      !=|3,23)
      HRR(1:LB,3,23)=CDy*(CDy*HRR(1:LB,12,1)+  &
                        2.D0*HRR(1:LB,23,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,6,1)+  &
                        4.D0*HRR(1:LB,13,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,3,1)+  &
                        2.D0*HRR(1:LB,7,1))+  &
                        HRR(1:LB,14,1))+  &
                        2.D0*HRR(1:LB,24,1))+  &
                        HRR(1:LB,39,1)
      !=|4,23)
      HRR(1:LB,4,23)=CDy*(CDy*HRR(1:LB,15,1)+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,8,1)+  &
                        4.D0*HRR(1:LB,16,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,4,1)+  &
                        2.D0*HRR(1:LB,9,1))+  &
                        HRR(1:LB,17,1))+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,44,1)
      !=|2,24)
      HRR(1:LB,2,24)=CDy*(CDy*(CDy*HRR(1:LB,5,1)+  &
                        3.D0*HRR(1:LB,12,1))+  &
                        3.D0*HRR(1:LB,23,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,2,1)+  &
                        3.D0*HRR(1:LB,6,1))+  &
                        3.D0*HRR(1:LB,13,1))+  &
                        HRR(1:LB,24,1))+  &
                        HRR(1:LB,39,1)
      !=|3,24)
      HRR(1:LB,3,24)=CDy*(CDy*(CDy*HRR(1:LB,6,1)+  &
                        3.D0*HRR(1:LB,13,1))+  &
                        3.D0*HRR(1:LB,24,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,3,1)+  &
                        3.D0*HRR(1:LB,7,1))+  &
                        3.D0*HRR(1:LB,14,1))+  &
                        HRR(1:LB,25,1))+  &
                        HRR(1:LB,40,1)
      !=|4,24)
      HRR(1:LB,4,24)=CDy*(CDy*(CDy*HRR(1:LB,8,1)+  &
                        3.D0*HRR(1:LB,16,1))+  &
                        3.D0*HRR(1:LB,28,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,4,1)+  &
                        3.D0*HRR(1:LB,9,1))+  &
                        3.D0*HRR(1:LB,17,1))+  &
                        HRR(1:LB,29,1))+  &
                        HRR(1:LB,45,1)
      !=|2,25)
      HRR(1:LB,2,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,2,1)+  &
                        4.D0*HRR(1:LB,6,1))+  &
                        6.D0*HRR(1:LB,13,1))+  &
                        4.D0*HRR(1:LB,24,1))+  &
                        HRR(1:LB,40,1)
      !=|3,25)
      HRR(1:LB,3,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,3,1)+  &
                        4.D0*HRR(1:LB,7,1))+  &
                        6.D0*HRR(1:LB,14,1))+  &
                        4.D0*HRR(1:LB,25,1))+  &
                        HRR(1:LB,41,1)
      !=|4,25)
      HRR(1:LB,4,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,4,1)+  &
                        4.D0*HRR(1:LB,9,1))+  &
                        6.D0*HRR(1:LB,17,1))+  &
                        4.D0*HRR(1:LB,29,1))+  &
                        HRR(1:LB,46,1)
      !=|2,26)
      HRR(1:LB,2,26)=CDz*HRR(1:LB,21,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,11,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,5,1)+  &
                        CDx*(CDz*HRR(1:LB,2,1)+  &
                        HRR(1:LB,8,1))+  &
                        3.D0*HRR(1:LB,15,1))+  &
                        3.D0*HRR(1:LB,26,1))+  &
                        HRR(1:LB,42,1)
      !=|3,26)
      HRR(1:LB,3,26)=CDz*HRR(1:LB,22,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,12,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,6,1)+  &
                        CDx*(CDz*HRR(1:LB,3,1)+  &
                        HRR(1:LB,9,1))+  &
                        3.D0*HRR(1:LB,16,1))+  &
                        3.D0*HRR(1:LB,27,1))+  &
                        HRR(1:LB,43,1)
      !=|4,26)
      HRR(1:LB,4,26)=CDz*HRR(1:LB,26,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,15,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,8,1)+  &
                        CDx*(CDz*HRR(1:LB,4,1)+  &
                        HRR(1:LB,10,1))+  &
                        3.D0*HRR(1:LB,18,1))+  &
                        3.D0*HRR(1:LB,30,1))+  &
                        HRR(1:LB,47,1)
      !=|2,27)
      HRR(1:LB,2,27)=CDz*HRR(1:LB,22,1)+  &
                        CDy*(CDz*HRR(1:LB,11,1)+  &
                        HRR(1:LB,26,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,12,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,5,1)+  &
                        2.D0*HRR(1:LB,15,1))+  &
                        CDx*(CDz*HRR(1:LB,6,1)+  &
                        CDy*(CDz*HRR(1:LB,2,1)+  &
                        HRR(1:LB,8,1))+  &
                        HRR(1:LB,16,1))+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        HRR(1:LB,43,1)
      !=|3,27)
      HRR(1:LB,3,27)=CDz*HRR(1:LB,23,1)+  &
                        CDy*(CDz*HRR(1:LB,12,1)+  &
                        HRR(1:LB,27,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,13,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,6,1)+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        CDx*(CDz*HRR(1:LB,7,1)+  &
                        CDy*(CDz*HRR(1:LB,3,1)+  &
                        HRR(1:LB,9,1))+  &
                        HRR(1:LB,17,1))+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,44,1)
      !=|4,27)
      HRR(1:LB,4,27)=CDz*HRR(1:LB,27,1)+  &
                        CDy*(CDz*HRR(1:LB,15,1)+  &
                        HRR(1:LB,30,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,16,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,8,1)+  &
                        2.D0*HRR(1:LB,18,1))+  &
                        CDx*(CDz*HRR(1:LB,9,1)+  &
                        CDy*(CDz*HRR(1:LB,4,1)+  &
                        HRR(1:LB,10,1))+  &
                        HRR(1:LB,19,1))+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,48,1)
      !=|2,28)
      HRR(1:LB,2,28)=CDz*HRR(1:LB,23,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,12,1)+  &
                        CDy*(CDz*HRR(1:LB,5,1)+  &
                        HRR(1:LB,15,1))+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        CDx*(CDz*HRR(1:LB,13,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,6,1)+  &
                        CDy*(CDz*HRR(1:LB,2,1)+  &
                        HRR(1:LB,8,1))+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        HRR(1:LB,28,1))+  &
                        HRR(1:LB,44,1)
      !=|3,28)
      HRR(1:LB,3,28)=CDz*HRR(1:LB,24,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,13,1)+  &
                        CDy*(CDz*HRR(1:LB,6,1)+  &
                        HRR(1:LB,16,1))+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        CDx*(CDz*HRR(1:LB,14,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,7,1)+  &
                        CDy*(CDz*HRR(1:LB,3,1)+  &
                        HRR(1:LB,9,1))+  &
                        2.D0*HRR(1:LB,17,1))+  &
                        HRR(1:LB,29,1))+  &
                        HRR(1:LB,45,1)
      !=|4,28)
      HRR(1:LB,4,28)=CDz*HRR(1:LB,28,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,16,1)+  &
                        CDy*(CDz*HRR(1:LB,8,1)+  &
                        HRR(1:LB,18,1))+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        CDx*(CDz*HRR(1:LB,17,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,9,1)+  &
                        CDy*(CDz*HRR(1:LB,4,1)+  &
                        HRR(1:LB,10,1))+  &
                        2.D0*HRR(1:LB,19,1))+  &
                        HRR(1:LB,32,1))+  &
                        HRR(1:LB,49,1)
      !=|2,29)
      HRR(1:LB,2,29)=CDz*HRR(1:LB,24,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,13,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,6,1)+  &
                        CDy*(CDz*HRR(1:LB,2,1)+  &
                        HRR(1:LB,8,1))+  &
                        3.D0*HRR(1:LB,16,1))+  &
                        3.D0*HRR(1:LB,28,1))+  &
                        HRR(1:LB,45,1)
      !=|3,29)
      HRR(1:LB,3,29)=CDz*HRR(1:LB,25,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,14,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,7,1)+  &
                        CDy*(CDz*HRR(1:LB,3,1)+  &
                        HRR(1:LB,9,1))+  &
                        3.D0*HRR(1:LB,17,1))+  &
                        3.D0*HRR(1:LB,29,1))+  &
                        HRR(1:LB,46,1)
      !=|4,29)
      HRR(1:LB,4,29)=CDz*HRR(1:LB,29,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,17,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,9,1)+  &
                        CDy*(CDz*HRR(1:LB,4,1)+  &
                        HRR(1:LB,10,1))+  &
                        3.D0*HRR(1:LB,19,1))+  &
                        3.D0*HRR(1:LB,32,1))+  &
                        HRR(1:LB,50,1)
      !=|2,30)
      HRR(1:LB,2,30)=CDz*(CDz*HRR(1:LB,11,1)+  &
                        2.D0*HRR(1:LB,26,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,5,1)+  &
                        4.D0*HRR(1:LB,15,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,2,1)+  &
                        2.D0*HRR(1:LB,8,1))+  &
                        HRR(1:LB,18,1))+  &
                        2.D0*HRR(1:LB,30,1))+  &
                        HRR(1:LB,47,1)
      !=|3,30)
      HRR(1:LB,3,30)=CDz*(CDz*HRR(1:LB,12,1)+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,6,1)+  &
                        4.D0*HRR(1:LB,16,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,3,1)+  &
                        2.D0*HRR(1:LB,9,1))+  &
                        HRR(1:LB,19,1))+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,48,1)
      !=|4,30)
      HRR(1:LB,4,30)=CDz*(CDz*HRR(1:LB,15,1)+  &
                        2.D0*HRR(1:LB,30,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,8,1)+  &
                        4.D0*HRR(1:LB,18,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,4,1)+  &
                        2.D0*HRR(1:LB,10,1))+  &
                        HRR(1:LB,20,1))+  &
                        2.D0*HRR(1:LB,33,1))+  &
                        HRR(1:LB,51,1)
      !=|2,31)
      HRR(1:LB,2,31)=CDz*(CDz*HRR(1:LB,12,1)+  &
                        2.D0*HRR(1:LB,27,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,5,1)+  &
                        2.D0*HRR(1:LB,15,1))+  &
                        HRR(1:LB,30,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,6,1)+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,2,1)+  &
                        2.D0*HRR(1:LB,8,1))+  &
                        HRR(1:LB,18,1))+  &
                        HRR(1:LB,31,1))+  &
                        HRR(1:LB,48,1)
      !=|3,31)
      HRR(1:LB,3,31)=CDz*(CDz*HRR(1:LB,13,1)+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,6,1)+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        HRR(1:LB,31,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,7,1)+  &
                        2.D0*HRR(1:LB,17,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,3,1)+  &
                        2.D0*HRR(1:LB,9,1))+  &
                        HRR(1:LB,19,1))+  &
                        HRR(1:LB,32,1))+  &
                        HRR(1:LB,49,1)
      !=|4,31)
      HRR(1:LB,4,31)=CDz*(CDz*HRR(1:LB,16,1)+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,8,1)+  &
                        2.D0*HRR(1:LB,18,1))+  &
                        HRR(1:LB,33,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,9,1)+  &
                        2.D0*HRR(1:LB,19,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,4,1)+  &
                        2.D0*HRR(1:LB,10,1))+  &
                        HRR(1:LB,20,1))+  &
                        HRR(1:LB,34,1))+  &
                        HRR(1:LB,52,1)
      !=|2,32)
      HRR(1:LB,2,32)=CDz*(CDz*HRR(1:LB,13,1)+  &
                        2.D0*HRR(1:LB,28,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,6,1)+  &
                        4.D0*HRR(1:LB,16,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,2,1)+  &
                        2.D0*HRR(1:LB,8,1))+  &
                        HRR(1:LB,18,1))+  &
                        2.D0*HRR(1:LB,31,1))+  &
                        HRR(1:LB,49,1)
      !=|3,32)
      HRR(1:LB,3,32)=CDz*(CDz*HRR(1:LB,14,1)+  &
                        2.D0*HRR(1:LB,29,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,7,1)+  &
                        4.D0*HRR(1:LB,17,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,3,1)+  &
                        2.D0*HRR(1:LB,9,1))+  &
                        HRR(1:LB,19,1))+  &
                        2.D0*HRR(1:LB,32,1))+  &
                        HRR(1:LB,50,1)
      !=|4,32)
      HRR(1:LB,4,32)=CDz*(CDz*HRR(1:LB,17,1)+  &
                        2.D0*HRR(1:LB,32,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,9,1)+  &
                        4.D0*HRR(1:LB,19,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,4,1)+  &
                        2.D0*HRR(1:LB,10,1))+  &
                        HRR(1:LB,20,1))+  &
                        2.D0*HRR(1:LB,34,1))+  &
                        HRR(1:LB,53,1)
      !=|2,33)
      HRR(1:LB,2,33)=CDz*(CDz*(CDz*HRR(1:LB,5,1)+  &
                        3.D0*HRR(1:LB,15,1))+  &
                        3.D0*HRR(1:LB,30,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,2,1)+  &
                        3.D0*HRR(1:LB,8,1))+  &
                        3.D0*HRR(1:LB,18,1))+  &
                        HRR(1:LB,33,1))+  &
                        HRR(1:LB,51,1)
      !=|3,33)
      HRR(1:LB,3,33)=CDz*(CDz*(CDz*HRR(1:LB,6,1)+  &
                        3.D0*HRR(1:LB,16,1))+  &
                        3.D0*HRR(1:LB,31,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,3,1)+  &
                        3.D0*HRR(1:LB,9,1))+  &
                        3.D0*HRR(1:LB,19,1))+  &
                        HRR(1:LB,34,1))+  &
                        HRR(1:LB,52,1)
      !=|4,33)
      HRR(1:LB,4,33)=CDz*(CDz*(CDz*HRR(1:LB,8,1)+  &
                        3.D0*HRR(1:LB,18,1))+  &
                        3.D0*HRR(1:LB,33,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,4,1)+  &
                        3.D0*HRR(1:LB,10,1))+  &
                        3.D0*HRR(1:LB,20,1))+  &
                        HRR(1:LB,35,1))+  &
                        HRR(1:LB,54,1)
      !=|2,34)
      HRR(1:LB,2,34)=CDz*(CDz*(CDz*HRR(1:LB,6,1)+  &
                        3.D0*HRR(1:LB,16,1))+  &
                        3.D0*HRR(1:LB,31,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,2,1)+  &
                        3.D0*HRR(1:LB,8,1))+  &
                        3.D0*HRR(1:LB,18,1))+  &
                        HRR(1:LB,33,1))+  &
                        HRR(1:LB,52,1)
      !=|3,34)
      HRR(1:LB,3,34)=CDz*(CDz*(CDz*HRR(1:LB,7,1)+  &
                        3.D0*HRR(1:LB,17,1))+  &
                        3.D0*HRR(1:LB,32,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,3,1)+  &
                        3.D0*HRR(1:LB,9,1))+  &
                        3.D0*HRR(1:LB,19,1))+  &
                        HRR(1:LB,34,1))+  &
                        HRR(1:LB,53,1)
      !=|4,34)
      HRR(1:LB,4,34)=CDz*(CDz*(CDz*HRR(1:LB,9,1)+  &
                        3.D0*HRR(1:LB,19,1))+  &
                        3.D0*HRR(1:LB,34,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,4,1)+  &
                        3.D0*HRR(1:LB,10,1))+  &
                        3.D0*HRR(1:LB,20,1))+  &
                        HRR(1:LB,35,1))+  &
                        HRR(1:LB,55,1)
      !=|2,35)
      HRR(1:LB,2,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,2,1)+  &
                        4.D0*HRR(1:LB,8,1))+  &
                        6.D0*HRR(1:LB,18,1))+  &
                        4.D0*HRR(1:LB,33,1))+  &
                        HRR(1:LB,54,1)
      !=|3,35)
      HRR(1:LB,3,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,3,1)+  &
                        4.D0*HRR(1:LB,9,1))+  &
                        6.D0*HRR(1:LB,19,1))+  &
                        4.D0*HRR(1:LB,34,1))+  &
                        HRR(1:LB,55,1)
      !=|4,35)
      HRR(1:LB,4,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,4,1)+  &
                        4.D0*HRR(1:LB,10,1))+  &
                        6.D0*HRR(1:LB,20,1))+  &
                        4.D0*HRR(1:LB,35,1))+  &
                        HRR(1:LB,56,1)
      !=|1,21)
      HRR(1:LB,1,21)=CDx*(CDx*(CDx*(CDx*HRR(1:LB,1,1)+  &
                        4.D0*HRR(1:LB,2,1))+  &
                        6.D0*HRR(1:LB,5,1))+  &
                        4.D0*HRR(1:LB,11,1))+  &
                        HRR(1:LB,21,1)
      !=|1,22)
      HRR(1:LB,1,22)=CDy*HRR(1:LB,11,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,5,1)+  &
                        CDx*(3.D0*CDy*HRR(1:LB,2,1)+  &
                        CDx*(CDy*HRR(1:LB,1,1)+  &
                        HRR(1:LB,3,1))+  &
                        3.D0*HRR(1:LB,6,1))+  &
                        3.D0*HRR(1:LB,12,1))+  &
                        HRR(1:LB,22,1)
      !=|1,23)
      HRR(1:LB,1,23)=CDy*(CDy*HRR(1:LB,5,1)+  &
                        2.D0*HRR(1:LB,12,1))+  &
                        CDx*(CDy*(2.D0*CDy*HRR(1:LB,2,1)+  &
                        4.D0*HRR(1:LB,6,1))+  &
                        CDx*(CDy*(CDy*HRR(1:LB,1,1)+  &
                        2.D0*HRR(1:LB,3,1))+  &
                        HRR(1:LB,7,1))+  &
                        2.D0*HRR(1:LB,13,1))+  &
                        HRR(1:LB,23,1)
      !=|1,24)
      HRR(1:LB,1,24)=CDy*(CDy*(CDy*HRR(1:LB,2,1)+  &
                        3.D0*HRR(1:LB,6,1))+  &
                        3.D0*HRR(1:LB,13,1))+  &
                        CDx*(CDy*(CDy*(CDy*HRR(1:LB,1,1)+  &
                        3.D0*HRR(1:LB,3,1))+  &
                        3.D0*HRR(1:LB,7,1))+  &
                        HRR(1:LB,14,1))+  &
                        HRR(1:LB,24,1)
      !=|1,25)
      HRR(1:LB,1,25)=CDy*(CDy*(CDy*(CDy*HRR(1:LB,1,1)+  &
                        4.D0*HRR(1:LB,3,1))+  &
                        6.D0*HRR(1:LB,7,1))+  &
                        4.D0*HRR(1:LB,14,1))+  &
                        HRR(1:LB,25,1)
      !=|1,26)
      HRR(1:LB,1,26)=CDz*HRR(1:LB,11,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,5,1)+  &
                        CDx*(3.D0*CDz*HRR(1:LB,2,1)+  &
                        CDx*(CDz*HRR(1:LB,1,1)+  &
                        HRR(1:LB,4,1))+  &
                        3.D0*HRR(1:LB,8,1))+  &
                        3.D0*HRR(1:LB,15,1))+  &
                        HRR(1:LB,26,1)
      !=|1,27)
      HRR(1:LB,1,27)=CDz*HRR(1:LB,12,1)+  &
                        CDy*(CDz*HRR(1:LB,5,1)+  &
                        HRR(1:LB,15,1))+  &
                        CDx*(2.D0*CDz*HRR(1:LB,6,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,2,1)+  &
                        2.D0*HRR(1:LB,8,1))+  &
                        CDx*(CDz*HRR(1:LB,3,1)+  &
                        CDy*(CDz*HRR(1:LB,1,1)+  &
                        HRR(1:LB,4,1))+  &
                        HRR(1:LB,9,1))+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        HRR(1:LB,27,1)
      !=|1,28)
      HRR(1:LB,1,28)=CDz*HRR(1:LB,13,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,6,1)+  &
                        CDy*(CDz*HRR(1:LB,2,1)+  &
                        HRR(1:LB,8,1))+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        CDx*(CDz*HRR(1:LB,7,1)+  &
                        CDy*(2.D0*CDz*HRR(1:LB,3,1)+  &
                        CDy*(CDz*HRR(1:LB,1,1)+  &
                        HRR(1:LB,4,1))+  &
                        2.D0*HRR(1:LB,9,1))+  &
                        HRR(1:LB,17,1))+  &
                        HRR(1:LB,28,1)
      !=|1,29)
      HRR(1:LB,1,29)=CDz*HRR(1:LB,14,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,7,1)+  &
                        CDy*(3.D0*CDz*HRR(1:LB,3,1)+  &
                        CDy*(CDz*HRR(1:LB,1,1)+  &
                        HRR(1:LB,4,1))+  &
                        3.D0*HRR(1:LB,9,1))+  &
                        3.D0*HRR(1:LB,17,1))+  &
                        HRR(1:LB,29,1)
      !=|1,30)
      HRR(1:LB,1,30)=CDz*(CDz*HRR(1:LB,5,1)+  &
                        2.D0*HRR(1:LB,15,1))+  &
                        CDx*(CDz*(2.D0*CDz*HRR(1:LB,2,1)+  &
                        4.D0*HRR(1:LB,8,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,1,1)+  &
                        2.D0*HRR(1:LB,4,1))+  &
                        HRR(1:LB,10,1))+  &
                        2.D0*HRR(1:LB,18,1))+  &
                        HRR(1:LB,30,1)
      !=|1,31)
      HRR(1:LB,1,31)=CDz*(CDz*HRR(1:LB,6,1)+  &
                        2.D0*HRR(1:LB,16,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,2,1)+  &
                        2.D0*HRR(1:LB,8,1))+  &
                        HRR(1:LB,18,1))+  &
                        CDx*(CDz*(CDz*HRR(1:LB,3,1)+  &
                        2.D0*HRR(1:LB,9,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,1,1)+  &
                        2.D0*HRR(1:LB,4,1))+  &
                        HRR(1:LB,10,1))+  &
                        HRR(1:LB,19,1))+  &
                        HRR(1:LB,31,1)
      !=|1,32)
      HRR(1:LB,1,32)=CDz*(CDz*HRR(1:LB,7,1)+  &
                        2.D0*HRR(1:LB,17,1))+  &
                        CDy*(CDz*(2.D0*CDz*HRR(1:LB,3,1)+  &
                        4.D0*HRR(1:LB,9,1))+  &
                        CDy*(CDz*(CDz*HRR(1:LB,1,1)+  &
                        2.D0*HRR(1:LB,4,1))+  &
                        HRR(1:LB,10,1))+  &
                        2.D0*HRR(1:LB,19,1))+  &
                        HRR(1:LB,32,1)
      !=|1,33)
      HRR(1:LB,1,33)=CDz*(CDz*(CDz*HRR(1:LB,2,1)+  &
                        3.D0*HRR(1:LB,8,1))+  &
                        3.D0*HRR(1:LB,18,1))+  &
                        CDx*(CDz*(CDz*(CDz*HRR(1:LB,1,1)+  &
                        3.D0*HRR(1:LB,4,1))+  &
                        3.D0*HRR(1:LB,10,1))+  &
                        HRR(1:LB,20,1))+  &
                        HRR(1:LB,33,1)
      !=|1,34)
      HRR(1:LB,1,34)=CDz*(CDz*(CDz*HRR(1:LB,3,1)+  &
                        3.D0*HRR(1:LB,9,1))+  &
                        3.D0*HRR(1:LB,19,1))+  &
                        CDy*(CDz*(CDz*(CDz*HRR(1:LB,1,1)+  &
                        3.D0*HRR(1:LB,4,1))+  &
                        3.D0*HRR(1:LB,10,1))+  &
                        HRR(1:LB,20,1))+  &
                        HRR(1:LB,34,1)
      !=|1,35)
      HRR(1:LB,1,35)=CDz*(CDz*(CDz*(CDz*HRR(1:LB,1,1)+  &
                        4.D0*HRR(1:LB,4,1))+  &
                        6.D0*HRR(1:LB,10,1))+  &
                        4.D0*HRR(1:LB,20,1))+  &
                        HRR(1:LB,35,1)
END SUBROUTINE KetHRR1515
