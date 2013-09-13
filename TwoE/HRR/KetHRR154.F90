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
   SUBROUTINE KetHRR154(LB,HRR)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB
      REAL(DOUBLE) :: HRR(1:LB,84,10)
      !=|21,5)
      HRR(1:LB,21,5)=CDx*(CDx*HRR(1:LB,21,1)+  &
                        2.D0*HRR(1:LB,36,1))+  &
                        HRR(1:LB,57,1)
      !=|22,5)
      HRR(1:LB,22,5)=CDx*(CDx*HRR(1:LB,22,1)+  &
                        2.D0*HRR(1:LB,37,1))+  &
                        HRR(1:LB,58,1)
      !=|23,5)
      HRR(1:LB,23,5)=CDx*(CDx*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,38,1))+  &
                        HRR(1:LB,59,1)
      !=|24,5)
      HRR(1:LB,24,5)=CDx*(CDx*HRR(1:LB,24,1)+  &
                        2.D0*HRR(1:LB,39,1))+  &
                        HRR(1:LB,60,1)
      !=|25,5)
      HRR(1:LB,25,5)=CDx*(CDx*HRR(1:LB,25,1)+  &
                        2.D0*HRR(1:LB,40,1))+  &
                        HRR(1:LB,61,1)
      !=|26,5)
      HRR(1:LB,26,5)=CDx*(CDx*HRR(1:LB,26,1)+  &
                        2.D0*HRR(1:LB,42,1))+  &
                        HRR(1:LB,64,1)
      !=|27,5)
      HRR(1:LB,27,5)=CDx*(CDx*HRR(1:LB,27,1)+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,65,1)
      !=|28,5)
      HRR(1:LB,28,5)=CDx*(CDx*HRR(1:LB,28,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,66,1)
      !=|29,5)
      HRR(1:LB,29,5)=CDx*(CDx*HRR(1:LB,29,1)+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,67,1)
      !=|30,5)
      HRR(1:LB,30,5)=CDx*(CDx*HRR(1:LB,30,1)+  &
                        2.D0*HRR(1:LB,47,1))+  &
                        HRR(1:LB,70,1)
      !=|31,5)
      HRR(1:LB,31,5)=CDx*(CDx*HRR(1:LB,31,1)+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,71,1)
      !=|32,5)
      HRR(1:LB,32,5)=CDx*(CDx*HRR(1:LB,32,1)+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,72,1)
      !=|33,5)
      HRR(1:LB,33,5)=CDx*(CDx*HRR(1:LB,33,1)+  &
                        2.D0*HRR(1:LB,51,1))+  &
                        HRR(1:LB,75,1)
      !=|34,5)
      HRR(1:LB,34,5)=CDx*(CDx*HRR(1:LB,34,1)+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,76,1)
      !=|35,5)
      HRR(1:LB,35,5)=CDx*(CDx*HRR(1:LB,35,1)+  &
                        2.D0*HRR(1:LB,54,1))+  &
                        HRR(1:LB,79,1)
      !=|21,6)
      HRR(1:LB,21,6)=CDy*HRR(1:LB,36,1)+  &
                        CDx*(CDy*HRR(1:LB,21,1)+  &
                        HRR(1:LB,37,1))+  &
                        HRR(1:LB,58,1)
      !=|22,6)
      HRR(1:LB,22,6)=CDy*HRR(1:LB,37,1)+  &
                        CDx*(CDy*HRR(1:LB,22,1)+  &
                        HRR(1:LB,38,1))+  &
                        HRR(1:LB,59,1)
      !=|23,6)
      HRR(1:LB,23,6)=CDy*HRR(1:LB,38,1)+  &
                        CDx*(CDy*HRR(1:LB,23,1)+  &
                        HRR(1:LB,39,1))+  &
                        HRR(1:LB,60,1)
      !=|24,6)
      HRR(1:LB,24,6)=CDy*HRR(1:LB,39,1)+  &
                        CDx*(CDy*HRR(1:LB,24,1)+  &
                        HRR(1:LB,40,1))+  &
                        HRR(1:LB,61,1)
      !=|25,6)
      HRR(1:LB,25,6)=CDy*HRR(1:LB,40,1)+  &
                        CDx*(CDy*HRR(1:LB,25,1)+  &
                        HRR(1:LB,41,1))+  &
                        HRR(1:LB,62,1)
      !=|26,6)
      HRR(1:LB,26,6)=CDy*HRR(1:LB,42,1)+  &
                        CDx*(CDy*HRR(1:LB,26,1)+  &
                        HRR(1:LB,43,1))+  &
                        HRR(1:LB,65,1)
      !=|27,6)
      HRR(1:LB,27,6)=CDy*HRR(1:LB,43,1)+  &
                        CDx*(CDy*HRR(1:LB,27,1)+  &
                        HRR(1:LB,44,1))+  &
                        HRR(1:LB,66,1)
      !=|28,6)
      HRR(1:LB,28,6)=CDy*HRR(1:LB,44,1)+  &
                        CDx*(CDy*HRR(1:LB,28,1)+  &
                        HRR(1:LB,45,1))+  &
                        HRR(1:LB,67,1)
      !=|29,6)
      HRR(1:LB,29,6)=CDy*HRR(1:LB,45,1)+  &
                        CDx*(CDy*HRR(1:LB,29,1)+  &
                        HRR(1:LB,46,1))+  &
                        HRR(1:LB,68,1)
      !=|30,6)
      HRR(1:LB,30,6)=CDy*HRR(1:LB,47,1)+  &
                        CDx*(CDy*HRR(1:LB,30,1)+  &
                        HRR(1:LB,48,1))+  &
                        HRR(1:LB,71,1)
      !=|31,6)
      HRR(1:LB,31,6)=CDy*HRR(1:LB,48,1)+  &
                        CDx*(CDy*HRR(1:LB,31,1)+  &
                        HRR(1:LB,49,1))+  &
                        HRR(1:LB,72,1)
      !=|32,6)
      HRR(1:LB,32,6)=CDy*HRR(1:LB,49,1)+  &
                        CDx*(CDy*HRR(1:LB,32,1)+  &
                        HRR(1:LB,50,1))+  &
                        HRR(1:LB,73,1)
      !=|33,6)
      HRR(1:LB,33,6)=CDy*HRR(1:LB,51,1)+  &
                        CDx*(CDy*HRR(1:LB,33,1)+  &
                        HRR(1:LB,52,1))+  &
                        HRR(1:LB,76,1)
      !=|34,6)
      HRR(1:LB,34,6)=CDy*HRR(1:LB,52,1)+  &
                        CDx*(CDy*HRR(1:LB,34,1)+  &
                        HRR(1:LB,53,1))+  &
                        HRR(1:LB,77,1)
      !=|35,6)
      HRR(1:LB,35,6)=CDy*HRR(1:LB,54,1)+  &
                        CDx*(CDy*HRR(1:LB,35,1)+  &
                        HRR(1:LB,55,1))+  &
                        HRR(1:LB,80,1)
      !=|21,7)
      HRR(1:LB,21,7)=CDy*(CDy*HRR(1:LB,21,1)+  &
                        2.D0*HRR(1:LB,37,1))+  &
                        HRR(1:LB,59,1)
      !=|22,7)
      HRR(1:LB,22,7)=CDy*(CDy*HRR(1:LB,22,1)+  &
                        2.D0*HRR(1:LB,38,1))+  &
                        HRR(1:LB,60,1)
      !=|23,7)
      HRR(1:LB,23,7)=CDy*(CDy*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,39,1))+  &
                        HRR(1:LB,61,1)
      !=|24,7)
      HRR(1:LB,24,7)=CDy*(CDy*HRR(1:LB,24,1)+  &
                        2.D0*HRR(1:LB,40,1))+  &
                        HRR(1:LB,62,1)
      !=|25,7)
      HRR(1:LB,25,7)=CDy*(CDy*HRR(1:LB,25,1)+  &
                        2.D0*HRR(1:LB,41,1))+  &
                        HRR(1:LB,63,1)
      !=|26,7)
      HRR(1:LB,26,7)=CDy*(CDy*HRR(1:LB,26,1)+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,66,1)
      !=|27,7)
      HRR(1:LB,27,7)=CDy*(CDy*HRR(1:LB,27,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,67,1)
      !=|28,7)
      HRR(1:LB,28,7)=CDy*(CDy*HRR(1:LB,28,1)+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,68,1)
      !=|29,7)
      HRR(1:LB,29,7)=CDy*(CDy*HRR(1:LB,29,1)+  &
                        2.D0*HRR(1:LB,46,1))+  &
                        HRR(1:LB,69,1)
      !=|30,7)
      HRR(1:LB,30,7)=CDy*(CDy*HRR(1:LB,30,1)+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,72,1)
      !=|31,7)
      HRR(1:LB,31,7)=CDy*(CDy*HRR(1:LB,31,1)+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,73,1)
      !=|32,7)
      HRR(1:LB,32,7)=CDy*(CDy*HRR(1:LB,32,1)+  &
                        2.D0*HRR(1:LB,50,1))+  &
                        HRR(1:LB,74,1)
      !=|33,7)
      HRR(1:LB,33,7)=CDy*(CDy*HRR(1:LB,33,1)+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,77,1)
      !=|34,7)
      HRR(1:LB,34,7)=CDy*(CDy*HRR(1:LB,34,1)+  &
                        2.D0*HRR(1:LB,53,1))+  &
                        HRR(1:LB,78,1)
      !=|35,7)
      HRR(1:LB,35,7)=CDy*(CDy*HRR(1:LB,35,1)+  &
                        2.D0*HRR(1:LB,55,1))+  &
                        HRR(1:LB,81,1)
      !=|21,8)
      HRR(1:LB,21,8)=CDz*HRR(1:LB,36,1)+  &
                        CDx*(CDz*HRR(1:LB,21,1)+  &
                        HRR(1:LB,42,1))+  &
                        HRR(1:LB,64,1)
      !=|22,8)
      HRR(1:LB,22,8)=CDz*HRR(1:LB,37,1)+  &
                        CDx*(CDz*HRR(1:LB,22,1)+  &
                        HRR(1:LB,43,1))+  &
                        HRR(1:LB,65,1)
      !=|23,8)
      HRR(1:LB,23,8)=CDz*HRR(1:LB,38,1)+  &
                        CDx*(CDz*HRR(1:LB,23,1)+  &
                        HRR(1:LB,44,1))+  &
                        HRR(1:LB,66,1)
      !=|24,8)
      HRR(1:LB,24,8)=CDz*HRR(1:LB,39,1)+  &
                        CDx*(CDz*HRR(1:LB,24,1)+  &
                        HRR(1:LB,45,1))+  &
                        HRR(1:LB,67,1)
      !=|25,8)
      HRR(1:LB,25,8)=CDz*HRR(1:LB,40,1)+  &
                        CDx*(CDz*HRR(1:LB,25,1)+  &
                        HRR(1:LB,46,1))+  &
                        HRR(1:LB,68,1)
      !=|26,8)
      HRR(1:LB,26,8)=CDz*HRR(1:LB,42,1)+  &
                        CDx*(CDz*HRR(1:LB,26,1)+  &
                        HRR(1:LB,47,1))+  &
                        HRR(1:LB,70,1)
      !=|27,8)
      HRR(1:LB,27,8)=CDz*HRR(1:LB,43,1)+  &
                        CDx*(CDz*HRR(1:LB,27,1)+  &
                        HRR(1:LB,48,1))+  &
                        HRR(1:LB,71,1)
      !=|28,8)
      HRR(1:LB,28,8)=CDz*HRR(1:LB,44,1)+  &
                        CDx*(CDz*HRR(1:LB,28,1)+  &
                        HRR(1:LB,49,1))+  &
                        HRR(1:LB,72,1)
      !=|29,8)
      HRR(1:LB,29,8)=CDz*HRR(1:LB,45,1)+  &
                        CDx*(CDz*HRR(1:LB,29,1)+  &
                        HRR(1:LB,50,1))+  &
                        HRR(1:LB,73,1)
      !=|30,8)
      HRR(1:LB,30,8)=CDz*HRR(1:LB,47,1)+  &
                        CDx*(CDz*HRR(1:LB,30,1)+  &
                        HRR(1:LB,51,1))+  &
                        HRR(1:LB,75,1)
      !=|31,8)
      HRR(1:LB,31,8)=CDz*HRR(1:LB,48,1)+  &
                        CDx*(CDz*HRR(1:LB,31,1)+  &
                        HRR(1:LB,52,1))+  &
                        HRR(1:LB,76,1)
      !=|32,8)
      HRR(1:LB,32,8)=CDz*HRR(1:LB,49,1)+  &
                        CDx*(CDz*HRR(1:LB,32,1)+  &
                        HRR(1:LB,53,1))+  &
                        HRR(1:LB,77,1)
      !=|33,8)
      HRR(1:LB,33,8)=CDz*HRR(1:LB,51,1)+  &
                        CDx*(CDz*HRR(1:LB,33,1)+  &
                        HRR(1:LB,54,1))+  &
                        HRR(1:LB,79,1)
      !=|34,8)
      HRR(1:LB,34,8)=CDz*HRR(1:LB,52,1)+  &
                        CDx*(CDz*HRR(1:LB,34,1)+  &
                        HRR(1:LB,55,1))+  &
                        HRR(1:LB,80,1)
      !=|35,8)
      HRR(1:LB,35,8)=CDz*HRR(1:LB,54,1)+  &
                        CDx*(CDz*HRR(1:LB,35,1)+  &
                        HRR(1:LB,56,1))+  &
                        HRR(1:LB,82,1)
      !=|21,9)
      HRR(1:LB,21,9)=CDz*HRR(1:LB,37,1)+  &
                        CDy*(CDz*HRR(1:LB,21,1)+  &
                        HRR(1:LB,42,1))+  &
                        HRR(1:LB,65,1)
      !=|22,9)
      HRR(1:LB,22,9)=CDz*HRR(1:LB,38,1)+  &
                        CDy*(CDz*HRR(1:LB,22,1)+  &
                        HRR(1:LB,43,1))+  &
                        HRR(1:LB,66,1)
      !=|23,9)
      HRR(1:LB,23,9)=CDz*HRR(1:LB,39,1)+  &
                        CDy*(CDz*HRR(1:LB,23,1)+  &
                        HRR(1:LB,44,1))+  &
                        HRR(1:LB,67,1)
      !=|24,9)
      HRR(1:LB,24,9)=CDz*HRR(1:LB,40,1)+  &
                        CDy*(CDz*HRR(1:LB,24,1)+  &
                        HRR(1:LB,45,1))+  &
                        HRR(1:LB,68,1)
      !=|25,9)
      HRR(1:LB,25,9)=CDz*HRR(1:LB,41,1)+  &
                        CDy*(CDz*HRR(1:LB,25,1)+  &
                        HRR(1:LB,46,1))+  &
                        HRR(1:LB,69,1)
      !=|26,9)
      HRR(1:LB,26,9)=CDz*HRR(1:LB,43,1)+  &
                        CDy*(CDz*HRR(1:LB,26,1)+  &
                        HRR(1:LB,47,1))+  &
                        HRR(1:LB,71,1)
      !=|27,9)
      HRR(1:LB,27,9)=CDz*HRR(1:LB,44,1)+  &
                        CDy*(CDz*HRR(1:LB,27,1)+  &
                        HRR(1:LB,48,1))+  &
                        HRR(1:LB,72,1)
      !=|28,9)
      HRR(1:LB,28,9)=CDz*HRR(1:LB,45,1)+  &
                        CDy*(CDz*HRR(1:LB,28,1)+  &
                        HRR(1:LB,49,1))+  &
                        HRR(1:LB,73,1)
      !=|29,9)
      HRR(1:LB,29,9)=CDz*HRR(1:LB,46,1)+  &
                        CDy*(CDz*HRR(1:LB,29,1)+  &
                        HRR(1:LB,50,1))+  &
                        HRR(1:LB,74,1)
      !=|30,9)
      HRR(1:LB,30,9)=CDz*HRR(1:LB,48,1)+  &
                        CDy*(CDz*HRR(1:LB,30,1)+  &
                        HRR(1:LB,51,1))+  &
                        HRR(1:LB,76,1)
      !=|31,9)
      HRR(1:LB,31,9)=CDz*HRR(1:LB,49,1)+  &
                        CDy*(CDz*HRR(1:LB,31,1)+  &
                        HRR(1:LB,52,1))+  &
                        HRR(1:LB,77,1)
      !=|32,9)
      HRR(1:LB,32,9)=CDz*HRR(1:LB,50,1)+  &
                        CDy*(CDz*HRR(1:LB,32,1)+  &
                        HRR(1:LB,53,1))+  &
                        HRR(1:LB,78,1)
      !=|33,9)
      HRR(1:LB,33,9)=CDz*HRR(1:LB,52,1)+  &
                        CDy*(CDz*HRR(1:LB,33,1)+  &
                        HRR(1:LB,54,1))+  &
                        HRR(1:LB,80,1)
      !=|34,9)
      HRR(1:LB,34,9)=CDz*HRR(1:LB,53,1)+  &
                        CDy*(CDz*HRR(1:LB,34,1)+  &
                        HRR(1:LB,55,1))+  &
                        HRR(1:LB,81,1)
      !=|35,9)
      HRR(1:LB,35,9)=CDz*HRR(1:LB,55,1)+  &
                        CDy*(CDz*HRR(1:LB,35,1)+  &
                        HRR(1:LB,56,1))+  &
                        HRR(1:LB,83,1)
      !=|21,10)
      HRR(1:LB,21,10)=CDz*(CDz*HRR(1:LB,21,1)+  &
                        2.D0*HRR(1:LB,42,1))+  &
                        HRR(1:LB,70,1)
      !=|22,10)
      HRR(1:LB,22,10)=CDz*(CDz*HRR(1:LB,22,1)+  &
                        2.D0*HRR(1:LB,43,1))+  &
                        HRR(1:LB,71,1)
      !=|23,10)
      HRR(1:LB,23,10)=CDz*(CDz*HRR(1:LB,23,1)+  &
                        2.D0*HRR(1:LB,44,1))+  &
                        HRR(1:LB,72,1)
      !=|24,10)
      HRR(1:LB,24,10)=CDz*(CDz*HRR(1:LB,24,1)+  &
                        2.D0*HRR(1:LB,45,1))+  &
                        HRR(1:LB,73,1)
      !=|25,10)
      HRR(1:LB,25,10)=CDz*(CDz*HRR(1:LB,25,1)+  &
                        2.D0*HRR(1:LB,46,1))+  &
                        HRR(1:LB,74,1)
      !=|26,10)
      HRR(1:LB,26,10)=CDz*(CDz*HRR(1:LB,26,1)+  &
                        2.D0*HRR(1:LB,47,1))+  &
                        HRR(1:LB,75,1)
      !=|27,10)
      HRR(1:LB,27,10)=CDz*(CDz*HRR(1:LB,27,1)+  &
                        2.D0*HRR(1:LB,48,1))+  &
                        HRR(1:LB,76,1)
      !=|28,10)
      HRR(1:LB,28,10)=CDz*(CDz*HRR(1:LB,28,1)+  &
                        2.D0*HRR(1:LB,49,1))+  &
                        HRR(1:LB,77,1)
      !=|29,10)
      HRR(1:LB,29,10)=CDz*(CDz*HRR(1:LB,29,1)+  &
                        2.D0*HRR(1:LB,50,1))+  &
                        HRR(1:LB,78,1)
      !=|30,10)
      HRR(1:LB,30,10)=CDz*(CDz*HRR(1:LB,30,1)+  &
                        2.D0*HRR(1:LB,51,1))+  &
                        HRR(1:LB,79,1)
      !=|31,10)
      HRR(1:LB,31,10)=CDz*(CDz*HRR(1:LB,31,1)+  &
                        2.D0*HRR(1:LB,52,1))+  &
                        HRR(1:LB,80,1)
      !=|32,10)
      HRR(1:LB,32,10)=CDz*(CDz*HRR(1:LB,32,1)+  &
                        2.D0*HRR(1:LB,53,1))+  &
                        HRR(1:LB,81,1)
      !=|33,10)
      HRR(1:LB,33,10)=CDz*(CDz*HRR(1:LB,33,1)+  &
                        2.D0*HRR(1:LB,54,1))+  &
                        HRR(1:LB,82,1)
      !=|34,10)
      HRR(1:LB,34,10)=CDz*(CDz*HRR(1:LB,34,1)+  &
                        2.D0*HRR(1:LB,55,1))+  &
                        HRR(1:LB,83,1)
      !=|35,10)
      HRR(1:LB,35,10)=CDz*(CDz*HRR(1:LB,35,1)+  &
                        2.D0*HRR(1:LB,56,1))+  &
                        HRR(1:LB,84,1)
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
      !=|21,2)
      HRR(1:LB,21,2)=CDx*HRR(1:LB,21,1)+  &
                        HRR(1:LB,36,1)
      !=|22,2)
      HRR(1:LB,22,2)=CDx*HRR(1:LB,22,1)+  &
                        HRR(1:LB,37,1)
      !=|23,2)
      HRR(1:LB,23,2)=CDx*HRR(1:LB,23,1)+  &
                        HRR(1:LB,38,1)
      !=|24,2)
      HRR(1:LB,24,2)=CDx*HRR(1:LB,24,1)+  &
                        HRR(1:LB,39,1)
      !=|25,2)
      HRR(1:LB,25,2)=CDx*HRR(1:LB,25,1)+  &
                        HRR(1:LB,40,1)
      !=|26,2)
      HRR(1:LB,26,2)=CDx*HRR(1:LB,26,1)+  &
                        HRR(1:LB,42,1)
      !=|27,2)
      HRR(1:LB,27,2)=CDx*HRR(1:LB,27,1)+  &
                        HRR(1:LB,43,1)
      !=|28,2)
      HRR(1:LB,28,2)=CDx*HRR(1:LB,28,1)+  &
                        HRR(1:LB,44,1)
      !=|29,2)
      HRR(1:LB,29,2)=CDx*HRR(1:LB,29,1)+  &
                        HRR(1:LB,45,1)
      !=|30,2)
      HRR(1:LB,30,2)=CDx*HRR(1:LB,30,1)+  &
                        HRR(1:LB,47,1)
      !=|31,2)
      HRR(1:LB,31,2)=CDx*HRR(1:LB,31,1)+  &
                        HRR(1:LB,48,1)
      !=|32,2)
      HRR(1:LB,32,2)=CDx*HRR(1:LB,32,1)+  &
                        HRR(1:LB,49,1)
      !=|33,2)
      HRR(1:LB,33,2)=CDx*HRR(1:LB,33,1)+  &
                        HRR(1:LB,51,1)
      !=|34,2)
      HRR(1:LB,34,2)=CDx*HRR(1:LB,34,1)+  &
                        HRR(1:LB,52,1)
      !=|35,2)
      HRR(1:LB,35,2)=CDx*HRR(1:LB,35,1)+  &
                        HRR(1:LB,54,1)
      !=|21,3)
      HRR(1:LB,21,3)=CDy*HRR(1:LB,21,1)+  &
                        HRR(1:LB,37,1)
      !=|22,3)
      HRR(1:LB,22,3)=CDy*HRR(1:LB,22,1)+  &
                        HRR(1:LB,38,1)
      !=|23,3)
      HRR(1:LB,23,3)=CDy*HRR(1:LB,23,1)+  &
                        HRR(1:LB,39,1)
      !=|24,3)
      HRR(1:LB,24,3)=CDy*HRR(1:LB,24,1)+  &
                        HRR(1:LB,40,1)
      !=|25,3)
      HRR(1:LB,25,3)=CDy*HRR(1:LB,25,1)+  &
                        HRR(1:LB,41,1)
      !=|26,3)
      HRR(1:LB,26,3)=CDy*HRR(1:LB,26,1)+  &
                        HRR(1:LB,43,1)
      !=|27,3)
      HRR(1:LB,27,3)=CDy*HRR(1:LB,27,1)+  &
                        HRR(1:LB,44,1)
      !=|28,3)
      HRR(1:LB,28,3)=CDy*HRR(1:LB,28,1)+  &
                        HRR(1:LB,45,1)
      !=|29,3)
      HRR(1:LB,29,3)=CDy*HRR(1:LB,29,1)+  &
                        HRR(1:LB,46,1)
      !=|30,3)
      HRR(1:LB,30,3)=CDy*HRR(1:LB,30,1)+  &
                        HRR(1:LB,48,1)
      !=|31,3)
      HRR(1:LB,31,3)=CDy*HRR(1:LB,31,1)+  &
                        HRR(1:LB,49,1)
      !=|32,3)
      HRR(1:LB,32,3)=CDy*HRR(1:LB,32,1)+  &
                        HRR(1:LB,50,1)
      !=|33,3)
      HRR(1:LB,33,3)=CDy*HRR(1:LB,33,1)+  &
                        HRR(1:LB,52,1)
      !=|34,3)
      HRR(1:LB,34,3)=CDy*HRR(1:LB,34,1)+  &
                        HRR(1:LB,53,1)
      !=|35,3)
      HRR(1:LB,35,3)=CDy*HRR(1:LB,35,1)+  &
                        HRR(1:LB,55,1)
      !=|21,4)
      HRR(1:LB,21,4)=CDz*HRR(1:LB,21,1)+  &
                        HRR(1:LB,42,1)
      !=|22,4)
      HRR(1:LB,22,4)=CDz*HRR(1:LB,22,1)+  &
                        HRR(1:LB,43,1)
      !=|23,4)
      HRR(1:LB,23,4)=CDz*HRR(1:LB,23,1)+  &
                        HRR(1:LB,44,1)
      !=|24,4)
      HRR(1:LB,24,4)=CDz*HRR(1:LB,24,1)+  &
                        HRR(1:LB,45,1)
      !=|25,4)
      HRR(1:LB,25,4)=CDz*HRR(1:LB,25,1)+  &
                        HRR(1:LB,46,1)
      !=|26,4)
      HRR(1:LB,26,4)=CDz*HRR(1:LB,26,1)+  &
                        HRR(1:LB,47,1)
      !=|27,4)
      HRR(1:LB,27,4)=CDz*HRR(1:LB,27,1)+  &
                        HRR(1:LB,48,1)
      !=|28,4)
      HRR(1:LB,28,4)=CDz*HRR(1:LB,28,1)+  &
                        HRR(1:LB,49,1)
      !=|29,4)
      HRR(1:LB,29,4)=CDz*HRR(1:LB,29,1)+  &
                        HRR(1:LB,50,1)
      !=|30,4)
      HRR(1:LB,30,4)=CDz*HRR(1:LB,30,1)+  &
                        HRR(1:LB,51,1)
      !=|31,4)
      HRR(1:LB,31,4)=CDz*HRR(1:LB,31,1)+  &
                        HRR(1:LB,52,1)
      !=|32,4)
      HRR(1:LB,32,4)=CDz*HRR(1:LB,32,1)+  &
                        HRR(1:LB,53,1)
      !=|33,4)
      HRR(1:LB,33,4)=CDz*HRR(1:LB,33,1)+  &
                        HRR(1:LB,54,1)
      !=|34,4)
      HRR(1:LB,34,4)=CDz*HRR(1:LB,34,1)+  &
                        HRR(1:LB,55,1)
      !=|35,4)
      HRR(1:LB,35,4)=CDz*HRR(1:LB,35,1)+  &
                        HRR(1:LB,56,1)
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
END SUBROUTINE KetHRR154
