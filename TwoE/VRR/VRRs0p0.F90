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
   SUBROUTINE VRRs0p0(LB,LK,VRR0,VRR1)
      USE DerivedTypes
      USE VScratchB
      USE GlobalScalars
      IMPLICIT REAL(DOUBLE) (W)
      INTEGER :: LB,LK
      REAL(DOUBLE), DIMENSION(1:LB,1:LK) :: VRR0,VRR1
      VRR0(1,2)=QCx*VRR0(1,1)+WQx*VRR1(1,1)
      VRR0(1,3)=QCy*VRR0(1,1)+WQy*VRR1(1,1)
      VRR0(1,4)=QCz*VRR0(1,1)+WQz*VRR1(1,1)
END SUBROUTINE VRRs0p0
SUBROUTINE MVRRs0p0(IXYZ,LBS,LKS,VS0,VS1,LBR,LKR,VR1)
USE DerivedTypes
USE VScratchB
USE GlobalScalars
IMPLICIT NONE
INTEGER IXYZ,LBS,LKS,LBR,LKR
REAL(DOUBLE) VS0(LBS,LKS),VS1(LBS,LKS),VR1(LBR,LKR)
SELECT CASE(IXYZ)
CASE(1)
VS0(1,2)=QCx*VS0(1,1)+WQx*VS1(1,1)-r1x2E*VR1(1,1)
VS0(1,3)=QCy*VS0(1,1)+WQy*VS1(1,1)
VS0(1,4)=QCz*VS0(1,1)+WQz*VS1(1,1)
CASE(2)
VS0(1,2)=QCx*VS0(1,1)+WQx*VS1(1,1)
VS0(1,3)=QCy*VS0(1,1)+WQy*VS1(1,1)-r1x2E*VR1(1,1)
VS0(1,4)=QCz*VS0(1,1)+WQz*VS1(1,1)
CASE(3)
VS0(1,2)=QCx*VS0(1,1)+WQx*VS1(1,1)
VS0(1,3)=QCy*VS0(1,1)+WQy*VS1(1,1)
VS0(1,4)=QCz*VS0(1,1)+WQz*VS1(1,1)-r1x2E*VR1(1,1)
CASE DEFAULT
WRITE(*,*) 'STOP IN MVRRs0p0'
STOP
END SELECT
END SUBROUTINE MVRRs0p0
