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

#include "MondoConfig.h"

MODULE Numerics
  USE DerivedTypes
  IMPLICIT NONE
  !
  ! Thresholds for linear scaling routines:
  !
  ! Sparse blocked matrix thresholds
  REAL(DOUBLE),DIMENSION(4) :: TrixNeglect=(/1.D-4, 1.D-5, 1.D-6,  1.D-7 /)
  ! HiCu threshold
#ifdef PARALLEL
  REAL(DOUBLE),DIMENSION(4) :: CubeNeglect=(/1.D-4, 1.D-6, 1.D-8,  1.D-10 /)
#else
  REAL(DOUBLE),DIMENSION(4) :: CubeNeglect=(/1.D-3, 1.D-5, 1.D-7,  1.D-9 /)
#endif
  ! QCTC and ONX threshold
  REAL(DOUBLE),DIMENSION(4) :: TwoENeglect=(/1.D-6, 1.D-8, 1.D-10, 1.D-12/)
  ! Distribution threshold
  REAL(DOUBLE),DIMENSION(4) :: DistNeglect=(/1.D-8, 1.D-10,1.D-12, 1.D-14/)
  !
  ! Convergence criteria for SCF and force routines:
  !
  ! Max error in the relative total energy
  REAL(DOUBLE),DIMENSION(4) :: ETol       =(/ 1.D-5, 1.D-7, 1.D-9,  1.D-11 /)
  ! Max element of the density matrix
  REAL(DOUBLE),DIMENSION(4) :: DTol       =(/ 1.D-2, 1.D-3, 1.D-4,  1.D-5  /)
  ! Max Cartesian gradient
  REAL(DOUBLE),DIMENSION(4) :: GTol       =(/ 5.D-2, 5.D-3, 5.D-4,  5.D-5  /)
  ! Max Cartesian displacement
  REAL(DOUBLE),DIMENSION(4) :: XTol       =(/ 1.D-1, 1.D-2, 1.D-3,  1.D-4  /)
  ! Max error in the response.
  REAL(DOUBLE),DIMENSION(4) :: RTol       =(/ 1.D-4, 1.D-5, 1.D-6,  1.D-7 /)
END MODULE Numerics
