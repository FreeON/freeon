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

! This module provides interfaces to file operation related C-wrapper programs.
!
! Author: Nicolas Bock <nbock@lanl.gov>

MODULE FileOperations

  IMPLICIT NONE

  INTERFACE FileCopyWrapper

    SUBROUTINE FileCopyWrapper (lenA, fileA, lenB, fileB)
      INTEGER :: lenA, lenB
      CHARACTER(LEN=*) :: fileA, fileB
    END SUBROUTINE FileCopyWrapper

  END INTERFACE FileCopyWrapper

  CONTAINS

  SUBROUTINE FileCopy (fileA, fileB)
    CHARACTER(LEN=*) :: fileA, fileB
    CALL FileCopyWrapper(LEN(TRIM(fileA)), TRIM(fileA), LEN(TRIM(fileB)), TRIM(fileB))
  END SUBROUTINE FileCopy

END MODULE FileOperations