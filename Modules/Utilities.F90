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

MODULE Utilities

  USE MondoLogger
  USE ParsingConstants
  USE Parse

  IMPLICIT NONE

  INTERFACE

    SUBROUTINE GetMemoryUsage_Wrapper (virtualMemorySize)
      INTEGER :: virtualMemorySize
    END SUBROUTINE GetMemoryUsage_Wrapper

    SUBROUTINE GetHostnameWrapper (hostname, maximumLength)
      CHARACTER(LEN = *) :: hostname
      INTEGER            :: maximumLength
    END SUBROUTINE GetHostnameWrapper

    SUBROUTINE GetStacksizeLimit (currentLimit, maximumLimit)
      INTEGER :: currentLimit, maximumLimit
    END SUBROUTINE GetStacksizeLimit

    SUBROUTINE FileCopyWrapper (lenA, fileA, lenB, fileB)
      INTEGER            :: lenA, lenB
      CHARACTER(LEN = *) :: fileA, fileB
    END SUBROUTINE FileCopyWrapper

    SUBROUTINE GetPWDWrapper (pwd, max_length)
      INTEGER            :: max_length
      CHARACTER(LEN = *) :: pwd
    END SUBROUTINE GetPWDWrapper

    SUBROUTINE TemporaryDirectory (path, length)
      CHARACTER(LEN = *), INTENT(INOUT) :: path
      INTEGER, INTENT(IN)               :: length
    END SUBROUTINE TemporaryDirectory

    FUNCTION GetPIDWrapper ()
      INTEGER :: GetPIDWrapper
    END FUNCTION GetPIDWrapper

  END INTERFACE

  ! The parameter sleeptime is in seconds. It can be fractions of a second, down
  ! to a microsecond.
  INTERFACE FreeONSleep

    SUBROUTINE FreeONSleep_integer (sleeptime)
      INTEGER, INTENT(IN) :: sleeptime
    END SUBROUTINE FreeONSleep_integer

    SUBROUTINE FreeONSleep_single (sleeptime)
      USE GlobalScalars
      REAL(SINGLE), INTENT(IN) :: sleeptime
    END SUBROUTINE FreeONSleep_single

  END INTERFACE FreeONSleep


CONTAINS

  SUBROUTINE NOP ()
    INTEGER :: i
    i = 0
  END SUBROUTINE NOP

  FUNCTION GetMemoryUsage ()
    INTEGER :: GetMemoryUsage
    INTEGER :: virtualMemorySize

    CALL GetMemoryUsage_Wrapper(virtualMemorySize)
    GetMemoryUsage = virtualMemorySize

    RETURN

  END FUNCTION GetMemoryUsage

  SUBROUTINE FileCopy (fileA, fileB)
    CHARACTER(LEN = *) :: fileA, fileB
    CALL MondoLog(DEBUG_MAXIMUM, "FileCopy", "copying "//TRIM(fileA)//" --> "//TRIM(FileB))
    CALL FileCopyWrapper(LEN(TRIM(fileA)), TRIM(fileA), LEN(TRIM(fileB)), TRIM(fileB))
  END SUBROUTINE FileCopy

  SUBROUTINE FileRemove (filename)
    CHARACTER(LEN = *), INTENT(IN) :: filename
    CALL MondoLog(DEBUG_MAXIMUM, "FileRemove", "removing "//TRIM(filename)// &
      " -> "//TRIM(EscapeFilename(filename)))
    CALL SYSTEM("rm -f "//TRIM(EscapeFilename(filename)))
  END SUBROUTINE

  SUBROUTINE GetPWD (pwd)
    CHARACTER(LEN = *) :: pwd
    CALL GetPWDWrapper(pwd, LEN(pwd))
  END SUBROUTINE GetPWD

END MODULE Utilities
