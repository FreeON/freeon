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

! This module provides logger functions to print "stuff" into some log file or
! the console.
!
! Author: Nicolas Bock <nbock@lanl.gov>

MODULE MondoLogger

  USE GlobalCharacters
  USE GlobalObjects
  USE ParsingConstants

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MondoLog, MondoLogPlain

CONTAINS

  ! Open the logfile.
  SUBROUTINE OpenLogfile(filename, fd)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER                      :: IOS, fd, fd_old
    LOGICAL                      :: isOpen, exists

    ! Does the file exist?
    INQUIRE(FILE = filename, OPENED = isOpen, EXIST = exists, &
      ERR = 11, IOSTAT = IOS, NUMBER = fd_old)

    IF(isOpen) THEN
      WRITE(*,"(A,I3,A)") "[MondoLog] logfile already open (fd = ", fd_old, ")"
      IF(fd /= fd_old) THEN
        WRITE(*,"(A)") "  WARNING: fd != fd_old"
      ENDIF
      CLOSE(fd_old)
      isOpen = .FALSE.
    ENDIF

    IF(exists.AND.(.NOT.isopen)) THEN

      ! Open existing file and position at the bottom (default)
      OPEN(UNIT = fd, FILE = filename, &
        ACCESS = "SEQUENTIAL", FORM = "FORMATTED", &
        POSITION = "APPEND", ERR = 12, IOSTAT = IOS, &
        STATUS = "OLD")

    ELSE

      ! Create a new file and open it
      WRITE(*,*) "[MondoLog] logfile does not exist, creating it"
      OPEN(UNIT = fd, FILE = filename, &
        ACCESS = "SEQUENTIAL", FORM = "FORMATTED", &
        ERR = 13, IOSTAT = IOS, STATUS = "NEW")

    ENDIF
    RETURN

11  WRITE(*,"(A)")    "[MondoLog.inquire] Fatal Error"
    WRITE(*,"(A,I3)") "  IOS  = ", IOS
    WRITE(*,"(A)")    "  file = "//TRIM(filename)
    STOP

12  WRITE(*,"(A)") "[MondoLog.append] Fatal Error"
    WRITE(*,"(A,I3)") "  IOS  = ", IOS
    WRITE(*,"(A)")    "  file = "//TRIM(filename)
    STOP

13  WRITE(*,"(A)") "[MondoLog.new] Fatal Error"
    WRITE(*,"(A,I3)") "  IOS  = ", IOS
    WRITE(*,"(A)")    "  file = "//TRIM(filename)
    STOP

  END SUBROUTINE OpenLogfile

  SUBROUTINE MondoLog(logLevel, tag, message)

    CHARACTER(LEN=*), INTENT(IN) :: message
    CHARACTER(LEN=*), INTENT(IN) :: tag
    INTEGER :: logLevel
    LOGICAL :: isOpen

    ! Check whether logLevel is sufficiently high.
    IF(logLevel <= PrintFlags%Key) THEN

      ! Open the logfile.
      CALL OpenLogfile(OutFile, 123)

      ! Write messsage.
      IF(LEN_TRIM(tag) > 0) THEN
        WRITE(123,"(A,A,A,A)") "[", tag, "] ", message
        WRITE(*,"(A,A,A,A)") "[", tag, "] ", message
      ELSE
        WRITE(123,"(A)") message
        WRITE(*,"(A)") message
      ENDIF

      ! Close logfile.
      CLOSE(123)

    ENDIF

  END SUBROUTINE MondoLog

  SUBROUTINE MondoLogPlain(message)

    CHARACTER(LEN=*), INTENT(IN) :: message

    CALL MondoLog(DEBUG_NONE, "", message)

  END SUBROUTINE MondoLogPlain

END MODULE MondoLogger