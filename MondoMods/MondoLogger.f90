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

  PUBLIC :: MondoLog, MondoLogPlain, ProcessName

CONTAINS

  FUNCTION ProcessName(Proc_O,Misc_O) RESULT (Tag)
    CHARACTER(LEN=*), OPTIONAL :: Proc_O
    CHARACTER(LEN=*), OPTIONAL :: Misc_O
    CHARACTER(LEN=32)          :: Tag
    CHARACTER(LEN=32)          :: Name
    CHARACTER(LEN=32),PARAMETER:: Blks='                        '
    CHARACTER(LEN=3), PARAMETER:: Colon =' : '
    CHARACTER(LEN=4), PARAMETER:: Colons=' :: '
    IF(PRESENT(Proc_O).AND.PRESENT(Misc_O))THEN
      Name=TRIM(ADJUSTL(TRIM(Proc_O)))//Colon//TRIM(Misc_O)
      Name(30:32)=":: "
      Tag=Name
    ELSEIF(PRESENT(Proc_O))THEN
      Name=TRIM(ADJUSTL(TRIM(Proc_O)))
      Name(30:32)=":: "
      Tag=Name
    ELSE
      Tag="" !Blks
    ENDIF
  END FUNCTION ProcessName

  ! Open the logfile.
  FUNCTION OpenLogfile(filename, fd) RESULT(fileOutput)
    CHARACTER(LEN=*), INTENT(IN)    :: filename
    INTEGER, INTENT(IN)             :: fd
    INTEGER                         :: IOS, fd_old
    LOGICAL                         :: isOpen, exists, fileOutput
    CHARACTER(LEN=INTERNAL_INT_LEN) :: fd_string

    ! Set the default result.
    fileOutput = .FALSE.

    ! Does the file exist?
    INQUIRE(FILE = filename, OPENED = isOpen, EXIST = exists, &
         ERR = 11, IOSTAT = IOS, NUMBER = fd_old)

    IF(isOpen) THEN
      ! The log file was opened somewhere else in the code and not closed
      ! properly. This indicates a possible bug.
      CLOSE(fd_old)
      isOpen = .FALSE.

      ! Log this also. Now that we have closed the log file we can use
      ! MondoLogger itself to properly log this event.
      WRITE(UNIT=fd_string, FMT="(I3)") fd_old
      fd_string = ADJUSTL(fd_string)
      CALL MondoLog(DEBUG_NONE, "OpenLogfile", "logfile already open (fd(old) = "//TRIM(fd_string)//")")
    ENDIF

    IF(exists.AND.(.NOT.isOpen)) THEN

      ! Open existing file and position at the bottom (default)
      OPEN(UNIT = fd, FILE = filename, &
           ACCESS = "SEQUENTIAL", FORM = "FORMATTED", &
           POSITION = "APPEND", ERR = 12, IOSTAT = IOS, &
           STATUS = "OLD")

    ELSE

      ! Create a new file and open it
      !WRITE(*,"(A)") "[MondoLog] logfile '"//TRIM(filename)//"' does not exist, creating it"
      OPEN(UNIT = fd, FILE = filename, &
           ACCESS = "SEQUENTIAL", FORM = "FORMATTED", &
           ERR = 13, IOSTAT = IOS, STATUS = "NEW")

    ENDIF
    fileOutput = .TRUE.
    RETURN

11  WRITE(*,"(A)")    "[MondoLog.inquire] Fatal Error"
    WRITE(*,"(A,I3)") "  IOS  = ", IOS
    WRITE(*,"(A)")    "  file = "//TRIM(filename)
    CALL Trap()

12  WRITE(*,"(A)") "[MondoLog.append] Fatal Error"
    WRITE(*,"(A,I3)") "  IOS  = ", IOS
    WRITE(*,"(A)")    "  file = "//TRIM(filename)
    CALL Trap()
    !WRITE(*,"(A)") "[MondoLog.new] Fatal Error"
    !WRITE(*,"(A,I3)") "  IOS  = ", IOS
    !WRITE(*,"(A)")    "  file = "//TRIM(filename)
13  RETURN

  END FUNCTION OpenLogfile

  SUBROUTINE MondoLog(logLevel, tag, message, file_O, line_O, NoIndent_O)

    CHARACTER(LEN=*), INTENT(IN)            :: message
    CHARACTER(LEN=*), INTENT(IN)            :: tag
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: file_O
    CHARACTER(LEN=2048)                     :: output
    CHARACTER(LEN=2048)                     :: line_string
    INTEGER, OPTIONAL, INTENT(IN)           :: line_O
    LOGICAL, OPTIONAL, INTENT(IN)           :: NoIndent_O
    INTEGER :: logLevel
    LOGICAL :: isOpen, fileOutput

    ! Check whether logLevel is sufficiently high.

    IF(logLevel <= PrintFlags%Key) THEN
      ! Open the logfile.
      fileOutput = OpenLogfile(OutFile, 123)

      ! Convert the line number into string.
      IF(PRESENT(line_O)) THEN
        WRITE(UNIT=line_string,FMT="(I20)") line_O
        line_string = ADJUSTL(line_string)
        IF(PRESENT(NoIndent_O))THEN
          output=line_string
        ELSE
          output=ProcessName(tag,line_string)
        ENDIF

      ELSEIF(PRESENT(file_O)) THEN
        IF(PRESENT(NoIndent_O))THEN
          output=line_string
        ELSE
          output=ProcessName(tag,file_O)
        ENDIF
      ELSE
        IF(PRESENT(NoIndent_O))THEN
          output=tag
        ELSE
          output=ProcessName(tag)
        ENDIF
      ENDIF

      IF(LEN(TRIM(output)) > 0) THEN
        output = TRIM(output)//' '//TRIM(message)
      ELSE
        output = TRIM(message)
      ENDIF

      IF(fileOutput) THEN
        WRITE(123,"(A)") TRIM(output)
      ENDIF

      WRITE(*,"(A)") TRIM(output)

      ! Close logfile.
      IF(fileOutput) CLOSE(123)

    ENDIF

  END SUBROUTINE MondoLog

  SUBROUTINE MondoLogPlain(message)

    CHARACTER(LEN=*), INTENT(IN) :: message

    CALL MondoLog(DEBUG_NONE, "",message, NoIndent_O=.TRUE.)

  END SUBROUTINE MondoLogPlain

END MODULE MondoLogger
