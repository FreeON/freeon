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

MODULE ParseCommands
  USE Parse
  USE InOut
#ifdef NAG
  USE F90_UNIX
#endif
  USE ControlStructures
  USE GlobalCharacters
  USE Utilities

  IMPLICIT NONE

CONTAINS
  !===============================================================================================
  ! PARSE THE COMMAND LINE AND GET RELATED ENV VARIABLES. CHECK EXISTENCE OF IN AND OUT FILES
  !===============================================================================================
  SUBROUTINE LoadCommands(N)
    TYPE(FileNames)     :: N
    TYPE(ARGMT)         :: Args
    CHARACTER(LEN=DCL)  :: PROCESS_ID
    INTEGER             :: DotDex
    LOGICAL             :: Exists
    INTEGER,EXTERNAL    :: GetPID

    ! Get command line arguments
    CALL Get(Args)
    IF(Args%NC==0) CALL MondoHalt(PRSE_ERROR,' No arguments to FreeON !')!

    ! Get current working directory.
    CALL GetPWD(N%M_PWD)

    ! Get environmental variables
    CALL GetEnv('MONDO_HOME',N%M_HOME)
    IF(LEN(TRIM(N%M_HOME)) == 0) THEN
      N%M_HOME = HAVE_MONDO_HOME
      CALL MondoLog(DEBUG_NONE,"FreeON", "env variable $(MONDO_HOME) not set. Using "//trim(N%M_HOME), "LoadCommand")
    ELSE
      CALL MondoLog(DEBUG_NONE, "FreeON", 'env variable $(MONDO_HOME) set to '//trim(N%M_HOME), "LoadCommand")
    ENDIF

    CALL GetEnv('MONDO_EXEC',N%M_EXEC)
    IF(LEN(TRIM(N%M_EXEC)) == 0) THEN
      N%M_EXEC = TRIM(N%M_HOME)//"/bin"
      CALL MondoLog(DEBUG_NONE, "FreeON", 'env variable $(MONDO_EXEC) not set. Using '//TRIM(N%M_EXEC), "LoadCommand")
    ELSE
      CALL MondoLog(DEBUG_NONE, "FreeON", 'env variable $(MONDO_EXEC) set to '//trim(N%M_EXEC), "LoadCommand")
    ENDIF

    CALL GetEnv('MONDO_SCRATCH',N%M_SCRATCH)
    IF(LEN(TRIM(N%M_SCRATCH)) == 0) THEN
      CALL MondoLog(DEBUG_NONE, "FreeON", 'env variable $(MONDO_SCRATCH) not set. Using '//TRIM(HAVE_MONDO_SCRATCH), "LoadCommand")
      N%M_SCRATCH = HAVE_MONDO_SCRATCH
      MONDO_SCRATCH = HAVE_MONDO_SCRATCH
    ELSE
      CALL MondoLog(DEBUG_NONE, "FreeON", 'env variable $(MONDO_SCRATCH) set to '//trim(N%M_SCRATCH), "LoadCommand")
    ENDIF

    ! Set path names etc
    N%M_PWD=TRIM(N%M_PWD)//'/'
    N%M_HOME=TRIM(N%M_HOME)//'/'

    ! Here is the scf name, tagged with the PID
    PROCESS_ID=IntToChar(GetPID())
    DotDex=INDEX(Args%C%C(1),'.')
    IF(DotDex==0) THEN
      CALL MondoLog(DEBUG_NONE, "FreeON", 'Parse error: no "." in input file name = <'//TRIM(N%IFile)//'>', "LoadCommand")
      STOP "Termination of FreeON"
    ENDIF

    CALL MondoLog(DEBUG_NONE, "FreeON", "running intut file "//TRIM(Args%C%C(1)))
    N%SCF_NAME=Args%C%C(1)(1:DotDex-1)//'_'//TRIM(PROCESS_ID)

    ! Come up with random scratch directory based on MONDO_SCRATCH.
    N%M_SCRATCH = TRIM(N%M_SCRATCH)//"/FreeON-scratch-"//TRIM(N%SCF_NAME)//"-XXXXXX"
    CALL TemporaryDirectory(N%M_SCRATCH, LEN(TRIM(N%M_SCRATCH)))
    N%M_SCRATCH=TRIM(N%M_SCRATCH)//'/'

    PWDName=TRIM(N%M_PWD)//TRIM(N%SCF_NAME)
    ScrName=TRIM(N%M_SCRATCH)//TRIM(N%SCF_NAME)

    ! Input file with full path
    N%IFile = TRIM(Args%C%C(1))
    IF(N%IFile(1:1) == '/') THEN
      ! This is an absolute path.
      N%IFile = TRIM(Args%C%C(1))
    ELSE
      ! Relative path, make it absolute.
      N%IFile=TRIM(N%M_PWD)//TRIM(Args%C%C(1))
    ENDIF

    ! Check to see that the input file exists
    INQUIRE(FILE=N%IFile,EXIST=Exists)
    IF(.NOT.Exists) THEN
      CALL MondoLog(DEBUG_NONE, "FreeON", 'Parse error: Input file "'//TRIM(N%IFile)//'" does not exist!', "LoadCommand")
      STOP "Termination of FreeON"
    ENDIF

    ! Determine if input file has a '.' in it.  If so, create SCF name from string up to the '.'
    DotDex=INDEX(Args%C%C(1),'.')
    IF(DotDex==0)THEN
      N%SCF_NAME=TRIM(Args%C%C(1))//'_'//TRIM(PROCESS_ID)
    ELSE
      N%SCF_NAME=Args%C%C(1)(1:DotDex-1)//'_'//TRIM(PROCESS_ID)
    ENDIF

    ! Out file with full path (demand an output file)
    ! Check to see that the output file does not exist (no overwrites allowed)
    !    INQUIRE(FILE=N%OFile,EXIST=Exists)
    !    IF(Exists) &
    !       CALL MondoHalt(PRSE_ERROR,' Ouput file: '//TRIM(N%OFile)//' already exists! ')

    ! Create user defined or implicit file names
    IF(Args%NC==1)THEN
      N%OFile=TRIM(PWDName)//OutF
      N%LFile=TRIM(PWDName)//LogF
      N%GFile=TRIM(PWDName)
    ELSEIF(Args%NC==2)THEN
      N%OFile=TRIM(N%M_PWD)//TRIM(Args%C%C(2))
      N%LFile=TRIM(PWDName)//LogF
      N%GFile=TRIM(PWDName)
    ELSEIF(Args%NC==3)THEN
      N%OFile=TRIM(N%M_PWD)//TRIM(Args%C%C(2))
      N%LFile=TRIM(N%M_PWD)//TRIM(Args%C%C(3))
      N%GFile=TRIM(PWDName)
    ELSEIF(Args%NC==4)THEN
      N%OFile=TRIM(N%M_PWD)//TRIM(Args%C%C(2))
      N%LFile=TRIM(N%M_PWD)//TRIM(Args%C%C(3))
      N%GFile=TRIM(N%M_PWD)//TRIM(Args%C%C(4))
    ENDIF
    N%HFile=TRIM(ScrName)//InfF

    ! If restart file set, it is set in ParseOptions
    N%RFile=""

  END SUBROUTINE LoadCommands
END MODULE ParseCommands
