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

#if defined(PARALLEL) || defined(PARALLEL_CLONES)
  USE MPI
#endif

  IMPLICIT NONE

CONTAINS

  !===============================================================================================
  ! PARSE THE COMMAND LINE AND GET RELATED ENV VARIABLES. CHECK EXISTENCE OF IN AND OUT FILES
  !===============================================================================================
  SUBROUTINE LoadCommands(N)
    TYPE(FileNames)     :: N
    TYPE(ARGMT)         :: Args
    CHARACTER(LEN=DCL)  :: PROCESS_ID
    INTEGER             :: i, indexBegin, indexEnd
    LOGICAL             :: Exists

    CALL MondoLog(DEBUG_NONE, "FreeON", "Version "//TRIM(PACKAGE_VERSION)//" starting...")

    ! Get command line arguments
    CALL Get(Args)
    IF(Args%NC==0) CALL MondoHalt(PRSE_ERROR,' No arguments to FreeON !')!

    ! Get current working directory.
    CALL GetPWD(N%M_PWD)

    ! Get environmental variables
    CALL GetEnv('FREEON_HOME',N%M_HOME)
    IF(LEN(TRIM(N%M_HOME)) == 0) THEN
      N%M_HOME = FREEON_HOME
      CALL MondoLog(DEBUG_NONE,"FreeON", "env variable $(FREEON_HOME) not set. Using "//trim(N%M_HOME), "LoadCommand")
    ELSE
      CALL MondoLog(DEBUG_NONE, "FreeON", 'env variable $(FREEON_HOME) set to '//trim(N%M_HOME), "LoadCommand")
    ENDIF

    CALL GetEnv('FREEON_BASISSETS',N%M_BASISSETS)
    IF(LEN(TRIM(N%M_BASISSETS)) == 0) THEN
      N%M_BASISSETS = FREEON_BASISSETS
      CALL MondoLog(DEBUG_NONE,"FreeON", "env variable $(FREEON_BASISSETS) not set. Using "//trim(N%M_BASISSETS), "LoadCommand")
    ELSE
      CALL MondoLog(DEBUG_NONE, "FreeON", 'env variable $(FREEON_BASISSETS) set to '//trim(N%M_BASISSETS), "LoadCommand")
    ENDIF

    CALL GetEnv('MONDO_EXEC',N%M_EXEC)
    IF(LEN(TRIM(N%M_EXEC)) == 0) THEN
      N%M_EXEC = TRIM(N%M_HOME)//"/bin"
      CALL MondoLog(DEBUG_NONE, "FreeON", 'env variable $(MONDO_EXEC) not set. Using '//TRIM(N%M_EXEC), "LoadCommand")
    ELSE
      CALL MondoLog(DEBUG_NONE, "FreeON", 'env variable $(MONDO_EXEC) set to '//trim(N%M_EXEC), "LoadCommand")
    ENDIF

    CALL GetEnv('FREEON_SCRATCH',N%M_SCRATCH)
    IF(LEN(TRIM(N%M_SCRATCH)) == 0) THEN
      CALL MondoLog(DEBUG_NONE, "FreeON", 'env variable $(FREEON_SCRATCH) not set. Using '//TRIM(HAVE_FREEON_SCRATCH), "LoadCommand")
      N%M_SCRATCH = HAVE_FREEON_SCRATCH
      FREEON_SCRATCH = HAVE_FREEON_SCRATCH
    ELSE
      CALL MondoLog(DEBUG_NONE, "FreeON", 'env variable $(FREEON_SCRATCH) set to '//trim(N%M_SCRATCH), "LoadCommand")
    ENDIF

    ! Set path names etc
    N%M_PWD=TRIM(N%M_PWD)//'/'
    N%M_HOME=TRIM(N%M_HOME)//'/'

    ! Here is the scf name, tagged with the PID.
    PROCESS_ID = IntToChar(GetPIDWrapper())

    CALL MondoLog(DEBUG_NONE, "FreeON", "input file "//TRIM(Args%C%C(1)), "LoadCommand")

    ! For the SCF name, we remove everything up to the last "/" in case the
    ! input file is a path.  If the filename contains a ".", we remove the last
    ! suffix. Other "." are not touched. Then we append the PID to the name.
    indexBegin = 1
    indexEnd = LEN(TRIM(Args%C%C(1)))
    N%SCF_NAME = TRIM(Args%C%C(1))

    IF(INDEX(N%SCF_NAME, "/", .TRUE.) > 0) THEN
      indexBegin = INDEX(Args%C%C(1), "/", .TRUE.)+1
    ENDIF

    IF(INDEX(N%SCF_NAME, ".", .TRUE.) > 0) THEN
      indexEnd = INDEX(Args%C%C(1), ".", .TRUE.)-1
    ENDIF

    N%SCF_NAME = TRIM(N%SCF_NAME(indexBegin:indexEnd))//'_'//TRIM(PROCESS_ID)
    CALL MondoLog(DEBUG_NONE, "FreeON", "SCF name "//TRIM(N%SCF_NAME), "LoadCommand")

    ! Come up with random scratch directory based on FREEON_SCRATCH.
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
      CALL MondoLog(DEBUG_NONE, "FreeON", 'Parse error: input file "'//TRIM(N%IFile)//'" does not exist!', "LoadCommand")
      STOP "Termination of FreeON"
    ENDIF

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
