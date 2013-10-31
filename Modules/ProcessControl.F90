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

MODULE ProcessControl
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE ParsingConstants
  USE MondoLogger

#if defined(PARALLEL) || defined(PARALLEL_CLONES)
  USE MPI
#endif

  IMPLICIT NONE

  !---------------------------------------------------------------------------------------------
  !  Status keys
  !
  INTEGER, PARAMETER     :: EXIT_ERROR=-120384
  INTEGER, PARAMETER     :: FORK_ERROR=-320498
  INTEGER, PARAMETER     :: DUMP_ERROR=-580234
  INTEGER, PARAMETER     :: SGNL_ERROR=-674034
  INTEGER, PARAMETER     :: MISC_ERROR=-843503
  INTEGER, PARAMETER     :: DRIV_ERROR=-704823
  INTEGER, PARAMETER     :: PRSE_ERROR=-803484
  INTEGER, PARAMETER     :: NEBS_ERROR=-924723
  INTEGER, PARAMETER     :: MPIS_ERROR=-975239
  INTEGER, PARAMETER     :: USUP_ERROR=-993942
  INTEGER, PARAMETER     :: INTC_ERROR=-135950
  INTEGER, PARAMETER     :: QMMM_ERROR=-458391
CONTAINS
  SUBROUTINE MondoHalt(IErr,Mssg)
    CHARACTER (LEN=*) :: Mssg
    INTEGER           :: IErr
    CHARACTER(LEN=10) :: stringIErr

    IF(IErr==EXIT_ERROR)THEN
      CALL MondoLogPlain(TRIM(Mssg)//' failed.')
    ELSEIF(IErr==DUMP_ERROR)THEN
      CALL MondoLogPlain(TRIM(Mssg)//' core dumped.')
    ELSEIF(IErr==SGNL_ERROR)THEN
      CALL MondoLogPlain(TRIM(Mssg)//' signaled.')
    ELSEIF(IErr==FORK_ERROR)THEN
      CALL MondoLogPlain('Unable to fork '//TRIM(Mssg))
    ELSEIF(IErr==DRIV_ERROR)THEN
      CALL MondoLogPlain('Error in SCFDriver: '//TRIM(Mssg))
    ELSEIF(IErr==PRSE_ERROR)THEN
      CALL MondoLogPlain('Error parsing: '//TRIM(Mssg))
    ELSEIF(IErr==MPIS_ERROR)THEN
      CALL MondoLogPlain('MPI error: '//TRIM(Mssg))
    ELSEIF(IErr==USUP_ERROR)THEN
      CALL MondoLogPlain('Unsupported feature: '//TRIM(Mssg))
    ELSEIF(IErr==NEBS_ERROR)THEN
      CALL MondoLogPlain('Error in NEB: '//TRIM(Mssg))
    ELSEIF(IErr==INTC_ERROR)THEN
      CALL MondoLogPlain('Error in IntCoo: '//TRIM(Mssg))
    ELSEIF(IErr==QMMM_ERROR)THEN
      CALL MondoLogPlain('Error in QMMM  : '//TRIM(Mssg))
    ELSE
      WRITE(stringIErr, *) IErr
      CALL MondoLogPlain("Unknown error (error code "//TRIM(stringIErr)//"): "//TRIM(Mssg))
    ENDIF
    WRITE(*,*)' '
    WRITE(*,*)' To have no errors '
    WRITE(*,*)' Would be life without meaning '
    WRITE(*,*)' No struggle, no joy '
    WRITE(*,*)' '
    STOP 'Termination of FreeON'
  END SUBROUTINE MondoHalt

  SUBROUTINE Halt(Strng)
    CHARACTER(LEN=*)           :: Strng
    CHARACTER(LEN=*),PARAMETER :: Motto='Frango ut patefaciam -- I break in order to reveal'
    CALL MondoLogPlain(Strng)
#if defined(PARALLEL) || defined(PARALLEL_CLONES)
    IF(InParallel)THEN
      WRITE(*,*)'Node that brought me down = ',MyId
      CALL Trap()
      CALL HaltMPI(Strng)
    ELSE
      WRITE(*,*)'GOING DOWN GOING DOWN GOING DOWN GOING DOWN GOING DOWN '
      WRITE(*,*)'<<'//TRIM(Strng)//'>>'
      WRITE(*,*)'GOING DOWN GOING DOWN GOING DOWN GOING DOWN GOING DOWN '
      CALL Trap()
    ENDIF
#else
    WRITE(*,*)'GOING DOWN GOING DOWN GOING DOWN GOING DOWN GOING DOWN '
    WRITE(*,*)'<<'//TRIM(Strng)//'>>'
    WRITE(*,*)'GOING DOWN GOING DOWN GOING DOWN GOING DOWN GOING DOWN '
    CALL Trap()
#endif
  END SUBROUTINE Halt

#if defined(PARALLEL) || defined(PARALLEL_CLONES)
  SUBROUTINE HaltMPI(InMsg,ErrorCode)
    CHARACTER(LEN=*), INTENT(IN)        :: InMsg
    INTEGER, OPTIONAL,INTENT(IN)        :: ErrorCode
    INTEGER                             :: IErr
    CHARACTER(LEN=MPI_MAX_ERROR_STRING) :: ErrMsg
    CHARACTER(LEN=2*DEFAULT_CHR_LEN)    :: LogMssg
    CHARACTER(LEN=INTERNAL_INT_LEN)     :: ChId
    INTEGER                             :: LenMsg
    WRITE(ChId,INTERNAL_INT_FMT)MyId
    ChId=ADJUSTL(ChId)
    LogMssg='GOING, GOING,... GONE DOWN: MyID = '//TRIM(ChId) &
         //Rtrn//'Message that killed me: '//TRIM(InMsg)
    IF(PRESENT(ErrorCode))THEN
      CALL MPI_Error_string(ErrorCode,ErrMsg,LenMsg,IErr)
      LogMssg=ADJUSTL(TRIM(LogMssg))//Rtrn &
           //'MPI ERROR = <'//TRIM(ErrMsg(1:LenMsg))//'>'
      ! CALL LogInTurn(LogMssg)
      CALL MondoLogPlain(LogMssg)
      CALL MPI_Abort(MPI_COMM_WORLD,ErrorCode,IErr)
    ELSE
      ! CALL LogInTurn(LogMssg)
      CALL MondoLogPlain(LogMssg)
      CALL MPI_Abort(MPI_COMM_WORLD,ErrorCode,IErr)
    ENDIF
  END SUBROUTINE HaltMPI

  SUBROUTINE LogInTurn(Mssg,IErr)
    CHARACTER(LEN=*),  INTENT(IN)  :: Mssg
    INTEGER, OPTIONAL, INTENT(OUT) :: IErr
    INTEGER                        :: I
    DO I=0,NPrc-1
      CALL MPI_BARRIER(MPI_COMM_WORLD,IErr)
      IF(MyId==I) CALL MondoLogPlain(Mssg)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IErr)
    ENDDO
  END SUBROUTINE LogInTurn
#endif

  SUBROUTINE Warn(Strng)
    CHARACTER (LEN=*) Strng
    CHARACTER (LEN=10*DEFAULT_CHR_LEN) :: Warning
    Warning='... ning Warning  Warning  Warning Wa ...'//Rtrn//TRIM(ADJUSTL(Strng))
    CALL MondoLogPlain(TRIM(Warning))
  END SUBROUTINE Warn

END MODULE ProcessControl
