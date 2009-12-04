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

MODULE Overlay
  USE InOut
  USE PunchHDF
  USE ControlStructures
  USE Utilities

  IMPLICIT NONE

CONTAINS

  ! Invoke the executable Ex first using the absolute M_EXEC path, and if that
  ! fails without any path, relying on the PATH environment variable and the
  ! shell.
  SUBROUTINE Invoke(Ex, N, S, M, TSSearchEndpoints_O)
    CHARACTER(LEN=*)     :: Ex
    TYPE(FileNames)      :: N
    TYPE(State)          :: S
    TYPE(Parallel)       :: M
    LOGICAL, OPTIONAL    :: TSSearchEndpoints_O

    CALL MondoLog(DEBUG_NONE, "Invoke", "invoking "//TRIM(Ex))
    IF(InvokeBackend(TRIM(N%M_EXEC)//"/"//Ex, N, S, M, TSSearchEndpoints_O) /= 0) THEN
      IF(InvokeBackend(Ex, N, S, M, TSSearchEndpoints_O) /= 0) THEN
        CALL MondoLog(DEBUG_NONE, "Invoke", "failed to spawn process")
        CALL Halt("This is fatal")
      ENDIF
    ENDIF

  END SUBROUTINE Invoke

  FUNCTION InvokeBackend(Ex, N, S, M, TSSearchEndpoints_O)
    INTEGER              :: InvokeBackend
    CHARACTER(LEN=*)     :: Ex
    TYPE(FileNames)      :: N
    TYPE(State)          :: S
    TYPE(Parallel)       :: M
    LOGICAL, OPTIONAL    :: TSSearchEndpoints_O

    INTEGER              :: I,J,K,L,NC,iCLUMP,IErr,NArg,MaxLen
    LOGICAL              :: ProgramFailed
    TYPE(CHR_VECT),SAVE  :: ArgV
    TYPE(INT_VECT),SAVE  :: IChr
    CHARACTER(LEN=2*DCL) :: CmndLine
    INTEGER              :: beginClump, endClump

#if (defined(PARALLEL) || defined(PARALLEL_CLONES)) && defined(MPI2)
    INTEGER                              :: SPAWN, INTRA_SPAWN, buffer_index
    CHARACTER, ALLOCATABLE, DIMENSION(:) :: message_buffer
    INTEGER, ALLOCATABLE, DIMENSION(:)   :: message_request
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: message_status
    LOGICAL, ALLOCATABLE, DIMENSION(:)   :: message_flag
#else
    INTERFACE
      FUNCTION Spawn(NC, MaxLen, IChr)
        INTEGER                         :: NC,MaxLen
        INTEGER, DIMENSION(1:NC*MaxLen) :: IChr
        INTEGER                         :: Spawn
      END FUNCTION Spawn
    END INTERFACE
#endif

    ! Set default return value (failure).
    InvokeBackend = -1

    IF(PRESENT(TSSearchEndpoints_O) .AND. TSSearchEndpoints_O) THEN
      ! Calculate the clones and the endpoint energies.
#ifdef PARALLEL
      CALL Halt("[FIXME] This is not implemented")
#else
      beginClump = 0
      endClump = M%Clumps+1
#endif
    ELSE
      ! Calculate the clones only.
      beginClump = 1
      endClump = M%Clumps
    ENDIF

    DO iCLUMP = beginClump, endClump

#if defined(PARALLEL) || defined(PARALLEL_CLONES)
      CALL MPIsArchive(N,M%NSpace,M%Clump%I(:,iCLUMP))
      CALL SetArgV(Ex,N,S,M,iCLUMP,NArg,ArgV)
#else
      CALL MPIsArchive(N,1,M%Clump%I(:,iCLUMP))
      CALL SetArgV(Ex,N,S,M,NArg,ArgV)
#endif

      ! This is the command line we are going to execute
      CmndLine=TRIM(ArgV%C(1))
      DO I=2,NArg
        CmndLine=TRIM(CmndLine)//Blnk//TRIM(ArgV%C(I))
      ENDDO

      ! Log this run
      CALL MondoLog(DEBUG_MAXIMUM, "Invoke", TRIM(CmndLine), "Clump "//TRIM(IntToChar(iCLUMP)))

#if (defined(PARALLEL) || defined(PARALLEL_CLONES)) && defined(MPI2)
      CALL MPI_COMM_SPAWN(ArgV%C(1),ArgV%C(2:NArg),M%Clump%I(3, iCLUMP),MPI_INFO_NULL,ROOT,MPI_COMM_SELF,SPAWN,MPI_ERRCODES_IGNORE,IErr)
      IF(IErr/=MPI_SUCCESS) THEN
        CALL MondoLog(DEBUG_NONE, "Invoke", "Could not spawn <"//TRIM(CmndLine)//">", "errorcode = "//TRIM(IntToChar(MPIS_ERROR)))
        CALL Delete(ArgV)
        RETURN
      ENDIF

      ! Merge the spawned and current local communicators
      CALL MPI_INTERCOMM_MERGE(SPAWN, .TRUE., INTRA_SPAWN, IErr)

      ! Wait for the kiddies to be done. We do that by waiting for them to send
      ! us a nice message instead of waiting at a barrier. This is to avoid
      ! spending CPU time in the front-end polling the barrier. Waiting for
      ! messages, we can do non-blocking receives and control the polling
      ! frequency ourselves.
      !
      ! First set up the receive buffers.
      ALLOCATE(message_buffer(M%Clump%I(3, iCLUMP)))
      ALLOCATE(message_request(M%Clump%I(3, iCLUMP)))
      ALLOCATE(message_flag(M%Clump%I(3, iCLUMP)))
      ALLOCATE(message_status(M%Clump%I(3, iCLUMP), MPI_STATUS_SIZE))

      DO I = 1, M%Clump%I(3, iCLUMP)
        !CALL MondoLog(DEBUG_NONE, "Invoke", "setting up receive buffers")
        CALL MPI_IRECV(message_buffer(I), 1, MPI_CHARACTER, MPI_ANY_SOURCE, 0, SPAWN, message_request(I), IErr)
        message_flag(I) = .FALSE.
      ENDDO

      ! Loop over children to check whether they are done and have sent a
      ! message.
      buffer_index = 0
      DO WHILE(buffer_index < M%Clump%I(3, iCLUMP))
        DO I = 1, M%Clump%I(3, iCLUMP)
          IF(.NOT. message_flag(I)) THEN
            !CALL MondoLog(DEBUG_NONE, "Invoke", "testing child "//TRIM(IntToChar(I)))
            CALL MPI_TEST(message_request(I), message_flag(I), message_status(I,:), IErr)
            IF(message_flag(I)) THEN
              !CALL MondoLog(DEBUG_NONE, "Invoke", "child "//TRIM(IntToChar(buffer_index+1))//" has finished")
              buffer_index = buffer_index+1
            ENDIF
          ENDIF
        ENDDO

        ! Sleep a little.
        !CALL MondoLog(DEBUG_NONE, "Invoke", "sleeping")
        CALL FreeONSleep(2)
        !CALL MondoLog(DEBUG_NONE, "Invoke", "done sleeping")
      ENDDO

      ! All children are done.
      !CALL MondoLog(DEBUG_NONE, "Invoke", "all children are done")

      ! Free the communicator.
      CALL MPI_COMM_FREE(SPAWN, IErr)
      CALL MPI_COMM_FREE(INTRA_SPAWN, IErr)

      ! Free the buffers.
      DEALLOCATE(message_buffer)
      DEALLOCATE(message_request)
      DEALLOCATE(message_flag)
      DEALLOCATE(message_status)
#else
      ! Create ASCII integer array to beat F9x/C incompatibility
      CALL CVToIV(NArg,ArgV,MaxLen,IChr)

      ! Spawn a sub process
      IErr=Spawn(NArg,MaxLen,IChr%I)
      CALL Delete(IChr)

      ! Bring this run down if not successful
      IF(IErr /= SUCCEED) THEN
        CALL MondoLog(DEBUG_NONE, "Invoke", "<"//TRIM(CmndLine)//">", "errorcode = "//TRIM(IntToChar(IErr)))
        CALL Delete(ArgV)
        RETURN
      ENDIF

      ! Double check success if a MONDO Exec ...
      HDF_CurrentID=OpenHDF(N%HFile)
      CALL Get(ProgramFailed,'ProgramFailed')
      CALL CloseHDF(HDF_CurrentID)
      IF(ProgramFailed) THEN
        CALL MondoHalt(-999,'<'//TRIM(CmndLine)//'>')
      ENDIF
#endif
      CALL Delete(ArgV)
    ENDDO

    ! Return success.
    InvokeBackend = 0
    RETURN

  END FUNCTION InvokeBackend

  !===============================================================
  ! CREATE A CHARACTER ARRAY OF NON-BLANK STRINGS THAT WILL
  ! BECOME THE ARGV ARRAY PASSED TO EXECVP BY SPAWN IF NOT MPI-2
  !===============================================================
#if defined(PARALLEL) || defined(PARALLEL_CLONES)
  SUBROUTINE SetArgV(Ex,N,S,M,cCLUMP,NArg,ArgV)
#else
  SUBROUTINE SetArgV(Ex,N,S,M,NArg,ArgV)
#endif
    CHARACTER(LEN=*)   :: Ex
    TYPE(FileNames)    :: N
    TYPE(State)        :: S
    TYPE(Parallel)     :: M
    INTEGER            :: I,K,NArg,cCLUMP,SNC,NewDex
    TYPE(CHR_VECT)     :: ArgT,ArgV

    ! Start...
    SNC=SIZE(S%Action%C)

#ifdef MPI2
    NArg=9+SNC
    CALL New(ArgT,NArg)
    ArgT%C(1)=Ex
    ArgT%C(2)=N%SCF_NAME
    NewDex = 2
    DO K=1,SNC
      ArgT%C(NewDex+K)=S%Action%C(K)
    ENDDO
    NewDex=2+SNC
    ArgT%C(NewDex+1)=IntToChar(S%Current%I(1))
    ArgT%C(NewDex+2)=IntToChar(S%Current%I(2))
    ArgT%C(NewDex+3)=IntToChar(S%Current%I(3))
    ArgT%C(NewDex+4)=IntToChar(S%Previous%I(1))
    ArgT%C(NewDex+5)=IntToChar(S%Previous%I(2))
    ArgT%C(NewDex+6)=IntToChar(S%Previous%I(3))
    ArgT%C(NewDex+7)=TRIM(N%M_SCRATCH)
#elif PARALLEL
    NArg=14+SNC
    CALL New(ArgT,NArg)
    ArgT%C(1)=M%Invoking
    ArgT%C(2)=M%ProcFlag
    ArgT%C(3)=IntToChar(M%Clump%I(3,cCLUMP))
    ArgT%C(4)=M%MachFlag
    ArgT%C(5)=M%MachFile
    ArgT%C(6)=Ex
    ArgT%C(7)=N%SCF_NAME
    NewDex = 7
    DO K=1,SNC
      ArgT%C(NewDex+K)=S%Action%C(K)
    ENDDO
    NewDex=7+SNC
    ArgT%C(NewDex+1)=IntToChar(S%Current%I(1))
    ArgT%C(NewDex+2)=IntToChar(S%Current%I(2))
    ArgT%C(NewDex+3)=IntToChar(S%Current%I(3))
    ArgT%C(NewDex+4)=IntToChar(S%Previous%I(1))
    ArgT%C(NewDex+5)=IntToChar(S%Previous%I(2))
    ArgT%C(NewDex+6)=IntToChar(S%Previous%I(3))
    ArgT%C(NewDex+7)=TRIM(N%M_SCRATCH)
#else
    NArg=9+SNC
    CALL New(ArgT,NArg)
    ArgT%C(1)=Ex
    ArgT%C(2)=N%SCF_NAME
    NewDex = 2
    DO K=1,SNC
      ArgT%C(NewDex+K)=S%Action%C(K)
    ENDDO
    NewDex=2+SNC
    ArgT%C(NewDex+1)=IntToChar(S%Current%I(1))
    ArgT%C(NewDex+2)=IntToChar(S%Current%I(2))
    ArgT%C(NewDex+3)=IntToChar(S%Current%I(3))
    ArgT%C(NewDex+4)=IntToChar(S%Previous%I(1))
    ArgT%C(NewDex+5)=IntToChar(S%Previous%I(2))
    ArgT%C(NewDex+6)=IntToChar(S%Previous%I(3))
    ArgT%C(NewDex+7)=TRIM(N%M_SCRATCH)
#endif
    K=NArg
    NArg=0
    DO I=1,K
      IF(ArgT%C(I)/="")THEN
        NArg=NArg+1
      ENDIF
    ENDDO
#ifdef MPI2
    ! Add space for a trailing blank
    NArg=NArg+1
#endif
    CALL New(ArgV,NArg)
    NArg=0
    DO I=1,K
      IF(ArgT%C(I)/="")THEN
        NArg=NArg+1
        ArgV%C(NArg)=ArgT%C(I)
      ENDIF
    ENDDO
#ifdef MPI2
    ! Here is the trailing blank MPI_COMMM_SPAWN wants
    NArg=NArg+1
    ArgV%C(NArg)= " "
#endif
    CALL Delete(ArgT)
  END SUBROUTINE SetArgV

  !===============================================================
  ! CREATE AN INTEGER ARRAY OF ASCII KEYS FROM AN ARRAY OF
  ! CHARACTER STRINGS; USE FOR PORTABLE F9x/C INTERFACE
  !===============================================================
  SUBROUTINE CVToIV(NArg,ArgV,MaxLen,IChr)
    TYPE(CHR_VECT) :: ArgV
    TYPE(INT_VECT) :: IChr
    INTEGER        :: I,J,K,L,NArg,MaxLen

    ! Max number of characters in an element of ArgV
    MaxLen=0
    DO I=1,NArg
      MaxLen=MAX(MaxLen,LEN(TRIM(ArgV%C(I))))
    ENDDO

    ! Integer array to hold ASCII char-code
    CALL New(IChr,NArg*MaxLen)

    ! Convert strings to ASCII keys
    K=0
    DO I=1,NArg
      L=LEN(TRIM(ArgV%C(I)))
      DO J=1,MaxLen
        K=K+1
        IF(J<=L)THEN
          IChr%I(K)=ICHAR(ArgV%C(I)(J:J))
        ELSE
          IChr%I(K)=IBlnk
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE CVToIV
END MODULE
