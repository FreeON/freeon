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
MODULE Overlay
  USE InOut
  USE PunchHDF
  USE ControlStructures
  IMPLICIT NONE
CONTAINS 
  !===============================================================
  !
  !===============================================================
  SUBROUTINE Invoke(Ex,N,S,M)
    CHARACTER(LEN=*)     :: Ex
    TYPE(FileNames)      :: N
    TYPE(State)          :: S
    TYPE(Parallel)       :: M
    INTEGER              :: I,J,K,L,NC,iCLUMP,IErr,NArg,MaxLen
    LOGICAL              :: ProgramFailed
    TYPE(CHR_VECT),SAVE  :: ArgV
    TYPE(INT_VECT),SAVE  :: IChr 
    CHARACTER(LEN=2*DCL) :: CmndLine
#if PARALLEL && MPI2
    INTEGER              :: SPAWN,ALL
#else
    INTERFACE 
       FUNCTION Spawn(NC,MaxLen,IChr)
         INTEGER                         :: NC,MaxLen
         INTEGER, DIMENSION(1:NC*MaxLen) :: IChr
         INTEGER                                     :: Spawn
       END FUNCTION Spawn
    END INTERFACE
#endif
    !------------------------------------------------------------!
    DO iCLUMP=1,M%Clumps
!       WRITE(*,*)'========================================================='
!       WRITE(*,*)' CLUMP = ',iCLUMP,' CLUMP = ',iCLUMP,' CLUMP = ',iCLUMP
!       WRITE(*,*)'========================================================='
#ifdef PARALLEL
       CALL MPIsArchive(N,M%NSpace,M%Clump%I(:,iCLUMP))
       CALL SetArgV(Ex,N,S,M,iCLUMP,NArg,ArgV)
#else 
       CALL MPIsArchive(N,1,M%Clump%I(:,iCLUMP))
       CALL SetArgV(Ex,N,S,M,NArg,ArgV)
#endif
       ! This is the command line we are going to execute 
       CmndLine=' '
       DO I=1,NArg
          CmndLine=TRIM(CmndLine)//Blnk//TRIM(ArgV%C(I))
          !          WRITE(*,*)I,' ARGVS <',TRIM(ArgV%C(I)),">"
       ENDDO
!       WRITE(*,*)iCLUMP,' COMMANDLINE = ',TRIM(CmndLine)
       ! Log this run
       CALL MondoLog(DEBUG_NONE, "Overlay:Invoke", TRIM(CmndLine))
#if MPI2
       CALL MPI_COMM_SPAWN(ArgV%C(1),ArgV%C(2:NArg),M%NProc,MPI_INFO_NULL, &
            ROOT,MPI_COMM_SELF,SPAWN,MPI_ERRCODES_IGNORE,IErr)
       IF(IErr/=MPI_SUCCESS)& 
          CALL MondoHalt(MPIS_ERROR,' Could not spawn <'//TRIM(CmndLine)//'>')
       ! Merge the spawned and current local communicators
       CALL MPI_INTERCOMM_MERGE(SPAWN,.TRUE.,ALL,IErr)
       ! Wait for the kiddies to be done
       CALL MPI_BARRIER(ALL,IErr)
#else
       ! Create ASCII integer array to beat F9x/C incompatibility
       CALL CVToIV(NArg,ArgV,MaxLen,IChr)
       ! Spawn a sub process 
       IErr=Spawn(NArg,MaxLen,IChr%I)
       ! Bring this run down if not successful
       IF(IErr/=SUCCEED)CALL MondoHalt(IErr,'<'//TRIM(CmndLine)//'>')
       ! Double check success if a MONDO Exec ...        
       HDF_CurrentID=OpenHDF(N%HFile)       
       CALL Get(ProgramFailed,'ProgramFailed')
       CALL CloseHDF(HDF_CurrentID)
       IF(ProgramFailed)CALL MondoHalt(-999,'<'//TRIM(CmndLine)//'>')
       CALL Delete(IChr)
#endif
       CALL Delete(ArgV)
    ENDDO
  END SUBROUTINE Invoke
  !===============================================================
  ! CREATE A CHARACTER ARRAY OF NON-BLANK STRINGS THAT WILL 
  ! BECOME THE ARGV ARRAY PASSED TO EXECVP BY SPAWN IF NOT MPI-2
  !===============================================================
#ifdef PARALLEL
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
    SNC=SIZE(S%Action%C)

#ifdef MPI2 
    NArg=8+SNC
    CALL New(ArgT,NArg)
    ArgT%C(1) =TRIM(N%M_EXEC)//'/'//Ex
    ArgT%C(2) =N%SCF_NAME
    DO K=1,SNC
      ArgT%C(2+K)=S%Action%C(K)
    ENDDO
    NewDex=2+SNC
    ArgT%C(NewDex+1)=IntToChar(S%Current%I(1))
    ArgT%C(NewDex+2)=IntToChar(S%Current%I(2))
    ArgT%C(NewDex+3)=IntToChar(S%Current%I(3))
    ArgT%C(NewDex+4)=IntToChar(S%Previous%I(1))
    ArgT%C(NewDex+5)=IntToChar(S%Previous%I(2))
    ArgT%C(NewDex+6)=IntToChar(S%Previous%I(3))
#elif PARALLEL
    NArg=13+SNC
    CALL New(ArgT,NArg)
    ArgT%C(1) =M%Invoking
    ArgT%C(2) =M%ProcFlag
    ArgT%C(3) =IntToChar(M%Clump%I(3,cCLUMP))
    ArgT%C(4) =M%MachFlag
    ArgT%C(5) =M%MachFile
    ArgT%C(6) =TRIM(N%M_EXEC)//'/'//Ex
    ArgT%C(7) =N%SCF_NAME
    DO K=1,SNC
       ArgT%C(7+K)=S%Action%C(K)
    ENDDO
    NewDex=7+SNC
    ArgT%C(NewDex+1)=IntToChar(S%Current%I(1))
    ArgT%C(NewDex+2)=IntToChar(S%Current%I(2))
    ArgT%C(NewDex+3)=IntToChar(S%Current%I(3))
    ArgT%C(NewDex+4)=IntToChar(S%Previous%I(1))
    ArgT%C(NewDex+5)=IntToChar(S%Previous%I(2))
    ArgT%C(NewDex+6)=IntToChar(S%Previous%I(3))
#else
    NArg=8+SNC
    CALL New(ArgT,NArg)
    ArgT%C(1) =TRIM(N%M_EXEC)//'/'//Ex
    ArgT%C(2) =N%SCF_NAME
    DO K=1,SNC
      ArgT%C(2+K)=S%Action%C(K)
    ENDDO
    NewDex=2+SNC
    ArgT%C(NewDex+1)=IntToChar(S%Current%I(1))
    ArgT%C(NewDex+2)=IntToChar(S%Current%I(2))
    ArgT%C(NewDex+3)=IntToChar(S%Current%I(3))
    ArgT%C(NewDex+4)=IntToChar(S%Previous%I(1))
    ArgT%C(NewDex+5)=IntToChar(S%Previous%I(2))
    ArgT%C(NewDex+6)=IntToChar(S%Previous%I(3))
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
    !------------------------------------------------------------!
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
