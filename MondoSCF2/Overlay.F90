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
    CHARACTER(LEN=*)   :: Ex
    TYPE(FileNames)    :: N
    TYPE(State)        :: S
    TYPE(Parallel)     :: M
    INTEGER            :: I,J,K,L,NC,iCLUMP,IErr,NArg,MaxLen
    LOGICAL            :: ProgramFailed
    TYPE(CHR_VECT)     :: ArgV
    TYPE(INT_VECT),SAVE :: IChr ! Save needed for NAG
    CHARACTER(LEN=2*DCL) :: CmndLine
    INTERFACE 
       FUNCTION Spawn(NC,MaxLen,IChr)
         INTEGER                         :: NC,MaxLen
         INTEGER, DIMENSION(1:NC*MaxLen) :: IChr
         INTEGER                                     :: Spawn
       END FUNCTION Spawn
    END INTERFACE
    !------------------------------------------------------------!
    DO iCLUMP=1,M%Clumps
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
       ENDDO
       ! WRITE(*,*)' COMMANDLINE  = ',TRIM(CmndLine)
       ! Log this run
       CALL Logger(CmndLine,.FALSE.)
       ! Create ASCII integer array to beat F9x/C incompatibility
       CALL CVToIV(NArg,ArgV,MaxLen,IChr)
#ifdef PARALLEL
       !        CALL Wait(1D0)
#endif
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
#ifdef PARALLEL 
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
    CALL New(ArgV,NArg)
    NArg=0
    DO I=1,K
       IF(ArgT%C(I)/="")THEN
          NArg=NArg+1
          ArgV%C(NArg)=ArgT%C(I)
       ENDIF
    ENDDO    
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
