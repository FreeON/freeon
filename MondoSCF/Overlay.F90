!  Authors: Matt Challacombe and Chee Kwan Gan
MODULE Overlay
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalCharacters
   USE ProcessControl
   USE Clock
#ifdef NAG
   USE F90_UNIX
#endif
   USE SCFLocals
   USE ParsingConstants
   USE Parse
   USE MemMan
   USE InOut
#ifdef PARALLEL
   USE MondoMPI
#endif
   IMPLICIT NONE
   CONTAINS 
!
      SUBROUTINE Invoke(Exec,CArgs_O,IArgs_O,AbsPath_O,MPIRun_O)
         CHARACTER(LEN=*),                      INTENT(IN) :: Exec
         CHARACTER(LEN=*),DIMENSION(:),OPTIONAL,INTENT(IN) :: CArgs_O
         INTEGER,         DIMENSION(:),OPTIONAL,INTENT(IN) :: IArgs_O
         LOGICAL,                      OPTIONAL            :: AbsPath_O,MPIRun_O
         INTEGER                                           :: I
         CHARACTER(LEN=4*DEFAULT_CHR_LEN)                  :: Command,CmndOpts

!         IF(PRESENT(MPIRun_O))THEN
!            CALL MWait(60.0D0) ! to live with myrinet ...
!         ENDIF
!-------------------------------------------------------------------------------
!        Find the character arguments if they exist
         CmndOpts=' '
         IF(PRESENT(CArgs_O))THEN
            DO I=1,SIZE(CArgs_O)
               CmndOpts=TRIM(CmndOpts)//Blnk//TRIM(CArgs_O(I))
            ENDDO
         ENDIF
!        Find the integer arguments if they exist
         IF(PRESENT(IArgs_O))THEN
            DO I=1,SIZE(IArgs_O)
               CmndOpts=TRIM(CmndOpts)//Blnk//TRIM(IntToChar(IArgs_O(I)))
            ENDDO
         ENDIF
!        Put the binary name together with the full path
         IF(PRESENT(AbsPath_O))THEN
            IF(AbsPath_O)THEN
               Command=Exec
            ELSE
               CALL MondoHalt(-99,'Bad logic in Invoke')
            ENDIF
         ELSE
            Command=TRIM(MONDO_Exec)//TRIM(Exec)
         ENDIF

#ifdef PARALLEL
         IF(PRESENT(MPIRun_O).OR.NPrc==1)THEN
#ifdef MPI2 
            CALL MPISpawn(Command,CmndOpts)
#else    
            CALL SerialSpawn(Command,CmndOpts,MPIRun_O=MPIRun_O,AbsPath_O=AbsPath_O)
#endif
         ELSE
#endif
            CALL SerialSpawn(Command,CmndOpts,AbsPath_O=AbsPath_O)
#ifdef PARALLEL 
         ENDIF
#endif
      END SUBROUTINE Invoke
!
      SUBROUTINE SerialSpawn(Command,CmndOpts,MPIRun_O,AbsPath_O)
         CHARACTER(LEN=*),INTENT(IN)      :: Command,CmndOpts
         TYPE(CHR_VECT)                   :: ArgV
         TYPE(INT_VECT)                   :: IChr
         INTEGER                          :: I,J,K,L,N,NArg,MaxLen,IErr
         CHARACTER(LEN=DEFAULT_CHR_LEN)   :: ENV_VAR
         CHARACTER(LEN=4*DEFAULT_CHR_LEN) :: CmndLine
         LOGICAL,     OPTIONAL            :: MPIRun_O,AbsPath_O
         LOGICAL                          :: ProgramFailed
         INTERFACE 
            FUNCTION Spawn(NC,MaxLen,IChr)
               INTEGER,                         INTENT(IN) :: NC,MaxLen
               INTEGER, DIMENSION(1:NC*MaxLen), INTENT(IN) :: IChr
               INTEGER                                     :: Spawn
            END FUNCTION Spawn
         END INTERFACE
!--------------------------------------------------------------------         
#if defined(PARALLEL) && !defined(MPI2) 
         IF(PRESENT(MPIRun_O))THEN
            IF(MPIRun_O)THEN
               IF(INDEX(MPI_INVOKE,'poe')/=0)THEN
                  CmndLine=TRIM(Command)//Blnk//TRIM(CmndOpts)
               ELSE
                  CmndLine=TRIM(MPI_INVOKE)//Blnk   &
                         //TRIM(MPI_FLAGS) //Blnk   & 
                         //TRIM(Command)   //Blnk   &
                         //TRIM(CmndOpts)              
               ENDIF
            ELSE
               CALL MondoHalt(-100,'Logic error 1 in Overlay.')
            ENDIF
         ELSE
#endif
            CmndLine=TRIM(Command)//Blnk//TRIM(CmndOpts)              
#if defined(PARALLEL) && !defined(MPI2)
         ENDIF
#endif
         CALL LineToChars(CmndLine,ArgV)
         NArg=SIZE(ArgV%C)
!        Exapnd environmental variables if any
         DO I=1,NArg
            IF(SCAN(ArgV%C(I),'$')/=0)THEN
               ArgV%C(I)=ADJUSTL(ArgV%C(I))  
               L=LEN(ArgV%C(I))
               CALL GETENV(TRIM(ArgV%C(I)(2:L)),ENV_VAR)
               ArgV%C(I)=ENV_VAR
            ENDIF
         ENDDO
         CmndLine=' '
         DO I=1,NArg
            CmndLine=TRIM(CmndLine)//Blnk//TRIM(ArgV%C(I))
         ENDDO
!        Max number of characters in an arg
         MaxLen=0
         DO I=1,NArg   
            MaxLen=MAX(MaxLen,LEN(TRIM(ArgV%C(I))))
         ENDDO
         CALL New(IChr,NArg*MaxLen)
!        Convert strings to ASCII keys
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
         IF(.NOT.PRESENT(AbsPath_O))THEN
            CALL OpenHDF(InfFile)
            CALL Put(.TRUE.,'ProgramFailed')
            CALL CloseHDF()
         ENDIF
!        Log this run
         CALL Logger(CmndLine,.FALSE.)
!        Fork a serial child
         IErr=Spawn(NArg,MaxLen,IChr%I)
!        Bring this run down if not successful
         IF(IErr/=SUCCEED)CALL MondoHalt(IErr,'<'//TRIM(CmndLine)//'>')
!        Double check success if a MONDO Exec ...        
         IF(.NOT.PRESENT(AbsPath_O))THEN
            CALL OpenHDF(InfFile)
            CALL Get(ProgramFailed,'ProgramFailed')
            CALL CloseHDF()
            IF(ProgramFailed)CALL MondoHalt(-999,'<'//TRIM(CmndLine)//'>')
         ENDIF
!        Tidy up ...
         CALL Delete(IChr)
         CALL Delete(ArgV)
      END SUBROUTINE SerialSpawn

#if defined(PARALLEL) && defined(MPI2) 
!--------------------------------------------------------------------------
      SUBROUTINE MPISpawn(Command,CmndOpts)
         CHARACTER(LEN=*),INTENT(IN)         :: Command,CmndOpts
         TYPE(CHR_VECT)                      :: ArgV
         TYPE(INT_VECT)                      :: ErrCodes
         INTEGER                             :: I,J,K,L,N,NArg,LenMsg,Status
         CHARACTER(LEN=2*DEFAULT_CHR_LEN)    :: CmndLine
         CHARACTER(LEN=MPI_MAX_ERROR_STRING) :: ErrMsg
         CHARACTER(LEN=2*DEFAULT_CHR_LEN)    :: LogMssg
         INTEGER                             :: NProcess,IErr,Child_Inter
!------------------------------------------------------------------------         
         Status=SUCCEED
!        CmndLine=TRIM(Command)//CmndOpts
         CmndLine=CmndOpts
!        CALL InitMPI()
         NProcess = MSize()
         IF(NProcess /= NPrc) THEN
           write(*,*) 'ERROR in MPISpawn(): NProcess not equal to NPrc'
         ENDIF
         IF(MyId==ROOT)THEN
            CALL LineToChars(CmndLine,ArgV,NULL_O=.TRUE.)
            NArg=SIZE(ArgV%C)
            CALL New(ErrCodes,NPrc)            
            CALL MPI_COMM_SPAWN(Command,ArgV%C,NPrc,MPI_INFO_NULL, &
                                ROOT,MPI_COMM_SELF,child_inter,ErrCodes%I,IErr)
            DO I=1,NPrc
               IF(ErrCodes%I(I)/=MPI_SUCCESS)THEN
                  CALL MPI_ERROR_STRING(ErrCodes%I(I),ErrMsg,LenMsg)
                  LogMssg='MPI ERROR = <'//TRIM(ErrMsg(1:LenMsg))//'>'
                  CALL Logger(LogMssg,.FALSE.)
                  Status=FAIL
               ENDIF
             ENDDO
!            Tidy up ...
             CALL Delete(ArgV)
             CALL Delete(ErrCodes)
         ENDIF
         IF(Status==FAIL)THEN
!           Bad run, call a halt :(
            CALL MondoHalt(MPIS_ERROR,'<'//TRIM(CmndLine)//'>')
         ELSE
!           Log succesful run :)
            CALL Logger(CmndLine,.FALSE.)
         ENDIF
      END SUBROUTINE MPISpawn
#endif
END MODULE
