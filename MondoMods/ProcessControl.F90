MODULE ProcessControl
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalCharacters
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
   INTEGER, PARAMETER     :: MPIS_ERROR=-975239
   INTEGER, PARAMETER     :: USUP_ERROR=-993942
   INTEGER, PARAMETER     :: INTC_ERROR=-135950
   INTEGER, PARAMETER     :: QMMM_ERROR=-458391
   CONTAINS 
      SUBROUTINE MondoHalt(IErr,Mssg)
         CHARACTER (LEN=*) :: Mssg
         INTEGER           :: IErr
         IF(IErr==EXIT_ERROR)THEN
            WRITE(*,*)  TRIM(Mssg)//' failed.'
            CALL Logger(TRIM(Mssg)//' failed.',.TRUE.)
         ELSEIF(IErr==DUMP_ERROR)THEN
            WRITE(*,*)  TRIM(Mssg)//' core dumped.'
            CALL Logger(TRIM(Mssg)//' core dumped.',.TRUE.)
         ELSEIF(IErr==SGNL_ERROR)THEN
            WRITE(*,*)  TRIM(Mssg)//' signaled.'
            CALL Logger(TRIM(Mssg)//' signaled.',.TRUE.)
         ELSEIF(IErr==FORK_ERROR)THEN
            WRITE(*,*)  'Unable to fork '//TRIM(Mssg)
            CALL Logger('Unable to fork '//TRIM(Mssg),.TRUE.)
         ELSEIF(IErr==DRIV_ERROR)THEN
            WRITE(*,*)  'Error in SCFDriver: '//TRIM(Mssg)
            CALL Logger('Error in SCFDriver: '//TRIM(Mssg),.TRUE.)
         ELSEIF(IErr==PRSE_ERROR)THEN
            WRITE(*,*)  'Error parsing: '//TRIM(Mssg)
            CALL Logger('Error parsing: '//TRIM(Mssg),.TRUE.)
         ELSEIF(IErr==MPIS_ERROR)THEN
            WRITE(*,*)  'MPI error: '//TRIM(Mssg)
            CALL Logger('MPI error: '//TRIM(Mssg),.TRUE.)
         ELSEIF(IErr==USUP_ERROR)THEN
            WRITE(*,*)  'Unsupported feature: '//TRIM(Mssg)
            CALL Logger('Unsupported feature: '//TRIM(Mssg),.TRUE.)
         ELSEIF(IErr==INTC_ERROR)THEN
            WRITE(*,*)  'Error in IntCoo: '//TRIM(Mssg)
            CALL Logger('Error in IntCoo: '//TRIM(Mssg),.TRUE.)
         ELSEIF(IErr==QMMM_ERROR)THEN
            WRITE(*,*)  'Error in QMMM  : '//TRIM(Mssg)
            CALL Logger('Error in QMMM  : '//TRIM(Mssg),.TRUE.)
         ELSE
            WRITE(*,*)  'Unknown error: '//TRIM(Mssg)
            CALL Logger('Unknown error: '//TRIM(Mssg),.TRUE.)
         ENDIF
         WRITE(*,*)' '
         WRITE(*,*)' To have no errors '
         WRITE(*,*)' Would be life without meaning '
         WRITE(*,*)' No struggle, no joy '
         WRITE(*,*)' '
         STOP 'Termination of MondoSCF'
      END SUBROUTINE MondoHalt

      SUBROUTINE Halt(Strng)
         CHARACTER (LEN=*) :: Strng
         CHARACTER(LEN=*),PARAMETER :: Motto='Frango ut patefaciam -- I break in order to reveal'
         CALL Logger(Strng,.TRUE.)
#ifdef PARALLEL
         IF(InParallel)THEN
            WRITE(*,*)' Node that brought me down = ',MyId
            CALL Trap()
            CALL HaltMPI(Strng)
         ELSE
            WRITE(*,*)' GOING DOWN GOING DOWN GOING DOWN GOING DOWN GOING DOWN '
            WRITE(*,*)'<<'//TRIM(Strng)//'>>'
            WRITE(*,*)' GOING DOWN GOING DOWN GOING DOWN GOING DOWN GOING DOWN '
            CALL Trap()
         ENDIF
#else
         WRITE(*,*)' GOING DOWN GOING DOWN GOING DOWN GOING DOWN GOING DOWN '
         WRITE(*,*)'<<'//TRIM(Strng)//'>>'
         WRITE(*,*)' GOING DOWN GOING DOWN GOING DOWN GOING DOWN GOING DOWN '
         CALL Trap()
#endif
      END SUBROUTINE Halt

#ifdef PARALLEL
      SUBROUTINE HaltMPI(InMsg,IErr)
         USE MPIInclude
         CHARACTER(LEN=*), INTENT(IN)        :: InMsg
         INTEGER, OPTIONAL,INTENT(IN)        :: IErr
         CHARACTER(LEN=MPI_MAX_ERROR_STRING) :: ErrMsg
         CHARACTER(LEN=2*DEFAULT_CHR_LEN)    :: LogMssg
         CHARACTER(LEN=INTERNAL_INT_LEN)     :: ChId
         INTEGER                             :: LenMsg
         WRITE(ChId,INTERNAL_INT_FMT)MyId
         ChId=ADJUSTL(ChId)
         LogMssg='GOING, GOING,... GONE DOWN: MyID = '//TRIM(ChId) &
                 //Rtrn//'Message that killed me: '//TRIM(InMsg)
         IF(PRESENT(IErr))THEN
            CALL MPI_ERROR_STRING(IErr,ErrMsg,LenMsg)
            LogMssg=ADJUSTL(TRIM(LogMssg))//Rtrn &
                  //'MPI ERROR = <'//TRIM(ErrMsg(1:LenMsg))//'>'
!            CALL LogInTurn(LogMssg)
            CALL Logger(LogMssg,.TRUE.)
            CALL MPI_ABORT(MPI_COMM_WORLD,IErr)
         ELSE
!            CALL LogInTurn(LogMssg)
            CALL Logger(LogMssg,.TRUE.)
            CALL MPI_ABORT(MPI_COMM_WORLD,0)
         ENDIF
      END SUBROUTINE HaltMPI

      SUBROUTINE LogInTurn(Mssg,IErr)
         USE MPIInclude
         CHARACTER(LEN=*), INTENT(IN) :: Mssg
         INTEGER, OPTIONAL,INTENT(IN) :: IErr
         INTEGER                      :: I
         DO I=0,NPrc-1
            CALL MPI_BARRIER(MPI_COMM_WORLD,IErr)
            IF(MyId==I)CALL Logger(Mssg,.TRUE.)
            CALL MPI_BARRIER(MPI_COMM_WORLD,IErr)
         ENDDO  
      END SUBROUTINE LogInTurn
#endif

      SUBROUTINE Warn(Strng)
         CHARACTER (LEN=*) Strng
         CHARACTER (LEN=DEFAULT_CHR_LEN) :: Warning
         Warning='... ning Warning  Warning  Warning Wa ...'//Rtrn//TRIM(ADJUSTL(Strng))
         CALL Logger(Warning,.TRUE.)
      END SUBROUTINE Warn

      SUBROUTINE Logger(Mssg,ToOutFile_O)
         CHARACTER(LEN=*)    :: Mssg
         LOGICAL,OPTIONAL    :: ToOutFile_O
         INTEGER             :: IOS
         LOGICAL             :: Opened, Exists, ToOutFile
         ToOutFile=.FALSE.
         IF(PRESENT(ToOutFile_O))ToOutFile=ToOutFile_O          
!------------------------------------------------------------------
!        Whats up with the log file?
!
         INQUIRE(FILE=TRIM(LogFile),OPENED=Opened, &
                 EXIST=Exists,ERR=11,IOSTAT=IOS)
!------------------------------------------------------------------
!        Open existing file and position at the bottom 
!
         IF(Exists.AND.(.NOT.Opened))THEN
            OPEN(UNIT=LgF,FILE=TRIM(LogFile), &
                 ACCESS='SEQUENTIAL', FORM='FORMATTED', &
                 POSITION='APPEND',ERR=11,IOSTAT=IOS,STATUS='OLD')
!------------------------------------------------------------------
!        Bad logic....
!
         ELSEIF(Exists.AND.Opened)THEN
           WRITE(*,*)' File '//TRIM(LogFile)//' already open'
!------------------------------------------------------------------
!        Create a new file and open it
!
         ELSE
           OPEN(UNIT=LgF,FILE=TRIM(LogFile), &
                ACCESS='SEQUENTIAL',FORM='FORMATTED', &
                ERR=11,IOSTAT=IOS,STATUS='NEW')         
         ENDIF
         WRITE(LgF,*)TRIM(Mssg)
         CLOSE(UNIT=LgF,STATUS='KEEP')

         IF(.NOT.ToOutFile)RETURN
!------------------------------------------------------------------
!        Whats up with the Out file?
!
         INQUIRE(FILE=TRIM(OutFile),OPENED=Opened, &
                 EXIST=Exists,ERR=12,IOSTAT=IOS)
!------------------------------------------------------------------
!        Open existing file and position at the bottom 
!
         IF(Exists.AND.(.NOT.Opened))THEN
            OPEN(UNIT=Out,FILE=TRIM(OutFile), &
                 ACCESS='SEQUENTIAL', FORM='FORMATTED', &
                 POSITION='APPEND',ERR=12,IOSTAT=IOS,STATUS='OLD')
!------------------------------------------------------------------
!        Bad Outic....
!
         ELSEIF(Exists.AND.Opened)THEN
           WRITE(*,*)' File '//TRIM(OutFile)//' already open'
!------------------------------------------------------------------
!        Create a new file and open it
!
         ELSE
           OPEN(UNIT=Out,FILE=TRIM(OutFile), &
                ACCESS='SEQUENTIAL',FORM='FORMATTED', &
                ERR=12,IOSTAT=IOS,STATUS='NEW')         
         ENDIF
         WRITE(Out,*)TRIM(Mssg)
         CLOSE(UNIT=Out,STATUS='KEEP')
         RETURN
      11 WRITE(*,*)TRIM(Mssg)
         WRITE(*,*)' Logging Error: IOS= ',IOS,' on file ',TRIM(LogFile)
         CALL Trap()
      12 WRITE(*,*)TRIM(Mssg)
         WRITE(*,*)' Logging Error: IOS= ',IOS,' on file ',TRIM(OutFile)
         CALL Trap()
      END SUBROUTINE Logger

END MODULE







