MODULE ParseCommands
  USE Parse
  USE InOut
#ifdef NAG
  USE F90_UNIX
#endif
  USE ControlStructures
  IMPLICIT NONE
CONTAINS
!===============================================================================================
! PARSE THE COMMAND LINE AND GET RELATED ENV VARIABLES. CHECK EXISTENCE OF IN AND OUT FILES
!===============================================================================================
  SUBROUTINE LoadCommands(N)
    TYPE(FileNames)     :: N
    TYPE(ARGMT)         :: Args
    CHARACTER(LEN=DCL)  :: PROCESS_ID,PWDName,SCRName
    INTEGER             :: DotDex
    LOGICAL             :: Exists
    INTEGER,EXTERNAL    :: GetPID
!-----------------------------------------------------------------------------------------------         
    ! Get command line arguments 
    CALL Get(Args)
    IF(Args%NC==0)CALL MondoHalt(PRSE_ERROR,' No arguments to MondoSCF !')!
    ! Get environmental variables
    CALL GetEnv('PWD',N%M_PWD)
    CALL GetEnv('MONDO_HOME',N%M_HOME)
    IF(LEN(TRIM(N%M_HOME))==0) &
       CALL MondoHalt(PRSE_ERROR,'env variable $(MONDO_HOME) not set.')
    CALL GetEnv('MONDO_SCRATCH',N%M_SCRATCH)
    IF(LEN(TRIM(N%M_SCRATCH))==0) &
       CALL MondoHalt(PRSE_ERROR,'env variable $(MONDO_SCRATCH) not set.')
    ! Set path names etc
    N%M_PWD=TRIM(N%M_PWD)//'/'
    N%M_HOME=TRIM(N%M_HOME)//'/'
    N%M_SCRATCH=TRIM(N%M_SCRATCH)//'/'
    ! Here is the scf name, tagged with the PID
    PROCESS_ID=IntToChar(GetPID())
    DotDex=INDEX(Args%C%C(1),'.')
    IF(DotDex==0) &
         CALL MondoHalt(PRSE_ERROR,' No "." in inut file name = <'//TRIM(N%IFile)//'>')
    N%SCF_NAME=Args%C%C(1)(1:DotDex-1)//'_'//TRIM(PROCESS_ID)
    PWDName=TRIM(N%M_PWD)//TRIM(N%SCF_NAME)
    SCRName=TRIM(N%M_SCRATCH)//TRIM(N%SCF_NAME)
    ! Input file with full path
    N%IFile=TRIM(N%M_PWD)//TRIM(Args%C%C(1))
    ! Check to see that the input file exists
    INQUIRE(FILE=N%IFile,EXIST=Exists)
    IF(.NOT.Exists) &
         CALL MondoHalt(PRSE_ERROR,' Input file: '//TRIM(N%IFile)//' does not exist! ')
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
       N%GFile=TRIM(PWDName)//GeoF
    ELSEIF(Args%NC==2)THEN
       N%OFile=TRIM(N%M_PWD)//TRIM(Args%C%C(2))         
       N%LFile=TRIM(PWDName)//LogF
       N%GFile=TRIM(PWDName)//GeoF
    ELSEIF(Args%NC==3)THEN
       N%OFile=TRIM(N%M_PWD)//TRIM(Args%C%C(2))         
       N%LFile=TRIM(N%M_PWD)//TRIM(Args%C%C(3))
       N%GFile=TRIM(PWDName)//GeoF
    ELSEIF(Args%NC==4)THEN
       N%OFile=TRIM(N%M_PWD)//TRIM(Args%C%C(2))         
       N%LFile=TRIM(N%M_PWD)//TRIM(Args%C%C(3))
       N%GFile=TRIM(N%M_PWD)//TRIM(Args%C%C(4))
    ENDIF
    N%HFile=TRIM(SCRName)//InfF
    ! If restart file set, it is set in ParseOptions
    N%RFile=""
    ! 
#if FULL_ON_FRONT_END_DEBUG    
    WRITE(*,*)' N%IFile = ',TRIM(N%IFile)
    WRITE(*,*)' N%OFile = ',TRIM(N%OFile)
    WRITE(*,*)' N%LFile = ',TRIM(N%LFile)
    WRITE(*,*)' N%HFile = ',TRIM(N%HFile)
    WRITE(*,*)' N%RFile = ',TRIM(N%RFile)
    WRITE(*,*)' N%GFile = ',TRIM(N%GFile)
#endif
  END SUBROUTINE LoadCommands
END MODULE ParseCommands
