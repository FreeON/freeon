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
    CHARACTER(LEN=DCL)  :: PROCESS_ID,PWDName
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
       CALL MondoHalt(PRSE_ERROR,' $(N%M_HOME) not set.')
    CALL GetEnv('MONDO_SCRATCH',N%M_SCRATCH)
    IF(LEN(TRIM(N%M_SCRATCH))==0) &
       CALL MondoHalt(PRSE_ERROR,' $(N%M_SCRATCH) not set.')
    ! Set path names etc
    N%M_PWD=TRIM(N%M_PWD)//'/'
    N%M_HOME=TRIM(N%M_HOME)//'/'
    N%M_SCRATCH=TRIM(N%M_SCRATCH)//'/'
    PROCESS_ID=IntToChar(GetPID())
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
    IF(Args%NC==1) &
         CALL MondoHalt(PRSE_ERROR,' Please specify an output file' )
    ! Out file with full path (demand an output file)
    N%OFile=TRIM(N%M_PWD)//TRIM(Args%C%C(2))         
    ! Check to see that the output file does not exist (no overwrites allowed)
    INQUIRE(FILE=N%IFile,EXIST=Exists)
    IF(Exists) &
       CALL MondoHalt(PRSE_ERROR,' Ouput file: '//TRIM(N%IFile)//' already exists! ')
    ! Create user defined or implicit file names
    IF(Args%NC==2)THEN
       N%LFile=TRIM(PWDName)//LogF
       N%GFile=TRIM(PWDName)//GeoF
    ELSEIF(Args%NC==3)THEN
       N%LFile=TRIM(N%M_PWD)//TRIM(Args%C%C(3))
       N%GFile=TRIM(PWDName)//GeoF
    ELSEIF(Args%NC==4)THEN
       N%LFile=TRIM(N%M_PWD)//TRIM(Args%C%C(3))
       N%GFile=TRIM(N%M_PWD)//TRIM(Args%C%C(4))
    ENDIF
    ! HDF archival and restart file 
    N%HFile=TRIM(N%M_SCRATCH)//TRIM(N%SCF_NAME)//InfF                
    ! clean up memory
    CALL Delete(Args)
    ! Initialize and open the new HDF file
    CALL InitHDF(N%HFile)
    N%NewFileID=OpenHDFFile(N%HFile)
    ! Set globally accesable HDFFileID 
    HDFFileID=N%NewFileID
    ! Put current working file names to HDF
    CALL Put(N%IFile,'inputfile')
    CALL Put(N%HFile,'infofile')
    CALL Put(N%LFile,'logfile')
    CALL Put(N%OFile,'outPutfile')
  END SUBROUTINE LoadCommands
END MODULE ParseCommands
