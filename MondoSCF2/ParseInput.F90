MODULE ParseInput
  USE Massage
  USE ParseBasis
  USE ParseOptions
  USE ParseCommands
  USE ParseDynamics
  USE ParsePeriodic
  USE ParseGeometries
  USE ControlStructures
  USE PrettyPrint
  USE ParseParallel
  USE ParseGeomOpt
  USE ParseExtraCoords
  USE ParseProperties, ONLY: LoadPropertyOptions
CONTAINS 
  !===============================================================
  ! 'NUFF SAID
  !===============================================================
  SUBROUTINE ParseTheInput(C)
    TYPE(Controls) :: C
    integer        :: iclone
    !-------------------------------------------------------------!
    ! Parse command line and load env and file names 
    CALL LoadCommands(C%Nams)
    ! Set global output and log file
    OutFile=C%Nams%OFile
    LogFile=C%Nams%LFile
    ! Allocate state variables
    CALL New(C%Stat%Current,3)
    CALL New(C%Stat%Previous,3)
    ! Parse generic options 
    CALL LoadOptions(C%Nams,C%Opts)
    ! Parse dynamics options
    IF(C%Opts%Grad==GRAD_DO_DYNAMICS) THEN
       CALL LoadDynamics(C%Nams,C%Opts,C%Geos,C%Dyns)
    ENDIF  
    ! Parse geometry or get from restart HDF 
    CALL LoadCoordinates(C%Nams,C%Opts,C%Dyns,C%Geos)
    ! Parse periodic info
    CALL LoadPeriodic(C%Nams,C%Opts,C%Geos,C%PBCs)
    ! Massage coodrinates, switch to AUs etc 
    CALL MassageCoordinates(C%Opts,C%Geos,C%PBCs)
    ! Load basis sets  
    CALL LoadBasisSets(C%Nams,C%Opts,C%Geos,C%Sets) 
    ! Parse in parallel info
    CALL LoadParallel(C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
    ! Load control of internal coord. optimizer
    CALL LoadGeomOpt(C%Nams,C%GOpt)
    ! Load constraints and extra internal coords
    CALL LoadExtraCoords(C%GOpt,C%Opts,C%Nams,C%Geos)
    ! Load CPSCF options.
    CALL LoadPropertyOptions(C%Nams,C%POpt)
    ! Check for Global conflicts.
    CALL ConflictCheck(C)
 END SUBROUTINE ParseTheInput
END MODULE ParseInput
