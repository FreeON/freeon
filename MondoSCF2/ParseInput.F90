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
#ifdef PARALLEL
  USE ParseParallel
#endif
CONTAINS 
  !===============================================================
  ! 'NUFF SAID
  !===============================================================
  SUBROUTINE ParseTheInput(C)
    TYPE(Controls) :: C
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
    CALL LoadDynamics(C%Nams,C%Opts,C%Geos,C%Dyns)
    ! Parse geometry or get from restart HDF 
    CALL LoadCoordinates(C%Nams,C%Opts,C%Dyns,C%Geos)
    ! Parse periodic info
    CALL LoadPeriodic(C%Nams,C%Geos,C%PBCs)
    ! Massage coodrinates, switch to AUs etc 
    CALL MassageCoordinates(C%Geos,C%PBCs)
    ! Load basis sets  
    CALL LoadBasisSets(C%Nams,C%Opts,C%Geos,C%Sets) 
#ifdef PARALLEL    
    ! Parse in parallel info
    CALL LoadParallel(C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
#endif
  END SUBROUTINE ParseTheInput
END MODULE ParseInput
