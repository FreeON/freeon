MODULE ParseInput
  USE ControlStructures
  USE ParseCommands
  USE ParseOptions
  USE ParseDynamics
  USE ParseGeometries
  USE ParseBasis
CONTAINS 
  !===============================================================
  ! 'NUFF SAID
  !===============================================================
  SUBROUTINE ParseTheInput(C)
    TYPE(Controls) :: C
    !-------------------------------------------------------------!
    ! Parse command line and load env and file names 
    CALL LoadCommands(C % Nams)
    ! Parse generic options 
    CALL LoadOptions(C % Nams, C % Opts)
    ! Parse dynamics options
    CALL LoadDynamics(C % Nams, C % Opts, C % Geos, C % Dyns)
    ! Parse geometry or get from restart HDF 
    CALL LoadGeometry(C % Nams, C % Opts, C % Geos)
    ! Parse periodic info
    CALL LoadPeriodic(C % Nams, C % Geos, C % PBCs)
    ! Massage coodrinates after periodic info is loaded
    CALL MassageCoordinates( C % Geos, C % PBCs)
    ! Load basis sets.  
    CALL LoadBasis( C % Nams, C % Geos, C % Sets ) 
  END SUBROUTINE ParseTheInput
END MODULE ParseInput
