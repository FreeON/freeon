MODULE ParseInput
  USE ControlStructures
  USE ParseCommands
  USE ParseOptions
  USE ParseDynamics
  USE ParseGeometries
  !USE ConflictDetect
  !USE ParseBasis
  !USE PrintParsed
CONTAINS 

  SUBROUTINE LoadInput(C)
    TYPE(Controls) :: C
    !-------------------------------------------------------------!
    ! Parse command line, load names and OPEN THE HDF.
    CALL LoadCommands(C % Nams)
    ! Parse generic options 
    CALL LoadOptions(C % Nams, C % Opts)
    ! Parse dynamics options
    CALL LoadDynamics(C % Nams, C % Opts, C % Geos, C % Dyns)
    ! Check for some logical errors
!    CALL ConflictCheck1(C % Opts, C % Dyns)
    ! Parse geometry or get from restart HDF 
    CALL LoadGeometry(C % Nams, C % Opts, C % Geos)
    ! Parse periodic info
    CALL LoadPeriodic(C % Nams, C % Geos, C % PBCs)
    ! Massage coodrinates after periodic info is loaded
    CALL MassageCoordinates( C % Geos, C % PBCs)
    ! Load basis sets.  
    CALL LoadBasis( C % Nams, C % Geos, C % Sets ) 
    ! >>>> NOTE BASIS SET ORDERING DEPENDS ON ORDERING OF 
    !  GEOMETRY.  NO ORDERING OF THE GEOMETRY CAN OCCUR PAST 
    !  THIS POINT UNLESS THE BASIS SET IS ALSO REORDERED <<<<

    ! Check for more logical errors
!    CALL ConflictCheck2(C % Opts, C % Dyns, C % Sets)

!    CALL PrintParsed()
  END SUBROUTINE LoadInput

END MODULE ParseInput
