MODULE GeometryKeys
  ! Coordinate input delemeters
  CHARACTER(LEN=15), PARAMETER :: GEOMETRY_BEGIN ='<BeginGeometry>'
  CHARACTER(LEN=14), PARAMETER :: GEOMETRY_NEXT  ='<NextGeometry>'
  CHARACTER(LEN=13), PARAMETER :: GEOMETRY_END   ='<EndGeometry>'  

  CHARACTER(LEN=15), PARAMETER :: PRODUCTS_BEGIN ='<BeginProducts>'
  CHARACTER(LEN=13), PARAMETER :: PRODUCTS_END   ='<EndProducts>'

  CHARACTER(LEN=16), PARAMETER :: REACTANTS_BEGIN='<BeginReactants>'
  CHARACTER(LEN=14), PARAMETER :: REACTANTS_END  ='<EndReactants>'
  ! Geometry=
  CHARACTER(LEN=8),  PARAMETER :: GEOMETRY          ='Geometry' 
  CHARACTER(LEN=8),  PARAMETER :: IN_AU             ='InAU'
  CHARACTER(LEN=6),  PARAMETER :: H_ORDER           ='HOrder'
  CHARACTER(LEN=6),  PARAMETER :: Z_ORDER           ='ZOrder'
  CHARACTER(LEN=11), PARAMETER :: RANDOM_ORDER      ='RandomOrder'
  CHARACTER(LEN=11), PARAMETER :: TRAVEL_ORDER      ='TravelOrder'
  CHARACTER(LEN=14), PARAMETER :: TABLETRAV_ORDER   ='TableTravOrder'
  CHARACTER(LEN=7),  PARAMETER :: NO_ORDER          ='NoOrder'
  ! Clones=
  CHARACTER(LEN=6),  PARAMETER :: CLONES            ='Clones'
  !------------------------------------------------------------------------------
  ! Charge=
  CHARACTER(LEN=6),  PARAMETER :: TOTAL_CHARGE      ='Charge' 
  !------------------------------------------------------------------------------
  ! Multiplicity=
  CHARACTER(LEN=12), PARAMETER :: MULTIPLICITY      ='Multiplicity' 
END MODULE GeometryKeys
