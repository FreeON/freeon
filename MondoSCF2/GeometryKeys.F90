MODULE GeometryKeys
  ! Coordinate input delemeters
  CHARACTER(LEN=*), PARAMETER :: GEOMETRY_BEGIN ='<BeginGeometry>'
  CHARACTER(LEN=*), PARAMETER :: GEOMETRY_NEXT  ='<NextGeometry>'
  CHARACTER(LEN=*), PARAMETER :: GEOMETRY_END   ='<EndGeometry>'  

  CHARACTER(LEN=*), PARAMETER :: PRODUCTS_BEGIN ='<BeginProducts>'
  CHARACTER(LEN=*), PARAMETER :: PRODUCTS_END   ='<EndProducts>'

  CHARACTER(LEN=*), PARAMETER :: REACTANTS_BEGIN='<BeginReactants>'
  CHARACTER(LEN=*), PARAMETER :: REACTANTS_END  ='<EndReactants>'
  ! Geometry=
  CHARACTER(LEN=*),  PARAMETER :: GEOMETRY          ='Geometry' 
  CHARACTER(LEN=*),  PARAMETER :: IN_AU             ='InAU'
  CHARACTER(LEN=*),  PARAMETER :: H_ORDER           ='HOrder'
  CHARACTER(LEN=*),  PARAMETER :: Z_ORDER           ='ZOrder'
  CHARACTER(LEN=*), PARAMETER :: RANDOM_ORDER      ='RandomOrder'
  CHARACTER(LEN=*), PARAMETER :: TRAVEL_ORDER      ='TravelOrder'
  CHARACTER(LEN=*), PARAMETER :: TABLETRAV_ORDER   ='TableTravOrder'
  CHARACTER(LEN=*),  PARAMETER :: NO_ORDER          ='NoOrder'
  ! Clones=
  CHARACTER(LEN=*),  PARAMETER :: CLONES            ='Clones'
  !------------------------------------------------------------------------------
  ! Charge=
  CHARACTER(LEN=*),  PARAMETER :: TOTAL_CHARGE      ='Charge' 
  !------------------------------------------------------------------------------
  ! Multiplicity=
  CHARACTER(LEN=*), PARAMETER :: MULTIPLICITY      ='Multiplicity' 
END MODULE GeometryKeys
