MODULE PoleNodeType
  USE Derivedtypes
  USE GlobalScalars   
  IMPLICIT NONE
  !=================================================================================================
  !  Hierarchical density node
  !=================================================================================================
  TYPE PoleNode
     LOGICAL                               :: Leaf     ! Is this a data containing node?
     INTEGER                               :: Bdex     ! Begining index of ORB list for this node
     INTEGER                               :: Edex     ! ENDign index of ORB list for this node
     INTEGER                               :: NQ       ! Number of centers
     INTEGER                               :: Ell      ! Ell type
     INTEGER                               :: EllCD    ! Maximium Ell of the Distributuions in the box
     REAL(DOUBLE)                          :: Zeta     ! Minimum exponent in this node
     REAL(DOUBLE)                          :: Strength ! Strength of the Pole
     REAL(DOUBLE)                          :: DMax2    ! (Max distance)^2 from node center to dist
     REAL(DOUBLE)                          :: WCoef    ! Weight for the Gaussian in the Box 
     TYPE(BBox)                            :: Box      ! Bounding Box of distribution (for PAC)
     TYPE(PoleNode),POINTER                :: Descend  ! Next node in tree descent
     TYPE(PoleNode),POINTER                :: Travrse  ! Next node in tree traversal
#ifdef POINTERS_IN_DERIVED_TYPES
     REAL(DOUBLE),DIMENSION(:),POINTER     :: S        ! Im component of the multipole tensor
     REAL(DOUBLE),DIMENSION(:),POINTER     :: C        ! Re component of the multipole tensor
     REAL(DOUBLE),DIMENSION(:),POINTER     :: Co       ! Coefficients of the HGTF density
#else
     REAL(DOUBLE),DIMENSION(:),ALLOCATABLE :: S        ! Im component of the multipole tensor
     REAL(DOUBLE),DIMENSION(:),ALLOCATABLE :: C        ! Re component of the multipole tensor
     REAL(DOUBLE),DIMENSION(:),ALLOCATABLE :: Co       ! Coefficients of the HGTF density
#endif
  END TYPE PoleNode
END MODULE PoleNodeType
