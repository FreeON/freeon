MODULE ERIGlobals
  USE DerivedTypes
  USE GlobalScalars
  IMPLICIT NONE
  REAL(DOUBLE),DIMENSION(1:HGLen) :: HGKet
  REAL(DOUBLE),DIMENSION(0:SPLen) :: SPKetC
  REAL(DOUBLE),DIMENSION(0:SPLen) :: SPKetS
!  REAL(DOUBLE),DIMENSION(500)     :: W
END MODULE ERIGlobals

MODULE PoleGlobals
  USE Derivedtypes
  USE GlobalScalars
  IMPLICIT NONE
  ! This is the dynamically set expansion order and tensor length for multipole expansions on the tree
  INTEGER                                      :: MaxPoleEll,LenPoleEll
  ! This is the dynamically set expansion order and tensor length for cell multipole expansions
  INTEGER                                      :: MaxPFFFEll,LenPFFFEll
  ! This is the fixed expansion order used for the Factorial function, scratch memory
  ! and other precomputed values
  INTEGER,PARAMETER                            :: FFELL=64
  INTEGER,PARAMETER                            :: FFELL2=2*FFELL
  INTEGER,PARAMETER                            :: FFLen=FFEll*(FFEll+3)/2
  INTEGER,PARAMETER                            :: FFLen2=FFEll2*(FFEll2+3)/2!
  ! Globally available precomputed factorial functions etc
  REAL(DOUBLE), DIMENSION(0:2*FFEll2)          :: Factorial
  REAL(DOUBLE), DIMENSION(0:FFEll2)            :: FactOlm0,FactMlm0
  REAL(DOUBLE), DIMENSION(0:FFLen2)            :: FactOlm2,FactMlm2
  ! The binomial distribution
  REAL(DOUBLE), DIMENSION(0:SPEll+1,0:FFELL)   :: FudgeFactorial
  ! Scratch space for tensor manipulations
  REAL(DOUBLE), DIMENSION(0:FFEll2)            :: Sine,Cosine,CoFact
  REAL(DOUBLE), DIMENSION(0:FFLen2)            :: ALegendreP,Spq,Cpq
  REAL(DOUBLE), DIMENSION(0:FFLen2,3)          :: DCpq,DSpq
END MODULE

MODULE RhoList

  USE Derivedtypes
  USE GlobalScalars
  IMPLICIT NONE
  TYPE HGLink
     INTEGER                   :: Ell      ! The highest angular symmetry in this cluster
     REAL(DOUBLE)              :: Zeta     ! The true zeta for the distributions
     REAL(DOUBLE),DIMENSION(3) :: Cent     ! This is the center of the distribution
     TYPE(DBL_VECT)            :: Coef     ! Coefficients of the HGTF density
     TYPE(HGLink),POINTER      :: Next     ! Next link in the chain
  END TYPE HGLink

END MODULE RhoList



MODULE PoleNodeType
  USE Derivedtypes
  USE GlobalScalars
  IMPLICIT NONE
  TYPE Pole
     INTEGER                               :: Ell      ! This is the ell of the multipole expansion
     INTEGER                               :: Complex  ! The distributional complexity represented by this expansion
     REAL(DOUBLE)                          :: Charge   ! Total nuclear charge used for wheighting translation distances
     REAL(DOUBLE)                          :: Delta    ! For debugging
     REAL(DOUBLE),DIMENSION(3)             :: Center   ! This is the center of the multipole expansion


#ifdef POINTERS_IN_DERIVED_TYPES
     REAL(DOUBLE),DIMENSION(:),POINTER     :: S        ! Im component of the multipole tensor
     REAL(DOUBLE),DIMENSION(:),POINTER     :: C        ! Re component of the multipole tensor
#else
     REAL(DOUBLE),DIMENSION(:),ALLOCATABLE :: S        ! Im component of the multipole tensor
     REAL(DOUBLE),DIMENSION(:),ALLOCATABLE :: C        ! Re component of the multipole tensor
#endif
  END TYPE Pole
!
  TYPE Herm
     INTEGER                               :: Ell      ! The highest angular symmetry in this cluster
     INTEGER                               :: Stack    ! Integer index of HG data in global stack array
#ifdef POINTERS_IN_DERIVED_TYPES
     INTEGER,DIMENSION(:),POINTER          :: NQ       ! Number of distributions per Ell
#else
     INTEGER,DIMENSION(:),ALLOCATABLE      :: NQ       ! Number of distributions per Ell
#endif
     TYPE(DBL_VECT),POINTER,DIMENSION(:)   :: Zeta     ! The true zeta for the distributions
     TYPE(DBL_VECT),POINTER,DIMENSION(:)   :: IHlf     ! Integral bound for use in Schwartz inequality
     TYPE(DBL_RNK2),POINTER,DIMENSION(:)   :: Cent     ! This is the center of the distribution
     TYPE(DBL_RNK2),POINTER,DIMENSION(:)   :: Coef     ! Coefficients of the HGTF density
  END TYPE Herm
  !
  TYPE ScalarHerm
     INTEGER      :: Ell      ! The highest angular symmetry in this cluster
     REAL(DOUBLE) :: Zeta     ! The true zeta for the distributions
     REAL(DOUBLE) :: IHlf     ! Integral bound for use in Schwartz inequality
     REAL(DOUBLE) :: Cent     ! This is the center of the distribution
#ifdef POINTERS_IN_DERIVED_TYPES
     REAL(DOUBLE),DIMENSION(:),POINTER     :: Coef     ! Coefficients of the HGTF density
#else
     REAL(DOUBLE),DIMENSION(:),ALLOCATABLE :: Coef     ! Coefficients of the HGTF density
#endif
  END TYPE ScalarHerm
!
  TYPE MAC
     REAL(DOUBLE),DIMENSION(0:7)           :: O
     REAL(DOUBLE)                          :: Delta    ! Max radius from box corners to expansion center
  END TYPE MAC
!
  TYPE PAC
     REAL(DOUBLE)                          :: Zeta     ! Effective Gaussian decay for PAC estimation
     REAL(DOUBLE)                          :: Wght     ! Effective weight for PAC estimation
  END TYPE PAC

  TYPE PAC8
     REAL(DOUBLE)                           :: Zeta     ! Effective Gaussian decay for PAC estimation
     REAL(DOUBLE),DIMENSION(-1:1,-1:1,-1:1) :: Wght     ! Effective weight for PAC estimation
  END TYPE PAC8
!
  TYPE QCPrim
     TYPE(PrimPair) :: Prim
     TYPE(MAC)      :: MAC
     TYPE(BBox)                            :: Box      !
!     TYPE(PAC)      :: PAC
!     REAL(DOUBLE)   :: IHalf
  END TYPE QCPrim

  TYPE PoleNode
     LOGICAL                               :: Leaf      ! Is this a leaf node (containing just one atom) ?
     INTEGER                               :: BdexE     ! Begining index of ORB charge list for this node
     INTEGER                               :: EdexE     ! Ending index of ORB charge list for this node
     INTEGER                               :: BdexN     ! Begining index of ORB nuclei list for this node
     INTEGER                               :: EdexN     ! Ending index of ORB nuclei for this node
     INTEGER                               :: NQ        ! Number of charges (distributions) bounded by this node
     INTEGER                               :: NAtms     ! Number of nuclear centers bounded by this node
     REAL(DOUBLE)                          :: IHalf     ! This is the max Schwartz factor (ab|ab)**0.5 used in the
                                                        ! locally direct part of the code
     REAL(DOUBLE)                          :: IHMin     ! This is the cooresponding min factor for the node
     TYPE(MAC)                             :: MAC
     TYPE(PAC)                             :: PAC
     TYPE(Pole)                            :: Pole
     TYPE(Herm)                            :: Herm
!
     TYPE(BBox)                            :: Box      ! Bounding Box of distribution (for PAC)
     TYPE(PoleNode),POINTER                :: Left
     TYPE(PoleNode),POINTER                :: Right
     TYPE(PoleNode),POINTER                :: Next
  END TYPE PoleNode

END MODULE PoleNodeType
