!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 1999, The University of California
!
!    DERIVED TYPES AND COMPOSED OBJECTS
!
MODULE DerivedTypes
   USE GlobalScalars
   USE GlobalCharacters
   IMPLICIT NONE
!==============================================================
!
!  FUNDAMENTAL ARRAY TYPES: VECTOR => RANK FOUR ARRAY
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!  INTEGER VECTOR
!
   TYPE INT_VECT
      INTEGER                      :: Alloc  !-- Allocation key (TRUE, TEMP, FALSE)
#ifdef POINTERS_IN_DERIVED_TYPES
      INTEGER,POINTER,DIMENSION(:) :: I      !-- Vector of integers
#else
      INTEGER,ALLOCATABLE,DIMENSION(:) :: I      !-- Vector of integers
#endif
   END TYPE                                      
!------------------------------------------------------------
!  INTEGER RANK2 ARRAY
!
   TYPE INT_RNK2
      INTEGER                          :: Alloc  !-- Is the array allocated yet?
#ifdef POINTERS_IN_DERIVED_TYPES
      INTEGER, POINTER, DIMENSION(:,:) :: I      !-- Rank 2 array of integers 
#else
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: I      !-- Rank 2 array of integers 
#endif
   END TYPE                                      
!------------------------------------------------------------
!  INTEGER RANK3 ARRAY
!
   TYPE INT_RNK3
      INTEGER                          :: Alloc  !-- Is the array allocated yet?
#ifdef POINTERS_IN_DERIVED_TYPES
      INTEGER,POINTER,DIMENSION(:,:,:) :: I      !-- Rank 3 array of integers 
#else
      INTEGER,ALLOCATABLE,DIMENSION(:,:,:) :: I      !-- Rank 3 array of integers 
#endif
   END TYPE                                      
!------------------------------------------------------------
!  INTEGER RANK4 ARRAY
!
   TYPE INT_RNK4
      INTEGER                            :: Alloc  !-- Allocation key
#ifdef POINTERS_IN_DERIVED_TYPES
      INTEGER,POINTER,DIMENSION(:,:,:,:) :: I      !-- Rank 4 array of integers 
#else
      INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:) :: I      !-- Rank 4 array of integers 
#endif
   END TYPE                                      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!  DOUBLE VECTOR
!
   TYPE DBL_VECT
      INTEGER                           :: Alloc  !-- Allocation key (TRUE, TEMP, FALSE)
#ifdef POINTERS_IN_DERIVED_TYPES
      REAL(DOUBLE),POINTER,DIMENSION(:) :: D      !-- Vector of doubles
#else
      REAL(DOUBLE),ALLOCATABLE,DIMENSION(:) :: D      !-- Vector of doubles
#endif
   END TYPE                                      
!------------------------------------------------------------
!  DOUBLE RANK2 ARRAY
!
   TYPE DBL_RNK2
      INTEGER                              :: Alloc  !-- Is the array allocated yet?
#ifdef POINTERS_IN_DERIVED_TYPES
      REAL(DOUBLE),POINTER,DIMENSION(:,:) :: D      !-- Rank 2 array of doubles 
#else
      REAL(DOUBLE),ALLOCATABLE, DIMENSION(:,:) :: D      !-- Rank 2 array of doubles 
#endif
   END TYPE                                      
!------------------------------------------------------------
!  DOUBLE RANK3 ARRAY
!
   TYPE DBL_RNK3
      INTEGER                               :: Alloc  !-- Is the array allocated yet?
#ifdef POINTERS_IN_DERIVED_TYPES
      REAL(DOUBLE),POINTER,DIMENSION(:,:,:) :: D      !-- Rank 3 array of doubles 
#else
      REAL(DOUBLE),ALLOCATABLE,DIMENSION(:,:,:) :: D      !-- Rank 3 array of doubles 
#endif
   END TYPE                                      
!------------------------------------------------------------
!  DOUBLE RANK4 ARRAY
!
   TYPE DBL_RNK4
      INTEGER                                 :: Alloc  !-- Is the array allocated yet?
#ifdef POINTERS_IN_DERIVED_TYPES
      REAL(DOUBLE),POINTER,DIMENSION(:,:,:,:) :: D      !-- Rank 4 array of doubles 
#else
      REAL(DOUBLE),ALLOCATABLE,DIMENSION(:,:,:,:) :: D      !-- Rank 4 array of doubles 
#endif
   END TYPE                                      
!------------------------------------------------------------
!  DOUBLE RANK6 ARRAY
!
   TYPE DBL_RNK6
      INTEGER                                     :: Alloc  !-- Is the array allocated yet?
#ifdef POINTERS_IN_DERIVED_TYPES
      REAL(DOUBLE),POINTER,DIMENSION(:,:,:,:,:,:) :: D      !-- Rank 6 array of doubles 
#else
      REAL(DOUBLE),ALLOCATABLE,DIMENSION(:,:,:,:,:,:) :: D      !-- Rank 6 array of doubles 
#endif
   END TYPE                                      
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!  CHARACTER STRING
!
   TYPE CHR_SCLR
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: C
   END TYPE                    
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!  CHARACTER STRING VECTOR
!
   TYPE CHR_VECT
      INTEGER                           :: Alloc  !-- Allocation key (TRUE, TEMP, FALSE)
#ifdef POINTERS_IN_DERIVED_TYPES
      CHARACTER(LEN=DEFAULT_CHR_LEN), &
                   POINTER,DIMENSION(:) :: C      !-- Vector of strings
#else
      CHARACTER(LEN=DEFAULT_CHR_LEN), &
                   ALLOCATABLE,DIMENSION(:) :: C      !-- Vector of strings
#endif
   END TYPE                                      
!
!==================================================================================
!
!  OBJECTS DERIVED THROUGH COMPOSITION: BASIS SETS, MATRICES, ETC. 
! 
!==================================================================================
!
#ifdef PARALLEL
!------------------------------------------------------------
!  DISTRIBUTED BLOCK COMPRESSED SPARSE ROW MATRIX
!
   TYPE DBCSR                                       
!     Local 
      INTEGER        :: Alloc  !-- Allocation key
      INTEGER        :: Node   !-- Node from which the data derived
      INTEGER        :: NAtms  !-- Number of atoms
      INTEGER        :: NBlks  !-- Number of non-zero blocks
      INTEGER        :: NNon0  !-- Number of non-zero matrix elements
      TYPE(INT_VECT) :: RowPt  !-- Row index  
      TYPE(INT_VECT) :: ColPt  !-- Coloumn index
      TYPE(INT_VECT) :: BlkPt  !-- Block index  
      TYPE(DBL_VECT) :: MTrix  !-- Blocks
!     Global
      INTEGER        :: GUpDate!-- Has the global information been uptdated?
      TYPE(INT_VECT) :: GRwPt  !-- Global row index  
      TYPE(INT_VECT) :: GClPt  !-- Global col index  
   END TYPE                                      
#endif
!------------------------------------------------------------
!  BLOCK COMPRESSED SPARSE ROW MATRIX
!
   TYPE BCSR                                       
      INTEGER        :: Alloc  !-- Allocation key
      INTEGER        :: NAtms  !-- Number of atoms
      INTEGER        :: NBlks  !-- Number of non-zero blocks
      INTEGER        :: NNon0  !-- Number of non-zero matrix elements
      TYPE(INT_VECT) :: RowPt  !-- Row index  
      TYPE(INT_VECT) :: ColPt  !-- Coloumn index
      TYPE(INT_VECT) :: BlkPt  !-- Block index  
      TYPE(DBL_VECT) :: MTrix  !-- Blocked matrices
   END TYPE                                      
!------------------------------------------------------------
!  Basis Sets
!
   TYPE BSET                                
      INTEGER          :: Alloc   !-- Allocation key
      CHARACTER(LEN= &
      BASESET_CHR_LEN) :: BName   !-- Basis set name
      INTEGER          :: BType   !-- Basis set index
      INTEGER          :: NAtms   !-- Number of atoms
      INTEGER          :: NBasF   !-- Number of basis functions
      INTEGER          :: NKind   !-- Number of atomic kinds (nuclear species)
      INTEGER          :: NCtrt   !-- Max number of contracted functions in a kind
      INTEGER          :: NPrim   !-- Max number of primitives in a contraction
      INTEGER          :: NASym   !-- Max value of angular symmetry (s=0, p=1, ... ) on a BF
      INTEGER          :: LMNLen  !-- Max basis function length 
      TYPE(INT_VECT)   :: Kinds   !-- Atomic kinds or species (Z numbers)
      TYPE(INT_VECT)   :: NCFnc   !-- Number of contracted functions per kind
      TYPE(INT_VECT)   :: BFKnd   !-- Number of basis functions per kind
      TYPE(INT_VECT)   :: LxDex   !-- Basis function index for X
      TYPE(INT_VECT)   :: LyDex   !-- Basis function index for Y
      TYPE(INT_VECT)   :: LzDex   !-- Basis function index for Z
      TYPE(INT_RNK2)   :: NPFnc   !-- Number of primitive functions 
      TYPE(INT_RNK2)   :: LStrt   !-- Starting basis function index per contraction, per kind
      TYPE(INT_RNK2)   :: LStop   !-- Stoping basis function index per contraction, per kind 
      TYPE(INT_RNK3)   :: ASymm   !-- Lo-hi angular symmetry per contraction, per kind
      TYPE(DBL_RNK3)   :: Expnt   !-- Exponent per primitive, per contraction, per kind
      TYPE(DBL_RNK4)   :: CCoef   !-- Contraction coefficient, per symmetry type,
                                  !   per primitive, per contraction, per kind
   END TYPE                                      
!----------------------------------------------------------------
!  Coordinates
!
   TYPE CRDS
!     Status variables
      INTEGER          :: Alloc     !-- Allocation key
      LOGICAL          :: InAU      !-- True if coordinates are in Atomic Units
      INTEGER          :: Ordrd     !-- Reordering key
      INTEGER          :: Confg     !-- Configuration number
!     Electronic coordinates
      INTEGER          :: NElec     !-- Number of electrons
      INTEGER          :: Multp     !-- Total spin multplicity
      REAL(DOUBLE)     :: TotCh     !-- Total charge
      INTEGER          :: NAlph     !-- Number of alpha electrons      
      INTEGER          :: NBeta     !-- Number of beta electrons
!     Misc
      REAL(DOUBLE)     :: ENucN     !-- Nuclear-nuclear repulsion energy
      TYPE(DBL_RNK2)   :: BndBox    !-- Bounding box of the system
#ifdef PERIODIC
!     Perodic Stuff
      LOGICAL              :: AtomW     !-- Wrap atoms back into box--BE CAREFUL
      LOGICAL,DIMENSION(3) :: AutoW     !-- Periodic in X, Y and or Z  direction
      LOGICAL              :: InVecForm !-- What form are the Lattice vectors in
      LOGICAL              :: InAtomCrd !-- Atomic or Fractional Coordinates
      LOGICAL              :: NoTransVec!-- Was the translate calculated or suppied
      TYPE(DBL_VECT)       :: TransVec  !-- Origin Translate Vector
      TYPE(DBL_RNK2)       :: BoxShape  !-- Box Shape Vectors
      TYPE(DBL_RNK2)       :: InvBoxSh  !-- Inverse of the Box Shape Vectors
      TYPE(DBL_RNK2)       :: BoxCarts  !-- Lattice coordinates 
#endif 
!     Atomic coordinates
      INTEGER          :: NAtms     !-- Number of atoms
      INTEGER          :: Nkind     !-- Number of atom kinds or types
      TYPE(INT_VECT)   :: AtNum     !-- Atomic number per atom      
      TYPE(INT_VECT)   :: AtTyp     !-- Atom type or kind per atom 
      TYPE(DBL_RNK2)   :: Carts     !-- Cartesian coordinates 
   END TYPE 
!-------------------------------------------------------------------------------------
!  Density in a Hermite Gaussian basis
!
  TYPE HGRho
     INTEGER          :: Alloc   !-- Allocation key
     INTEGER          :: AllocRE !-- Allocation key for RhoEst
     INTEGER          :: NExpt   !-- Number of exponents
     INTEGER          :: NDist   !-- Number of distributions 
     INTEGER          :: NCoef   !-- Number of coefficients
     TYPE(INT_VECT)   :: NQ      !-- Number of distributions per exponent (NExpt)
     TYPE(INT_VECT)   :: Lndx    !-- Max angular symmetry per exponent (NExpt)
     TYPE(INT_VECT)   :: OffQ    !-- Distribution offset (NExpt)
     TYPE(INT_VECT)   :: OffR    !-- Coefficient offset  (NExpt)
     TYPE(DBL_VECT)   :: Expt    !-- Exponents (NExpt)
     TYPE(DBL_VECT)   :: Qx      !-- x-coordinates of distribution center (NDist)
     TYPE(DBL_VECT)   :: Qy      !-- y-coordinates of distribution center (NDist)
     TYPE(DBL_VECT)   :: Qz      !-- z-coordinates of distribution center (NDist)
     TYPE(DBL_VECT)   :: Est     !-- Value of the Norm <rho_exp(r)|rho_exp(r)> (NDist)
     TYPE(DBL_VECT)   :: Co      !-- Density coefficients (NCoef)
  END TYPE HGRho
!--------------------------------------------------------------------------
!  Bounding Box Type
!
   TYPE BBox
      INTEGER                            :: Tier     ! Level of this box
      INTEGER                            :: Number   ! Box number
      REAL(DOUBLE),   DIMENSION(3,2)     :: BndBox   ! Upper and lower limits of the box
      REAL(DOUBLE),   DIMENSION(3)       :: Center   ! Center of the box
      REAL(DOUBLE),   DIMENSION(3)       :: Half     ! Half each dimensions length
   END TYPE 
!--------------------------------------------------------------------------
! Atom Pair Type
!
  TYPE AtomPair
     REAL(DOUBLE),DIMENSION(3)      :: A,B
     REAL(DOUBLE)                   :: AB2
     INTEGER                        :: KA
     INTEGER                        :: KB
     INTEGER                        :: NA
     INTEGER                        :: NB
     LOGICAL                        :: SameAtom
  ENDTYPE AtomPair
!--------------------------------------------------------------------------
! Primative Atom Pair Type
!
  TYPE PrimPair
     INTEGER                        :: Ell  
     REAL(DOUBLE)                   :: Z
     REAL(DOUBLE),DIMENSION(3)      :: P
     REAL(DOUBLE),DIMENSION(MaxLen) :: Ket
     REAL(DOUBLE)                   :: Extent
     TYPE(BBox)                     :: Box
  ENDTYPE PrimPair
#ifdef PERIODIC
!-----------------------------------------------------------
! The Set of Cells needed to sum over: CellSet
!
  TYPE CellSet
     INTEGER                        :: Alloc
     INTEGER                        :: NCells
     TYPE(DBL_RNK2)                 :: CellCarts
  ENDTYPE CellSet
#endif
!------------------------------------------------------------
!  Numerical thresholds
!
   TYPE TOLS                                       
      REAL(DOUBLE) :: Cube  !-- Cubature 
      REAL(DOUBLE) :: Dist  !-- Distribution threshold
      REAL(DOUBLE) :: TwoE  !-- Two electron integral threshold
      REAL(DOUBLE) :: Trix  !-- Matrix threshold
   END TYPE 
!------------------------------------------------------------
!  Debuging flags
!
   TYPE DEBG                                       
      INTEGER :: Key  !-- Debug level
      INTEGER :: Mat  !-- Debug matrices      
      INTEGER :: Set  !-- Debug basis set
      INTEGER :: Int  !-- Debug integrals
      INTEGER :: Rho  !-- Debug density
      INTEGER :: Fmt  !-- Debug formating
   END TYPE                                      
!------------------------------------------------------------
!  Timing/performance statistics
!
   TYPE TIME                                       
      REAL(DOUBLE) :: CPUS   !-- Accumulated CPU sec 
      REAL(DOUBLE) :: Wall   !-- Accumulated wall sec
      REAL(DOUBLE) :: CStrt  !-- Start of cpu  sec accumulation
      REAL(DOUBLE) :: WStrt  !-- Start of wall sec accumulation
      REAL(DOUBLE) :: FLOP   !-- Floating point opperations
   END TYPE                                      
!------------------------------------------------------------
!  Memory statistics
!
   TYPE MEMS
      INTEGER :: Allocs    !-- Number of ALLOCATEs untill now
      INTEGER :: DeAllocs  !-- Number of DEALLOCATEs untill now
      INTEGER :: MemTab    !-- Current total number of bytes allocated
      INTEGER :: MaxMem    !-- Max number of total bytes allocated 
      INTEGER :: MaxAlloc  !-- Max number of bytes allocated by a call to New
   END TYPE                                      
!------------------------------------------------------------
!  Argument lists
!
   TYPE ARGMT
      INTEGER        :: Alloc     !-- Allocation key
      INTEGER        :: NI        !-- Length of integer vector
      INTEGER        :: NC        !-- Length of character vector
      TYPE(INT_VECT) :: I         !-- Vector of integer arguments
      TYPE(CHR_VECT) :: C         !-- Vector of character arguments
   END TYPE
!------------------------------------------------------------
!  MPI Derived data type pointers
!
   TYPE MPI_INDX
      INTEGER        :: Alloc     !-- Allocation key
      INTEGER        :: Type      !-- Derived type
      INTEGER        :: NBlks     !-- Number of blocks in type
      TYPE(INT_VECT) :: Disp      !-- Vector of displacements
      TYPE(INT_VECT) :: Blks      !-- Vector of block sizes
   END TYPE
!
END MODULE



