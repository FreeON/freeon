!-----------------------------------------------------------------
!    DERIVED TYPES AND COMPOSED OBJECTS
!    Author:  Matt Challacombe
!-----------------------------------------------------------------
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
!--------------------------------------------------------------------
!
   TYPE CHR10_VECT
      INTEGER                           :: Alloc  !-- Allocation key (TRUE, TEMP, FALSE)
#ifdef POINTERS_IN_DERIVED_TYPES
      CHARACTER(LEN=10), &
                   POINTER,DIMENSION(:) :: C      !-- Vector of strings
#else
      CHARACTER(LEN=10), &
                   ALLOCATABLE,DIMENSION(:) :: C      !-- Vector of strings
#endif
   END TYPE
!
!----------------------------------------------------------------------
!
   TYPE LOG_VECT
      INTEGER                           :: Alloc  !-- Allocation key (TRUE, TEMP, FALSE)
#ifdef POINTERS_IN_DERIVED_TYPES
      LOGICAL,POINTER,DIMENSION(:) :: L      !-- Vector of logicals
#else
      LOGICAL,ALLOCATABLE,DIMENSION(:) :: L      !-- Vector of logicals
#endif
   END TYPE
!
!==================================================================================
!
!  OBJECTS DERIVED THROUGH COMPOSITION: BASIS SETS, MATRICES, ETC.
!
!==================================================================================
!------------------------------------------------------------
!  DISTRIBUTED BLOCK COMPRESSED SPARSE ROW MATRIX
!
   TYPE DBCSR
!     Local
      INTEGER        :: NSMat  !-- Number of spin matrices
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
!------------------------------------------------------------
!  BLOCK COMPRESSED SPARSE ROW MATRIX
!
   TYPE BCSR
      INTEGER        :: NSMat  !-- Number of spin matrices
      !TYPE(INT_VECT) :: SInfo  !--
      INTEGER        :: Alloc  !-- Allocation key
      INTEGER        :: NAtms  !-- Number of atoms
      INTEGER        :: NBlks  !-- Number of non-zero blocks
      INTEGER        :: NNon0  !-- Number of non-zero matrix elements
      TYPE(INT_VECT) :: RowPt  !-- Row index
      TYPE(INT_VECT) :: ColPt  !-- Coloumn index
      TYPE(INT_VECT) :: BlkPt  !-- Block index
      TYPE(DBL_VECT) :: MTrix  !-- Blocked matrices
   END TYPE BCSR
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
      TYPE(CHR10_VECT) :: AtNam   !-- Atomname
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
      LOGICAL          :: HasECPs
      INTEGER          :: MxProjL
      INTEGER          :: Typ1Fnk
      INTEGER          :: Typ2Fnk
      TYPE(DBL_VECT)   :: NCoreEl
      ! Type one ECP data
      TYPE(INT_VECT)   :: NTyp1PF
      TYPE(INT_RNK2)   :: Typ1Ell
      TYPE(DBL_RNK2)   :: Typ1Exp
      TYPE(DBL_RNK2)   :: Typ1CCo
      ! Type two ECP data
      TYPE(INT_VECT)   :: ProjEll
      TYPE(INT_RNK2)   :: NTyp2PF
      TYPE(INT_RNK3)   :: Typ2Ell
      TYPE(DBL_RNK3)   :: Typ2Exp
      TYPE(DBL_RNK3)   :: Typ2CCo
   END TYPE
!------------------------------------------------------------------------------------
! The Set of Cells needed to sum over: CellSet
!
  TYPE CellSet
     INTEGER                        :: Alloc
     INTEGER                        :: NCells
     REAL(DOUBLE)                   :: Radius
     TYPE(DBL_RNK2)                 :: CellCarts
  ENDTYPE CellSet
!------------------------------------------------------------------------------------
! Periodic Information
!
  TYPE PBCInfo
     INTEGER                     :: Alloc      !-- Allocation key
     INTEGER                     :: Dimen      !-- Dimension of the System
     INTEGER                     :: PFFMaxEll  !-- Maxium angular symmetry in the PFF contribution
     INTEGER                     :: PFFWelSep  !-- Well seperated criteria in the PFF, using strict FMM definition
!!     Depricated:
!!     LOGICAL                     :: AtomW      !-- Wrap atoms back into box--BE CAREFUL
     LOGICAL                     :: InVecForm  !-- What form are the Lattice vectors in
     LOGICAL                     :: InAtomCrd  !-- Atomic or Fractional Coordinates
     LOGICAL                     :: Translate  !-- Should the Atomic Coordinated be Translated
     REAL(DOUBLE)                :: CellVolume !-- Cell Volume
     REAL(DOUBLE)                :: Epsilon    !-- Epsilon at Infinity (Metal == Infinity)
     REAL(DOUBLE)                :: DipoleFAC  !-- Normalization of the Dipole Term
     REAL(DOUBLE)                :: QupoleFAC  !-- Normalization of Quadrupole Term
     TYPE(INT_VECT)              :: AutoW      !-- Periodic in X, Y and or Z  direction
     TYPE(INT_VECT)              :: SuperCell  !-- SuperCell of Primative Cell
     TYPE(DBL_VECT)              :: CellCenter !-- Center of Cell
     TYPE(DBL_VECT)              :: TransVec   !-- Origin Translate Vector
     TYPE(DBL_RNK2)              :: BoxShape   !-- Box Shape Vectors
     TYPE(DBL_RNK2)              :: InvBoxSh   !-- Inverse of the Box Shape Vectors
     TYPE(DBL_RNK2)              :: LatFrc     !-- Lattice Forces
  END TYPE PBCInfo
!----------------------------------------------------------------------
!
   TYPE BONDDATA
     INTEGER        :: Alloc     !-- Allocation key
     INTEGER        :: N
     TYPE(INT_RNK2) :: IJ
     TYPE(DBL_VECT) :: Length
     TYPE(CHR10_VECT) :: Type
     TYPE(INT_VECT) :: HBExtraSN !serial # of third atom in HBond
     TYPE(INT_VECT) :: HBExtraNC !nuclear charge of third at in HBond
     TYPE(INT_VECT) :: LonelyAtom   ! 1 for lonely, 0 otherwise
   END TYPE BONDDATA
!
!----------------------------------------------------------------------
!
   TYPE ATOMBONDS
     INTEGER        :: Alloc     !-- Allocation key
     INTEGER        :: N1
     INTEGER        :: N2
     TYPE(INT_VECT) :: Count
     TYPE(INT_RNK2) :: Bonds
     TYPE(INT_RNK2) :: Atoms
   END TYPE ATOMBONDS
!
   TYPE INTC
     INTEGER          :: Alloc     !-- Allocation key
     INTEGER          :: N
     TYPE(CHR10_VECT) :: Def
     TYPE(INT_RNK2)   :: Atoms
     TYPE(INT_RNK2)   :: Cells
     TYPE(DBL_VECT)   :: Value
     TYPE(LOG_VECT)   :: Constraint
     TYPE(DBL_VECT)   :: ConstrValue
     TYPE(LOG_VECT)   :: Active
     TYPE(DBL_VECT)   :: PredVal
     TYPE(DBL_VECT)   :: PredGrad
     TYPE(DBL_VECT)   :: InvHess
   END TYPE INTC
!
!----------------------------------------------------------------------
!
   TYPE PBCFits
     INTEGER        :: Alloc     !-- Allocation key
     INTEGER        :: MaxMem
     INTEGER        :: ActMem
     TYPE(DBL_VECT) :: AWeights
     TYPE(DBL_RNK2) :: PBCValues
     TYPE(DBL_RNK2) :: PBCGrads
   END TYPE PBCFits
!
!----------------------------------------------------------------------
!
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
      REAL(DOUBLE)     :: ETotal    !-- Total SCF Energy at this geometry
      TYPE(DBL_VECT)   :: ETotalPerSCF  !-- Total SCF energy per SCF cycle
      REAL(DOUBLE)     :: GradRMS   !-- RMS error in gradient at this geometry
      REAL(DOUBLE)     :: GradMax   !-- Max error in gradient at this geometry
      LOGICAL          :: Unstable  !-- SCF is unstable at this geometry
      TYPE(DBL_RNK2)   :: BndBox    !-- Bounding box of the system
!     Perodic Stuff
      TYPE(PBCInfo)    :: PBC       !-- Periodic Information
      TYPE(CellSet)    :: OvCells   !-- Overlap (basis set) cell set
      TYPE(CellSet)    :: InCells   !-- Inner (Coulomb) cell set
!     Atomic coordinates
      INTEGER          :: NAtms     !-- Number of atoms
      INTEGER          :: NKind     !-- Number of atom kinds or types
      TYPE(DBL_VECT)   :: AtNum     !-- Atomic number per atom
      TYPE(INT_VECT)   :: AtTyp     !-- Atom type or kind per atom
      TYPE(CHR10_VECT) :: AtNam     !-- Atomname
!!!!      TYPE(CHR10_VECT) :: AtMMTyp   !-- Molecular Mechanics atomtype
      TYPE(DBL_VECT)   :: AtMss     !-- Atomic Mass per Atom
      TYPE(INT_VECT)   :: CConstrain!-- Atom type or kind per atom
      TYPE(INT_VECT)   :: DoFreq    !-- Atom included in the frequency calculation
      TYPE(DBL_RNK2)   :: Carts     !-- Cartesian coordinates
      TYPE(DBL_RNK2)   :: BoxCarts  !-- Fractional coordinates
      TYPE(DBL_RNK2)   :: Velocity  !-- velocities
      TYPE(DBL_RNK2)   :: Gradients !-- Gradients
      TYPE(DBL_RNK2)   :: Fext      !-- An External Force
      ! Used in QUICCA:
      TYPE(DBL_RNK2)   :: Displ     !-- Displaced Cartesian geometry
      TYPE(DBL_RNK2)   :: PBCDispl  !-- Displaced Box Shape
      LOGICAL          :: LatticeOnly
      INTEGER          :: AltCount
   END TYPE
!-------------------------------------------------------------------------------------
!  Cartesian Multipoles of the Density
!
  TYPE CMPoles
     INTEGER          :: Alloc   !-- Allocation key
     REAL(DOUBLE)     :: MPole   !-- Monopole
     TYPE(DBL_VECT)   :: DPole   !-- Dipoles (Dim=3)
     TYPE(DBL_VECT)   :: QPole   !-- Quadrupoles (Dim=6)
  END TYPE CMPoles
!-------------------------------------------------------------------------------------
!  Density in a Hermite Gaussian basis
!
  TYPE HGRho
     INTEGER          :: Alloc   !-- Allocation key
     INTEGER          :: NSDen   !-- Number of Spin densities
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
     TYPE(DBL_VECT)   :: Co      !-- Density coefficients (NCoef)
  END TYPE HGRho
!-------------------------------------------------------------------------------------
!  ONX distribution buffers
!
  TYPE DBuf
     INTEGER          :: Alloc   !-- Allocation key
     INTEGER          :: NShells
     INTEGER          :: NPrim
     INTEGER          :: MInfo
     INTEGER          :: MAXDis
     INTEGER          :: MAXPrm
     INTEGER          :: MAXD    !-- Max dist in a l/c set: TBufC,TBufP
     INTEGER          :: MAXT    !-- Max angular symmetry types: TCode,TCPop
     INTEGER          :: MAXK    !-- Max contraction types: CCode,TCPop
     INTEGER          :: MAXC    !-- Contracted data set size: TBufC
     INTEGER          :: MAXP    !-- Primitive data set size: TBufP
     INTEGER          :: LenCC
     INTEGER          :: LenTC
     TYPE(INT_VECT)   :: TCode   !-- Angular symmetry types
     TYPE(INT_VECT)   :: CCode   !-- Contraction lengths
     TYPE(INT_RNK2)   :: TCPop   !-- Is this TCode and CCode populated?
     TYPE(DBL_RNK2)   :: TBufC   !-- Temp. contracted dis. buffer for 2-e estimates
     TYPE(DBL_RNK3)   :: TBufP   !-- Temp. primitive dis. buffer for 2-e estimates
     TYPE(INT_RNK4)   :: DisPtr  !-- Pointers to locations in DisBuf and PrmBuf
     TYPE(DBL_VECT)   :: DisBuf  !-- Contracted dis. buffer for 2-e computation
     TYPE(DBL_VECT)   :: PrmBuf  !-- Primitive dis. buffer for 2-e computation
  END TYPE DBuf
!-------------------------------------------------------------------------------------
!  ONX SL info
!
  TYPE DSL
     INTEGER          :: Alloc   !-- Allocation key
     INTEGER          :: MAXSL   !-- Maximum size of the SL buffer
     TYPE(INT_VECT)   :: SLDis   !-- Pointer to distribution info
     TYPE(INT_VECT)   :: SLPrm   !-- Pointer to primitive info
     TYPE(INT_VECT)   :: SLKey   !-- Basis set offset key (for 2-e symmetry)
  END TYPE DSL
!-------------------------------------------------------------------------------------
!  ONX integral buffers
!
  TYPE IBuf
     INTEGER          :: Alloc   !-- Allocation key
     INTEGER          :: NPrim
     INTEGER          :: MAXI    !-- Size of W1 and W2
     INTEGER          :: Lval    !-- Angular symmetry of Gamma table in memory
     INTEGER          :: MAXL    !-- Max L in the Gamma tables
     INTEGER          :: MaxInts !-- Maximum number of two-e in a vector loop
     TYPE(DBL_VECT)   :: W1      !-- Two-e scratch space
     TYPE(DBL_VECT)   :: W2      !-- Two-e scratch space
     TYPE(DBL_RNK2)   :: CB      !-- Bra contraction coefficients
     TYPE(DBL_RNK3)   :: CK      !-- Ket contraction coefficients
     TYPE(DBL_RNK2)   :: WR
     TYPE(DBL_RNK2)   :: WZ
  END TYPE IBuf
!-------------------------------------------------------------------------------------
!  ONX integral drivers
!
  TYPE IDrv
     INTEGER          :: Alloc   !-- Allocation key
     INTEGER          :: LngVRR
     INTEGER          :: LngLoc
     INTEGER          :: LngDrv
     INTEGER          :: LngCC
     INTEGER          :: id,is,nr,ns
     TYPE(INT_VECT)   :: VLOC
     TYPE(INT_VECT)   :: CDrv
     TYPE(INT_RNK2)   :: SLOC
  END TYPE IDrv
!-------------------------------------------------------------------------------------
!  ONX integral space
!
  TYPE ISpc
    INTEGER           :: L1,L2,L3,L4
    INTEGER           :: NB1,NB2
    INTEGER           :: NK1,NK2
    INTEGER           :: NVRR
  END TYPE ISpc
!-------------------------------------------------------------------------------------
!  ONX gradient driver
!
  TYPE GradD
     INTEGER          :: Alloc   !-- Allocation key
     INTEGER          :: LG1,LG2,LG3,LG4,LG5
     INTEGER          :: LenDa,LenDc,LenDb,LenDn
     INTEGER          :: NCON
     INTEGER          :: NLOCB1,NLOCB2,NLOCB3,NLOCK2,NLOCK3
     TYPE(INT_RNK2)   :: GDrv1,GDrv2,GDrv3,GDrv4,GDrv5
  END TYPE GradD
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
! Primitive Pair Type
!
  TYPE PrimPair
     INTEGER                   :: Ell,KA,KB,CFA,CFB,PFA,PFB
     REAL(DOUBLE)              :: ZA,ZB,Zeta,Xi,AB2
     REAL(DOUBLE),DIMENSION(3) :: A,B,P,Pw
  END TYPE
!--------------------------------------------------------------------------
! Atom Pair Type
!
  TYPE AtomPair
     REAL(DOUBLE),DIMENSION(3) :: A,B
     REAL(DOUBLE)              :: AB2
     INTEGER                   :: KA,KB,NA,NB
     LOGICAL                   :: SameAtom
  ENDTYPE AtomPair
!------------------------------------------------------------
!  Numerical thresholds
!
   TYPE TOLS
      REAL(DOUBLE) :: Cube  !-- Cubature
      REAL(DOUBLE) :: Dist  !-- Distribution threshold
      REAL(DOUBLE) :: TwoE  !-- Two electron integral threshold
      REAL(DOUBLE) :: Trix  !-- Matrix threshold
      REAL(DOUBLE) :: ETol  !-- Relative error in total energy sought
      REAL(DOUBLE) :: DTol  !-- Max difference in density matrix sought
   END TYPE
!------------------------------------------------------------
!  Debuging flags
!
   TYPE DEBG
      INTEGER :: Key        !-- Debug level
      INTEGER :: Chk        !-- Debug check sums
      INTEGER :: Mat        !-- Debug matrices
      INTEGER :: Set        !-- Debug basis set
      INTEGER :: Int        !-- Debug integrals
      INTEGER :: Rho        !-- Debug density
      INTEGER :: Fmt        !-- Debug formating
      INTEGER :: GeOp       !-- Debug geometry optimizer
      INTEGER :: GeOpSpInfo !-- Debug sparsity info in geometry optimizer
      INTEGER :: MM         !-- Debug molecular mechanics energies and forces
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
!---------------------------------------------------------------
!
   TYPE BMATR
     INTEGER        :: Alloc     !-- Allocation key
     TYPE(INT_RNK2) :: IB
     TYPE(DBL_RNK2) :: B
     TYPE(INT_VECT) :: BLI
     TYPE(DBL_RNK2) :: BL
   END TYPE BMATR
!
   TYPE Cholesky
     INTEGER        :: Alloc     !-- Allocation key
     TYPE(DBL_VECT) :: GcScale
     TYPE(INT_VECT) :: Perm
     TYPE(INT_VECT) :: IPerm
     TYPE(INT_VECT) :: ChRowPt
     TYPE(INT_VECT) :: ChColPt
     TYPE(DBL_VECT) :: ChDiag
     TYPE(DBL_VECT) :: ChFact
   END TYPE
!
!----------------------------------------------------------------------
!
   TYPE Sp1x1
     INTEGER        :: Alloc     !-- Allocation key
     INTEGER        :: NRow
     INTEGER        :: NZ
     TYPE(INT_VECT) :: IA
     TYPE(INT_VECT) :: JA
     TYPE(DBL_VECT) :: AN
   END TYPE Sp1x1
!
!----------------------------------------------------------------------
!
   TYPE TOPMATS
     TYPE(INT_VECT) :: ITot12,JTot12
     TYPE(INT_VECT) :: IExcl,JExcl
   END TYPE TOPMATS
!
!----------------------------------------------------------------------
!
   TYPE IntCBox
     INTEGER        :: Alloc     !-- Allocation key
     INTEGER        :: N,NX,NY,NZ,NatmsLoc
     TYPE(INT_VECT) :: I
     TYPE(INT_VECT) :: J
   END TYPE IntCBox
!
!----------------------------------------------------------------------
!
   TYPE OUTPDATA
     INTEGER        :: Alloc     !-- Allocation key
     INTEGER        :: N
     TYPE(INT_RNK2) :: IJKL
     TYPE(CHR10_VECT) :: Type
   END TYPE OUTPDATA
!
!----------------------------------------------------------------------
!
   TYPE ANGLEDATA
     INTEGER        :: Alloc     !-- Allocation key
     INTEGER        :: N
     TYPE(INT_RNK2) :: IJK
     TYPE(CHR10_VECT) :: Type
   END TYPE ANGLEDATA
!
!----------------------------------------------------------------------
!
END MODULE

