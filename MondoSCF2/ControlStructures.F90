MODULE ControlStructures
  USE GlobalScalars
  USE DerivedTypes
  USE BasisSetParameters
!
  INTEGER, PARAMETER                :: MaxSets=6
!
  TYPE FileNames
     INTEGER                        :: NewFileID
     INTEGER                        :: OldFileID
     CHARACTER(LEN=DCL)             :: M_PWD
     CHARACTER(LEN=DCL)             :: M_HOME
     CHARACTER(LEN=DCL)             :: M_SCRATCH
     CHARACTER(LEN=DCL)             :: SCF_NAME
     CHARACTER(LEN=DCL)             :: IFile
     CHARACTER(LEN=DCL)             :: OFile
     CHARACTER(LEN=DCL)             :: LFile
     CHARACTER(LEN=DCL)             :: GFile
     CHARACTER(LEN=DCL)             :: HFile
     CHARACTER(LEN=DCL)             :: RFile
  END TYPE FileNames
!
  TYPE Options
     INTEGER                        :: NMthds
     INTEGER                        :: NModls
     INTEGER                        :: NThrsh
     INTEGER                        :: NSteps
     INTEGER                        :: Guess
     INTEGER                        :: Grad
     INTEGER                        :: Coordinates
     INTEGER,   DIMENSION(MaxSets)  :: Methods
     INTEGER,   DIMENSION(MaxSets)  :: Models
     INTEGER,   DIMENSION(MaxSets)  :: AccuracyLevels
     LOGICAL                        :: OneBase,DoGDIIS
     TYPE(INT_VECT)                 :: RestartState
     TYPE(TOLS),DIMENSION(MaxSets)  :: Thresholds
     TYPE(DEBG)                     :: PFlags
  END TYPE Options

  TYPE Parallel
     INTEGER                        :: NProc
     INTEGER                        :: NTime
     INTEGER                        :: NSpace
     CHARACTER(LEN=DCL)             :: MPIRun
  END TYPE Parallel

  TYPE Dynamics
     INTEGER                         :: MDAlgorithm
     INTEGER                         :: CRDfreq
     INTEGER                         :: VELfreq
     INTEGER                         :: ENEfreq
     INTEGER                         :: RESfreq
     INTEGER                         :: MAX_STEPS
     INTEGER                         :: TRANSfreq
     INTEGER                         :: ROTATfreq 
     REAL(DOUBLE)                    :: DT
     REAL(DOUBLE)                    :: TEMP0
     REAL(DOUBLE)                    :: TEMP
     REAL(DOUBLE)                    :: PRES
     REAL(DOUBLE)                    :: TTAU
     REAL(DOUBLE)                    :: PTAU
     LOGICAL                         :: REM_TRANS
     LOGICAL                         :: REM_ROTAT
     LOGICAL                         :: RESTRT
     LOGICAL                         :: AtomWrap
     LOGICAL                         :: CLOBBER
     CHARACTER(LEN=DCL)              :: RESTRT_IN
     CHARACTER(LEN=DCL)              :: RESTRT_OUT
     CHARACTER(LEN=DCL)              :: CRD_NAME
     CHARACTER(LEN=DCL)              :: VEL_NAME
     CHARACTER(LEN=DCL)              :: ENE_NAME
  END TYPE Dynamics

  TYPE Geometries
     INTEGER                         :: Klones
     TYPE(CRDS),POINTER,DIMENSION(:) :: Klone
  END TYPE Geometries

  TYPE State
     CHARACTER                       :: Action
     CHARACTER                       :: SubAction
     TYPE(INT_VECT)                  :: Current
     TYPE(INT_VECT)                  :: Previous
  END TYPE State

  TYPE BasisSets
     INTEGER                           :: NBSets
     INTEGER,  DIMENSION(MaxSets)      :: NExpt
     CHARACTER(LEN=BASESET_CHR_LEN),&
                DIMENSION(MaxSets)     :: BName  
     TYPE(BSET),POINTER,            &
                DIMENSION(:,:)         :: BSets
     TYPE(INT_VECT),POINTER,        &
                DIMENSION(:,:)         :: OffS,BSiz,LnDex
     TYPE(DBL_VECT),POINTER,        &
                DIMENSION(:,:)         :: DExpt
  END TYPE BasisSets

#ifdef PERIODIC
  TYPE Periodics
     INTEGER                           :: Dimen      !-- Dimension of the System
     LOGICAL                           :: AtomW      !-- Wrap atoms back into box--BE CAREFUL
     LOGICAL                           :: PFFOvRide  !-- Override of Automatic PFF stuff     
     LOGICAL                           :: InVecForm  !-- What form are the Lattice vectors in
     LOGICAL                           :: InAtomCrd  !-- Atomic or Fractional Coordinates
     LOGICAL                           :: Translate  !-- Should the Atomic Coordinated be Translated 
     LOGICAL                           :: Trans_COM  !-- Weither to Translate to The center of The Box
     LOGICAL,DIMENSION(3)              :: AutoW      !-- Periodic in X, Y and or Z  direction
     REAL(DOUBLE)                      :: Epsilon    !-- Epsilon at Infinity (Metal == Infinity)
  END TYPE Periodics
#endif

  TYPE Controls
     TYPE(FileNames)  :: Nams
     TYPE(Options)    :: Opts
     TYPE(Dynamics)   :: Dyns
     TYPE(Geometries) :: Geos
#ifdef PERIODIC
     TYPE(Periodics)  :: PBCs
#endif
     TYPE(BasisSets)  :: Sets
!     TYPE(MMechanics) :: QMMM
#ifdef PARALLEL
     TYPE(Parallel)   :: MPIs
#endif
     TYPE(State)      :: Stat
  END TYPE Controls

  INTEGER                     :: NLoc
  INTEGER, DIMENSION(MaxSets) :: Location
  CHARACTER(LEN=DCL)          :: Mssg

END MODULE ControlStructures
