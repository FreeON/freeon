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
     INTEGER                        :: GradOpt
     INTEGER                        :: Coordinates
     INTEGER,   DIMENSION(MaxSets)  :: Methods
     INTEGER,   DIMENSION(MaxSets)  :: Models
     INTEGER,   DIMENSION(MaxSets)  :: AccuracyLevels
     LOGICAL                        :: OneBase,DoGDIIS
     TYPE(INT_VECT)                 :: RestartState
     TYPE(TOLS),DIMENSION(MaxSets)  :: Thresholds
     TYPE(DEBG)                     :: PFlags
  END TYPE Options

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
     INTEGER                         :: Images
     TYPE(INT_VECT)                  :: State
     TYPE(CRDS),POINTER,DIMENSION(:) :: Image
  END TYPE Geometries

#ifdef PERIODIC
  ! << CJ README >><< CJ README >><< CJ README >><< CJ README >><< CJ README >><< CJ README >>
  !
  ! >>>>>>>>>THESE SORT OF ITEMS DO NOT NEED TO BE ATTATCHED TO EACH AND EVERY GEOMETRY.<<<<<<  
  ! >>>>>>>>>THEY SHOULD BE DEFINED HERE ONCE AND FOR ALL FOR THE ENTIRE SIMMULATION. <<<<<<<< 
  !
  ! >>>>>>>>> ALSO, DOESNT LOOK LIKE BOX CARTS OR BOX VECT ARE USED FOR ANYTHING. <<<<<<<<<<<<
  ! >>>>>>>>> DO WE REALLY NEED TO KEEP THEM AROUND ? <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !
  !<< CJ README >><< CJ README >><< CJ README >><< CJ README >><< CJ README >><< CJ README >>


  !
  !<<CJ, HUGH AND GRAME README>> <<CJ, HUGH AND GRAME README>> <<CJ, HUGH AND GRAME README>> 
  !
  !    >>>>>>>>>DO WE WANT TO KEEP THE UNIT CELL CONSTANT FOR ALL IMAGES?<<<<<<<<<<<<<<<
  !
  !<<CJ, HUGH AND GRAME README>> <<CJ, HUGH AND GRAME README>> <<CJ, HUGH AND GRAME README>> 
  !
  TYPE Periodics
     INTEGER                         :: Dimen      !-- Dimension of the System
     LOGICAL                         :: AtomW      !-- Wrap atoms back into box--BE CAREFUL
     LOGICAL                         :: PFFOvRide  !-- Override of Automatic PFF stuff     
     LOGICAL                         :: InVecForm  !-- What form are the Lattice vectors in
     LOGICAL                         :: InAtomCrd  !-- Atomic or Fractional Coordinates
     LOGICAL                         :: Translate  !-- Should the Atomic Coordinated be Translated 
     LOGICAL                         :: Trans_COM  !-- Weither to Translate to The center of The Box
     LOGICAL,DIMENSION(3)            :: AutoW      !-- Periodic in X, Y and or Z  direction
     REAL(DOUBLE)                    :: Epsilon    !-- Epsilon at Infinity (Metal == Infinity)
  END TYPE Periodics
#endif

  TYPE BasisSets
     INTEGER                           ::NBSets
     CHARACTER(LEN=BASESET_CHR_LEN),&
                  DIMENSION(MaxSets)   :: BName  
     TYPE(INT_VECT),POINTER,        &
                DIMENSION(:,:)         :: OffS,BSiz
     TYPE(BSET),POINTER,            &
                DIMENSION(:,:)         :: BSets
  END TYPE BasisSets

  TYPE Controls
     TYPE(FileNames)  :: Nams
     TYPE(Options)    :: Opts
     TYPE(Dynamics)   :: Dyns
     TYPE(Geometries) :: Geos

     TYPE(Periodics)  :: PBCs
     TYPE(BasisSets)  :: Sets
!     TYPE(Parallel)   :: MPIs
!     TYPE(MMechanics) :: QMMM
  END TYPE Controls
  

  INTEGER                     :: NLoc
  INTEGER, DIMENSION(MaxSets) :: Location
  CHARACTER(LEN=DCL)          :: Mssg

END MODULE ControlStructures
