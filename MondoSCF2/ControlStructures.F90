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
     CHARACTER(LEN=DCL)             :: M_EXEC
     CHARACTER(LEN=DCL)             :: M_SCRATCH
     CHARACTER(LEN=DCL)             :: SCF_NAME
     CHARACTER(LEN=DCL)             :: IFile        !-- Input file
     CHARACTER(LEN=DCL)             :: OFile        !-- Output file
     CHARACTER(LEN=DCL)             :: LFile        !-- Log file
     CHARACTER(LEN=DCL)             :: GFile        !-- Geometry file
     CHARACTER(LEN=DCL)             :: HFile        !-- HDF5 file
     CHARACTER(LEN=DCL)             :: RFile        !-- Restart HDF5 file
  END TYPE FileNames

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
     LOGICAL                        :: DoGDIIS,SteepStep
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
     INTEGER                         :: Clones
     TYPE(CRDS),POINTER,DIMENSION(:) :: Clone
  END TYPE Geometries

  TYPE BasisSets
     INTEGER                           :: NBSets
     INTEGER,  DIMENSION(MaxSets)      :: NExpt
     INTEGER,  DIMENSION(MaxSets)      :: MxAts
     INTEGER,  DIMENSION(MaxSets)      :: MxN0s
     INTEGER,  DIMENSION(MaxSets)      :: MxBlk
     CHARACTER(LEN=BASESET_CHR_LEN),&
          DIMENSION(MaxSets)     :: BName  
     TYPE(BSET),POINTER,            &
          DIMENSION(:,:)         :: BSets
     TYPE(INT_VECT),POINTER,        &
          DIMENSION(:,:)         :: OffS,BSiz,LnDex
     TYPE(DBL_VECT),POINTER,        &
          DIMENSION(:,:)         :: DExpt
     REAL(DOUBLE),POINTER,        &
          DIMENSION(:,:)         :: AtomPairThresh,PrimPairThresh 
  END TYPE BasisSets

#ifdef PERIODIC
  TYPE Periodics
     INTEGER                           :: Dimen      !-- Dimension of the System
     INTEGER                           :: PFFMAXLAY
     INTEGER                           :: PFFMAXELL  
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

  TYPE Parallel
#ifdef PARALLEL
     CHARACTER(LEN=DCL)                :: Invoking
     CHARACTER(LEN=DCL)                :: ProcFlag
     CHARACTER(LEN=DCL)                :: MachFlag
     CHARACTER(LEN=DCL)                :: MachFile
     INTEGER                           :: NProc
     INTEGER                           :: NSpace
     INTEGER, DIMENSION(MaxSets)       :: MxAtsNode,MxBlkNode,MxN0sNode
     TYPE(INT_VECT),POINTER,        &
          DIMENSION(:,:)         :: Beg,End,GLO
#endif     
     INTEGER                           :: Clumps
     TYPE(INT_RNK2)                    :: Clump
  END TYPE Parallel

  TYPE State
     CHARACTER(LEN=DCL)              :: Action
     CHARACTER(LEN=DCL)              :: SubAction
     TYPE(INT_VECT)                  :: Current,Previous
  END TYPE State

  TYPE OptimizerType
     INTEGER                            :: AccL
     CHARACTER(LEN=DEFAULT_CHR_LEN)     :: CoordType
     LOGICAL                            :: DoInternals
     LOGICAL                            :: DoRotOff
     LOGICAL                            :: DoTranslOff
     INTEGER                            :: ActStep   
     INTEGER                            :: ReDefIntC
     INTEGER                            :: MaxGeOpSteps
     INTEGER                            :: BlkGeomSize 
     INTEGER                            :: IMaxGrad
     INTEGER                            :: IMaxGradNoConstr
     INTEGER                            :: MaxIt_GrdTrf
     INTEGER                            :: MaxIt_CooTrf
     INTEGER                            :: NConstr 
     INTEGER                            :: NCartConstr 
     INTEGER                            :: GDIISMinDomCount
     INTEGER                            :: GDIISInit
     INTEGER                            :: GDIISMaxMem
     INTEGER                            :: LSStepMax
     REAL(DOUBLE)                       :: AINVThrsh    
     REAL(DOUBLE)                       :: BMatThrsh   
     REAL(DOUBLE)                       :: StreHessian  
     REAL(DOUBLE)                       :: BendHessian 
     REAL(DOUBLE)                       :: LinBHessian 
     REAL(DOUBLE)                       :: OutPHessian 
     REAL(DOUBLE)                       :: TorsHessian 
     REAL(DOUBLE)                       :: StpDescInvH 
     REAL(DOUBLE)                       :: MaxGrad      
     REAL(DOUBLE)                       :: RMSGrad      
     REAL(DOUBLE)                       :: OldRMSGrad      
     REAL(DOUBLE)                       :: MaxGradNoConstr
     REAL(DOUBLE)                       :: RMSGradNoConstr
     REAL(DOUBLE)                       :: MaxStreDispl 
     REAL(DOUBLE)                       :: MaxBendDispl 
     REAL(DOUBLE)                       :: MaxLinBDispl 
     REAL(DOUBLE)                       :: MaxOutPDispl 
     REAL(DOUBLE)                       :: MaxTorsDispl 
     REAL(DOUBLE)                       :: RMSIntDispl   
     REAL(DOUBLE)                       :: GradCrit     
     REAL(DOUBLE)                       :: StreConvCrit
     REAL(DOUBLE)                       :: BendConvCrit
     REAL(DOUBLE)                       :: OutPConvCrit
     REAL(DOUBLE)                       :: LinBConvCrit
     REAL(DOUBLE)                       :: TorsConvCrit
     REAL(DOUBLE)                       :: AINVThresh    
     REAL(DOUBLE)                       :: GrdTrfCrit    
     REAL(DOUBLE)                       :: MaxGradDiff   
     REAL(DOUBLE)                       :: CooTrfCrit    
     REAL(DOUBLE)                       :: RMSCrit    
     REAL(DOUBLE)                       :: MaxCartDiff   
     REAL(DOUBLE)                       :: DistRefresh   
     REAL(DOUBLE)                       :: ConstrMax     
     REAL(DOUBLE)                       :: ConstrMaxCrit 
     REAL(DOUBLE)                       :: GDIISBandWidth 
     REAL(DOUBLE)                       :: GDIISMetric
     LOGICAL                            :: NoGDIIS    
     LOGICAL                            :: GeOpConvgd    
     LOGICAL                            :: GDIISOn
     LOGICAL                            :: GDIISMetricOn
     CHARACTER(LEN=DEFAULT_CHR_LEN)     :: GDIISCoordType
  END TYPE OptimizerType

  TYPE Controls
     TYPE(FileNames)  :: Nams
     TYPE(Options)    :: Opts
     TYPE(OptimizerType)  :: Mizr
     TYPE(Dynamics)   :: Dyns
     TYPE(Geometries) :: Geos
#ifdef PERIODIC
     TYPE(Periodics)  :: PBCs
#endif
     TYPE(BasisSets)  :: Sets
     TYPE(Parallel)   :: MPIs
     TYPE(State)      :: Stat
  END TYPE Controls

  INTEGER                     :: NLoc
  INTEGER, DIMENSION(MaxSets) :: Location
  CHARACTER(LEN=DCL)          :: Mssg

END MODULE ControlStructures
