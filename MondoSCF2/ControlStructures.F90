MODULE ControlStructures
  USE GlobalScalars
  USE DerivedTypes
  USE BasisSetParameters
  INTEGER, PARAMETER                :: MaxSets=6
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
     CHARACTER(LEN=DCL)             :: ReactantsFile!-- HDF5 with optimized reactants state
     CHARACTER(LEN=DCL)             :: ProductsFile !-- HDF5 with optimized products state
  END TYPE FileNames

  TYPE Options
     INTEGER                        :: NMthds
     INTEGER                        :: NModls
     INTEGER                        :: NThrsh
     INTEGER                        :: NSteps
     INTEGER                        :: Guess
     INTEGER                        :: Grad
     INTEGER                        :: EndPts
     INTEGER                        :: Coordinates
     INTEGER,   DIMENSION(MaxSets)  :: Methods
     INTEGER,   DIMENSION(MaxSets)  :: Models
     INTEGER,   DIMENSION(MaxSets)  :: AccuracyLevels
     LOGICAL                        :: DoGDIIS,SteepStep
     TYPE(INT_VECT)                 :: RestartState,ProductsState,ReactantsState
     TYPE(TOLS),DIMENSION(MaxSets)  :: Thresholds
     TYPE(DEBG)                     :: PFlags
     REAL(DOUBLE)                   :: NEBSpring
     LOGICAL                        :: NEBClimb      
     CHARACTER(LEN=3)               :: GeomPrint
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
     TYPE(CHR_VECT)                    :: Action ! Array of character options that tell a program what to do
     TYPE(INT_VECT)                    :: Current,Previous ! Tracking the SCF, basis and geometry states
  END TYPE State  
 
  TYPE GDIIS
    LOGICAL                            :: NoGDIIS    
    LOGICAL                            :: On
    INTEGER                            :: Init
    INTEGER                            :: MaxMem
  END TYPE GDIIS
  !
  TYPE Constr
     INTEGER                            :: NConstr
     INTEGER                            :: NCartConstr
     REAL(DOUBLE)                       :: ConstrMax
     REAL(DOUBLE)                       :: ConstrMaxCrit
     LOGICAL                            :: DoFixMM
     LOGICAL                            :: TSSearch
  END TYPE Constr
  !
  TYPE BackTrf
     INTEGER                            :: MaxIt_CooTrf
     REAL(DOUBLE)                       :: CooTrfCrit
     REAL(DOUBLE)                       :: RMSCrit
     REAL(DOUBLE)                       :: MaxCartDiff
     REAL(DOUBLE)                       :: DistRefresh
  END TYPE BackTrf
  !
  TYPE GrdTrf
     INTEGER                            :: MaxIt_GrdTrf
     REAL(DOUBLE)                       :: GrdTrfCrit
     REAL(DOUBLE)                       :: MaxGradDiff
  END TYPE GrdTrf
  !
  TYPE GConvCrit 
     REAL(DOUBLE)                       :: Grad
     REAL(DOUBLE)                       :: Stre
     REAL(DOUBLE)                       :: Bend
     REAL(DOUBLE)                       :: OutP
     REAL(DOUBLE)                       :: LinB
     REAL(DOUBLE)                       :: Tors
     INTEGER                            :: MaxGeOpSteps
     LOGICAL                            :: DoBackTr
  END TYPE GConvCrit 
  !
  TYPE GOptStat 
     INTEGER                            :: ActStep
     REAL(DOUBLE)                       :: MaxStreDispl
     REAL(DOUBLE)                       :: MaxBendDispl
     REAL(DOUBLE)                       :: MaxLinBDispl
     REAL(DOUBLE)                       :: MaxOutPDispl
     REAL(DOUBLE)                       :: MaxTorsDispl
     REAL(DOUBLE)                       :: MaxLDispl
     INTEGER                            :: IMaxGrad
     INTEGER                            :: IMaxCGrad
     INTEGER                            :: IMaxLGrad
     REAL(DOUBLE)                       :: MaxGrad
     REAL(DOUBLE)                       :: MaxCGrad
     REAL(DOUBLE)                       :: MaxLGrad
     REAL(DOUBLE)                       :: MaxDMult
     REAL(DOUBLE)                       :: RMSGrad
     INTEGER                            :: IMaxGradNoConstr
     REAL(DOUBLE)                       :: MaxGradNoConstr
     REAL(DOUBLE)                       :: RMSGradNoConstr
     LOGICAL                            :: GeOpConvgd
     REAL(DOUBLE)                       :: RMSIntDispl
  END TYPE GOptStat 
  !
  TYPE Hessian 
     REAL(DOUBLE)                       :: Stre
     REAL(DOUBLE)                       :: Bend
     REAL(DOUBLE)                       :: LinB
     REAL(DOUBLE)                       :: OutP
     REAL(DOUBLE)                       :: Tors
     REAL(DOUBLE)                       :: StpDescInvH
  END TYPE Hessian 
  !
  TYPE StepSize 
     REAL(DOUBLE)                       :: Stre
     REAL(DOUBLE)                       :: Bend
     REAL(DOUBLE)                       :: LinB
     REAL(DOUBLE)                       :: OutP
     REAL(DOUBLE)                       :: Tors
     REAL(DOUBLE)                       :: Cart
  END TYPE StepSize 
  !
  TYPE TrfCtrl
     LOGICAL                            :: DoFullTrf
     LOGICAL                            :: DoClssTrf
     LOGICAL                            :: DoInternals
     LOGICAL                            :: DoNewChol 
     LOGICAL                            :: DoRotOff
     LOGICAL                            :: DoTranslOff
     LOGICAL                            :: Linearity
     INTEGER,DIMENSION(3)               :: ThreeAt
     INTEGER,DIMENSION(3)               :: ThreeAt_2
     REAL(DOUBLE),DIMENSION(3)          :: TranslAt1
     REAL(DOUBLE),DIMENSION(3,3)        :: RotAt2ToX
     REAL(DOUBLE),DIMENSION(3,3)        :: RotAt3ToXY
     REAL(DOUBLE),DIMENSION(3)          :: TranslAt1_2
     REAL(DOUBLE),DIMENSION(3,3)        :: RotAt2ToX_2
     REAL(DOUBLE),DIMENSION(3,3)        :: RotAt3ToXY_2
  END TYPE TrfCtrl
  !
  TYPE CoordCtrl
     INTEGER                            :: RefreshIn
     INTEGER                            :: Refresh
     REAL(DOUBLE)                       :: VDWFact
     LOGICAL                            :: DoSelect
     CHARACTER(LEN=DEFAULT_CHR_LEN)     :: CoordType
     INTEGER                            :: NCov  
     INTEGER                            :: NExtra
     INTEGER                            :: NStre
     INTEGER                            :: NBend
     INTEGER                            :: NLinB
     INTEGER                            :: NOutP
     INTEGER                            :: NTors
     REAL(DOUBLE)                       :: LinCrit
     REAL(DOUBLE)                       :: TorsLinCrit
     REAL(DOUBLE)                       :: OutPCrit
  END TYPE CoordCtrl

  TYPE GeomOpt
     INTEGER                         :: Optimizer
     LOGICAL                         :: RestartBas
     TYPE(CoordCtrl)                 :: CoordCtrl
     TYPE(TrfCtrl)                   :: TrfCtrl
     TYPE(StepSize)                  :: StepSize
     TYPE(Hessian)                   :: Hessian
     TYPE(GOptStat)                  :: GOptStat
     TYPE(GConvCrit)                 :: GConvCrit
     TYPE(GrdTrf)                    :: GrdTrf
     TYPE(BackTrf)                   :: BackTrf
     TYPE(Constr)                    :: Constr
     TYPE(GDIIS)                     :: GDIIS
     TYPE(INTC)                      :: ExtIntCs
  END TYPE GeomOpt

  TYPE RespOpts
     LOGICAL               :: StcAlpha
     LOGICAL, DIMENSION(3) :: AlphaAxis
     LOGICAL               :: StcBeta
     LOGICAL               :: StcGamma
  END TYPE RespOpts

  TYPE PropOpts
     TYPE(RespOpts) :: Resp
  END TYPE PropOpts

  TYPE Controls
     TYPE(FileNames)  :: Nams
     TYPE(Options)    :: Opts
     TYPE(GeomOpt)    :: GOpt
     TYPE(Dynamics)   :: Dyns
     TYPE(Geometries) :: Geos
     TYPE(Periodics)  :: PBCs
     TYPE(BasisSets)  :: Sets
     TYPE(Parallel)   :: MPIs
     TYPE(State)      :: Stat
     TYPE(PropOpts)   :: POpt
  END TYPE Controls

  INTEGER                     :: NLoc
  INTEGER, DIMENSION(MaxSets) :: Location
  CHARACTER(LEN=DCL)          :: Mssg

END MODULE ControlStructures
