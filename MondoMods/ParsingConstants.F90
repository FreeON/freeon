MODULE ParsingConstants
   USE GlobalScalars
   IMPLICIT NONE
!------------------------------------------------  
!  Parsing constants
   CHARACTER(LEN=26), PARAMETER   :: Upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   CHARACTER(LEN=26), PARAMETER   :: Lower='abcdefghijklmnopqrstuvwxyz'
   CHARACTER(LEN=23), PARAMETER   :: Numbers='0123456789-+.0123456789'
   CHARACTER(LEN=5),  PARAMETER   :: Special='$/_#*'
   CHARACTER(LEN=54), PARAMETER   :: Characters=Special//Lower//Numbers
   CHARACTER(LEN=9),  PARAMETER   :: Delimiters='[(|")],= '
   CHARACTER(LEN=4) , PARAMETER   :: Stars='****'
   CHARACTER(LEN=1) , PARAMETER   :: Space=' '
!  Parsing keys for <Options.Program=> or <Options.GLOBAL_DEBUG=>
   CHARACTER(LEN=11), PARAMETER :: GLOBAL_DEBUG ='DebugAll'
   CHARACTER(LEN=8),  PARAMETER :: DBG_NONE     ='NoDebug'
   CHARACTER(LEN=8),  PARAMETER :: DBG_MINIMUM  ='MinDebug'
   CHARACTER(LEN=8),  PARAMETER :: DBG_MEDIUM   ='MedDebug'
   CHARACTER(LEN=8),  PARAMETER :: DBG_MAXIMUM  ='MaxDebug'
   CHARACTER(LEN=8),  PARAMETER :: DBG_MATRICES ='MatDebug'
   CHARACTER(LEN=8),  PARAMETER :: DBG_PRT_SETS ='BasisSet'
   CHARACTER(LEN=8),  PARAMETER :: DBG_PRT_INTS ='PrintInt'
   CHARACTER(LEN=8),  PARAMETER :: DBG_PRT_RHO  ='PrintRho'
   CHARACTER(LEN=8),  PARAMETER :: PLT_MATRICES ='PlotMats'
   CHARACTER(LEN=8),  PARAMETER :: DBG_MMA_STYLE='MmaStyle'
   CHARACTER(LEN=8),  PARAMETER :: DBG_DBL_STYLE='DblStyle'
   CHARACTER(LEN=8),  PARAMETER :: DBG_FLT_STYLE='FltStyle'
!------------------------------------------------- 
   INTEGER,           PARAMETER :: DEBUG_NONE    =038108  ! Nothing
   INTEGER,           PARAMETER :: DEBUG_MINIMUM =138408  ! Timing and memory (Default)
   INTEGER,           PARAMETER :: DEBUG_MEDIUM  =285082  ! + check sums, sparsity, thresholds
   INTEGER,           PARAMETER :: DEBUG_MAXIMUM =420942  ! + intermediate values
   INTEGER,           PARAMETER :: DEBUG_MATRICES=382842  ! Print Matrices
   INTEGER,           PARAMETER :: DEBUG_DENSITY =390485  ! Print Density
   INTEGER,           PARAMETER :: DEBUG_BASISSET=409843  ! Print Basis sets
   INTEGER,           PARAMETER :: DEBUG_INTEGRAL=485583  ! Print Integrals
   INTEGER,           PARAMETER :: PLOT_MATRICES =608948  ! Plot Matrices
   INTEGER,           PARAMETER :: DEBUG_MMASTYLE=848423  ! Print in Mathematica style
   INTEGER,           PARAMETER :: DEBUG_FLTSTYLE=480484  ! Print in float style
   INTEGER,           PARAMETER :: DEBUG_DBLSTYLE=504843  ! Print in scientific (D) style
END MODULE
