MODULE ParsingConstants
   USE GlobalScalars
   IMPLICIT NONE
!------------------------------------------------  
!  Parsing constants
   CHARACTER(LEN=*), PARAMETER   :: Upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   CHARACTER(LEN=*), PARAMETER   :: Lower='abcdefghijklmnopqrstuvwxyz'
   CHARACTER(LEN=*), PARAMETER   :: Numbers='0123456789-+.0123456789'
   CHARACTER(LEN=*),  PARAMETER   :: Special='$/_#*='
   CHARACTER(LEN=*), PARAMETER   :: Characters=Special//Lower//Numbers
   CHARACTER(LEN=*),  PARAMETER   :: Delimiters='[(|")],= '
   CHARACTER(LEN=*) , PARAMETER   :: Stars='****'
   CHARACTER(LEN=*) , PARAMETER   :: Space=' '
!  Parsing keys for <Options.Program=> or <Options.GLOBAL_DEBUG=>
   CHARACTER(LEN=*), PARAMETER :: GLOBAL_DEBUG ='DebugAll'
   CHARACTER(LEN=*),  PARAMETER :: DBG_NONE     ='NoDebug'
   CHARACTER(LEN=*),  PARAMETER :: DBG_MINIMUM  ='MinDebug'
   CHARACTER(LEN=*),  PARAMETER :: DBG_MEDIUM   ='MedDebug'
   CHARACTER(LEN=*),  PARAMETER :: DBG_MAXIMUM  ='MaxDebug'
   CHARACTER(LEN=*),  PARAMETER :: DBG_MATRICES ='MatDebug'
   CHARACTER(LEN=*),  PARAMETER :: DBG_CHKSUMS  ='CheckSums'
   CHARACTER(LEN=*),  PARAMETER :: DBG_PRT_SETS ='BasisSet'
   CHARACTER(LEN=*),  PARAMETER :: DBG_PRT_INTS ='PrintInt'
   CHARACTER(LEN=*),  PARAMETER :: DBG_PRT_RHO  ='PrintRho'
   CHARACTER(LEN=*),  PARAMETER :: PLT_MATRICES ='PlotMats'
   CHARACTER(LEN=*),  PARAMETER :: DBG_MMA_STYLE='MmaStyle'
   CHARACTER(LEN=*),  PARAMETER :: DBG_DBL_STYLE='DblStyle'
   CHARACTER(LEN=*),  PARAMETER :: DBG_FLT_STYLE='FltStyle'
   CHARACTER(LEN=*),  PARAMETER :: DBG_GEOP_MIN ='MinGeOp'
   CHARACTER(LEN=*),  PARAMETER :: DBG_GEOP_MAX ='MaxGeOp'
   CHARACTER(LEN=*),  PARAMETER :: DBG_PRT_MM   ='MMDebug'

   INTEGER, PARAMETER     :: SFC_NONE   =3480481
   INTEGER, PARAMETER     :: SFC_PEANO  =5308208
   INTEGER, PARAMETER     :: SFC_HILBERT=4808471
   INTEGER, PARAMETER     :: SFC_RANDOM =5058108
   INTEGER, PARAMETER     :: SFC_TRAVEL =2505811
   INTEGER, PARAMETER     :: SFC_TableTrav =2505812
!------------------------------------------------- 
   INTEGER,           PARAMETER :: DEBUG_NONE    =038108  ! Nothing
   INTEGER,           PARAMETER :: DEBUG_MINIMUM =138408  ! Timing and memory (Default)
   INTEGER,           PARAMETER :: DEBUG_MEDIUM  =285082  ! + sparsity, thresholds
   INTEGER,           PARAMETER :: DEBUG_CHKSUMS =414234  ! + check sums
   INTEGER,           PARAMETER :: DEBUG_MAXIMUM =420942  ! + intermediate values
   INTEGER,           PARAMETER :: DEBUG_MATRICES=482842  ! Print Matrices
   INTEGER,           PARAMETER :: DEBUG_DENSITY =490485  ! Print Density
   INTEGER,           PARAMETER :: DEBUG_BASISSET=509843  ! Print Basis sets
   INTEGER,           PARAMETER :: DEBUG_INTEGRAL=585583  ! Print Integrals
   INTEGER,           PARAMETER :: PLOT_MATRICES =608948  ! Plot Matrices
   INTEGER,           PARAMETER :: DEBUG_MMASTYLE=848423  ! Print in Mathematica style
   INTEGER,           PARAMETER :: DEBUG_FLTSTYLE=480484  ! Print in float style
   INTEGER,           PARAMETER :: DEBUG_DBLSTYLE=504843  ! Print in scientific (D) style
   INTEGER,           PARAMETER :: DEBUG_GEOP_MIN=568356  ! Print geometry optimization data 
   INTEGER,           PARAMETER :: DEBUG_GEOP_MAX=568357  ! Print geometry optimization data 
   INTEGER,           PARAMETER :: DEBUG_MM      =129462  ! Print molecular mechanics related data
END MODULE
