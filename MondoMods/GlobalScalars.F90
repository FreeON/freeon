MODULE GlobalScalars
   IMPLICIT NONE
!-------------------------------------------------------------------------------
!  Default system types
!
   INTEGER, PARAMETER :: INTEGR=KIND(0)              !--Integer*4 or Integer*8
   INTEGER, PARAMETER :: SINGLE=KIND(0.0)            !--Real*4
   INTEGER, PARAMETER :: DOUBLE=KIND(0.D0)           !--Real*8
!-------------------------------------------------------------------------------
!  Defined precision types
!
   INTEGER, PARAMETER :: INT1=SELECTED_INT_KIND(2)  !--Integer*1
   INTEGER, PARAMETER :: INT2=SELECTED_INT_KIND(4)  !--Integer*2
   INTEGER, PARAMETER :: INT4=SELECTED_INT_KIND(9)  !--Integer*4
   INTEGER, PARAMETER :: INT8=SELECTED_INT_KIND(18) !--Integer*8
!-------------------------------------------------------------------------------
!  Pi and related numbers
!
   REAL(DOUBLE), PARAMETER :: Pi        =3.1415926535897932385D0  ! Pi
   REAL(DOUBLE), PARAMETER :: Pi2       =9.8696044010893586188D0  ! Pi^2
   REAL(DOUBLE), PARAMETER :: Pi3       =3.1006276680299820175D1  ! Pi^3
   REAL(DOUBLE), PARAMETER :: SqrtPi    =1.7724538509055160273D0  ! Sqrt(Pi)
   REAL(DOUBLE), PARAMETER :: TwoPi5x2  =3.4986836655249725693D1  ! 2*Pi^(5/2)
   REAL(DOUBLE), PARAMETER :: Sqrt2Pi5x2=2.4739429451193148050D1  ! Sqrt(2)*Pi^(5/2)
   REAL(DOUBLE), PARAMETER :: DegToRad  =1.7453292519943295769D-2 ! Pi/180 
!-------------------------------------------------------------------------------
!  Whole numbers
!
   REAL(DOUBLE), PARAMETER :: Zero=0.0D0, One  =1.0D0, &
                              Two =2.0D0, Three=3.0D0, &
                              Four=4.0D0, Five =5.0D0, &
                              Six =6.0D0
!-------------------------------------------------------------------------------
!  Max Angular Symmetries (d; 4+1=5, f;5+1=6, etc [+1 for gradients])
!   
   INTEGER,      PARAMETER :: HGEll=5
   INTEGER,      PARAMETER :: SPEll=5
   INTEGER,      PARAMETER :: FFEll=16
!
   INTEGER,      PARAMETER :: SPEll2=2*SPEll
   INTEGER,      PARAMETER :: FFEll2=2*FFEll
!
   INTEGER,      PARAMETER :: HGLen=(HGEll+1)*(HGEll+2)*(HGEll+3)/6
   INTEGER,      PARAMETER :: SPLen=SPEll*(SPEll+3)/2                ! note, poles start from 0.
   INTEGER,      PARAMETER :: SPLen2=(2*SPEll)*((2*SPEll)+3)/2
   INTEGER,      PARAMETER :: FFLen=FFEll*(FFEll+3)/2                ! note, poles start from 0.
   INTEGER,      PARAMETER :: FFLen2=(2*FFEll)*((2*FFEll)+3)/2
!-------------------------------------------------------------------------------
!  Other numbers
!
   INTEGER,      PARAMETER :: BIG_INT      =2**28               ! bigest integer*4
   REAL(DOUBLE), PARAMETER :: Half         =One/Two             ! 1/2
   REAL(DOUBLE), PARAMETER :: ThreeHalves  =Three/Two           ! 3/2
   REAL(DOUBLE), PARAMETER :: FiveHalves   =Five/Two            ! 5/2
   REAL(DOUBLE), PARAMETER :: FiveFourths  =Five/Four           ! 5/4
   REAL(DOUBLE), PARAMETER :: AngstromsToAU=1.889725988578923D0 ! Angstronms -> AU
   REAL(DOUBLE), PARAMETER :: BIG_DBL      =HUGE(One)           ! bigest machine rep double
   REAL(DOUBLE), PARAMETER :: NuclearExpnt =1.D16               ! Exponent for nuclear delta 
!-------------------------------------------------------------------------------
!  Status keys
!
   INTEGER, PARAMETER      :: SUCCEED=0
   INTEGER, PARAMETER      :: FAIL=-1
!-------------------------------------------------  
!  Generic true false keys
!
   INTEGER, PARAMETER      :: STATUS_TRUE =90575
   INTEGER, PARAMETER      :: STATUS_FALSE=60672
!-------------------------------------------------  
!  ROOT (MPI and Tree)
!
   INTEGER, PARAMETER      :: ROOT=0
!-------------------------------------------------  IO Unit numbers 
   INTEGER, PARAMETER      :: Geo=22                 ! Unit for geometries
   INTEGER, PARAMETER      :: Bas=33                 ! Unit for basis file IO
   INTEGER, PARAMETER      :: Seq=44                 ! Unit for sequential, binary IO
   INTEGER, PARAMETER      :: Tmp=55                 ! Unit for temp/scratch files
   INTEGER, PARAMETER      :: Plt=66                 ! Unit for EPS plotting
   INTEGER, PARAMETER      :: Inp=77                 ! Unit for ASCI input files
   INTEGER, PARAMETER      :: Out=88                 ! Unit for ASCI output files
   INTEGER, PARAMETER      :: LgF=99                 ! Unit for ASCI log files
!-------------------------------------------------  Matrix dimensions
   INTEGER, SAVE           :: NBasF,NEl,NAtoms
#ifdef PARALLEL
   INTEGER, SAVE           :: MaxAtmsNode,MaxBlksNode,MaxNon0Node               
#endif
   INTEGER, SAVE           :: MaxAtms,MaxBlks,MaxNon0,MaxBlkSize
!-------------------------------------------------  
!  Memory managment
! 
   INTEGER, SAVE           :: SizeOfInt,SizeOfDbl
   REAL(DOUBLE), SAVE      :: IntToMB,DblToMB
#ifdef PARALLEL
!-------------------------------------------------  
!  MPI Scalars, default values
!
   INTEGER, SAVE           :: MyID=ROOT
   INTEGER, SAVE           :: NPrc=1
   LOGICAL, SAVE           :: InParallel=.FALSE.
   INTEGER, PARAMETER      :: MaxProc=1024
#endif
END MODULE GlobalScalars
