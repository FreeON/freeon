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
   INTEGER,      PARAMETER :: SPEll2=2*SPEll
!
   INTEGER,      PARAMETER :: HGLen=(HGEll+1)*(HGEll+2)*(HGEll+3)/6
   INTEGER,      PARAMETER :: SPLen=SPEll*(SPEll+3)/2                ! note, poles start from 0.
   INTEGER,      PARAMETER :: SPLen2=(2*SPEll)*((2*SPEll)+3)/2
!-------------------------------------------------------------------------------
! MD
!
    REAL(DOUBLE), PARAMETER :: KelvinToHartrees=3.166815208D-6  ! Boltz's Const. in Hartrees/Kelvin (NIST)
    REAL(DOUBLE), PARAMETER :: HartreesToKelvin=3.157746614D+5           
    REAL(DOUBLE), PARAMETER :: InternalTimeToSeconds=1.032749873D-15
    REAL(DOUBLE), PARAMETER :: SecondsToInternalTime=0.968288669D+15
    REAL(DOUBLE), PARAMETER :: BohrsToAngstroms=0.5291772083D0
!-------------------------------------------------------------------------------
!  Other numbers
!
   INTEGER,      PARAMETER :: BIG_INT      =2**28               ! bigest integer*4
   REAL(DOUBLE), PARAMETER :: Half         =One/Two             ! 1/2
   REAL(DOUBLE), PARAMETER :: ThreeHalves  =Three/Two           ! 3/2
   REAL(DOUBLE), PARAMETER :: FiveHalves   =Five/Two            ! 5/2
   REAL(DOUBLE), PARAMETER :: FiveFourths  =Five/Four           ! 5/4
   REAL(DOUBLE), PARAMETER :: AngstromsToAU=1.889725988578923D0 ! Angstronms -> AU
!#ifdef MMech
! data from NIST home page
   REAL(DOUBLE), PARAMETER :: JToHartree=2.29371276D17          ! Joul -> Hartree
   REAL(DOUBLE), PARAMETER :: C_Avogadro=6.02214199D23          ! Avogadro const 
   REAL(DOUBLE), PARAMETER :: e2PerAngstroemToKJPerMol=1389.3548461690 ! e2/A -> KJ/mol
   REAL(DOUBLE), PARAMETER :: KJPerMolPerAngstToHPerBohr=.00020155297074504836D0 ! KJ/mol/Angstroem -> Hartree/Bohr
   REAL(DOUBLE), PARAMETER :: HToJoule=4.3597482D-18 ! Hartree -> Joule 
!#endif
   REAL(DOUBLE), PARAMETER :: BIG_DBL      =HUGE(One)           ! bigest machine rep double
   REAL(DOUBLE), PARAMETER :: SMALL_DBL    =TINY(One)           ! smallest machine rep double
   REAL(DOUBLE), PARAMETER :: NuclearExpnt =1.D16               ! Exponent for nuclear delta 
   REAL(DOUBLE), PARAMETER :: LinCrit =1.D0  ! criterium for linearity of an angle in degree
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
   INTEGER, PARAMETER      :: INP_UNIT=115  
   INTEGER, PARAMETER      :: CRD_UNIT=110           
   INTEGER, PARAMETER      :: VEL_UNIT=111 
   INTEGER, PARAMETER      :: ENE_UNIT=112
   INTEGER, PARAMETER      :: RES_IN_UNIT=113
   INTEGER, PARAMETER      :: RES_OUT_UNIT=114
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
!-------------------------------------------------
!
! Van der Waals radii as listed in JMol program
! Only for 70 atoms, upto atomic number 83
! 0.0 values indicate elements for which this value was 
! not available, for these cases 2.1 A is recommended to 
! be used in bonding schemes.
! These data serve mainly for setting up the bonding scheme
! including VdW bonds.
!
   REAL(DOUBLE), DIMENSION(1:120) :: VDWRadii
   DATA VDWRadii(  1: 10) /1.2D0,  1.4D0,  1.82D0,1.372D0,0.795D0,1.7D0,  1.55D0, 1.52D0, 1.47D0, 1.54D0/
   DATA VDWRadii( 11: 20) /2.27D0, 1.73D0, 1.7D0, 2.1D0,  1.8D0,  1.8D0,  1.75D0, 1.88D0, 2.75D0, 2.45D0/
   DATA VDWRadii( 21: 30) /1.37D0, 1.37D0, 1.37D0,1.37D0, 1.37D0, 1.456D0,0.88D0, 0.69D0, 0.72D0, 0.74D0/
   DATA VDWRadii( 31: 40) /1.37D0, 1.95D0, 1.85D0,1.9D0,  1.85D0, 2.02D0, 1.58D0, 2.151D0,1.801D0,1.602D0/
   DATA VDWRadii( 41: 50) /1.468D0,1.526D0,1.36D0,1.339D0,1.345D0,1.376D0,1.27D0, 1.424D0,1.663D0,2.1D0/   
   DATA VDWRadii( 51: 60) /2.05D0, 2.06D0, 1.98D0,2.0D0,  1.84D0, 2.243D0,1.877D0,2.1D0,  2.1D0,  2.1D0/   
   DATA VDWRadii( 61: 70) /2.1D0,  2.1D0,  2.1D0, 2.1D0,  2.1D0,  2.1D0,  2.1D0,  2.1D0,  2.1D0,  2.1D0/   
   DATA VDWRadii( 71: 80) /2.17D0, 1.58D0,1.467D0,1.534D0,1.375D0,1.353D0,1.357D0,1.75D0,1.66D0, 1.55D0/  
   DATA VDWRadii( 81: 90) /1.96D0, 2.02D0, 2.15D0, 2.1D0,  2.1D0,  2.1D0,  2.1D0,  2.1D0,  2.1D0,  2.1D0/
   DATA VDWRadii( 91:100) /2.1D0,  2.1D0,  2.1D0, 2.1D0,  2.1D0,  2.1D0,  2.1D0,  2.1D0,  2.1D0,  2.1D0/   
   DATA VDWRadii(101:110) /2.1D0,  2.1D0,  2.1D0, 2.1D0,  2.1D0,  2.1D0,  2.1D0,  2.1D0,  2.1D0,  2.1D0/   
   DATA VDWRadii(111:120) /2.1D0,  2.1D0,  2.1D0, 2.1D0,  2.1D0,  2.1D0,  2.1D0,  2.1D0,  2.1D0,  2.1D0/   
!
!-------------------------------------------------
! Slater radii for bonding scheme recognition
! See J.C. Slater, JCP, V41,N10,p3199,y1964
! For nobel gas atoms and other missing elements
! the 0.75*VDWRadii value is recommended
!
   REAL(DOUBLE),DIMENSION(1:120) :: SLRadii
   DATA SLRadii(  1) /0.25D0/
   DATA SLRadii(  3: 9) /1.45D0, 1.05D0, 0.85D0, 0.70D0, 0.65D0, 0.60D0, 0.5D0 /
   DATA SLRadii( 11: 17) /1.8D0, 1.5D0, 1.25D0, 1.1D0, 1.0D0, 1.0D0, 1.0D0/
   DATA SLRadii( 19: 35) /2.2D0, 1.8D0, 1.6D0, 1.4D0, 1.35D0, 1.4D0, 1.4D0, &
                          1.4D0, 1.35D0, 1.35D0, 1.35D0, 1.35D0, 1.30D0, 1.25D0, &
                          1.15D0, 1.15D0, 1.15D0/
   DATA SLRadii( 37: 53) /2.35D0, 2.00D0, 1.80D0, 1.55D0, 1.45D0, 1.45D0, 1.35D0, &
                          1.30D0, 1.35D0, 1.40D0, 1.60D0, 1.55D0, 1.55D0, 1.45D0, & 
                          1.45D0, 1.40D0, 1.40D0/
   DATA SLRadii( 55: 84) /2.60D0, 2.15D0, 1.95D0,&
                          1.85D0, 1.85D0, 1.85D0, 1.85D0, 1.85D0, 1.85D0,&
                          1.80D0, &
                          1.75D0, 1.75D0, 1.75D0, 1.75D0, 1.75D0, 1.75D0, 1.75D0,&
                          1.55D0, 1.45D0, 1.35D0, 1.35D0, 1.30D0,&
                          1.35D0, 1.35D0, 1.35D0, 1.50D0, 1.90D0, 1.9D0, 1.60D0, 1.90D0 /
   DATA SLRadii( 88: 95) /2.15D0, 1.95D0, 1.80D0, 1.80D0, 1.75D0, 1.75D0, 1.75D0, 1.75D0/
!
   DATA SLRadii(2)   /1.05D0/ !!! He
   DATA SLRadii(10)  /1.155D0/ !!! Ne
   DATA SLRadii(18)  /1.410D0/ !!! Ar
   DATA SLRadii(36)  /1.515D0/ !!! Kr
   DATA SLRadii(54)  /1.500D0/ !!! Xe
   DATA SLRadii(86)  /1.550D0/ !!! Rn
!
!  DATA SLRadii(82)  /1.900D0/ !!! Pb set to Tl
   DATA SLRadii(85)  /1.550D0/ !!! At
   DATA SLRadii(87)  /1.550D0/ !!! Fr
!
   DATA SLRadii(96:120) /1.550D0,1.550D0,1.550D0,1.550D0,1.550D0,&
                         1.550D0,1.550D0,1.550D0,1.550D0,1.550D0,&
                         1.550D0,1.550D0,1.550D0,1.550D0,1.550D0,&
                         1.550D0,1.550D0,1.550D0,1.550D0,1.550D0,&
                         1.550D0,1.550D0,1.550D0,1.550D0,1.550D0/ 
!
! SCF global arrays 
!
   INTEGER,DIMENSION(3),SAVE :: Current
   INTEGER,DIMENSION(3),SAVE :: Previous
!
! Whether Periodic Boundary Condition is On
!
   LOGICAL :: PBC_On


#ifdef PARALLEL
!-------------------------------------------------  
!  MPI Scalars, default values
!
   INTEGER, SAVE           :: MyID=ROOT
   INTEGER, SAVE           :: MyKlone=1
   INTEGER, SAVE           :: NPrc=1
   LOGICAL, SAVE           :: InParallel=.FALSE.
   INTEGER, PARAMETER      :: MaxProc=1024
#endif
END MODULE GlobalScalars
