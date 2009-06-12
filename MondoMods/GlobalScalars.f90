!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------

!#define OLD_STYLE

#if defined OLD_STYLE
#warning using old constants (buggy)
#endif

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
  !  Whole numbers
  !
  REAL(DOUBLE), PARAMETER :: &
    Zero  = 0.0D0, &
    One   = 1.0D0, &
    Two   = 2.0D0, &
    Three = 3.0D0, &
    Four  = 4.0D0, &
    Five  = 5.0D0, &
    Six   = 6.0D0, &
    Seven = 7.0D0, &
    Eight = 8.0D0

  !-------------------------------------------------------------------------------
  !  Pi and related numbers
  !
#if defined OLD_STYLE
  REAL(DOUBLE), PARAMETER :: Pi         = 3.1415926535897932385D0  ! Pi
  REAL(DOUBLE), PARAMETER :: Pi2        = 9.8696044010893586188D0  ! Pi^2
  REAL(DOUBLE), PARAMETER :: Pi3        = 3.1006276680299820175D1  ! Pi^3
  REAL(DOUBLE), PARAMETER :: SqrtPi     = 1.7724538509055160273D0  ! Sqrt(Pi)
  REAL(DOUBLE), PARAMETER :: TwoPi      = 6.2831853071795864770D0  ! 2*Pi
  REAL(DOUBLE), PARAMETER :: TwoPi5x2   = 3.4986836655249725693D1  ! 2*Pi^(5/2)
  REAL(DOUBLE), PARAMETER :: Sqrt2Pi5x2 = 2.4739429451193148050D1  ! Sqrt(2)*Pi^(5/2)
  REAL(DOUBLE), PARAMETER :: DegToRad   = 1.7453292519943295769D-2 ! Pi/180
#else
  REAL(DOUBLE), PARAMETER :: Pi         = 3.1415926535897932385D0  ! Pi
  REAL(DOUBLE), PARAMETER :: Pi2        = Pi*Pi                    ! Pi^2
  REAL(DOUBLE), PARAMETER :: Pi3        = Pi*Pi*Pi                 ! Pi^3
  REAL(DOUBLE), PARAMETER :: SqrtPi     = 1.7724538509055160273D0  ! Sqrt(Pi)
  REAL(DOUBLE), PARAMETER :: TwoPi      = Two*Pi                   ! 2*Pi
  REAL(DOUBLE), PARAMETER :: TwoPi5x2   = 3.4986836655249725693D1  ! 2*Pi^(5/2)
  REAL(DOUBLE), PARAMETER :: Sqrt2Pi5x2 = 2.4739429451193148050D1  ! Sqrt(2)*Pi^(5/2)
  REAL(DOUBLE), PARAMETER :: DegToRad   = Pi/180.0D0               ! Pi/180
  REAL(DOUBLE), PARAMETER :: RadToDeg   = 1/DegToRad               ! 1/(Pi/180)
#endif
  !-------------------------------------------------------------------------------
  !  Max Angular Symmetries (d; 4+1=5, f;5+1=6, etc [+1 for gradients])
  !
  INTEGER,      PARAMETER :: BFEll=4          ! Max ell for a basis function+1; 4 == f+1 (+1 for gradients)
  INTEGER,      PARAMETER :: BFLen=(BFEll+1)*(BFEll+2)*(BFEll+3)/6
  INTEGER,      PARAMETER :: PrjEll=3         ! Max ell for ECP related projection opperators
  INTEGER,      PARAMETER :: ECPEll=2         ! Max radial exponent in Gaussian expansions of the ECP
  INTEGER,      PARAMETER :: HGEll=2*BFEll-1  ! Max ell for a distribution (bf product) and its derivative
  INTEGER,      PARAMETER :: SPEll=HGEll      !
  INTEGER,      PARAMETER :: SPEll2=2*SPEll   !
  INTEGER,      PARAMETER :: HGLen=(HGEll+1)*(HGEll+2)*(HGEll+3)/6
  INTEGER,      PARAMETER :: SPLen=SPEll*(SPEll+3)/2                ! note, poles start from 0.
  INTEGER,      PARAMETER :: SPLen2=(2*SPEll)*((2*SPEll)+3)/2

  !-------------------------------------------------------------------------------
  ! MD and others
  !
#if defined OLD_STYLE
  REAL(DOUBLE), PARAMETER :: KelvinToHartrees=3.166815208D-6    ! Boltz's Const. in Hartrees/Kelvin (NIST)
  REAL(DOUBLE), PARAMETER :: HartreesToKelvin=3.157746614D+5
  REAL(DOUBLE), PARAMETER :: InternalTimeToSeconds=1.032749873D-15
  REAL(DOUBLE), PARAMETER :: InternalTimeToFemtoseconds=1.032749873D0
  REAL(DOUBLE), PARAMETER :: SecondsToInternalTime=0.968288669D+15
  REAL(DOUBLE), PARAMETER :: FemtosecondsToInternalTime=0.968288669D0
  REAL(DOUBLE), PARAMETER :: BohrsToAngstroms=0.5291772083D0    ! AU ->  Angstronms
  REAL(DOUBLE), PARAMETER :: AngstromsToAU=1.889725988578923D0  ! Angstronms -> AU
  REAL(DOUBLE), PARAMETER :: GPaToAU=3.398928928849693861282D-5 ! GPa -> AU
#else
  REAL(DOUBLE), PARAMETER :: KelvinToHartrees=3.166815208D-6    ! Boltz's Const. in Hartrees/Kelvin (NIST)
  REAL(DOUBLE), PARAMETER :: HartreesToKelvin=One/KelvinToHartrees
  REAL(DOUBLE), PARAMETER :: InternalTimeToSeconds=1.032749873D-15
  REAL(DOUBLE), PARAMETER :: InternalTimeToFemtoseconds=1.032749873D0
  REAL(DOUBLE), PARAMETER :: SecondsToInternalTime=One/InternalTimeToSeconds
  REAL(DOUBLE), PARAMETER :: FemtosecondsToInternalTime=One/InternalTimeToFemtoseconds
  REAL(DOUBLE), PARAMETER :: AUToAngstroms=0.52917720859D0      ! AU ->  Angstroms
  REAL(DOUBLE), PARAMETER :: AngstromsToAU=One/AUToAngstroms    ! Angstroms -> AU
  REAL(DOUBLE), PARAMETER :: GPaToAU=3.398928928849693861282D-5 ! GPa -> AU
#endif
  ! These constants not yet verified:
  REAL(DOUBLE), PARAMETER :: au2eV=27.21139613182D0           ! au to eV
  REAL(DOUBLE), PARAMETER :: eV2au=1/au2eV                    ! eV to au
  REAL(DOUBLE), PARAMETER :: eVK=11676.33D0                   ! eV to K
  REAL(DOUBLE), PARAMETER :: eVcm=8065.73D0                   ! eV to cm-1
  REAL(DOUBLE), PARAMETER :: eVkkal=23.06035D0                ! eV to kkal
  REAL(DOUBLE), PARAMETER :: auDeb=4.803242                   ! exA to Debye
  REAL(DOUBLE), PARAMETER :: amuToKg = 1.66053886D-27         ! atomic mass units into kg.

  ! Symplectic expansion coefficient 4th order by McLachlan and Atela,
  ! Nonlinearity, vol 5, 541 (1992)
  REAL(DOUBLE), DIMENSION(4) :: Symplectic_4th_Order_a, Symplectic_4th_Order_b
  DATA Symplectic_4th_Order_a / &
     0.5153528374311229364D0, &
    -0.085782019412973646D0,  &
     0.4415830236164665242D0, &
     0.1288461583653841854D0 /
  DATA Symplectic_4th_Order_b / &
     0.1344961992774310892D0, &
    -0.2248198030794208058D0, &
     0.7563200005156682911D0, &
     0.3340036032863214255D0 /
  !-------------------------------------------------------------------------------
  !  Other numbers
  !
  INTEGER,      PARAMETER :: BIG_INT     = 2**28               ! bigest integer*4
  REAL(DOUBLE), PARAMETER :: Half        = One/Two             ! 1/2
  REAL(DOUBLE), PARAMETER :: OneThird    = One/Three           ! 1/3
  REAL(DOUBLE), PARAMETER :: TwoThird    = Two/Three           ! 2/3
  REAL(DOUBLE), PARAMETER :: ThreeHalves = Three/Two           ! 3/2
  REAL(DOUBLE), PARAMETER :: FiveHalves  = Five/Two            ! 5/2
  REAL(DOUBLE), PARAMETER :: FiveFourths = Five/Four           ! 5/4
  !#ifdef MMech
  ! data from NIST home page
  REAL(DOUBLE), PARAMETER :: JToHartree = 2.29371276D17          ! Joul -> Hartree
  REAL(DOUBLE), PARAMETER :: C_Avogadro = 6.02214199D23          ! Avogadro const
  REAL(DOUBLE), PARAMETER :: e2PerAngstroemToKJPerMol = 1389.3548461690 ! e2/A -> KJ/mol
  REAL(DOUBLE), PARAMETER :: KJPerMolPerAngstToHPerBohr = .00020155297074504836D0 ! KJ/mol/Angstroem -> Hartree/Bohr
  REAL(DOUBLE), PARAMETER :: HToJoule = 4.3597482D-18 ! Hartree -> Joule
  !#endif
  REAL(DOUBLE), PARAMETER :: BIG_DBL      = HUGE(One)           ! bigest machine rep double
  REAL(DOUBLE), PARAMETER :: SMALL_DBL    = TINY(One)           ! smallest machine rep double
  REAL(DOUBLE), PARAMETER :: NuclearExpnt = 1.D16               ! Exponent for nuclear delta

  !-------------------------------------------------------------------------------
  !  Status keys
  !
  INTEGER, PARAMETER      :: SUCCEED = 0
  INTEGER, PARAMETER      :: FAIL = -1

  !-------------------------------------------------
  !  Generic true false keys
  !
  INTEGER, PARAMETER      :: STATUS_TRUE  = 905788
  INTEGER, PARAMETER      :: STATUS_FALSE = 60672
  !-------------------------------------------------
  !  ROOT (MPI and Tree)
  !
  INTEGER, PARAMETER      :: ROOT=0
  !-------------------------------------------------  IO Unit numbers
  INTEGER, PARAMETER      :: Geo=22                 ! Unit for geometries
  INTEGER, PARAMETER      :: Bas=33                 ! Unit for basis file IO
  INTEGER, PARAMETER      :: GBas=34                ! Unit for Ghost basis file IO
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
  INTEGER, SAVE           :: NBasF,NEl,NAlph,NBeta,NAtoms
  REAL(DOUBLE), SAVE      :: TotCh
  INTEGER, SAVE           :: MaxAtmsNode,MaxBlksNode,MaxNon0Node
  INTEGER, SAVE           :: MaxAtms,MaxBlks,MaxNon0,MaxBlkSize
  INTEGER, PARAMETER      :: DIIS_MAX_MATRIX_SIZE=25
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
  ! Hydrogen modified from 0.25 to 0.30
  ! Oxigen modified from 0.6 to 0.61
  !
  REAL(DOUBLE),DIMENSION(1:120) :: SLRadii
  DATA SLRadii(  1) /0.30D0/
  DATA SLRadii(  3: 9) /1.45D0, 1.05D0, 0.85D0, 0.70D0, 0.65D0, 0.61D0, 0.5D0 /
  DATA SLRadii( 11: 17) /1.8D0, 1.50D0, 1.25D0, 1.1D0, 1.1D0, 1.0D0, 1.0D0/
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
  REAL(DOUBLE),DIMENSION(1:3)     :: Lindh_K
  DATA Lindh_K(1:3)   /0.45D0, 0.15D0, 0.005D0/
  REAL(DOUBLE),DIMENSION(1:3,1:3) :: Lindh_Alpha
  DATA Lindh_Alpha(1,1:3) /1.0000D0, 0.3949D0, 0.3949D0/
  DATA Lindh_Alpha(2,1:3) /0.3949D0, 0.2800D0, 0.2800D0/
  DATA Lindh_Alpha(3,1:3) /0.3949D0, 0.2800D0, 0.2800D0/
  !
  REAL(DOUBLE),DIMENSION(1:3,1:3) :: Lindh_R
  DATA Lindh_R(1,1:3) /1.3500D0, 2.1000D0, 2.5300D0/
  DATA Lindh_R(2,1:3) /2.1000D0, 2.8700D0, 3.4000D0/
  DATA Lindh_R(3,1:3) /2.5300D0, 3.4000D0, 3.4000D0/
  !these are the elements considered for H-bridges
  INTEGER,DIMENSION(10)  :: HBondList
  DATA HBondList(1:9) /5,7,8,9,16,34,17,35,53/
  ! Tabulated statistical Q-test values for 90% confidence
  REAL(DOUBLE),DIMENSION(10) :: QTest90
  DATA QTest90(1:10) /Zero,Zero,0.94D0,0.76D0,0.64D0,0.56D0,0.51D0,0.47D0,0.44D0,0.41D0/
  ! Diagonals of tabulated statistical F-test values for 95% confidence
  ! CAUTION! Values out of the 10x10 table may contain some inaccurate
  ! diagonals
  REAL(DOUBLE),DIMENSION(30) :: FTest95D
  DATA FTest95D( 1:10) /  Zero,19.00D0,9.28D0,6.39D0,5.05D0,4.28D0,3.79D0,3.44D0,3.18D0,2.98D0/
  DATA FTest95D(11:20) /2.79D0,2.69D0,2.53D0,2.46D0,2.40D0,2.28D0,2.23D0,2.19D0,2.16D0,2.12D0/
  DATA FTest95D(21:30) /2.10D0,2.07D0,2.05D0,2.02D0,1.99D0,1.96D0,1.93D0,1.90D0,1.87D0,1.84D0/
  !
  ! Maximum number of internal coords per atom
  ! this is needed to estimate storage of IntCs in HDF
  !
  INTEGER,PARAMETER :: IntCPerAtom=30
  !
  ! Maximum array size for remembered PBC lattice optimization points
  !
  INTEGER,PARAMETER  :: LattMaxMem=10
  !
  ! SCF global arrays
  !
  INTEGER,DIMENSION(3),SAVE :: Current
  INTEGER,DIMENSION(3),SAVE :: Previous
  !
  ! Whether Periodic Boundary Condition is On
  !
  LOGICAL :: PBC_On
  !
  !  MPI Scalars, default values
  !
  INTEGER, SAVE           :: MyID=ROOT
  INTEGER, SAVE           :: NPrc
  LOGICAL, SAVE           :: InParallel=.FALSE.
  INTEGER, PARAMETER      :: MaxProc=1024
  !
  !  Send in the clones...
  !
  INTEGER, SAVE           :: NClones
  INTEGER, SAVE           :: MyClone
END MODULE GlobalScalars
