MODULE ParsingKeys
   USE GlobalScalars
   USE GlobalCharacters
   IMPLICIT NONE
!------------------------------------------------  
!  General keys for parsing input
!
!  Title
   CHARACTER(LEN=12), PARAMETER :: BEGIN_TITLE   ='<BeginTitle>'
   CHARACTER(LEN=10), PARAMETER :: END_TITLE     ='<EndTitle>'
!  Generic options
   CHARACTER(LEN=14), PARAMETER :: BEGIN_OPTIONS ='<BeginOptions>'
   CHARACTER(LEN=12), PARAMETER :: END_OPTIONS   ='<EndOptions>'
#ifdef PERIODIC
!  Periodic
   CHARACTER(LEN=19), PARAMETER :: BEGIN_PERIODIC='<BeginPeriodic>'
   CHARACTER(LEN=17), PARAMETER :: END_PERIODIC  ='<EndPeriodic>'
#endif
!  Geometry
   CHARACTER(LEN=15), PARAMETER :: BEGIN_GEOMETRY                 ='<BeginGeometry>'
   CHARACTER(LEN=14), PARAMETER :: NEXT_GEOMETRY_MONDO_DEFAULT    ='<NextGeometry>'
   CHARACTER(LEN=13), PARAMETER :: END_GEOMETRY                   ='<EndGeometry>'
   CHARACTER(LEN=23), PARAMETER :: BEGIN_NEXT_GEOMETRY_MSI_ARCHIVE='input file for discover'
   CHARACTER(LEN=3),  PARAMETER :: END_NEXT_GEOMETRY_MSI_ARCHIVE  ='end'
!=====================================================================================
!  Parsing keys between <BeingOptions> and <EndOptions>:
!
!  <Options.Charge=>
   CHARACTER(LEN=6),  PARAMETER :: TOTAL_CHARGE ='Charge' 
!  <Options.Multiplicity=>
   CHARACTER(LEN=12),  PARAMETER :: MULTIPLICITY ='Multiplicity' 
!  <Options.SCF>=
   CHARACTER(LEN=9),  PARAMETER :: SCF_OPTION='SCFMethod' 
   CHARACTER(LEN=4),  PARAMETER :: SCF_SDMM  ='SDMM' 
   CHARACTER(LEN=2),  PARAMETER :: SCF_RHHF  ='RH'
   CHARACTER(LEN=6),  PARAMETER :: SCF_INKF  ='InkFok'
   CHARACTER(LEN=4),  PARAMETER :: SCF_DIIS  ='DIIS'
!  SCF keys
   INTEGER, PARAMETER :: RH_R_SCF  =10982348 ! Restricted Roothaan-Hall 
   INTEGER, PARAMETER :: RH_U_SCF  =43240823 ! Unrestricted Roothaan-Hall 
   INTEGER, PARAMETER :: SDMM_R_SCF=58638502 ! Restricted Simplified Density Matrix Minimization 
   INTEGER, PARAMETER :: SDMM_U_SCF=92874018 ! Unrestricted Simplified Density Matrix Minimization 
   INTEGER, PARAMETER :: EDMM_R_SCF=48530687 ! Restricted Extrapolated Density Matrix Minimization 
   INTEGER, PARAMETER :: EDMM_U_SCF=40834070 ! Unrestricted ExtrapolatedDensity Matrix Minimization 
!  <Options.Guess=>
   CHARACTER(LEN=5),  PARAMETER :: GUESS_OPTION   ='Guess'
   CHARACTER(LEN=4),  PARAMETER :: GUESS_CORE     ='Core'
   CHARACTER(LEN=8),  PARAMETER :: GUESS_SUPER    ='SuperPos'
   CHARACTER(LEN=7),  PARAMETER :: GUESS_RESTART  ='Restart'
   CHARACTER(LEN=7),  PARAMETER :: RESTART_INFO   ='InfoFile'
!  Guess keys
   INTEGER, PARAMETER :: GUESS_EQ_CORE=14334 ! Density matrix from core Hamiltonian
   INTEGER, PARAMETER :: GUESS_EQ_SUPR=40823 ! Density matrix from superposition of 
!  <Options.Geometry>=
   CHARACTER(LEN=8),  PARAMETER :: GEOMETRY     ='Geometry' 
   CHARACTER(LEN=8),  PARAMETER :: IN_AU        ='InAu'
   CHARACTER(LEN=6),  PARAMETER :: H_ORDER      ='HOrder'
   CHARACTER(LEN=6),  PARAMETER :: Z_ORDER      ='ZOrder'
   CHARACTER(LEN=11), PARAMETER :: RANDOM_ORDER ='RandomOrder'
   CHARACTER(LEN=7),  PARAMETER :: NO_ORDER     ='NoOrder'
   CHARACTER(LEN=9),  PARAMETER :: MSI_FORMAT   ='MSIFormat'
   CHARACTER(LEN=10), PARAMETER :: XMOL_FORMAT  ='XMolFormat'
!--------------------------------------------------------
!  Options involving force evaluations
!
   INTEGER, PARAMETER           :: GRAD_NO_GRAD   = 1000001 ! Do no gradeints evaluation
!  <Options.Grad>=
   CHARACTER(LEN=4),  PARAMETER :: GRADIENTS      ='Grad'
   CHARACTER(LEN=5 ), PARAMETER :: FORCE          ='Force'
   INTEGER, PARAMETER           :: GRAD_ONE_FORCE = 1084814 ! Perform one force evaluation
!  <Options.Opt=>
   CHARACTER(LEN=3),  PARAMETER :: OPTIMIZATION   ='Opt'
   CHARACTER(LEN=5),  PARAMETER :: OPT_QUNEW      ='QuNew'
   CHARACTER(LEN=2),  PARAMETER :: OPT_TSTATE     ='TS'
   INTEGER, PARAMETER           :: GRAD_QNEW_OPT  = 3489343 ! Perform quasi-newton geometry optimization
   INTEGER, PARAMETER           :: GRAD_TS_SEARCH = 3577711 ! Perform a gradients only transition state search.
!  <Options.MD=>
   CHARACTER(LEN=2),  PARAMETER :: DYNAMICS       ='MD'
   CHARACTER(LEN=6),  PARAMETER :: MD_VERLET      ='Verlet'
   CHARACTER(LEN=5),  PARAMETER :: MAX_STEPS      ='Steps'
   CHARACTER(LEN=8),  PARAMETER :: MD_TIME_STEP   ='TimeStep' 
   INTEGER, PARAMETER           :: GRAD_VERLET_MD = 6413123 ! Perform verlet dynamics
#ifdef PERIODIC

!  Parsing keys for <Options.Periodic=>
   CHARACTER(LEN=9),  PARAMETER :: PBOUNDRY     ='Periodic' 
!
   CHARACTER(LEN=3),  PARAMETER :: PBC_OFF      ='Off'
   CHARACTER(LEN=4),  PARAMETER :: PBC_X        ='On-X'
   CHARACTER(LEN=4),  PARAMETER :: PBC_Y        ='On-Y'
   CHARACTER(LEN=4),  PARAMETER :: PBC_Z        ='On-Z'
   CHARACTER(LEN=5),  PARAMETER :: PBC_XY       ='On-XY'
   CHARACTER(LEN=5),  PARAMETER :: PBC_XZ       ='On-XZ'
   CHARACTER(LEN=5),  PARAMETER :: PBC_YZ       ='On-YZ'
   CHARACTER(LEN=6),  PARAMETER :: PBC_XYZ      ='On-XYZ'
!
   CHARACTER(LEN=8),  PARAMETER :: ATOMW_ON     ='AtomWrap'
   CHARACTER(LEN=10), PARAMETER :: ATOMW_OFF    ='NoAtomWrap'
   CHARACTER(LEN=9),  PARAMETER :: LVF_VEC      ='VecFormat'
   CHARACTER(LEN=9),  PARAMETER :: LVF_ANG      ='AngFormat'
!
   CHARACTER(LEN=2),  PARAMETER :: TRAN_VEC     ='tv'
   CHARACTER(LEN=2),  PARAMETER :: ALAT_VEC     ='av'
   CHARACTER(LEN=2),  PARAMETER :: BLAT_VEC     ='bv'
   CHARACTER(LEN=2),  PARAMETER :: CLAT_VEC     ='cv'
!
   CHARACTER(LEN=10), PARAMETER :: CRT_ATOM     ='AtomCoord'
   CHARACTER(LEN=10), PARAMETER :: CRT_FRAC     ='FracCoord'
#endif
!-------------------------------------------------  
!  Parsing keys for <Options.BasisSets=>
!
   CHARACTER(LEN=9),  PARAMETER :: BASIS_SETS    ='BasisSets' 
!---------------------------------------------------------
!  Parsing keys for <Options.Accuracy=>
!
   CHARACTER(LEN=8),  PARAMETER :: ACCURACY_OPTION   ='Accuracy' 
   CHARACTER(LEN=6),  PARAMETER :: ACCURACY_CHEEZY   ='Loose'
   CHARACTER(LEN=4),  PARAMETER :: ACCURACY_GOOD     ='Good'
   CHARACTER(LEN=5),  PARAMETER :: ACCURACY_TIGHT    ='Tight'
   CHARACTER(LEN=9),  PARAMETER :: ACCURACY_RETENTIVE='VeryTight'
#ifdef PARALLEL
!-------------------------------------------------  
!  Parsing keys for <Options.MPIRun=>
!
   CHARACTER(LEN=6),  PARAMETER :: MPI_OPTION    ='MPIRun'
#endif 
END MODULE
