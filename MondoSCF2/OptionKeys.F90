MODULE OptionKeys
  IMPLICIT NONE
  !  Generic options
  CHARACTER(LEN=*), PARAMETER :: OPTIONS_BEGIN      ='<BeginOptions>'
  CHARACTER(LEN=*), PARAMETER :: OPTIONS_END        ='<EndOptions>'
  ! Options:SCF=
  CHARACTER(LEN=*),  PARAMETER :: SCF_OPTION         ='SCFMethod' 
  ! Restricted Roothaan-Hall 
  CHARACTER(LEN=*),  PARAMETER :: SCF_RHHF           ='RH'
  INTEGER,           PARAMETER :: RH_R_SCF           =10982348 
  ! Restricted Simplified Density Matrix Minimization 
  CHARACTER(LEN=*),  PARAMETER :: SCF_SDMM           ='SDMM' 
  INTEGER,           PARAMETER :: SDMM_R_SCF         =58638502 
  ! Restricted Palser-Manulopolus Purification 
  CHARACTER(LEN=*),  PARAMETER :: SCF_PM             ='PM' 
  INTEGER,           PARAMETER :: PM_R_SCF           =58638503 
  ! Restricted quadratic Trace Correcting Purification 
  CHARACTER(LEN=*),  PARAMETER :: SCF_SP2            ='TC2'
  INTEGER,           PARAMETER :: SP2_R_SCF          =58638504 
  ! Restricted quartic Trace Correcting Purification
  CHARACTER(LEN=*),  PARAMETER :: SCF_SP4            ='TC4'
  INTEGER,           PARAMETER :: SP4_R_SCF          =58638505 
  ! Restricted quartic Trace Re-Setting Purification
  CHARACTER(LEN=*),  PARAMETER :: SCF_TS4            ='TRS4'
  INTEGER,           PARAMETER :: TS4_R_SCF          =58638506 
  !---------------------------------------------------------
  ! Options:Accuracy=
  CHARACTER(LEN=*),  PARAMETER :: ACCURACY_OPTION    ='Accuracy' 
  CHARACTER(LEN=*),  PARAMETER :: ACCURACY_CHEEZY    ='Loose'
  CHARACTER(LEN=*),  PARAMETER :: ACCURACY_GOOD      ='Good'
  CHARACTER(LEN=*),  PARAMETER :: ACCURACY_TIGHT     ='Tight'
  CHARACTER(LEN=*),  PARAMETER :: ACCURACY_RETENTIVE ='VeryTight'
  !------------------------------------------------------------------------------
  ! Options:Guess=
  CHARACTER(LEN=*),  PARAMETER :: GUESS_OPTION       ='Guess'
  ! Density matrix from superposition of atomic STO-NG densities
  CHARACTER(LEN=*),  PARAMETER :: GUESS_SUPER        ='SuperPos'
  INTEGER,           PARAMETER :: GUESS_EQ_SUPR      =40823 
  ! Density matrix from core Hamiltonian
  CHARACTER(LEN=*),  PARAMETER :: GUESS_CORE         ='Core'
  INTEGER,           PARAMETER :: GUESS_EQ_CORE      =14334 
  ! Possibly restart from HDFFile density matrix  
  CHARACTER(LEN=*),  PARAMETER :: GUESS_RESTART      ='Restart'
  INTEGER,           PARAMETER :: GUESS_EQ_RESTART   =34344
  ! These are Guess options (hardwired) for solving CPSCF equations 
  ! the polarizability
  CHARACTER(LEN=*),  PARAMETER :: CPSCF_OPTION       ='CPSCF'   ! unused for now
  CHARACTER(LEN=*),  PARAMETER :: GUESS_DIPOLE       ='Dipole'  ! unused for now
  INTEGER,           PARAMETER :: GUESS_EQ_DIPOLE    = 2834032  ! hard wired for now
  CHARACTER(LEN=*),  PARAMETER :: RESTART_INFO       ='HDFFile'
  !------------------------------------------------------------------------------
  ! Options:Restart=
  CHARACTER(LEN=*),  PARAMETER :: RESTART_OPTION     ='Restart'   
  CHARACTER(LEN=*),  PARAMETER :: RESTART_NEWGEOM    ='ParseGeom'
  CHARACTER(LEN=*), PARAMETER :: RESTART_TC2PERT    ='TC2Perturb'
  CHARACTER(LEN=*), PARAMETER :: TC2PERT_OLDFOCK    ='OldFockMat'
  !------------------------------------------------------------------------------
  ! Options:InkFok=
  CHARACTER(LEN=*),  PARAMETER :: INKFOCK_OPTION     ='InkFok'
  CHARACTER(LEN=*),  PARAMETER :: INKFOCK_ON         ='On'
  !------------------------------------------------------------------------------
  ! Options:Grad=  
  CHARACTER(LEN=*),   PARAMETER :: GRADIENTS         ='Grad'
  ! Do no gradeint evaluation
  INTEGER,            PARAMETER :: GRAD_NO_GRAD      = 1000001 
  ! Perform one force evaluation, with print out of the forces
  CHARACTER(LEN=*),   PARAMETER :: GRAD_FORCE        ='OneForce'
  INTEGER,            PARAMETER :: GRAD_ONE_FORCE    = 1084814 
  ! Geometry optimization
  CHARACTER(LEN=*),   PARAMETER :: GRAD_OPTIMIZE     ='Optimize'
  INTEGER,            PARAMETER :: GRAD_GO_DOWNHILL  = 3489343 
  ! Molecular dynamics 
  CHARACTER(LEN=*),   PARAMETER :: GRAD_DYNAMICS     = 'MD'      
  INTEGER,            PARAMETER :: GRAD_DO_DYNAMICS  = 6413123 
  ! Transition state
  CHARACTER(LEN=*),   PARAMETER :: GRAD_TS_SEARCH    ='TSSearch'
  INTEGER,            PARAMETER :: GRAD_TS_SEARCH_NEB= 3577711 
  ! GDIIS 
  CHARACTER(LEN=*),   PARAMETER :: GRAD_GDIIS        ='GDIIS'
  ! Approximate diagonal Hessian in internal coordinates
  CHARACTER(LEN=*),  PARAMETER :: GRAD_APPRX_HESS   ='ApproxHessian'
  ! Coordinate types for gradient operations
  CHARACTER(LEN=*),   PARAMETER  :: GRAD_INTERNALS   ='PrimInt'
  INTEGER,            PARAMETER  :: GRAD_INTS_OPT    =83458086
  CHARACTER(LEN=*),   PARAMETER  :: GRAD_CARTESIAN   ='Cartesian'
  INTEGER,            PARAMETER  :: GRAD_CART_OPT    =34243945
  !------------------------------------------------------------------------------
  ! Options:NEB=
  CHARACTER(LEN=*),  PARAMETER :: NEB_OPTION       ='NEB'
  ! Climbing Image
  CHARACTER(LEN=*),   PARAMETER :: NEB_CLIMB         ='NEBClimb'
  ! Start NEB from reactants and products HDF files
  CHARACTER(LEN=*),  PARAMETER :: NEB_READ_HDF       ='HDF'
  INTEGER,           PARAMETER :: ENDPOINTS_FROM_HDF =456675
  ! HDF files to read enpoint info from 
  CHARACTER(LEN=*),  PARAMETER :: NEB_REACTANTS_HDF  ='Reactants'
  CHARACTER(LEN=*),  PARAMETER :: NEB_PRODUCTS_HDF   ='Products'
  ! Read NEB reactants and products geometries from input file
  CHARACTER(LEN=*),  PARAMETER :: NEB_READ_INP       ='Inp'
  INTEGER,           PARAMETER :: ENDPOINTS_FROM_INP =134455
  !
  CHARACTER(LEN=*),   PARAMETER :: NEB_SPRING        ='NEBSpring'
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  ! Options: Output=
  CHARACTER(LEN=*),  PARAMETER :: OUTPUT_OPTION       ='Output'
  !
  CHARACTER(LEN=*),  PARAMETER :: OUTPUT_PDB          ='PDB'
  INTEGER,           PARAMETER :: PDB_FILE	      = -40
  !
  CHARACTER(LEN=*),  PARAMETER :: OUTPUT_XYZ          ='XYZ'
  INTEGER,           PARAMETER :: XYZ_FILE	      = -41
  !
  CHARACTER(LEN=*),  PARAMETER :: OUTPUT_XCD          ='XSF'
  INTEGER,           PARAMETER :: XSF_FILE	      = -42
END MODULE OptionKeys
