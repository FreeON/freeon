MODULE OptionKeys
  IMPLICIT NONE
  !  Generic options
  CHARACTER(LEN=14), PARAMETER :: OPTIONS_BEGIN     ='<BeginOptions>'
  CHARACTER(LEN=12), PARAMETER :: OPTIONS_END       ='<EndOptions>'
  ! Options:SCF=
  CHARACTER(LEN=9),  PARAMETER :: SCF_OPTION        ='SCFMethod' 
  ! Restricted Roothaan-Hall 
  CHARACTER(LEN=2),  PARAMETER :: SCF_RHHF          ='RH'
  INTEGER,           PARAMETER :: RH_R_SCF          =10982348 
  ! Restricted Simplified Density Matrix Minimization 
  CHARACTER(LEN=4),  PARAMETER :: SCF_SDMM          ='SDMM' 
  INTEGER,           PARAMETER :: SDMM_R_SCF        =58638502 
  ! Restricted Palser-Manulopolus Purification 
  CHARACTER(LEN=2),  PARAMETER :: SCF_PM            ='PM' 
  INTEGER,           PARAMETER :: PM_R_SCF          =58638503 
  ! Restricted quadratic Trace Correcting Purification 
  CHARACTER(LEN=3),  PARAMETER :: SCF_SP2           ='TC2'
  INTEGER,           PARAMETER :: SP2_R_SCF         =58638504 
  ! Restricted quartic Trace Correcting Purification
  CHARACTER(LEN=3),  PARAMETER :: SCF_SP4           ='TC4'
  INTEGER,           PARAMETER :: SP4_R_SCF         =58638505 
  ! Restricted quartic Trace Re-Setting Purification 
  CHARACTER(LEN=4),  PARAMETER :: SCF_TS4           ='TRS4'
  INTEGER,           PARAMETER :: TS4_R_SCF         =58638506 
  !---------------------------------------------------------
  ! Options:Accuracy=
  CHARACTER(LEN=8),  PARAMETER :: ACCURACY_OPTION   ='Accuracy' 
  CHARACTER(LEN=6),  PARAMETER :: ACCURACY_CHEEZY   ='Loose'
  CHARACTER(LEN=4),  PARAMETER :: ACCURACY_GOOD     ='Good'
  CHARACTER(LEN=5),  PARAMETER :: ACCURACY_TIGHT    ='Tight'
  CHARACTER(LEN=9),  PARAMETER :: ACCURACY_RETENTIVE='VeryTight'
  !------------------------------------------------------------------------------
  ! Options:Guess=
  CHARACTER(LEN=5),  PARAMETER :: GUESS_OPTION      ='Guess'
  ! Density matrix from superposition of atomic STO-NG densities
  CHARACTER(LEN=8),  PARAMETER :: GUESS_SUPER       ='SuperPos'
  INTEGER,           PARAMETER :: GUESS_EQ_SUPR     =40823 
  ! Density matrix from core Hamiltonian
  CHARACTER(LEN=4),  PARAMETER :: GUESS_CORE        ='Core'
  INTEGER,           PARAMETER :: GUESS_EQ_CORE     =14334 
  ! Possibly restart from HDFFile density matrix  
  CHARACTER(LEN=7),  PARAMETER :: GUESS_RESTART     ='Restart'
  INTEGER,           PARAMETER :: GUESS_EQ_RESTART  =34344
  CHARACTER(LEN=7),  PARAMETER :: RESTART_INFO      ='HDFFile'
  !------------------------------------------------------------------------------
  ! Options:InkFok=
  CHARACTER(LEN=6),  PARAMETER :: INKFOCK_OPTION    ='InkFok'
  CHARACTER(LEN=2),  PARAMETER :: INKFOCK_ON        ='On'
  !------------------------------------------------------------------------------
  ! Options:Grad=  
  CHARACTER(LEN=4),   PARAMETER :: GRADIENTS         ='Grad'
  ! Do no gradeint evaluation
  INTEGER,            PARAMETER :: GRAD_NO_GRAD      = 1000001 
  ! Perform one force evaluation, with print out of the forces
  CHARACTER(LEN=8),   PARAMETER :: GRAD_FORCE        ='OneForce'
  INTEGER,            PARAMETER :: GRAD_ONE_FORCE    = 1084814 
  ! Perform quasi-newton geometry optimization for each basis set in turn
  CHARACTER(LEN=5),   PARAMETER :: GRAD_QUNEW        ='QuNew'
  INTEGER,            PARAMETER :: GRAD_QNEW_OPT     = 3489343 
  ! Optimizer type is set to Steepest Descent
  CHARACTER(LEN=7),   PARAMETER :: GRAD_SDESCENT     ='StpDesc'
  INTEGER,            PARAMETER :: GRAD_SDESCENT_OPT = 3876123 
  ! Perform optimization/dynamics on all basis sets in sequence
  CHARACTER(LEN=8),   PARAMETER :: GRAD_ALL_BASIS    ='AllBasis'
  ! Constant energy molecular dynamics with the velocity verlet integrator
  CHARACTER(LEN=2),   PARAMETER :: GRAD_DYNAMICS     = 'MD'      
  INTEGER,            PARAMETER :: GRAD_DO_DYNAMICS  = 6413123 
  ! Perform a gradients only transition state search using nugged elastic band
  ! See below for <NEBand>  key words and keys
  CHARACTER(LEN=8),   PARAMETER :: GRAD_TS_SEARCH    ='TSSearch'
  INTEGER,            PARAMETER :: GRAD_TS_SEARCH_NEB= 3577711 
  ! To use GDIIS or not 
  CHARACTER(LEN=8),   PARAMETER :: GRAD_GDIIS        ='GDIIS'
  ! Coordinate type to use for gradient operations
  CHARACTER(LEN=9),   PARAMETER  :: GRAD_INTERNALS   ='Internals'
  INTEGER,            PARAMETER  :: GRAD_INTS_OPT    =83458086
  CHARACTER(LEN=9),   PARAMETER  :: GRAD_CARTESIAN   ='Cartesian'
  INTEGER,            PARAMETER  :: GRAD_CART_OPT    =34243945
  !------------------------------------------------------------------------------
END MODULE OptionKeys
