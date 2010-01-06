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
MODULE OptionKeys
  IMPLICIT NONE
  !  Generic options
  CHARACTER(LEN=*), PARAMETER :: OPTIONS_BEGIN      ='<BeginOptions>'
  CHARACTER(LEN=*), PARAMETER :: OPTIONS_END        ='<EndOptions>'
  ! Options:SCF=
  CHARACTER(LEN=*), PARAMETER :: SCF_OPTION         ='SCFMethod'
  ! Restricted Roothaan-Hall
  CHARACTER(LEN=*), PARAMETER :: SCF_RHHF           ='RH'
  INTEGER,          PARAMETER :: RH_R_SCF           =10982348
  ! Restricted Simplified Density Matrix Minimization
  CHARACTER(LEN=*), PARAMETER :: SCF_SDMM           ='SDMM'
  INTEGER,          PARAMETER :: SDMM_R_SCF         =58638502
  ! Restricted Palser-Manulopolus Purification
  CHARACTER(LEN=*), PARAMETER :: SCF_PM             ='PM'
  INTEGER,          PARAMETER :: PM_R_SCF           =58638503
  ! Restricted quadratic Trace Correcting Purification
  CHARACTER(LEN=*), PARAMETER :: SCF_SP2            ='TC2'
  INTEGER,          PARAMETER :: SP2_R_SCF          =58638504
  ! Restricted quartic Trace Correcting Purification
  CHARACTER(LEN=*), PARAMETER :: SCF_SP4            ='TC4'
  INTEGER,          PARAMETER :: SP4_R_SCF          =58638505
  ! Restricted quartic Trace Re-Setting Purification
  CHARACTER(LEN=*), PARAMETER :: SCF_TS4            ='TRS4'
  INTEGER,          PARAMETER :: TS4_R_SCF          =58638506
  !------------------------------------------------------------------------------
  ! Options:Convergence=
  !
  CHARACTER(LEN=*), PARAMETER :: CNVRGNC_OPTION ='SCFConvergence'
  !
  CHARACTER(LEN=*), PARAMETER :: CNVRGNC_OVRIDE ='OverideSCFConvergence'
  ! Not yet in the basin of convergence,
  ! because SCF just started, basis set switch, etc
  INTEGER,          PARAMETER :: NOTYET_CNVRGNC =0
  ! DIIS
  CHARACTER(LEN=*), PARAMETER :: CNVRGNC_DIIS   ='DIIS'
  INTEGER,          PARAMETER :: DIIS_CNVRGNC   =42424242
  ! DIIS with incremental Fock builds
  CHARACTER(LEN=*), PARAMETER :: CNVRGNC_INCF   ='DIIS-INC'
  INTEGER,          PARAMETER :: INCF_CNVRGNC   =42734534
  ! ODA
  CHARACTER(LEN=*), PARAMETER :: CNVRGNC_ODA    ='ODA'
  INTEGER,          PARAMETER :: ODA_CNVRGNC    =43434343
  ! ODA then DIIS
  CHARACTER(LEN=*), PARAMETER :: CNVRGNC_ODMIX  ='ODA-DIIS'
  INTEGER,          PARAMETER :: ODMIX_CNVRGNC  =44444444
  ! DIIS then ODA on Fail
  CHARACTER(LEN=*), PARAMETER :: CNVRGNC_DOMIX  ='DIIS-ODA'
  INTEGER,          PARAMETER :: DOMIX_CNVRGNC  =45454545
  ! Stanton Mixing
  CHARACTER(LEN=*), PARAMETER :: CNVRGNC_SMIX   ='SMIX'
  INTEGER,          PARAMETER :: SMIX_CNVRGNC   =46464646
  !------------------------------------------------------------------------------
  ! Options:Accuracy=
  CHARACTER(LEN=*), PARAMETER :: ACCURACY_OPTION    ='Accuracy'
  CHARACTER(LEN=*), PARAMETER :: ACCURACY_CHEEZY    ='Loose'
  CHARACTER(LEN=*), PARAMETER :: ACCURACY_GOOD      ='Good'
  CHARACTER(LEN=*), PARAMETER :: ACCURACY_TIGHT     ='Tight'
  CHARACTER(LEN=*), PARAMETER :: ACCURACY_RETENTIVE ='VeryTight'
  !------------------------------------------------------------------------------
  ! Options:Guess=
  CHARACTER(LEN=*), PARAMETER :: GUESS_OPTION       ='Guess'
  ! Density matrix from superposition of atomic STO-NG densities
  CHARACTER(LEN=*), PARAMETER :: GUESS_SUPER        ='SuperPos'
  INTEGER,          PARAMETER :: GUESS_EQ_SUPR      =40823
  ! Density matrix from core Hamiltonian
  CHARACTER(LEN=*), PARAMETER :: GUESS_CORE         ='Core'
  INTEGER,          PARAMETER :: GUESS_EQ_CORE      =14334
  ! Possibly restart from HDFFile density matrix
  CHARACTER(LEN=*), PARAMETER :: GUESS_RESTART      ='Restart'
  INTEGER,          PARAMETER :: GUESS_EQ_RESTART   =34344
  CHARACTER(LEN=*), PARAMETER :: GUESS_NUGUESS      ='ReGuess'
  INTEGER,          PARAMETER :: GUESS_EQ_NUGUESS   =45524
  CHARACTER(LEN=*), PARAMETER :: GUESS_NEWGEOM      ='ReParse'
  INTEGER,          PARAMETER :: GUESS_EQ_NEWGEOM   =53223
  !------------------------------------------------------------------------------
  ! These are Guess options (hardwired) for solving CPSCF equations
  ! the polarizability
  CHARACTER(LEN=*), PARAMETER :: CPSCF_OPTION       ='CPSCF'   ! unused for now
  CHARACTER(LEN=*), PARAMETER :: GUESS_DIPOLE       ='Dipole'  ! unused for now
  INTEGER,          PARAMETER :: GUESS_EQ_DIPOLE    = 2834032  ! hard wired for now
  INTEGER,          PARAMETER :: GUESS_EQ_NOGUESS   = 4523123  ! hard wired for now
  CHARACTER(LEN=*), PARAMETER :: RESTART_INFO       ='HDFFile'
  !------------------------------------------------------------------------------
  ! Options:InkFok=
  CHARACTER(LEN=*), PARAMETER :: INKFOCK_OPTION     ='InkFok'
  CHARACTER(LEN=*), PARAMETER :: INKFOCK_ON         ='On'
  !------------------------------------------------------------------------------
  ! Options:Grad=
  CHARACTER(LEN=*),  PARAMETER  :: GRADIENTS         ='Grad'
  ! Do no gradeint evaluation
  INTEGER,           PARAMETER  :: GRAD_NO_GRAD      = 1000001
  ! Maximum number of geometry steps.
  CHARACTER(LEN=*), PARAMETER   :: OPT_NSTEPS        = "MaxSteps"
  ! Perform one force evaluation, with print out of the forces
  CHARACTER(LEN=*),  PARAMETER  :: GRAD_FORCE        ='OneForce'
  INTEGER,           PARAMETER  :: GRAD_ONE_FORCE    = 1084814
  ! Scan from Geometry 1 to Geometry 2
  CHARACTER(LEN=*),  PARAMETER  :: GRAD_SCAN         ='Scan'
  INTEGER,           PARAMETER  :: GRAD_DO_SCAN      = 1666666
  ! Geometry optimization
  CHARACTER(LEN=*),  PARAMETER  :: GRAD_OPTIMIZE     ='Optimize'
  INTEGER,           PARAMETER  :: GRAD_GO_DOWNHILL  = 3489343
  ! Molecular dynamics
  CHARACTER(LEN=*),  PARAMETER  :: GRAD_DYNAMICS     = 'MolecularDynamics'
  INTEGER,           PARAMETER  :: GRAD_DO_DYNAMICS  = 6413123
  ! Hybrid Mondte-Carlo
  CHARACTER(LEN=*),  PARAMETER  :: GRAD_HYBRIDMC     = 'HybridMonteCarlo'
  INTEGER,           PARAMETER  :: GRAD_DO_HYBRIDMC  = 7776665
  ! Transition state
  CHARACTER(LEN=*),  PARAMETER  :: GRAD_TS_SEARCH    ='TSSearch'
  INTEGER,           PARAMETER  :: GRAD_TS_SEARCH_NEB= 3577711
  ! GDIIS
  CHARACTER(LEN=*),  PARAMETER  :: GRAD_GDIIS        ='GDIIS'
  ! Approximate diagonal Hessian in internal coordinates
  CHARACTER(LEN=*), PARAMETER   :: GRAD_APPRX_HESS   ='ApproxHessian'
  ! Coordinate types for gradient operations
  CHARACTER(LEN=*),  PARAMETER  :: GRAD_INTERNALS    ='PrimInt'
  INTEGER,           PARAMETER  :: GRAD_INTS_OPT     =83458086
  CHARACTER(LEN=*),  PARAMETER  :: GRAD_CARTESIAN    ='Cartesian'
  INTEGER,           PARAMETER  :: GRAD_CART_OPT     =34243945
  ! Different cartesian optimizer types.
  CHARACTER(LEN=*),  PARAMETER  :: GRAD_OPTIMIZE_SD                 = "SteepestDescent"
  INTEGER,           PARAMETER  :: GRAD_OPTIMIZE_SD_VALUE           = 1
  CHARACTER(LEN=*),  PARAMETER  :: GRAD_OPTIMIZE_CG                 = "CG"
  INTEGER,           PARAMETER  :: GRAD_OPTIMIZE_CG_VALUE           = 2
  CHARACTER(LEN=*),  PARAMETER  :: GRAD_OPTIMIZE_GLOBAL_CG          = "Global-CG"
  INTEGER,           PARAMETER  :: GRAD_OPTIMIZE_GLOBAL_CG_VALUE    = 3
  CHARACTER(LEN=*),  PARAMETER  :: GRAD_OPTIMIZE_LBFGS              = "LBFGS"
  INTEGER,           PARAMETER  :: GRAD_OPTIMIZE_LBFGS_VALUE        = 4
  CHARACTER(LEN=*),  PARAMETER  :: GRAD_OPTIMIZE_GLOBAL_LBFGS       = "Global-LBFGS"
  INTEGER,           PARAMETER  :: GRAD_OPTIMIZE_GLOBAL_LBFGS_VALUE = 5
  ! Hessian and frequencies
  CHARACTER(LEN=*),  PARAMETER  :: GRAD_NHESSIAN     ='NumFreq'
  INTEGER,           PARAMETER  :: GRAD_DO_NHESSIAN  =13458986
  !------------------------------------------------------------------------------
  ! Options:NEB=
  CHARACTER(LEN=*), PARAMETER :: NEB_OPTION           = 'NEB'
  ! Climbing Image
  CHARACTER(LEN=*), PARAMETER :: NEB_CLIMB            = 'NEBClimb'
  ! Double nudging
  CHARACTER(LEN=*), PARAMETER :: NEB_DOUBLE_NUDGE     = "NEBDoubleNudge"
  ! Steepest Descent parameters.
  CHARACTER(LEN=*), PARAMETER :: NEB_STEEP_ALPHA      = "NEBSteepAlpha"
  CHARACTER(LEN=*), PARAMETER :: NEB_STEEP_MAX_MOVE   = "NEBSteepMaxMove"
  ! Conjugate Gradient parameters.
  CHARACTER(LEN=*), PARAMETER :: CG_STEPSIZE          = "ConjugateGradientdR"
  CHARACTER(LEN=*), PARAMETER :: CG_MAX_MOVE          = "ConjugateGradientMaxMove"
  ! Start NEB from reactants and products HDF files
  CHARACTER(LEN=*), PARAMETER :: NEB_READ_HDF         = 'HDF'
  INTEGER,          PARAMETER :: ENDPOINTS_FROM_HDF   = 456675
  ! HDF files to read enpoint info from
  CHARACTER(LEN=*), PARAMETER :: NEB_REACTANTS_HDF    = 'Reactants'
  CHARACTER(LEN=*), PARAMETER :: NEB_PRODUCTS_HDF     = 'Products'
  ! Read NEB reactants and products geometries from input file
  CHARACTER(LEN=*), PARAMETER :: NEB_READ_INP         = 'Inp'
  INTEGER,          PARAMETER :: ENDPOINTS_FROM_INP   = 134455
  !
  CHARACTER(LEN=*),  PARAMETER :: NEB_SPRING          = 'NEBSpring'
  CHARACTER(LEN=*),  PARAMETER :: NEB_REACTANT_ENERGY = 'NEBReactantEnergy'
  CHARACTER(LEN=*),  PARAMETER :: NEB_PRODUCT_ENERGY  = 'NEBProductEnergy'

  !------------------------------------------------------------------------------
  ! Options:
  CHARACTER(LEN=*), PARAMETER :: Op_MinSCF ='MinSCF'
  CHARACTER(LEN=*), PARAMETER :: Op_MaxSCF ='MaxSCF'
  CHARACTER(LEN=*), PARAMETER :: RQICycles ='MaxRQI'
  CHARACTER(LEN=*), PARAMETER :: RQIGuess  ='RQIGuess'

  !------------------------------------------------------------------------------
  ! Option: misc
  CHARACTER(LEN=*), PARAMETER :: Op_Pressure ='Pressure'

  !------------------------------------------------------------------------------
  ! Options: Output=
  CHARACTER(LEN=*), PARAMETER :: OUTPUT_OPTION       ='Output'

  CHARACTER(LEN=*), PARAMETER :: OUTPUT_PDB          ='PDB'
  INTEGER,          PARAMETER :: PDB_FILE            = -40

  CHARACTER(LEN=*), PARAMETER :: OUTPUT_XYZ          ='XYZ'
  INTEGER,          PARAMETER :: XYZ_FILE            = -41

  CHARACTER(LEN=*), PARAMETER :: OUTPUT_XCD          ='XSF'
  INTEGER,          PARAMETER :: XSF_FILE            = -42

  CHARACTER(LEN=*), PARAMETER :: OUTPUT_CIF          ='CIF'
  INTEGER,          PARAMETER :: CIF_FILE            = -43

  CHARACTER(LEN=*), PARAMETER :: RECYCLE_HDF_OPTION    = "RecycleHDF"
  CHARACTER(LEN=*), PARAMETER :: RANDOM_SCRATCH_OPTION = "RandomScratch"
  CHARACTER(LEN=*), PARAMETER :: CLEAN_SCRATCH_OPTION  = "CleanScratch"

  ! For debugging, use a Lennard-Jones potential. This can be useful to test the
  ! optimizer for instance.
  CHARACTER(LEN=*), PARAMETER :: USE_LENNARD_JONES = "LennardJones"

END MODULE OptionKeys
