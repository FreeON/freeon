!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
MODULE ParsingKeys
   USE GlobalScalars
   USE GlobalCharacters
   IMPLICIT NONE
!------------------------------------------------  
!  General keys for parsing input
!
   CHARACTER(LEN=15), PARAMETER :: BEGIN_GEOMETRY='<BeginGeometry>'
   CHARACTER(LEN=14), PARAMETER :: NEXT_GEOMETRY_MONDO_DEFAULT = &
                                                  '<NextGeometry>'
   CHARACTER(LEN=23), PARAMETER :: BEGIN_NEXT_GEOMETRY_MSI_ARCHIVE   = &
                                                  'input file for discover'
   CHARACTER(LEN=3),  PARAMETER :: END_NEXT_GEOMETRY_MSI_ARCHIVE   =   &
                                                  'end'
   CHARACTER(LEN=13), PARAMETER :: END_GEOMETRY  ='<EndGeometry>'
   CHARACTER(LEN=12), PARAMETER :: BEGIN_TITLE   ='<BeginTitle>'
   CHARACTER(LEN=10), PARAMETER :: END_TITLE     ='<EndTitle>'
   CHARACTER(LEN=14), PARAMETER :: BEGIN_OPTIONS ='<BeginOptions>'
   CHARACTER(LEN=12), PARAMETER :: END_OPTIONS   ='<EndOptions>'
   CHARACTER(LEN=15), PARAMETER :: BEGIN_INFNAME ='<BeginInfoFile>'
   CHARACTER(LEN=13), PARAMETER :: END_INFNAME   ='<EndInfoFile>'
#ifdef PERIODIC
   CHARACTER(LEN=19), PARAMETER :: BEGIN_PERIODIC='<BeginPeriodic>'
   CHARACTER(LEN=17), PARAMETER :: END_PERIODIC  ='<EndPeriodic>'
#endif
!-------------------------------------------------  
!  Parsing keys for <Options.Charge=>
!
   CHARACTER(LEN=6),  PARAMETER :: TOTAL_CHARGE ='Charge' 
!-------------------------------------------------  
!  Parsing keys for <Options.Multiplicity=>
!
   CHARACTER(LEN=12),  PARAMETER :: MULTIPLICITY ='Multiplicity' 
!------------------------------------------------- -------------------------------------------------  
!  Parsing keys for <Options.SCF=>
   CHARACTER(LEN=3),  PARAMETER :: SCF_OPTION='SCFMethod' 
   CHARACTER(LEN=4),  PARAMETER :: SCF_SDMM  ='SDMM' 
   CHARACTER(LEN=2),  PARAMETER :: SCF_RHHF  ='RH'
   CHARACTER(LEN=6),  PARAMETER :: SCF_INKF  ='InkFok'
   CHARACTER(LEN=4),  PARAMETER :: SCF_DIIS  ='DIIS'
!  SCF method keys
   INTEGER, PARAMETER :: RH_R_SCF  =10982348 ! Restricted Roothaan-Hall 
   INTEGER, PARAMETER :: RH_U_SCF  =43240823 ! Unrestricted Roothaan-Hall 
   INTEGER, PARAMETER :: SDMM_R_SCF=58638502 ! Restricted Simplified Density Matrix Minimization 
   INTEGER, PARAMETER :: SDMM_U_SCF=92874018 ! Unrestricted Simplified Density Matrix Minimization 
   INTEGER, PARAMETER :: EDMM_R_SCF=48530687 ! Restricted Extrapolated Density Matrix Minimization 
   INTEGER, PARAMETER :: EDMM_U_SCF=40834070 ! Unrestricted ExtrapolatedDensity Matrix Minimization 
!  SCF guess keys
   INTEGER, PARAMETER :: GUESS_EQ_CORE=14334 ! Density matrix from core Hamiltonian
   INTEGER, PARAMETER :: GUESS_EQ_SUPR=40823 ! Density matrix from superposition of 
                                             ! spherically averaged atomic density matrices
!-----------------------------------------------------------------------------------------------------  
!  Parsing keys for <Options.Guess=>
   CHARACTER(LEN=5),  PARAMETER :: GUESS_OPTION   ='Guess'
   CHARACTER(LEN=4),  PARAMETER :: GUESS_CORE     ='Core'
   CHARACTER(LEN=8),  PARAMETER :: GUESS_SUPER    ='SuperPos'
   CHARACTER(LEN=7),  PARAMETER :: GUESS_RESTART  ='Restart'
   CHARACTER(LEN=7),  PARAMETER :: RESTART_INFO   ='InfoFile'
!---------------------------------------------------------
!  Parsing keys for <Options.Geometry=>
!
   CHARACTER(LEN=8),  PARAMETER :: GEOMETRY     ='Geometry' 
   CHARACTER(LEN=8),  PARAMETER :: IN_AU        ='InAu'
   CHARACTER(LEN=6),  PARAMETER :: H_ORDER      ='HOrder'
   CHARACTER(LEN=6),  PARAMETER :: Z_ORDER      ='ZOrder'
   CHARACTER(LEN=11), PARAMETER :: RANDOM_ORDER ='RandomOrder'
   CHARACTER(LEN=7),  PARAMETER :: NO_ORDER     ='NoOrder'
   CHARACTER(LEN=9),  PARAMETER :: MSI_FORMAT   ='MSIFormat'
   CHARACTER(LEN=10), PARAMETER :: XMOL_FORMAT  ='XMolFormat'
!-----------------------------------------------------------------------------------------------------  
!  Parsing keys for <Options.Force=>
   CHARACTER(LEN=11), PARAMETER :: FACTION      ='ForceAction'   
   CHARACTER(LEN=17), PARAMETER :: MOLDYN       ='MolecularDynamics'
   CHARACTER(LEN=20), PARAMETER :: GEOOPT       ='GeometryOptimization'
   CHARACTER(LEN=5 ), PARAMETER :: FORCE        ='Force'
!-----------------------------------------------------------------------------------------------------  
!  Parsing keys for <Options.Force.MOLDYN=>
   CHARACTER(LEN=11), PARAMETER :: MOLDYN_NS     ='MD_NumSteps' 
   CHARACTER(LEN=11), PARAMETER :: MOLDYN_TS     ='MD_TimeStep' 
   CHARACTER(LEN=11), PARAMETER :: MOLDYN_VS     ='MD_VelScale' 
!---------------------------------------------------------------
!  Parsing keys for <Options.Periodic=>
!
#ifdef PERIODIC
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
