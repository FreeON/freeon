MODULE GeomOptKeys
   USE GlobalScalars
   USE GlobalCharacters
   IMPLICIT NONE
!------------------------------------------------  
!
   CHARACTER(LEN=3),  PARAMETER :: OPTIMIZATION     ='Opt'
   CHARACTER(LEN=5),  PARAMETER :: OPT_QUNEW        ='QuNew'
   CHARACTER(LEN=7),  PARAMETER :: OPT_StpDesc      ='StpDesc'
   CHARACTER(LEN=8),  PARAMETER :: OPT_DiagHess     ='DiagHess'
   CHARACTER(LEN=6),  PARAMETER :: OPT_BiSect       ='BiSect'
   CHARACTER(LEN=8),  PARAMETER :: OPT_CartDIIS     ='CartDIIS'
   CHARACTER(LEN=7),  PARAMETER :: OPT_IntDIIS      ='IntDIIS'
   CHARACTER(LEN=7),  PARAMETER :: OPT_NoGDIIS      ='NoGDIIS'
   CHARACTER(LEN=10), PARAMETER :: OPT_GradNorm     ='DoGradNorm'
   CHARACTER(LEN=8),  PARAMETER :: OPT_Pictures     ='Pictures'
   CHARACTER(LEN=8),  PARAMETER :: OPT_ExplLatt     ='ExplLatt'
   CHARACTER(LEN=9),  PARAMETER :: OPT_Alternate    ='Alternate'
   CHARACTER(LEN=15), PARAMETER :: OPT_FixedAtomsFirst='FixedAtomsFirst'
   CHARACTER(LEN=9),  PARAMETER :: OPT_DoThreeAt    ='DoThreeAt' 
   CHARACTER(LEN=8),  PARAMETER :: OPT_NoBackTr     ='NoBackTr'    
   CHARACTER(LEN=12), PARAMETER :: OPT_DoAtomBackTr ='DoAtomBackTr'
   CHARACTER(LEN=12), PARAMETER :: OPT_DoLattBackTr ='DoLattBackTr'
   CHARACTER(LEN=8),  PARAMETER :: OPT_NoRotOff     ='NoRotOff'
   CHARACTER(LEN=11), PARAMETER :: OPT_NoTranslOff  ='NoTranslOff'
!  Perform quasi-newton geometry optimization for each basis set in turn
   INTEGER, PARAMETER           :: GRAD_QNEW_OPT    = 3489343 
!  Optimizer type is set to Steepest Descent
   INTEGER, PARAMETER           :: GRAD_STPDESC_OPT = 3876123 
!  Optimizer type is set to Diagonal Hessian 
   INTEGER, PARAMETER           :: GRAD_DIAGHESS_OPT = 8942901 
   INTEGER, PARAMETER           :: GRAD_BISECT_OPT   = 8942902 
!--------------------------------------------------------------------------------------
!
    CHARACTER(LEN=9),  PARAMETER :: COORDTYPE='CoordType'
    CHARACTER(LEN=7),  PARAMETER :: CoordType_PrimInt='PrimInt'
    CHARACTER(LEN=9),  PARAMETER :: CoordType_Cartesian='Cartesian'
    CHARACTER(LEN=7),  PARAMETER :: CoordType_FullTrf='FullTrf'
    CHARACTER(LEN=9),  PARAMETER :: CoordType_DoNewChol='DoNewChol'
    CHARACTER(LEN=9),  PARAMETER :: CoordType_DoClssTrf='DoClssTrf'
    CHARACTER(LEN=7),  PARAMETER :: CoordType_DoFixMM='DoFixMM'
    CHARACTER(LEN=8),  PARAMETER :: MaxAtoms='MaxAtoms'
    CHARACTER(LEN=10), PARAMETER :: MaxLattice='MaxLattice'
!-------------------------------------------------
!  Parsing keys for <Options.StpDescInvH=>
!
    CHARACTER(LEN=11),  PARAMETER :: STPDESCINVH='StpDescInvH'
!-------------------------------------------------
!  Parsing keys for <Options.VDWFACT=>
!
    CHARACTER(LEN=7),  PARAMETER :: VDWFACT='VDWFact'
!-------------------------------------------------
!  Parsing keys for <Options.Refresh=>
!
    CHARACTER(LEN=11),  PARAMETER :: INTCREFRESH='IntCRefresh'
!-------------------------------------------------
!  Parsing keys for <Options.MaxAngle=> and <Options.MaxStre=>
!
    CHARACTER(LEN=8),  PARAMETER :: MaxAngle='MaxAngle'
    CHARACTER(LEN=7),  PARAMETER :: MaxStre='MaxStre'
!---------------------------------------------------------
END MODULE GeomOptKeys
