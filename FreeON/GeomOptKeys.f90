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
MODULE GeomOptKeys
   USE GlobalScalars
   USE GlobalCharacters
   IMPLICIT NONE
!------------------------------------------------
!
   CHARACTER(LEN=3),  PARAMETER :: OPTIMIZATION       ='Opt'
   CHARACTER(LEN=5),  PARAMETER :: OPT_QUNEW          ='QuNew'
   CHARACTER(LEN=7),  PARAMETER :: OPT_StpDesc        ='StpDesc'
   CHARACTER(LEN=8),  PARAMETER :: OPT_DiagHess       ='DiagHess'
   CHARACTER(LEN=6),  PARAMETER :: OPT_BiSect         ='BiSect'
   CHARACTER(LEN=8),  PARAMETER :: OPT_CartDIIS       ='CartDIIS'
   CHARACTER(LEN=7),  PARAMETER :: OPT_IntDIIS        ='IntDIIS'
   CHARACTER(LEN=7),  PARAMETER :: OPT_NoGDIIS        ='NoGDIIS'
   CHARACTER(LEN=10), PARAMETER :: OPT_GradNorm       ='DoGradNorm'
   CHARACTER(LEN=8),  PARAMETER :: OPT_Pictures       ='Pictures'
   CHARACTER(LEN=9),  PARAMETER :: OPT_PrtBackTr      ='PrtBackTr'
   CHARACTER(LEN=7),  PARAMETER :: OPT_NoBTRep        ='NoBTRep'
   CHARACTER(LEN=8),  PARAMETER :: OPT_ExplLatt       ='ExplLatt'
   CHARACTER(LEN=9),  PARAMETER :: OPT_DoQFilter      ='DoQFilter'
   CHARACTER(LEN=9),  PARAMETER :: OPT_Alternate      ='Alternate'
   CHARACTER(LEN=12), PARAMETER :: OPT_LatticeStart   ='LatticeStart'
   CHARACTER(LEN=8),  PARAMETER :: OPT_RatioABC       ='RatioABC'
   CHARACTER(LEN=14), PARAMETER :: OPT_RatioAlpBetGam ='RatioAlpBetGam'
   CHARACTER(LEN=9),  PARAMETER :: OPT_DoThreeAt      ='DoThreeAt'
   CHARACTER(LEN=8),  PARAMETER :: OPT_NoBackTr       ='NoBackTr'
   CHARACTER(LEN=12), PARAMETER :: OPT_DoAtomBackTr   ='DoAtomBackTr'
   CHARACTER(LEN=12), PARAMETER :: OPT_DoLattBackTr   ='DoLattBackTr'
   CHARACTER(LEN=8),  PARAMETER :: OPT_NoRotOff       ='NoRotOff'
   CHARACTER(LEN=11), PARAMETER :: OPT_NoTranslOff    ='NoTranslOff'
   CHARACTER(LEN=10), PARAMETER :: OPT_NonCovTors     ='NonCovTors'
   CHARACTER(LEN=10), PARAMETER :: OPT_NonCovBend     ='NonCovBend'
   CHARACTER(LEN=9),  PARAMETER :: OPT_HBondOnly      ='HBondOnly'
   CHARACTER(LEN=14), PARAMETER :: OPT_NoFragmConnect ='NoFragmConnect'
!  Perform quasi-newton geometry optimization for each basis set in turn
   INTEGER, PARAMETER           :: GRAD_QNEW_OPT      = 3489343
!  Optimizer type is set to Steepest Descent
   INTEGER, PARAMETER           :: GRAD_STPDESC_OPT   = 3876123
!  Optimizer type is set to Diagonal Hessian
   INTEGER, PARAMETER           :: GRAD_DIAGHESS_OPT  = 8942901
   INTEGER, PARAMETER           :: GRAD_BISECT_OPT    = 8942902
!--------------------------------------------------------------------------------------
!
    CHARACTER(LEN=9),  PARAMETER :: COORDTYPE           ='CoordType'
    CHARACTER(LEN=7),  PARAMETER :: CoordType_PrimInt   ='PrimInt'
    CHARACTER(LEN=9),  PARAMETER :: CoordType_Cartesian ='Cartesian'
    CHARACTER(LEN=7),  PARAMETER :: CoordType_FullTrf   ='FullTrf'
    CHARACTER(LEN=9),  PARAMETER :: CoordType_DoNewChol ='DoNewChol'
    CHARACTER(LEN=9),  PARAMETER :: CoordType_DoClssTrf ='DoClssTrf'
    CHARACTER(LEN=7),  PARAMETER :: CoordType_DoFixMM   ='DoFixMM'
    CHARACTER(LEN=12), PARAMETER :: MaxAtomSteps        ='MaxAtomSteps'
    CHARACTER(LEN=15), PARAMETER :: MaxLatticeSteps     ='MaxLatticeSteps'
!-------------------------------------------------
!  Parsing keys for <Options.StpDescInvH=>
!
    CHARACTER(LEN=11),  PARAMETER :: STPDESCINVH ='StpDescInvH'
!-------------------------------------------------
!  Parsing keys for <Options.VDWFACT=>
!
    CHARACTER(LEN=7),  PARAMETER :: VDWFACT ='VDWFact'
!-------------------------------------------------
!  Parsing keys for <Options.Refresh=>
!
    CHARACTER(LEN=11),  PARAMETER :: INTCREFRESH ='IntCRefresh'
!-------------------------------------------------
!  Parsing keys for <Options.MaxAngle=> and <Options.MaxStre=>
!
    CHARACTER(LEN=8),  PARAMETER :: MaxAngle ='MaxAngle'
    CHARACTER(LEN=7),  PARAMETER :: MaxStre  ='MaxStre'
!---------------------------------------------------------
END MODULE GeomOptKeys
