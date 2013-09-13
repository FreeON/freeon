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
MODULE ParsingConstants
  USE GlobalScalars

  IMPLICIT NONE

  !------------------------------------------------
  ! Parsing constants
  CHARACTER(LEN=*), PARAMETER :: Upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  CHARACTER(LEN=*), PARAMETER :: Lower='abcdefghijklmnopqrstuvwxyz'
  CHARACTER(LEN=*), PARAMETER :: Numbers='0123456789-+.0123456789'
  CHARACTER(LEN=*), PARAMETER :: Special='$/_#*='
  CHARACTER(LEN=*), PARAMETER :: Characters=Special//Lower//Numbers
  CHARACTER(LEN=*), PARAMETER :: Delimiters='[(|")],= '
  CHARACTER(LEN=*), PARAMETER :: Stars='****'
  CHARACTER(LEN=*), PARAMETER :: Space=' '

  ! Parsing keys for <Options.Program=> or <Options.GLOBAL_DEBUG=>
  CHARACTER(LEN=*), PARAMETER :: GLOBAL_DEBUG  = 'DebugAll'
  CHARACTER(LEN=*), PARAMETER :: DBG_NONE      = 'NoDebug'
  CHARACTER(LEN=*), PARAMETER :: DBG_MINIMUM   = 'MinDebug'
  CHARACTER(LEN=*), PARAMETER :: DBG_MEDIUM    = 'MedDebug'
  CHARACTER(LEN=*), PARAMETER :: DBG_MAXIMUM   = 'MaxDebug'
  CHARACTER(LEN=*), PARAMETER :: DBG_MATRICES  = 'MatDebug'
  CHARACTER(LEN=*), PARAMETER :: DBG_CHKSUMS   = 'CheckSums'
  CHARACTER(LEN=*), PARAMETER :: DBG_PRT_SETS  = 'BasisSet'
  CHARACTER(LEN=*), PARAMETER :: DBG_PRT_INTS  = 'PrintInt'
  CHARACTER(LEN=*), PARAMETER :: DBG_PRT_RHO   = 'PrintRho'
  CHARACTER(LEN=*), PARAMETER :: PLT_MATRICES  = 'PlotMats'
  CHARACTER(LEN=*), PARAMETER :: DBG_MMA_STYLE = 'MmaStyle'
  CHARACTER(LEN=*), PARAMETER :: DBG_MM_STYLE  = 'MMStyle'
  CHARACTER(LEN=*), PARAMETER :: DBG_DBL_STYLE = 'DblStyle'
  CHARACTER(LEN=*), PARAMETER :: DBG_FLT_STYLE = 'FltStyle'
  CHARACTER(LEN=*), PARAMETER :: DBG_GEOP_MIN  = 'MinGeOp'
  CHARACTER(LEN=*), PARAMETER :: DBG_GEOP_MAX  = 'MaxGeOp'
  CHARACTER(LEN=*), PARAMETER :: DBG_PRT_MM    = 'MMDebug'
  CHARACTER(LEN=*), PARAMETER :: DBG_PRT_FRC   = 'FrcDebug'

  INTEGER, PARAMETER :: SFC_NONE       = 3480481
  INTEGER, PARAMETER :: SFC_PEANO      = 5308208
  INTEGER, PARAMETER :: SFC_HILBERT    = 4808471
  INTEGER, PARAMETER :: SFC_RANDOM     = 5058108
  INTEGER, PARAMETER :: SFC_TRAVEL     = 2505811
  INTEGER, PARAMETER :: SFC_TableTrav  = 2505812

  INTEGER, PARAMETER :: DEBUG_NONE     = 0       ! Nothing
  INTEGER, PARAMETER :: DEBUG_MINIMUM  = 1       ! Timing and memory (Default)
  INTEGER, PARAMETER :: DEBUG_MEDIUM   = 2       ! + sparsity, thresholds
  INTEGER, PARAMETER :: DEBUG_CHKSUMS  = 3       ! + check sums
  INTEGER, PARAMETER :: DEBUG_MAXIMUM  = 4       ! + intermediate values

  INTEGER, PARAMETER :: DEBUG_MATRICES = 482842  ! Print Matrices
  INTEGER, PARAMETER :: DEBUG_DENSITY  = 490485  ! Print Density
  INTEGER, PARAMETER :: DEBUG_BASISSET = 509843  ! Print Basis sets
  INTEGER, PARAMETER :: DEBUG_INTEGRAL = 585583  ! Print Integrals
  INTEGER, PARAMETER :: PLOT_MATRICES  = 608948  ! Plot Matrices
  INTEGER, PARAMETER :: DEBUG_MMASTYLE = 848423  ! Print in Mathematica style
  INTEGER, PARAMETER :: DEBUG_MMSTYLE  = 848424  ! Print in MatrixMarket style
  INTEGER, PARAMETER :: DEBUG_FLTSTYLE = 480484  ! Print in float style
  INTEGER, PARAMETER :: DEBUG_DBLSTYLE = 504843  ! Print in scientific (D) style
  INTEGER, PARAMETER :: DEBUG_GEOP_MIN = 568356  ! Print geometry optimization data
  INTEGER, PARAMETER :: DEBUG_GEOP_MAX = 568357  ! Print geometry optimization data
  INTEGER, PARAMETER :: DEBUG_MM       = 129462  ! Print molecular mechanics related data
  INTEGER, PARAMETER :: DEBUG_FRC      = 680976  ! Print Forceses
END MODULE ParsingConstants
