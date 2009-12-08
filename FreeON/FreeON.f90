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
!==============================================================

#include "MondoConfig.h"

PROGRAM FreeON
  USE SCFs
  USE Macros
  USE Response
  USE PunchHDF
  USE Optimizer
  USE ParseInput
  USE ZippyQuote
  USE MDynamics
  USE MonteCarlo
  USE MondoLogger
  USE RayleighQuotientIteration

#if defined(PARALLEL) || defined(PARALLEL_CLONES)
  USE MondoMPI
#endif

  IMPLICIT NONE

  TYPE(Controls) :: C

  CALL Init(PerfMon)
  CALL Init(MemStats)

#if (defined(PARALLEL) || defined(PARALLEL_CLONES)) && defined(MPI2)
  CALL InitMPI()
  InParallel = .FALSE.
  IF(MyID == 0) THEN
#endif

    ! Load the input file.
    CALL ParseTheInput(C)

    ! Initialize controls
    CALL InitGlobal(C)

    ! Much ado about gradients
    SELECT CASE(C%Opts%Grad)
    CASE(GRAD_NO_GRAD)
      CALL SinglePoints(C)
      CALL CPSCF(C)
      CALL TDSCF(C)
    CASE(GRAD_GO_DOWNHILL)
      CALL Descender(C)
    CASE(GRAD_TS_SEARCH_NEB)
      CALL Descender(C)
    CASE(GRAD_DO_DYNAMICS)
      CALL MD(C)
    CASE(GRAD_DO_HYBRIDMC)
      CALL HybridMC(C)
    CASE(GRAD_DO_SCAN)
      CALL Halt('SCAN Not Implimented')
      ! CALL ScanGeom(C)
    CASE(GRAD_ONE_FORCE)
      CALL SinglePoints(C)
      CALL Force(C%Sets%NBSets,1,C%Nams,C%Opts,C%Stat,C%Geos,C%Sets,C%MPIs)
    CASE(GRAD_DO_NHESSIAN)
      CALL NHessian(C)
    END SELECT
    ! Something surreal to celebrate this run
    CALL ZippySez(C)
    !--------------------------------------------------------
    CALL TimeStamp('Successful FreeON run',.FALSE.)
    !--------------------------------------------------------
#if (defined(PARALLEL) || defined(PARALLEL_CLONES)) && defined(MPI2)
  ENDIF

  CALL AlignFrontends()
#endif

END PROGRAM FreeON
