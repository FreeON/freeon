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
MODULE ParseInput
  USE Massage
  USE ParseBasis
  USE ParseOptions
  USE ParseCommands
  USE ParseDynamics
  USE ParsePeriodic
  USE ParseGeometries
  USE ControlStructures
  USE PrettyPrint
  USE ParseParallel
  USE ParseGeomOpt
  USE ParseExtraCoords
  USE ParseProperties, ONLY: LoadPropertyOptions
  USE PrintParsed
  USE MondoLogger

  IMPLICIT NONE

CONTAINS
  !===============================================================
  ! 'NUFF SAID
  !===============================================================
  SUBROUTINE ParseTheInput(C)
    TYPE(Controls) :: C
    INTEGER        :: iCLONE

    ! Parse command line and load env and file names
    CALL LoadCommands(C%Nams)

    ! Set global output and log file
    OutFile=C%Nams%OFile
    LogFile=C%Nams%LFile

    ! Allocate state variables
    CALL New(C%Stat%Current,3)
    CALL New(C%Stat%Previous,3)

    ! Print startup
    CALL StartUpHeader(C%Nams)

    ! Parse generic options
    CALL LoadOptions(C%Nams,C%Opts)

    ! Parse dynamics options
    IF(C%Opts%Grad==GRAD_DO_DYNAMICS .OR. C%Opts%Grad==GRAD_DO_HYBRIDMC ) THEN
      CALL MondoLog(DEBUG_MAXIMUM, "ParseTheInput", "loading dynamics")
      CALL LoadDynamics(C%Nams,C%Opts,C%Dyns)
    ELSE
      C%Dyns%DoingMD = .FALSE.
    ENDIF

    ! Parse geometry or get from restart HDF
    CALL LoadCoordinates(C%Nams,C%Opts,C%Dyns,C%Geos)

    ! Parse periodic info
    CALL LoadPeriodic(C%Nams,C%Opts,C%Geos,C%PBCs)

    ! Massage coordinates, switch to AUs etc
    CALL MassageCoordinates(C%Opts,C%Geos,C%PBCs)

    ! Load basis sets
    CALL LoadBasisSets(C%Nams,C%Opts,C%Geos,C%Sets)

    ! Parse in parallel info
    CALL LoadParallel(C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)

    ! Load control of internal coord. optimizer
    CALL LoadGeomOpt(C%Nams,C%GOpt,C%Geos%Clone(1)%PBC%Dimen)

    ! Load constraints and extra internal coords
    CALL LoadExtraCoords(C%GOpt,C%Opts,C%Nams,C%Geos)

    ! Load CPSCF options.
    CALL LoadPropertyOptions(C%Nams,C%POpt)

    ! Check for Global conflicts.
    CALL ConflictCheck(C)

  END SUBROUTINE ParseTheInput
END MODULE ParseInput
