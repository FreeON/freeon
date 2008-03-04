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

#include <MondoConfig.h>

MODULE PrintParsed
  USE InOut
#ifdef NAG
  USE F90_UNIX
#endif
  USE PrettyPrint
  USE ControlStructures
  USE Utilities

  IMPLICIT NONE

CONTAINS
  !===============================================================================================
  ! PRINT THE MONDOSCF BANNER AND THE INPUT FILE TO THE OUTPUT FILE
  !===============================================================================================
  SUBROUTINE PrintsStartUp(N)
    TYPE(FileNames)    :: N
    CHARACTER(LEN=DCL) :: OFile,M_HOST,M_MACH,M_SYST,M_VRSN,M_PLAT
    CHARACTER(LEN=500) :: Line
    CHARACTER(LEN=256) :: hostname
    INTEGER            :: I
    INTEGER            :: HDF5_majnum, HDF5_minnum, HDF5_relnum
    INTEGER            :: StackCurrent, StackMax
    !-----------------------------------------------------------------------------------------------
    CALL GetEnv('MONDO_HOST',M_HOST)
    IF(LEN(TRIM(M_HOST)) == 0) M_HOST = HAVE_MONDO_HOST
    CALL GetEnv('MONDO_MACH',M_MACH)
    IF(LEN(TRIM(M_MACH)) == 0) M_MACH = HAVE_MONDO_MACH
    CALL GetEnv('MONDO_SYST',M_SYST)
    IF(LEN(TRIM(M_SYST)) == 0) M_SYST = HAVE_MONDO_SYST
    CALL GetEnv('MONDO_VRSN',M_VRSN)
    IF(LEN(TRIM(M_VRSN)) == 0) M_VRSN = HAVE_MONDO_VRSN
    CALL GetEnv('MONDO_PLAT',M_PLAT)
    IF(LEN(TRIM(M_PLAT)) == 0) M_PLAT = HAVE_MONDO_PLAT

    ! Print the Input File to the Output File
    OPEN(UNIT=Inp,FILE=N%IFile,STATUS='Old')
    CALL MondoLogPlain("--------- input file = "//TRIM(N%IFile))
    DO
      READ(Inp, '(A500)', END=101) Line
      CALL MondoLogPlain(TRIM(Line))
    ENDDO
101 CONTINUE
    CLOSE(Inp)
    CALL MondoLogPlain("--------- end of input file")

    ! Print MondoSCF banner and authorship
    CALL MondoLogPlain(' __    __                 _       ____________ ______ ')
    CALL MondoLogPlain('|  \  /  |               | |     /       /    |      |')
    CALL MondoLogPlain("|   \/   | ___  _ __   __| | ___/   ____/   __|  |==='")
    CALL MondoLogPlain("|        |/ _ \| '_ \ / _  |/ _ \____  \   (__|  ____|")
    CALL MondoLogPlain('|  |\/|  | (_) | | | | (_| | (_) )     /\     |  |    ')
    CALL MondoLogPlain('|__|  |__|\___/|_| |_|\____|\___/_____/  \____|__|    ')
    CALL MondoLogPlain("")
    CALL MondoLogPlain("Version "//VERSION)
    CALL MondoLogPlain("local version: "//LOCAL_VERSION)
    CALL MondoLogPlain("")
    CALL MondoLogPlain('A program suite for O(N) SCF theory and ab initio MD ')
    CALL MondoLogPlain("")
    CALL MondoLogPlain("Authors:")
    CALL MondoLogPlain('Matt Challacombe, C.J. Tymczak, Karoly Nemeth,       ')
    CALL MondoLogPlain('Valery Weber, Chee Kwan Gan, Eric Schwegler,         ')
    CALL MondoLogPlain('Nicolas Bock, Anders Niklasson, and Graeme Henkelman ')
    CALL MondoLogPlain("")
    CALL MondoLogPlain('Los Alamos National Laboratory                       ')
    CALL MondoLogPlain('LA-CC-04-086 (formerly 01-2)                         ')
    CALL MondoLogPlain('Copyright 2001, University of California.            ')

    ! Get information about hdf version used.
    CALL HDF5Version(HDF5_majnum, HDF5_minnum, HDF5_relnum)

    ! Get information about stack.
    CALL GetStacksizeLimit(StackCurrent, StackMax)

    ! Get hostname.
    CALL GetHostnameWrapper(hostname, 256)

    ! Write information on host, platform, etc
    CALL MondoLogPlain("")
    CALL MondoLogPlain("*** Summary ***")
    CALL MondoLogPlain("")
    CALL MondoLogPlain("running on "//TRIM(hostname))
    CALL MondoLogPlain("")
    CALL MondoLogPlain("** Build Information **")
    CALL MondoLogPlain("")
    CALL MondoLogPlain("compiled for "//TRIM(M_PLAT))
    CALL MondoLogPlain("on "//TRIM(M_HOST)//", a "//TRIM(M_MACH)//" machine")
    CALL MondoLogPlain("running "//TRIM(M_SYST)//" "//TRIM(M_VRSN))
    CALL MondoLogPlain("using HDF5 library version " &
      //TRIM(IntToChar(HDF5_majnum))//"." &
      //TRIM(IntToChar(HDF5_minnum))//"." &
      //TRIM(IntToChar(HDF5_relnum)))
    IF(StackCurrent < 0) THEN
      CALL MondoLogPlain("current stacksize limit: unlimited")
    ELSE
      CALL MondoLogPlain("current stacksize limit: "//TRIM(FltToShrtChar(StackCurrent/1024.0D0/1024.0D0))//" MB")
      CALL MondoLogPlain("")
      CALL MondoLogPlain("Notice:")
      CALL MondoLogPlain("In case you see unexplained segmentation violations, try")
      CALL MondoLogPlain("raising the stacksize limit.")
      CALL MondoLogPlain("")
    ENDIF
    CALL MondoLogPlain("current HDF file = "//TRIM(N%HFile))
    CALL MondoLogPlain("")

    ! Print out a timestamp.
    CALL TimeStamp("Starting MondoSCF")
    CALL MondoLogPlain("")

  END SUBROUTINE PrintsStartUp

END MODULE PrintParsed
