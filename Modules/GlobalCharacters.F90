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
MODULE GlobalCharacters
   USE GlobalScalars
   IMPLICIT NONE
!------------------------------------------------
!  Lengths and formats for default character string stuff
!
   INTEGER, PARAMETER          :: DEFAULT_CHR_LEN=256
   INTEGER, PARAMETER          :: DCL=DEFAULT_CHR_LEN
   CHARACTER(LEN=*), PARAMETER :: DEFAULT_CHR_FMT='(A256)'
!-------------------------------------------------
!  Lengths and formats for internal IO
!
   INTEGER, PARAMETER          :: INTERNAL_INT_LEN=22
   INTEGER, PARAMETER          :: INTERNAL_DBL_LEN=22
   INTEGER, PARAMETER          :: INTERNAL_FLT_LEN=22
   CHARACTER(LEN=*), PARAMETER :: INTERNAL_INT_FMT='(I22)'
   CHARACTER(LEN=*),PARAMETER  :: INTERNAL_DBL_FMT='(D22.16)'
   CHARACTER(LEN=*), PARAMETER :: INTERNAL_FLT_FMT='(F22.16)'
!   CHARACTER(LEN=*), PARAMETER :: INTERNAL_CHR_FMT='(A)'
!-------------------------------------------------
!  Lengths and formats for basis sets
!
   INTEGER, PARAMETER          :: BASESET_CHR_LEN=16
!-------------------------------------------------
!  Lengths and formats for string output
!
   CHARACTER(LEN=*), PARAMETER :: OUTPUT_STR_FMT='(2x,A)'
!-------------------------------------------------
!  Environmental variables
!
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: MONDO_PWD
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: FREEON_HOME
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: MONDO_SCRATCH
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: MONDO_EXEC
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: MONDO_HOST
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: MONDO_MACH
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: MONDO_SYST
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: MONDO_VRSN
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: MONDO_PLAT
!--------------------------------------------------------
!  File names
!
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: ScrName
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: PWDName
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: InpFile
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: OutFile
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: InfFile
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: GeoFile
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: LogFile
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: BasFile
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: Restart
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: H5File
!-------------------------------------------------
!  File postfixes
!
   CHARACTER(LEN=*), PARAMETER :: InpF='.inp'
   CHARACTER(LEN=*), PARAMETER :: OutF='.out'
   CHARACTER(LEN=*), PARAMETER :: InfF='.hdf'
   CHARACTER(LEN=*), PARAMETER :: GeoF='.pdb'
   CHARACTER(LEN=*), PARAMETER :: LogF='.log'
   CHARACTER(LEN=*), PARAMETER :: BasF='.bas'
!-------------------------------------------------
!  Title variable
!
   INTEGER, PARAMETER          :: TITLE_LINES=10
   CHARACTER(LEN=TITLE_LINES*DEFAULT_CHR_LEN)  :: SCFTitle
!-------------------------------------------------
!  Misc. character variables
!
   CHARACTER(LEN=*), PARAMETER :: Rtrn=CHAR(10)
   CHARACTER(LEN=*), PARAMETER :: BakSlash=CHAR(92)
   CHARACTER(LEN=*), PARAMETER :: LeftParenStar='(* '//Rtrn
   CHARACTER(LEN=*), PARAMETER :: RightParenStar=Rtrn//' *)'
   CHARACTER(LEN=*), PARAMETER :: Blanks= &
     '                                                                '
   CHARACTER(LEN=*), PARAMETER :: Blnk=' '
   INTEGER,          PARAMETER :: IBlnk=ICHAR(Blnk)
   CHARACTER(LEN=*), PARAMETER :: Delta='/'//BakSlash
   INTEGER                     :: RecycleHDF

!
! SCF global characters
!
   CHARACTER(LEN=3)  :: SCFCycl
   CHARACTER(LEN=3)  :: PrvCycl
   CHARACTER(LEN=3)  :: CurCycl
   CHARACTER(LEN=3)  :: NxtCycl
   CHARACTER(LEN=20) :: SCFActn
   CHARACTER(LEN=3)  :: CurBase
   CHARACTER(LEN=3)  :: PrvBase
   CHARACTER(LEN=6)  :: CurGeom
   CHARACTER(LEN=6)  :: PrvGeom
   CHARACTER(LEN=6)  :: NxtGeom
   CHARACTER(LEN=6)  :: CurClone
!
END MODULE
