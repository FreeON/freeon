!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 1999, The University of California
!
MODULE GlobalCharacters
   USE GlobalScalars
   IMPLICIT NONE
!------------------------------------------------  
!  Lengths and formats for default character string stuff 
!
   INTEGER, PARAMETER          :: DEFAULT_CHR_LEN=128
   CHARACTER(LEN=6), PARAMETER :: DEFAULT_CHR_FMT='(A128)'
!-------------------------------------------------  
!  Lengths and formats for internal IO
!
   INTEGER, PARAMETER          :: INTERNAL_INT_LEN=22
   INTEGER, PARAMETER          :: INTERNAL_DBL_LEN=22
   INTEGER, PARAMETER          :: INTERNAL_FLT_LEN=22
   CHARACTER(LEN=7), PARAMETER :: INTERNAL_INT_FMT='(I22)'
   CHARACTER(LEN=10),PARAMETER :: INTERNAL_DBL_FMT='(D22.16)'
   CHARACTER(LEN=9), PARAMETER :: INTERNAL_FLT_FMT='(F22.16)'
!   CHARACTER(LEN=5), PARAMETER :: INTERNAL_CHR_FMT='(A)'
!-------------------------------------------------  
!  Lengths and formats for basis sets
!
   INTEGER, PARAMETER          :: BASESET_CHR_LEN=16   
!-------------------------------------------------  
!  Lengths and formats for string output 
!
   CHARACTER(LEN=6), PARAMETER :: OUTPUT_STR_FMT='(2x,A)'
!-------------------------------------------------  
!  File names
! 
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: InpFile
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: OutFile
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: InfFile
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: GeoFile
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: LogFile
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: BasFile
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: SCFName
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: HDFFile
!-------------------------------------------------  
!  File postfixes 
!
   CHARACTER(LEN=4), PARAMETER :: InpF='.inp'
   CHARACTER(LEN=4), PARAMETER :: OutF='.out'
   CHARACTER(LEN=4), PARAMETER :: InfF='.inf'
!- - -- -  - 
!  change to once F77->F90 transition is final
!   CHARACTER(LEN=4), PARAMETER :: InfF='.scf' 
!
   CHARACTER(LEN=4), PARAMETER :: GeoF='.geo'
   CHARACTER(LEN=4), PARAMETER :: LogF='.log'
   CHARACTER(LEN=4), PARAMETER :: BasF='.bas'
!-------------------------------------------------  
!  Title variable
!
   INTEGER, PARAMETER          :: TITLE_LINES=10
   CHARACTER(LEN=TITLE_LINES &
            *DEFAULT_CHR_LEN)  :: SCFTitle
!-------------------------------------------------
!  Misc. character variables
!
   CHARACTER(LEN=1), PARAMETER :: Rtrn=CHAR(10)
   CHARACTER(LEN=1), PARAMETER :: BakSlash=CHAR(92)
   CHARACTER(LEN=5), PARAMETER :: LeftParenStar='(* '//Rtrn
   CHARACTER(LEN=5), PARAMETER :: RightParenStar=Rtrn//' *)'
   CHARACTER(LEN=64), PARAMETER:: Blanks= &
     '                                                               '
   CHARACTER(LEN=1), PARAMETER :: Blnk=' '
   INTEGER,          PARAMETER :: IBlnk=ICHAR(Blnk)
   CHARACTER(LEN=2), PARAMETER :: Delta='/'//BakSlash
END MODULE

