!===============================================================================
!
!       Droits de reproduction et de diffusion réservés. © 2000 CEA/CNRS.
!              (Laboratoire de Dynamique Moléculaire/IBS/DSV) 2000
!
!===============================================================================
!
!                Copyright © 2000 CEA/CNRS. All Rights Reserved.
!              (Laboratoire de Dynamique Moléculaire/IBS/DSV) 2000
!
!===============================================================================
!                               The Parsing Module
!===============================================================================
!
! . The current parsing unit:
!
!   UNIT                           The file unit number.
!
! . The current command line:
!
!   CURLIN                         The current line.
!   LINE                           Input line containing the words.
!   LININD                         The line index array for words.
!   NLINES                         The number of lines on a command line.
!   NWORDS                         The number of words on the line.
!   QREAD                          Flags to indicate whether a word has been read.
!   WSTART, WSTOP                  The word starting and stopping places.
!
! . The current command word:
!
!   WORD                           Input word.
!   WRDLEN                         Length of the input word.
!
! . Stream stack:
!
!   NUNITS                         The number of units on the stack.
!   UNITS                          The unit stack.
!
! . Functions:
!
!   EMPTY                          Check to see if the line is empty.
!
! . Subroutines:
!
!   CHECK_UNPARSED_CHARACTERS      Check for unparsed characters.
!   GET_INTEGER                    Get an integer from the command line.
!   GET_LINE                       Read in the next command line.
!   GET_RANGE_INTEGER              Get an integer range.
!   GET_REAL                       Get a real from the command line.
!   GET_STRING                     Get a string from the command line.
!   GET_TAGGED_INTEGER             Get an integer after a tag.
!   GET_TAGGED_LOGICAL             Get a logical after a tag.
!   GET_TAGGED_REAL                Get a real after a tag.
!   GET_TAGGED_REAL_VECTOR         Get a real vector after a tag.
!   GET_TAGGED_STRING              Get a string after a tag.
!   GET_WORD                       Get the next word on the command line.
!   PARSE_ERROR                    There has been a parsing error.
!   POP_UNIT                       Pop the unit for parsing.
!   PUSH_UNIT                      Push the unit for parsing.
!
! . Notes:
!
!   The unit number for the parsing procedures is set using PUSH_UNIT. The
!   units being read are stored on a small stack. The unit is reinitialized
!   using POP_UNIT.
!
!   Lines are read in using GET_LINE. Some processing (although not much) is
!   done. In particular, all characters are converted to upper case and some
!   special characters are converted to blanks. No conversion is done within
!   strings (which are indicated by double quotes ""). Any characters after
!   and including a "!" are treated as comments and are ignored. Errors will
!   occur if the subroutine has problems reading the command line (including
!   an end of file) or if there are unparsed characters on the previous command
!   line. No continuation lines are allowed although multiple lines can occur
!   on a single line each separated by a ";". Note that "!" overrides a ";"
!   although not a double quotes.
!
!   The subroutine GET_WORD picks off strings from the line read and processed
!   by GET_LINE. A word is defined by any string of characters separated from a
!   subsequent string by a space, a comma or an equals sign.
!
!   GET_INTEGER, GET_REAL and GET_STRING are used to return the next word on
!   the command line as an integer, a real or a string.
!
!   PARSE_ERROR is the parsing error subroutine. It calls the normal error
!   subroutine but prints out the current command line first. It is private
!   to the module.
!
!   EMPTY is a logical function that can be used to see if there are extra
!   characters on the command line.
!
!   CHECK_UNPARSED_CHARACTERS checks to see if the current command line is
!   empty. If it is not an error is produced. It is used only within the
!   module.
!
!===============================================================================
MODULE PARSING

! . Module declarations.
USE DEFINITIONS,  ONLY : DP, LINE_LENGTH, MAX_UNITS
USE IO_UNITS,     ONLY : OUTPUT
USE STATUS,       ONLY : ERROR
USE STRING,       ONLY : DECODE_INTEGER, DECODE_LOGICAL, DECODE_REAL, TO_UPPER_CASE

IMPLICIT NONE
PUBLIC
PRIVATE :: CURLIN, LINE, LININD, NLINES, NUNITS, NWORDS, QREAD, UNIT, UNITS, WSTART, WSTOP
SAVE

! . Module parameters.
INTEGER, PARAMETER :: MAX_WORDS = ( LINE_LENGTH + 3 ) / 2

! . Module scalars.
CHARACTER ( LEN = LINE_LENGTH ) :: LINE, WORD
INTEGER                         :: CURLIN = 0, NLINES = 0, NUNITS = 0, NWORDS = 0, UNIT = -1, WRDLEN = 0

! . Module arrays.
INTEGER, DIMENSION(1:MAX_UNITS)   :: UNITS  = 0
INTEGER, DIMENSION(1:MAX_WORDS)   :: WSTART = 0, WSTOP = 0
INTEGER, DIMENSION(1:MAX_WORDS+1) :: LININD = 0
LOGICAL, DIMENSION(1:MAX_WORDS)   :: QREAD  = .FALSE.

!===============================================================================
CONTAINS
!===============================================================================

   !-----------------------------------
   SUBROUTINE CHECK_UNPARSED_CHARACTERS
   !-----------------------------------

   ! . Check to see if characters remain on the line.
   IF ( .NOT. EMPTY ( ) ) THEN
      CALL PARSE_ERROR ( "CHECK_UNPARSED_CHARACTERS", "There are unparsed characters on the line." )
   END IF

   END SUBROUTINE CHECK_UNPARSED_CHARACTERS

   !-----------------
   FUNCTION EMPTY ( )
   !-----------------

   ! . Function declarations.
   LOGICAL :: EMPTY

   ! . Local scalars.
   INTEGER :: N, START, STOP

   ! . Initialize EMPTY.
   EMPTY = .TRUE.

   ! . Check CURLIN.
   IF ( CURLIN <= 0 ) RETURN

   ! . Get the indices for the current line.
   START = LININD(CURLIN) + 1
   STOP  = LININD(CURLIN + 1)
   N     = STOP - START + 1

   ! . Set the logical flag.
   EMPTY = ( COUNT ( QREAD(START:STOP) ) == N )

   END FUNCTION EMPTY

   !--------------------------------
   SUBROUTINE GET_INTEGER ( NUMBER )
   !--------------------------------

   ! . Argument declarations.
   INTEGER, INTENT(OUT) :: NUMBER

   ! . Get the word containing the integer.
   CALL GET_WORD

   ! . Decode the integer.
   NUMBER = DECODE_INTEGER ( WORD(1:WRDLEN) )

   END SUBROUTINE GET_INTEGER

   !----------------------------------
   SUBROUTINE GET_LINE ( END_OF_FILE )
   !----------------------------------

   ! . Scalar arguments.
   LOGICAL, INTENT(OUT), OPTIONAL :: END_OF_FILE

   ! . Local scalars.
   CHARACTER ( LEN = 1 ) :: CHR
   INTEGER               :: I, IOSTAT, LINLEN
   LOGICAL               :: QBLANK, QSEMICOLON, QSEPAR, QSTRING

   ! . Local parameters.
   CHARACTER ( LEN = 1 ), PARAMETER :: TAB = CHAR ( 9 )

   ! . Initialize the end of file argument.
   IF ( PRESENT ( END_OF_FILE ) ) END_OF_FILE = .FALSE.

   ! . Check for unparsed characters.
   CALL CHECK_UNPARSED_CHARACTERS

   ! . There are existing processed lines in LINE.
   IF ( CURLIN < NLINES ) THEN
      CURLIN = CURLIN + 1
      RETURN
   END IF

   ! . Top of the loop for reading lines.
   10 CONTINUE

   ! . Read the line.
   READ ( UNIT, "(A)", END = 20, IOSTAT = IOSTAT ) LINE

   ! . Check for an error.
   IF ( IOSTAT /= 0 ) CALL ERROR ( "GET_LINE", "I/O Error.", IOSTAT )

   ! . Initialization.
   LINLEN     = 0
   QBLANK     = .FALSE.
   QSEMICOLON = .FALSE.
   QSEPAR     = .FALSE.
   QSTRING    = .FALSE.

   ! . Loop over the characters in the line.
   DO I = 1,LINE_LENGTH

      ! . Get the character.
      CHR = LINE(I:I)

      ! . Check for a double quotes.
      IF ( CHR == '"' ) THEN
         QBLANK  = .FALSE.
         QSTRING = .NOT. QSTRING
      END IF

      ! . Do further processing for characters outside a string.
      IF ( .NOT. QSTRING ) THEN

         ! . Substitute various miscellaneous characters for blanks.
         IF ( ( CHR == TAB ) .OR. ( CHR == "=" ) .OR. ( CHR == "," ) ) CHR = " "

         ! . Substitute a blank for a semicolon.
         IF ( CHR == ";" ) THEN
            QSEMICOLON = .TRUE.
            CHR = " "
         END IF

         ! . A comment.
         IF ( CHR == "!" ) THEN
            EXIT

         ! . A blank.
         ELSE IF ( CHR == " " ) THEN
            QBLANK = .TRUE.
            QSEPAR = .TRUE.

         ! . Remaining characters.
         ELSE
            QBLANK = .FALSE.
            CHR = TO_UPPER_CASE ( CHR )
         END IF

      END IF

      ! . There is a non-blank character.
      IF ( .NOT. QBLANK ) THEN

         ! . Check for a previous semicolon or blank separator.
         IF ( QSEPAR ) THEN
            IF ( LINLEN > 0 ) THEN
               LINLEN = LINLEN + 1
               IF ( QSEMICOLON ) THEN
                  LINE(LINLEN:LINLEN) = ";"
               ELSE
                  LINE(LINLEN:LINLEN) = " "
               END IF
            END IF
            QSEMICOLON = .FALSE.
            QSEPAR     = .FALSE.
         END IF

         ! . Write out the character.
         LINLEN = LINLEN + 1
         LINE(LINLEN:LINLEN) = CHR
 
      END IF
   END DO

   ! . Remove a trailing blank.
   IF ( LINLEN > 0 ) THEN
      IF ( LINE(LINLEN:LINLEN) == " " ) LINLEN = LINLEN - 1
   END IF

   ! . Read another line if there are no characters on the line.
   IF ( LINLEN <= 0 ) GO TO 10

   ! . Initialize the line variables for the line.
   CURLIN = 1
   NLINES = 1
   LININD = 0

   ! . Initialize the word variables for the line.
   NWORDS = 0
   QREAD  = .FALSE.

   ! . Initialization.
   QSTRING = .FALSE.

   ! . Fill the data for the first word.
   WSTART(1) = 1

   ! . Fill the word index arrays.
   DO I = 1,LINLEN

      ! . Check for a string.
      IF ( LINE(I:I) == '"' ) QSTRING = .NOT. QSTRING

      ! . There is no string.
      IF ( .NOT. QSTRING ) THEN

         ! . There is a blank or semicolon.
         IF ( ( LINE(I:I) == " " ) .OR. ( LINE(I:I) == ";" ) ) THEN

            ! . The current word has finished.
            NWORDS = NWORDS + 1
            WSTOP(NWORDS)    = I - 1
            WSTART(NWORDS+1) = I + 1

            ! . There is a semicolon.
            IF ( LINE(I:I) == ";" ) THEN

               ! . Increment the number of lines.
               NLINES = NLINES + 1

               ! . Set the index for the line.
               LININD(NLINES) = NWORDS

            END IF
         END IF
      END IF
   END DO

   ! . Fill the data for the last word.
   NWORDS = NWORDS + 1
   WSTOP(NWORDS) = LINLEN

   ! . Fill the last element of LININD.
   LININD(NLINES+1) = NWORDS

   ! . Finished so return.
   RETURN

   ! . End of file condition.
   20 CONTINUE

   ! . Check for the presence of the argument.
   IF ( PRESENT ( END_OF_FILE ) ) THEN
      END_OF_FILE = .TRUE.
   ELSE
      CALL ERROR ( "GET_LINE", "End of file reached." )
   END IF

   END SUBROUTINE GET_LINE

   !------------------------------------------------------------------------
   SUBROUTINE GET_RANGE_INTEGER ( START, STOP, DEFAULT_START, DEFAULT_STOP )
   !------------------------------------------------------------------------

   ! . Argument declarations.
   INTEGER, INTENT(IN)  :: DEFAULT_START, DEFAULT_STOP
   INTEGER, INTENT(OUT) :: START, STOP

   ! . Local scalars.
   INTEGER :: COLON

   ! . Get the word containing the range.
   CALL GET_WORD

   ! . There is a word.
   IF ( WRDLEN > 0 ) THEN

      ! . Check for a colon.
      COLON = INDEX ( WORD, ":" )

      ! . There is no colon.
      IF ( COLON <= 0 ) THEN

         ! . START and STOP are the same.
         START = DECODE_INTEGER ( WORD(1:WRDLEN) )
         STOP  = START

      ! . There is a colon.
      ELSE

         ! . Get the first integer.
         IF ( COLON > 1 ) THEN
            START = DECODE_INTEGER ( WORD(1:(COLON-1)) )
         ! . Use the default value for START.
         ELSE
            START = DEFAULT_START
         END IF

         ! . Get the last integer.
         IF ( COLON < WRDLEN ) THEN
            STOP = DECODE_INTEGER ( WORD((COLON+1):WRDLEN) )
         ! . Use the default value for STOP.
         ELSE
            STOP = DEFAULT_STOP
         END IF

      END IF

   ! . There is no word.
   ELSE

      ! . Return nothing.
      START = 1
      STOP  = 0

   END IF

   END SUBROUTINE GET_RANGE_INTEGER

   !-----------------------------
   SUBROUTINE GET_REAL ( NUMBER )
   !-----------------------------

   ! . Argument declarations.
   REAL ( KIND = DP ), INTENT(OUT) :: NUMBER

   ! . Get the word containing the real.
   CALL GET_WORD

   ! . Decode the real.
   NUMBER = DECODE_REAL ( WORD(1:WRDLEN) )

   END SUBROUTINE GET_REAL

   !--------------------
   SUBROUTINE GET_STRING
   !--------------------

   ! . Get the word containing the string.
   CALL GET_WORD

   ! . There is a word.
   IF ( WRDLEN > 0 ) THEN

      ! . Check that the word is a string.
      IF ( ( WORD(1:1) /= '"' ) .OR. ( WORD(WRDLEN:WRDLEN) /= '"' ) ) THEN
         CALL PARSE_ERROR ( "GET_STRING", "Invalid string constant." )
         WRDLEN = 0

      ! . Decode the string.
      ELSE
         WORD(1:WRDLEN-2) = WORD(2:WRDLEN-1)
         WRDLEN = WRDLEN - 2
      END IF

   END IF

   END SUBROUTINE GET_STRING

   !-----------------------------------------------------
   SUBROUTINE GET_TAGGED_INTEGER ( NUMBER, TAG, DEFAULT )
   !-----------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN)  :: TAG
   INTEGER,               INTENT(IN)  :: DEFAULT
   INTEGER,               INTENT(OUT) :: NUMBER

   ! . Local scalars.
   INTEGER :: IWORD

   ! . Initialize NUMBER.
   NUMBER = DEFAULT

   ! . Loop over the unread words to see if a tag exists.
   DO IWORD = (LININD(CURLIN)+1),LININD(CURLIN+1)

      ! . A word exists that has not been read.
      IF ( .NOT. QREAD(IWORD) ) THEN

         ! . The word matches the tag.
         IF ( LINE(WSTART(IWORD):WSTOP(IWORD)) == TAG ) THEN

            ! . Interpret the next word as an integer.
            IF ( IWORD < LININD(CURLIN+1) ) THEN
               NUMBER = DECODE_INTEGER ( LINE(WSTART(IWORD+1):WSTOP(IWORD+1)) )
            ELSE
               CALL PARSE_ERROR ( "GET_TAGGED_INTEGER", "Missing integer." )
            END IF

            ! . Flag the two words as read and return.
            QREAD(IWORD:IWORD+1) = .TRUE. ; RETURN

         END IF
      END IF
   END DO

   END SUBROUTINE GET_TAGGED_INTEGER

   !----------------------------------------------------
   SUBROUTINE GET_TAGGED_LOGICAL ( QFLAG, TAG, DEFAULT )
   !----------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN)  :: TAG
   LOGICAL,               INTENT(IN)  :: DEFAULT
   LOGICAL,               INTENT(OUT) :: QFLAG

   ! . Local scalars.
   INTEGER :: IWORD

   ! . Initialize QFLAG.
   QFLAG = DEFAULT

   ! . Loop over the unread words to see if a tag exists.
   DO IWORD = (LININD(CURLIN)+1),LININD(CURLIN+1)

      ! . A word exists that has not been read.
      IF ( .NOT. QREAD(IWORD) ) THEN

         ! . The word matches the tag.
         IF ( LINE(WSTART(IWORD):WSTOP(IWORD)) == TAG ) THEN

            ! . Interpret the next word as a logical.
            IF ( IWORD < LININD(CURLIN+1) ) THEN
               QFLAG = DECODE_LOGICAL ( LINE(WSTART(IWORD+1):WSTOP(IWORD+1)) )
            ELSE
               CALL PARSE_ERROR ( "GET_TAGGED_LOGICAL", "Missing logical." )
            END IF

            ! . Flag the two words as read and return.
            QREAD(IWORD:IWORD+1) = .TRUE. ; RETURN

         END IF
      END IF
   END DO

   END SUBROUTINE GET_TAGGED_LOGICAL

   !--------------------------------------------------
   SUBROUTINE GET_TAGGED_REAL ( NUMBER, TAG, DEFAULT )
   !--------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN)  :: TAG
   REAL ( KIND = DP ),    INTENT(IN)  :: DEFAULT
   REAL ( KIND = DP ),    INTENT(OUT) :: NUMBER

   ! . Local scalars.
   INTEGER :: IWORD

   ! . Initialize NUMBER.
   NUMBER = DEFAULT

   ! . Loop over the unread words to see if a tag exists.
   DO IWORD = (LININD(CURLIN)+1),LININD(CURLIN+1)

      ! . A word exists that has not been read.
      IF ( .NOT. QREAD(IWORD) ) THEN

         ! . The word matches the tag.
         IF ( LINE(WSTART(IWORD):WSTOP(IWORD)) == TAG ) THEN

            ! . Interpret the next word as a real.
            IF ( IWORD < LININD(CURLIN+1) ) THEN
               NUMBER = DECODE_REAL ( LINE(WSTART(IWORD+1):WSTOP(IWORD+1)) )
            ELSE
               CALL PARSE_ERROR ( "GET_TAGGED_REAL", "Missing real." )
            END IF

            ! . Flag the two words as read and return.
            QREAD(IWORD:IWORD+1) = .TRUE. ; RETURN

         END IF
      END IF
   END DO

   END SUBROUTINE GET_TAGGED_REAL

   !---------------------------------------------------------
   SUBROUTINE GET_TAGGED_REAL_VECTOR ( NUMBER, TAG, DEFAULT )
   !---------------------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN)  :: TAG
   REAL ( KIND = DP ),    INTENT(IN)  :: DEFAULT

   ! . Array arguments.
   REAL ( KIND = DP ), DIMENSION(:), INTENT(OUT) :: NUMBER

   ! . Local scalars.
   INTEGER :: I, IWORD, N

   ! . Get the size of NUMBER.
   N = SIZE ( NUMBER )

   ! . Check N.
   IF ( N <= 0 ) RETURN

   ! . Initialize NUMBER.
   NUMBER = DEFAULT

   ! . Loop over the unread words to see if a tag exists.
   DO IWORD = (LININD(CURLIN)+1),LININD(CURLIN+1)

      ! . A word exists that has not been read.
      IF ( .NOT. QREAD(IWORD) ) THEN

         ! . The word matches the tag.
         IF ( LINE(WSTART(IWORD):WSTOP(IWORD)) == TAG ) THEN

            ! . Interpret the next N words as reals.
            DO I = (IWORD+1),(IWORD+N)
               IF ( I <= LININD(CURLIN+1) ) THEN
                  NUMBER(I-IWORD) = DECODE_REAL ( LINE(WSTART(I):WSTOP(I)) )
               ELSE
                  CALL PARSE_ERROR ( "GET_TAGGED_REAL_VECTOR", "Missing real." )
               END IF
            END DO

            ! . Flag the N+1 words as read and return.
            QREAD(IWORD:IWORD+N) = .TRUE. ; RETURN

         END IF
      END IF
   END DO

   END SUBROUTINE GET_TAGGED_REAL_VECTOR

   !-------------------------------------------
   SUBROUTINE GET_TAGGED_STRING ( STRING, TAG )
   !-------------------------------------------

   ! . Scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN)  :: TAG
   CHARACTER ( LEN = * ), INTENT(OUT) :: STRING

   ! . Local scalars.
   INTEGER :: IWORD

   ! . Initialize STRING.
   STRING = " "

   ! . Loop over the unread words to see if a tag exists.
   DO IWORD = (LININD(CURLIN)+1),LININD(CURLIN+1)

      ! . A word exists that has not been read.
      IF ( .NOT. QREAD(IWORD) ) THEN

         ! . The word matches the tag.
         IF ( LINE(WSTART(IWORD):WSTOP(IWORD)) == TAG ) THEN

            ! . Interpret the next word as a string.
            IF ( IWORD < LININD(CURLIN+1) ) THEN

               ! . Check that the next word is a string.
               IF ( ( LINE(WSTART(IWORD+1):WSTART(IWORD+1)) /= '"' ) .OR. &
                    ( LINE(WSTOP (IWORD+1):WSTOP (IWORD+1)) /= '"' ) ) THEN
                  CALL PARSE_ERROR ( "GET_TAGGED_STRING", "Invalid string constant." )

               ! . Decode the string.
               ELSE
                  STRING = LINE(WSTART(IWORD+1)+1:WSTOP(IWORD+1)-1)
               END IF

            ELSE
               CALL PARSE_ERROR ( "GET_TAGGED_STRING", "Missing string." )
            END IF

            ! . Flag the two words as read and return.
            QREAD(IWORD:IWORD+1) = .TRUE. ; RETURN

         END IF
      END IF
   END DO

   END SUBROUTINE GET_TAGGED_STRING

   !------------------
   SUBROUTINE GET_WORD
   !------------------

   ! . Local scalars.
   INTEGER :: IWORD, START, STOP

   ! . Initialize WORD and WRDLEN.
   WORD   = " "
   WRDLEN = 0

   ! . Loop over the words to find the next available word.
   DO IWORD = (LININD(CURLIN)+1),LININD(CURLIN+1)

      ! . A word exists that has not been read.
      IF ( .NOT. QREAD(IWORD) ) THEN

         ! . Set the word flag.
         QREAD(IWORD) = .TRUE.

         ! . Get the starting and stopping points for the word and its length.
         START  = WSTART(IWORD)
         STOP   = WSTOP(IWORD)
         WRDLEN = STOP - START + 1

         ! . Fill the word variable.
         WORD(1:WRDLEN) = LINE(START:STOP)

         ! . Exit.
         EXIT

      END IF
   END DO

   END SUBROUTINE GET_WORD

   !------------------------------------------
   SUBROUTINE PARSE_ERROR ( ROUTINE, MESSAGE )
   !------------------------------------------

   ! . Scalar argument declarations.
   CHARACTER ( LEN = * ), INTENT(IN) :: MESSAGE, ROUTINE

   ! . Local scalars.
   INTEGER :: START, STOP

   ! . Get the indices for the line.
   START = WSTART(LININD(CURLIN)+1)
   STOP  = WSTOP(LININD(CURLIN+1))

   ! . Write out the first word of the current command line.
   WRITE ( OUTPUT, "(/A)" ) "The Current Command Line: "//LINE(START:STOP)

   ! . Call the standard error routine.
   CALL ERROR ( ROUTINE, MESSAGE )

   END SUBROUTINE PARSE_ERROR

   !------------------
   SUBROUTINE POP_UNIT
   !------------------

   ! . Check for unparsed characters.
   CALL CHECK_UNPARSED_CHARACTERS

   ! . Pop the old unit from the stack.
   IF ( NUNITS > 1 ) THEN
      UNIT = UNITS(NUNITS-1)
   ! . Reinitialize the unit number.
   ELSE
      UNIT = -1
   END IF

   ! . Decrement NUNITS.
   NUNITS = MAX ( NUNITS - 1, 0 )

   END SUBROUTINE POP_UNIT

   !------------------------------
   SUBROUTINE PUSH_UNIT ( STREAM )
   !------------------------------

   ! . Scalar argument declarations.
   INTEGER, INTENT(IN) :: STREAM

   ! . Save the old unit.
   IF ( NUNITS > 0 ) UNITS(NUNITS) = UNIT

   ! . Increment the unit number.
   NUNITS = NUNITS + 1

   ! . Assign the unit number.
   UNIT = STREAM

   END SUBROUTINE PUSH_UNIT

END MODULE PARSING
