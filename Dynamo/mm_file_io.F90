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
!                    The Molecular Mechanics File I/O Module
!===============================================================================
!
! . Subroutines (public):
!
!   MM_FILE_PROCESS                Process a formatted MM definition file.
!   MM_FILE_READ                   Read in an unformatted MM file.
!   MM_FILE_WRITE                  Write out an unformatted MM file.
!
! . Subroutines (private):
!
!   READ_LINK_TYPE                 Read a link type.
!   READ_PARAMETER_DATA            Read the parameter data.
!   READ_RESIDUE_TYPE              Read a residue type.
!   READ_TYPE_DATA                 Read the type data.
!   READ_VARIANT_TYPE              Read a variant type.
!
!   WRITE_LINK_TYPE                Write a link type.
!   WRITE_PARAMETER_DATA           Write the parameter data.
!   WRITE_RESIDUE_TYPE             Write a residue type.
!   WRITE_TYPE_DATA                Write the type data.
!   WRITE_VARIANT_TYPE             Write a variant type.
!
!===============================================================================
MODULE MM_FILE_IO

! . Module declarations.
USE CONSTANTS,   ONLY : KCAL_TO_KJ, TO_RADIANS
USE DEFINITIONS, ONLY : DP, FORCE_FIELD, MAX_RECORD_LENGTH, VERSION
USE FILES
USE IO_UNITS
USE PARSING
USE STATUS,      ONLY : ERROR

USE MM_FILE_DATA
USE SEQUENCE,    ONLY : RESIDUE_NAME_LENGTH

IMPLICIT NONE
PRIVATE
PUBLIC :: MM_FILE_PROCESS, MM_FILE_READ, MM_FILE_WRITE
SAVE

!===============================================================================
CONTAINS
!===============================================================================

   !-----------------------------------------------
   SUBROUTINE MM_FILE_PROCESS ( FILE_OUT, FILE_IN )
   !-----------------------------------------------

   ! . Essential scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE_OUT

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: FILE_IN

   ! . Local scalars.
   INTEGER            :: IOSTAT, TUNIT, UNIT
   REAL ( KIND = DP ) :: VERNUM

   ! . MM file section flags.
   LOGICAL :: QELECT, QLJ, QLINKS, QPARAM, QRESID, QUNITS, QVARIA

   ! . Parameter section flags.
   LOGICAL :: QANGL, QBOND, QDIHE, QIMPR

   ! . Units option flags.
   LOGICAL :: QKCAL

   !-------------------------------------------------------------------------
   ! . Get the unit number for the MM file.
   !-------------------------------------------------------------------------
   ! . Check to see if the FILE argument is present.
   IF ( PRESENT ( FILE_IN ) ) THEN

      ! . Get the next unit number.
      UNIT = NEXT_UNIT ( )

      ! . Open the file.
      OPEN ( UNIT, FILE = FILE_IN, ACTION = "READ", STATUS = "OLD", IOSTAT = IOSTAT )

      ! . Check for an error.
      IF ( IOSTAT /= 0 ) CALL ERROR ( "MM_FILE_PROCESS", "Input file I/O error.", IOSTAT )
   
   ! . The unit for reading is the input stream.
   ELSE
      UNIT = INPUT
   END IF

   ! . Set the parsing unit.
   CALL PUSH_UNIT ( UNIT )

   ! . Open a temporary scratch file.
   TUNIT = NEXT_UNIT ( )
   OPEN ( TUNIT, FORM = "UNFORMATTED", RECL = MAX_RECORD_LENGTH, STATUS = "SCRATCH", IOSTAT = IOSTAT )

   ! . Check for an error.
   IF ( IOSTAT /= 0 ) CALL ERROR ( "MM_FILE_PROCESS", "Scratch file I/O error.", IOSTAT )

   !-------------------------------------------------------------------------
   ! . Read the data from the file.
   !-------------------------------------------------------------------------
   ! . Initialize the keyword flags.
   QELECT = .FALSE.
   QLJ    = .FALSE.
   QLINKS = .FALSE.
   QPARAM = .FALSE.
   QRESID = .FALSE.
   QUNITS = .FALSE.
   QVARIA = .FALSE.

   ! . Initialize the PARAMETER keyword flags.
   QANGL = .FALSE.
   QBOND = .FALSE.
   QDIHE = .FALSE.
   QIMPR = .FALSE.  

   ! . Initialize the UNITS option flags.
   QKCAL = .FALSE.

   ! . Check the file header.
   CALL GET_LINE
   CALL GET_WORD
   IF ( WORD(1:WRDLEN) /= "MM_DEFINITIONS" ) THEN
      CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Invalid file header." )
   END IF

   ! . Get the force field.
   CALL GET_WORD
   IF ( WORD(1:WRDLEN) /= FORCE_FIELD ) THEN
      CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Invalid FORCE FIELD." )
   END IF

   ! . Get the version number of the program.
   CALL GET_REAL ( VERNUM )
   IF ( VERNUM > VERSION ) THEN
      CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Invalid file version number." )
   END IF

   ! . Read in the type data.
   CALL GET_LINE
   CALL GET_WORD
   IF ( WORD(1:WRDLEN) == "TYPES" ) THEN
      CALL READ_TYPES
   ELSE
      CALL PARSE_ERROR ( "MM_FILE_PROCESS", "TYPES section not present." )
   END IF

   ! . Top of the loop for reading.
   10 CONTINUE

   ! . Get the next line in the file.
   CALL GET_LINE

   ! . Get the next word.
   CALL GET_WORD

   ! . Check the word.
   SELECT CASE ( WORD(1:WRDLEN) )
   CASE ( "ELECTROSTATICS" ) ; CALL READ_ELECTROSTATICS
   CASE ( "END"            ) ; GO TO 20
   CASE ( "LENNARD_JONES"  ) ; CALL READ_LENNARD_JONES
   CASE ( "LINKS"          ) ; CALL READ_LINKS
   CASE ( "PARAMETERS"     ) ; CALL READ_PARAMETERS
   CASE ( "RESIDUES"       ) ; CALL READ_RESIDUES
   CASE ( "UNITS"          ) ; CALL READ_UNITS
   CASE ( "VARIANTS"       ) ; CALL READ_VARIANTS
   CASE DEFAULT ; CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Unrecognized section keyword." )
   END SELECT

   ! . Read some more data.
   GO TO 10

   ! . There is no more data to read.
   20 CONTINUE

   ! . Check the section logical flags.
   IF ( .NOT. QRESID ) CALL ERROR ( "MM_FILE_PROCESS", "RESIDUES section missing." )

   ! . Check that the residue and variant data are consistent.
   CALL CHECK_RESIDUE_AND_VARIANT_DATA

   ! . Process the units commands.
   IF ( QKCAL ) CALL CHANGE_KCAL_TO_KJ

   !-------------------------------------------------------------------------
   ! . Finish up.
   !-------------------------------------------------------------------------
   ! . Close the files.
   CLOSE (  UNIT )
   CLOSE ( TUNIT )

   ! . Write out some information.
   IF ( UNIT == INPUT ) THEN
      WRITE ( OUTPUT, "(/A)" ) "MM file data read from input stream."
   ELSE
      WRITE ( OUTPUT, "(/A)" ) "MM file data read from " // FILE_IN
   END IF

   ! . Write out a summary of the MM file data.
   CALL MM_FILE_SUMMARY

   ! . Write out the MM file data to a file.
   CALL MM_FILE_WRITE ( FILE_OUT )

   ! . Free the MM file data structure.
   CALL MM_FILE_INITIALIZE

   ! . Reset the parsing unit.
   CALL POP_UNIT

   !============================================================================
   CONTAINS
   !============================================================================

      !---------------------------
      SUBROUTINE CHANGE_KCAL_TO_KJ
      !---------------------------

      ! . Change the force constant data.
      ANGLPS(1:NANGLPS)%FC = KCAL_TO_KJ * ANGLPS(1:NANGLPS)%FC
      BONDPS(1:NBONDPS)%FC = KCAL_TO_KJ * BONDPS(1:NBONDPS)%FC
      DIHEPS(1:NDIHEPS)%V0 = KCAL_TO_KJ * DIHEPS(1:NDIHEPS)%V0
      DIHEPS(1:NDIHEPS)%V1 = KCAL_TO_KJ * DIHEPS(1:NDIHEPS)%V1
      DIHEPS(1:NDIHEPS)%V2 = KCAL_TO_KJ * DIHEPS(1:NDIHEPS)%V2
      DIHEPS(1:NDIHEPS)%V3 = KCAL_TO_KJ * DIHEPS(1:NDIHEPS)%V3
      IMPRPS(1:NIMPRPS)%V0 = KCAL_TO_KJ * IMPRPS(1:NIMPRPS)%V0
      IMPRPS(1:NIMPRPS)%V1 = KCAL_TO_KJ * IMPRPS(1:NIMPRPS)%V1
      IMPRPS(1:NIMPRPS)%V2 = KCAL_TO_KJ * IMPRPS(1:NIMPRPS)%V2
      IMPRPS(1:NIMPRPS)%V3 = KCAL_TO_KJ * IMPRPS(1:NIMPRPS)%V3

      ! . Change the type data.
      TYPES(1:NTYPES)%EPSILON = KCAL_TO_KJ * TYPES(1:NTYPES)%EPSILON

      END SUBROUTINE CHANGE_KCAL_TO_KJ

      !----------------------------------------
      SUBROUTINE CHECK_RESIDUE_AND_VARIANT_DATA
      !----------------------------------------

      ! . Local scalars.
      INTEGER :: R, V

      ! . There are no variants.
      IF ( NVARIANTS <= 0 ) RETURN

      ! . Loop over the variants.
      DO V = 1,NVARIANTS

         ! . Check for a specific variant.
         IF ( VARIANTS(V)%RESIDUE_NAME /= " " ) THEN

            ! . Loop over the residues.
            DO R = 1,NRESIDUES
               IF ( VARIANTS(V)%RESIDUE_NAME == RESIDUES(R)%NAME ) GO TO 10
            END DO

            ! . There is an error.
            CALL ERROR ( "MM_FILE_PROCESS", "Unknown residue name for VARIANT declaration." )

            ! . Successful exit.
            10 CONTINUE

         END IF
      END DO

      END SUBROUTINE CHECK_RESIDUE_AND_VARIANT_DATA

      !---------------------------
      FUNCTION MATCH_TYPE ( NAME )
      !---------------------------

      ! . Function declarations.
      INTEGER :: MATCH_TYPE

      ! . Scalar argument declarations.
      CHARACTER ( LEN = * ), INTENT(IN) :: NAME

      ! . Local scalars.
      INTEGER :: I

      ! . Initialize MATCH_TYPE.
      MATCH_TYPE = -1

      ! . Check the length of the name.
      IF ( LEN ( NAME ) > TYPE_NAME_LENGTH ) RETURN

      ! . Loop over the types.
      DO I = 1,NTYPES
         IF ( NAME == TYPES(I)%NAME ) THEN
            MATCH_TYPE = I ; RETURN
         END IF
      END DO

      ! . Check for type X.
      IF ( NAME == DUMMY_TYPE ) MATCH_TYPE = DUMMY_CODE

      END FUNCTION MATCH_TYPE

      !---------------------
      SUBROUTINE READ_ANGLES
      !---------------------

      ! . Local scalars.
      INTEGER            :: I, J, K, P, Q
      REAL ( KIND = DP ) :: EQ, FC

      ! . Check and set the angles flag.
      IF ( QANGL) CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate ANGLES keyword." )
      QANGL = .TRUE.

      ! . Initialize the number of angle parameters.
      NANGLPS = 0

      ! . Rewind the scratch file.
      REWIND ( TUNIT )

      ! . Top of the loop for reading.
      10 CONTINUE

      ! . Get the next line.
      CALL GET_LINE
      CALL GET_WORD

      ! . Read in the data for the next type.
      IF ( WORD(1:WRDLEN) /= "END" ) THEN

         ! . Get the type of the first atom.
         I = MATCH_TYPE ( WORD(1:WRDLEN) )

         ! . Get the types of the second and third atoms.
         CALL GET_WORD
         J = MATCH_TYPE ( WORD(1:WRDLEN) )
         CALL GET_WORD
         K = MATCH_TYPE ( WORD(1:WRDLEN) )

         ! . Check for an error.
         IF ( ( I <= 0 ) .OR. ( J <= 0 ) .OR. ( K <= 0 ) ) THEN
            CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Invalid angle parameter types." )
         END IF

         ! . Read in the remaining data.
         CALL GET_REAL ( FC )
         CALL GET_REAL ( EQ )

         ! . Convert the angle to radians.
         EQ = TO_RADIANS * EQ

         ! . Increment the number of angle parameters.
         NANGLPS = NANGLPS + 1

         ! . Write out the angle data.
         WRITE ( TUNIT ) MIN ( I, K ), J, MAX ( I, K ), FC, EQ

         ! . Read the next line.
         GO TO 10

      END IF

      ! . Allocate space for the angle data.
      ALLOCATE ( ANGLPS(1:NANGLPS) )

      ! . Fill the ANGLPS array using the data on the file.
      REWIND ( TUNIT )
      DO P = 1,NANGLPS

         ! . Read in the data.
         READ ( TUNIT ) ANGLPS(P)%TYPE1, ANGLPS(P)%TYPE2, ANGLPS(P)%TYPE3, ANGLPS(P)%FC, ANGLPS(P)%EQ

         ! . Check for a duplicate angle parameter.
         DO Q = 1,(P-1)
            IF ( ( ANGLPS(P)%TYPE1 == ANGLPS(Q)%TYPE1 ) .AND. &
                 ( ANGLPS(P)%TYPE2 == ANGLPS(Q)%TYPE2 ) .AND. &
                 ( ANGLPS(P)%TYPE3 == ANGLPS(Q)%TYPE3 ) ) THEN
               CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate angle parameter." )
            END IF
         END DO

      END DO

      END SUBROUTINE READ_ANGLES

      !--------------------
      SUBROUTINE READ_BONDS
      !--------------------

      ! . Local scalars.
      INTEGER            :: I, J, P, Q
      REAL ( KIND = DP ) :: EQ, FC

      ! . Check and set the bonds flag.
      IF ( QBOND ) CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate BONDS keyword." )
      QBOND = .TRUE.

      ! . Initialize the number of bond parameters.
      NBONDPS = 0

      ! . Rewind the scratch file.
      REWIND ( TUNIT )

      ! . Top of the loop for reading.
      10 CONTINUE

      ! . Get the next line.
      CALL GET_LINE
      CALL GET_WORD

      ! . Read in the data for the next type.
      IF ( WORD(1:WRDLEN) /= "END" ) THEN

         ! . Find the type of the first atom.
         I = MATCH_TYPE ( WORD(1:WRDLEN) )

         ! . Find the type of the second atom.
         CALL GET_WORD
         J = MATCH_TYPE ( WORD(1:WRDLEN) )

         ! . Check the types.
         IF ( ( I <= 0 ) .OR. ( J <= 0 ) ) THEN
            CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Invalid bond parameter type." )
         END IF

         ! . Read in the remaining data.
         CALL GET_REAL ( FC )
         CALL GET_REAL ( EQ )

         ! . Increment the number of bond parameters.
         NBONDPS = NBONDPS + 1

         ! . Write out the bond data.
         WRITE ( TUNIT ) MIN ( I, J ), MAX ( I, J ), FC, EQ

         ! . Read the next line.
         GO TO 10

      END IF

      ! . Allocate space for the bond data.
      ALLOCATE ( BONDPS(1:NBONDPS) )

      ! . Fill the BONDPS array using the data on the file.
      REWIND ( TUNIT )
      DO P = 1,NBONDPS

         ! . Read in the data.
         READ ( TUNIT ) BONDPS(P)%TYPE1, BONDPS(P)%TYPE2, BONDPS(P)%FC, BONDPS(P)%EQ

         ! . Check for a duplicate bond parameter.
         DO Q = 1,(P-1)
            IF ( ( BONDPS(P)%TYPE1 == BONDPS(Q)%TYPE1 ) .AND. &
                 ( BONDPS(P)%TYPE2 == BONDPS(Q)%TYPE2 ) ) THEN
               CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate bond parameter." )
            END IF
         END DO

      END DO

      END SUBROUTINE READ_BONDS

      !------------------------
      SUBROUTINE READ_DIHEDRALS
      !------------------------

      ! . Local scalars.
      INTEGER            :: I, J, K, L, P, Q, TMP
      REAL ( KIND = DP ) :: V0, V1, V2, V3

      ! . Check and set the dihedrals flag.
      IF ( QDIHE ) CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate DIHEDRALS keyword." )
      QDIHE = .TRUE.

      ! . Initialize the number of dihedral parameters.
      NDIHEPS = 0

      ! . Rewind the scratch file.
      REWIND ( TUNIT )

      ! . Top of the loop for reading.
      10 CONTINUE

      ! . Get the next line.
      CALL GET_LINE
      CALL GET_WORD

      ! . Read in the data for the next type.
      IF ( WORD(1:WRDLEN) /= "END" ) THEN

         ! . Get the type of the first atom.
         I = MATCH_TYPE ( WORD(1:WRDLEN) )

         ! . Get the types of the remaining atoms.
         CALL GET_WORD
         J = MATCH_TYPE ( WORD(1:WRDLEN) )
         CALL GET_WORD
         K = MATCH_TYPE ( WORD(1:WRDLEN) )
         CALL GET_WORD
         L = MATCH_TYPE ( WORD(1:WRDLEN) )

         ! . Check for an error.
         IF ( ( I < 0 ) .OR. ( J <=0 ) .OR. ( K <=0 ) .OR. ( L < 0 ) ) THEN
            CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Invalid dihedral parameter types." )
         END IF

         ! . Reorder the data.
         IF ( J > K ) THEN
            TMP = J ; J = K ; K = TMP
            TMP = I ; I = L ; L = TMP
         ELSE IF ( J == K ) THEN
            IF ( I > L ) THEN
               TMP = I ; I = L ; L = TMP
            END IF
         END IF

         ! . Read in the remaining data.
         CALL GET_REAL ( V0 )
         CALL GET_REAL ( V1 )
         CALL GET_REAL ( V2 )
         CALL GET_REAL ( V3 )

         ! . Scale all the energy terms by a half.
         V0 = 0.5_DP * V0 ; V1 = 0.5_DP * V1 ; V2 = 0.5_DP * V2 ; V3 = 0.5_DP * V3

         ! . Increment the number of dihedral parameters.
         NDIHEPS = NDIHEPS + 1

         ! . Write out the dihedral data.
         WRITE ( TUNIT ) I, J, K, L, V0, V1, V2, V3

         ! . Read the next line.
         GO TO 10

      END IF

      ! . Allocate space for the dihedral data.
      ALLOCATE ( DIHEPS(1:NDIHEPS) )

      ! . Fill the DIHEPS array using the data on the file.
      REWIND ( TUNIT )
      DO P = 1,NDIHEPS

         ! . Read in the data.
         READ ( TUNIT ) DIHEPS(P)%TYPE1, DIHEPS(P)%TYPE2, DIHEPS(P)%TYPE3, DIHEPS(P)%TYPE4, &
                        DIHEPS(P)%V0,    DIHEPS(P)%V1,    DIHEPS(P)%V2,    DIHEPS(P)%V3

         ! . Check for a duplicate dihedral parameter.
         DO Q = 1,(P-1)
            IF ( ( DIHEPS(P)%TYPE1 == DIHEPS(Q)%TYPE1 ) .AND. &
                 ( DIHEPS(P)%TYPE2 == DIHEPS(Q)%TYPE2 ) .AND. &
                 ( DIHEPS(P)%TYPE3 == DIHEPS(Q)%TYPE3 ) .AND. &
                 ( DIHEPS(P)%TYPE4 == DIHEPS(Q)%TYPE4 ) ) THEN
               CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate dihedral parameter." )
            END IF
         END DO
      END DO

      END SUBROUTINE READ_DIHEDRALS

      !-----------------------------
      SUBROUTINE READ_ELECTROSTATICS
      !-----------------------------

      ! . Check and set the electrostatics flag.
      IF ( QELECT ) CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate ELECTROSTATICS keyword." )
      QELECT = .TRUE.

      ! . Read the next word.
      CALL GET_WORD

      ! . Check the word.
      SELECT CASE ( WORD(1:WRDLEN) )
      CASE ( "SCALE" ) ; CALL GET_REAL ( REF_SCALE_EL )
      CASE DEFAULT ; CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Unrecognized ELECTROSTATICS keyword." )
      END SELECT

      END SUBROUTINE READ_ELECTROSTATICS

      !------------------------
      SUBROUTINE READ_IMPROPERS
      !------------------------

      ! . Local scalars.
      INTEGER            :: I, J, K, L, P, Q
      REAL ( KIND = DP ) :: V0, V1, V2, V3

      ! . Check and set the impropers flag.
      IF ( QIMPR ) CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate IMPROPERS keyword." )
      QIMPR = .TRUE.

      ! . Initialize the number of improper parameters.
      NIMPRPS = 0

      ! . Rewind the scratch file.
      REWIND ( TUNIT )

      ! . Top of the loop for reading.
      10 CONTINUE

      ! . Get the next line.
      CALL GET_LINE
      CALL GET_WORD

      ! . Read in the data for the next type.
      IF ( WORD(1:WRDLEN) /= "END" ) THEN

         ! . Get the type of the first atom.
         I = MATCH_TYPE ( WORD(1:WRDLEN) )

         ! . Get the types of the remaining atoms.
         CALL GET_WORD
         J = MATCH_TYPE ( WORD(1:WRDLEN) )
         CALL GET_WORD
         K = MATCH_TYPE ( WORD(1:WRDLEN) )
         CALL GET_WORD
         L = MATCH_TYPE ( WORD(1:WRDLEN) )

         ! . Check for an error.
         IF ( ( I < 0 ) .OR. ( J < 0 ) .OR. ( K <= 0 ) .OR. ( L <= 0 ) ) THEN
            CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Invalid improper parameter types." )
         END IF

         ! . Read in the remaining data.
         CALL GET_REAL ( V0 )
         CALL GET_REAL ( V1 )
         CALL GET_REAL ( V2 )
         CALL GET_REAL ( V3 )

         ! . Scale all the energy terms by a half.
         V0 = 0.5_DP * V0 ; V1 = 0.5_DP * V1 ; V2 = 0.5_DP * V2 ; V3 = 0.5_DP * V3

         ! . Increment the number of dihedral parameters.
         NIMPRPS = NIMPRPS + 1

         ! . Write out the improper data.
         WRITE ( TUNIT ) I, J, K, L, V0, V1, V2, V3

         ! . Read the next line.
         GO TO 10

      END IF

      ! . Allocate space for the dihedral data.
      ALLOCATE ( IMPRPS(1:NIMPRPS) )

      ! . Fill the IMPRPS array using the data on the file.
      REWIND ( TUNIT )
      DO P = 1,NIMPRPS

         ! . Read in the data.
         READ ( TUNIT ) IMPRPS(P)%TYPE1, IMPRPS(P)%TYPE2, IMPRPS(P)%TYPE3, IMPRPS(P)%TYPE4, &
                        IMPRPS(P)%V0,    IMPRPS(P)%V1,    IMPRPS(P)%V2,    IMPRPS(P)%V3

         ! . Check for a duplicate dihedral parameter.
         DO Q = 1,(P-1)
            IF ( ( IMPRPS(P)%TYPE1 == IMPRPS(Q)%TYPE1 ) .AND. &
                 ( IMPRPS(P)%TYPE2 == IMPRPS(Q)%TYPE2 ) .AND. &
                 ( IMPRPS(P)%TYPE3 == IMPRPS(Q)%TYPE3 ) .AND. &
                 ( IMPRPS(P)%TYPE4 == IMPRPS(Q)%TYPE4 ) ) THEN
               CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate improper parameter." )
            END IF
         END DO
      END DO

      END SUBROUTINE READ_IMPROPERS

      !----------------------------
      SUBROUTINE READ_LENNARD_JONES
      !----------------------------

      ! . Check and set the Lennard-Jones flag.
      IF ( QLJ ) CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate LENNARD_JONES keyword." )
      QLJ = .TRUE.

      ! . Read the next word.
      CALL GET_WORD

      ! . Check the word.
      SELECT CASE ( WORD(1:WRDLEN) )
      CASE ( "SCALE" ) ; CALL GET_REAL ( REF_SCALE_LJ )
      CASE DEFAULT ; CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Unrecognized LENNARD_JONES keyword." )
      END SELECT

      END SUBROUTINE READ_LENNARD_JONES

      !--------------------
      SUBROUTINE READ_LINKS
      !--------------------

      ! . Local scalars.
      INTEGER         :: I, J
      TYPE(LINK_TYPE) :: LTEMP

      ! . Check and set the variants flag.
      IF ( QLINKS ) CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate LINKS keyword." )
      QLINKS = .TRUE.

      ! . Rewind the scratch file.
      REWIND ( TUNIT )

      ! . Initialize the number of links.
      NRESLINK = 0

      ! . Get the first line after the section header.
      CALL GET_LINE
      CALL GET_WORD
      IF ( WORD(1:WRDLEN) /= "LINK" ) CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Missing LINK keyword." )

      ! . Top of the loop for reading link information.
      10 CONTINUE

      ! . Get the link name.
      CALL GET_WORD
      IF ( WRDLEN > RESIDUE_NAME_LENGTH ) THEN
         CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Link name too long." )
      ELSE
         LTEMP%NAME = WORD(1:WRDLEN)
      END IF

      ! . Read the variant information.
      CALL READ_SINGLE_VARIANT ( LTEMP%RESIDUE1 )
      CALL READ_SINGLE_VARIANT ( LTEMP%RESIDUE2 )

      ! . Get the next word.
      CALL GET_LINE
      CALL GET_WORD

      ! . Write out the residue information to the temporary file.
      IF ( ( WORD(1:WRDLEN) == "END" ) .OR. ( WORD(1:WRDLEN) == "LINK" ) ) THEN

         ! . Increment the number of links.
         NRESLINK = NRESLINK + 1

         ! . Write out the link data.
         CALL WRITE_LINK_TYPE ( TUNIT, LTEMP )

         ! . Initialize LTEMP.
         CALL INITIALIZE_LINK_TYPE ( LTEMP )

         ! . Read in further link data.
         IF ( WORD(1:WRDLEN) == "LINK" ) GO TO 10

      ELSE
         CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Missing LINK or END keyword." )
      END IF

      ! . Rewind the scratch file.
      REWIND ( TUNIT )

      ! . Allocate space for the link information.
      ALLOCATE ( RESLINK(1:NRESLINK) )

      ! . Loop over the links.
      DO I = 1,NRESLINK

         ! . Read in the link data.
         CALL READ_LINK_TYPE ( TUNIT, RESLINK(I) )

         ! . Check for duplicate names.
         DO J = 1,(I-1)
            IF ( RESLINK(I)%NAME == RESLINK(J)%NAME ) THEN
               CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate LINK names." )
            END IF
         END DO

      END DO

      END SUBROUTINE READ_LINKS

      !-------------------------
      SUBROUTINE READ_PARAMETERS
      !-------------------------

      ! . Check and set the parameters flag.
      IF ( QPARAM ) CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate PARAMETERS keyword." )
      QPARAM = .TRUE.

      ! . Top of the loop for reading.
      10 CONTINUE

      ! . Get the next line in the file.
      CALL GET_LINE

      ! . Get the next word.
      CALL GET_WORD

      ! . Check the word.
      SELECT CASE ( WORD(1:WRDLEN) )
      CASE ( "ANGLES"    ) ; CALL READ_ANGLES
      CASE ( "BONDS"     ) ; CALL READ_BONDS
      CASE ( "END"       ) ; RETURN
      CASE ( "DIHEDRALS" ) ; CALL READ_DIHEDRALS
      CASE ( "IMPROPERS" ) ; CALL READ_IMPROPERS
      CASE DEFAULT ; CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Unrecognized PARAMETER keyword." )
      END SELECT

      ! . Read some more data.
      GO TO 10

      END SUBROUTINE READ_PARAMETERS

      !-----------------------
      SUBROUTINE READ_RESIDUES
      !-----------------------

      ! . Local scalars.
      INTEGER            :: I, J, P, TYPE
      TYPE(RESIDUE_TYPE) :: RTEMP

      ! . Check and set the residues flag.
      IF ( QRESID ) CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate RESIDUES keyword." )
      QRESID = .TRUE.

      ! . Rewind the scratch file.
      REWIND ( TUNIT )

      ! . Initialize the number of residues.
      NRESIDUES = 0

      ! . Get the first line after the section header.
      CALL GET_LINE
      CALL GET_WORD
      IF ( WORD(1:WRDLEN) /= "RESIDUE" ) CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Missing RESIDUE keyword." )

      ! . Top of the loop for reading residue information.
      10 CONTINUE

      ! . Get the residue name.
      CALL GET_WORD
      IF ( WRDLEN > RESIDUE_NAME_LENGTH ) THEN
         CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Residue name too long." )
      ELSE
         RTEMP%NAME = WORD(1:WRDLEN)
      END IF

      ! . Get the number of atoms, bonds and impropers for the residue.
      CALL GET_LINE
      CALL GET_INTEGER ( RTEMP%MM_NATOMS     )
      CALL GET_INTEGER ( RTEMP%NBONDS     )
      CALL GET_INTEGER ( RTEMP%NIMPROPERS )

      ! . Allocate space for the RTEMP arrays.
      CALL ALLOCATE_RESIDUE_TYPE ( RTEMP )

      ! . Read in the atom type data.
      DO I = 1,RTEMP%MM_NATOMS

         ! . Read in the next line.
         CALL GET_LINE

         ! . Get the atom name.
         CALL GET_WORD
         IF ( WRDLEN > ATOM_NAME_LENGTH ) CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Atom name too long." )
         RTEMP%NAMES(I) = WORD(1:WRDLEN)

         ! . Check for a duplicate atom name.
         DO J = 1,(I-1)
            IF ( RTEMP%NAMES(I) == RTEMP%NAMES(J) ) THEN
               CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate atom name in residue." )
            END IF
         END DO

         ! . Get the atom type.
         CALL GET_WORD
         TYPE = MATCH_TYPE ( WORD(1:WRDLEN) )
         IF ( TYPE <= 0 ) CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Unknown type for atom in residue." )
         RTEMP%TYPES(I) = TYPE

         ! . Get the charge.
         CALL GET_REAL ( RTEMP%CHARGES(I) )

      END DO

      ! . Read in the bond data.
      DO P = 1,RTEMP%NBONDS

         ! . Read in and check the atom names.
         CALL GET_LINE
         DO I = 1,2
            CALL GET_WORD
            RTEMP%BONDS(I,P) = WORD(1:WRDLEN)
            IF ( .NOT. ANY ( WORD(1:WRDLEN) == RTEMP%NAMES(1:RTEMP%MM_NATOMS) ) ) THEN
               IF ( ( WORD(1:WRDLEN) /= LINK_LATERAL ) .AND. &
                    ( WORD(1:WRDLEN) /= LINK_MINUS   ) .AND. &
                    ( WORD(1:WRDLEN) /= LINK_PLUS    ) ) THEN
                  CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Unknown atom name in bond." )
               END IF
            END IF
         END DO

         ! . Check for duplicate names.
         IF ( RTEMP%BONDS(1,P) == RTEMP%BONDS(2,P) ) THEN
            CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Atom names same in bond." )
         END IF

      END DO

      ! . Read in the improper data.
      DO P = 1,RTEMP%NIMPROPERS

         ! . Read in and check the atom names.
         CALL GET_LINE
         DO I = 1,4
            CALL GET_WORD
            RTEMP%IMPROPERS(I,P) = WORD(1:WRDLEN)
            IF ( .NOT. ANY ( WORD(1:WRDLEN) == RTEMP%NAMES(1:RTEMP%MM_NATOMS) ) ) THEN
               IF ( ( WORD(1:WRDLEN) /= LINK_LATERAL ) .AND. &
                    ( WORD(1:WRDLEN) /= LINK_MINUS   ) .AND. &
                    ( WORD(1:WRDLEN) /= LINK_PLUS    ) ) THEN
                  CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Unknown atom name in improper." )
               END IF
            END IF
         END DO

         ! . Check for duplicate names.
         IF ( ( RTEMP%IMPROPERS(1,P) == RTEMP%IMPROPERS(2,P) ) .OR. &
              ( RTEMP%IMPROPERS(1,P) == RTEMP%IMPROPERS(3,P) ) .OR. &
              ( RTEMP%IMPROPERS(1,P) == RTEMP%IMPROPERS(4,P) ) .OR. &
              ( RTEMP%IMPROPERS(2,P) == RTEMP%IMPROPERS(3,P) ) .OR. &
              ( RTEMP%IMPROPERS(2,P) == RTEMP%IMPROPERS(4,P) ) .OR. &
              ( RTEMP%IMPROPERS(3,P) == RTEMP%IMPROPERS(4,P) ) ) THEN
            CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Atom names same in improper." )
         END IF

      END DO

      ! . Get the next word.
      CALL GET_LINE
      CALL GET_WORD

      ! . Write out the residue information to the temporary file.
      IF ( ( WORD(1:WRDLEN) == "END" ) .OR. ( WORD(1:WRDLEN) == "RESIDUE" ) ) THEN

         ! . Increment the number of residues.
         NRESIDUES = NRESIDUES + 1

         ! . Write out the residue data.
         CALL WRITE_RESIDUE_TYPE ( TUNIT, RTEMP )

         ! . Initialize RTEMP.
         CALL INITIALIZE_RESIDUE_TYPE ( RTEMP )

         ! . Read in further residue data.
         IF ( WORD(1:WRDLEN) == "RESIDUE" ) GO TO 10

      ELSE
         CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Missing RESIDUE or END keyword." )
      END IF

      ! . Rewind the scratch file.
      REWIND ( TUNIT )

      ! . Allocate space for the residue information.
      ALLOCATE ( RESIDUES(1:NRESIDUES) )

      ! . Loop over the residues.
      DO I = 1,NRESIDUES

         ! . Read in the residue data.
         CALL READ_RESIDUE_TYPE ( TUNIT, RESIDUES(I) )

         ! . Check for duplicate names.
         DO J = 1,(I-1)
            IF ( RESIDUES(I)%NAME == RESIDUES(J)%NAME ) THEN
               CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate RESIDUE names." )
            END IF
         END DO

      END DO

      END SUBROUTINE READ_RESIDUES

      !---------------------------------------
      SUBROUTINE READ_SINGLE_VARIANT ( VTEMP )
      !---------------------------------------

      ! . Scalar argument declarations.
      TYPE(VARIANT_TYPE), INTENT(INOUT) :: VTEMP

      ! . Local scalars.
      INTEGER :: I, J, P, TYPE

      ! . Get the number of deletes, adds, types, charges, bonds and impropers.
      CALL GET_LINE
      CALL GET_INTEGER ( VTEMP%NDELETES   )
      CALL GET_INTEGER ( VTEMP%NADDS      )
      CALL GET_INTEGER ( VTEMP%NCHARGES   )
      CALL GET_INTEGER ( VTEMP%NBONDS     )
      CALL GET_INTEGER ( VTEMP%NIMPROPERS )

      ! . Allocate space for the VTEMP arrays.
      CALL ALLOCATE_VARIANT_TYPE ( VTEMP )

      ! . Read in the atom deletion data.
      DO I = 1,VTEMP%NDELETES

         ! . Read in the name of the atom to delete.
         CALL GET_LINE
         CALL GET_WORD
         VTEMP%DELATM(I) = WORD(1:WRDLEN)

         ! . Check for a duplicate atom name.
         DO J = 1,(I-1)
            IF ( VTEMP%DELATM(I) == VTEMP%DELATM(J) ) THEN
               CALL PARSE_ERROR ( "MM_FILE_PROCESS", "An atom has been deleted twice." )
            END IF
         END DO

      END DO

      ! . Read in the data for the atoms to add.
      DO I = 1,VTEMP%NADDS

         ! . Read in the name of the atom to add.
         CALL GET_LINE
         CALL GET_WORD
         VTEMP%ADDATM(I) = WORD(1:WRDLEN)

         ! . Check for a duplicate atom name.
         DO J = 1,(I-1)
            IF ( VTEMP%ADDATM(I) == VTEMP%ADDATM(J) ) THEN
               CALL PARSE_ERROR ( "MM_FILE_PROCESS", "An atom has been added twice." )
            END IF
         END DO

         ! . Read in the atom type.
         CALL GET_WORD
         TYPE = MATCH_TYPE ( WORD(1:WRDLEN) )
         IF ( TYPE > 0 ) THEN
            VTEMP%ADDTYP(I) = TYPE
         ELSE
            CALL PARSE_ERROR ( "MM_FILE_PROCESS", "The type of the atom to be added is unknown." )
         END IF

         ! . Get the charge.
         CALL GET_REAL ( VTEMP%ADDCHG(I) )

      END DO

      ! . Read in the data for the atoms whose charges are to be changed.
      DO I = 1,VTEMP%NCHARGES

         ! . Read in the name of the atom to add.
         CALL GET_LINE
         CALL GET_WORD
         VTEMP%CHGATM(I) = WORD(1:WRDLEN)

         ! . Check for a duplicate atom name.
         DO J = 1,(I-1)
            IF ( VTEMP%CHGATM(I) == VTEMP%CHGATM(J) ) THEN
               CALL PARSE_ERROR ( "MM_FILE_PROCESS", "An atom's charge has been changed twice." )
            END IF
         END DO

         ! . Get the charge.
         CALL GET_REAL ( VTEMP%CHGCHG(I) )

      END DO

      ! . Read in the bond data.
      DO P = 1,VTEMP%NBONDS

         ! . Read in and check the atom names.
         CALL GET_LINE
         DO I = 1,2
            CALL GET_WORD
            VTEMP%BONDS(I,P) = WORD(1:WRDLEN)
         END DO

         ! . Check for duplicate names.
         IF ( VTEMP%BONDS(1,P) == VTEMP%BONDS(2,P) ) THEN
            CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Atom names same in variant bond." )
         END IF

      END DO

      ! . Read in the improper data.
      DO P = 1,VTEMP%NIMPROPERS

         ! . Read in and check the atom names.
         CALL GET_LINE
         DO I = 1,4
            CALL GET_WORD
            VTEMP%IMPROPERS(I,P) = WORD(1:WRDLEN)
         END DO

         ! . Check for duplicate names.
         IF ( ( VTEMP%IMPROPERS(1,P) == VTEMP%IMPROPERS(2,P) ) .OR. &
              ( VTEMP%IMPROPERS(1,P) == VTEMP%IMPROPERS(3,P) ) .OR. &
              ( VTEMP%IMPROPERS(1,P) == VTEMP%IMPROPERS(4,P) ) .OR. &
              ( VTEMP%IMPROPERS(2,P) == VTEMP%IMPROPERS(3,P) ) .OR. &
              ( VTEMP%IMPROPERS(2,P) == VTEMP%IMPROPERS(4,P) ) .OR. &
              ( VTEMP%IMPROPERS(3,P) == VTEMP%IMPROPERS(4,P) ) ) THEN
            CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Atom names same in variant improper." )
         END IF

      END DO

      END SUBROUTINE READ_SINGLE_VARIANT

      !--------------------
      SUBROUTINE READ_TYPES
      !--------------------

      ! . Local scalars.
      CHARACTER ( LEN = TYPE_NAME_LENGTH ) :: TNAM
      INTEGER                              :: I, J, NUMBER
      REAL ( KIND = DP )                   :: TEPS, TSIG

      ! . Initialize the number of types.
      NTYPES = 0

      ! . Rewind the scratch file.
      REWIND ( TUNIT )

      ! . Top of the loop for reading.
      10 CONTINUE

      ! . Get the next line and its first word.
      CALL GET_LINE
      CALL GET_WORD

      ! . Read in the data for the next type.
      IF ( WORD(1:WRDLEN) /= "END" ) THEN

         ! . Check the word length.
         IF ( WRDLEN > TYPE_NAME_LENGTH ) THEN
            CALL PARSE_ERROR ( "MM_FILE_PROCESS", "TYPE name too long." )
         ELSE
            TNAM = WORD(1:WRDLEN)
         END IF

         ! . Increment the number of types.
         NTYPES = NTYPES + 1

         ! . Read in the remaining type data.
         CALL GET_INTEGER ( NUMBER )
         CALL GET_REAL    ( TSIG   )
         CALL GET_REAL    ( TEPS   )

         ! . Write out the type data.
         WRITE ( TUNIT ) TNAM, NUMBER, TSIG, TEPS

         ! . Read the next line.
         GO TO 10

      END IF

      ! . Allocate space for the type data.
      ALLOCATE ( TYPES(1:NTYPES) )

      ! . Fill the TYPES array using the data on the file.
      REWIND ( TUNIT )
      DO I = 1,NTYPES
         READ ( TUNIT ) TNAM, NUMBER, TSIG, TEPS
         TYPES(I)%NAME    = TNAM
         TYPES(I)%NUMBER  = NUMBER
         TYPES(I)%EPSILON = TEPS
         TYPES(I)%SIGMA   = TSIG
      END DO

      ! . Check the type names.
      DO I = 1,NTYPES

         ! . Get the type name.
         TNAM = TYPES(I)%NAME

         ! . Check for type X.
         IF ( TNAM == DUMMY_TYPE ) CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Type name `X' is reserved." ) 

         ! . Check for a duplication.
         DO J = 1,(I-1)
            IF ( TNAM == TYPES(J)%NAME ) THEN
               CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate TYPE name." )
            END IF 
         END DO
      END DO

      END SUBROUTINE READ_TYPES

      !--------------------
      SUBROUTINE READ_UNITS
      !--------------------

      ! . Check and set the electrostatics flag.
      IF ( QUNITS ) CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate UNITS keyword." )
      QUNITS = .TRUE.

      ! . Read the next word.
      CALL GET_WORD

      ! . Check the word.
      SELECT CASE ( WORD(1:WRDLEN) )
      CASE ( "KCAL/MOLE" ) ; QKCAL = .TRUE.
      CASE DEFAULT ; CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Unrecognized UNITS keyword." )
      END SELECT

      END SUBROUTINE READ_UNITS

      !-----------------------
      SUBROUTINE READ_VARIANTS
      !-----------------------

      ! . Local scalars.
      INTEGER            :: I, J
      TYPE(VARIANT_TYPE) :: VTEMP

      ! . Check and set the variants flag.
      IF ( QVARIA ) CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate VARIANTS keyword." )
      QVARIA = .TRUE.

      ! . Rewind the scratch file.
      REWIND ( TUNIT )

      ! . Initialize the number of variants.
      NVARIANTS = 0

      ! . Get the first line after the section header.
      CALL GET_LINE
      CALL GET_WORD
      IF ( WORD(1:WRDLEN) /= "VARIANT" ) CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Missing VARIANT keyword." )

      ! . Top of the loop for reading variant information.
      10 CONTINUE

      ! . Get the variant name.
      CALL GET_WORD
      IF ( WRDLEN > RESIDUE_NAME_LENGTH ) THEN
         CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Variant name too long." )
      ELSE
         VTEMP%NAME = WORD(1:WRDLEN)
      END IF

      ! . Initialize the residue name and exclusive flags.
      VTEMP%RESIDUE_NAME = " "

      ! . Get the next word.
      CALL GET_WORD
      IF ( WRDLEN > 0 ) VTEMP%RESIDUE_NAME = WORD(1:WRDLEN)

      ! . Read the variant information.
      CALL READ_SINGLE_VARIANT ( VTEMP )

      ! . Get the next word.
      CALL GET_LINE
      CALL GET_WORD

      ! . Write out the residue information to the temporary file.
      IF ( ( WORD(1:WRDLEN) == "END" ) .OR. ( WORD(1:WRDLEN) == "VARIANT" ) ) THEN

         ! . Increment the number of residues.
         NVARIANTS = NVARIANTS + 1

         ! . Write out the residue data.
         CALL WRITE_VARIANT_TYPE ( TUNIT, VTEMP )

         ! . Initialize VTEMP.
         CALL INITIALIZE_VARIANT_TYPE ( VTEMP )

         ! . Read in further variant data.
         IF ( WORD(1:WRDLEN) == "VARIANT" ) GO TO 10

      ELSE
         CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Missing VARIANT or END keyword." )
      END IF

      ! . Rewind the scratch file.
      REWIND ( TUNIT )

      ! . Allocate space for the variant information.
      ALLOCATE ( VARIANTS(1:NVARIANTS) )

      ! . Loop over the variants.
      DO I = 1,NVARIANTS

         ! . Read in the variant data.
         CALL READ_VARIANT_TYPE ( TUNIT, VARIANTS(I) )

         ! . Check for duplicate names.
         DO J = 1,(I-1)
            IF ( ( VARIANTS(I)%NAME         == VARIANTS(J)%NAME         ) .AND. &
                 ( VARIANTS(I)%RESIDUE_NAME == VARIANTS(J)%RESIDUE_NAME ) ) THEN
               CALL PARSE_ERROR ( "MM_FILE_PROCESS", "Duplicate VARIANT names." )
            END IF
         END DO

      END DO

      END SUBROUTINE READ_VARIANTS

   END SUBROUTINE MM_FILE_PROCESS

   !-------------------------------
   SUBROUTINE MM_FILE_READ ( FILE )
   !-------------------------------

   ! . Essential scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE

   ! . Local scalars.
   CHARACTER ( LEN = LEN ( FORCE_FIELD ) ) :: HDRSTR
   INTEGER                                 :: I, IOSTAT, UNIT
   REAL ( KIND = DP )                      :: VERNUM

   ! . Open the file.
   UNIT = NEXT_UNIT ( )
   OPEN ( UNIT, ACTION = "READ", FILE = FILE, FORM = "UNFORMATTED", RECL = MAX_RECORD_LENGTH, STATUS = "UNKNOWN", &
                IOSTAT = IOSTAT )

   ! . Check for an error.
   IF ( IOSTAT /= 0 ) CALL ERROR ( "MM_FILE_READ", "I/O Error.", IOSTAT )

   ! . Read in the header.
   READ ( UNIT ) HDRSTR, VERNUM
   IF ( ( HDRSTR /= FORCE_FIELD ) .OR. ( VERNUM > VERSION ) ) CALL ERROR ( "MM_FILE_READ", "Invalid MM file header." )

   ! . Read in the counters.
   READ ( UNIT ) NANGLPS, NBONDPS, NDIHEPS, NIMPRPS, NRESLINK, NRESIDUES, NTYPES, NVARIANTS

   ! . Create the data structure.
   CALL MM_FILE_ALLOCATE

   ! . Read in the electrostatics options.
   READ ( UNIT ) REF_SCALE_EL, REF_SCALE_LJ

   ! . Read in the parameter data.
   CALL READ_PARAMETER_DATA ( UNIT )

   ! . Read in the type data.
   CALL READ_TYPE_DATA ( UNIT )

   ! . Read in the link data.
   DO I = 1,NRESLINK
      CALL READ_LINK_TYPE ( UNIT, RESLINK(I) )
   END DO

   ! . Read in the residue data.
   DO I = 1,NRESIDUES
      CALL READ_RESIDUE_TYPE ( UNIT, RESIDUES(I) )
   END DO

   ! . Read in the variant data.
   DO I = 1,NVARIANTS
      CALL READ_VARIANT_TYPE ( UNIT, VARIANTS(I) )
   END DO

   ! . Close the file.
   CLOSE ( UNIT )

   ! . Write out some information.
   WRITE ( OUTPUT, "(/A)" ) "MM file data read in binary format from " // FILE

   END SUBROUTINE MM_FILE_READ

   !--------------------------------
   SUBROUTINE MM_FILE_WRITE ( FILE )
   !--------------------------------

   ! . Essential scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: FILE

   ! . Local scalars.
   INTEGER :: I, IOSTAT, UNIT

   ! . Open the file.
   UNIT = NEXT_UNIT ( )
   OPEN ( UNIT, ACTION = "WRITE", FILE = FILE, FORM = "UNFORMATTED", RECL = MAX_RECORD_LENGTH, STATUS = "UNKNOWN", &
                IOSTAT = IOSTAT )

   ! . Check for an error.
   IF ( IOSTAT /= 0 ) CALL ERROR ( "MM_FILE_WRITE", "I/O Error.", IOSTAT )

   ! . Write out the header.
   WRITE ( UNIT ) FORCE_FIELD, VERSION

   ! . Write out the counters.
   WRITE ( UNIT ) NANGLPS, NBONDPS, NDIHEPS, NIMPRPS, NRESLINK, NRESIDUES, NTYPES, NVARIANTS

   ! . Write out the electrostatics options.
   WRITE ( UNIT ) REF_SCALE_EL, REF_SCALE_LJ

   ! . Write out the parameter data.
   CALL WRITE_PARAMETER_DATA ( UNIT )

   ! . Write out the type data.
   CALL WRITE_TYPE_DATA ( UNIT )

   ! . Write out the link data.
   DO I = 1,NRESLINK
      CALL WRITE_LINK_TYPE ( UNIT, RESLINK(I) )
   END DO

   ! . Write out the residue data.
   DO I = 1,NRESIDUES
      CALL WRITE_RESIDUE_TYPE ( UNIT, RESIDUES(I) )
   END DO

   ! . Write out the variant data.
   DO I = 1,NVARIANTS
      CALL WRITE_VARIANT_TYPE ( UNIT, VARIANTS(I) )
   END DO

   ! . Close the file.
   CLOSE ( UNIT )

   ! . Write out some information.
   WRITE ( OUTPUT, "(/A)" ) "MM file data written in binary format to " // FILE

   END SUBROUTINE MM_FILE_WRITE

   !============================================================================
   ! . Other routines.
   !============================================================================

   !----------------------------------------
   SUBROUTINE READ_LINK_TYPE ( UNIT, LTYPE )
   !----------------------------------------
   
   ! . Scalar arguments.
   INTEGER,         INTENT(IN)  :: UNIT
   TYPE(LINK_TYPE), INTENT(OUT) :: LTYPE

   ! . Read in the link variables.
   READ ( UNIT ) LTYPE%NAME

   ! . Read in the link variants.
   CALL READ_VARIANT_TYPE ( UNIT, LTYPE%RESIDUE1 )
   CALL READ_VARIANT_TYPE ( UNIT, LTYPE%RESIDUE2 )

   END SUBROUTINE READ_LINK_TYPE

   !--------------------------------------
   SUBROUTINE READ_PARAMETER_DATA ( UNIT )
   !--------------------------------------
   
   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: UNIT

   ! . Read in the parameter counters.
   READ ( UNIT ) NANGLPS, NBONDPS, NDIHEPS, NIMPRPS

   ! . Read in the parameter arrays.
   IF ( NANGLPS > 0 ) READ ( UNIT ) ANGLPS
   IF ( NBONDPS > 0 ) READ ( UNIT ) BONDPS
   IF ( NDIHEPS > 0 ) READ ( UNIT ) DIHEPS
   IF ( NIMPRPS > 0 ) READ ( UNIT ) IMPRPS

   END SUBROUTINE READ_PARAMETER_DATA

   !-------------------------------------------
   SUBROUTINE READ_RESIDUE_TYPE ( UNIT, RTYPE )
   !-------------------------------------------
   
   ! . Scalar arguments.
   INTEGER,            INTENT(IN)  :: UNIT
   TYPE(RESIDUE_TYPE), INTENT(OUT) :: RTYPE

   ! . Local scalars.
   INTEGER :: I, J

   ! . Read in the residue scalars.
   READ ( UNIT ) RTYPE%NAME, RTYPE%MM_NATOMS, RTYPE%NBONDS, RTYPE%NIMPROPERS

   ! . Allocate space for the type.
   CALL ALLOCATE_RESIDUE_TYPE ( RTYPE )

   ! . Read in the residue arrays.
   IF ( RTYPE%MM_NATOMS     > 0 ) READ ( UNIT ) ( RTYPE%NAMES(I), RTYPE%TYPES(I), RTYPE%CHARGES(I), I = 1,RTYPE%MM_NATOMS )
   IF ( RTYPE%NBONDS     > 0 ) READ ( UNIT ) ( ( RTYPE%BONDS(J,I), J = 1,2 ), I = 1,RTYPE%NBONDS )
   IF ( RTYPE%NIMPROPERS > 0 ) READ ( UNIT ) ( ( RTYPE%IMPROPERS(J,I), J = 1,4 ), I = 1,RTYPE%NIMPROPERS )

   END SUBROUTINE READ_RESIDUE_TYPE

   !---------------------------------
   SUBROUTINE READ_TYPE_DATA ( UNIT )
   !---------------------------------
   
   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: UNIT

   ! . Read in the type counter.
   READ ( UNIT ) NTYPES

   ! . Read in the type array.
   IF ( NTYPES > 0 ) READ ( UNIT ) TYPES

   END SUBROUTINE READ_TYPE_DATA

   !-------------------------------------------
   SUBROUTINE READ_VARIANT_TYPE ( UNIT, VTYPE )
   !-------------------------------------------
   
   ! . Scalar arguments.
   INTEGER,            INTENT(IN)  :: UNIT
   TYPE(VARIANT_TYPE), INTENT(OUT) :: VTYPE

   ! . Local scalars.
   INTEGER :: I, J

   ! . Read in the variant scalars.
   READ ( UNIT ) VTYPE%NAME,     VTYPE%RESIDUE_NAME, VTYPE%NADDS,      VTYPE%NBONDS, &
                 VTYPE%NCHARGES, VTYPE%NDELETES,     VTYPE%NIMPROPERS

   ! . Allocate space for the type.
   CALL ALLOCATE_VARIANT_TYPE ( VTYPE )

   ! . Read in the variant arrays.
   IF ( VTYPE%NADDS      > 0 ) READ ( UNIT ) ( VTYPE%ADDATM(I), VTYPE%ADDCHG(I), VTYPE%ADDTYP(I), I = 1,VTYPE%NADDS )
   IF ( VTYPE%NBONDS     > 0 ) READ ( UNIT ) ( ( VTYPE%BONDS(J,I), J = 1,2 ), I = 1,VTYPE%NBONDS )
   IF ( VTYPE%NCHARGES   > 0 ) READ ( UNIT ) ( VTYPE%CHGATM(I), VTYPE%CHGCHG(I), I = 1,VTYPE%NCHARGES )
   IF ( VTYPE%NDELETES   > 0 ) READ ( UNIT ) ( VTYPE%DELATM(I), I = 1,VTYPE%NDELETES )
   IF ( VTYPE%NIMPROPERS > 0 ) READ ( UNIT ) ( ( VTYPE%IMPROPERS(J,I), J = 1,4 ), I = 1,VTYPE%NIMPROPERS )

   END SUBROUTINE READ_VARIANT_TYPE

   !-----------------------------------------
   SUBROUTINE WRITE_LINK_TYPE ( UNIT, LTYPE )
   !-----------------------------------------
   
   ! . Scalar arguments.
   INTEGER,         INTENT(IN) :: UNIT
   TYPE(LINK_TYPE), INTENT(IN) :: LTYPE

   ! . Write out the link variables.
   WRITE ( UNIT ) LTYPE%NAME

   ! . Write out the link variants.
   CALL WRITE_VARIANT_TYPE ( UNIT, LTYPE%RESIDUE1 )
   CALL WRITE_VARIANT_TYPE ( UNIT, LTYPE%RESIDUE2 )

   END SUBROUTINE WRITE_LINK_TYPE

   !---------------------------------------
   SUBROUTINE WRITE_PARAMETER_DATA ( UNIT )
   !---------------------------------------
   
   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: UNIT

   ! . Write out the parameter counters.
   WRITE ( UNIT ) NANGLPS, NBONDPS, NDIHEPS, NIMPRPS

   ! . Write out the parameter arrays.
   IF ( NANGLPS > 0 ) WRITE ( UNIT ) ANGLPS
   IF ( NBONDPS > 0 ) WRITE ( UNIT ) BONDPS
   IF ( NDIHEPS > 0 ) WRITE ( UNIT ) DIHEPS
   IF ( NIMPRPS > 0 ) WRITE ( UNIT ) IMPRPS

   END SUBROUTINE WRITE_PARAMETER_DATA

   !--------------------------------------------
   SUBROUTINE WRITE_RESIDUE_TYPE ( UNIT, RTYPE )
   !--------------------------------------------
   
   ! . Scalar arguments.
   INTEGER,            INTENT(IN) :: UNIT
   TYPE(RESIDUE_TYPE), INTENT(IN) :: RTYPE

   ! . Local scalars.
   INTEGER :: I, J

   ! . Write out the residue scalars.
   WRITE ( UNIT ) RTYPE%NAME, RTYPE%MM_NATOMS, RTYPE%NBONDS, RTYPE%NIMPROPERS

   ! . Write out the residue arrays.
   IF ( RTYPE%MM_NATOMS     > 0 ) WRITE ( UNIT ) ( RTYPE%NAMES(I), RTYPE%TYPES(I), RTYPE%CHARGES(I), I = 1,RTYPE%MM_NATOMS )
   IF ( RTYPE%NBONDS     > 0 ) WRITE ( UNIT ) ( ( RTYPE%BONDS(J,I), J = 1,2 ), I = 1,RTYPE%NBONDS )
   IF ( RTYPE%NIMPROPERS > 0 ) WRITE ( UNIT ) ( ( RTYPE%IMPROPERS(J,I), J = 1,4 ), I = 1,RTYPE%NIMPROPERS )

   END SUBROUTINE WRITE_RESIDUE_TYPE

   !----------------------------------
   SUBROUTINE WRITE_TYPE_DATA ( UNIT )
   !----------------------------------
   
   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: UNIT

   ! . Write out the type counter.
   WRITE ( UNIT ) NTYPES

   ! . Write out the type array.
   IF ( NTYPES > 0 ) WRITE ( UNIT ) TYPES

   END SUBROUTINE WRITE_TYPE_DATA

   !--------------------------------------------
   SUBROUTINE WRITE_VARIANT_TYPE ( UNIT, VTYPE )
   !--------------------------------------------
   
   ! . Scalar arguments.
   INTEGER,            INTENT(IN) :: UNIT
   TYPE(VARIANT_TYPE), INTENT(IN) :: VTYPE

   ! . Local scalars.
   INTEGER :: I, J

   ! . Write out the variant scalars.
   WRITE ( UNIT ) VTYPE%NAME,     VTYPE%RESIDUE_NAME, VTYPE%NADDS,      VTYPE%NBONDS, &
                  VTYPE%NCHARGES, VTYPE%NDELETES,     VTYPE%NIMPROPERS

   ! . Write out the variant arrays.
   IF ( VTYPE%NADDS      > 0 ) WRITE ( UNIT ) ( VTYPE%ADDATM(I), VTYPE%ADDCHG(I), VTYPE%ADDTYP(I), I = 1,VTYPE%NADDS )
   IF ( VTYPE%NBONDS     > 0 ) WRITE ( UNIT ) ( ( VTYPE%BONDS(J,I), J = 1,2 ), I = 1,VTYPE%NBONDS )
   IF ( VTYPE%NCHARGES   > 0 ) WRITE ( UNIT ) ( VTYPE%CHGATM(I), VTYPE%CHGCHG(I), I = 1,VTYPE%NCHARGES )
   IF ( VTYPE%NDELETES   > 0 ) WRITE ( UNIT ) ( VTYPE%DELATM(I), I = 1,VTYPE%NDELETES )
   IF ( VTYPE%NIMPROPERS > 0 ) WRITE ( UNIT ) ( ( VTYPE%IMPROPERS(J,I), J = 1,4 ), I = 1,VTYPE%NIMPROPERS )

   END SUBROUTINE WRITE_VARIANT_TYPE

END MODULE MM_FILE_IO
