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
!                         The Coordinate File I/O Module
!===============================================================================
!
! . Subroutines:
!
!   COORDINATES_DEFINE             Define a structure from a coordinate file.
!   COORDINATES_READ               Read the coordinates from a coordinate file.
!   COORDINATES_WRITE              Write a coordinate file.
!
! . Notes:
!
!   COORDINATES_DEFINE is used to define a system (its composition, sequence
!   and symmetry) from a coordinate file. All existing atoms, sequence and
!   symmmetry data is lost.
!
!   COORDINATES_READ and COORDINATES_WRITE read and write coordinate files
!   respectively. COORDINATES_READ just reads in new coordinate data, the
!   atoms, sequence and symmetry variables already having been defined.
!   The sequence of residues and subsystems must be the same for the
!   subroutine to work, but the atoms within each residue can be in an
!   arbitrary order (and can even be missing). Whereas COORDINATES_READ
!   will only change ATMCRD (at most), COORDINATES_WRITE does not change
!   any atoms, sequence or symmetry data.
!
!   For all three of the above subroutines, the argument FILE gives the file name
!   to read from or write to. If it is not present the standard input or output
!   streams are used (whichever is appropriate). For COORDINATES_READ and
!   COORDINATES_WRITE there is also a second argument which is an optional
!   array that can contain the data to read into or write out from. This is
!   useful when dealing with other atom data (such as forces or velocities).
!
!===============================================================================
MODULE COORDINATE_IO

! . Module declarations.
USE CONSTANTS,   ONLY : UNDEFINED
USE DEFINITIONS, ONLY : DP
USE ELEMENTS,    ONLY : MASS
USE FILES,       ONLY : NEXT_UNIT
USE IO_UNITS,    ONLY : INPUT, OUTPUT
USE PARSING
USE STATUS,      ONLY : ERROR

USE ATOMS
USE SEQUENCE
USE SYMMETRY

IMPLICIT NONE
PUBLIC

!==============================================================================
CONTAINS
!==============================================================================

   !-------------------------------------
   SUBROUTINE COORDINATES_DEFINE ( FILE )
   !-------------------------------------

   ! . Optional scalar arugments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: FILE

   ! . Local scalars.
   INTEGER :: I, IATOM, IOSTAT, IRES, ISUB, JATOM, JSUB, NATM, NATMIN, NRES, NRESIN, NUMATM, NUMRES, NSUBIN, UNIT

   ! . Check to see if the FILE argument is present.
   IF ( PRESENT ( FILE ) ) THEN

      ! . Get the next unit number.
      UNIT = NEXT_UNIT ( )

      ! . Open the file.
      OPEN ( UNIT, FILE = FILE, ACTION = "READ", STATUS = "OLD", IOSTAT = IOSTAT )

      ! . If there has been an error return.
      IF ( IOSTAT /= 0 ) CALL ERROR ( "COORDINATES_DEFINE", "I/O Error.", IOSTAT )

   ! . The unit for reading is the input stream.
   ELSE
      UNIT = INPUT
   END IF

   ! . Set the parsing unit.
   CALL PUSH_UNIT ( UNIT )

   ! . Initialize the atom, sequence and symmetry data structures.
   CALL    ATOMS_INITIALIZE
   CALL SEQUENCE_INITIALIZE
   CALL SYMMETRY_INITIALIZE

   ! . Read in the first line of the file.
   CALL GET_LINE
   CALL GET_INTEGER ( NATMIN )
   CALL GET_INTEGER ( NRESIN )
   CALL GET_INTEGER ( NSUBIN )

   ! . Check the counters.
   IF ( ( NATMIN <= 0 ) .OR. ( NRESIN <= 0 ) .OR. ( NSUBIN <= 0 ) ) THEN
      CALL PARSE_ERROR ( "COORDINATES_DEFINE", "Invalid atoms, residue or subsystem counter." )
   END IF

   ! . Allocate the atom and sequence data arrays.
   CALL ATOMS_ALLOCATE ( NATMIN )
   CALL SEQUENCE_ALLOCATE ( NRESIN, NSUBIN )

   ! . Read in any symmetry definitions.
   CALL SYMMETRY_RECORD_READ

   ! . Initialization.
   NUMATM = 0
   NUMRES = 0

   ! . Loop over the subsystems.
   DO ISUB = 1,NSUBSYS

      ! . Get the first line of the subsystem.
      IF ( ISUB > 1 ) THEN
         CALL GET_LINE
         CALL GET_WORD
      END IF
      IF ( WORD(1:WRDLEN) /= "SUBSYSTEM" ) THEN
         CALL PARSE_ERROR ( "COORDINATES_DEFINE", "SUBSYSTEM label invalid." )
      END IF
      CALL GET_INTEGER ( I )
      CALL GET_WORD
      IF ( ( WRDLEN <= 0 ) .OR. ( WRDLEN > SUBSYSTEM_NAME_LENGTH ) ) THEN
         CALL PARSE_ERROR ( "COORDINATES_DEFINE", "SUBSYSTEM name length invalid." )
      END IF

      ! . Check the subsystem names.
      DO JSUB = 1,(ISUB-1)
         ! . Compare the names.
         IF ( SUBNAM(JSUB) == WORD(1:WRDLEN) ) THEN
            CALL PARSE_ERROR ( "COORDINATES_DEFINE", "Two SUBSYSTEM names are the same." )
         END IF
      END DO

      ! . Save the subsystem name.
      SUBNAM(ISUB) = WORD(1:WRDLEN)

      ! . Fill the subsystem arrays.
      SUBIND(ISUB) = NUMRES

      ! . Get the number of residues in the subsystem.
      CALL GET_LINE
      CALL GET_INTEGER ( NRES )
      IF ( ( NRES <= 0 ) .OR. ( ( NUMRES+NRES ) > NRESID ) ) THEN
         CALL PARSE_ERROR ( "COORDINATES_DEFINE", "Invalid number of residues in a subsystem." )
      END IF

      ! . Loop over the residues in the subsystem.
      DO IRES = 1,NRES

         ! . Increment the NUMRES counter.
         NUMRES = NUMRES + 1

         ! . Get the first line of the residue.
         CALL GET_LINE
         CALL GET_WORD
         IF ( WORD(1:WRDLEN) /= "RESIDUE" ) THEN
            CALL PARSE_ERROR ( "COORDINATES_DEFINE", "RESIDUE label invalid." )
         END IF
         CALL GET_INTEGER ( I )
         IF ( I /= IRES ) THEN
            CALL PARSE_ERROR ( "COORDINATES_DEFINE", "RESIDUE number invalid." )
         END IF
         CALL GET_WORD
         IF ( ( WRDLEN <= 0 ) .OR. ( WRDLEN > RESIDUE_NAME_LENGTH ) ) THEN
            CALL PARSE_ERROR ( "COORDINATES_DEFINE", "RESIDUE name length invalid." )
         END IF

         ! . Save the residue name.
         RESNAM(NUMRES) = WORD(1:WRDLEN)

         ! . Fill the residue arrays.
         RESIND(NUMRES) = NUMATM

         ! . Get the number of atoms in the residue.
         CALL GET_LINE
         CALL GET_INTEGER ( NATM )
         IF ( ( NATM <= 0 ).OR.( (NUMATM+NATM) > MM_NATOMS ) ) THEN
            CALL PARSE_ERROR ( "COORDINATES_DEFINE", "Invalid number of atoms in a residue." )
         END IF

         ! . Loop over the atoms in the residue.
         DO IATOM = 1,NATM

            ! . Increment the atom counter.
            NUMATM = NUMATM + 1

            ! . Parse an atom line ( the order, name and number ).
            CALL GET_LINE
            CALL GET_INTEGER ( I )
            CALL GET_WORD
            IF ( WRDLEN <= ATOM_NAME_LENGTH ) THEN
               ATMNAM(NUMATM) = WORD(1:WRDLEN)
            ELSE
               CALL PARSE_ERROR ( "COORDINATES_DEFINE", "ATOM name length invalid." )
            ENDIF

            ! . Check the atom names.
            DO JATOM = (RESIND(NUMRES)+1),(NUMATM-1)
               ! . Compare the names.
               IF ( ATMNAM(JATOM) == WORD(1:WRDLEN) ) THEN
                  CALL PARSE_ERROR ( "COORDINATES_DEFINE", "Two ATOM names are the same." )
               END IF
            END DO

            ! . Read in the remaining data.
            CALL GET_INTEGER ( ATMNUM(NUMATM) )
            CALL GET_REAL    ( ATMCRD(1,NUMATM) )
            CALL GET_REAL    ( ATMCRD(2,NUMATM) )
            CALL GET_REAL    ( ATMCRD(3,NUMATM) )

         END DO
      END DO
   END DO

   ! . Check NUMATM and NUMRES.
   IF ( ( NUMATM /= MM_NATOMS ) .OR. ( NUMRES /= NRESID ) ) THEN
      CALL PARSE_ERROR ( "COORDINATES_DEFINE", "ATOM or RESIDUE number mismatch." )
   END IF

   ! . Fill the last RESIND and SUBIND elements.
   RESIND(NRESID+1)  = NUMATM
   SUBIND(NSUBSYS+1) = NUMRES

   ! . Fill the atom mass array.
   DO IATOM = 1,MM_NATOMS
      IF ( ATMNUM(IATOM) > 0 ) THEN
         ATMMAS(IATOM) = MASS(ATMNUM(IATOM))
      ELSE
         ATMMAS(IATOM) = 0.0_DP
      END IF
   END DO

   ! . Close the input file if necessary.
   IF ( UNIT /= INPUT ) CLOSE ( UNIT )

   ! . Reset the parsing unit.
   CALL POP_UNIT

   ! . Do some printing.
   IF ( PRESENT ( FILE ) ) THEN
      WRITE ( OUTPUT, "(/A)" ) "System read from " // FILE
   ELSE
      WRITE ( OUTPUT, "(/A)" ) "System read from the input stream."
   END IF

   END SUBROUTINE COORDINATES_DEFINE

   !----------------------------------------------------
   SUBROUTINE COORDINATES_READ ( FILE, DATA, SELECTION )
   !----------------------------------------------------

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: FILE

   ! . Optional array arguments.
   LOGICAL,            DIMENSION(1:MM_NATOMS),     INTENT(IN),    OPTIONAL :: SELECTION
   REAL ( KIND = DP ), DIMENSION(1:3,1:MM_NATOMS), INTENT(INOUT), OPTIONAL :: DATA

   ! . Local scalars.
   CHARACTER ( LEN = 16 ) :: TAG
   INTEGER                :: I, IATOM, IOSTAT, IRES, ISUB, JATOM, LTAG, NATM, NATMIN, NRES, NRESIN, NSUB, NSUBIN, UNIT

   ! . Local arrays.
   LOGICAL,            DIMENSION(1:MM_NATOMS)     :: QATOM
   LOGICAL,            DIMENSION(1:NRESID)     :: QRESIDUE
   LOGICAL,            DIMENSION(1:NSUBSYS)    :: QSUBSYSTEM
   REAL ( KIND = DP ), DIMENSION(1:3,1:MM_NATOMS) :: ATMDAT

   ! . No residues have been defined.
   IF ( NRESID <= 0 ) RETURN

   ! . Check to see if the FILE argument is present.
   IF ( PRESENT ( FILE ) ) THEN

      ! . Get the next unit number.
      UNIT = NEXT_UNIT ( )

      ! . Open the file.
      OPEN ( UNIT, FILE = FILE, ACTION = "READ", STATUS = "OLD", IOSTAT = IOSTAT )

      ! . Check for an error.
      IF ( IOSTAT /= 0 ) CALL ERROR ( "COORDINATES_READ", "I/O Error.", IOSTAT )

   ! . The unit for reading is the input stream.
   ELSE
      UNIT = INPUT
   END IF

   ! . Check to see if the SELECTION argument is present.
   IF ( PRESENT ( SELECTION ) ) THEN

      ! . Save the SELECTION array in QATOM.
      QATOM = SELECTION

      ! . Determine which residues are active.
      DO IRES = 1,NRESID
         QRESIDUE(IRES) = ANY ( QATOM(RESIND(IRES)+1:RESIND(IRES+1)) )
      END DO

      ! . Determine which subsystems are active.
      DO ISUB = 1,NSUBSYS
         QSUBSYSTEM(ISUB) = ANY ( QRESIDUE(SUBIND(ISUB)+1:SUBIND(ISUB+1)) )
      END DO

   ! . All data is to be read.
   ELSE
      QATOM = .TRUE. ; QRESIDUE = .TRUE. ; QSUBSYSTEM = .TRUE.
   END IF

   ! . Determine the total number of atoms, residues and subsystems.
   NATM = COUNT ( QATOM ) ; NRES = COUNT ( QRESIDUE ) ; NSUB = COUNT ( QSUBSYSTEM )

   ! . Initialize the ATMDAT array.
   ATMDAT = UNDEFINED

   ! . Set the parsing unit.
   CALL PUSH_UNIT ( UNIT )

   ! . Read in the first line of the file.
   CALL GET_LINE
   CALL GET_INTEGER ( NATMIN )
   CALL GET_INTEGER ( NRESIN )
   CALL GET_INTEGER ( NSUBIN )

   ! . Check the counters.
   IF ( ( NATMIN > MM_NATOMS ) .OR. ( NRESIN /= NRES ) .OR. ( NSUBIN /= NSUB ) ) THEN
      CALL PARSE_ERROR ( "COORDINATES_READ", "Invalid atoms, residue or subsystem counter." )
   END IF

   ! . Read in any symmetry definitions.
   CALL SYMMETRY_RECORD_READ

   ! . Initialize some counters.
   NATM = 0 ; NSUB = 0

   ! . Loop over the subsystems.
   DO ISUB = 1,NSUBSYS

      ! . Check to see if the subsystem is to be included.
      IF ( .NOT. QSUBSYSTEM(ISUB) ) CYCLE

      ! . Increment the subsystem number.
      NSUB = NSUB + 1

      ! . Get the first line of the subsystem.
      IF ( NSUB > 1 ) THEN
         CALL GET_LINE
         CALL GET_WORD
      END IF
      IF ( WORD(1:WRDLEN) /= "SUBSYSTEM" ) THEN
         CALL PARSE_ERROR ( "COORDINATES_READ", "SUBSYSTEM label invalid." )
      END IF
      CALL GET_INTEGER ( I )
      CALL GET_WORD
      IF ( ( WRDLEN <= 0 ) .OR. ( WRDLEN > SUBSYSTEM_NAME_LENGTH ) ) THEN
         CALL PARSE_ERROR ( "COORDINATES_READ", "SUBSYSTEM name length invalid." )
      END IF

      ! . Check the subsystem name.
      IF ( SUBNAM(ISUB) /= WORD(1:WRDLEN) ) THEN
         CALL PARSE_ERROR ( "COORDINATES_READ", "SUBSYSTEM name mismatch." )
      END IF

      ! . Get the number of residues in the subsystem.
      CALL GET_LINE
      CALL GET_INTEGER ( NRES )
      IF ( NRES /= COUNT ( QRESIDUE(SUBIND(ISUB)+1:SUBIND(ISUB+1)) ) ) THEN
         CALL PARSE_ERROR ( "COORDINATES_READ", "Invalid number of residues in a subsystem." )
      END IF

      ! . Loop over the residues in the subsystem.
      DO IRES = (SUBIND(ISUB)+1),SUBIND(ISUB+1)

         ! . Check to see if the residue is to be included.
         IF ( .NOT. QRESIDUE(IRES) ) CYCLE

         ! . Get the first line of the residue.
         CALL GET_LINE
         CALL GET_WORD
         IF ( WORD(1:WRDLEN) /= "RESIDUE" ) THEN
            CALL PARSE_ERROR ( "COORDINATES_READ", "RESIDUE label invalid." )
         END IF
         CALL GET_INTEGER ( I )
         IF ( I /= ( IRES - SUBIND(ISUB) ) ) THEN
            CALL PARSE_ERROR ( "COORDINATES_READ", "RESIDUE number invalid." )
         END IF
         CALL GET_WORD
         IF ( ( WRDLEN <= 0 ) .OR. ( WRDLEN > RESIDUE_NAME_LENGTH ) ) THEN
            CALL PARSE_ERROR ( "COORDINATES_READ", "RESIDUE name length invalid." )
         END IF

         ! . Check the residue name.
         IF ( RESNAM(IRES) /= WORD(1:WRDLEN) ) THEN
            CALL PARSE_ERROR ( "COORDINATES_READ", "RESIDUE name mismatch." )
         END IF

         ! . Get the number of atoms in the residue.
         CALL GET_LINE
         CALL GET_INTEGER ( NATM )
         IF ( NATM > ( RESIND(IRES+1) - RESIND(IRES) ) ) THEN
            CALL PARSE_ERROR ( "COORDINATES_READ", "Too many atoms in a residue." )
         END IF

         ! . Read in the atoms in the residue.
         DO IATOM = 1,NATM

            ! . Parse an atom line ( Atom's order, name and number ).
            CALL GET_LINE
            CALL GET_INTEGER ( I )
            CALL GET_WORD
            IF ( WRDLEN > ATOM_NAME_LENGTH ) THEN
               CALL PARSE_ERROR ( "COORDINATES_READ", "ATOM name length invalid." )
            ENDIF

            ! . Loop over the atoms in the definition.
            DO JATOM = (RESIND(IRES)+1),RESIND(IRES+1)

               ! . Check the name.
               IF ( ATMNAM(JATOM) == WORD(1:WRDLEN) ) THEN

                  ! . Process the line only for selected atoms.
                  IF ( QATOM(JATOM) ) THEN

                     ! . Parse the remaining data.
                     CALL GET_INTEGER ( I )
                     IF ( ATMNUM(JATOM) /= I ) THEN
                        CALL PARSE_ERROR ( "COORDINATES_READ", "ATOM number mismatch." )
                     END IF

                     ! . Get the coordinates.
                     CALL GET_REAL ( ATMDAT(1,JATOM) )
                     CALL GET_REAL ( ATMDAT(2,JATOM) )
                     CALL GET_REAL ( ATMDAT(3,JATOM) )

                  END IF

                  ! . Skip out of the loop.
                  GO TO 10

               END IF
            END DO

            ! . No match was found.
            CALL PARSE_ERROR ( "COORDINATES_READ", "Unknown ATOM name in residue." )

            ! . End of the loop.
            10 CONTINUE

         END DO
      END DO
   END DO

   ! . Close the input file if necessary.
   IF ( UNIT /= INPUT ) CLOSE ( UNIT )

   ! . Reset the parsing unit.
   CALL POP_UNIT

   ! . Put data in DATA if the DATA argument is present.
   IF ( PRESENT ( DATA ) ) THEN

      ! . Transfer data only for the selected atoms.
      DO IATOM = 1,MM_NATOMS
         IF ( QATOM(IATOM) ) DATA(1:3,IATOM) = ATMDAT(1:3,IATOM)
      END DO
      TAG = "Data"

   ! . By default put the data in ATMCRD.
   ELSE

      ! . Save the existing ATMCRD if necessary.
      IF ( ( SIZE ( ATMCRD, 1 ) == 3 ) .AND. ( SIZE ( ATMCRD, 2 ) == MM_NATOMS ) ) THEN
         DO IATOM = 1,MM_NATOMS
            IF ( QATOM(IATOM) ) ATMCRD(1:3,IATOM) = ATMDAT(1:3,IATOM)
         END DO
      ! . Reallocate ATMCRD.
      ELSE
         IF ( ALLOCATED ( ATMCRD ) ) DEALLOCATE ( ATMCRD )
         ALLOCATE ( ATMCRD(1:3,1:MM_NATOMS) )
         ATMCRD = ATMDAT
      END IF
      TAG = "Coordinates"

   END IF

   ! . Get the length of the tag.
   LTAG = LEN_TRIM ( TAG )

   ! . Do some printing.
   IF ( PRESENT ( FILE ) ) THEN
      WRITE ( OUTPUT, "(/A)" ) TAG(1:LTAG)//" read from " // FILE
   ELSE
      WRITE ( OUTPUT, "(/A)" ) TAG(1:LTAG)//" read from the input stream."
   END IF

   END SUBROUTINE COORDINATES_READ

   !-----------------------------------------------------
   SUBROUTINE COORDINATES_WRITE ( FILE, DATA, SELECTION )
   !-----------------------------------------------------

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: FILE

   ! . Optional array arguments.
   LOGICAL,            DIMENSION(1:MM_NATOMS),     INTENT(IN), OPTIONAL :: SELECTION
   REAL ( KIND = DP ), DIMENSION(1:3,1:MM_NATOMS), INTENT(IN), OPTIONAL :: DATA

   ! . Local scalars.
   CHARACTER ( LEN = 16 ) :: TAG
   INTEGER                :: IATOM, IOSTAT, IRES, ISUB, LTAG, NATM, NLAST, NRES, NSUB, UNIT

   ! . Local arrays.
   LOGICAL,            DIMENSION(1:MM_NATOMS)     :: QATOM
   LOGICAL,            DIMENSION(1:NRESID)     :: QRESIDUE
   LOGICAL,            DIMENSION(1:NSUBSYS)    :: QSUBSYSTEM
   REAL ( KIND = DP ), DIMENSION(1:3,1:MM_NATOMS) :: ATMDAT

   ! . Return if no residues have been defined.
   IF ( NRESID <= 0 ) RETURN

   ! . Write out DATA if the DATA argument is present.
   IF ( PRESENT ( DATA ) ) THEN
      ATMDAT = DATA
      TAG    = "Data"
   ! . By default write out the contents of ATMCRD.
   ELSE
      ATMDAT = ATMCRD
      TAG    = "Coordinates"
   END IF

   ! . Get the length of the tag.
   LTAG = LEN_TRIM ( TAG )

   ! . Check to see if the FILE argument is present.
   IF ( PRESENT ( FILE ) ) THEN

      ! . Get the next unit number.
      UNIT = NEXT_UNIT ( )

      ! . Open the file.
      OPEN ( UNIT, FILE = FILE, ACTION = "WRITE", STATUS = "REPLACE", IOSTAT = IOSTAT )

      ! . If there has been an error return.
      IF ( IOSTAT /= 0 ) CALL ERROR ( "COORDINATES_WRITE", "I/O Error.", IOSTAT )

      ! . Do some printing.
      WRITE ( OUTPUT, "(/A)" ) TAG(1:LTAG)//" written to " // FILE

   ! . The unit for writing is the output stream.
   ELSE

      ! . Assign the unit number.
      UNIT = OUTPUT

      ! . Do some printing.
      WRITE ( OUTPUT, "(/A/)" ) TAG(1:LTAG)//" written to the output stream."

   END IF

   ! . Check to see if the selection array is present.
   IF ( PRESENT ( SELECTION ) ) THEN

      ! . Save the SELECTION array in QATOM.
      QATOM = SELECTION

      ! . Determine which residues are active.
      DO IRES = 1,NRESID
         QRESIDUE(IRES) = ANY ( QATOM(RESIND(IRES)+1:RESIND(IRES+1)) )
      END DO

      ! . Determine which subsystems are active.
      DO ISUB = 1,NSUBSYS
         QSUBSYSTEM(ISUB) = ANY ( QRESIDUE(SUBIND(ISUB)+1:SUBIND(ISUB+1)) )
      END DO

   ! . All data is to be printed.
   ELSE
      QATOM = .TRUE. ; QRESIDUE = .TRUE. ; QSUBSYSTEM = .TRUE.
   END IF

   ! . Determine the total number of atoms, residues and subsystems.
   NATM = COUNT ( QATOM ) ; NRES = COUNT ( QRESIDUE ) ; NSUB = COUNT ( QSUBSYSTEM )

   ! . Write out the separator.
   WRITE ( UNIT, "('!',79('='))" )

   ! . Write out the number of atoms, residues and subsystems.
   WRITE ( UNIT, "(3I6,A)" ) NATM, NRES, NSUB, " ! # of atoms, residues and subsystems."

   ! . Write out any symmetry information.
   CALL SYMMETRY_RECORD_WRITE ( UNIT )

   ! . Initialize some counters.
   NATM = 0 ; NSUB = 0

   ! . Loop over the subsystems.
   DO ISUB = 1,NSUBSYS

      ! . Check to see if the subsystem is to be included.
      IF ( .NOT. QSUBSYSTEM(ISUB) ) CYCLE

      ! . Increment the subsystem number.
      NSUB = NSUB + 1

      ! . Find the number of the last active residue in the subsystem.
      DO IRES = SUBIND(ISUB+1),(SUBIND(ISUB)+1),-1
         IF ( QRESIDUE(IRES) ) THEN
            NLAST = IRES ; EXIT
         END IF
      END DO

      ! . Write out a subsystem header.
      WRITE ( UNIT, "('!',79('='))" )
      WRITE ( UNIT, "(A,I6,2X,A)" ) "Subsystem", NSUB, SUBNAM(ISUB)
      WRITE ( UNIT, "(I6,A)" ) COUNT ( QRESIDUE(SUBIND(ISUB)+1:SUBIND(ISUB+1)) ), " ! # of residues."
      WRITE ( UNIT, "('!',79('='))" )

      ! . Loop over the residues.
      DO IRES = (SUBIND(ISUB)+1),SUBIND(ISUB+1)

         ! . Check to see if the residue is to be included.
         IF ( .NOT. QRESIDUE(IRES) ) CYCLE

         ! . Write out a residue header.
         WRITE ( UNIT, "(A,I6,2X,A)" ) "Residue", ( IRES - SUBIND(ISUB) ), RESNAM(IRES)
         WRITE ( UNIT, "(I6,A)" ) COUNT ( QATOM(RESIND(IRES)+1:RESIND(IRES+1)) ), " ! # of atoms."

         ! . Loop over the atoms.
         DO IATOM = (RESIND(IRES)+1),RESIND(IRES+1)

            ! . Check to see if the atom is to be included.
            IF ( .NOT. QATOM(IATOM) ) CYCLE

            ! . Increment the atom number.
            NATM = NATM + 1

            ! . Write out the atom data.
            WRITE ( UNIT, "(I6,3X,A8,3X,I4,2X,3F18.10)" ) NATM, ATMNAM(IATOM), ATMNUM(IATOM), ATMDAT(1:3,IATOM)

         END DO

         ! . Write out a comment separator.
         IF ( IRES /= NLAST ) WRITE ( UNIT, "('!',79('-'))" )

      END DO
   END DO

   ! . Terminate the file.
   WRITE ( UNIT, "('!',79('='))" )

   ! . Close the output file if necessary.
   IF ( UNIT /= OUTPUT ) CLOSE ( UNIT )

   END SUBROUTINE COORDINATES_WRITE

END MODULE COORDINATE_IO
