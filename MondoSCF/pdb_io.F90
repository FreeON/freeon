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
!                       The PDB Coordinate File I/O Module
!===============================================================================
!
! . Subroutines:
!
!   PDB_DEFINE                     Define a system from a PDB file.
!   PDB_READ                       Read the coordinates from a PDB file.
!   PDB_WRITE                      Write a PDB file.
!
! . Notes:
!
!   The  PDB format that is recognized is relatively simple as many of the
!   special cards are ignored. Only "ATOM" and "HETATM" cards are read. Also
!   PDB files only recognize residue names of 3 characters and atoms names of
!   4 characters. In addition, subsystem names are not used although in the
!   format for this module they are added on at the end of the line (columns
!   73 to 76). All files should terminate with an "END".
!
!   The use of the subroutines is basically similar to that of the subroutines
!   in the module COORDINATE_IO. However, the optional arguments DATA are not
!   found (ATMCRD is always used as the repository or the source of the data)
!   and no symmetry information is processed. The other problem is that as
!   the element type for an atom is not indicated in a PDB file, the ATMMAS
!   and ATMNUM arrays are NOT filled by either PDB_DEFINE or PDB_READ.
!
!===============================================================================
MODULE PDB_IO

! . Module declarations.
USE CONSTANTS,   ONLY : UNDEFINED
USE DEFINITIONS, ONLY : DP
USE FILES,       ONLY : NEXT_UNIT
USE IO_UNITS,    ONLY : INPUT, OUTPUT
USE STATUS,      ONLY : ERROR

USE ATOMS
USE SEQUENCE
USE SYMMETRY

IMPLICIT NONE
PUBLIC

! . Module parameters.
CHARACTER ( LEN = 46 ), PARAMETER :: CHARMM_FORMAT = "(A6,I5,1X,A4,1X,A4,2X,I4,3X,3F8.3,2F6.2,6X,A4)", &
                                     PDB_FORMAT    = "(A6,I5,1X,A4,A,A3,2X,I4,4X,3F8.3,2F6.2,6X,A4)"
!                                     PDB_FORMAT    = "(A6,I5,1X,A4,1X,A3,2X,I4,4X,3F8.3,2F6.2,6X,A4)"


! . Module scalars.
CHARACTER ( LEN = 46 ) :: FORMAT
CHARACTER ( LEN = 76 ) :: LINE

!==============================================================================
CONTAINS
!==============================================================================

   !-------------------------------------
   SUBROUTINE PDB_DEFINE ( FILE, CHARMM )
   !-------------------------------------

   ! . Optional scalar arugments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: FILE
   LOGICAL,               INTENT(IN), OPTIONAL :: CHARMM

   ! . Local scalars.
   CHARACTER ( LEN = 4 ) :: ANAME, RNAME, SNAME, SNAME_OLD
   CHARACTER ( LEN = 6 ) :: DUMMY
   CHARACTER ( LEN = 1 ) :: ACR   
   INTEGER               :: I, IATOM, IOSTAT, IRES, IRES_OLD, ISUB, ISUB_OLD, NATMIN, NRESIN, NSUBIN, TUNIT, UNIT
   REAL ( KIND = DP )    :: DUM1, DUM2, X, Y, Z

   !---------------------------------------------------------------------------
   ! . Initialization.
   !---------------------------------------------------------------------------
   ! . Check to see if the FILE argument is present.
   IF ( PRESENT ( FILE ) ) THEN

      ! . Get the next unit number.
      UNIT = NEXT_UNIT ( )

      ! . Open the file.
      OPEN ( UNIT, FILE = FILE, ACTION = "READ", STATUS = "OLD", IOSTAT = IOSTAT )

      ! . If there has been an error.
      IF ( IOSTAT /= 0 ) CALL ERROR ( "PDB_DEFINE", "I/O Error.", IOSTAT )

   ! . The unit for reading is the input stream.
   ELSE
      UNIT = INPUT
   END IF

   ! . Initialize the atom and sequence data structures.
   CALL ATOMS_INITIALIZE
   CALL SEQUENCE_INITIALIZE

   ! . Get the next unit number.
   TUNIT = NEXT_UNIT ( )

   ! . Open the file.
   OPEN ( TUNIT, ACTION = "READWRITE", FORM = "UNFORMATTED", STATUS = "SCRATCH", IOSTAT = IOSTAT )

   ! . If there has been an error.
   IF ( IOSTAT /= 0 ) CALL ERROR ( "PDB_DEFINE", "I/O Error.", IOSTAT )

   ! . Set the format.
   FORMAT = PDB_FORMAT
   IF ( PRESENT ( CHARMM ) ) THEN
      IF ( CHARMM ) FORMAT = CHARMM_FORMAT
   END IF

   !---------------------------------------------------------------------------
   ! . Read in the file and write out the necessary data to a scratch file.
   !---------------------------------------------------------------------------
   ! . Initialization.
   IRES_OLD  = 0
   NATMIN    = 0
   NRESIN    = 0
   NSUBIN    = 0
   SNAME_OLD = "????"

   ! . Top of the loop for reading lines.
   10 CONTINUE

   ! . Read in a line from the file.
   READ ( UNIT, "(A)", END = 20, IOSTAT = IOSTAT ) LINE

   ! . The line is an atom line.
   IF ( ( LINE(1:4) == "ATOM" ) .OR. ( LINE(1:6) == "HETATM" ) ) THEN

      ! . Read in the data on the line.
      READ ( LINE, FORMAT, IOSTAT = IOSTAT ) DUMMY, IATOM, ANAME, ACR, RNAME, IRES, X, Y, Z, DUM1, DUM2, SNAME
!     READ ( LINE, FORMAT, IOSTAT = IOSTAT ) DUMMY, IATOM, ANAME, RNAME, IRES, X, Y, Z, DUM1, DUM2, SNAME

    IF(ACR == ' ' .OR. ACR == 'A') THEN
      ! . Check for an error.
      IF ( IOSTAT /= 0 ) THEN
         WRITE ( OUTPUT, "(/A)" ) "Line = ", LINE(1:LEN_TRIM(LINE))
         CALL ERROR ( "PDB_DEFINE", "Error decoding line." )
      END IF

      ! . Reset the names.
      ANAME = ADJUSTL ( ANAME )
      RNAME = ADJUSTL ( RNAME )
      SNAME = ADJUSTL ( SNAME )

      ! . Increment the number of atoms.
      NATMIN = NATMIN + 1

      ! . Check the subsystem name.
      IF ( SNAME /= SNAME_OLD ) THEN
         NSUBIN    = NSUBIN + 1
         SNAME_OLD = SNAME
      END IF

      ! . Check the residue number.
      IF ( ( SNAME /= SNAME_OLD ) .OR. ( IRES /= IRES_OLD ) ) THEN
         IRES_OLD  = IRES
         NRESIN    = NRESIN + 1
      END IF

      ! . Write out the data to the scratch file.
      WRITE ( TUNIT ) NATMIN, IRES, NSUBIN, ANAME, RNAME, SNAME, X, Y, Z

      ! . Read in the next line.
    END IF !!! ON ACR

      GO TO 10

   ! . The line is a termination line.
   ELSE IF ( LINE(1:3) == "END" ) THEN
      CONTINUE

   ! . Ignore all other lines.
   ELSE
      GO TO 10
   END IF

   ! . Check for an error.
   IF ( IOSTAT /= 0 ) CALL ERROR ( "PDB_DEFINE", "I/O Error.", IOSTAT )

   ! . End of file.
   20 CONTINUE

   !---------------------------------------------------------------------------
   ! . Reprocess the data on the scratch file.
   !---------------------------------------------------------------------------
   ! . Allocate the atom and sequence data arrays.
   CALL ATOMS_ALLOCATE ( NATMIN )
   CALL SEQUENCE_ALLOCATE ( NRESIN, NSUBIN )

   ! . Rewind the scratch file.
   REWIND ( TUNIT )

   ! . Initialize some counters.
   IRES_OLD = 0
   ISUB_OLD = 0

   ! . Loop over the atom data in the scratch file.
   DO IATOM = 1,MM_NATOMS

      ! .Read in data from the scratch file.
      READ ( TUNIT ) I, IRES, ISUB, ANAME, RNAME, SNAME, X, Y, Z

      ! . Check for a new subsystem.
      IF ( ISUB /= ISUB_OLD ) THEN

         ! . Check the subsystem names.
         DO I = 1,ISUB_OLD
            IF ( SUBNAM(I) == SNAME ) THEN
               CALL ERROR ( "PDB_DEFINE", "Duplicate subsystem names: "//SNAME//"." )
            END IF
         END DO

         ! . Assign the new subsystem.
         SUBIND(ISUB) = IRES - 1
         SUBNAM(ISUB) = SNAME

         ! . Reset the subsystem counter.
         ISUB_OLD = ISUB

      END IF

      ! . Check for a new residue.
      IF ( ( ISUB /= ISUB_OLD ) .OR. ( IRES /= IRES_OLD ) ) THEN

         ! . Assign the new residue.
         RESIND(IRES) = IATOM - 1
         RESNAM(IRES) = RNAME

         ! . Reset the residue counter.
         IRES_OLD = IRES

      END IF

      ! . Check the atom names.
      DO I = (RESIND(IRES)+1),(IATOM-1)
         IF ( ATMNAM(I) == ANAME ) THEN
            CALL ERROR ( "PDB_DEFINE", "Duplicate atom name "//ANAME//" in residue "//RNAME//"." )
         END IF
      END DO

      ! . Assign the atom data.
      ATMNAM(IATOM)   = ANAME
      ATMCRD(1,IATOM) = X
      ATMCRD(2,IATOM) = Y
      ATMCRD(3,IATOM) = Z

   END DO

   ! . Fill the last RESIND and SUBIND elements.
   RESIND(NRESID+1)  = MM_NATOMS
   SUBIND(NSUBSYS+1) = NRESID

   !---------------------------------------------------------------------------
   ! . Finish up.
   !---------------------------------------------------------------------------
   ! . Initialize the atom mass and atom number arrays.
   ATMMAS = 0.0_DP
   ATMNUM = -1

   ! . Close the input file if necessary.
   IF ( UNIT /= INPUT ) CLOSE ( UNIT, IOSTAT = IOSTAT )

   ! . Close the scratch file.
   CLOSE ( TUNIT, IOSTAT = IOSTAT )

   ! . Check for an error.
   IF ( IOSTAT /= 0 ) CALL ERROR ( "PDB_DEFINE", "I/O Error.", IOSTAT )

   ! . Do some printing.
   IF ( PRESENT ( FILE ) ) THEN
      WRITE ( OUTPUT, "(/A)" ) "System read in PDB format from " // FILE
   ELSE
      WRITE ( OUTPUT, "(/A)" ) "System read in PDB format from the input stream."
   END IF

   END SUBROUTINE PDB_DEFINE

   !-----------------------------------
   SUBROUTINE PDB_READ ( FILE, CHARMM )
   !-----------------------------------

   ! . Optional scalar arugments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: FILE
   LOGICAL,               INTENT(IN), OPTIONAL :: CHARMM

   ! . Local scalars.
   CHARACTER ( LEN = 4 ) :: ANAME, RNAME, SNAME
   CHARACTER ( LEN = 6 ) :: DUMMY
   CHARACTER ( LEN = 1 ) :: ACR   
   INTEGER               :: I, IATOM, IOSTAT, IRES, ISUB, LENNAM, UNIT
   REAL ( KIND = DP )    :: DUM1, DUM2, X, Y, Z

   ! . Return if no residues have been defined.
   IF ( NRESID <= 0 ) RETURN

   ! . Check to see if the FILE argument is present.
   IF ( PRESENT ( FILE ) ) THEN

      ! . Get the next unit number.
      UNIT = NEXT_UNIT ( )

      ! . Open the file.
      OPEN ( UNIT, FILE = FILE, ACTION = "READ", STATUS = "OLD", IOSTAT = IOSTAT )

      ! . Check for an error.
      IF ( IOSTAT /= 0 ) CALL ERROR ( "PDB_READ", "I/O Error.", IOSTAT )

   ! . The unit for reading is the input stream.
   ELSE
      UNIT = INPUT
   END IF

   ! . Set the format and the name length.
   FORMAT = PDB_FORMAT
   LENNAM = 3
   IF ( PRESENT ( CHARMM ) ) THEN
      IF ( CHARMM ) THEN
         FORMAT = CHARMM_FORMAT
         LENNAM = 4
      END IF
   END IF

   ! . Reallocate and initialize the ATMCRD array.
   IF ( ALLOCATED ( ATMCRD ) ) DEALLOCATE ( ATMCRD )
   ALLOCATE ( ATMCRD(1:3,1:MM_NATOMS) )
   ATMCRD = UNDEFINED

   ! . Top of the loop for reading lines.
   10 CONTINUE

   ! . Read in a line from the file.
   READ ( UNIT, "(A)", END = 30, IOSTAT = IOSTAT ) LINE

   ! . The line is an atom line.
   IF ( ( LINE(1:4) == "ATOM" ) .OR. ( LINE(1:6) == "HETATM" ) ) THEN

      ! . Read in the data on the line.
!     READ ( LINE, FORMAT, IOSTAT = IOSTAT ) DUMMY, IATOM, ANAME, RNAME, IRES, X, Y, Z, DUM1, DUM2, SNAME
      READ ( LINE, FORMAT, IOSTAT = IOSTAT ) DUMMY, IATOM, ANAME, ACR, RNAME, IRES, X, Y, Z, DUM1, DUM2, SNAME

      ! . Check for an error.
      IF ( IOSTAT /= 0 ) THEN
         WRITE ( OUTPUT, "(/A)" ) "Line = ", LINE(1:LEN_TRIM(LINE))
         CALL ERROR ( "PDB_READ", "Error decoding line." )
      END IF

      ! . Reset the names.
      ANAME = ADJUSTL ( ANAME )
      RNAME = ADJUSTL ( RNAME )
      SNAME = ADJUSTL ( SNAME )

      ! . Locate the subsystem.
      DO I = 1,NSUBSYS
         IF ( SNAME == SUBNAM(I)(1:4) ) THEN
            ISUB = I ; GO TO 20
         END IF
      END DO

      ! . No subsystem was located.
      CALL ERROR ( "PDB_READ", "Subsystem "//SNAME//" not found." )

      ! . The subsystem was found.
      20 CONTINUE

      ! . Increment IRES.
      IRES = IRES + SUBIND(ISUB)

      ! . Check the residue.
      IF ( RESNAM(IRES)(1:LENNAM) /= RNAME ) THEN
         CALL ERROR ( "PDB_READ", "Residue mismatch: "//RESNAM(IRES)(1:LENNAM)//" and "//RNAME//"." )
      END IF

      ! . Loop over the atoms in the residue.
      DO I = (RESIND(IRES)+1),RESIND(IRES+1)

         ! . The atom has been found.
         IF ( ( ATMNAM(I)(1:4) == ANAME ) .AND. ( ALL ( ATMCRD(1:3,I) == UNDEFINED ) ) ) THEN
            ATMCRD(1,I) = X
            ATMCRD(2,I) = Y
            ATMCRD(3,I) = Z
            GO TO 10
         END IF

      END DO

      ! . No atom was found.
      CALL ERROR ( "PDB_READ", "Atom "//ANAME//" not found in residue "//RNAME//"." )

   ! . The line is a termination line.
   ELSE IF ( LINE(1:3) == "END" ) THEN
      CONTINUE

   ! . Ignore all other lines.
   ELSE
      GO TO 10
   END IF

   ! . End of file.
   30 CONTINUE

   ! . Close the input file if necessary.
   IF ( UNIT /= INPUT ) CLOSE ( UNIT, IOSTAT = IOSTAT )

   ! . Check for an error.
   IF ( IOSTAT /= 0 ) CALL ERROR ( "PDB_READ", "I/O Error.", IOSTAT )

   ! . Do some printing.
   IF ( PRESENT ( FILE ) ) THEN
      WRITE ( OUTPUT, "(/A)" ) "Coordinates read in PDB format from " // FILE
   ELSE
      WRITE ( OUTPUT, "(/A)" ) "Coordinates read in PDB format from the input stream."
   END IF

   END SUBROUTINE PDB_READ

   !------------------------------------
   SUBROUTINE PDB_WRITE ( FILE, CHARMM )
   !------------------------------------

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: FILE
   LOGICAL,               INTENT(IN), OPTIONAL :: CHARMM

   ! . Local scalars.
   CHARACTER ( LEN = 4 ) :: ANAME
   CHARACTER ( LEN = 1 ) :: ACR   
   INTEGER               :: IATOM, IOSTAT, IRES, ISUB, LENNAM, NATM, NRES, UNIT

   ! . Return if no residues have been defined.
   IF ( NRESID <= 0 ) RETURN

   ! . Check to see if the FILE argument is present.
   IF ( PRESENT ( FILE ) ) THEN

      ! . Get the next unit number.
      UNIT = NEXT_UNIT ( )

      ! . Open the file.
      OPEN ( UNIT, FILE = FILE, ACTION = "WRITE", STATUS = "REPLACE", IOSTAT = IOSTAT )

      ! . If there has been an error return.
      IF ( IOSTAT /= 0 ) CALL ERROR ( "PDB_WRITE", "I/O Error.", IOSTAT )

      ! . Do some printing.
      WRITE ( OUTPUT, "(/A)" ) "Coordinates written in PDB format to " // FILE

   ! . The unit for writing is the output stream.
   ELSE

      ! . Assign the unit number.
      UNIT = OUTPUT

      ! . Print out a blank line.
      WRITE ( OUTPUT, "(/A/)" ) "Coordinates written in PDB format to the output stream."

   END IF

   ! . Set the format and the name length.
   FORMAT = PDB_FORMAT
   LENNAM = 3
   IF ( PRESENT ( CHARMM ) ) THEN
      IF ( CHARMM ) THEN
         FORMAT = CHARMM_FORMAT
         LENNAM = 4
      END IF
   END IF

   ! . Loop over the subsystems.
   DO ISUB = 1,NSUBSYS

      ! . Get the number of residues.
      NRES = SUBIND(ISUB+1) - SUBIND(ISUB)

      ! . Loop over the residues.
      DO IRES = (SUBIND(ISUB)+1),SUBIND(ISUB+1)

         ! . Get the number of atoms.
         NATM = RESIND(IRES+1) - RESIND(IRES)

         ! . Loop over the atoms.
         DO IATOM = (RESIND(IRES)+1),RESIND(IRES+1)

            ! . Get the atom name.
            IF ( ATMNAM(IATOM)(4:4) == " " ) THEN
               ANAME = " " // ATMNAM(IATOM)(1:3)
            ELSE
               ANAME = ATMNAM(IATOM)(1:4)
            END IF

            ! . Write out the atom data.
!           WRITE ( UNIT, FORMAT ) "ATOM  ", IATOM, ANAME, RESNAM(IRES)(1:LENNAM), IRES - SUBIND(ISUB), &
            WRITE ( UNIT, FORMAT ) "ATOM  ", IATOM, ANAME, ACR, RESNAM(IRES)(1:LENNAM), IRES - SUBIND(ISUB), &
                                                  ATMCRD(1:3,IATOM), 0.0_DP, 0.0_DP, SUBNAM(ISUB)(1:4)

         END DO
      END DO
   END DO

   ! . Terminate the file.
   WRITE ( UNIT, "(A3)" ) "END"

   ! . Close the output file if necessary and check for any errors.
   IF ( UNIT /= OUTPUT ) THEN
      CLOSE ( UNIT, IOSTAT = IOSTAT )
      IF ( IOSTAT /= 0 ) CALL ERROR ( "PDB_WRITE", "I/O Error.", IOSTAT )
   END IF

   END SUBROUTINE PDB_WRITE

END MODULE PDB_IO
