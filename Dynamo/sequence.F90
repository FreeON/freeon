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
!                              The Sequence Module
!===============================================================================
!
! . Parameter data:
!
!   RESIDUE_NAME_LENGTH            The maximum length of a residue name.
!   SUBSYSTEM_NAME_LENGTH          The maximum length of a subsystem name.
!
! . Scalar data:
!
!   NRESID                         The number of residues.
!   NSUBSYS                        The number of subsystems.
!
! . Array data:
!
!   RESIND                         The residue index array (to atoms).
!   RESNAM                         The residue names.
!
!   SUBIND                         The subsystem index array (to residues).
!   SUBNAM                         The subsystem names.
!
! . Subroutines:
!
!   SEQUENCE_ALLOCATE              Create the sequence data structure.
!   SEQUENCE_INITIALIZE            Initialize the sequence data structure.
!   SEQUENCE_PRINT                 Print the sequence.
!   SEQUENCE_SUMMARY               Print a summary of the sequence.
!   SEQUENCE_WRITE                 Write out a sequence file.
!
! . Notes:
!
!   The sequence module stores data about the residues and subsystems into
!   which a system is divided. In addition to the number and names of the
!   residues and subsystems there are index arrays.
!
!   The index arrays are structured as follows. For residue i, the element
!   RESIND(i)+1 gives the first atom of the residue and RESIND(i+1) gives
!   the last element. The total number of atoms in the residue is, thus,
!   RESIND(i+1) - RESIND(i). Note that the total length of the array is
!   NRESID + 1, not NRESID. The SUBIND array has the same structure, except
!   that it gives an index into the residue array for each subsystem.
!
!===============================================================================
MODULE SEQUENCE

! . Module declarations.
USE DEFINITIONS, ONLY : LINE_LENGTH
USE FILES,       ONLY : NEXT_UNIT
USE IO_UNITS,    ONLY : OUTPUT
USE STATUS,      ONLY : ERROR

IMPLICIT NONE
PUBLIC
SAVE

! . Parameter declarations.
INTEGER, PARAMETER ::   RESIDUE_NAME_LENGTH = 32, &
                      SUBSYSTEM_NAME_LENGTH = 32

! . Scalar data.
INTEGER :: NRESID = 0, NSUBSYS = 0

! . Character array data.
CHARACTER ( LEN =   RESIDUE_NAME_LENGTH ), ALLOCATABLE, DIMENSION(:) :: RESNAM
CHARACTER ( LEN = SUBSYSTEM_NAME_LENGTH ), ALLOCATABLE, DIMENSION(:) :: SUBNAM

! . Integer array data.
INTEGER, ALLOCATABLE, DIMENSION(:) :: RESIND, SUBIND

!==============================================================================
CONTAINS
!==============================================================================

   !------------------------------------------
   SUBROUTINE SEQUENCE_ALLOCATE ( NRES, NSUB )
   !------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: NRES, NSUB

   ! . Initialize the scalar variables.
   NRESID  = NRES   
   NSUBSYS = NSUB

   ! . Allocate all the arrays.
   ALLOCATE ( RESIND(1:NRESID+1), RESNAM(1:NRESID), SUBIND(1:NSUBSYS+1), SUBNAM(1:NSUBSYS) )

   END SUBROUTINE SEQUENCE_ALLOCATE

   !-----------------------------
   SUBROUTINE SEQUENCE_INITIALIZE
   !-----------------------------

   ! . Initialize the scalar variables.
   NRESID  = 0   
   NSUBSYS = 0

   ! . Deallocate all the arrays.
   IF ( ALLOCATED ( RESIND ) ) DEALLOCATE ( RESIND )
   IF ( ALLOCATED ( RESNAM ) ) DEALLOCATE ( RESNAM )
   IF ( ALLOCATED ( SUBIND ) ) DEALLOCATE ( SUBIND )
   IF ( ALLOCATED ( SUBNAM ) ) DEALLOCATE ( SUBNAM )

   END SUBROUTINE SEQUENCE_INITIALIZE

   !------------------------
   SUBROUTINE SEQUENCE_PRINT
   !------------------------

   ! . Local scalars.
   INTEGER :: I, ISUB, LENGTH

   ! . Check that a sequence is present.
   IF ( ( NRESID > 0 ) .AND. ( NSUBSYS > 0 ) ) THEN

      ! . Write out the header.
      WRITE ( OUTPUT, "(/32('-'),A,31('-'))" ) " System Sequence "

      ! . Loop over the subsystems.
      DO ISUB = 1,NSUBSYS

         ! . Write out the subsystem name.
         LENGTH = LEN_TRIM ( SUBNAM(ISUB) )
         WRITE ( OUTPUT, "(A,I4,2X,A)" ) "Subsystem ", ISUB, SUBNAM(ISUB)(1:LENGTH)//":"

         ! . Write out the residue names.
         WRITE ( OUTPUT, "(2(8X,A))" ) ( RESNAM(I), I = (SUBIND(ISUB)+1),SUBIND(ISUB+1) )
      END DO

      ! . Write out the terminator.
      WRITE ( OUTPUT, "(80('-'))" )

   END IF

   END SUBROUTINE SEQUENCE_PRINT

   !--------------------------
   SUBROUTINE SEQUENCE_SUMMARY
   !--------------------------

   ! . Local scalars.
   INTEGER :: I, ISUB, J, LENGTH, N, NRES

   ! . Local arrays.
   CHARACTER ( LEN = RESIDUE_NAME_LENGTH ), ALLOCATABLE, DIMENSION(:) :: NAME
   INTEGER,                                 ALLOCATABLE, DIMENSION(:) :: FREQUENCY

   ! . Check that a sequence is present.
   IF ( ( NRESID > 0 ) .AND. ( NSUBSYS > 0 ) ) THEN

      ! . Write out the header.
      WRITE ( OUTPUT, "(/24('-'),A,24('-'))"  ) " Summary of the System Sequence "

      ! . Loop over the subsystems.
      DO ISUB = 1,NSUBSYS

         ! . Write out the subsystem name.
         LENGTH = LEN_TRIM ( SUBNAM(ISUB) )
         WRITE ( OUTPUT, "(A,I4,2X,A)" ) "Subsystem ", ISUB, SUBNAM(ISUB)(1:LENGTH)//":"

         ! . Calculate the number of residues in the subsystem.
         NRES = SUBIND(ISUB+1) - SUBIND(ISUB)

         ! . Allocate and initialize some temporary arrays.
         ALLOCATE ( FREQUENCY(1:NRES), NAME(1:NRES) )
         FREQUENCY = 0

         ! . Loop over the residues in the subsystem.
         N = 0
         DO I = (SUBIND(ISUB)+1),SUBIND(ISUB+1)

            ! . Search for a name match with previously stored residues.
            DO J = 1,N
               IF ( RESNAM(I) == NAME(J) ) THEN
                  FREQUENCY(J) = FREQUENCY(J) + 1
                  GO TO 10
               END IF
            END DO

            ! . There is no match.
            N = N + 1
            FREQUENCY(N) = 1
            NAME(N)      = RESNAM(I)

            ! . End of the loop.
            10 CONTINUE

         END DO

         ! . Write out the residue names and their frequencies.
         WRITE ( OUTPUT, "(2(2X,A,I6))" ) ( NAME(I), FREQUENCY(I), I = 1,N )

         ! . Deallocate the temporary arrays.
         DEALLOCATE ( FREQUENCY, NAME )

      END DO

      ! . Write out the terminator.
      WRITE ( OUTPUT, "(80('-'))" )

   END IF

   END SUBROUTINE SEQUENCE_SUMMARY

   !---------------------------------
   SUBROUTINE SEQUENCE_WRITE ( FILE )
   !---------------------------------

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: FILE

   ! . Local scalars.
   CHARACTER ( LEN = LINE_LENGTH ) :: LINE
   INTEGER                         :: II, IOSTAT, IRES, ISUB, LENGTH, N, START, STOP, UNIT

   ! . Check to see if the FILE argument is present.
   IF ( PRESENT ( FILE ) ) THEN

      ! . Get the next unit number.
      UNIT = NEXT_UNIT ( )

      ! . Open the file.
      OPEN ( UNIT, FILE = FILE, ACTION = "WRITE", STATUS = "REPLACE", IOSTAT = IOSTAT )

      ! . If there has been an error return.
      IF ( IOSTAT /= 0 ) CALL ERROR ( "SEQUENCE_WRITE", "I/O Error.", IOSTAT )

      ! . Do some printing.
      WRITE ( OUTPUT, "(/A)" ) "Sequence written to "//FILE

   ! . The unit for writing is the output stream.
   ELSE

      ! . Assign the unit number.
      UNIT = OUTPUT

      ! . Do some printing.
      WRITE ( OUTPUT, "(/A/)" ) "Sequence written to the output stream."

   END IF

   ! . Write out the sequence file header.
   WRITE ( UNIT, "(A)" ) "Sequence"

   ! . Write out the number of subsystems.
   WRITE ( UNIT, "(I6)" ) NSUBSYS

   ! . Loop over the subsystems.
   DO ISUB = 1,NSUBSYS

      ! . Write out the subsystem header.
      WRITE ( UNIT, "(/A)" ) "Subsystem "//TRIM ( SUBNAM(ISUB) )

      ! . Write out the number of residues in the subsystem.
      WRITE ( UNIT, "(I6)" ) ( SUBIND(ISUB+1) - SUBIND(ISUB) )

      ! . Find the longest residue name.
      LENGTH = 0
      DO IRES = (SUBIND(ISUB)+1),SUBIND(ISUB+1)
         LENGTH = MAX ( LENGTH, LEN_TRIM ( RESNAM(IRES) ) )
      END DO

      ! . Calculate the number of residue names that can be fit onto a line.
      N = ( LINE_LENGTH + 3 ) / ( LENGTH + 3 )

      ! . Top of the loop over residues.
      START = SUBIND(ISUB) + 1
      10 CONTINUE

         ! . Calculate the stopping position.
	 STOP = MIN ( ( START + N - 1 ), SUBIND(ISUB+1) )

         ! . Blank out the line.
	 LINE = REPEAT ( " ", LINE_LENGTH )

         ! . Fill the line.
	 DO IRES = START,STOP
	    II = ( IRES - START ) * ( LENGTH + 3 )
	    LINE(II+1:II+LENGTH) = RESNAM(IRES)(1:LENGTH)
	    IF ( IRES /= STOP ) LINE(II+LENGTH+1:II+LENGTH+3) = " ; "
	 END DO

         ! . Write out the lines.
         WRITE ( UNIT, "(A)" ) TRIM ( LINE )

         ! . Reset START.
	 START = STOP + 1

      ! . Check for termination.
      IF ( START <= SUBIND(ISUB+1) ) GO TO 10

      ! . Write out the subsystem terminator.
      WRITE ( UNIT, "(A)" ) "End"

   END DO

   ! . Write out the sequence terminator.
   WRITE ( UNIT, "(/A)" ) "End"

   ! . Close the file.
   IF ( UNIT /= OUTPUT ) CLOSE ( UNIT )

   END SUBROUTINE SEQUENCE_WRITE

END MODULE SEQUENCE
