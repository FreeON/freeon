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
!                              The Symmetry Module
!===============================================================================
!
! . Parameter data:
!
!   BOX_NAME_LENGTH                 The length of the box type name.
!
! . Scalar data:
!
!   BOX_TYPE                        The type of periodic box.
!   QBOX                            The flag to indicate if a box exists.
!
! . Array data:
!
!   BOXL                            The sides of the rectangular box.
!
! . Subroutines:
!
!   SYMMETRY_CUBIC_BOX              Define a cubic periodic box.
!   SYMMETRY_INITIALIZE             Initialize the symmetry variables.
!   SYMMETRY_ORTHORHOMBIC_BOX       Define an orthorhombic periodic box.
!   SYMMETRY_RECORD_READ            Read a symmetry record.
!   SYMMETRY_RECORD_WRITE           Write a symmetry record.
!   SYMMETRY_SUMMARY                Write out a summary of the symmetry data.
!
! . Notes:
!
!   The symmetry module defines the symmetry of a system. In this version
!   there are only a limited number of options. It is only possible to use
!   translational symmetry and, at present, only two possibilities are
!   allowed - cubic and orthorhombic systems. These are particularly useful
!   for generating periodic systems for use with periodic boundary conditions.
!
!===============================================================================
MODULE SYMMETRY

! . Module declarations.
USE DEFINITIONS, ONLY : DP
USE IO_UNITS,    ONLY : OUTPUT
USE PARSING
USE STATUS,      ONLY : ERROR

IMPLICIT NONE
PUBLIC
SAVE

! . Parameter data.
INTEGER, PARAMETER :: BOX_NAME_LENGTH = 16

! . Scalar data.
CHARACTER ( LEN = BOX_NAME_LENGTH ) :: BOX_TYPE = " "
LOGICAL                             :: QBOX     = .FALSE.

! . Array data.
REAL ( KIND = DP ), DIMENSION(1:3) :: BOXL = 0.0_DP

!==============================================================================
CONTAINS
!==============================================================================

   !----------------------------------
   SUBROUTINE SYMMETRY_CUBIC_BOX ( A )
   !----------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN) :: A

   ! . Initialize the symmetry variables.
   CALL SYMMETRY_INITIALIZE

   ! . Set the box parameters.
   BOXL = ABS ( A ) ; BOX_TYPE = "CUBIC" ; QBOX = .TRUE.

   ! . Write out the symmetry data.
   CALL SYMMETRY_SUMMARY

   END SUBROUTINE SYMMETRY_CUBIC_BOX

   !-----------------------------
   SUBROUTINE SYMMETRY_INITIALIZE
   !-----------------------------

   ! . Reset the box parameters.
   BOXL = 0.0_DP ; BOX_TYPE = " " ; QBOX = .FALSE.

   END SUBROUTINE SYMMETRY_INITIALIZE

   !-----------------------------------------------
   SUBROUTINE SYMMETRY_ORTHORHOMBIC_BOX ( A, B, C )
   !-----------------------------------------------

   ! . Scalar arguments.
   REAL ( KIND = DP ), INTENT(IN) :: A, B, C


   ! . Initialize the symmetry variables.
   CALL SYMMETRY_INITIALIZE

   ! . Set the box parameters.
   BOXL(1) = ABS ( A ) ; BOXL(2) = ABS ( B ) ; BOXL(3) = ABS ( C ) ; BOX_TYPE = "ORTHORHOMBIC" ; QBOX = .TRUE.

   ! . Write out the data.
   CALL SYMMETRY_SUMMARY

   END SUBROUTINE SYMMETRY_ORTHORHOMBIC_BOX

   !------------------------------------------
   SUBROUTINE SYMMETRY_RECORD_READ ( COMPACT )
   !------------------------------------------

   ! . Scalar arguments.
   LOGICAL, INTENT(IN), OPTIONAL :: COMPACT

   ! . Local scalars.
   CHARACTER ( LEN = LEN ( BOX_TYPE ) ) :: TMPBOX
   INTEGER                              :: I, ICARD
   LOGICAL                              :: QCOMPACT

   ! . Get the first word of the next line.
   CALL GET_LINE
   CALL GET_WORD

   ! . Check for symmetry information.
   IF ( WORD(1:WRDLEN) == "SYMMETRY" ) THEN

      ! . Check the format.
      IF ( PRESENT ( COMPACT ) ) THEN
         QCOMPACT = COMPACT
      ELSE
         QCOMPACT = .FALSE.
      END IF

      ! . Full format.
      IF ( .NOT. QCOMPACT ) THEN

         ! . Get the rest of the information on the line.
         CALL GET_INTEGER ( ICARD )

         ! . Check the number of cards.
         IF ( ICARD /= 1 ) CALL PARSE_ERROR ( "SYMMETRY_RECORD_READ", "Invalid number of symmetry cards." )

         ! . Get the first word on the next line.
         CALL GET_LINE

      END IF

      ! . Interpret the next word on the line.
      CALL GET_WORD ; TMPBOX = WORD(1:WRDLEN)

      ! . Branch on the box type.
      SELECT CASE ( TMPBOX )
      CASE ( "CUBIC           " ) ; CALL GET_REAL ( BOXL(1) ) ; BOXL(2:3) = BOXL(1)
      CASE ( "ORTHORHOMBIC    " ) ; DO I = 1,3 ; CALL GET_REAL ( BOXL(I) ) ; END DO
      CASE DEFAULT ; CALL ERROR ( "SYMMETRY_RECORD_READ", "Unknown periodic box type." )
      END SELECT

      ! . Set some symmetry variables.
      BOX_TYPE = TMPBOX
      QBOX     = .TRUE.

      ! . Reposition the file correctly by getting the first word on the next line.
      CALL GET_LINE
      CALL GET_WORD

      ! . Print a summary of the symmetry data.
      CALL SYMMETRY_SUMMARY

   END IF

   END SUBROUTINE SYMMETRY_RECORD_READ

   !-------------------------------------------------
   SUBROUTINE SYMMETRY_RECORD_WRITE ( UNIT, COMPACT )
   !-------------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN)           :: UNIT
   LOGICAL, INTENT(IN), OPTIONAL :: COMPACT

   ! . Local scalars.
   INTEGER :: ICARD
   LOGICAL :: QCOMPACT

   ! . Check the format.
   IF ( PRESENT ( COMPACT ) ) THEN
      QCOMPACT = COMPACT
   ELSE
      QCOMPACT = .FALSE.
   END IF

   ! . Compact format.
   IF ( QCOMPACT ) THEN

      ! . There is symmetry data.
      IF ( QBOX ) THEN

         ! . Write out the periodic box information.
         SELECT CASE ( BOX_TYPE )
         CASE ( "CUBIC           " ) ; WRITE ( UNIT, "(A,F18.10)"  ) "SYMMETRY "//BOX_TYPE, BOXL(1)
         CASE ( "ORTHORHOMBIC    " ) ; WRITE ( UNIT, "(A,3F18.10)" ) "SYMMETRY "//BOX_TYPE, BOXL(1:3)
         CASE DEFAULT ; CALL ERROR ( "SYMMETRY_RECORD_WRITE", "Unknown periodic box type." )
         END SELECT

      ! . There is no symmetry data so write a blank line.
      ELSE
         WRITE ( UNIT, "(1X)" )
      END IF

   ! . Full format.
   ELSE

      ! . There is symmetry data.
      IF ( QBOX ) THEN

         ! . Initialize the number of cards.
         ICARD = 1

         ! . Write out the symmetry header.
         WRITE ( UNIT, "('!',79('='))" )
         WRITE ( UNIT, "(A,I6)" ) "Symmetry", ICARD

         ! . Write out the periodic box information.
         SELECT CASE ( BOX_TYPE )
         CASE ( "CUBIC           " ) ; WRITE ( UNIT, "(A,F18.10)"  ) BOX_TYPE, BOXL(1)
         CASE ( "ORTHORHOMBIC    " ) ; WRITE ( UNIT, "(A,3F18.10)" ) BOX_TYPE, BOXL(1:3)
         CASE DEFAULT ; CALL ERROR ( "SYMMETRY_RECORD_WRITE", "Unknown periodic box type." )
         END SELECT

      END IF
   END IF

   END SUBROUTINE SYMMETRY_RECORD_WRITE

   !--------------------------
   SUBROUTINE SYMMETRY_SUMMARY
   !--------------------------

   ! . Local scalars.
   REAL ( KIND = DP ) :: VOLUME

   ! . Check that there is symmetry.
   IF ( QBOX ) THEN

      ! . Calculate the volume.
      VOLUME = PRODUCT ( BOXL )

      ! . Branch on the box type.
      SELECT CASE ( BOX_TYPE )
      CASE ( "CUBIC           " )
         WRITE ( OUTPUT, "(/27('-'),A,26('-'))"  ) " Summary of Cubic Box Data "
         WRITE ( OUTPUT, "(A,F14.4,2X,A,F14.4)" ) "Box Length             = ", BOXL(1), &
                                                  "Box Volume             = ", VOLUME
      CASE ( "ORTHORHOMBIC    " )
         WRITE ( OUTPUT, "(/23('-'),A,23('-'))"  ) " Summary of Orthorhombic Box Data "
         WRITE ( OUTPUT, "(A,F14.4,2X,A,F14.4)" ) "Box Length - X         = ", BOXL(1), &
                                                  "Box Length - Y         = ", BOXL(2)
         WRITE ( OUTPUT, "(A,F14.4,2X,A,F14.4)" ) "Box Length - Z         = ", BOXL(3), &
                                                  "Box Volume             = ", VOLUME
      END SELECT

      ! . Write out the terminator.
      WRITE ( OUTPUT, "(80('-'))" )

   ! . There is no symmetry.
   ELSE

      WRITE ( OUTPUT, "(/A)"  ) "No symmetry is defined."

   END IF

   END SUBROUTINE SYMMETRY_SUMMARY

END MODULE SYMMETRY
