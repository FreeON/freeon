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
!                                The Files Module
!===============================================================================
!
! . Functions:
!
!   NEXT_UNIT                      Get the next available FORTRAN stream number.
!
! . Notes:
!
!   The Fortran 90 INQUIRE statement is used to find the next unit which is not
!   open. The input and output streams are always assumed to be assigned and
!   are not checked.
!
!===============================================================================
MODULE FILES

! . Module declarations.
USE DEFINITIONS, ONLY : MAX_UNITS
USE IO_UNITS,    ONLY : INPUT, OUTPUT

IMPLICIT NONE
PUBLIC

!==============================================================================
CONTAINS
!==============================================================================

   !---------------------
   FUNCTION NEXT_UNIT ( )
   !---------------------

   ! . Function declarations.
   INTEGER :: NEXT_UNIT

   ! . Local scalars.
   INTEGER :: I, UNIT
   LOGICAL :: QOPEN

   ! . Initialize the unit number.
   UNIT = -1

   ! . Loop over the number of units available.
   DO I = 1,MAX_UNITS

      ! . Check all units except for the standard input and output streams.
      IF ( ( I /= INPUT ) .AND. ( I /= OUTPUT ) ) THEN

         ! . Inquire if the unit is connected.
         INQUIRE ( I, OPENED = QOPEN )

         ! . If the unit is unconnected assign the unit and exit the loop.
         IF ( .NOT. QOPEN ) THEN
            UNIT = I ; EXIT
         END IF

      END IF

   END DO

   ! . Assign a value to the function.
   NEXT_UNIT = UNIT

   END FUNCTION NEXT_UNIT

END MODULE FILES
