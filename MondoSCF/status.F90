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
!                               The Status Module
!===============================================================================
!
! . Subroutines:
!
!   ERROR                    Write out an error message and stop the execution
!                            of the program.
!
!===============================================================================
MODULE STATUS

! . Module declarations.
USE IO_UNITS, ONLY : OUTPUT

IMPLICIT NONE
PUBLIC

!==============================================================================
CONTAINS
!==============================================================================

   !------------------------------------------
   SUBROUTINE ERROR ( ROUTINE, MESSAGE, CODE )
   !------------------------------------------

   ! . Essential scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: MESSAGE, ROUTINE

   ! . Optional scalar arguments.
   INTEGER, INTENT(IN), OPTIONAL :: CODE

   ! . Local scalars.
   INTEGER :: LENM, LENR

   ! . Find the lengths of the character strings.
   LENM = LEN_TRIM ( MESSAGE )
   LENR = LEN_TRIM ( ROUTINE )

   ! . Write out the message.
   WRITE ( OUTPUT, "(/A)" ) "Error in " // ROUTINE(1:LENR) // ":"
   WRITE ( OUTPUT, "(A)"  ) MESSAGE(1:LENM)
   IF ( PRESENT ( CODE ) ) WRITE ( OUTPUT, "(A,I6)" ) "Error code = ", CODE

   ! . Terminate the program.
   WRITE ( OUTPUT, "(/A)" ) "Program terminating."
   STOP
   
   END SUBROUTINE ERROR

END MODULE STATUS
