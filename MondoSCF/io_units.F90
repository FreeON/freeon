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
!                              The I/O Units Module
!===============================================================================
!
! . I/O Channels:
!
!   INPUT                          The input  stream unit number.
!   OUTPUT                         The output stream unit number.
!
! . Notes:
!
!   The I/O unit values are parameters and so cannot be changed.
!
!===============================================================================
MODULE IO_UNITS

! . Module declarations.
USE DEFINITIONS, ONLY : DEFAULT_INPUT, DEFAULT_OUTPUT

IMPLICIT NONE
PUBLIC

! . Integer data.
INTEGER, PARAMETER :: INPUT = DEFAULT_INPUT, OUTPUT = DEFAULT_OUTPUT

END MODULE IO_UNITS
