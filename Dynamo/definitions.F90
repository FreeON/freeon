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
!                             The Definitions Module
!===============================================================================
!
! . Parameter data:
!
!   DEFAULT_INPUT                   The default input stream.
!   DEFAULT_OUTPUT                  The default output stream.
!   DP                              The Default Precision for reals.
!   FORCE_FIELD                     The module library force field.
!   LINE_LENGTH                     The maximum line length.
!   MAX_RECORD_LENGTH               The maximum binary file record length.
!   MAX_UNITS                       The maximum number of files.
!   SP                              Single Precision.
!   VERSION                         The version number of the program.
!
! . Notes:
!
!   This version has been tested with the NAG compiler on PC/Linux and HP/HP-UX
!   machines and the DEC compiler on DEC Alpha systems.
!
!   A 32 bit precision is sufficient for integers. For real numbers a 64 bit
!   precision is used throughout.
!
!===============================================================================
MODULE DEFINITIONS

! . F90 module declarations.
!USE F90_KIND, ONLY : DOUBLE, SINGLE ! NAG compiler.

IMPLICIT NONE
PUBLIC

! . The default program input stream.
INTEGER, PARAMETER :: DEFAULT_INPUT = 5

! . The default program output stream.
INTEGER, PARAMETER :: DEFAULT_OUTPUT = 6

! . The default precision parameter for real numbers.
!INTEGER, PARAMETER :: DP = DOUBLE ! NAG compiler.
INTEGER, PARAMETER :: DP = KIND(0.D0) ! 8      ! Cray T3E, DEC, PG and SGI compilers.

! . The module library force field.
CHARACTER ( LEN = 16 ), PARAMETER :: FORCE_FIELD = "OPLS_AA"

! . The input line length.
INTEGER, PARAMETER :: LINE_LENGTH = 132

! . The maximum record length for a binary file.
INTEGER, PARAMETER :: MAX_RECORD_LENGTH = 2**20 ! Linux.

! . The maximum FORTRAN stream number.
INTEGER, PARAMETER :: MAX_UNITS = 100

! . The single precision parameter for real numbers.
!INTEGER, PARAMETER :: SP = SINGLE ! NAG compiler.
INTEGER, PARAMETER :: SP = KIND(0.0) !4      ! Cray T3E, DEC, PG and SGI compilers.

! . The version number of the module library.
REAL ( KIND = DP ), PARAMETER :: VERSION = 2.0_DP

END MODULE DEFINITIONS
