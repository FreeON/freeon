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
!                                The Atoms Module
!===============================================================================
!
! . Parameter data:
!
!   ATOM_NAME_LENGTH                The maximum length of an atom name.
!
! . Scalar data:
!
!   MM_NATOMS                          The total number of atoms.
!   NFIXED                          The number of fixed atoms.
!   NFREE                           The number of free atoms (MM_NATOMS - NFIXED).
!
!   MM_NATOMSMM                        The number of MM atoms.
!   MM_NATOMSQM                        The number of QM atoms.
!
! . Array data:
!
!   ATMCRD                          The coordinates of the atoms.
!   ATMFIX                          The array to indicate which atoms are fixed.
!   ATMIND                          The free atom index array.
!   ATMMAS                          The atom masses.
!   ATMNAM                          The atom names.
!   ATMNUM                          The atom numbers.
!   ATMQMI                          The indices of the QM atoms.
!
! . Subroutines:
!
!   ATOMS_ALLOCATE                  Allocate the atoms variables.
!   ATOMS_FIX                       Fill the ATMFIX array.
!   ATOMS_FORMULA                   Write out the formula for the structure.
!   ATOMS_INITIALIZE                Initialize the atoms variables.
!   ATOMS_SUMMARY                   Print a summary of the atom data.
!
! . Notes:
!
!   The atoms module stores the data for the atoms in the system. In addition
!   to the number of atoms in the system, for each atom there are its mass,
!   its name and its coordinates. All data is public and, hence, accessible.
!
!   All atoms are supposed to be MM by default.
!
!===============================================================================
MODULE ATOMS

! . Module declarations.
USE CONSTANTS,   ONLY : UNDEFINED
USE DEFINITIONS, ONLY : DP
USE ELEMENTS,    ONLY : NELEMENTS, SYMBOL
USE IO_UNITS,    ONLY : OUTPUT
!USE DYNAMO_VARIABLES

IMPLICIT NONE
PUBLIC
SAVE

! . Parameter declarations.
INTEGER, PARAMETER :: ATOM_NAME_LENGTH = 8

! . Scalar data.
INTEGER :: MM_NATOMS = 0, MM_NATOMSMM = 0, MM_NATOMSQM = 0, NFIXED = 0, NFREE = 0

! . Character array data.
CHARACTER ( LEN = ATOM_NAME_LENGTH ), ALLOCATABLE, DIMENSION(:) :: ATMNAM

! . Integer array data.
INTEGER, ALLOCATABLE, DIMENSION(:) :: ATMIND, ATMNUM, ATMQMI

! . Logical array data.
LOGICAL, ALLOCATABLE, DIMENSION(:) :: ATMFIX

! . Real array data.
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:)   :: ATMMAS
REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: ATMCRD

!==============================================================================
CONTAINS
!==============================================================================

   !------------------------------
   SUBROUTINE ATOMS_ALLOCATE ( N )
   !------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: N

   ! . Local scalars.
   INTEGER :: I

   ! . Initialize the scalar variables.
   MM_NATOMS   = N
   MM_NATOMSMM = N
   MM_NATOMSQM = 0
   NFIXED   = 0
   NFREE    = N

   ! . Allocate all the arrays.
   ALLOCATE ( ATMCRD(1:3,1:MM_NATOMS), ATMFIX(1:MM_NATOMS), ATMIND(1:MM_NATOMS), ATMMAS(1:MM_NATOMS), ATMNAM(1:MM_NATOMS), &
              ATMNUM(1:MM_NATOMS),     ATMQMI(1:MM_NATOMS) )

   ! . Initialize the ATMCRD, ATMFIX and ATMQMI arrays.
   ATMCRD = UNDEFINED
   ATMFIX = .FALSE.
   ATMQMI = 0

   ! . Initialize ATMIND.
   DO I = 1,MM_NATOMS
      ATMIND(I) = I
   END DO

   END SUBROUTINE ATOMS_ALLOCATE

   !----------------------------
   SUBROUTINE ATOMS_FIX ( QFIX )
   !----------------------------

   ! . Array arguments.
   LOGICAL, DIMENSION(1:MM_NATOMS), INTENT(IN) :: QFIX

   ! . Local scalars.
   INTEGER :: I, IFREE

   ! . There are atoms.
   IF ( MM_NATOMS > 0 ) THEN

      ! . Copy the QFIX array to ATMFIX.
      ATMFIX = QFIX

      ! . Set the number of fixed atoms.
      NFIXED = COUNT ( ATMFIX )

      ! . Set the number of free atoms
      NFREE = MM_NATOMS - NFIXED

      ! . Fill the ATMIND array.
      IFREE = 0
      DO I = 1,MM_NATOMS
         IF ( ATMFIX(I) ) THEN
	    ATMIND(I) = 0
	 ELSE
	    IFREE = IFREE + 1
	    ATMIND(I) = IFREE
	 END IF
      END DO

   ! . There are no atoms.
   ELSE
      NFIXED = 0
   END IF

   ! . Write out some data.
   WRITE ( OUTPUT, "(/A,I6,A)" ) "The coordinates of ", NFIXED, " atoms have been fixed."

   END SUBROUTINE ATOMS_FIX

   !-----------------------
   SUBROUTINE ATOMS_FORMULA
   !-----------------------

   ! . Local scalars.
   INTEGER :: I, N, NOTHER

   ! . Local arrays.
   INTEGER, ALLOCATABLE, DIMENSION(:) :: FREQUENCY, INDEX

   ! . Check that there are atoms.
   IF ( MM_NATOMS > 0 ) THEN

      ! . Allocate and initialize some temporary arrays.
      ALLOCATE ( FREQUENCY(1:NELEMENTS), INDEX(1:NELEMENTS) )
      FREQUENCY = 0

      ! . Loop over the atoms.
      NOTHER = 0
      DO I = 1,MM_NATOMS
         IF ( ( ATMNUM(I) <= 0 ) .OR. ( ATMNUM(I) > NELEMENTS ) ) THEN
            NOTHER = NOTHER + 1
         ELSE
            FREQUENCY(ATMNUM(I)) = FREQUENCY(ATMNUM(I)) + 1
         END IF
      END DO

      ! . Contract the frequency array and create the index array.
      N = 0
      DO I = 1,NELEMENTS
         IF ( FREQUENCY(I) > 0 ) THEN
            N = N + 1
            FREQUENCY(N) = FREQUENCY(I)
            INDEX(N)     = I
         END IF
      END DO

      ! . Write out the formula.
      WRITE ( OUTPUT, "(/33('-'),A,33('-'))"  ) " Atom Formula "
      WRITE ( OUTPUT, "(6(A2,I8,4X))" ) ( SYMBOL(INDEX(I)), FREQUENCY(I), I = 1,N )
      IF ( NOTHER > 0 ) THEN
         WRITE ( OUTPUT, "(/A,I12)" ) "Number of Unknown Atoms = ", NOTHER
      END IF
      WRITE ( OUTPUT, "(80('-'))" )

      ! . Deallocate the temporary arrays.
      DEALLOCATE ( FREQUENCY, INDEX )

   END IF

   END SUBROUTINE ATOMS_FORMULA

   !--------------------------
   SUBROUTINE ATOMS_INITIALIZE
   !--------------------------

   ! . Initialize the scalar variables.
   MM_NATOMS   = 0
   MM_NATOMSMM = 0
   MM_NATOMSQM = 0
   NFIXED   = 0
   NFREE    = 0

   ! . Deallocate all the arrays.
   IF ( ALLOCATED ( ATMCRD ) ) DEALLOCATE ( ATMCRD )
   IF ( ALLOCATED ( ATMFIX ) ) DEALLOCATE ( ATMFIX )
   IF ( ALLOCATED ( ATMIND ) ) DEALLOCATE ( ATMIND )
   IF ( ALLOCATED ( ATMMAS ) ) DEALLOCATE ( ATMMAS )
   IF ( ALLOCATED ( ATMNAM ) ) DEALLOCATE ( ATMNAM )
   IF ( ALLOCATED ( ATMNUM ) ) DEALLOCATE ( ATMNUM )
   IF ( ALLOCATED ( ATMQMI ) ) DEALLOCATE ( ATMQMI )

   END SUBROUTINE ATOMS_INITIALIZE

   !-----------------------
   SUBROUTINE ATOMS_SUMMARY
   !-----------------------

   ! . Local scalars.
   INTEGER            :: NHEAVY, NHYDRO, NOTHER
   REAL ( KIND = DP ) :: TOTMAS

   ! . Check that there are some atoms.
   IF ( MM_NATOMS > 0 ) THEN

      ! . Determine the number of atoms in each class and the total mass.
      NHYDRO = COUNT ( ATMNUM(1:MM_NATOMS) == 1 )
      NOTHER = COUNT ( ATMNUM(1:MM_NATOMS) <= 0 ) + COUNT ( ATMNUM(1:MM_NATOMS) > NELEMENTS )
      NHEAVY = MM_NATOMS - NHYDRO - NOTHER
      TOTMAS = SUM ( ATMMAS(1:MM_NATOMS) )

      ! . Write out the data.
      WRITE ( OUTPUT, "(/29('-'),A,29('-'))"  ) " Summary of Atom Data "
      WRITE ( OUTPUT, "(A,I14,2X,A,I14)"   ) "Number of Atoms        = ", MM_NATOMS, &
                                             "Number of Hydrogens    = ", NHYDRO
      WRITE ( OUTPUT, "(A,I14,2X,A,I14)"   ) "Number of Heavy Atoms  = ", NHEAVY, &
                                             "Number of Unknowns     = ", NOTHER
      WRITE ( OUTPUT, "(A,I14,2X,A,F14.1)" ) "Number of Fixed Atoms  = ", NFIXED, &
                                             "Total Mass (a.m.u.)    = ", TOTMAS
      WRITE ( OUTPUT, "(80('-'))" )

   END IF

   END SUBROUTINE ATOMS_SUMMARY

END MODULE ATOMS
