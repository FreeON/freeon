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
!                            The Connectivity Module
!===============================================================================
!
! . Subroutines:
!
!   CONNECTIVITY_ANGLES            Generate the angle connectivity list.
!   CONNECTIVITY_BONDS             Generate the bond connectivity list.
!   CONNECTIVITY_DIHEDRALS         Generate the dihedral connectivity list.
!   CONVERT_ALL_TO_1234            Convert all arrays to 1-2,1-3,1-4 lists.
!   CONVERT_ALL_TO_14              Convert all arrays to 1-4 lists.
!   CONVERT_BONDS_TO_12            Convert BONDS to 1-2 lists.
!   CONVERT_12_TO_BONDS            Convert 1-2 to BONDS lists.
!   GENERATE_12                    Generate 1-2 lists.
!   ORDER_BONDS                    Order a bond list.
!   REMOVE_BONDS_REDUNDANCIES      Remove redundancies from a bond list.
!
!===============================================================================
MODULE CONNECTIVITY

! . Module declarations.
USE DEFINITIONS, ONLY : DP
USE ELEMENTS,    ONLY : RADII
USE IO_UNITS,    ONLY : OUTPUT
USE SORT_dyn

USE ATOMS
USE SEQUENCE

IMPLICIT NONE
PUBLIC
PRIVATE :: CONVERT_BONDS_TO_12, CONVERT_12_TO_BONDS, REMOVE_BONDS_REDUNDANCIES, REMOVE_123_REDUNDANCIES

!==============================================================================
CONTAINS
!==============================================================================

   !-----------------------------------------------
   SUBROUTINE CONNECTIVITY_ANGLES ( ANGLES, BONDS )
   !-----------------------------------------------

   ! . Array argument declarations.
   INTEGER, DIMENSION(:,:), POINTER    :: ANGLES
   INTEGER, DIMENSION(:,:), INTENT(IN) :: BONDS

   ! . Local scalars.
   INTEGER :: I, IBOND, J, JBOND, K, MAXANG1, MAXANG2, NANGLE1, NANGLE2, NANGLES, NBONDS

   ! . Local arrays.
   INTEGER, ALLOCATABLE, DIMENSION(:)            :: LISTI, LISTJ
   INTEGER,              DIMENSION(:,:), POINTER :: ANGLE1, ANGLE2

   ! . Get the number of bonds.
   NBONDS = SIZE ( BONDS, 2 )

   ! . Check the array dimensions.
   IF ( ( SIZE ( BONDS, 1 ) /= 2 ) .OR. ( NBONDS <= 1 ) ) THEN
      ALLOCATE ( ANGLES(1:3,1:0) ) ; RETURN
   END IF

   ! . Allocate some temporary arrays.
   ALLOCATE ( LISTI(1:MM_NATOMS+1), LISTJ(1:NBONDS) )

   ! . Convert the BONDS list to a 1-2 list.
   CALL CONVERT_BONDS_TO_12 ( BONDS, LISTI, LISTJ )

   ! . Guess the maximum number of angles.
   MAXANG1 = 4 * NBONDS
   MAXANG2 = 4 * NBONDS

   ! . Allocate the temporary angle arrays.
   ALLOCATE ( ANGLE1(1:3,1:MAXANG1), ANGLE2(1:3,1:MAXANG2) )

   ! . Initialize the number of angles.
   NANGLE1 = 0
   NANGLE2 = 0

   ! . Loop over the atoms.
   DO I = 1,MM_NATOMS

      ! . Loop over the bonds for atom I.
      DO IBOND = (LISTI(I)+1),LISTI(I+1)

         ! . Get the atom index.
         J = LISTJ(IBOND)

         ! . Loop over the bonds for atom J with atoms of lower index.
         DO K = (I+1),(J-1)
            DO JBOND = (LISTI(K)+1),LISTI(K+1)
               IF ( LISTJ(JBOND) == J ) THEN

                  ! . Check NANGLE1.
                  IF ( NANGLE1 == MAXANG1 ) CALL RESIZE ( NANGLE1, MAXANG1, ANGLE1 )

                  ! . Save the angle.
                  NANGLE1 = NANGLE1 + 1
                  ANGLE1(1,NANGLE1) = I
                  ANGLE1(2,NANGLE1) = J
                  ANGLE1(3,NANGLE1) = K

               END IF
            END DO
         END DO

         ! . Loop over the bonds for atom J with atoms of higher index.
         DO JBOND = (LISTI(J)+1),LISTI(J+1)

           ! . Get the atom index.
           K = LISTJ(JBOND)

           ! . Check NANGLE1.
           IF ( NANGLE1 == MAXANG1 ) CALL RESIZE ( NANGLE1, MAXANG1, ANGLE1 )

           ! . Save the angle.
           NANGLE1 = NANGLE1 + 1
           ANGLE1(1,NANGLE1) = I
           ANGLE1(2,NANGLE1) = J
           ANGLE1(3,NANGLE1) = K

         END DO

         ! . Loop over the remaining bonds for atom I.
         DO JBOND = (IBOND+1),LISTI(I+1)

           ! . Get the atom index.
           K = LISTJ(JBOND)

           ! . Check NANGLE2.
           IF ( NANGLE2 == MAXANG2 ) CALL RESIZE ( NANGLE2, MAXANG2, ANGLE2 )

           ! . Save the angle.
           NANGLE2 = NANGLE2 + 1
           ANGLE2(1,NANGLE2) = J
           ANGLE2(2,NANGLE2) = I
           ANGLE2(3,NANGLE2) = K

         END DO
      END DO
   END DO

   ! . Deallocate some temporary arrays.
   DEALLOCATE ( LISTI, LISTJ )

   ! . Calculate the total number of angles.
   NANGLES = NANGLE1 + NANGLE2

   ! . Allocate space for the ANGLES array.
   ALLOCATE ( ANGLES(1:3,1:NANGLES) )

   ! . Save the angles.
   ANGLES(1:3,        1:NANGLE1) = ANGLE1(1:3,1:NANGLE1)
   ANGLES(1:3,NANGLE1+1:NANGLES) = ANGLE2(1:3,1:NANGLE2)

   ! . Do some printing.
   WRITE ( OUTPUT, "(/24('-'),A,24('-'))" ) " Angle Connectivity Information "
   WRITE ( OUTPUT, "(A,I14,2X,A,I14)" ) "Number of Angles       = ", NANGLES, &
                                        "Number of Connections  = ", NBONDS
   WRITE ( OUTPUT, "(80('-'))" )

   ! . Deallocate the remaining temporary arrays.
   DEALLOCATE ( ANGLE1, ANGLE2 )

   !===========================================================================
   CONTAINS
   !===========================================================================

      !------------------------------------------
      SUBROUTINE RESIZE ( NANGLE, MAXANG, ANGLE )
      !------------------------------------------

      ! . Scalar argument declarations.
      INTEGER, INTENT(INOUT) :: MAXANG
      INTEGER, INTENT(IN)    :: NANGLE

      ! . Array argument declarations.
      INTEGER, DIMENSION(:,:), POINTER :: ANGLE

      ! . Save the existing lists.
      ALLOCATE ( ANGLES(1:3,1:MAXANG) )
      ANGLES(1:3,1:MAXANG) = ANGLE(1:3,1:MAXANG)

      ! . Increment MAXANG.
      MAXANG = MAXANG + NBONDS

      ! . Reallocate the temporary lists and resave the data.
      DEALLOCATE ( ANGLE )
      ALLOCATE ( ANGLE(1:3,1:MAXANG) )
      ANGLE(1:3,1:NANGLE) = ANGLES(1:3,1:NANGLE)

      ! . Deallocate ANGLES.
      DEALLOCATE ( ANGLES )

      END SUBROUTINE RESIZE

   END SUBROUTINE CONNECTIVITY_ANGLES

   !----------------------------------------------
   SUBROUTINE CONNECTIVITY_BONDS ( BONDS, BUFFER )
   !----------------------------------------------

   ! . Scalar argument declarations.
   REAL ( KIND = DP ), INTENT(IN)  :: BUFFER

   ! . Array argument declarations.
   INTEGER, DIMENSION(:,:), POINTER :: BONDS

   ! . Local scalars.
   INTEGER :: NCONN

   ! . Local arrays.
   INTEGER, ALLOCATABLE, DIMENSION(:)          :: LISTI
   INTEGER,              DIMENSION(:), POINTER :: LISTJ

   ! . Allocate some temporary arrays.
   ALLOCATE ( LISTI(1:MM_NATOMS+1) )

   ! . Calculate the 1-2 lists.
   CALL GENERATE_12 ( LISTI, LISTJ, BUFFER )

   ! . Find out the number of connections.
   NCONN = SIZE ( LISTJ )

   ! . Allocate BONDS.
   ALLOCATE ( BONDS(1:2,1:NCONN) )

   ! . Convert between the lists.
   CALL CONVERT_12_TO_BONDS ( BONDS, LISTI, LISTJ )

   ! . Deallocate the temporary arrays.
   DEALLOCATE ( LISTI, LISTJ )

   END SUBROUTINE CONNECTIVITY_BONDS

   !-------------------------------------------------------------
   SUBROUTINE CONNECTIVITY_DIHEDRALS ( DIHEDRALS, BONDS, ANGLES )
   !-------------------------------------------------------------

   ! . Array argument declarations.
   INTEGER, DIMENSION(:,:), POINTER    :: DIHEDRALS
   INTEGER, DIMENSION(:,:), INTENT(IN) :: ANGLES, BONDS

   ! . Local scalars.
   INTEGER :: I, IANGLE, IBOND, J, K, L, MAXDIH1, MAXDIH2, NANGLES, NANGLE1, &
              NBONDS, NDIHEDRALS, NDIHEDRAL1, NDIHEDRAL2

   ! . Local arrays.
   INTEGER, ALLOCATABLE, DIMENSION(:)            :: LISTI, LISTJ
   INTEGER,              DIMENSION(:,:), POINTER :: DIHEDRAL1, DIHEDRAL2

   ! . Get the number of bonds and angles.
   NBONDS  = SIZE ( BONDS,  2 )
   NANGLES = SIZE ( ANGLES, 2 )

   ! . Check the number of bonds and angles.
   IF ( ( NBONDS <= 1 ) .OR. ( NANGLES <= 0 ) ) THEN
      ALLOCATE ( DIHEDRALS(1:4,1:0) ) ; RETURN
   END IF

   ! . Check the array dimensions.
   IF ( ( SIZE ( BONDS, 1 ) /= 2 ) .OR. ( SIZE ( ANGLES, 1 ) /= 3 ) ) RETURN

   ! . Allocate some temporary arrays.
   ALLOCATE ( LISTI(1:MM_NATOMS+1), LISTJ(1:NBONDS) )

   ! . Convert the BONDS list to a 1-2 list.
   CALL CONVERT_BONDS_TO_12 ( BONDS, LISTI, LISTJ )

   ! . Guess the maximum number of dihedrals.
   MAXDIH1 = 6 * NANGLES
   MAXDIH2 = 6 * NANGLES

   ! . Allocate the temporary dihedral arrays.
   ALLOCATE ( DIHEDRAL1(1:4,1:MAXDIH1), DIHEDRAL2(1:4,1:MAXDIH2) )

   ! . Initialize the number of dihedrals.
   NDIHEDRAL1 = 0
   NDIHEDRAL2 = 0

   ! . Find out how many angles of type 1 there are.
   DO I = 1,NANGLES
      IF ( ANGLES(1,I) > ANGLES(2,I) ) THEN
         NANGLE1 = I - 1 ; EXIT
      END IF
   END DO

   ! . Loop over the angles of type 1.
   DO IANGLE = 1,NANGLE1

      ! . Get the atom indices.
      I = ANGLES(1,IANGLE)
      J = ANGLES(2,IANGLE)
      K = ANGLES(3,IANGLE)

      ! . Loop over the bonds of lower index for atom K (I<J<K).
      IF ( K > J ) THEN

         ! . Loop over all lower atoms.
         DO L = 1,(K-1)

            ! . Skip if L is the same as I or J.
            IF ( ( L == I ) .OR. ( L == J ) ) CYCLE

            ! . Loop over the bonds for the atom.
            DO IBOND = (LISTI(L)+1),LISTI(L+1)
               IF ( LISTJ(IBOND) == K ) THEN

                  ! . Check NDIHEDRAL1.
                  IF ( NDIHEDRAL1 == MAXDIH1 ) CALL RESIZE ( NDIHEDRAL1, MAXDIH1, DIHEDRAL1 )

                  ! . Save the dihedral.
                  NDIHEDRAL1 = NDIHEDRAL1 + 1
                  DIHEDRAL1(1,NDIHEDRAL1) = I
                  DIHEDRAL1(2,NDIHEDRAL1) = J
                  DIHEDRAL1(3,NDIHEDRAL1) = K
                  DIHEDRAL1(4,NDIHEDRAL1) = L

               END IF
            END DO
         END DO
      END IF

      ! . Loop over the bonds of higher index for atom K.
      DO IBOND = (LISTI(K)+1),LISTI(K+1)

         ! . Get the atom index.
         L = LISTJ(IBOND)

         ! . Skip if L is the same as I or J.
         IF ( ( L == I ) .OR. ( L == J ) ) CYCLE

         ! . Check NDIHEDRAL1.
         IF ( NDIHEDRAL1 == MAXDIH1 ) CALL RESIZE ( NDIHEDRAL1, MAXDIH1, DIHEDRAL1 )

         ! . Save the dihedral.
         NDIHEDRAL1 = NDIHEDRAL1 + 1
         DIHEDRAL1(1,NDIHEDRAL1) = I
         DIHEDRAL1(2,NDIHEDRAL1) = J
         DIHEDRAL1(3,NDIHEDRAL1) = K
         DIHEDRAL1(4,NDIHEDRAL1) = L

      END DO

      ! . Loop over the bonds for atom I.
      DO IBOND = (LISTI(I)+1),LISTI(I+1)

         ! . Get the atom index.
         L = LISTJ(IBOND)

         ! . Skip if L is the same as J or K.
         IF ( ( L == J ) .OR. ( L == K ) ) CYCLE

         ! . Check NDIHEDRAL2.
         IF ( NDIHEDRAL2 == MAXDIH2 ) CALL RESIZE ( NDIHEDRAL2, MAXDIH2, DIHEDRAL2 )

         ! . Save the dihedral.
         NDIHEDRAL2 = NDIHEDRAL2 + 1
         DIHEDRAL2(1,NDIHEDRAL2) = L
         DIHEDRAL2(2,NDIHEDRAL2) = I
         DIHEDRAL2(3,NDIHEDRAL2) = J
         DIHEDRAL2(4,NDIHEDRAL2) = K

      END DO
   END DO

   ! . Deallocate some temporary arrays.
   DEALLOCATE ( LISTI, LISTJ )

   ! . Calculate the total number of dihedrals.
   NDIHEDRALS = NDIHEDRAL1 + NDIHEDRAL2

   ! . Allocate space for the DIHEDRALS array.
   ALLOCATE ( DIHEDRALS(1:4,1:NDIHEDRALS) )

   ! . Save the dihedrals.
   DIHEDRALS(1:4,           1:NDIHEDRAL1) = DIHEDRAL1(1:4,1:NDIHEDRAL1)
   DIHEDRALS(1:4,NDIHEDRAL1+1:NDIHEDRALS) = DIHEDRAL2(1:4,1:NDIHEDRAL2)

   ! . Do some printing.
   WRITE ( OUTPUT, "(/23('-'),A,22('-'))" ) " Dihedral Connectivity Information "
   WRITE ( OUTPUT, "(A,I14,2X,A,I14)" ) "Number of Dihedrals    = ", NDIHEDRALS, &
                                        "Number of Angles       = ", NANGLES
   WRITE ( OUTPUT, "(A,I14,2X,A,I14)" ) "Number of Connections  = ", NBONDS
   WRITE ( OUTPUT, "(80('-'))" )

   ! . Deallocate the remaining temporary arrays.
   DEALLOCATE ( DIHEDRAL1, DIHEDRAL2 )

   !===========================================================================
   CONTAINS
   !===========================================================================

      !------------------------------------------------
      SUBROUTINE RESIZE ( NDIHEDRAL, MAXDIH, DIHEDRAL )
      !------------------------------------------------

      ! . Scalar argument declarations.
      INTEGER, INTENT(INOUT) :: MAXDIH
      INTEGER, INTENT(IN)    :: NDIHEDRAL

      ! . Array argument declarations.
      INTEGER, DIMENSION(:,:), POINTER :: DIHEDRAL

      ! . Save the existing lists.
      ALLOCATE ( DIHEDRALS(1:4,1:MAXDIH) )
      DIHEDRALS(1:4,1:MAXDIH) = DIHEDRAL(1:4,1:MAXDIH)

      ! . Increment MAXDIH.
      MAXDIH = MAXDIH + NANGLES

      ! . Reallocate the temporary lists and resave the data.
      DEALLOCATE ( DIHEDRAL )
      ALLOCATE ( DIHEDRAL(1:4,1:MAXDIH) )
      DIHEDRAL(1:4,1:NDIHEDRAL) = DIHEDRALS(1:4,1:NDIHEDRAL)

      ! . Deallocate DIHEDRALS.
      DEALLOCATE ( DIHEDRALS )

      END SUBROUTINE RESIZE

   END SUBROUTINE CONNECTIVITY_DIHEDRALS

   !------------------------------------------------------------------------
   SUBROUTINE CONVERT_ALL_TO_1234 ( BONDS, ANGLES, DIHEDRALS, LISTI, LISTJ )
   !------------------------------------------------------------------------

   ! . Array argument declarations.
   INTEGER, DIMENSION(:,:),        INTENT(IN)  :: ANGLES, BONDS, DIHEDRALS
   INTEGER, DIMENSION(1:MM_NATOMS+1), INTENT(OUT) :: LISTI
   INTEGER, DIMENSION(:),              POINTER :: LISTJ

   ! . Local scalars.
   INTEGER :: NANGL, NBOND, NDIHE, NDIM

   ! . Local arrays.
   INTEGER, DIMENSION(:,:), POINTER :: TMPLIS

   ! . Get the array dimensions.
   NANGL = SIZE ( ANGLES,    2 )
   NBOND = SIZE ( BONDS,     2 )
   NDIHE = SIZE ( DIHEDRALS, 2 )

   ! . Check the array dimensions.
   IF ( ( SIZE ( BONDS, 1 ) /= 2 ) .OR. ( SIZE ( ANGLES, 1 ) /= 3 ) .OR. ( SIZE ( DIHEDRALS, 1 ) /= 4 ) ) THEN
      LISTI = - HUGE ( 0 ) ; ALLOCATE ( LISTJ(1:0) ) ; RETURN
   END IF

   ! . Calculate the total length of the temporary array.
   NDIM = NANGL + NBOND + NDIHE

   ! . Allocate the temporary array.
   ALLOCATE ( TMPLIS(1:2,1:NDIM) )

   ! . Copy BONDS to TMPLIS.
   TMPLIS(1:2,1:NBOND) = BONDS(1:2,1:NBOND)

   ! . Copy the first and third columns of ANGLES to TMPLIS.
   TMPLIS(1,NBOND+1:NBOND+NANGL) = ANGLES(1,1:NANGL)
   TMPLIS(2,NBOND+1:NBOND+NANGL) = ANGLES(3,1:NANGL)

   ! . Copy the first and fourth columns of DIHEDRALS to TMPLIS.
   TMPLIS(1,NDIM-NDIHE+1:NDIM) = DIHEDRALS(1,1:NDIHE)
   TMPLIS(2,NDIM-NDIHE+1:NDIM) = DIHEDRALS(4,1:NDIHE)

   ! . Reorder the temporary dihedral array.
   CALL ORDER_BONDS ( TMPLIS )

   ! . Remove any redundancies in the list.
   CALL REMOVE_BONDS_REDUNDANCIES ( TMPLIS )

   ! . Get the new dimension.
   NDIM = SIZE ( TMPLIS, 2 )

   ! . Allocate space for LISTJ.
   ALLOCATE ( LISTJ(1:NDIM) )

   ! . Convert the reordered array to 1-2 form.
   CALL CONVERT_BONDS_TO_12 ( TMPLIS, LISTI, LISTJ )

   ! . Deallocate the temporary array.
   DEALLOCATE ( TMPLIS )

   END SUBROUTINE CONVERT_ALL_TO_1234

   !-------------------------------------------------------------------------
   SUBROUTINE CONVERT_ALL_TO_14 ( BONDS, ANGLES, I1234, J1234, LISTI, LISTJ )
   !-------------------------------------------------------------------------

   ! . Array argument declarations.
   INTEGER, DIMENSION(:,:),        INTENT(IN)  :: ANGLES, BONDS
   INTEGER, DIMENSION(1:MM_NATOMS+1), INTENT(IN)  :: I1234
   INTEGER, DIMENSION(:),          INTENT(IN)  :: J1234
   INTEGER, DIMENSION(1:MM_NATOMS+1), INTENT(OUT) :: LISTI
   INTEGER, DIMENSION(:),              POINTER :: LISTJ

   ! . Local scalars.
   INTEGER :: NANGL, NBOND, NDIM

   ! . Local arrays.
   INTEGER, DIMENSION(:),   POINTER :: I123, J123
   INTEGER, DIMENSION(:,:), POINTER :: TMPLIS

   ! . Get the array dimensions.
   NANGL = SIZE ( ANGLES, 2 )
   NBOND = SIZE ( BONDS,  2 )

   ! . Check the array dimensions.
   IF ( ( SIZE ( BONDS, 1 ) /= 2 ) .OR. ( SIZE ( ANGLES, 1 ) /= 3 ) ) THEN
      LISTI = - HUGE ( 0 ) ; ALLOCATE ( LISTJ(1:0) ) ; RETURN
   END IF

   ! . Calculate the total length of the temporary array.
   NDIM = NANGL + NBOND

   ! . Allocate the temporary array.
   ALLOCATE ( TMPLIS(1:2,1:NDIM) )

   ! . Copy BONDS to TMPLIS.
   TMPLIS(1:2,1:NBOND) = BONDS(1:2,1:NBOND)

   ! . Copy the first and third columns of ANGLES to TMPLIS.
   TMPLIS(1,NBOND+1:NBOND+NANGL) = ANGLES(1,1:NANGL)
   TMPLIS(2,NBOND+1:NBOND+NANGL) = ANGLES(3,1:NANGL)

   ! . Reorder the temporary dihedral array.
   CALL ORDER_BONDS ( TMPLIS )

   ! . Remove any redundancies in the list.
   CALL REMOVE_BONDS_REDUNDANCIES ( TMPLIS )

   ! . Get the new dimension.
   NDIM = SIZE ( TMPLIS, 2 )

   ! . Allocate space for I123 and J123.
   ALLOCATE ( I123(1:MM_NATOMS+1), J123(1:NDIM) )

   ! . Convert the reordered array to 1-2 form.
   CALL CONVERT_BONDS_TO_12 ( TMPLIS, I123, J123 )

   ! . Deallocate the temporary array.
   DEALLOCATE ( TMPLIS )

   ! . Remove the 1-2 and 1-3 interactions from the 1-2-3-4 lists.
   CALL REMOVE_123_REDUNDANCIES ( I123, J123, I1234, J1234, LISTI, LISTJ )

   ! . Deallocate the 1-2-3 lists.
   DEALLOCATE ( I123, J123 )

   END SUBROUTINE CONVERT_ALL_TO_14

   !-----------------------------------------------------
   SUBROUTINE CONVERT_BONDS_TO_12 ( BONDS, LISTI, LISTJ )
   !-----------------------------------------------------

   ! . Array argument declarations.
   INTEGER, DIMENSION(:,:),        INTENT(IN)  :: BONDS
   INTEGER, DIMENSION(1:MM_NATOMS+1), INTENT(OUT) :: LISTI
   INTEGER, DIMENSION(:),          INTENT(OUT) :: LISTJ

   ! . Local scalars.
   INTEGER :: I, IBOND, N, NBONDS, NOLD

   ! . Get the dimension of BONDS.
   NBONDS = SIZE ( BONDS, 2 )

   ! . Check the array dimensions.
   IF ( ( SIZE ( BONDS, 1 ) /= 2 ) .OR. ( SIZE ( LISTJ ) /= NBONDS ) ) THEN
      LISTI = - HUGE ( 0 ) ; RETURN
   END IF

   ! . Initialize LISTI.
   LISTI = 0

   ! . Loop over the connections to obtain the atom frequencies.
   DO IBOND = 1,NBONDS

      ! . Get the atom index.
      I = BONDS(1,IBOND)

      ! . Update the frequency.
      LISTI(I) = LISTI(I) + 1

   END DO   

   ! . Calculate the cumulative frequencies.
   N = 0
   DO I = 1,MM_NATOMS

      ! . Save N.
      NOLD = N

      ! . Increment the frequency.
      N = N + LISTI(I)

      ! . Fill the first element of LISTI.
      LISTI(I) = NOLD

   END DO

   ! . Set the last element of LISTI.
   LISTI(MM_NATOMS+1) = N

   ! . Copy the second column of BONDS to LISTJ.
   LISTJ(1:NBONDS) = BONDS(2,1:NBONDS)

   END SUBROUTINE CONVERT_BONDS_TO_12

   !-----------------------------------------------------
   SUBROUTINE CONVERT_12_TO_BONDS ( BONDS, LISTI, LISTJ )
   !-----------------------------------------------------

   ! . Array argument declarations.
   INTEGER, DIMENSION(:,:),        INTENT(OUT) :: BONDS
   INTEGER, DIMENSION(1:MM_NATOMS+1), INTENT(IN)  :: LISTI
   INTEGER, DIMENSION(:),          INTENT(IN)  :: LISTJ

   ! . Local scalars.
   INTEGER :: I, J, NBONDS

   ! . Get the dimension of LIST.
   NBONDS = SIZE ( BONDS, 2 )

   ! . Check the array dimensions and the last entry in LISTI.
   IF ( ( SIZE ( BONDS, 1 ) /=      2 ) .OR. &
        ( SIZE ( LISTJ )    /= NBONDS ) .OR. &
        ( LISTI(MM_NATOMS+1)   /= NBONDS ) ) THEN
      BONDS = - HUGE ( 0 ) ; RETURN
   END IF

   ! . Loop over the atoms.
   DO I = 1,MM_NATOMS

      ! . Loop over the interactions.
      DO J = (LISTI(I)+1),LISTI(I+1)

         ! . Fill the elements of BONDS.
         BONDS(1,J) = I
         BONDS(2,J) = LISTJ(J)

      END DO

   END DO

   END SUBROUTINE CONVERT_12_TO_BONDS

   !----------------------------------------------
   SUBROUTINE GENERATE_12 ( LISTI, LISTJ, BUFFER )
   !----------------------------------------------

   ! . Scalar argument declarations.
   REAL ( KIND = DP ), INTENT(IN) :: BUFFER

   ! . Array argument declarations.
   INTEGER, DIMENSION(1:MM_NATOMS+1), INTENT(OUT) :: LISTI
   INTEGER, DIMENSION(:),          POINTER     :: LISTJ

   ! . Local scalars.
   INTEGER            :: I, IINT, IRES, J, JRES, MAXCON, NCONN, NINT, START
   REAL ( KIND = DP ) :: CUTSQ, RADFAC, RIJ2

   ! . Local allocatable arrays.
   INTEGER,            ALLOCATABLE, DIMENSION(:)   :: INDEX, TEMPJ
   REAL ( KIND = DP ), ALLOCATABLE, DIMENSION(:,:) :: CENTER, EXTENT

   ! . Other local arrays.
   REAL ( KIND = DP ), DIMENSION(1:3) :: DR, RMAX, RMIN

   ! . Do nothing if there are no residues.
   IF ( NRESID <= 0 ) RETURN

   !---------------------------------------------------------------------------
   ! . Calculate the centers of the residues and their extents.
   !---------------------------------------------------------------------------
   ! . Guess the maximum number of connections.
   MAXCON = 4 * MM_NATOMS

   ! . Determine the cutoff distance.
   CUTSQ = - HUGE ( 0.0_DP )
   DO I = 1,MM_NATOMS
      CUTSQ = MAX ( CUTSQ, RADII(ATMNUM(I)) )
   END DO
   CUTSQ = ( BUFFER + 2.0_DP * CUTSQ )**2

   ! . Allocate the temporary arrays.
   ALLOCATE ( CENTER(1:3,1:NRESID), EXTENT(1:3,1:NRESID), INDEX(1:MM_NATOMS), TEMPJ(1:MAXCON) )

   !---------------------------------------------------------------------------
   ! . Calculate the centers of the residues and their extents.
   !---------------------------------------------------------------------------
   ! . Loop over the residues.
   DO IRES = 1,NRESID

      ! . Initialize the maximum and minimum arrays.
      RMAX = - HUGE ( 0.0_DP )
      RMIN =   HUGE ( 0.0_DP )

      ! . Find the maximum and minimum coordinates.
      DO I = (RESIND(IRES)+1),RESIND(IRES+1)
         RMAX = MAX ( ATMCRD(1:3,I), RMAX )
         RMIN = MIN ( ATMCRD(1:3,I), RMIN )
      END DO

      ! . Calculate the center and the extent.
      CENTER(1:3,IRES) = 0.5_DP * ( RMAX + RMIN )
      EXTENT(1:3,IRES) = 0.5_DP * ( RMAX - RMIN )

   END DO

   !---------------------------------------------------------------------------
   ! . Create the lists.
   !---------------------------------------------------------------------------
   ! . Initialize the number of connections.
   NCONN = 0

   ! . Outer loop over residues.
   DO IRES = 1,NRESID

      ! . Initialize the number of interactions.
      NINT = 0

      ! . Inner loop over residues to determine which are in range.
      DO JRES = IRES,NRESID

         ! . Calculate the distance differences in each dimension.
         DR = MAX ( ( ABS ( CENTER(1:3,IRES) - CENTER(1:3,JRES) ) - &
                            EXTENT(1:3,IRES) - EXTENT(1:3,JRES) ), 0.0_DP )

         ! . Calculate the distance squared.
         RIJ2 = DOT_PRODUCT ( DR, DR )

         ! . Check to see if the residue is in range.
         IF ( RIJ2 <= CUTSQ ) THEN

            ! . Loop over the atoms in JRES and fill the interaction array.
            DO J = (RESIND(JRES)+1),RESIND(JRES+1)
               NINT = NINT+1
               INDEX(NINT) = J
            END DO

         END IF
      END DO

      ! . Calculate the starting point of the loop.
      START = RESIND(IRES)+1

      ! . Loop over the atoms in IRES.
      DO I = START,RESIND(IRES+1)

         ! . Fill the LISTI element for the atom.
         LISTI(I) = NCONN

         ! . Calculate a local radius factor.
         RADFAC = RADII(ATMNUM(I)) + BUFFER

         ! . Loop over the atoms in the interaction array.
         DO IINT = I-START+2,NINT

            ! . Get the index of the second atom.
            J = INDEX(IINT)

            ! . Calculate the distance between the particles.
            DR   = ATMCRD(1:3,I) - ATMCRD(1:3,J)
            RIJ2 = DOT_PRODUCT ( DR, DR )

            ! . Check the distance.
            IF ( RIJ2 < ( RADFAC + RADII(ATMNUM(J)) )**2 ) THEN

               ! . Check the number of interactions.
               IF ( NCONN == MAXCON ) THEN

                  ! . Save the existing lists.
                  ALLOCATE ( LISTJ(1:NCONN) )
                  LISTJ(1:NCONN) = TEMPJ(1:NCONN)

                  ! . Increment MAXCON.
                  MAXCON = MAXCON + MM_NATOMS

                  ! . Reallocate the temporary lists and resave the data.
                  DEALLOCATE ( TEMPJ )
                  ALLOCATE ( TEMPJ(1:MAXCON) )
                  TEMPJ(1:NCONN) = LISTJ(1:NCONN)

                  ! . Deallocate LISTJ.
                  DEALLOCATE ( LISTJ )

               END IF

               ! . Increment the number of interactions.
               NCONN = NCONN + 1

               ! . Add the interaction to the list.
               TEMPJ(NCONN) = J

            END IF
         END DO
      END DO
   END DO

   ! . Fill the last LISTI element.
   LISTI(MM_NATOMS+1) = NCONN

   ! . Deallocate some arrays.
   DEALLOCATE ( CENTER, EXTENT, INDEX )

   !---------------------------------------------------------------------------
   ! . Finish up.
   !---------------------------------------------------------------------------
   ! . Allocate the LISTJ array.
   ALLOCATE ( LISTJ(1:NCONN) )

   ! . Save the connectivity data.
   LISTJ(1:NCONN) = TEMPJ(1:NCONN)

   ! . Do some printing.
   WRITE ( OUTPUT, "(/25('-'),A,24('-'))" ) " Bond Connectivity Information "
   WRITE ( OUTPUT, "(A,I14,2X,A,F14.1)" ) "Number of Connections  = ", NCONN, &
                                          "Buffer Distance        = ", BUFFER
   WRITE ( OUTPUT, "(80('-'))" )

   ! . Deallocate the remaining temporary arrays.
   DEALLOCATE ( TEMPJ )

   END SUBROUTINE GENERATE_12

   !-------------------------------
   SUBROUTINE ORDER_BONDS ( BONDS )
   !-------------------------------

   ! . Array argument declarations.
   INTEGER, DIMENSION(:,:), INTENT(INOUT)  :: BONDS

   ! . Local scalars.
   INTEGER :: B1, B1OLD, B2, IBOND, NBONDS, START, STOP

   ! . Local arrays.
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: TMPBND

   ! . Get the dimension of BONDS.
   NBONDS = SIZE ( BONDS, 2 )

   ! . Check the array dimensions.
   IF ( ( SIZE ( BONDS, 1 ) /= 2 ) .OR. ( NBONDS <= 0 ) ) RETURN

   ! . Allocate the sorting arrays.
   ALLOCATE ( TMPBND(1:NBONDS,1:2) )

   ! . Loop over the bonds in the list putting the atom with the lower index first.
   DO IBOND = 1,NBONDS
      B1 = BONDS(1,IBOND)
      B2 = BONDS(2,IBOND)
      TMPBND(IBOND,1) = MIN ( B1, B2 )
      TMPBND(IBOND,2) = MAX ( B1, B2 )
   END DO

   ! . Reorder the two arrays (based upon the first index).
   CALL SORT_INTEGER ( TMPBND(1:NBONDS,1), TMPBND(1:NBONDS,2) )

   ! . Initialization.
   B1OLD = TMPBND(1,1)
   START = 1

   ! . Reorder the bonds for each atom.
   DO IBOND = 2,NBONDS

      ! . Get the index of the first atom.
      B1 = TMPBND(IBOND,1)

      ! . The atom is different.
      IF ( B1 /= B1OLD ) THEN

         ! . Get the stopping point.
         STOP = IBOND - 1

         ! . Do the reordering.
         CALL SORT_INTEGER ( TMPBND(START:STOP,2) )

         ! . Reset the starting point.
         START = IBOND

      END IF

      ! . Reset B1OLD.
      B1OLD = B1

   END DO

   ! . Reorder the last group.
   IF ( START < NBONDS ) THEN
      CALL SORT_INTEGER ( TMPBND(START:NBONDS,2) )
   END IF

   ! . Save the sorted bond array.
   DO IBOND = 1,NBONDS
      BONDS(1,IBOND) = TMPBND(IBOND,1)
      BONDS(2,IBOND) = TMPBND(IBOND,2)
   END DO

   ! . Deallocate the temporary arrays.
   DEALLOCATE ( TMPBND )

   END SUBROUTINE ORDER_BONDS

   !---------------------------------------------
   SUBROUTINE REMOVE_BONDS_REDUNDANCIES ( BONDS )
   !---------------------------------------------

   ! . Array argument declarations.
   INTEGER, DIMENSION(:,:), POINTER :: BONDS

   ! . Local scalars.
   INTEGER :: I, NBND, NBONDS

   ! . Local arrays.
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: TMPBND

   ! . Initialization.
   NBND   = 1
   NBONDS = SIZE ( BONDS, 2 )

   ! . Loop over the bonds.
   DO I = 2,NBONDS

      ! . Copy the bond if there is no match.
      IF ( .NOT. ( ( BONDS(1,I) == BONDS(1,I-1) ) .AND. ( BONDS(2,I) == BONDS(2,I-1) ) ) ) THEN
         NBND = NBND + 1
         BONDS(1,NBND) = BONDS(1,I)
         BONDS(2,NBND) = BONDS(2,I)
      END IF

   END DO

   ! . Resize BONDS if there are fewer bonds.
   IF ( NBND < NBONDS ) THEN

      ! . Save BONDS.
      ALLOCATE ( TMPBND(1:2,1:NBND) ) ; TMPBND = BONDS(1:2,1:NBND)

      ! . Resize and fill BONDS.
      DEALLOCATE ( BONDS ) ; ALLOCATE ( BONDS(1:2,1:NBND) ) ; BONDS = TMPBND

      ! . Deallocate the temporary array.
      DEALLOCATE ( TMPBND )

   END IF

   END SUBROUTINE REMOVE_BONDS_REDUNDANCIES

   !----------------------------------------------------------------------------
   SUBROUTINE REMOVE_123_REDUNDANCIES ( I123, J123, I1234, J1234, LISTI, LISTJ )
   !----------------------------------------------------------------------------

   ! . Array argument declarations.
   INTEGER, DIMENSION(1:MM_NATOMS+1), INTENT(IN)  :: I123, I1234
   INTEGER, DIMENSION(:),          INTENT(IN)  :: J123, J1234
   INTEGER, DIMENSION(1:MM_NATOMS+1), INTENT(OUT) :: LISTI
   INTEGER, DIMENSION(:),              POINTER :: LISTJ

   ! . Local scalars.
   INTEGER :: I, IINT, J, N

   ! . Local arrays.
   INTEGER, ALLOCATABLE, DIMENSION(:) :: TMPLIS

   ! . Get and check the list size.
   N = SIZE ( J1234 )

   ! . Allocate space for the lists.
   ALLOCATE ( TMPLIS(1:N) )

   ! . Initialize the number of interactions.
   N = 0

   ! . Outer loop over atoms.
   DO I = 1,MM_NATOMS

      ! . Set the index for the atom.
      LISTI(I) = N

      ! . Loop over the interactions.
      DO IINT = (I1234(I)+1),I1234(I+1)

         ! . Get the second atom of the interaction.
         J = J1234(IINT)

         ! . Include the interaction if there is no match.
         IF ( .NOT. ( ANY ( J123(I123(I)+1:I123(I+1)) == J ) ) ) THEN
            N = N + 1
            TMPLIS(N) = J
         END IF

      END DO
   END DO

   ! . Set the last index.
   LISTI(MM_NATOMS+1) = N

   ! . Allocate and save the lists.
   ALLOCATE ( LISTJ(1:N) ) ; LISTJ = TMPLIS(1:N)

   ! . Deallocate the temporary lists.
   DEALLOCATE ( TMPLIS )

   END SUBROUTINE REMOVE_123_REDUNDANCIES

END MODULE CONNECTIVITY
