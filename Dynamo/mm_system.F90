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
!                              The MM System Module
!===============================================================================
!
! . Subroutines (public):
!
!   MM_SYSTEM_CONSTRUCT            Construct the MM system.
!
! . Subroutines (private):
!
!   FILL_ATOMS_AND_MM_TERMS        Fill the atoms and MM terms data structures.
!   GENERATE_STRUCTURE             Generate the atom, bond and improper lists.
!   READ_SEQUENCE_FILE             Read in the sequence file.
!
!===============================================================================
MODULE MM_SYSTEM

! . Module declarations.
USE DEFINITIONS, ONLY : DP, MAX_RECORD_LENGTH
USE ELEMENTS,    ONLY : MASS
USE FILES,       ONLY : NEXT_UNIT
USE IO_UNITS
USE PARSING
USE STATUS,      ONLY : ERROR

USE ATOMS
USE CONNECTIVITY
USE MM_FILE_DATA
USE MM_FILE_IO,  ONLY : MM_FILE_READ
USE MM_TERMS
USE SEQUENCE
USE SYMMETRY

IMPLICIT NONE
PRIVATE
PUBLIC :: MM_SYSTEM_CONSTRUCT

! . Scalars for first processing step.
INTEGER :: TOTLNK, TOTVAR

! . Arrays for first processing step.
INTEGER, ALLOCATABLE, DIMENSION(:)   :: LNKIND, RINDEX, VARIND, VARRES
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: LNKNUM

! . Scalars for second processing step.
INTEGER :: TOTATOM, TOTBOND, TOTIMPR

! . Arrays for second processing step.
CHARACTER ( LEN = ATOM_NAME_LENGTH ), ALLOCATABLE, DIMENSION(:)   :: TMPNAM
INTEGER,                              ALLOCATABLE, DIMENSION(:)   :: TMPTYP
INTEGER,                              ALLOCATABLE, DIMENSION(:,:) :: TMPBND, TMPIMP
REAL ( KIND = DP ),                   ALLOCATABLE, DIMENSION(:)   :: TMPCHG

!===============================================================================
CONTAINS
!===============================================================================

   !--------------------------------------------------
   SUBROUTINE MM_SYSTEM_CONSTRUCT ( MMFILE, SEQUENCE )
   !--------------------------------------------------

   ! . Essential scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN) :: MMFILE

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: SEQUENCE

   ! . Local scalars.
   INTEGER :: IOSTAT, TUNIT

   ! . Initialize all relevant data structures.
   CALL ATOMS_INITIALIZE
   CALL MM_TERMS_INITIALIZE
   CALL SEQUENCE_INITIALIZE
   CALL SYMMETRY_INITIALIZE

   ! . Open a temporary scratch file.
   TUNIT = NEXT_UNIT ( )
   OPEN ( TUNIT, FORM = "UNFORMATTED", RECL = MAX_RECORD_LENGTH, STATUS = "SCRATCH", IOSTAT = IOSTAT )

   ! . Check for an error.
   IF ( IOSTAT /= 0 ) CALL ERROR ( "MM_SYSTEM_CONSTRUCT", "I/O Error.", IOSTAT )

   ! . Read the MM file data.
   CALL MM_FILE_READ ( MMFILE )

   ! . Read the sequence file.
   CALL READ_SEQUENCE_FILE ( TUNIT, SEQUENCE )

   ! . Check the residue, variant and link data.
   CALL GENERATE_STRUCTURE ( TUNIT )

   ! . Fill the ATOMS and MM_TERMS data structures.
   CALL FILL_ATOMS_AND_MM_TERMS ( TUNIT )

   ! . Close the temporary file.
   CLOSE ( TUNIT )

   ! . Free the MM file data structure.
   CALL MM_FILE_INITIALIZE

   ! . Write out a summary of the atoms, sequence and MM terms data.
   CALL    ATOMS_SUMMARY
   CALL MM_TERMS_SUMMARY
   CALL SEQUENCE_SUMMARY

   END SUBROUTINE MM_SYSTEM_CONSTRUCT

   !-------------------------------------------
   SUBROUTINE FILL_ATOMS_AND_MM_TERMS ( TUNIT )
   !-------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: TUNIT

   ! . Local scalars.
   INTEGER :: I, MM, TOTANGL, TOTDIHE, TOTEXCL, TOTEX14
   LOGICAL :: QMISSING

   ! . Local pointers.
   INTEGER, DIMENSION(:),   POINTER :: TMPEXC, TMPE14
   INTEGER, DIMENSION(:,:), POINTER :: TMPANG, TMPDIH

   ! . Initialization.
   QMISSING = .FALSE.

   ! . Initialize the pointers.
   NULLIFY ( TMPANG, TMPDIH, TMPEXC, TMPE14 )

   ! . Allocate the atoms data structure.
   CALL ATOMS_ALLOCATE ( TOTATOM )

   ! . Fill the atoms data structure elements.
   DO I = 1,MM_NATOMS

      ! . Save the atom name.
      ATMNAM(I) = TMPNAM(I)

      ! . Get the atom type.
      MM = TMPTYP(I)

      ! . Fill the data in ATOMS.
      ATMNUM(I) = TYPES(MM)%NUMBER
      IF ( ATMNUM(I) > 0 ) THEN
         ATMMAS(I) = MASS(ATMNUM(I))
      ELSE
         ATMMAS(I) = 0.0_DP
      END IF

   END DO

   ! . Put the bond array in the correct order.
   CALL ORDER_BONDS ( TMPBND )

   ! . Calculate the angle lists.
   CALL CONNECTIVITY_ANGLES ( TMPANG, TMPBND )

   ! . Calculate the dihedral lists.
   CALL CONNECTIVITY_DIHEDRALS ( TMPDIH, TMPBND, TMPANG )

   ! . Get the total number of angles and dihedrals.
   TOTANGL = SIZE ( TMPANG, 2 )
   TOTDIHE = SIZE ( TMPDIH, 2 )

   ! . Rewind the scratch file.
   REWIND ( TUNIT )

   ! . Determine the number of bonds.
   CALL DETERMINE_BONDS

   ! . Determine the number of angles.
   CALL DETERMINE_ANGLES

   ! . Determine the number of dihedrals.
   CALL DETERMINE_DIHEDRALS

   ! . Determine the number of impropers.
   CALL DETERMINE_IMPROPERS

   ! . Allocate the MM terms data structure.
   CALL MM_TERMS_ALLOCATE ( TOTBOND, TOTANGL, TOTDIHE, TOTIMPR )

   ! . Fill the MM terms atom data.
   DO I = 1,MM_NATOMS

      ! . Get the atom type.
      MM = TMPTYP(I)

      ! . Fill the data elements.
      ATMCHG(I)   = TMPCHG(I)
      ATMEPS(I)   = 2.0_DP * SQRT ( TYPES(MM)%EPSILON )
      ATMSIG(I)   =          SQRT ( TYPES(MM)%SIGMA   )
      ATMTYP(I)   =                 TYPES(MM)%NAME

   END DO

   ! . Set the non-bond scale factors.
   SCALE_EL14 = REF_SCALE_EL
   SCALE_LJ14 = REF_SCALE_LJ

   ! . Calculate the 1-4 arrays.
   ATMCHG14 = SQRT ( SCALE_EL14 ) * ATMCHG
   ATMEPS14 = SQRT ( SCALE_LJ14 ) * ATMEPS

   ! . Rewind the scratch file.
   REWIND ( TUNIT )

   ! . Fill the bond arrays.
   DO I = 1,NBONDS
      READ ( TUNIT ) BONDS(I)%I, BONDS(I)%J, BONDS(I)%EQ, BONDS(I)%FC
   END DO

   ! . Fill the angle arrays.
   DO I = 1,NANGLES
      READ ( TUNIT ) ANGLES(I)%I, ANGLES(I)%J, ANGLES(I)%K, ANGLES(I)%EQ, ANGLES(I)%FC
   END DO

   ! . Fill the dihedral arrays.
   DO I = 1,NDIHEDRALS
      READ ( TUNIT ) DIHEDRALS(I)%I, DIHEDRALS(I)%J, DIHEDRALS(I)%K, DIHEDRALS(I)%L, &
                     DIHEDRALS(I)%FC, DIHEDRALS(I)%PERIOD, DIHEDRALS(I)%PHASE
   END DO

   ! . Fill the improper arrays.
   DO I = 1,NIMPROPERS
      READ ( TUNIT ) IMPROPERS(I)%I, IMPROPERS(I)%J, IMPROPERS(I)%K, IMPROPERS(I)%L, &
                     IMPROPERS(I)%FC, IMPROPERS(I)%PERIOD, IMPROPERS(I)%PHASE
   END DO

   ! . Calculate the 1-2-3-4 exclusion lists.
   CALL CONVERT_ALL_TO_1234 ( TMPBND, TMPANG, TMPDIH, ATMEXCI, TMPEXC )

   ! . Calculate the 1-4 exclusion lists.
   CALL CONVERT_ALL_TO_14 ( TMPBND, TMPANG, ATMEXCI, TMPEXC, ATME14I, TMPE14 )

   ! . Get the total number of exclusions.
   TOTEXCL = SIZE ( TMPEXC )
   TOTEX14 = SIZE ( TMPE14 )

   ! . Allocate the exclusion arrays.
   ALLOCATE ( ATMEXCJ(1:TOTEXCL), ATME14J(1:TOTEX14) )

   ! . Copy over the exclusion arrays.
   ATMEXCJ = TMPEXC
   ATME14J = TMPE14

   ! . Deallocate the temporary arrays.
   DEALLOCATE ( TMPANG, TMPBND, TMPCHG, TMPDIH, TMPEXC, TMPE14, TMPIMP, TMPNAM, TMPTYP )

   ! . Check for a missing parameter.
   IF ( QMISSING ) CALL ERROR ( "FILL_ATOMS_AND_MM_TERMS", "Missing energy function parameters." )

   !============================================================================
   CONTAINS
   !============================================================================

      !--------------------------
      SUBROUTINE DETERMINE_ANGLES
      !--------------------------

      ! . Local scalars.
      INTEGER :: A1, A2, A3, I, MM, N, NANGL, T1, T2, T3

      ! . Store the number of angles.
      N     = 0
      NANGL = TOTANGL

      ! . Loop over the angles.
      DO I = 1,NANGL

         ! . Get the atom indices.
         A1 = TMPANG(1,I)
         A2 = TMPANG(2,I)         
         A3 = TMPANG(3,I)         

         ! . Get the types of the atoms.
         T1 = TMPTYP(A1)
         T2 = TMPTYP(A2)
         T3 = TMPTYP(A3)

         ! . Loop over the parameters that are available.
         DO MM = 1,NANGLPS
            ! . Check for a match.
            IF ( ( T1 == ANGLPS(MM)%TYPE1 .AND. T2 == ANGLPS(MM)%TYPE2 .AND. T3 == ANGLPS(MM)%TYPE3 ) .OR. &
                 ( T3 == ANGLPS(MM)%TYPE1 .AND. T2 == ANGLPS(MM)%TYPE2 .AND. T1 == ANGLPS(MM)%TYPE3 ) ) THEN
               ! . Save the angle if the force constant is non-zero.
               IF ( ANGLPS(MM)%FC /= 0.0_DP ) THEN
                  N = N + 1
                  WRITE ( TUNIT ) A1, A2, A3, ANGLPS(MM)%EQ, ANGLPS(MM)%FC
               END IF
               GO TO 10
            END IF
         END DO

         ! . No match.
         QMISSING = .TRUE.
         WRITE ( OUTPUT, "(A,3(2X,A))" ) "Missing angle parameter for:", TYPES(T1)%NAME, TYPES(T2)%NAME, TYPES(T3)%NAME

         ! . End of the loop.
         10 CONTINUE

      END DO

      ! . Set the number of angles.
      TOTANGL = N

      END SUBROUTINE DETERMINE_ANGLES

      !-------------------------
      SUBROUTINE DETERMINE_BONDS
      !-------------------------

      ! . Local scalars.
      INTEGER :: B1, B2, I, MM, N, NBND, T1, T2

      ! . Store the number of bonds.
      N    = 0
      NBND = TOTBOND

      ! . Loop over the bonds.
      DO I = 1,NBND

         ! . Get the atom indices.
         B1 = TMPBND(1,I)
         B2 = TMPBND(2,I)         

         ! . Get the types of the atoms.
         T1 = TMPTYP(B1)
         T2 = TMPTYP(B2)

         ! . Loop over the parameters that are available.
         DO MM = 1,NBONDPS
            ! . Check for a match.
            IF ( ( T1 == BONDPS(MM)%TYPE1 .AND. T2 == BONDPS(MM)%TYPE2 ) .OR. &
                 ( T2 == BONDPS(MM)%TYPE1 .AND. T1 == BONDPS(MM)%TYPE2 ) ) THEN

               ! . Save the bond if the force constant is non-zero.
               IF ( BONDPS(MM)%FC /= 0.0_DP ) THEN
                  N = N + 1
                  WRITE ( TUNIT ) B1, B2, BONDPS(MM)%EQ, BONDPS(MM)%FC
               END IF
               GO TO 10
            END IF
         END DO

         ! . No match.
         QMISSING = .TRUE.
         WRITE ( OUTPUT, "(A,2(2X,A))" ) "Missing bond parameter for:", TYPES(T1)%NAME, TYPES(T2)%NAME

         ! . End of the loop.
         10 CONTINUE

      END DO

      ! . Set the number of bonds.
      TOTBOND = N

      END SUBROUTINE DETERMINE_BONDS

      !-----------------------------
      SUBROUTINE DETERMINE_DIHEDRALS
      !-----------------------------

      ! . Local scalars.
      INTEGER :: D1, D2, D3, D4, I, MM, MM1234, MMX23X, N, NDIHE, TMP, T1, T2, T3, T4

      ! . Store the number of dihedrals.
      N     = 0
      NDIHE = TOTDIHE

      ! . Loop over the dihedrals.
      DO I = 1,NDIHE

         ! . Get the atom indices.
         D1 = TMPDIH(1,I)
         D2 = TMPDIH(2,I)         
         D3 = TMPDIH(3,I)         
         D4 = TMPDIH(4,I)         

         ! . Get the types of the atoms.
         T1 = TMPTYP(D1)
         T2 = TMPTYP(D2)
         T3 = TMPTYP(D3)
         T4 = TMPTYP(D4)

         ! . Order the types.
         IF ( T2 > T3 ) THEN
            TMP = T2 ; T2 = T3 ; T3 = TMP
            TMP = T1 ; T1 = T4 ; T4 = TMP
         ELSE IF ( T2 == T3 ) THEN
            IF ( T1 > T4 ) THEN
               TMP = T1 ; T1 = T4 ; T4 = TMP
            END IF
         END IF

         ! . Initialization.
         MM1234 = 0
         MMX23X = 0

         ! . Loop over the parameters that are available.
         DO MM = 1,NDIHEPS

            ! . Check for a match.
            IF ( T2 == DIHEPS(MM)%TYPE2 .AND. T3 == DIHEPS(MM)%TYPE3 ) THEN

               ! . There is a specific match.
               IF ( T1 == DIHEPS(MM)%TYPE1 .AND. T4 == DIHEPS(MM)%TYPE4 ) THEN
                  MM1234 = MM
               ! . There is a non-specific match.
               ELSE IF ( DIHEPS(MM)%TYPE1 == DUMMY_CODE .AND. DIHEPS(MM)%TYPE4 == DUMMY_CODE ) THEN
                  MMX23X = MM
               END IF

            END IF

         END DO

         ! . There was no specific match.
         IF ( MM1234 == 0 ) THEN
            ! . There is a non-specific match.
            IF ( MMX23X > 0 ) THEN
               MM = MMX23X
            ! . There is no match.
            ELSE
               MM = 0
               QMISSING = .TRUE.
               WRITE ( OUTPUT, "(A,4(2X,A))" ) "Missing dihedral parameter for:", &
                      TYPES(T1)%NAME, TYPES(T2)%NAME, TYPES(T3)%NAME, TYPES(T4)%NAME
            END IF
         ! . There was a specific match.
         ELSE
            MM = MM1234
         END IF

         ! . Write out the parameter data.
         IF ( MM > 0 ) THEN

            ! . V0.
            IF ( DIHEPS(MM)%V0 /= 0.0_DP ) THEN
               N = N + 1
               WRITE ( TUNIT ) D1, D2, D3, D4, DIHEPS(MM)%V0, 0,  1.0_DP
            END IF

            ! . V1.
            IF ( DIHEPS(MM)%V1 /= 0.0_DP ) THEN
               N = N + 1
               WRITE ( TUNIT ) D1, D2, D3, D4, DIHEPS(MM)%V1, 1,  1.0_DP
            END IF

            ! . V2.
            IF ( DIHEPS(MM)%V2 /= 0.0_DP ) THEN
               N = N + 1
               WRITE ( TUNIT ) D1, D2, D3, D4, DIHEPS(MM)%V2, 2, -1.0_DP
            END IF

            ! . V3.
            IF ( DIHEPS(MM)%V3 /= 0.0_DP ) THEN
               N = N + 1
               WRITE ( TUNIT ) D1, D2, D3, D4, DIHEPS(MM)%V3, 3,  1.0_DP
            END IF
         END IF

      END DO

      ! . Set the number of dihedrals.
      TOTDIHE = N

      END SUBROUTINE DETERMINE_DIHEDRALS

      !-----------------------------
      SUBROUTINE DETERMINE_IMPROPERS
      !-----------------------------

      ! . Local scalars.
      INTEGER :: I, I1, I2, I3, I4, MM, MM1234, MMX234, MMXX34, N, NIMPR, T1, T2, T3, T4

      ! . Store the number of impropers.
      N     = 0
      NIMPR = TOTIMPR

      ! . Loop over the impropers.
      DO I = 1,NIMPR

         ! . Get the atom indices.
         I1 = TMPIMP(1,I)
         I2 = TMPIMP(2,I)         
         I3 = TMPIMP(3,I)         
         I4 = TMPIMP(4,I)         

         ! . Get the types of the atoms.
         T1 = TMPTYP(I1)
         T2 = TMPTYP(I2)
         T3 = TMPTYP(I3)
         T4 = TMPTYP(I4)

         ! . Initialization.
         MM1234 = 0
         MMX234 = 0
         MMXX34 = 0

         ! . Loop over the parameters that are available.
         DO MM = 1,NIMPRPS
            ! . Check for a match.
            IF ( T3 == IMPRPS(MM)%TYPE3 .AND. T4 == IMPRPS(MM)%TYPE4 ) THEN
               IF ( T2 == IMPRPS(MM)%TYPE2 ) THEN
                  IF ( T1 == IMPRPS(MM)%TYPE1 ) THEN
                    MM1234 = MM
                  ELSE IF ( IMPRPS(MM)%TYPE1 == DUMMY_CODE ) THEN
                    MMX234 = MM
                  END IF
               ELSE IF ( IMPRPS(MM)%TYPE2 == DUMMY_CODE ) THEN
                  IF ( IMPRPS(MM)%TYPE1 == DUMMY_CODE ) MMXX34 = MM
               END IF
            END IF
         END DO

         ! . There is a specific match.
         IF ( MM1234 > 0 ) THEN
            MM = MM1234
         ELSE IF ( MMX234 > 0 ) THEN
            MM = MMX234
         ELSE IF ( MMXX34 > 0 ) THEN
            MM = MMXX34
         ! . There is no match.
         ELSE
            MM = 0
            QMISSING = .TRUE.
            WRITE ( OUTPUT, "(A,4(2X,A))" ) "Missing improper parameter for:", &
                   TYPES(T1)%NAME, TYPES(T2)%NAME, TYPES(T3)%NAME, TYPES(T4)%NAME
         END IF

         ! . Write out the parameter data.
         IF ( MM > 0 ) THEN

            ! . V0.
            IF ( IMPRPS(MM)%V0 /= 0.0_DP ) THEN
               N = N + 1
               WRITE ( TUNIT ) I1, I2, I3, I4, IMPRPS(MM)%V0, 0,  1.0_DP
            END IF

            ! . V1.
            IF ( IMPRPS(MM)%V1 /= 0.0_DP ) THEN
               N = N + 1
               WRITE ( TUNIT ) I1, I2, I3, I4, IMPRPS(MM)%V1, 1,  1.0_DP
            END IF

            ! . V2.
            IF ( IMPRPS(MM)%V2 /= 0.0_DP ) THEN
               N = N + 1
               WRITE ( TUNIT ) I1, I2, I3, I4, IMPRPS(MM)%V2, 2, -1.0_DP
            END IF

            ! . V3.
            IF ( IMPRPS(MM)%V3 /= 0.0_DP ) THEN
               N = N + 1
               WRITE ( TUNIT ) I1, I2, I3, I4, IMPRPS(MM)%V3, 3,  1.0_DP
            END IF
         END IF

      END DO

      ! . Set the number of impropers.
      TOTIMPR = N

      END SUBROUTINE DETERMINE_IMPROPERS

   END SUBROUTINE FILL_ATOMS_AND_MM_TERMS

   !--------------------------------------
   SUBROUTINE GENERATE_STRUCTURE ( TUNIT )
   !--------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: TUNIT

   ! . Local parameters.
   INTEGER, PARAMETER :: CODE_LATERAL = -1, CODE_MINUS = -2, CODE_PLUS = -3

   ! . Local type declarations.
   TYPE RESIDUE_INFORMATION
      INTEGER :: MM_NATOMS, NBONDS, NIMPRS
      CHARACTER ( LEN = ATOM_NAME_LENGTH ), DIMENSION(:),   POINTER :: NAMES
      CHARACTER ( LEN = ATOM_NAME_LENGTH ), DIMENSION(:,:), POINTER :: BONDS, IMPRS
      INTEGER,                              DIMENSION(:),   POINTER :: BINDX, IINDX, TYPES
      REAL ( KIND = DP ),                   DIMENSION(:),   POINTER :: CHRGS
   END TYPE RESIDUE_INFORMATION

   ! . Local scalars.
   INTEGER                   :: I, INDATM, INDBND, INDIMP, J, K, MM, NATM, NBND, NIMP, POS_NEW, POS_OLD
   TYPE(RESIDUE_INFORMATION) :: RESINF

   ! . Local arrays.
   INTEGER, ALLOCATABLE, DIMENSION(:)   :: LNKHGH, LNKLOW, NEG_LINK
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: BNDATM, IMPATM

   !-------------------------------------------------------------------------
   ! . Process the information in the temporary arrays.
   !-------------------------------------------------------------------------
   ! . Rewind the scratch file.
   REWIND ( TUNIT )

   ! . Allocate some temporary arrays.
   ALLOCATE ( LNKHGH(1:TOTLNK), LNKLOW(1:TOTLNK), NEG_LINK(1:NRESID+1) )

   ! . Initialize some counters.
   TOTATOM = 0
   TOTBOND = 0
   TOTIMPR = 0

   ! . Initialize the atom indices for the between residue links.
   POS_NEW  = -1
   NEG_LINK = -1

   ! . Loop over the residues.
   DO I = 1,NRESID

      ! . Fill the RESIND element for the residue.
      RESIND(I) = TOTATOM

      ! . Get the MM index for the residue.
      MM = RINDEX(I)

      ! . Get the number of atoms, bonds and impropers for the residue.
      NATM = RESIDUES(MM)%MM_NATOMS
      NBND = RESIDUES(MM)%NBONDS
      NIMP = RESIDUES(MM)%NIMPROPERS

      ! . Allocate and fill the RESINF type for the residue.
      CALL RESIDUE_ALLOCATE ( NATM, NBND, NIMP, RESINF )
      CALL RESIDUE_FILL ( RESINF, RESIDUES(MM) )

      ! . Check for a variant.
      IF ( ANY ( VARRES(1:TOTVAR) == I ) ) THEN

         ! . Loop over the variants.
         DO J = 1,TOTVAR
            IF ( VARRES(J) == I ) THEN
               MM = VARIND(J)
               CALL RESIDUE_MODIFY ( RESINF, VARIANTS(MM), 0 )
            END IF
         END DO            

      END IF

      ! . Check for a link.
      IF ( ANY ( LNKNUM(1:2,1:TOTLNK) == I ) ) THEN

         ! . Loop over the links.
         DO J = 1,TOTLNK
            IF ( LNKNUM(1,J) == I ) THEN
               MM = LNKIND(J)
               CALL RESIDUE_MODIFY ( RESINF, RESLINK(MM)%RESIDUE1, J )
            ELSE IF ( LNKNUM(2,J) == I ) THEN
               MM = LNKIND(J)
               CALL RESIDUE_MODIFY ( RESINF, RESLINK(MM)%RESIDUE2, J )
            END IF
         END DO            

      END IF

      ! . Assign some counters.
      NATM = RESINF%MM_NATOMS
      NBND = RESINF%NBONDS
      NIMP = RESINF%NIMPRS

      ! . Allocate space for the bond and improper indices.
      ALLOCATE ( BNDATM(1:2,1:NBND), IMPATM(1:4,1:NIMP) )

      ! . Transform the name arrays to the index arrays.
      CALL NAMES_TO_INDICES ( TOTATOM, 2, NBND, RESINF%BONDS, RESINF%NAMES, BNDATM )
      CALL NAMES_TO_INDICES ( TOTATOM, 4, NIMP, RESINF%IMPRS, RESINF%NAMES, IMPATM )

      ! . Save POS_OLD.
      POS_OLD = POS_NEW

      ! . Produce the final bond list.
      CALL RESIDUE_BONDS ( I, POS_NEW, NEG_LINK(I), NBND, RESINF%BINDX, BNDATM, LNKNUM, LNKHGH, LNKLOW )

      ! . Produce an intermediate improper list.
      CALL RESIDUE_IMPROPERS ( I, POS_OLD, NIMP, RESINF%IINDX, IMPATM, LNKNUM, LNKLOW )

      ! . Write out the residue information to the scratch file.
      WRITE ( TUNIT ) NATM, NBND, NIMP

      ! . Write out the atom arrays.
      IF ( NATM > 0 ) THEN
         WRITE ( TUNIT ) RESINF%CHRGS(1:NATM)
         WRITE ( TUNIT ) RESINF%NAMES(1:NATM)
         WRITE ( TUNIT ) RESINF%TYPES(1:NATM)
      END IF

      ! . Write out the bond array.
      IF ( NBND > 0 ) WRITE ( TUNIT ) BNDATM(1:2,1:NBND)

      ! . Fill the improper array.
      IF ( NIMP > 0 ) WRITE ( TUNIT ) IMPATM(1:4,1:NIMP)

      ! . Reinitialize RESINF.
      CALL RESIDUE_INITIALIZE ( RESINF )

      ! . Deallocate the bond and improper arrays.
      DEALLOCATE ( BNDATM, IMPATM )

      ! . Increment the atom, bond and improper counters.
      TOTATOM = TOTATOM + NATM
      TOTBOND = TOTBOND + NBND
      TOTIMPR = TOTIMPR + NIMP

   END DO

   ! . Fill the last RESIND element.
   RESIND(NRESID+1) = TOTATOM

   ! . Deallocate some temporary arrays.
   DEALLOCATE ( LNKIND, LNKLOW, LNKNUM, RINDEX, VARIND, VARRES )

   ! . Check the value of POS_OLD.
   IF ( POS_NEW > 0 ) THEN
      CALL ERROR ( "GENERATE_STRUCTURE", "The last residue cannot have a `+R' bond specification." )
   END IF

   !-------------------------------------------------------------------------
   ! . Fill the temporary atom, bond and improper arrays.
   !-------------------------------------------------------------------------
   ! . Allocate the temporary arrays.
   ALLOCATE ( TMPBND(1:2,1:TOTBOND), TMPCHG(1:TOTATOM), TMPIMP(1:4,1:TOTIMPR), &
              TMPNAM(1:TOTATOM), TMPTYP(1:TOTATOM) )

   ! . Rewind the scratch file.
   REWIND ( TUNIT )

   ! . Initialization.
   INDATM = 0
   INDBND = 0
   INDIMP = 0

   ! . Loop over the residues.
   DO I = 1,NRESID

      ! . Get the counters for the residue.
      READ ( TUNIT ) NATM, NBND, NIMP

      ! . Fill the atom arrays.
      IF ( NATM > 0 ) THEN
         READ ( TUNIT ) TMPCHG(INDATM+1:INDATM+NATM)
         READ ( TUNIT ) TMPNAM(INDATM+1:INDATM+NATM)
         READ ( TUNIT ) TMPTYP(INDATM+1:INDATM+NATM)
      END IF

      ! . Fill the bond array.
      IF ( NBND > 0 ) READ ( TUNIT ) ( ( TMPBND(J,K), J = 1,2 ), K = INDBND+1,INDBND+NBND )

      ! . Fill the improper array.
      IF ( NIMP > 0 ) READ ( TUNIT ) ( ( TMPIMP(J,K), J = 1,4 ), K = INDIMP+1,INDIMP+NIMP )

      ! . Finish processing of the improper lists.
      CALL RESIDUE_IMPROPERS_LAST ( I, NEG_LINK, INDIMP+1, INDIMP+NIMP, TMPIMP, LNKHGH )

      ! . Increment the index counters.
      INDATM = INDATM + NATM
      INDBND = INDBND + NBND
      INDIMP = INDIMP + NIMP

   END DO


   ! . Deallocate the remaining temporary arrays.
   DEALLOCATE ( LNKHGH, NEG_LINK )

   !============================================================================
   CONTAINS
   !============================================================================

      !--------------------------------------------------------------------------
      SUBROUTINE NAMES_TO_INDICES ( TOTATOM, M, N, NAME_IN, NAME_REF, INDEX_OUT )
      !--------------------------------------------------------------------------

      ! . Scalar arguments.
      INTEGER, INTENT(IN) :: M, N, TOTATOM

      ! . Array arguments.
!      CHARACTER ( LEN = ATOM_NAME_LENGTH ), DIMENSION(:),   INTENT(IN)  :: NAME_REF
!      CHARACTER ( LEN = ATOM_NAME_LENGTH ), DIMENSION(:,:), INTENT(IN)  :: NAME_IN
      CHARACTER ( LEN = * ), DIMENSION(:),   INTENT(IN)  :: NAME_REF
      CHARACTER ( LEN = * ), DIMENSION(:,:), INTENT(IN)  :: NAME_IN
      INTEGER,               DIMENSION(:,:), INTENT(OUT) :: INDEX_OUT

      ! . Local scalars.
      INTEGER :: I, J, K, NREF, TMP

      ! . Initialization.
      NREF = SIZE ( NAME_REF )

      ! . Loop over the elements in the array.
      DO I = 1,N
         DO J = 1,M

            ! . Find the atom index.
            TMP = 0
            SELECT CASE ( NAME_IN(J,I) )
            CASE ( LINK_LATERAL ) ; TMP = CODE_LATERAL
            CASE ( LINK_MINUS   ) ; TMP = CODE_MINUS
            CASE ( LINK_PLUS    ) ; TMP = CODE_PLUS
            CASE DEFAULT
               DO K = 1,NREF
                  IF ( NAME_IN(J,I) == NAME_REF(K) ) THEN
                     TMP = TOTATOM + K ; EXIT
                  END IF
               END DO
            END SELECT

            ! . Assign the atom index.
            IF ( TMP == 0 ) THEN
               CALL ERROR ( "GENERATE_STRUCTURE", "Unknown atom name: "//NAME_IN(J,I)//"." )
            ELSE
               INDEX_OUT(J,I) = TMP
            END IF

         END DO
      END DO

      END SUBROUTINE NAMES_TO_INDICES

      !-------------------------------------------------------
      SUBROUTINE RESIDUE_ALLOCATE ( NATM, NBND, NIMP, RESINF )
      !-------------------------------------------------------

      ! . Scalar arguments.
      INTEGER,                   INTENT(IN)  :: NATM, NBND, NIMP
      TYPE(RESIDUE_INFORMATION), INTENT(OUT) :: RESINF

      ! . Initialize the RESINF scalars.
      RESINF%MM_NATOMS = NATM
      RESINF%NBONDS = NBND
      RESINF%NIMPRS = NIMP

      ! . Initialize the pointers.
      NULLIFY ( RESINF%BINDX, RESINF%BONDS, RESINF%CHRGS, RESINF%IINDX, RESINF%IMPRS, RESINF%NAMES, RESINF%TYPES )

      ! . Allocate the arrays in RESINF.
      ALLOCATE ( RESINF%CHRGS(1:NATM), RESINF%NAMES(1:NATM), RESINF%TYPES(1:NATM) )
      ALLOCATE ( RESINF%BONDS(1:2,1:NBND), RESINF%BINDX(1:NBND) )
      ALLOCATE ( RESINF%IMPRS(1:4,1:NIMP), RESINF%IINDX(1:NIMP) )

      END SUBROUTINE RESIDUE_ALLOCATE

      !---------------------------------------------------------------------------------------------------
      SUBROUTINE RESIDUE_BONDS ( CURRES, POS_LINK, NEG_LINK, NBND, BINDX, BNDATM, LNKNUM, LNKHGH, LNKLOW )
      !---------------------------------------------------------------------------------------------------

      ! . Scalar arguments.
      INTEGER, INTENT(IN)    :: CURRES
      INTEGER, INTENT(INOUT) :: NBND, NEG_LINK, POS_LINK

      ! . Array arguments.
      INTEGER, DIMENSION(:),   INTENT(IN)    :: BINDX
      INTEGER, DIMENSION(:,:), INTENT(IN)    :: LNKNUM
      INTEGER, DIMENSION(:),   INTENT(INOUT) :: LNKHGH, LNKLOW
      INTEGER, DIMENSION(:,:), INTENT(INOUT) :: BNDATM

      ! . Local scalars.
      INTEGER :: ATM_NEG, ATM_POS, B1, B2, I, NBNDTMP

      ! . Initialization.
      ATM_NEG = -1
      ATM_POS = -1

      ! . Store the number of bonds counter.
      NBNDTMP = NBND
      NBND    = 0

      ! . Loop over the bonds.
      DO I = 1,NBNDTMP

         ! . Get the bond indices.
         B1 = MIN ( BNDATM(1,I), BNDATM(2,I) )
         B2 = MAX ( BNDATM(1,I), BNDATM(2,I) )

         ! . The bond is normal.
         IF ( B1 > 0 ) THEN

            NBND = NBND + 1
            BNDATM(1,NBND) = B1
            BNDATM(2,NBND) = B2

         ! . The bond is to another residue.
         ELSE

            ! . There is a second link atom.
            IF ( B2 <= 0 ) CALL ERROR ( "GENERATE_STRUCTURE", "A bond contains no non-link atoms." )

            ! . Check the type of connection.
            SELECT CASE ( B1 )
            ! . Bond to a distant residue.
            CASE ( CODE_LATERAL )

               ! . The current residue has the highest index in the link.
               IF ( CURRES == MAXVAL ( LNKNUM(1:2,BINDX(I)) ) ) THEN

                  ! . Store the bond.
                  IF ( LNKLOW(BINDX(I)) > 0 ) THEN

                     ! . Store the bond.
                     NBND = NBND + 1
                     BNDATM(1,NBND) = LNKLOW(BINDX(I))
                     BNDATM(2,NBND) = B2

                     ! . Store the high index.
                     LNKHGH(BINDX(I)) = B1

                  ! . There is no atom specified for the previous residue of the link.
                  ELSE
                     CALL ERROR ( "GENERATE_STRUCTURE", "There is a missing `*R' specification in a link." )
                  END IF

               ! . Store the link.
               ELSE
                  LNKLOW(BINDX(I)) = B2
               END IF

            ! . Bond to the previous residue.
            CASE ( CODE_MINUS )

               ! . Check ATM_NEG.
               IF ( ATM_NEG > 0 ) THEN
                  CALL ERROR ( "GENERATE_STRUCTURE", "Only one connection to the previous residue is allowed." )
               ELSE
                  ATM_NEG = B2
               END IF

               ! . Save the bond.
               IF ( POS_LINK > 0 ) THEN
                  NBND = NBND + 1
                  BNDATM(1,NBND) = POS_LINK
                  BNDATM(2,NBND) = B2
               ! . The previous residue did not have a `+R'.
               ELSE
                  CALL ERROR ( "GENERATE_STRUCTURE", "Invalid bond to a previous residue." )
               END IF

            ! . Bond to the next residue.
            CASE ( CODE_PLUS )

               ! . Check ATM_POS.
               IF ( ATM_POS > 0 ) THEN
                  CALL ERROR ( "GENERATE_STRUCTURE", "Only one connection to the next residue is allowed." )
               ELSE
                  ATM_POS = B2
               END IF

            ! . Unknown index.
            CASE DEFAULT
               CALL ERROR ( "GENERATE_STRUCTURE", "Invalid bond atom index." )
            END SELECT
         END IF

      END DO

      ! . Set the NEG_LINK and POS_LINK counters.
      NEG_LINK = ATM_NEG
      POS_LINK = ATM_POS

      END SUBROUTINE RESIDUE_BONDS

      !------------------------------------------------
      SUBROUTINE RESIDUE_COPY ( RESINF_IN, RESINF_OUT )
      !------------------------------------------------

      ! . Scalar arguments.
      TYPE(RESIDUE_INFORMATION), INTENT(IN)    :: RESINF_IN
      TYPE(RESIDUE_INFORMATION), INTENT(INOUT) :: RESINF_OUT

      ! . Local scalars.
      INTEGER :: NATM, NBND, NIMP

      ! . Determine the number of elements to copy.
      NATM = MIN ( RESINF_IN%MM_NATOMS, RESINF_OUT%MM_NATOMS )
      NBND = MIN ( RESINF_IN%NBONDS, RESINF_OUT%NBONDS )
      NIMP = MIN ( RESINF_IN%NIMPRS, RESINF_OUT%NIMPRS )

      ! . Copy the arrays.
      IF ( NATM > 0 ) THEN
         RESINF_OUT%CHRGS(1:NATM) = RESINF_IN%CHRGS(1:NATM)
         RESINF_OUT%NAMES(1:NATM) = RESINF_IN%NAMES(1:NATM)
         RESINF_OUT%TYPES(1:NATM) = RESINF_IN%TYPES(1:NATM)
      END IF

      IF ( NBND > 0 ) THEN
         RESINF_OUT%BINDX(    1:NBND) = RESINF_IN%BINDX(    1:NBND)
         RESINF_OUT%BONDS(1:2,1:NBND) = RESINF_IN%BONDS(1:2,1:NBND)
      END IF

      IF ( NIMP > 0 ) THEN
         RESINF_OUT%IINDX(    1:NIMP) = RESINF_IN%IINDX(    1:NIMP)
         RESINF_OUT%IMPRS(1:4,1:NIMP) = RESINF_IN%IMPRS(1:4,1:NIMP)
      END IF

      END SUBROUTINE RESIDUE_COPY

      !------------------------------------------
      SUBROUTINE RESIDUE_FILL ( RESINF, RESIDUE )
      !------------------------------------------

      ! . Scalar arguments.
      TYPE(RESIDUE_INFORMATION), INTENT(INOUT) :: RESINF
      TYPE(RESIDUE_TYPE),        INTENT(IN)    :: RESIDUE

      ! . Copy the arrays.
      IF ( RESINF%MM_NATOMS > 0 ) THEN
         RESINF%CHRGS(1:RESINF%MM_NATOMS) = RESIDUE%CHARGES(1:RESINF%MM_NATOMS)
         RESINF%NAMES(1:RESINF%MM_NATOMS) = RESIDUE%NAMES(1:RESINF%MM_NATOMS)
         RESINF%TYPES(1:RESINF%MM_NATOMS) = RESIDUE%TYPES(1:RESINF%MM_NATOMS)
      END IF

      IF ( RESINF%NBONDS > 0 ) THEN
         RESINF%BONDS(1:2,1:RESINF%NBONDS) = RESIDUE%BONDS(1:2,1:RESINF%NBONDS)
         RESINF%BINDX(1:RESINF%NBONDS) = 0
      END IF

      IF ( RESINF%NIMPRS > 0 ) THEN
         RESINF%IMPRS(1:4,1:RESINF%NIMPRS) = RESIDUE%IMPROPERS(1:4,1:RESINF%NIMPRS)
         RESINF%IINDX(1:RESINF%NIMPRS) = 0
      END IF

      END SUBROUTINE RESIDUE_FILL

      !-------------------------------------------------------------------------------------
      SUBROUTINE RESIDUE_IMPROPERS ( CURRES, POS_OLD, NIMP, IINDX, IMPATM, LNKNUM, LNKLOW )
      !-------------------------------------------------------------------------------------

      ! . Scalar arguments.
      INTEGER, INTENT(IN) :: CURRES,NIMP, POS_OLD

      ! . Array arguments.
      INTEGER, DIMENSION(:),   INTENT(IN)    :: IINDX, LNKLOW
      INTEGER, DIMENSION(:,:), INTENT(IN)    :: LNKNUM
      INTEGER, DIMENSION(:,:), INTENT(INOUT) :: IMPATM

      ! . Local scalars.
      INTEGER :: I, J

      ! . Loop over the impropers.
      DO I = 1,NIMP

         ! . Loop over the atoms in the impropers.
         DO J = 1,4

            ! . Branch on the type of atom.
            SELECT CASE ( IMPATM(J,I) )
            ! . Bond to a distant residue.
            CASE ( CODE_LATERAL )

               ! . The current residue has the highest index in the link.
               IF ( CURRES == MAXVAL ( LNKNUM(1:2,IINDX(I)) ) ) THEN

                  ! . Store the atom index.
                  IF ( LNKLOW(IINDX(I)) > 0 ) THEN

                     IMPATM(J,I) = LNKLOW(IINDX(I))

                  ! . There is no atom specified for the previous residue of the link.
                  ELSE
                     CALL ERROR ( "GENERATE_STRUCTURE", "There is a missing `*R' specification in a link." )
                  END IF

               ! . Reset the atom index.
               ELSE
                  IMPATM(J,I) = - IINDX(I)
               END IF

            ! . Bond to the previous residue.
            CASE ( CODE_MINUS )

               ! . Store the atom index.
               IF ( POS_OLD > 0 ) THEN
                  IMPATM(J,I) = POS_OLD
               ELSE
                  CALL ERROR ( "GENERATE_STRUCTURE", "Improper specified to previous residue without bond to `-R'." )
               END IF

            ! . Bond to the next residue.
            CASE ( CODE_PLUS ) ; IMPATM(J,I) = 0
            END SELECT

         END DO
      END DO

      END SUBROUTINE RESIDUE_IMPROPERS

      !----------------------------------------------------------------------------------
      SUBROUTINE RESIDUE_IMPROPERS_LAST ( CURRES, NEG_LINK, START, STOP, TMPIMP, LNKHGH )
      !----------------------------------------------------------------------------------

      ! . Scalar arguments.
      INTEGER, INTENT(IN) :: CURRES, START, STOP

      ! . Array arguments.
      INTEGER, DIMENSION(:),   INTENT(IN)    :: LNKHGH, NEG_LINK
      INTEGER, DIMENSION(:,:), INTENT(INOUT) :: TMPIMP

      ! . Local scalars.
      INTEGER :: I, INDEX, J

      ! . Loop over the impropers.
      DO I = START,STOP

         ! . Loop over the atoms in the impropers.
         DO J = 1,4

            ! . There is a bond to the next residue.
            IF ( TMPIMP(J,I) == 0 ) THEN

               IF ( NEG_LINK(CURRES+1) > 0 ) THEN
                  TMPIMP(J,I) = NEG_LINK(CURRES+1)
               ELSE
                  CALL ERROR ( "GENERATE_STRUCTURE", "Improper specified without corresponding `+R' bond." )
               END IF

            ! . There is a link.
            ELSE IF ( TMPIMP(J,I) < 0 ) THEN

               INDEX = ABS ( TMPIMP(J,I) )
               IF ( LNKHGH(INDEX) > 0 ) THEN
                  TMPIMP(J,I) = LNKHGH(INDEX)
               ELSE
                  CALL ERROR ( "GENERATE_STRUCTURE", "There is a missing `*R' specification in a link." )
               END IF
            END IF
         END DO
      END DO

      END SUBROUTINE RESIDUE_IMPROPERS_LAST

      !---------------------------------------
      SUBROUTINE RESIDUE_INITIALIZE ( RESINF )
      !---------------------------------------

      ! . Scalar arguments.
      TYPE(RESIDUE_INFORMATION), INTENT(INOUT) :: RESINF

      ! . Deallocate the arrays in RESINF.
      DEALLOCATE ( RESINF%CHRGS, RESINF%NAMES, RESINF%TYPES )
      DEALLOCATE ( RESINF%BONDS, RESINF%BINDX )
      DEALLOCATE ( RESINF%IMPRS, RESINF%IINDX )

      ! . Initialize the RESINF scalars.
      RESINF%MM_NATOMS = 0
      RESINF%NBONDS = 0
      RESINF%NIMPRS = 0

      END SUBROUTINE RESIDUE_INITIALIZE

      !---------------------------------------------------
      SUBROUTINE RESIDUE_MODIFY ( RESINF, VARIANT, INDEX )
      !---------------------------------------------------

      ! . Scalar arguments.
      INTEGER,                   INTENT(IN)    :: INDEX
      TYPE(RESIDUE_INFORMATION), INTENT(INOUT) :: RESINF
      TYPE(VARIANT_TYPE),        INTENT(IN)    :: VARIANT

      ! . Local scalars.
      INTEGER                   :: I, J, NATMTMP, NBNDTMP, NIMPTMP, START
      LOGICAL                   :: QFOUND
      TYPE(RESIDUE_INFORMATION) :: RESTMP

      ! . Fill some local counters.
      NATM = RESINF%MM_NATOMS
      NBND = RESINF%NBONDS
      NIMP = RESINF%NIMPRS

      ! . Calculate the number of extra atoms, bonds and impropers to add.
      NATMTMP = MAX ( ( NATM + VARIANT%NADDS - VARIANT%NDELETES ), NATM )
      NBNDTMP = NBND + VARIANT%NBONDS
      NIMPTMP = NIMP + VARIANT%NIMPROPERS

      ! . Allocate space for the new residue information.
      CALL RESIDUE_ALLOCATE ( NATMTMP, NBNDTMP, NIMPTMP, RESTMP )

      ! . Copy over the data.
      CALL RESIDUE_COPY ( RESINF, RESTMP )

      ! . Initialize START.
      START = NATM + 1

      ! . Check to see if atoms are to be deleted.
      IF ( VARIANT%NDELETES > 0 ) THEN

         ! . Loop over the atoms to delete.
         DO I = 1,VARIANT%NDELETES

            ! . Delete the atom from the list.
            NATMTMP = NATM
            QFOUND  = .FALSE.
            NATM    = 0
            DO J = 1,NATMTMP
               IF ( VARIANT%DELATM(I) == RESTMP%NAMES(J) ) THEN
                  QFOUND = .TRUE.
                  START  = MIN ( START, J )
               ELSE
                  NATM = NATM + 1
                  RESTMP%CHRGS(NATM) = RESTMP%CHRGS(J)
                  RESTMP%NAMES(NATM) = RESTMP%NAMES(J)
                  RESTMP%TYPES(NATM) = RESTMP%TYPES(J)
               END IF
            END DO

            ! . Unknown atom.
            IF ( .NOT. QFOUND ) CALL ERROR ( "GENERATE_STRUCTURE", "Unknown atom to delete: " // &
                                                         VARIANT%DELATM(I)(1:LEN_TRIM(VARIANT%DELATM(I)))//"." )

            ! . Remove bonds and impropers concerning the deleted atom.
            NBNDTMP = NBND
            NBND = 0
            DO J = 1,NBNDTMP
               IF ( .NOT. ANY ( RESTMP%BONDS(1:2,J) == VARIANT%DELATM(I) ) ) THEN
                  NBND = NBND + 1
                  RESTMP%BINDX(    NBND) = RESTMP%BINDX(    J)
                  RESTMP%BONDS(1:2,NBND) = RESTMP%BONDS(1:2,J)
               END IF
            END DO

            NIMPTMP = NIMP
            NIMP = 0
            DO J = 1,NIMPTMP
               IF ( .NOT. ANY ( RESTMP%IMPRS(1:4,J) == VARIANT%DELATM(I) ) ) THEN
                  NIMP = NIMP + 1
                  RESTMP%IINDX(    NIMP) = RESTMP%IINDX(    J)
                  RESTMP%IMPRS(1:4,NIMP) = RESTMP%IMPRS(1:4,J)
               END IF
            END DO

         END DO

      END IF

      ! . Check to see if atoms are to be added.
      IF ( VARIANT%NADDS > 0 ) THEN

         ! . Make a space for the data.
         IF ( START <= NATM ) THEN
            DO I = NATM,START,-1
               RESTMP%CHRGS(I+VARIANT%NADDS) = RESTMP%CHRGS(I)
               RESTMP%NAMES(I+VARIANT%NADDS) = RESTMP%NAMES(I)
               RESTMP%TYPES(I+VARIANT%NADDS) = RESTMP%TYPES(I)
            END DO    
         END IF

         ! . Add the data.
         RESTMP%CHRGS(START:VARIANT%NADDS+START-1) = VARIANT%ADDCHG(1:VARIANT%NADDS)
         RESTMP%NAMES(START:VARIANT%NADDS+START-1) = VARIANT%ADDATM(1:VARIANT%NADDS)
         RESTMP%TYPES(START:VARIANT%NADDS+START-1) = VARIANT%ADDTYP(1:VARIANT%NADDS)

         ! . Increment NATM.
         NATM = NATM + VARIANT%NADDS

      END IF

      ! . Check if any charges are to be changed.
      IF ( VARIANT%NCHARGES > 0 ) THEN

         ! . Loop over the charges.
         DO I = 1,VARIANT%NCHARGES

            ! . Loop over the atoms.
            DO J = 1,NATM
               IF ( VARIANT%CHGATM(I) == RESTMP%NAMES(J) ) THEN
                  RESTMP%CHRGS(J) = VARIANT%CHGCHG(I)
                  GO TO 20
               END IF
            END DO

            ! . Type name not found.
            CALL ERROR ( "GENERATE_STRUCTURE", "Unknown atom to change charge: " // &
                                            VARIANT%CHGATM(I)(1:LEN_TRIM(VARIANT%CHGATM(I)))//"." )

            ! . End of loop.
            20 CONTINUE

         END DO

      END IF

      ! . Check if any bonds are to be added.
      IF ( VARIANT%NBONDS > 0 ) THEN

         ! . Add the bonds.
         RESTMP%BINDX(    NBND+1:NBND+VARIANT%NBONDS) = INDEX
         RESTMP%BONDS(1:2,NBND+1:NBND+VARIANT%NBONDS) = VARIANT%BONDS(1:2,1:VARIANT%NBONDS)

         ! . Increment NBND.
         NBND = NBND + VARIANT%NBONDS

      END IF

      ! . Check if any impropers are to be added.
      IF ( VARIANT%NIMPROPERS > 0 ) THEN

         ! . Add the impropers.
         RESTMP%IINDX(    NIMP+1:NIMP+VARIANT%NIMPROPERS) = INDEX
         RESTMP%IMPRS(1:4,NIMP+1:NIMP+VARIANT%NIMPROPERS) = VARIANT%IMPROPERS(1:4,1:VARIANT%NIMPROPERS)

         ! . Increment NIMP.
         NIMP = NIMP + VARIANT%NIMPROPERS

      END IF

      ! . Initialize the residue.
      CALL RESIDUE_INITIALIZE ( RESINF )

      ! . Reallocate space for the residue.
      CALL RESIDUE_ALLOCATE ( NATM, NBND, NIMP, RESINF )

      ! . Copy back the data.
      CALL RESIDUE_COPY ( RESTMP, RESINF )

      ! . Initialize the temporary residue type.
      CALL RESIDUE_INITIALIZE ( RESTMP )

      END SUBROUTINE RESIDUE_MODIFY

   END SUBROUTINE GENERATE_STRUCTURE

   !--------------------------------------------
   SUBROUTINE READ_SEQUENCE_FILE ( TUNIT, FILE )
   !--------------------------------------------

   ! . Scalar arguments.
   INTEGER, INTENT(IN) :: TUNIT

   ! . Optional scalar arguments.
   CHARACTER ( LEN = * ), INTENT(IN), OPTIONAL :: FILE

   ! . Local scalars.
   CHARACTER ( LEN = RESIDUE_NAME_LENGTH ) :: LNKNAM, TMPNAM, VNAME, VRNAM
   INTEGER                                 :: I, IOSTAT, IRES, ISUB, J, K, LINDEX, MM, N, NRESIN, &
                                              NRESTOT, NSUBIN, NVAR, START, STOP, UNIT, VINDEX,   &
                                              VNUM, VGENERIC, VSPECIFIC

   ! . Local arrays.
   CHARACTER ( LEN = RESIDUE_NAME_LENGTH ), DIMENSION(1:2) :: LNKRES, LNKSEG
   INTEGER,                                 DIMENSION(1:2) :: LNUM

   ! . Local allocatable arrays.
   CHARACTER ( LEN =   RESIDUE_NAME_LENGTH ), ALLOCATABLE, DIMENSION(:) :: RNAME
   CHARACTER ( LEN = SUBSYSTEM_NAME_LENGTH ), ALLOCATABLE, DIMENSION(:) :: SNAME
   INTEGER,                                   ALLOCATABLE, DIMENSION(:) :: SINDEX, SUBVAR

   !-------------------------------------------------------------------------
   ! . Get the unit number for the sequence information.
   !-------------------------------------------------------------------------
   ! . Check to see if the sequence file argument is present.
   IF ( PRESENT ( FILE ) ) THEN

      ! . Get the next unit number.
      UNIT = NEXT_UNIT ( )

      ! . Open the file.
      OPEN ( UNIT, FILE = FILE, ACTION = "READ", STATUS = "OLD", IOSTAT = IOSTAT )

      ! . Check for an error.
      IF ( IOSTAT /= 0 ) CALL ERROR ( "READ_SEQUENCE_FILE", "I/O Error.", IOSTAT )
   
   ! . The unit for reading is the input stream.
   ELSE
      UNIT = INPUT
   END IF

   ! . Set the parsing unit.
   CALL PUSH_UNIT ( UNIT )

   ! . Rewind the scratch file.
   REWIND ( TUNIT )

   !-------------------------------------------------------------------------
   ! . Read in the subsystem residue and variant data from the sequence file.
   !-------------------------------------------------------------------------
   ! . Check the sequence header.
   CALL GET_LINE
   CALL GET_WORD
   IF ( WORD(1:WRDLEN) /= "SEQUENCE" ) THEN
      CALL PARSE_ERROR ( "READ_SEQUENCE_FILE", "SEQUENCE label invalid." )
   END IF

   ! . Read in the number of subsystems.
   CALL GET_LINE
   CALL GET_INTEGER ( NSUBIN )

   ! . Allocate space for the subsystem arrays.
   ALLOCATE ( SINDEX(1:NSUBIN+1), SNAME(1:NSUBIN), SUBVAR(1:NSUBIN+1) )
   
   ! . Initialize some counters.
   NRESTOT = 0

   ! . Loop over the subsystems.
   DO ISUB = 1,NSUBIN

      ! . Read in the subsystem header.
      CALL GET_LINE
      CALL GET_WORD
      IF ( WORD(1:WRDLEN) /= "SUBSYSTEM" ) THEN
         CALL PARSE_ERROR ( "READ_SEQUENCE_FILE", "SUBSYSTEM label invalid." )
      END IF

      ! . Get and check the subsystem name.
      CALL GET_WORD
      IF ( ANY ( SNAME(1:ISUB-1) == WORD(1:WRDLEN) ) ) THEN
         CALL PARSE_ERROR ( "READ_SEQUENCE_FILE", "Duplicate subsystem names." )
      END IF

      ! . Store the name.
      SNAME(ISUB) = WORD(1:WRDLEN)

      ! . Fill the subsystem index array.
      SINDEX(ISUB) = NRESTOT
   
      ! . Read in the number of residues in the subsystem.
      CALL GET_LINE
      CALL GET_INTEGER ( NRESIN )

      ! . Increment the total number of residues.
      NRESTOT = NRESTOT + NRESIN

      ! . Allocate space for the residue names.
      ALLOCATE ( RINDEX(1:NRESIN), RNAME(1:NRESIN) )

      ! . Top of the loop to read in the residue names.
      IRES = 0
      5 CONTINUE

         ! . Get the residue name and its frequency of occurrence.
         CALL GET_LINE
         CALL GET_WORD ; TMPNAM = WORD(1:WRDLEN)
         IF ( EMPTY ( ) ) THEN
            N = 1
         ELSE
            CALL GET_INTEGER ( N )
         END IF

         ! . Check the number of residues.
         IF ( ( N <= 0 ) .OR. ( IRES + N ) > NRESIN ) THEN
            CALL PARSE_ERROR ( "READ_SEQUENCE_FILE", "Subsystem residue count error." )
         END IF

         ! . Loop over the possible residue names.
         DO MM = 1,NRESIDUES
            ! . The residue names match.
            IF ( TMPNAM == RESIDUES(MM)%NAME ) THEN
               DO I = IRES+1,IRES+N
                  RINDEX(I) = MM
                  RNAME(I)  = TMPNAM
               END DO
               GO TO 10
            END IF
         END DO

         ! . No name has been found.
         CALL PARSE_ERROR ( "READ_SEQUENCE_FILE", "Unknown residue name: "//WORD(1:WRDLEN)//"." )

         ! . End of loop.
         10 CONTINUE

         ! . Increment IRES and decide whether to go back to the top of the loop.
         IRES = IRES + N
         IF ( IRES < NRESIN ) GO TO 5

      ! . Write out the residue indices and names to the scratch file.
      WRITE ( TUNIT ) RINDEX(1:NRESIN), RNAME(1:NRESIN)

      ! . Read in the next word.
      NVAR = 0
      20 CONTINUE
      CALL GET_LINE
      CALL GET_WORD

      ! . There are variants present.
      IF ( WORD(1:WRDLEN) == "VARIANT" ) THEN

         ! . Increment the number of variants.
         NVAR = NVAR + 1

         ! . Read in the details of the variant.
         CALL GET_WORD
         VNAME = WORD(1:WRDLEN)
         CALL GET_WORD
         VRNAM = WORD(1:WRDLEN)
         CALL GET_INTEGER ( VNUM )

         ! . Check the residue name and number.
         IF ( ( VNUM < 1 ) .OR. ( VNUM > NRESIN ) .OR. ( VRNAM /= RNAME(VNUM) ) ) THEN
            CALL PARSE_ERROR ( "READ_SEQUENCE_FILE", "Invalid residue name or number in VARIANT statement." )
         END IF

         ! . Search for a specific or a generic variant name.
         VSPECIFIC = 0
         VGENERIC  = 0
         DO MM = 1,NVARIANTS
            IF ( ( VNAME == VARIANTS(MM)%NAME ) .AND. ( VRNAM == VARIANTS(MM)%RESIDUE_NAME ) ) THEN
               VSPECIFIC = MM
            ELSE IF ( ( VNAME == VARIANTS(MM)%NAME ) .AND. ( VARIANTS(MM)%RESIDUE_NAME == " " ) ) THEN
               VGENERIC = MM
            END IF
         END DO

         ! . Check the names.
         IF ( VSPECIFIC > 0 ) THEN
            VINDEX = VSPECIFIC
         ELSE IF ( VGENERIC > 0 ) THEN
            VINDEX = VGENERIC
         ELSE
            CALL PARSE_ERROR ( "READ_SEQUENCE_FILE", "Unknown variant name: "//VNAME(1:LEN_TRIM(VNAME))//"." )
         END IF

         ! . Write out the variant to the temporary file.
         WRITE ( TUNIT ) VINDEX, VNUM

         ! . Read in the next line.
         GO TO 20

      ! . The subsystem has finished.
      ELSE IF ( WORD(1:WRDLEN) == "END" ) THEN

         ! . Deallocate the residue arrays.
         DEALLOCATE ( RINDEX, RNAME )

         ! . Save the number of variants for the subsystem.
         SUBVAR(ISUB) = NVAR

      ! . Unrecognized command.
      ELSE
         CALL PARSE_ERROR ( "READ_SEQUENCE_FILE", "Invalid subsystem block terminator." )
      END IF

   END DO

   ! . Fill the final element in the subsystem index array.
   SINDEX(NSUBIN+1) = NRESTOT

   !----------------------------------------------------------------------------
   ! . Fill some parts of the SEQUENCE data structure.
   !----------------------------------------------------------------------------
   ! . Allocate the data structures.
   CALL SEQUENCE_ALLOCATE ( NRESTOT, NSUBIN )

   ! . Fill the subsystem arrays.
   SUBIND = SINDEX
   SUBNAM = SNAME

   ! . Calculate the total number of variants.
   TOTVAR = SUM ( SUBVAR(1:NSUBSYS) )

   ! . Allocate some temporary arrays.
   ALLOCATE ( RINDEX(1:NRESID), VARIND(1:TOTVAR), VARRES(1:TOTVAR) )

   ! . Rewind the scratch file.
   REWIND ( TUNIT )

   ! . Initialize some counters.
   NVAR = 0

   ! . Read in the residue and variant data from the scratch file.
   DO ISUB = 1,NSUBSYS

      ! . Find the starting and stopping points for the read.
      START = SUBIND(ISUB)+1
      STOP  = SUBIND(ISUB+1)

      ! . Read in the residue data.
      READ ( TUNIT ) RINDEX(START:STOP), RESNAM(START:STOP)

      ! . Loop over the variants for the subsystem.
      DO I = 1,SUBVAR(ISUB)

         ! . Increment NVAR.
         NVAR = NVAR + 1

         ! . Read in the variant data.
         READ ( TUNIT ) VARIND(NVAR), VNUM

         ! . Increment VARNUM.
         VARRES(NVAR) = VNUM + START - 1

      END DO

   END DO

   ! . Deallocate some temporary arrays.
   DEALLOCATE ( SINDEX, SNAME, SUBVAR )

   !-------------------------------------------------------------------------
   ! . Read in the link data from the sequence file.
   !-------------------------------------------------------------------------
   ! . Rewind the scratch file.
   REWIND ( TUNIT )

   ! . Read the next line in the file.
   TOTLNK = 0
   30 CONTINUE
   CALL GET_LINE
   CALL GET_WORD

   ! . There are links present.
   IF ( WORD(1:WRDLEN) == "LINK" ) THEN

      ! . Increment the number of links.
      TOTLNK = TOTLNK + 1

      ! . Read in the details of the link
      CALL GET_WORD
      LNKNAM = WORD(1:WRDLEN)
      CALL GET_WORD
      LNKSEG(1) = WORD(1:WRDLEN)
      CALL GET_WORD
      LNKRES(1) = WORD(1:WRDLEN)
      CALL GET_INTEGER ( LNUM(1) )
      CALL GET_WORD
      LNKSEG(2) = WORD(1:WRDLEN)
      CALL GET_WORD
      LNKRES(2) = WORD(1:WRDLEN)
      CALL GET_INTEGER ( LNUM(2) )

      ! . Loop over the possible link names.
      DO MM = 1,NRESLINK
         IF ( LNKNAM == RESLINK(MM)%NAME ) THEN
            LINDEX = MM
            GO TO 40
         END IF
      END DO

      ! . No name has been found.
      CALL PARSE_ERROR ( "READ_SEQUENCE_FILE", "Unknown link name: "//LNKNAM(1:LEN_TRIM(LNKNAM))//"." )

      ! . End of loop.
      40 CONTINUE

      ! . Check the subsystem and residue data.
      DO J = 1,2

         ! . Loop over the possible subsystem names.
         DO K = 1,NSUBSYS
            IF ( LNKSEG(J) == SUBNAM(K) ) THEN
               ISUB = K
               GO TO 50
            END IF
         END DO

         ! . No name has been found.
         CALL PARSE_ERROR ( "READ_SEQUENCE_FILE", "Unknown link subsystem name: "//LNKSEG(J)(1:LEN_TRIM(LNKSEG(J)))//"." )

         ! . End of loop.
         50 CONTINUE

         ! . Increment LNKNUM.
         LNUM(J) = LNUM(J) + SUBIND(ISUB)

         ! . Check the residue name and number.
         IF ( ( LNUM(J) <= SUBIND(ISUB) ) .OR. ( LNUM(J) > SUBIND(ISUB+1) ) .OR. &
              ( LNKRES(J) /= RESNAM(LNUM(J)) ) ) THEN
            CALL PARSE_ERROR ( "READ_SEQUENCE_FILE", "Invalid residue name or number in LINK statement." )
         END IF

      END DO

      ! . Write out the link data.
      WRITE ( TUNIT ) LINDEX, LNUM(1), LNUM(2)

      ! . Read in the next line.
      GO TO 30

   ! . The sequence has finished.
   ELSE IF ( WORD(1:WRDLEN) == "END" ) THEN
   
      ! . Close the sequence file.
      CLOSE ( UNIT )

      ! . Reset the parsing unit.
      CALL POP_UNIT

   ! . Unrecognized command.
   ELSE
      CALL PARSE_ERROR ( "READ_SEQUENCE_FILE", "Invalid sequence block terminator." )
   END IF

   ! . Allocate some temporary arrays.
   ALLOCATE ( LNKIND(1:TOTLNK), LNKNUM(1:2,1:TOTLNK) )

   ! . Rewind the scratch file.
   REWIND ( TUNIT )

   ! . Read in the link data.
   DO I = 1,TOTLNK
      READ ( TUNIT ) LNKIND(I), LNKNUM(1,I), LNKNUM(2,I)
   END DO

   !----------------------------------------------------------------------------
   ! . Finish up.
   !----------------------------------------------------------------------------
   ! . Write out some information.
   IF ( UNIT == INPUT ) THEN
      WRITE ( OUTPUT, "(/A)" ) "Sequence data read from input stream."
   ELSE
      WRITE ( OUTPUT, "(/A)" ) "Sequence data read from " // FILE
   END IF

   END SUBROUTINE READ_SEQUENCE_FILE

END MODULE MM_SYSTEM
