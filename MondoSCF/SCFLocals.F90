MODULE SCFLocals
   USE DerivedTypes
   USE GlobalCharacters
   IMPLICIT NONE
!-------------------------------------------------  
!  SCF Dimensions
   INTEGER, PARAMETER         :: MaxSets=10
   INTEGER, PARAMETER         :: MaxSCFs=50
!  Clean up control
   LOGICAL, PARAMETER         :: TidyFiles=.TRUE.
!------------------------------------------------------------------------------------------------
!  SCF convergence keys
!
   INTEGER, PARAMETER :: SCF_OSCILATOR=10302 ! SCF is out of control (oscilating)
   INTEGER, PARAMETER :: SCF_STAGNATED=30875 ! SCF stagnated
   INTEGER, PARAMETER :: SCF_CONVERGED=42341 ! SCF converged to requested accuracy
!------------------------------------------------------------------------------------------------  
!  Object to control the SCF
   TYPE SCFControls
      CHARACTER(LEN=DEFAULT_CHR_LEN)     :: Info   ! Info file 
      CHARACTER(LEN=DEFAULT_CHR_LEN)     :: Name   ! SCF name 
      TYPE(INT_VECT)                     :: Previous
      TYPE(INT_VECT)                     :: Current
      INTEGER                            :: NSet   ! Number of basis sets 
      INTEGER                            :: NGeom  ! Number of configurations 
      INTEGER                            :: NCGC   ! Number of CG cycles
      LOGICAL                            :: SuperP ! Start from the superposition of atomic densitites
      LOGICAL                            :: InkFok ! Incremental Fock builds with difference densities
      LOGICAL                            :: DIIS   ! Use direct inversion in iterative subspace
      LOGICAL                            :: Rest   ! Restart from a specified Info file
      INTEGER                            :: Fail   ! SCF convergence key
      CHARACTER(LEN=DEFAULT_CHR_LEN), &
                   DIMENSION(MaxSets)    :: BName  ! Basis set name
      INTEGER,     DIMENSION(MaxSets)    :: Method ! SCF method
      INTEGER,     DIMENSION(MaxSets)    :: Model  ! Model chemistry
      INTEGER,     DIMENSION(MaxSets)    :: AccL   ! Accuracy level to use 
      INTEGER,     DIMENSION(MaxSets)    :: NCyc   ! Number of SCF cycles taken for each set
      LOGICAL,     DIMENSION(MaxSets)    :: Cnvrgd ! Convergence flag
      REAL(DOUBLE),DIMENSION(MaxSets)    :: EErr   ! Relative error in total energy
      REAL(DOUBLE),DIMENSION(0:MaxSCFs,5):: Stats  ! Statistics     
   END TYPE
!------------------------------------------------------------------------------------------------  
!  Thresholds (Loose,Good,Tight,VeryTight):

   REAL(DOUBLE),DIMENSION(4) :: CubeNeglect=(/1.D-4, 1.D-5, 1.D-6,  1.D-8 /)
   REAL(DOUBLE),DIMENSION(4) :: TrixNeglect=(/1.D-2, 1.D-3, 1.D-4,  1.D-6 /)
   REAL(DOUBLE),DIMENSION(4) :: TwoENeglect=(/1.D-5, 1.D-6, 1.D-7,  1.D-9 /)
   REAL(DOUBLE),DIMENSION(4) :: DistNeglect=(/1.D-6, 1.D-8, 1.D-10, 1.D-12/)
!
   REAL(DOUBLE),DIMENSION(4) :: ETol       =(/1.D-6, 1.D-8, 1.D-10, 1.D-12/)
   REAL(DOUBLE),DIMENSION(4) :: DTol       =(/1.D-2, 1.D-3, 1.D-4,   1.D-6/)
!-----------------------------------------------------------------------------------------------
!  Asymptotic dimensioning parameters for memory limits (needs work)
   REAL(DOUBLE),DIMENSION(4) :: BandWidth=(/7.D2 ,7.D2 ,7.D2 ,7.D2 /)
   REAL(DOUBLE),DIMENSION(4) :: BWDecay  =(/1.D-6,1.D-6,1.D-6,1.D-6/)
!-------------------------------------------------  
!  Global variables for file control
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: MONDO_PWD 
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: MONDO_HOME   
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: MONDO_SCRATCH
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: MONDO_EXEC
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: SCF_NAME
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: PROCESS_ID 
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: HOST_NAME
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: MACHINE_ID
#ifdef PARALLEL
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: MPI_INVOKE
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: MPI_FLAGS
#endif
!  Globa 
END MODULE

