MODULE SCFLocals
   USE DerivedTypes
   USE GlobalCharacters
   USE Parse
   IMPLICIT NONE
!--------------------------------------------------
!  Global status variables
!
   INTEGER,DIMENSION(3) :: Current
   INTEGER,DIMENSION(3) :: Previous
   INTEGER              :: PCyc,CCyc,NCyc   
   CHARACTER(LEN=3)     :: PrvCycl,CurCycl,SCFCycl,NxtCycl
   INTEGER              :: PBas,CBas
   CHARACTER(LEN=3)     :: PrvBase,CurBase
   INTEGER              :: PGeo,CGeo,NGeo
   CHARACTER(LEN=3)     :: PrvGeom,CurGeom,NxtGeom
   CHARACTER(LEN=20)    :: SCFActn
!-----------------------------------------------------------
   INTEGER,               PARAMETER :: DCL=DEFAULT_CHR_LEN
   CHARACTER(LEN=DCL),    PARAMETER :: False  = '.FALSE.'   ! False flag
   CHARACTER(LEN=DCL),    PARAMETER :: True   = '.TRUE.'    ! True flag
   CHARACTER(LEN=DCL),    PARAMETER :: Direct = 'Direct'    ! Do direct SCF
   CHARACTER(LEN=DCL),    PARAMETER :: Restart= 'Restart'   ! Restart from InfoFile
   CHARACTER(LEN=DCL),    PARAMETER :: Switch = 'Switch'    ! Do a basis set switch
   CHARACTER(LEN=DCL),    PARAMETER :: InkFok = 'InkFok'    ! Do incremental Fock build
   CHARACTER(LEN=DCL),    PARAMETER :: Core   = 'Core'      ! Use core Hamiltontian guess
   CHARACTER(LEN=DCL),    PARAMETER :: Minimal= 'Minimal'   ! Use superposition of AO DMs guess
   CHARACTER(LEN=DCL),    PARAMETER :: Blank  =' '
   CHARACTER(LEN=DCL),DIMENSION(9)  :: CtrlVect
!-------------------------------------------------  
!  SCF Dimensions
   INTEGER, PARAMETER         :: MaxSets=10
   INTEGER, PARAMETER         :: MaxSCFs=64
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
      INTEGER,DIMENSION(3)               :: Previous
      INTEGER,DIMENSION(3)               :: Current
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
!
      INTEGER                            :: Guess
!
      INTEGER                            :: Grad
      REAL(DOUBLE),DIMENSION(2)          :: MDVar
!
   END TYPE
!------------------------------------------------------------------------------------------------  
!  Thresholds (Loose ~4 digits, Good ~6 digits, Tight ~8 digits, VeryTight ~10 digits):
!
   REAL(DOUBLE),DIMENSION(4) :: TrixNeglect=(/1.D-4, 1.D-5, 1.D-6,  1.D-7 /)
   REAL(DOUBLE),DIMENSION(4) :: CubeNeglect=(/1.D-4, 1.D-6, 1.D-8,  1.D-10 /)
   REAL(DOUBLE),DIMENSION(4) :: TwoENeglect=(/1.D-6, 1.D-8, 1.D-10, 1.D-12/)
   REAL(DOUBLE),DIMENSION(4) :: DistNeglect=(/1.D-8, 1.D-10,1.D-12, 1.D-14/)
   REAL(DOUBLE),DIMENSION(4) :: ETol       =(/1.D-5, 1.D-7, 1.D-9,  1.D-11/)
   REAL(DOUBLE),DIMENSION(4) :: DTol       =(/1.D-2, 1.D-3, 1.D-4,  1.D-6 /)
!-----------------------------------------------------------------------------------------------
!  Asymptotic dimensioning parameters for memory limits of BCSR and DBCSR matrices
!
   REAL(DOUBLE),DIMENSION(4) :: BandWidth=(/ 1.D3, 1.D3,1.3D3,1.6D3/)
   REAL(DOUBLE),DIMENSION(4) :: BWDecay  =(/1.D-4,1.D-4,1.D-3,1.D-2/)
!-----------------------------------------------------------------------------------------------
!  Global variables for file control
!
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: SCF_NAME
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: PROCESS_ID 
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: HOST_NAME
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: MACHINE_ID
#if defined(PARALLEL) && !defined(MPI2) 
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: MPI_INVOKE
   CHARACTER(LEN=DEFAULT_CHR_LEN),SAVE :: MPI_FLAGS
#endif
  CONTAINS
     SUBROUTINE SetGlobalCtrlIndecies(Ctrl)
        TYPE(SCFControls),        INTENT(IN) :: Ctrl
        Previous=Ctrl%Previous
        Current=Ctrl%Current
!       Set SCF cycle indecies
        PCyc=Previous(1)
        CCyc=Current(1)
        NCyc=CCyc+1
        PrvCycl=TRIM(IntToChar(PCyc))
        CurCycl=TRIM(IntToChar(CCyc))
        SCFCycl=TRIM(IntToChar(CCyc))
        NxtCycl=TRIM(IntToChar(NCyc))
!       Set basis set indecies
        PBas=Previous(2)
        CBas=Current(2)
        PrvBase=TRIM(IntToChar(PBas))
        CurBase=TRIM(IntToChar(CBas))
!       Set geometry indecies
        PGeo=Previous(3)
        CGeo=Current(3)
        NGeo=CGeo+1
        PrvGeom=TRIM(IntToChar(PGeo))
        CurGeom=TRIM(IntToChar(CGeo))
        NxtGeom=TRIM(IntToChar(NGeo))
     END SUBROUTINE SetGlobalCtrlIndecies
!
     FUNCTION SetCtrlVect(Ctrl,Actn1_O,Actn2_O) RESULT(CVect)
        TYPE(SCFControls),        INTENT(IN) :: Ctrl
        CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Actn1_O,Actn2_O
        CHARACTER(LEN=DCL),DIMENSION(9)      :: CVect
        CVect(1)=Ctrl%Name
        CVect(2)=" "
        CVect(3)=" "
        IF(PRESENT(Actn1_O))CVect(2)=Actn1_O
        IF(PRESENT(Actn2_O))CVect(3)=Actn2_O
        CVect(4:9)=(/(IntToChar(Ctrl%Current(1))), &
                     (IntToChar(Ctrl%Current(2))), &
                     (IntToChar(Ctrl%Current(3))), &
                     (IntToChar(Ctrl%Previous(1))),&
                     (IntToChar(Ctrl%Previous(2))),&
                     (IntToChar(Ctrl%Previous(3))) /)              
     END FUNCTION SetCtrlVect
!
END MODULE

