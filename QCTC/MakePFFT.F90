#ifdef PARALLEL_CLONES
PROGRAM MakePFFT
#ifdef PERIODIC
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE ProcessControl
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE BraBloks
  USE PFFTen
  USE AtomPairs
  IMPLICIT NONE
  TYPE(TIME)                     :: TimePFFT
  TYPE(CRDS)                     :: GM
  TYPE(DBL_VECT)                 :: TenC,TenS
  TYPE(CellSet)                  :: CS
  TYPE(ARGMT)                    :: Args
  INTEGER                        :: MaxEll
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
  CHARACTER(LEN=8),PARAMETER     :: Prog='MakePFFT'
  !-------------------------------------------------------------------------------- 
  ! Start up macro
  CALL StartUp(Args,Prog)
  ! Get the geometry
  CALL Get(GM,Tag_O=CurGeom)
  IF(GM%PBC%PFFMaxEll>FFELL) THEN
     CALL MondoHalt(0,'MaxELL > FFELL Halting in CalculatePFFT')
  ELSEIF(GM%PBC%Dimen>0)THEN
     ! Set up the multipole arrays
     CALL MultipoleSetUp()
     ! Allocate the tensors
     MaxEll=GM%PBC%PFFMaxEll
     CALL New(TenC,LSP(2*MaxEll),0)
     CALL New(TenS,LSP(2*MaxEll),0)
     ! Calculate the tensors ...
     CALL CalculatePFFT(MaxEll,GM,Args,CS_OUT,TenC,TenS)
     ! ... and put them to HDF
     CALL Put(TenC,'PFFTensorC')
     CALL Put(TenS,'PFFTensorS')
     ! Delete
     CALL Delete(GM)
     CALL Delete(Args)
     CALL Delete(TenC)
     CALL Delete(TenS)
  ENDIF
  ! didn't count flops, any accumulation is residual
  ! from matrix routines
  PerfMon%FLOP=Zero 
  ! Shutdown 
  CALL ShutDown(Prog)
#endif
END PROGRAM MakePFFT
#else
PROGRAM MakePFFT
#ifdef PERIODIC
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE ProcessControl
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE BraBloks
  USE PFFTen
  USE AtomPairs
  IMPLICIT NONE
!
  TYPE(TIME)                     :: TimePFFT
  TYPE(CRDS)                     :: GM
  TYPE(ARGMT)                    :: Args
  CHARACTER(LEN=8),PARAMETER     :: Prog='MakePFFT'
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
  INTEGER                        :: MaxEll, I
!-------------------------------------------------------------------------------- 
! Start up macro
  CALL StartUp(Args,Prog)
! Get basis set and geometry and Set up and Store the outer cell sum
#ifdef MMech
  IF(HasMM()) THEN
    CALL Get(GM,Tag_O='GM_MM'//CurGeom)
    CALL SetCellNumber(GM,1.D-10)
  ELSE
    CALL Get(GM,Tag_O=CurGeom)
    CALL SetCellNumber(GM)
  ENDIF
#else
  CALL Get(GM,Tag_O=CurGeom)
  CALL SetCellNumber(GM,CS_OUT)
#endif 
  CALL Put(CS_OUT,'CS_OUT',Tag_O=CurBase) 
! Set Up the Multipoles  
  CALL MultipoleSetUp()
! Allocate memory
  MaxEll = GM%PBC%PFFMaxEll
  IF(MaxELL > FFELL) THEN
     CALL MondoHalt(0,'MaxELL > FFELL Halting in CalculatePFFT')
  ENDIF
! Calculate and Store the Tensors  
  CALL CalculatePFFT(MaxEll,GM,Args)
! Output
  CALL PPrint(CS_IN, 'inner sum',Prog,Unit_O=6)
  CALL PPrint(CS_OUT,'outer sum',Prog,Unit_O=6)
! Delete
  CALL Delete(GM)
  CALL Delete(Args)
  CALL Delete(TensorC)
  CALL Delete(TensorS)
! didn't count flops, any accumulation is residual
! from matrix routines
  PerfMon%FLOP=Zero 
! Shutdown 
  CALL ShutDown(Prog)
#endif
END PROGRAM MakePFFT
#endif
