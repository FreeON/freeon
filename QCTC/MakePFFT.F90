!    FAST O(N lg N) COMPUTATION OF THE COULOMB MATRIX
!    Authors:  Matt Challacombe and CJ Tymczak
!==============================================================================
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
  INTEGER                        :: MaxEll
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
  CALL SetCellNumber(GM)
#endif 
  CALL Put_CellSet(CS_OUT,'CS_OUT'//CurBase//CurGeom) 
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
  CALL PPrint(CS_IN, 'inner sum',Prog)
  CALL PPrint(CS_OUT,'outer sum',Prog)
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
