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
!
  TYPE(TIME)                     :: TimePFFT
  TYPE(BSET)                     :: BS
  TYPE(CRDS)                     :: GM
  TYPE(ARGMT)                    :: Args
  CHARACTER(LEN=8),PARAMETER     :: Prog='MakePFFT'
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
  INTEGER                        :: MaxEll
!-------------------------------------------------------------------------------- 
! Start up macro
  CALL StartUp(Args,Prog)
! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
! Set Up the Multipoles  
  CALL MultipoleSetUp(FFEll2)
! Allocate memory
  MaxEll = GM%PBC%PFFMaxEll+BS%NASym+1
  IF(MaxELL > FFELL) THEN
     CALL MondoHalt(0,'MaxELL > FFELL Halting in CalculatePFFT')
  ENDIF
  CALL New(TensorC,LSP(MaxEll),0)
  CALL New(TensorS,LSP(MaxEll),0)
! Calculate and Store the Tensors  
  CALL CalculatePFFT(MaxEll,GM,BS,Args)
! Delete
  CALL Delete(BS)
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
