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
#ifdef MMech
  IF(MMOnly()) THEN
    CALL Get(GM,Tag_O='GM_MM'//CurGeom)
    MaxEll = GM%PBC%PFFMaxEll+1
  ELSE IF(HasMM()) THEN
    CALL Get(BS,Tag_O=CurBase)
    CALL Get(GM,Tag_O='GM_MM'//CurGeom)
    MaxEll = GM%PBC%PFFMaxEll+BS%NASym+1
  ELSE
    CALL Get(BS,Tag_O=CurBase)
    CALL Get(GM,Tag_O=CurGeom)
    MaxEll = GM%PBC%PFFMaxEll+BS%NASym+1
  ENDIF
#else
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
  MaxEll = GM%PBC%PFFMaxEll+BS%NASym+1
#endif
! Set Up the Multipoles  
  CALL MultipoleSetUp(FFEll2)
! Allocate memory
  IF(MaxELL > FFELL) THEN
     CALL MondoHalt(0,'MaxELL > FFELL Halting in CalculatePFFT')
  ENDIF
  CALL New(TensorC,LSP(2*MaxEll),0)
  CALL New(TensorS,LSP(2*MaxEll),0)
! Calculate and Store the Tensors  
  CALL CalculatePFFT(MaxEll,GM,Args)
! Delete
#ifdef MMech
  IF(HasQM()) THEN
    CALL Delete(BS)
  ENDIF
#else
  CALL Delete(BS)
#endif
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
