PROGRAM MakePFFT
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
!    Set up the multipole arrays
     CALL MultipoleSetUp()
!    Allocate the tensors
     MaxEll=GM%PBC%PFFMaxEll
     CALL New(TenC,LSP(2*MaxEll),0)
     CALL New(TenS,LSP(2*MaxEll),0)
!    Calculate the tensors ...
     CALL CalculatePFFT(MaxEll,GM,Args,CS_IN,TenC,TenS)
!    Put them to HDF
     CALL Put(TenC,'PFFTensorC')
     CALL Put(TenS,'PFFTensorS')
!    Delete
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
!
END PROGRAM MakePFFT
