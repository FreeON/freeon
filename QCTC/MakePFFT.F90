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
  TYPE(DBL_RNK2)                 :: BoxShape
  TYPE(DBL_RNK3)                 :: dTenC,dTenS
  TYPE(CellSet)                  :: CS
  TYPE(ARGMT)                    :: Args
  INTEGER                        :: MaxEll,I,J,K
  REAL(DOUBLE)                   :: DDelta 
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
  CHARACTER(LEN=8),PARAMETER     :: Prog='MakePFFT'
!-------------------------------------------------------------------------------- 
! Start up macro
  CALL StartUp(Args,Prog)
! Get the geometry
  CALL Get(GM,Tag_O=CurGeom)
  IF(GM%PBC%Dimen>0)THEN
     IF(GM%PBC%PFFMaxEll>FFELL) THEN
        CALL MondoHalt(0,'MaxELL > FFELL Halting in CalculatePFFT')
     ENDIF
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
!
!    Do the Derivaltive Arrays
!
!    Allocate the derivative tensors
     MaxEll=GM%PBC%PFFMaxEll
     CALL New(dTenC,(/LSP(2*MaxEll),3,3/),(/0,1,1/))
     CALL New(dTenS,(/LSP(2*MaxEll),3,3/),(/0,1,1/))
     CALL New(BoxShape,(/3,3/))
     CALL New(CS,CS_IN%NCells)
     DDelta = 1.D-8
!    Initialize CS
     BoxShape%D=GM%PBC%BoxShape%D
     dTenC%D=Zero
     dTenS%D=Zero
     DO K=1,CS_IN%NCells
        CS_IN%CellCarts%D(:,K) = AtomToFrac(GM,CS_IN%CellCarts%D(:,K))
     ENDDO
!    Calculate the tensors numerically
     DO I=1,3
        DO J=1,3
           IF(GM%PBC%AutoW%I(I)==1 .AND. GM%PBC%AutoW%I(J)==1) THEN
              GM%PBC%BoxShape%D(I,J) = BoxShape%D(I,J) + DDelta
              GM%PBC%InvBoxSh%D      = InverseMatrix(GM%PBC%BoxShape%D)
              DO K=1,CS_IN%NCells
                 CS%CellCarts%D(:,K) = FracToAtom(GM,CS_IN%CellCarts%D(:,K))
              ENDDO
!
              CALL CalculatePFFT(MaxEll,GM,Args,CS,TenC,TenS)
              dTenC%D(:,I,J) = TenC%D(:)
              dTenS%D(:,I,J) = TenS%D(:)
!
              GM%PBC%BoxShape%D(I,J) = BoxShape%D(I,J) - DDelta  
              GM%PBC%InvBoxSh%D      = InverseMatrix(GM%PBC%BoxShape%D)
              DO K=1,CS_IN%NCells
                 CS%CellCarts%D(:,K) = FracToAtom(GM,CS_IN%CellCarts%D(:,K))
              ENDDO
!
              CALL CalculatePFFT(MaxEll,GM,Args,CS,TenC,TenS)
              dTenC%D(:,I,J) = (TenC%D(:)-dTenC%D(:,I,J))/(Two*DDelta)
              dTenS%D(:,I,J) = (TenS%D(:)-dTenS%D(:,I,J))/(Two*DDelta)
           ENDIF
        ENDDO
     ENDDO
!    Put them to HDF
     CALL Put(dTenC,'dPFFTensorC')
     CALL Put(dTenS,'dPFFTensorS')
!    Delete
     CALL Delete(TenC)
     CALL Delete(TenS)  
     CALL Delete(dTenC)
     CALL Delete(dTenS)
  ENDIF
! Delete
  CALL Delete(GM)
  CALL Delete(Args)
! didn't count flops, any accumulation is residual
! from matrix routines
  PerfMon%FLOP=Zero 
! Shutdown 
  CALL ShutDown(Prog)
!
END PROGRAM MakePFFT
