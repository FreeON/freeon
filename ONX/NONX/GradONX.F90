PROGRAM GradONX
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE Clock
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE ONXParameters
  USE ContractionScaling
  USE ONXMemory
  USE MatFilter
  USE InitExchangeMatrix
  IMPLICIT NONE
  TYPE(BCSR)          :: D
  TYPE(DBL_RNK2)      :: XFrc
  TYPE(BSET)          :: BS
  TYPE(CRDS)          :: GM
  TYPE(ARGMT)         :: Args
  TYPE(DBuf)          :: DB        ! distribution buffers
  TYPE(IBuf)          :: IB        ! 2-e eval buffers
  TYPE(DSL)           :: SB        ! distribution pointers
  TYPE(IDrv)          :: Drv       ! VRR/contraction drivers
  TYPE(INT_VECT)      :: NameBuf
  TYPE(INT_RNK2)      :: SubInd
  TYPE(INT_RNK2)      :: BfnInd
!--------------------------------------------------------------------------------
! Large buffers to hold the sorted distribution and multipole data.
!--------------------------------------------------------------------------------
  TYPE(DBL_VECT)         :: PrmBuf,DisBuf,CBuf,SBuf
  TYPE(DBL_RNK2)         :: VecBuf
  TYPE(INT_VECT)         :: BfnBuf
  TYPE(INT_RNK2)         :: iSA,iSrt
!--------------------------------------------------------------------------------
! Misc. variables and parameters...
!--------------------------------------------------------------------------------
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: InFile
  CHARACTER(LEN=6),PARAMETER     :: Prog='XForce'
!--------------------------------------------------------------------------------
  CALL StartUp(Args,Prog)
  InFile=TRIM(SCFName)//'_Cyc'//TRIM(IntToChar(Args%i%i(1)))
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
  CALL Get(BSiz,'atsiz',Tag_O=CurBase)
  CALL Get(OffS,'atoff',Tag_O=CurBase)
  CALL Get(NBasF,'nbasf',Tag_O=CurBase)
  CALL Get(D,TrixFile('D',Args,0)) !InFile,'.D')
  CALL TrnMatBlk(BS,GM,D)
  Gradient=.TRUE.
  CALL New(NameBuf,NAtoms)
  CALL New(BfnInd,(/BS%NAtms,BS%NCtrt/))
  CALL New(XFrc,(/3,NAtoms/))
  XFrc%D=0.0D0
  CALL New(SubInd,(/3,NBasF/))
  CALL InitSubInd(BS,GM,SubInd)
!--------------------------------------------------------------------------------
! Compute and sort the distribution buffers
!--------------------------------------------------------------------------------
  CALL RangeOfDensity(D,NameBuf,BfnInd,DB,BS,GM)
  DO WHILE (ErrorCode/=eAOK) 
    CALL MemInit(DB,IB,SB,Drv,BS,BS)
    CALL DisOrderGrad(BS,GM,DB,IB,SB,Drv,NameBuf)
  END DO
!--------------------------------------------------------------------------------
! All set to compute the exchange gradient
!--------------------------------------------------------------------------------
  CALL ComputeXForce(BS,GM,D,XFrc,DB,IB,SB,Drv,SubInd,BfnInd)
!--------------------------------------------------------------------------------
! Clean up...
!--------------------------------------------------------------------------------
  CALL Delete(D)
  CALL Delete(DB)
  CALL Delete(IB)
  CALL Delete(SB)
  CALL Delete(Drv)
  CALL Delete(SubInd)
  CALL Delete(NameBuf)
  CALL Delete(BfnInd)
  CALL Delete(BS)
  CALL Delete(GM)
  call halt('enough for now')
  CALL ShutDown(Prog)
END PROGRAM GradONX

