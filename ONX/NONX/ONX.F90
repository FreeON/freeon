PROGRAM ONX
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
  USE ONXMemory
  USE ContractionScaling
  USE MatFilter
  USE InitExchangeMatrix
#ifdef PARALLEL
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef PARALLEL
  TYPE(DBCSR)         :: D
  TYPE(DBCSR)         :: K
  TYPE(BCSR)          :: KTotal
#else
  TYPE(BCSR)          :: D
  TYPE(BCSR)          :: K
  INTEGER             :: MyID=0
#endif
  TYPE(BSET)          :: BSc
  TYPE(CRDS)          :: GMc
  TYPE(BSET)          :: BSp
  TYPE(CRDS)          :: GMp
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
  INTEGER                :: iSwitch
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: InFile
  CHARACTER(LEN=3),PARAMETER     :: Prog='ONX'
!--------------------------------------------------------------------------------
  CALL StartUp(Args,Prog)
  InFile=TRIM(ScrName)//'_Cyc'//TRIM(IntToChar(Args%i%i(1)))
  CALL Get(BSc,Tag_O=CurBase)
  CALL Get(GMc,Tag_O=CurGeom)
  CALL Get(BSp,Tag_O=PrvBase)
  CALL Get(GMp,Tag_O=PrvGeom)
  CALL New(NameBuf,NAtoms)
  CALL Get(BSiz,'atsiz',Tag_O=PrvBase)
  CALL Get(OffS,'atoff',Tag_O=PrvBase)
  CALL Get(NBasF,'nbasf',Tag_O=PrvBase)
!----------------------------------------------
! Get the Density Matrix
!----------------------------------------------
  IF(SCFActn=='InkFok')THEN
     CALL Get(D,TrixFile('DeltaD',Args,0))
  ELSEIF(SCFActn=='BasisSetSwitch')THEN
     CALL Get(D,TrixFile('D',Args,-1))
  ELSE
     CALL Get(D,TrixFile('D',Args,0))
  ENDIF
!
  CALL TrnMatBlk(BSp,GMp,D)
  CALL Get(BSiz,'atsiz',Tag_O=CurBase)
  CALL Get(OffS,'atoff',Tag_O=CurBase)
  CALL Get(NBasF,'nbasf',Tag_O=CurBase)
  CALL New(BfnInd,(/BSp%NAtms,BSp%NCtrt/))
!--------------------------------------------------------------------------------
! Compute and sort the distribution buffers
!--------------------------------------------------------------------------------
  CALL RangeOfDensity(D,NameBuf,BfnInd,DB,BSp,GMp)
  1000 DO WHILE (ErrorCode/=eAOK) 
    CALL MemInit(DB,IB,SB,Drv,BSc,BSp)
    CALL DisOrder(BSc,GMc,BSp,GMp,DB,IB,SB,Drv,NameBuf)
  END DO
!--------------------------------------------------------------------------------
! Allocate space for the exchange matrix. The routines below make sure 
! that there is *always* enough space allocated for the exchange matrix. 
! If this becomes too much to hold in memory then the variable ONXRange 
! should be set to some hard cut-off, and KThresh in CalcK will need to 
! be set to 1 instead of 0 (and possibly other things need to be changed...)
!--------------------------------------------------------------------------------
  CALL RangeOfExchange(BSc,GMc,BSp,GMp,D,NameBuf)
  CALL New(K,(/NRows+1,NCols,NElem/))
  CALL New(SubInd,(/3,NBasF/))
  CALL InitK(BSc,GMc,K,NameBuf)
  CALL InitSubInd(BSc,GMc,SubInd)
!--------------------------------------------------------------------------------
! All set to compute the exchange matrix
!--------------------------------------------------------------------------------
  CALL ComputeKg(BSc,GMc,BSp,GMp,D,K,DB,IB,SB,Drv,SubInd,BfnInd)
  IF (ErrorCode/=eAOK) THEN
    CALL Delete(K)
    CALL Delete(SubInd)
    GOTO 1000
  END IF
  CALL ComputeKe(BSc,GMc,BSp,GMp,D,K,DB,IB,SB,Drv,SubInd,BfnInd)
  IF (ErrorCode/=eAOK) THEN
    CALL Delete(K)
    CALL Delete(SubInd)
    GOTO 1000
  END IF
!--------------------------------------------------------------------------------
! Free up some space that we dont need anymore.
!--------------------------------------------------------------------------------
  CALL Delete(D)
  CALL Delete(DB)
  CALL Delete(IB)
  CALL Delete(SB)
  CALL Delete(Drv)
  CALL Delete(SubInd)
!--------------------------------------------------------------------------------
! Collect the distributed exchange matrices on the root node.
!--------------------------------------------------------------------------------
  CALL Fillout_BCSR(BSc,GMc,K)
  CALL TrnMatBlk(BSc,GMc,K)
  CALL ONXFilter(BSc,GMc,K,NameBuf,Thresholds%Trix)
  CALL Put(K,TrixFile('K',Args,0))!InFile,'.K')
  CALL Put(K%NBlks,'nki')
  CALL Put(K%NNon0,'nkm')

  CALL PPrint(K,'K')
  CALL PChkSum(K,'K',Prog)
  CALL PPrint(K,'K')
  CALL Plot(K,'K')
!--------------------------------------------------------------------------------
! Clean up...
!--------------------------------------------------------------------------------
  CALL Delete(K)
  CALL Delete(NameBuf)
  CALL Delete(BfnInd)
  CALL Delete(BSc)
  CALL Delete(GMc)
  CALL Delete(BSp)
  CALL Delete(GMp)
  CALL ShutDown(Prog)
END PROGRAM ONX

