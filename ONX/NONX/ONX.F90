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
  TYPE(DBuf)          :: DB
  TYPE(IBuf)          :: IB
  TYPE(INT_VECT)      :: NameBuf
  TYPE(INT_RNK2)      :: SubInd
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
  CHARACTER(LEN=5),PARAMETER     :: Prog='PONX'
!--------------------------------------------------------------------------------
  CALL StartUp(Args,Prog)
  InFile=TRIM(SCFName)//'_Cyc'//TRIM(IntToChar(Args%i%i(1)))
  CALL Get(BSc,Tag_O=CurBase)
  CALL Get(GMc,Tag_O=CurGeom)
  CALL Get(BSp,Tag_O=PrvBase)
  CALL Get(GMp,Tag_O=PrvGeom)
  CALL New(NameBuf,NAtoms)
  CALL Get(BSiz,'atsiz',Tag_O=PrvBase)
  CALL Get(OffS,'atoff',Tag_O=PrvBase)
  CALL Get(NBasF,'nbasf',Tag_O=PrvBase)

  CALL Get(D,TrixFile('D',Args,0)) !InFile,'.D')
  CALL TrnMatBlk(BSp,GMp,D)

  CALL Get(BSiz,'atsiz',Tag_O=CurBase)
  CALL Get(OffS,'atoff',Tag_O=CurBase)
  CALL Get(NBasF,'nbasf',Tag_O=CurBase)
!--------------------------------------------------------------------------------
! Compute and sort the distribution buffers
!--------------------------------------------------------------------------------
  write(*,*) "RangeOfDensity"
  CALL RangeOfDensity(D,NameBuf)

  DO WHILE (ErrorCode/=0) 
    CALL MemInit(DB,IB,BSc,BSp)
    write(*,*) "DisOrder"
    CALL DisOrder(BSc,GMc,BSp,GMp,DB,IB,NameBuf) 
  END DO
!--------------------------------------------------------------------------------
! Allocate space for the exchange matrix. The routines below make sure 
! that there is *always* enough space allocated for the exchange matrix. 
! If this becomes too much to hold in memory then the variable ONXRange 
! should be set to some hard cut-off, and KThresh in CalcK will need to 
! be set to 1 instead of 0 (and possibly other things need to be changed...)
!--------------------------------------------------------------------------------
!  CALL RangeOfExchange(BSc,GMc,BSp,GMp,D,NameBuf)
!  CALL New(K,(/NRows+1,NCols,NElem/))
!  CALL New(SubInd,(/3,NBasF/))
!  CALL InitK(BSc,GMc,K,NameBuf,SubInd)
!--------------------------------------------------------------------------------
! All set to loop over Ltypes and contraction types
!--------------------------------------------------------------------------------
!  CALL Looper(BSc,GMc,BSp,GMp,D,K,PrmBuf,DisBuf,CBuf,SBuf,VecBuf,BfnBuf, &
!               iSA,iSrt,SubInd)
!--------------------------------------------------------------------------------
! Free up some space.
!--------------------------------------------------------------------------------
!  CALL Delete(D)
!  CALL Delete(PrmBuf)
!  CALL Delete(DisBuf)
!  CALL Delete(CBuf)
!  CALL Delete(SBuf)
!  CALL Delete(VecBuf)
!  CALL Delete(BfnBuf)
!  CALL Delete(iSA)
!  CALL Delete(iSrt)
!  CALL Delete(SubInd)
!--------------------------------------------------------------------------------
! Collect the distributed exchange matrices on the root node.
!--------------------------------------------------------------------------------
!  CALL Fillout_BCSR(BSc,GMc,K)
!  CALL TrnMatBlk(BSc,GMc,K)
!  CALL ONXFilter(BSc,GMc,K,NameBuf,Thresholds%Trix)
!  CALL Put(K,TrixFile('K',Args,0))!InFile,'.K')
!  CALL Put(K%NBlks,'nki')
!  CALL Put(K%NNon0,'nkm')
!
!  CALL PPrint(K,'K')
!  IF(PrintFlags%Key>=DEBUG_MEDIUM) THEN
!    CALL PChkSum(K,'K',Prog)
!  END IF
!  IF(PrintFlags%Mat==DEBUG_MATRICES)THEN
!     CALL PPrint(K,'K')
!  ELSEIF(PrintFlags%Mat==PLOT_MATRICES)THEN
!     CALL Plot(K,'K')
!  ENDIF
!--------------------------------------------------------------------------------
! Clean up...
!--------------------------------------------------------------------------------
!  CALL Delete(K)
!  CALL Delete(NameBuf)
!  CALL Delete(BSc)
!  CALL Delete(GMc)
!  CALL Delete(BSp)
!  CALL Delete(GMp)
  CALL ShutDown(Prog)
END PROGRAM ONX

