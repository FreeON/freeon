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
  TYPE(BUFL)          :: LenF
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
!  CALL PPrint(D,'D used in ONX')
  CALL TrnMatBlk(BSp,GMp,D)

  CALL Get(BSiz,'atsiz',Tag_O=CurBase)
  CALL Get(OffS,'atoff',Tag_O=CurBase)
  CALL Get(NBasF,'nbasf',Tag_O=CurBase)

!--------------------------------------------------------------------------------
! Compute and sort the distribution buffers, save them to disk.
!--------------------------------------------------------------------------------
  CALL RangeOfDensity(D,NameBuf)
  CALL ONXBuf(BSc,GMc,BSp,GMp,LenF,NameBuf)
!--------------------------------------------------------------------------------
! Allocate the exact space needed for the distribution buffers.
!--------------------------------------------------------------------------------
  CALL New(PrmBuf,LenF%PrmL)
  CALL New(DisBuf,LenF%DisL)
  CALL New(CBuf,LenF%CosL)
  CALL New(SBuf,LenF%SinL)
  CALL New(VecBuf,(/3,LenF%VecL/))
  CALL New(BfnBuf,LenF%BfnL)
  CALL New(iSA,(/5,LenF%iSAL/))
  CALL New(iSrt,(/BSc%NCtrt,NAtoms/))
!--------------------------------------------------------------------------------
! Read in the buffers, then erase the files from disk.
!--------------------------------------------------------------------------------
  CALL ONXIO(MyID,FRead,FilePrm,TypeR,LenF%PrmL,PrmBuf%D(1))
  CALL ONXIO(MyID,FRead,FileDis,TypeR,LenF%DisL,DisBuf%D(1))
  CALL ONXIO(MyID,FRead,FileVec,TypeR,3*LenF%VecL,VecBuf%D(1,1))
  CALL ONXIO(MyID,FRead,FileCos,TypeR,LenF%CosL,CBuf%D(1))
  CALL ONXIO(MyID,FRead,FileSin,TypeR,LenF%SinL,SBuf%D(1))
  CALL ONXIO(MyID,FRead,FileiSA,TypeI,5*LenF%iSAL,iSA%I(1,1))
  CALL ONXIO(MyID,FRead,FileBfn,TypeI,LenF%BfnL,BfnBuf%I(1))
  CALL ONXIO(MyID,FRead,FileSrt,TypeI,LenF%SrtL,iSrt%I(1,1))
  CALL ONXIO(MyID,FErase,FilePrm,0,0,0)
  CALL ONXIO(MyID,FErase,FileDis,0,0,0)
  CALL ONXIO(MyID,FErase,FileVec,0,0,0)
  CALL ONXIO(MyID,FErase,FileCos,0,0,0)
  CALL ONXIO(MyID,FErase,FileSin,0,0,0)
  CALL ONXIO(MyID,FErase,FileiSA,0,0,0)
  CALL ONXIO(MyID,FErase,FileBfn,0,0,0)
  CALL ONXIO(MyID,FErase,FileSrt,0,0,0)
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
  CALL InitK(BSc,GMc,K,NameBuf,SubInd)
!--------------------------------------------------------------------------------
! All set to compute the exchange matrix! 
!--------------------------------------------------------------------------------
  CALL ONXLoop(BSc,GMc,BSp,GMp,D,K,PrmBuf,DisBuf,CBuf,SBuf,VecBuf,BfnBuf, &
               iSA,iSrt,SubInd,LenF)
!--------------------------------------------------------------------------------
! Free up some space.
!--------------------------------------------------------------------------------
  CALL Delete(D)
  CALL Delete(PrmBuf)
  CALL Delete(DisBuf)
  CALL Delete(CBuf)
  CALL Delete(SBuf)
  CALL Delete(VecBuf)
  CALL Delete(BfnBuf)
  CALL Delete(iSA)
  CALL Delete(iSrt)
  CALL Delete(SubInd)
!--------------------------------------------------------------------------------
! Collect the distributed exchange matrices on the root node.
!--------------------------------------------------------------------------------
#ifdef PARALLEL
  CALL ONXFilter(BS,GM,K,NameBuf,Thresholds%Trix*1.0D-2)
  CALL New(KTotal,(/NAtoms,MaxBlks,MaxNon0/))
  CALL InitK(BS,GM,KTotal)
  CALL GatherK(BS,GM,NameBuf,KTotal,K)
  CALL Delete(K)
  CALL Fillout_BCSR(BS,GM,KTotal)
  CALL TrnMatBlk(BS,GM,KTotal)
  CALL ONXFilter(BS,GM,KTotal)
  CALL PPrint(KTotal,'KTotal after filter, from PONX')
  CALL Put(KTotal,TrixFile('K',Args,0)) !InFile,'.K')
  CALL Put(KTotal%NBlks,'nki')
  CALL Put(KTotal%NNon0,'nkm')
  IF(PrintFlags%Key>=DEBUG_MEDIUM) THEN
    CALL PChkSum(KTotal,'K',Prog)
  END IF
  IF(PrintFlags%Mat==DEBUG_MATRICES)THEN
     CALL PPrint(KTotal,'K')
  ELSEIF(PrintFlags%Mat==PLOT_MATRICES)THEN
     CALL Plot(KTotal,'K')
  ENDIF
  CALL Delete(KTotal)
#else
  CALL Fillout_BCSR(BSc,GMc,K)
  CALL TrnMatBlk(BSc,GMc,K)
  CALL ONXFilter(BSc,GMc,K,NameBuf,Thresholds%Trix)
  CALL Put(K,TrixFile('K',Args,0))!InFile,'.K')
  CALL Put(K%NBlks,'nki')
  CALL Put(K%NNon0,'nkm')

  CALL PPrint(K,'K')
  IF(PrintFlags%Key>=DEBUG_MEDIUM) THEN
    CALL PChkSum(K,'K',Prog)
  END IF
  IF(PrintFlags%Mat==DEBUG_MATRICES)THEN
     CALL PPrint(K,'K')
  ELSEIF(PrintFlags%Mat==PLOT_MATRICES)THEN
     CALL Plot(K,'K')
  ENDIF
#endif
!--------------------------------------------------------------------------------
! Clean up...
!--------------------------------------------------------------------------------
  CALL Delete(K)
  CALL Delete(NameBuf)
  CALL Delete(BSc)
  CALL Delete(GMc)
  CALL Delete(BSp)
  CALL Delete(GMp)
  CALL ShutDown(Prog)
END PROGRAM ONX

