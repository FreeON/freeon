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
  USE DOrder
  USE ONXParameters
  USE ONXMemory
  USE ContractionScaling
  USE MatFilter
  USE InitExchangeMatrix
#ifdef PARALLEL_ONX
  USE MondoMPI
#endif
  IMPLICIT NONE
#ifdef PARALLEL_ONX
  TYPE(DBCSR)         :: D
  TYPE(DBCSR)         :: K,T1,T2
  TYPE(BCSR)          :: KTotal
#else
  TYPE(BCSR)          :: D
  TYPE(BCSR)          :: K,T1,T2
!  INTEGER             :: MyID=0
#endif
  TYPE(BSET)          :: BSc
  TYPE(CRDS)          :: GMc
  TYPE(BSET)          :: BSp
  TYPE(CRDS)          :: GMp
  TYPE(ARGMT)         :: Args
#ifdef PERIODIC
  TYPE(DBuf)          :: DB1,DB2   ! distribution buffers
#else
  TYPE(DBuf)          :: DB1       ! distribution buffers
#endif
  TYPE(IBuf)          :: IB        ! 2-e eval buffers
  TYPE(DSL)           :: SB        ! distribution pointers
  TYPE(ISpc)          :: IS
  TYPE(IDrv)          :: Drv       ! VRR/contraction drivers
  TYPE(INT_VECT)      :: NameBuf,Stat
  TYPE(INT_RNK2)      :: SubInd
  TYPE(INT_RNK2)      :: BfnInd
!--------------------------------------------------------------------------------
! Large buffers to hold the sorted distribution and multipole data.f
!--------------------------------------------------------------------------------
  TYPE(DBL_VECT)         :: PrmBuf,DisBuf,CBuf,SBuf
  TYPE(DBL_RNK2)         :: VecBuf
  TYPE(INT_VECT)         :: BfnBuf
  TYPE(INT_RNK2)         :: iSA,iSrt
!--------------------------------------------------------------------------------
! Misc. variables and parameters...
!--------------------------------------------------------------------------------
  INTEGER                :: iSwitch,ErrorCodeTmp,OldFileID
#ifdef PERIODIC
  INTEGER                :: NC1,NC2
#endif
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: InFile,RestartHDF
  CHARACTER(LEN=3),PARAMETER     :: Prog='ONX'
!--------------------------------------------------------------------------------
  CALL StartUp(Args,Prog)
#ifdef PARALLEL_CLONES
#else
  CALL Get(CS_OUT,'CS_OUT',Tag_O=CurBase)
#endif
  InFile=TRIM(ScrName)//'_Cyc'//TRIM(IntToChar(Args%i%i(1)))
  IF(SCFActn=='Restart')THEN
#ifdef PARALLEL_CLONES
     ! Close current group and HDF
     CALL CloseHDFGroup(H5GroupID)
     CALL CloseHDF(HDFFileID)
     ! Open old group and HDF
     HDF_CurrentID=OpenHDF(Restart)
     OldFileID=HDF_CurrentID
     CALL New(Stat,3)
     CALL Get(Stat,'current_state')
     PrvCycl=TRIM(IntToChar(Stat%I(1)))
     PrvBase=TRIM(IntToChar(Stat%I(2)))
     PrvGeom=TRIM(IntToChar(Stat%I(3)))
     HDF_CurrentID=OpenHDFGroup(HDF_CurrentID,"Clone #"//TRIM(IntToChar(MyClone)))
     CALL Get(BSp,Tag_O=PrvBase)
!    Get the previous geometry, ASSUMING that 
!    we are not extrapolating the DM
     CALL Get(GMp,Tag_O=PrvGeom)
     CALL Get(BSiz,'atsiz',Tag_O=PrvBase)
     CALL Get(OffS,'atoff',Tag_O=PrvBase)
     CALL Get(NBasF,'nbasf',Tag_O=PrvBase)
     ! Close the old hdf up 
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(OldFileID)
     ! Reopen current group and HDF
     HDFFileID=OpenHDF(H5File)
     H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
     HDF_CurrentID=H5GroupID
#else
!    Get the old information
     CALL Get(RestartHDF,'OldInfo')
     CALL CloseHDF(HDF_CurrentID)
     HDF_CurrentID=OpenHDF(RestartHDF)
     CALL New(Stat,3)
     CALL Get(Stat,'current')
     PrvCycl=TRIM(IntToChar(Stat%I(1)))
     PrvBase=TRIM(IntToChar(Stat%I(2)))
     PrvGeom=TRIM(IntToChar(Stat%I(3)))
     CALL Get(BSp,Tag_O=PrvBase)
!    Get the previous geometry, ASSUMING that 
!    we are not extrapolating the DM
     CALL Get(GMp,Tag_O=PrvGeom)
     CALL Get(BSiz,'atsiz',Tag_O=PrvBase)
     CALL Get(OffS,'atoff',Tag_O=PrvBase)
     CALL Get(NBasF,'nbasf',Tag_O=PrvBase)
     CALL CloseHDF(HDF_CurrentID)
     HDF_CurrentID=OpenHDF(InfFile)     
#endif
  ELSE
     CALL Get(BSp,Tag_O=PrvBase)
     CALL Get(GMp,Tag_O=PrvGeom)
     CALL Get(BSiz,'atsiz',Tag_O=PrvBase)
     CALL Get(OffS,'atoff',Tag_O=PrvBase)
     CALL Get(NBasF,'nbasf',Tag_O=PrvBase)
  ENDIF
  CALL Get(BSc,Tag_O=CurBase)
  CALL Get(GMc,Tag_O=CurGeom)
  CALL New(NameBuf,NAtoms)
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
  CALL RangeOfDensity(D,NameBuf,BfnInd,DB1,BSp,GMp)
  CALL RangeOfDensity(D,NameBuf,BfnInd,DB2,BSp,GMp)

     PBC_h=Zero
     PBC_g=Zero
1000 CONTINUE 
!    Have to bump both buffers 1 and 2, otherwise if buffer2 is to
!    small, buffer 1 gets bumped infinately untill you run out of memory!
     ErrorCodeTmp=ErrorCode
     DO WHILE (ErrorCode/=eAOK) 
        CALL MemInit(DB1,IB,SB,Drv,BSc,BSp)
        CALL DisOrder(BSc,GMc,BSp,GMp,DB1,IB,SB,Drv,NameBuf)        
     ENDDO
     ErrorCode=ErrorCodeTmp
     DO WHILE (ErrorCode/=eAOK) 
        CALL MemInit(DB2,IB,SB,Drv,BSc,BSp)
        CALL DisOrder(BSc,GMc,BSp,GMp,DB2,IB,SB,Drv,NameBuf)        
     ENDDO

!     ErrorCode=eInit
!     CALL MemInit(DB2,IB,SB,Drv,BSc,BSp)
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
!
#ifdef PERIODIC
! Periodic double sum over R and Rprime
  DO NC1=1,CS_OUT%NCells
     PBC_h=CS_OUT%CellCarts%D(:,NC1)
     PBC_g=Zero
     CALL DisOrder(BSc,GMc,BSp,GMp,DB1,IB,SB,Drv,NameBuf)
     DO NC2=1,CS_OUT%NCells
        PBC_g=CS_OUT%CellCarts%D(:,NC2)
        CALL DisOrder(BSc,GMc,BSp,GMp,DB2,IB,SB,Drv,NameBuf)
!--------------------------------------------------------------------------------
!       All set to compute the exchange matrix
!--------------------------------------------------------------------------------
        CALL ComputeKg(BSc,GMc,BSp,GMp,D,K,DB1,DB2,IB,SB,IS,Drv,SubInd,BfnInd)
        IF(ErrorCode/=eAOK) THEN
           CALL Delete(K)
           CALL Delete(SubInd)
           GOTO 1000
        END IF
        CALL ComputeKe(BSc,GMc,BSp,GMp,D,K,DB1,DB2,IB,SB,IS,Drv,SubInd,BfnInd)
#else
        CALL ComputeKg(BSc,GMc,BSp,GMp,D,K,DB1,DB1,IB,SB,IS,Drv,SubInd,BfnInd)
        IF(ErrorCode/=eAOK) THEN
           CALL Delete(K)
           CALL Delete(SubInd)
           GOTO 1000
        END IF
        CALL ComputeKe(BSc,GMc,BSp,GMp,D,K,DB1,DB1,IB,SB,IS,Drv,SubInd,BfnInd)
#endif
        IF(ErrorCode/=eAOK) THEN
          CALL Delete(K)
          CALL Delete(SubInd)
          GOTO 1000
        ENDIF
#ifdef PERIODIC
     ENDDO
  ENDDO
  CALL Delete(DB2)
#endif
!--------------------------------------------------------------------------------
! Free up some space that we dont need anymore.
!--------------------------------------------------------------------------------
  CALL Delete(D)
  CALL Delete(DB1)
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
  K%NAtms=NAtoms 
! Add in correction if incremental K build
  IF(SCFActn=='InkFok')THEN
     CALL New(T1)
     CALL New(T2)
     CALL Get(T1,TrixFile('K',Args,-1))
     CALL Add(K,T1,T2)
     CALL Filter(K,T2)
     CALL Delete(T1)
     CALL Delete(T2)
  ENDIF
  CALL Put(K,TrixFile('K',Args,0))
  CALL PChkSum(K,'Kx['//TRIM(SCFCycl)//']',Prog)
  CALL PPrint( K,'Kx['//TRIM(SCFCycl)//']')
!  CALL PPrint( K,'Kx['//TRIM(SCFCycl)//']',Unit_O=6)
!  CALL Plot(   K,'Kx['//TRIM(SCFCycl)//']')
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

