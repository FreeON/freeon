PROGRAM ONX
!H=================================================================================
!H PROGRAM ONX
!H
!H  OPTIONS:
!H  DEBUGING: Use -.... to print some stuff.
!H  INFO    : Use -DONX_INFO to print some stuff.
!H
!H Comment:
!H
!H=================================================================================
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
#ifdef PERIODIC  !per
  USE AtomPairs  !per
#endif           !per
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
  USE ONXComptK , ONLY: ComputeK
  USE ONXDOrder , ONLY: DisOrder,PBC
  USE ONXMemInit, ONLY: MemInit
  USE ONXInit   , ONLY: InitSubInd,InitBfnInd,InitK
#ifdef PARALLEL_ONX
  USE ONXFillOut, ONLY: FillOutFASTMAT
  USE ONXRng    , ONLY: RangeOfExchangeFASTMAT
#else
  USE ONXFillOut, ONLY: FillOutBCSR
  USE ONXRng    , ONLY: RangeOfExchangeBCSR
#endif
!old #ifdef PARALLEL_ONX
#ifdef PARALLEL_ONX !per
  USE ONXGet    , ONLY: Get_Essential_RowCol
  USE PartDrv   , ONLY: PDrv_Initialize,PDrv_Finalize
  USE MondoMPI
  USE FastMatrices
#endif
  IMPLICIT NONE
#ifdef PARALLEL_ONX
  TYPE(FastMat),POINTER          :: KFastMat
  TYPE(FastMat),POINTER          :: DFastMat
  TYPE(DBCSR)                    :: D
  TYPE(BCSR )                    :: K
#else
  TYPE(BCSR)                     :: D
  TYPE(BCSR)                     :: K,T1,T2
#endif
  TYPE(BSET)                     :: BSc
  TYPE(CRDS)                     :: GMc
  TYPE(BSET)                     :: BSp
  TYPE(CRDS)                     :: GMp
  TYPE(ARGMT)                    :: Args
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!old #ifdef PARALLEL_ONX
  TYPE(DBuf)                     :: DBC        ! distribution buffers for indeces in C
  TYPE(DBuf)                     :: DBD        ! distribution buffers for indeces in D
!old#else
!old  TYPE(DBuf)                 :: DB         ! distribution buffers
!old#endif
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  TYPE(IBuf)                     :: IB         ! 2-e eval buffers
  TYPE(DSL)                      :: SB         ! distribution pointers
  TYPE(ISpc)                     :: IS
  TYPE(IDrv)                     :: Drv        ! VRR/contraction drivers
  TYPE(INT_VECT)                 :: Stat
  TYPE(INT_RNK2)                 :: SubInd
  TYPE(INT_VECT)                 :: BfnInd
#ifdef PARALLEL_ONX
  TYPE(INT_VECT)                 :: RowPt,ColPt
  INTEGER                        :: NBrRow,NBrCol
#endif
  INTEGER                        :: ErrorCodeTmp
!--------------------------------------------------------------------------------
! Large buffers to hold the sorted distribution and multipole data.
!--------------------------------------------------------------------------------
  TYPE(DBL_VECT)                 :: PrmBuf,DisBuf,CBuf,SBuf
  TYPE(DBL_RNK2)                 :: VecBuf
  TYPE(INT_VECT)                 :: BfnBuf
  TYPE(INT_RNK2)                 :: iSA,iSrt
!--------------------------------------------------------------------------------
! Misc. variables and parameters...
!--------------------------------------------------------------------------------
  INTEGER                        :: iSwitch,OldFileID
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: InFile,RestartHDF
  CHARACTER(LEN=3),PARAMETER     :: Prog='ONX'
  TYPE(DBL_RNK2)                 :: B      
  REAL(DOUBLE)                   :: Time1,Time2,tmpT,xTotNERIs,TmK,TmDO,TmRE,TmFO,TmTM
#ifdef PARALLEL_ONX
  INTEGER :: IErr
  REAL(DOUBLE)                   :: TmBegK,TmEndK,TmKT,TmBegKT,TmEndKT
  TYPE(DBL_VECT)                 :: TmKArr
#endif
!old#ifdef PERIODIC                                 !per
  INTEGER                        :: I,NCC,NCD       !per
  TYPE(CellSet)                  :: CSTemp          !per
!new REAL(DOUBLE), DIMENSION(3) :: PBC              !per
!old#endif                                          !per
!--------------------------------------------------------------------------------
! vw comments:
! o clean up the code.
! o Optimise the fastmat in the innermost loop.
!--------------------------------------------------------------------------------
!
#ifdef PARALLEL_ONX
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
#else
  CALL StartUp(Args,Prog)
#endif

  InFile=TRIM(ScrName)//'_Cyc'//TRIM(IntToChar(Args%i%i(1)))
  IF(SCFActn=='Restart')THEN
!toremove #ifdef PARALLEL_CLONES
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

     !CALL Get(GMp,Tag_O=PrvGeom) !old
     CALL Get(GMp,Tag_O=CurGeom) !per

     CALL Get(BSiz ,'atsiz',Tag_O=PrvBase)
     CALL Get(OffS ,'atoff',Tag_O=PrvBase)
     CALL Get(NBasF,'nbasf',Tag_O=PrvBase)
     !  CALL Get(CS_OUT,'CS_OUT',Tag_O=CurBase)
     ! Close the old hdf up 
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(OldFileID)
     ! Reopen current group and HDF
     HDFFileID=OpenHDF(H5File)
     H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
     HDF_CurrentID=H5GroupID
!toremove#else
!toremove!    Get the old information
!toremove     CALL Get(RestartHDF,'OldInfo')
!toremove     CALL CloseHDF(HDF_CurrentID)
!toremove     HDF_CurrentID=OpenHDF(RestartHDF)
!toremove     CALL New(Stat,3)
!toremove     CALL Get(Stat,'current')
!toremove     PrvCycl=TRIM(IntToChar(Stat%I(1)))
!toremove     PrvBase=TRIM(IntToChar(Stat%I(2)))
!toremove     PrvGeom=TRIM(IntToChar(Stat%I(3)))
!toremove     CALL Get(BSp,Tag_O=PrvBase)
!toremove!    Get the previous geometry, ASSUMING that 
!toremove!    we are not extrapolating the DM
!toremove     CALL Get(GMp,Tag_O=PrvGeom)
!toremove     CALL Get(BSiz ,'atsiz',Tag_O=PrvBase)
!toremove     CALL Get(OffS ,'atoff',Tag_O=PrvBase)
!toremove     CALL Get(NBasF,'nbasf',Tag_O=PrvBase)
!toremove     CALL CloseHDF(HDF_CurrentID)
!toremove     HDF_CurrentID=OpenHDF(InfFile)     
!toremove#endif
  ELSE
     CALL Get(BSp,Tag_O=PrvBase)
!    Get the current geometry here...
     CALL Get(GMp,Tag_O=CurGeom)
     CALL Get(BSiz ,'atsiz',Tag_O=PrvBase)
     CALL Get(OffS ,'atoff',Tag_O=PrvBase)
     CALL Get(NBasF,'nbasf',Tag_O=PrvBase)
  ENDIF
  CALL Get(BSc,Tag_O=CurBase)
  CALL Get(GMc,Tag_O=CurGeom)

!----------------------------------------------
! Get the Density Matrix
!----------------------------------------------
  IF(SCFActn=='InkFok')THEN
#ifndef PARALLEL_ONX
     CALL Get(D,TrixFile('DeltaD',Args,0))
#else
     CALL PDrv_Initialize(DFastMat,TrixFile('DeltaD',Args,0),'ONXPart',Args)
#endif
  ELSEIF(SCFActn=='BasisSetSwitch')THEN
#ifndef PARALLEL_ONX
     CALL Get(D,TrixFile('D',Args,-1))
#else
     CALL PDrv_Initialize(DFastMat,TrixFile('D',Args,-1),'ONXPart',Args)
#endif
  ELSE
#ifndef PARALLEL_ONX
     CALL Get(D,TrixFile('D',Args,0))
#else
     CALL PDrv_Initialize(DFastMat,TrixFile('D',Args,0),'ONXPart',Args)
#endif
  ENDIF
!
! Start total timing.
#ifdef PARALLEL_ONX
  TmBegKT = MPI_WTIME()
#endif

#ifdef PARALLEL_ONX
  CALL TrnMatBlk(BSp,GMp,DFastMat)
#else
  CALL TrnMatBlk(BSp,GMp,D       )
#endif
  CALL Get(BSiz ,'atsiz',Tag_O=CurBase)
  CALL Get(OffS ,'atoff',Tag_O=CurBase)
  CALL Get(NBasF,'nbasf',Tag_O=CurBase)
  CALL New(BfnInd,NAtoms)
  !
!--------------------------------------------------------------------------------
! Compute and sort the distribution buffers
!--------------------------------------------------------------------------------
  PBC(:) = Zero !per
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef PARALLEL_ONX
  time1 = MPI_WTIME()
  CALL Get_Essential_RowCol(DFastMat,RowPt,NbrRow,ColPt,NbrCol)
  CALL InitBfnInd(DBC,BSp,GMp,RowPt,NbrRow)
  CALL InitBfnInd(DBD,BSp,GMp,ColPt,NbrCol,BfnInd)
  time2 = MPI_WTIME()
  !write(*,*) 'List',time2-time1,MyID
#else
  CALL CPU_TIME(time1)
!old  CALL InitBfnInd(DB,BSp,GMp,BfnInd)
  CALL InitBfnInd(DBC,BSp,GMp,BfnInd) !per
  CALL InitBfnInd(DBD,BSp,GMp,BfnInd) !per
  CALL CPU_TIME(time2)
  !write(*,*) 'List',time2-time1
#endif
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
#ifdef PARALLEL_ONX
  time1 = MPI_WTIME()
!1000 CONTINUE
  ! Have to bump both buffers 1 and 2, otherwise if buffer2 is to
  ! small, buffer 1 gets bumped infinately untill you run out of memory!
  ErrorCodeTmp = ErrorCode
  DO WHILE (ErrorCode/=eAOK)
     CALL MemInit(DBC,IB,SB,Drv,BSc,BSp)
     CALL DisOrder(BSc,GMc,BSp,GMp,DBC,IB,SB,Drv,RowPt,NBrRow) !To add PBC,BufN,BufT,SchT
  ENDDO
  ErrorCode = ErrorCodeTmp
  DO WHILE (ErrorCode/=eAOK)
     CALL MemInit(DBD,IB,SB,Drv,BSc,BSp)
     CALL DisOrder(BSc,GMc,BSp,GMp,DBD,IB,SB,Drv,ColPt,NBrCol) !To add PBC,BufN,BufT,SchT
  ENDDO
  time2 = MPI_WTIME()
  !write(*,*) 'DisOrder',time2-time1,MyID
#else
!old  CALL CPU_TIME(time1)
!old1000 DO WHILE (ErrorCode/=eAOK)
!old     CALL MemInit(DB,IB,SB,Drv,BSc,BSp)
!old     CALL DisOrder(BSc,GMc,BSp,GMp,DB,IB,SB,Drv)
!old  ENDDO
!old  CALL CPU_TIME(time2)
!old  ! write(*,*) 'DisOrder',time2-time1
1000 CONTINUE
  CALL CPU_TIME(time1)
!    Have to bump both buffers 1 and 2, otherwise if buffer2 is to
!    small, buffer 1 gets bumped infinately untill you run out of memory!
  ErrorCodeTmp=ErrorCode
  DO WHILE (ErrorCode/=eAOK) 
     CALL MemInit(DBC,IB,SB,Drv,BSc,BSp)
     CALL DisOrder(BSc,GMc,BSp,GMp,DBC,IB,SB,Drv)         !To add PBC,BufN,BufT,SchT
  ENDDO
  ErrorCode=ErrorCodeTmp
  DO WHILE (ErrorCode/=eAOK) 
     CALL MemInit(DBD,IB,SB,Drv,BSc,BSp)
     CALL DisOrder(BSc,GMc,BSp,GMp,DBD,IB,SB,Drv)         !To add PBC,BufN,BufT,SchT
  ENDDO
  CALL CPU_TIME(time2)
  ! write(*,*) 'DisOrder',time2-time1
#endif
  !
!--------------------------------------------------------------------------------
! Allocate space for the exchange matrix. The routines below make sure 
! that there is *always* enough space allocated for the exchange matrix. 
! If this becomes too much to hold in memory then the variable ONXRange 
! should be set to some hard cut-off, and KThresh in CalcK will need to 
! be set to 1 instead of 0 (and possibly other things need to be changed...)
!--------------------------------------------------------------------------------
  TmRE = Zero
#ifdef PARALLEL_ONX
  time1 = MPI_WTIME()
  CALL RangeOfExchangeFASTMAT(BSc,GMc,BSp,GMp,DFastMat)
  time2 = MPI_WTIME()
  !write(*,*) 'Range of K',time2-time1,MyID
#else
  CALL CPU_TIME(time1)
  CALL RangeOfExchangeBCSR(BSc,GMc,BSp,GMp,D)
  CALL CPU_TIME(time2)
  !write(*,*) 'Range of K',time2-time1
#endif
  TmRE = time2-time1
  !
!------------------------------------------------------------------------------- 
!     Initialize the matrix and associated indecies
!------------------------------------------------------------------------------- 
#ifdef PARALLEL_ONX
  CALL New_FASTMAT(KFastMat,0,(/0,0/))
#else
  CALL New(K,(/NRows+1,NCols,NElem/))
  CALL SetEq(K%MTrix,Zero)
  CALL InitK(BSc,GMc,K)
#endif
  !
  CALL New(SubInd,(/3,NBasF/))
  CALL InitSubInd(BSc,GMc,SubInd)
  !
#ifndef PARALLEL_ONX
  IF(ErrorCode/=eAOK) THEN
     !CALL Delete_FASTMAT(KFastMat)
     !#else
     CALL Delete(K)
     CALL Delete(SubInd)
     GOTO 1000
  ENDIF
#endif
  !
!--------------------------------------------------------------------------------
! All set to compute the exchange matrix
!--------------------------------------------------------------------------------
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#ifdef PARALLEL_ONX
  TmK  = Zero
  TmDO = Zero
#endif
  xTotNERIs = Zero                                                                            !per
  !
  ! Set up third sum                                                                          !per
  CALL New_CellSet_Cube(CSTemp,GMc%PBC%AutoW,GMc%PBC%BoxShape,(/1,1,1/))                      !per
  !
  ! Periodic double sum over R and Rprime                                                     !per
  DO NCC = 1,CS_OUT%NCells                                                                    !per
     PBC(:) = CS_OUT%CellCarts%D(:,NCC)                                                       !per
#ifdef PARALLEL_ONX
     time1 = MPI_WTIME()
     CALL DisOrder(BSc,GMc,BSp,GMp,DBC,IB,SB,Drv,RowPt,NBrRow)    !To add PBC,BufN,BufT,SchT  !per
     time2 = MPI_WTIME()
#else
     CALL CPU_TIME(time1)
     CALL DisOrder(BSc,GMc,BSp,GMp,DBC,IB,SB,Drv)                 !To add PBC,BufN,BufT,SchT  !per
     CALL CPU_TIME(time2)
#endif
     TmDO = TmDO+time2-time1
     DO NCD = 1,CS_OUT%NCells                                                                 !per
        PBC(:) = CS_OUT%CellCarts%D(:,NCD)                                                    !per
#ifdef PARALLEL_ONX
        time1 = MPI_WTIME()
        CALL DisOrder(BSc,GMc,BSp,GMp,DBD,IB,SB,Drv,ColPt,NBrCol) !To add PBC,BufN,BufT,SchT  !per
        time2 = MPI_WTIME()
#else
        CALL CPU_TIME(time1)
        CALL DisOrder(BSc,GMc,BSp,GMp,DBD,IB,SB,Drv)              !To add PBC,BufN,BufT,SchT  !per
        CALL CPU_TIME(time2)
#endif
        TmDO = TmDO+time2-time1
!--------------------------------------------------------------------------------             !per
!       All set to compute the exchange matrix                                                !per
!--------------------------------------------------------------------------------             !per
#ifdef PARALLEL_ONX
        time1 = MPI_WTIME()                                                                   !per
        CALL ComputeK(BSc,GMc,BSp,GMp,DFastMat,KFastMat,DBC,DBD,IB,SB,IS,Drv,SubInd,BfnInd)   !per
        time2 = MPI_WTIME()                                                                   !per
        TmK = TmK+time2-time1                                                                 !per
#else                                                                                         !per
        CALL ComputeK(BSc,GMc,BSp,GMp,D       ,K       ,DBC,DBD,IB,SB,IS,Drv,SubInd,BfnInd)   !per
        IF(ErrorCode/=eAOK) THEN                                                              !per
           CALL Delete(K)                                                                     !per
           CALL Delete(SubInd)                                                                !per
           GOTO 1000                                                                          !per
        ENDIF                                                                                 !per
#endif                                                                                        !per
        xTotNERIs = xTotNERIs+xNERIs                                                          !per
     ENDDO                                                                                    !per
  ENDDO                                                                                       !per
  !
  ! Collect number of integrals.                                                              !per
!#ifdef PARALLEL_ONX                                                                           !per
!  xTotNERIs = Reduce(xTotNERIs)                                                               !per
!  TmDO = Reduce(TmDO)
!  TmRE = Reduce(TmRE)
!#endif                                                                                        !per

!old#ifdef PARALLEL_ONX
!old  time1 = MPI_WTIME()
!old  TmBegK = MPI_WTIME()
!old  CALL ComputeK(BSc,GMc,BSp,GMp,DFastMat,KFastMat,DBC,DBD,IB,SB,IS,Drv,SubInd,BfnInd)
!old  TmEndK = MPI_WTIME()
!old  TmK = TmEndK-TmBegK
!old  time2 = MPI_WTIME()
!old  write(*,*) 'K Build',time2-time1,MyID
!old  ! Collect number of integrals.
!old  xTotNERIs = Zero
!old  xTotNERIs = Reduce(xNERIs)
!old#else
!old  CALL CPU_TIME(time1)
!old  CALL ComputeK(BSc,GMc,BSp,GMp,D,K,DB,IB,SB,IS,Drv,SubInd,BfnInd)
!old  CALL CPU_TIME(time2)
!old  write(*,*) 'K Build',time2-time1
!old#endif

#ifndef PARALLEL_ONX
  IF(ErrorCode/=eAOK) THEN
     CALL Delete(K)
     CALL Delete(SubInd)
     GOTO 1000
  ENDIF
#endif
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !
!--------------------------------------------------------------------------------
! Redistribute partition informations.
!--------------------------------------------------------------------------------
#ifdef PARALLEL_ONX
  CALL PDrv_Finalize(DFastMat,CollectInPar_O=.TRUE.)
#endif
  !
!--------------------------------------------------------------------------------
! Free up some space that we dont need anymore.
!--------------------------------------------------------------------------------
!old#ifdef PARALLEL_ONX
!old  CALL Delete_FastMat1(DFastMat)
!old  CALL Delete(RowPt )
!old  CALL Delete(ColPt )
!old  CALL Delete(DBC   )
!old  CALL Delete(DBD   )
!old#else
!old  CALL Delete(D     )
!old  CALL Delete(DB    )
!old#endif

#ifdef PARALLEL_ONX                   !per
  CALL Delete_FastMat1(DFastMat)      !per
  CALL Delete(RowPt )                 !per
  CALL Delete(ColPt )                 !per
#else                                 !per
  CALL Delete(D     )                 !per
#endif                                !per
  CALL Delete(DBC   )                 !per
  CALL Delete(DBD   )                 !per

  CALL Delete(IB    )
  CALL Delete(SB    )
  CALL Delete(Drv   )
  CALL Delete(SubInd)
  !
!--------------------------------------------------------------------------------
! Collect the distributed exchange matrices on the root node.
!--------------------------------------------------------------------------------
  TmFO = Zero
#ifdef PARALLEL_ONX
  time1 = MPI_WTIME()
  CALL FillOutFastMat(BSc,GMc,KFastMat)
  time2 = MPI_WTIME()
  !write(*,*) 'Symmetrized',time2-time1,MyID
#else
  CALL CPU_TIME(time1)
  CALL FillOutBCSR(BSc,GMc,K)
  CALL CPU_TIME(time2)
  !write(*,*) 'Symmetrized',time2-time1
#endif
  TmFO = time2-time1
  !
  !
  TmTM = Zero
#ifdef PARALLEL_ONX
  time1 = MPI_WTIME()
  CALL TrnMatBlk(BSc,GMc,KFastMat)
  time2 = MPI_WTIME()
  !write(*,*) 'Normalization2',time2-time1,MyID
#else
  CALL CPU_TIME(time1)
  CALL TrnMatBlk(BSc,GMc,K)
  CALL CPU_TIME(time2)
  !write(*,*) 'Normalization2',time2-time1
#endif  
  TmTM = time2-time1
  !
  ! Collect number of integrals.
#ifdef PARALLEL_ONX
  xTotNERIs = Reduce(xTotNERIs)
  TmDO = Reduce(TmDO)
  TmRE = Reduce(TmRE)
  TmFO = Reduce(TmFO)
  TmTM = Reduce(TmTM)
#endif
  !
  !
#ifdef PARALLEL_ONX
  !
  ! End Total Timing
  TmEndKT = MPI_WTIME()
  TmKT = TmEndKT-TmBegKT

  IF(SCFActn == 'InkFok') STOP 'InkFok in PARALLEL ONX is not supported.'
  !time1 = MPI_WTIME()
  ! Collect the data on the root.
  CALL Redistribute_FASTMAT(KFastMat)
  !time2 = MPI_WTIME()
  !write(*,*) 'Redistribute',MyID,time2-time1
  !time1 = MPI_WTIME()
  CALL Set_BCSR_EQ_DFASTMAT(K,KFastMat)
  !time2 = MPI_WTIME()
  !write(*,*) 'Set_BCSR_Eq_FASTMAT',MyID,time2-time1
  CALL Delete_FastMat1(KFastMat)
  !
#ifdef ONX_INFO
  ! Imbalance stuff.
  CALL New(TmKArr,NPrc)
  CALL MPI_Gather(TmK,1,MPI_DOUBLE_PRECISION,TmKArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  IF(MyID == ROOT) THEN
     CALL PImbalance(TmKArr,NPrc,Prog_O='ComputeK')
  ENDIF
  CALL Delete(TmKArr)
  !
  ! Imbalance stuff.
  CALL New(TmKArr,NPrc)
  CALL MPI_Gather(TmKT,1,MPI_DOUBLE_PRECISION,TmKArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  IF(MyID == ROOT) THEN
     CALL PImbalance(TmKArr,NPrc,Prog_O='ONX')
  ENDIF
  CALL Delete(TmKArr)
#endif
  !
#else
  ! Add in correction if incremental K build
  IF(SCFActn == 'InkFok')THEN
     CALL New(T1)
     CALL New(T2)
     CALL Get(T1,TrixFile('K',Args,-1))
     CALL Add(K,T1,T2)
     CALL Filter(K,T2)
     CALL Delete(T1)
     CALL Delete(T2)
  ENDIF
#endif
  !
  !
#ifdef ONX_INFO
#ifdef PARALLEL_ONX
  IF(MyID.EQ.ROOT) WRITE(*,*) 'NbrTot ERI = ',xTotNERIs
  IF(MyID.EQ.ROOT) WRITE(*,*) 'TmDO = ',TmDO
#else
  WRITE(*,*) 'ONX: NbrTot ERI = ',xTotNERIs
  WRITE(*,*) 'ONX: TmDO = ',TmDO
  WRITE(*,*) 'ONX: TmRE = ',TmRE
  WRITE(*,*) 'ONX: TmFO = ',TmFO
  WRITE(*,*) 'ONX: TmTM = ',TmTM
#endif
#endif
  !
  ! Save on disc.
  CALL Put(K,TrixFile('K',Args,0))
  CALL PChkSum(K,'Kx['//TRIM(SCFCycl)//']',Prog)
  CALL PPrint( K,'Kx['//TRIM(SCFCycl)//']')
  CALL Plot(   K,'Kx['//TRIM(SCFCycl)//']')
!--------------------------------------------------------------------------------
! Clean up...
!--------------------------------------------------------------------------------
  CALL Delete(K     )
  CALL Delete(BfnInd)
  CALL Delete(BSc   )
  CALL Delete(GMc   )
  CALL Delete(BSp   )
  CALL Delete(GMp   )
!  IF(myid==0) write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  CALL ShutDown(Prog)
  !
END PROGRAM ONX
























