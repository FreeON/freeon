PROGRAM ONX
!H=================================================================================
!H PROGRAM ONX
!H
!H  OPTIONS:
!H  DEBUGING: Use -DONX_DBUG to print some stuff. Cannot be used yet
!H  INFO    : Use -DONX_INFO to print some stuff.
!H
!H Comment:
!H
!H=================================================================================
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE AtomPairs  !per
  USE Clock
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE ONXParameters
  USE ONXMemory
  USE ONXCtrSclg, ONLY: TrnMatBlk   !old USE ContractionScaling
  USE ONXComptK , ONLY: ComputeK
  USE ONXDOrder , ONLY: DisOrder
  USE ONXMemInit, ONLY: MemInit
  USE ONXInit   , ONLY: InitSubInd,InitBfnInd,InitK
#ifdef PARALLEL
  USE ONXFillOut, ONLY: FillOutFASTMAT
  USE ONXRng    , ONLY: RangeOfExchangeFASTMAT
#else
  USE ONXFillOut, ONLY: FillOutBCSR
  USE ONXRng    , ONLY: RangeOfExchangeBCSR
#endif
#ifdef PARALLEL
  USE ONXGet    , ONLY: Get_Essential_RowCol
  USE PartDrv   , ONLY: PDrv_Initialize,PDrv_Finalize
  USE MondoMPI
  USE FastMatrices
#endif
  IMPLICIT NONE
#ifdef PARALLEL
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
  TYPE(DBuf)                     :: DBC        ! distribution buffers for indeces in C
  TYPE(DBuf)                     :: DBD        ! distribution buffers for indeces in D
  TYPE(IBuf)                     :: IB         ! 2-e eval buffers
  TYPE(DSL)                      :: SB         ! distribution pointers
  TYPE(ISpc)                     :: IS
  TYPE(IDrv)                     :: Drv        ! VRR/contraction drivers
  TYPE(INT_VECT)                 :: Stat
  TYPE(INT_RNK2)                 :: SubInd
  TYPE(INT_VECT)                 :: BfnInd
#ifdef PARALLEL
  TYPE(INT_VECT)                 :: RowPt,ColPt
  INTEGER                        :: NBrRow,NBrCol
#endif
  INTEGER                        :: ErrorCodeTmp
!--------------------------------------------------------------------------------
! Misc. variables and parameters...
!--------------------------------------------------------------------------------
  INTEGER                        :: iSwitch,OldFileID
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: InFile,RestartHDF
  CHARACTER(LEN=3),PARAMETER     :: Prog='ONX'
  TYPE(DBL_VECT)                 :: PBC
  TYPE(DBL_RNK3)                 :: SchT
  TYPE(INT_RNK3)                 :: BufT
  TYPE(INT_RNK2)                 :: BufN
  REAL(DOUBLE)                   :: Time1,Time2,xTotNERIs,TmK,TmDO,TmRE,TmFO,TmTM
#ifdef PARALLEL
  INTEGER                        :: IErr
  REAL(DOUBLE)                   :: TmBegK,TmEndK,TmKT,TmBegKT,TmEndKT
  TYPE(DBL_VECT)                 :: TmKArr,NERIsArr,TmDOArr,TmREArr,TmFOArr,TmTMArr,TmKTArr
#endif
  INTEGER                        :: I,NCC,NCD       !per
  CHARACTER(100)                 :: User
!--------------------------------------------------------------------------------
! vw comments:
! o clean up the code.
! o Optimise the fastmat in the innermost loop.
!--------------------------------------------------------------------------------
!
#ifdef PARALLEL
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
#else
  CALL StartUp(Args,Prog)
#endif
!
!CALL GetEnv('USER',User)
!IF(TRIM(User).EQ.'tymczak') WRITE(*,*)'    CJ, The ONX Ghost is back!'
!
  InFile=TRIM(ScrName)//'_Cyc'//TRIM(IntToChar(Args%i%i(1)))
  IF(SCFActn=='Restart')THEN
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
     CALL Delete(Stat) !NEW_29.07.03
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
#ifndef PARALLEL
     CALL Get(D,TrixFile('DeltaD',Args,0))
#else
     CALL PDrv_Initialize(DFastMat,TrixFile('DeltaD',Args,0),'ONXPart',Args)
#endif
  ELSEIF(SCFActn=='BasisSetSwitch')THEN
#ifndef PARALLEL
     CALL Get(D,TrixFile('D',Args,-1))
#else
     CALL PDrv_Initialize(DFastMat,TrixFile('D',Args,-1),'ONXPart',Args)
#endif
  ELSE
#ifndef PARALLEL
     CALL Get(D,TrixFile('D',Args,0))
#else
     CALL PDrv_Initialize(DFastMat,TrixFile('D',Args,0),'ONXPart',Args)
#endif
  ENDIF
  !
  ! Start total timing.
#ifdef PARALLEL
  TmBegKT = MPI_WTIME()
#endif
  !
#ifdef PARALLEL
  CALL TrnMatBlk(BSp,GMp,DFastMat)
#else
  CALL TrnMatBlk(BSp,GMp,D       )
#endif
  CALL Get(BSiz ,'atsiz',Tag_O=CurBase)
  CALL Get(OffS ,'atoff',Tag_O=CurBase)
  CALL Get(NBasF,'nbasf',Tag_O=CurBase)
  CALL New(BfnInd,NAtoms)
  !
  CALL New(PBC,3)
  CALL SetEq(PBC,Zero)
  !
!--------------------------------------------------------------------------------
! Compute and sort the distribution buffers
!--------------------------------------------------------------------------------
#ifdef PARALLEL
  time1 = MPI_WTIME()
  CALL Get_Essential_RowCol(DFastMat,RowPt,NbrRow,ColPt,NbrCol)
  CALL InitBfnInd(DBC,BSp,GMp,RowPt,NbrRow)
  CALL InitBfnInd(DBD,BSp,GMp,ColPt,NbrCol,BfnInd)
  time2 = MPI_WTIME()
  !write(*,*) 'List',time2-time1,MyID
#else
  CALL CPU_TIME(time1)
  CALL InitBfnInd(DBC,BSp,GMp,BfnInd)
  CALL InitBfnInd(DBD,BSp,GMp,BfnInd)
  CALL CPU_TIME(time2)
  !write(*,*) 'List',time2-time1
#endif
!
#ifdef PARALLEL
  time1 = MPI_WTIME()
  ! Have to bump both buffers 1 and 2, otherwise if buffer2 is to
  ! small, buffer 1 gets bumped infinately untill you run out of memory!
  ErrorCodeTmp = ErrorCode
  DO WHILE (ErrorCode/=eAOK)
     CALL MemInit(DBC,IB,SB,Drv,BSc,BSp,SchT,BufT,BufN)
     CALL DisOrder(PBC,BSc,GMc,BSp,GMp,DBC,IB,SB,Drv,SchT,BufT,BufN,RowPt,NBrRow)
  ENDDO
  ErrorCode = ErrorCodeTmp
  DO WHILE (ErrorCode/=eAOK)
     CALL MemInit(DBD,IB,SB,Drv,BSc,BSp,SchT,BufT,BufN)
     CALL DisOrder(PBC,BSc,GMc,BSp,GMp,DBD,IB,SB,Drv,SchT,BufT,BufN,ColPt,NBrCol)
  ENDDO
  time2 = MPI_WTIME()
  !write(*,*) 'DisOrder',time2-time1,MyID
#else
1000 CONTINUE
  CALL CPU_TIME(time1)
!    Have to bump both buffers 1 and 2, otherwise if buffer2 is to
!    small, buffer 1 gets bumped infinately untill you run out of memory!
  ErrorCodeTmp=ErrorCode
  DO WHILE (ErrorCode/=eAOK) 
     CALL MemInit(DBC,IB,SB,Drv,BSc,BSp,SchT,BufT,BufN)
     CALL DisOrder(PBC,BSc,GMc,BSp,GMp,DBC,IB,SB,Drv,SchT,BufT,BufN)
  ENDDO
!new if(CS_OUT%NCells.NE.1) THEN
  ErrorCode=ErrorCodeTmp
  DO WHILE (ErrorCode/=eAOK) 
     CALL MemInit(DBD,IB,SB,Drv,BSc,BSp,SchT,BufT,BufN)
     CALL DisOrder(PBC,BSc,GMc,BSp,GMp,DBD,IB,SB,Drv,SchT,BufT,BufN)
  ENDDO
!new endif
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
#ifdef PARALLEL
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
#ifdef PARALLEL
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
#ifndef PARALLEL
  IF(ErrorCode/=eAOK) THEN
     CALL Delete(K)
     CALL Delete(SubInd)
     GOTO 1000
  ENDIF
#endif
  !
!--------------------------------------------------------------------------------
! All set to compute the exchange matrix
!--------------------------------------------------------------------------------
#ifdef PARALLEL
  TmK  = Zero
  TmDO = Zero
#endif
  xTotNERIs = Zero                                                                            !per
  ! Periodic double sum over R and Rprime                                                     !per
  DO NCC = 1,CS_OUT%NCells                                                                    !per
     PBC%D(:) = CS_OUT%CellCarts%D(:,NCC)                                                     !per
#ifdef PARALLEL
     time1 = MPI_WTIME()
     CALL DisOrder(PBC,BSc,GMc,BSp,GMp,DBC,IB,SB,Drv,SchT,BufT,BufN,RowPt,NBrRow)
     time2 = MPI_WTIME()
#else
     CALL CPU_TIME(time1)
!new if(CS_OUT%NCells.NE.1) THEN
     CALL DisOrder(PBC,BSc,GMc,BSp,GMp,DBC,IB,SB,Drv,SchT,BufT,BufN)
!new endif
     CALL CPU_TIME(time2)
#endif
     TmDO = TmDO+time2-time1
     IF(DBC%LenTC.EQ.0) CYCLE
     DO NCD = 1,CS_OUT%NCells                                                                 !per
        PBC%D(:) = CS_OUT%CellCarts%D(:,NCD)                                                  !per
#ifdef PARALLEL
        time1 = MPI_WTIME()
        CALL DisOrder(PBC,BSc,GMc,BSp,GMp,DBD,IB,SB,Drv,SchT,BufT,BufN,ColPt,NBrCol)
        time2 = MPI_WTIME()
#else
        CALL CPU_TIME(time1)
!new if(CS_OUT%NCells.NE.1) THEN
        CALL DisOrder(PBC,BSc,GMc,BSp,GMp,DBD,IB,SB,Drv,SchT,BufT,BufN)
!new endif
        CALL CPU_TIME(time2)
#endif
        TmDO = TmDO+time2-time1
        IF(DBD%LenTC.EQ.0) CYCLE
!--------------------------------------------------------------------------------
!       All set to compute the exchange matrix
!--------------------------------------------------------------------------------
#ifdef PARALLEL
        time1 = MPI_WTIME()
        CALL ComputeK(BSc,GMc,BSp,GMp,DFastMat,KFastMat,DBC,DBD,IB,SB,IS,Drv,SubInd,BfnInd)
        time2 = MPI_WTIME()
        TmK = TmK+time2-time1
#else
!new if(CS_OUT%NCells.EQ.1) THEN
!       CALL ComputeK(BSc,GMc,BSp,GMp,D       ,K       ,DBC,DBC,IB,SB,IS,Drv,SubInd,BfnInd)
!new then
        CALL ComputeK(BSc,GMc,BSp,GMp,D       ,K       ,DBC,DBD,IB,SB,IS,Drv,SubInd,BfnInd)
!new endif
        IF(ErrorCode/=eAOK) THEN
           CALL Delete(K)
           CALL Delete(SubInd)
           GOTO 1000
        ENDIF
#endif
        xTotNERIs = xTotNERIs+xNERIs
     ENDDO
  ENDDO
  !
#ifndef PARALLEL
  IF(ErrorCode/=eAOK) THEN
     CALL Delete(K)
     CALL Delete(SubInd)
     GOTO 1000
  ENDIF
#endif
  !
!--------------------------------------------------------------------------------
! Redistribute partition informations.
!--------------------------------------------------------------------------------
#ifdef PARALLEL
  CALL PDrv_Finalize(DFastMat,CollectInPar_O=.TRUE.)
#endif
  !
!--------------------------------------------------------------------------------
! Free up some space that we dont need anymore.
!--------------------------------------------------------------------------------
#ifdef PARALLEL
  CALL Delete_FastMat1(DFastMat)
  CALL Delete(RowPt )
  CALL Delete(ColPt )
#else
  CALL Delete(D     )
#endif
  CALL Delete(DBC   )
  CALL Delete(DBD   )
  CALL Delete(PBC   )
  !
  CALL Delete(IB    )
  CALL Delete(SB    )
  CALL Delete(Drv   )
  CALL Delete(SubInd)
  !
!--------------------------------------------------------------------------------
! Collect the distributed exchange matrices on the root node.
!--------------------------------------------------------------------------------
  TmFO = Zero
#ifdef PARALLEL
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
#ifdef PARALLEL
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
#ifdef PARALLEL
  !
  ! End Total Timing
  TmEndKT = MPI_WTIME()
  TmKT = TmEndKT-TmBegKT
  !
  IF(MyID.EQ.ROOT) THEN
     CALL New(TmKArr  ,NPrc)
     CALL New(TmKTArr ,NPrc)
     CALL New(NERIsArr,NPrc)
     CALL New(TmDOArr ,NPrc)
     CALL New(TmREArr ,NPrc)
     CALL New(TmFOArr ,NPrc)
     CALL New(TmTMArr ,NPrc)
  ENDIF
  CALL MPI_Gather(TmK      ,1,MPI_DOUBLE_PRECISION,TmKArr%D(1)  ,1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmKT     ,1,MPI_DOUBLE_PRECISION,TmKTArr%D(1) ,1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(xTotNERIs,1,MPI_DOUBLE_PRECISION,NERIsArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmDO     ,1,MPI_DOUBLE_PRECISION,TmDOArr%D(1) ,1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmRE     ,1,MPI_DOUBLE_PRECISION,TmREArr%D(1) ,1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmFO     ,1,MPI_DOUBLE_PRECISION,TmFOArr%D(1) ,1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmTM     ,1,MPI_DOUBLE_PRECISION,TmTMArr%D(1) ,1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
#endif
  !
  !
#ifdef PARALLEL
  !
  IF(SCFActn == 'InkFok') CALL Halt('InkFok in PARALLEL ONX is not supported.')
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
#ifdef PARALLEL
  !
  IF(MyID.EQ.ROOT) THEN
     ! Imbalance stuff.
     CALL PImbalance(TmKArr ,NPrc,Prog_O='ComputeK')
     CALL PImbalance(TmKTArr,NPrc,Prog_O='ONX'     )
     !
     WRITE(*,1001) SUM(TmKArr%D  )/DBLE(NPrc),MINVAL(TmKArr%D  ),MAXVAL(TmKArr%D  )
     WRITE(*,1002) SUM(NERIsArr%D)/DBLE(NPrc),MINVAL(NERIsArr%D),MAXVAL(NERIsArr%D)
     WRITE(*,1003) SUM(TmDOArr%D )/DBLE(NPrc),MINVAL(TmDOArr%D ),MAXVAL(TmDOArr%D )
     WRITE(*,1004) SUM(TmREArr%D )/DBLE(NPrc),MINVAL(TmREArr%D ),MAXVAL(TmREArr%D )
     WRITE(*,1005) SUM(TmFOArr%D )/DBLE(NPrc),MINVAL(TmFOArr%D ),MAXVAL(TmFOArr%D )
     WRITE(*,1006) SUM(TmTMArr%D )/DBLE(NPrc),MINVAL(TmTMArr%D ),MAXVAL(TmTMArr%D )
     !
     CALL Delete(TmKArr  )
     CALL Delete(TmKTArr )
     CALL Delete(NERIsArr)
     CALL Delete(TmDOArr )
     CALL Delete(TmREArr )
     CALL Delete(TmFOArr )
     CALL Delete(TmTMArr )
     !
1001 FORMAT(' ONX: Ave TmK  = ',F15.2,', Min TmK  = ',F15.2,', Max TmK  = ',F15.2)
1002 FORMAT(' ONX: Tot ERI  = ',F15.2,', Min ERI  = ',F15.2,', Max ERI  = ',F15.2)
1003 FORMAT(' ONX: Ave TmDO = ',F15.2,', Min TmDO = ',F15.2,', Max TmDO = ',F15.2)
1004 FORMAT(' ONX: Ave TmRE = ',F15.2,', Min TmRE = ',F15.2,', Max TmRE = ',F15.2)
1005 FORMAT(' ONX: Ave TmFO = ',F15.2,', Min TmFO = ',F15.2,', Max TmFO = ',F15.2)
1006 FORMAT(' ONX: Ave TmTM = ',F15.2,', Min TmTM = ',F15.2,', Max TmTM = ',F15.2)
  ENDIF
  !
#else
  !
  WRITE(*,1001) TmK
  WRITE(*,1002) xTotNERIs
  WRITE(*,1003) TmDO
  WRITE(*,1004) TmRE
  WRITE(*,1005) TmFO
  WRITE(*,1006) TmTM
  !
1001 FORMAT(' ONX: Tot TmK  = ',F15.2)
1002 FORMAT(' ONX: Tot ERI  = ',F15.2)
1003 FORMAT(' ONX: Tot TmDO = ',F15.2)
1004 FORMAT(' ONX: Tot TmRE = ',F15.2)
1005 FORMAT(' ONX: Tot TmFO = ',F15.2)
1006 FORMAT(' ONX: Tot TmTM = ',F15.2)
  !
#endif !PARALLEL
#endif !ONX_INFO
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
  !
  CALL ShutDown(Prog)
  !
END PROGRAM ONX
