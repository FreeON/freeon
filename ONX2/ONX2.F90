PROGRAM ONX2
  !
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  !
  !Old
  USE ONXParameters
  USE ONXInit   , ONLY: InitK
  USE ONXCtrSclg, ONLY: TrnMatBlk
#ifdef PARALLEL
  USE ONXRng    , ONLY: RangeOfExchangeFASTMAT
  USE ONXFillOut, ONLY: FillOutFASTMAT
#else
  USE ONXRng    , ONLY: RangeOfExchangeBCSR
  USE ONXFillOut, ONLY: FillOutBCSR
#endif

#ifdef PARALLEL
  !USE ONXGet    , ONLY: Get_Essential_RowCol
  USE PartDrv   , ONLY: PDrv_Initialize,PDrv_Finalize
  USE MondoMPI
  USE FastMatrices
#endif
  !
  !New
  USE ONX2DataType
  USE ONX2List
  USE ONX2ComputK
  !
  IMPLICIT NONE
  !
#ifdef PARALLEL
  TYPE(FASTMAT),POINTER          :: KxFastMat
  TYPE(FASTMAT),POINTER          :: DFastMat
  TYPE(DBCSR)                    :: D
  TYPE(BCSR )                    :: Kx
#else
  TYPE(BCSR)                     :: D
  TYPE(BCSR)                     :: Kx,T1,T2
#endif
  TYPE(BSET)                     :: BSc
  TYPE(CRDS)                     :: GMc
  TYPE(BSET)                     :: BSp
  TYPE(CRDS)                     :: GMp
  TYPE(ARGMT)                    :: Args
  TYPE(INT_VECT)                 :: Stat
!--------------------------------------------------------------------------------
! Misc. variables and parameters...
!--------------------------------------------------------------------------------
  INTEGER                        :: OldFileID
  REAL(DOUBLE)                   :: time1,time2
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: InFile,RestartHDF
  CHARACTER(LEN=*),PARAMETER     :: Prog='ONX2'
!--------------------------------------------------------------------------------
  INTEGER :: IErr
  REAL(DOUBLE)   :: TmML,TmRE,TmKx,TmFO,TmTM
  TYPE(DBL_VECT) :: TmKxArr,TmMLArr,TmREArr,TmFOArr,TmTMArr!,TmKTArr,NERIsArr
!New
  TYPE(CList2), DIMENSION(:), POINTER :: ListC,ListD

!
#ifdef PARALLEL
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
#else
  CALL StartUp(Args,Prog)
#endif
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
     ! Get the previous geometry, ASSUMING that 
     ! we are not extrapolating the DM
     CALL Get(GMp,Tag_O=CurGeom)
     CALL Get(BSiz ,'atsiz',Tag_O=PrvBase)
     CALL Get(OffS ,'atoff',Tag_O=PrvBase)
     CALL Get(NBasF,'nbasf',Tag_O=PrvBase)
     ! Close the old hdf up 
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(OldFileID)
     ! Reopen current group and HDF
     HDFFileID=OpenHDF(H5File)
     H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
     HDF_CurrentID=H5GroupID
     CALL Delete(Stat)
  ELSE
     CALL Get(BSp,Tag_O=PrvBase)
     ! Get the current geometry here...
     CALL Get(GMp,Tag_O=CurGeom)
     CALL Get(BSiz ,'atsiz',Tag_O=PrvBase)
     CALL Get(OffS ,'atoff',Tag_O=PrvBase)
     CALL Get(NBasF,'nbasf',Tag_O=PrvBase)
  ENDIF
  CALL Get(BSc,Tag_O=CurBase)
  CALL Get(GMc,Tag_O=CurGeom)
  !
  !
!!$  SELECT CASE(SCFActn)
!!$  !IF(SCFActn=='StartResponse'.OR.SCFActn=='FockPrimeBuild')THEN
!!$  CASE('StartResponse','FockPrimeBuild')
!!$#ifndef PARALLEL
!!$     CALL Get(D,TrixFile('DPrime'//Args%C%C(4),Args,0))
!!$#else
!!$     CALL PDrv_Initialize(DFastMat,TrixFile('DPrime'//Args%C%C(4),Args,0),'ONXPart',Args)
!!$#endif
!!$  !ELSEIF(SCFActn=='InkFok')THEN
!!$  CASE('InkFok')
!!$#ifndef PARALLEL
!!$     CALL Get(D,TrixFile('DeltaD',Args,0))
!!$#else
!!$     CALL PDrv_Initialize(DFastMat,TrixFile('DeltaD',Args,0),'ONXPart',Args)
!!$#endif
!!$  !ELSEIF(SCFActn=='BasisSetSwitch')THEN
!!$  CASE('BasisSetSwitch')
!!$#ifndef PARALLEL
!!$     CALL Get(D,TrixFile('D',Args,-1))
!!$#else
!!$     CALL PDrv_Initialize(DFastMat,TrixFile('D',Args,-1),'ONXPart',Args)
!!$#endif
!!$  !ELSE
!!$  CASE DEFAULT
#ifdef PARALLEL
     CALL PDrv_Initialize(DFastMat,TrixFile('D',Args,0),'ONXPart',Args)
#else
     CALL Get(D,TrixFile('D',Args,0))
#endif
!!$  !ENDIF
!!$  END SELECT
  !
  !
#ifdef PARALLEL
  CALL TrnMatBlk(BSp,GMp,DFastMat)
#else
  CALL TrnMatBlk(BSp,GMp,D       )
#endif
  !
  !
#ifdef PARALLEL
  IF(MyID.EQ.ROOT) &
#endif
  WRITE(*,*) 'WE ARE IN ONX2 WE ARE IN ONX2 WE ARE IN ONX2 WE ARE IN ONX2 WE ARE IN ONX2'
  !
  !
  !WRITE(*,*) 'allocate List'
#ifdef PARALLEL
  Time1 = MPI_WTIME()
  CALL AllocList(ListC,NAtoms) !!!!!!!!!!!!!!!!!! Add atom1-atom2 or sthg like that !!!!!!!!!!!!!!!!!!
  CALL AllocList(ListD,NAtoms) !!!!!!!!!!!!!!!!!! Add atom1-atom2 or sthg like that !!!!!!!!!!!!!!!!!!
  Time2 = MPI_WTIME()
#else
  CALL CPU_TIME(Time1)
  CALL AllocList(ListC,NAtoms)
  CALL GetBufferSize(GMc,BSc)
  CALL CPU_TIME(Time2)
#endif
  !WRITE(*,*) 'allocate List: ok',Time2-Time1
  !
  !------------------------------------------------
  ! Make the distribution list(s).
  !WRITE(*,*) 'make List'
#ifdef PARALLEL
  Time1 = MPI_WTIME()
  CALL MakeList(ListC,GMc,BSc,CS_OUT) !!!!!!!!!!!!!!!!!! Add atom list !!!!!!!!!!!!!!!!!!
  CALL MakeList(ListD,GMc,BSc,CS_OUT) !!!!!!!!!!!!!!!!!! Add atom list !!!!!!!!!!!!!!!!!!
  Time2 = MPI_WTIME()
#else
  CALL CPU_TIME(Time1)
  CALL MakeList(ListC,GMc,BSc,CS_OUT)
  CALL CPU_TIME(Time2)
#endif
  TmML = Time2-Time1
  !WRITE(*,*) 'make List: ok',Time2-Time1
  !
  !------------------------------------------------
  !
#ifdef ONX2_DBUG
  WRITE(*,*) 'Print List'
#ifdef PARALLEL
  CALL PrintList(ListC)
  CALL PrintList(ListD)
#else
  CALL PrintList(ListC)
#endif
  WRITE(*,*) 'Print List:ok'
#endif
  !
  !
  !------------------------------------------------
  !
#ifdef PARALLEL
  time1 = MPI_WTIME()
  CALL RangeOfExchangeFASTMAT(BSc,GMc,BSp,GMp,DFastMat)
  time2 = MPI_WTIME()
#else
  CALL CPU_TIME(time1)
  CALL RangeOfExchangeBCSR(BSc,GMc,BSp,GMp,D)
  CALL CPU_TIME(time2)
#endif
  TmRE = time2-time1
  !
  !
  !------------------------------------------------
  !
#ifdef PARALLEL
  CALL New_FASTMAT(KxFastMat,0,(/0,0/))
#else
  !write(*,*) 'ONXRange',ONXRange
  CALL New(Kx,(/NRows+1,NCols,NElem/))
  CALL SetEq(Kx%MTrix,Zero)
  CALL InitK(BSc,GMc,Kx)
#endif
  !
  !------------------------------------------------
  ! Compute Exchange matrix.
  !WRITE(*,*) 'Compute Kx'
#ifdef PARALLEL
  Time1 = MPI_WTIME()
  CALL ComputK(DFastMat,KxFastMat,ListC,ListD,GMc,BSc,CS_OUT)
  Time2 = MPI_WTIME()
#else
  CALL CPU_TIME(Time1)
  CALL ComputK(D,Kx,ListC,ListC,GMc,BSc,CS_OUT)
  CALL CPU_TIME(Time2)
#endif
  TmKx = Time2-Time1
  !WRITE(*,*) 'Compute Kx:ok',Time2-Time1
  !
  !------------------------------------------------
  ! Free up some space. Deallocate the list(s).
  !WRITE(*,*) 'deallocate List'
#ifdef PARALLEL
  Time1 = MPI_WTIME()
  CALL DeAllocList(ListC)
  CALL DeAllocList(ListD)
  Time2 = MPI_WTIME()
#else
  CALL CPU_TIME(Time1)
  CALL DeAllocList(ListC)
  CALL CPU_TIME(Time2)
#endif
  !WRITE(*,*) 'deallocate List:ok',Time2-Time1
  !
  !------------------------------------------------
  ! 
#ifdef PARALLEL
  time1 = MPI_WTIME()
  CALL FillOutFastMat(BSc,GMc,KxFastMat)
  time2 = MPI_WTIME()
#else
  CALL CPU_TIME(time1)
  CALL FillOutBCSR(BSc,GMc,Kx)
  CALL CPU_TIME(time2)
#endif
  TmFO = time2-time1
  !
  !------------------------------------------------
  ! Normilization.
#ifdef PARALLEL
  time1 = MPI_WTIME()
  CALL TrnMatBlk(BSc,GMc,KxFastMat)
  time2 = MPI_WTIME()
#else
  CALL CPU_TIME(time1)
  CALL TrnMatBlk(BSc,GMc,Kx       )
  CALL CPU_TIME(time2)
#endif
  TmTM = time2-time1
  !
  !
  !------------------------------------------------
  ! Redistribute partition informations.
#ifdef PARALLEL
  CALL PDrv_Finalize(DFastMat,CollectInPar_O=.TRUE.)
  CALL Delete_FastMat1(DFastMat)
#else
  CALL Delete(D)
#endif
  !
  !
  !------------------------------------------------
  !
#ifdef PARALLEL
  !
  ! End Total Timing
  !TmEndKT = MPI_WTIME()
  !TmKT = TmEndKT-TmBegKT
  !
  CALL New(TmKxArr ,NPrc)
  CALL New(TmMLArr ,NPrc)
  CALL New(TmREArr ,NPrc)
  CALL New(TmFOArr ,NPrc)
  CALL New(TmTMArr ,NPrc)
  !CALL New(TmKTArr ,NPrc)
  !CALL New(NERIsArr,NPrc)
  !
  CALL MPI_Gather(TmKx,1,MPI_DOUBLE_PRECISION,TmKxArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmML,1,MPI_DOUBLE_PRECISION,TmMLArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmRE,1,MPI_DOUBLE_PRECISION,TmREArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmFO,1,MPI_DOUBLE_PRECISION,TmFOArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmTM,1,MPI_DOUBLE_PRECISION,TmTMArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  !CALL MPI_Gather(TmKT     ,1,MPI_DOUBLE_PRECISION,TmKTArr%D(1) ,1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  !CALL MPI_Gather(xTotNERIs,1,MPI_DOUBLE_PRECISION,NERIsArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
#endif
  !
  !------------------------------------------------
  !
#ifdef PARALLEL
  !
  IF(SCFActn == 'InkFok') CALL Halt('InkFok in PARALLEL ONX is not supported.')
  ! Collect the data on the root.
  CALL Redistribute_FASTMAT(KxFastMat)
  CALL Set_BCSR_EQ_DFASTMAT(Kx,KxFastMat)
  CALL Delete_FastMat1(KxFastMat)
  !
#else
  ! Add in correction if incremental K build
  IF(SCFActn == 'InkFok')THEN
     CALL New(T1)
     CALL New(T2)
     CALL Get(T1,TrixFile('K',Args,-1))
     CALL Add(Kx,T1,T2)
     CALL Filter(Kx,T2)
     CALL Delete(T1)
     CALL Delete(T2)
  ENDIF
#endif
  !
  !------------------------------------------------
  ! Save on disc.
  CALL Put(Kx,TrixFile('K',Args,0))
  CALL PChkSum(Kx,'Kx['//TRIM(SCFCycl)//']',Prog)
  CALL PPrint( Kx,'Kx['//TRIM(SCFCycl)//']')
  CALL Plot(   Kx,'Kx['//TRIM(SCFCycl)//']')
  !



#ifdef PARALLEL
  !
  IF(MyID.EQ.ROOT) THEN
     ! Imbalance stuff.
     CALL PImbalance(TmKxArr ,NPrc,Prog_O='ComputeK')
     !CALL PImbalance(TmKTArr,NPrc,Prog_O='ONX'     )
     !
     WRITE(*,1001) SUM(TmKxArr%D)/DBLE(NPrc),MINVAL(TmKxArr%D  ),MAXVAL(TmKxArr%D  )
     !WRITE(*,1002) SUM(NERIsArr%D)           ,MINVAL(NERIsArr%D),MAXVAL(NERIsArr%D)
     WRITE(*,1003) SUM(TmMLArr%D )/DBLE(NPrc),MINVAL(TmMLArr%D ),MAXVAL(TmMLArr%D )
     WRITE(*,1004) SUM(TmREArr%D )/DBLE(NPrc),MINVAL(TmREArr%D ),MAXVAL(TmREArr%D )
     WRITE(*,1005) SUM(TmFOArr%D )/DBLE(NPrc),MINVAL(TmFOArr%D ),MAXVAL(TmFOArr%D )
     WRITE(*,1006) SUM(TmTMArr%D )/DBLE(NPrc),MINVAL(TmTMArr%D ),MAXVAL(TmTMArr%D )
     !
1001 FORMAT(' ONX: Ave TmKx = ',F15.2,', Min TmKx = ',F15.2,', Max TmKx = ',F15.2)
!1002 FORMAT(' ONX: Tot ERI  = ',F15.2,', Min ERI  = ',F15.2,', Max ERI  = ',F15.2)
1003 FORMAT(' ONX: Ave TmML = ',F15.2,', Min TmML = ',F15.2,', Max TmML = ',F15.2)
1004 FORMAT(' ONX: Ave TmRE = ',F15.2,', Min TmRE = ',F15.2,', Max TmRE = ',F15.2)
1005 FORMAT(' ONX: Ave TmFO = ',F15.2,', Min TmFO = ',F15.2,', Max TmFO = ',F15.2)
1006 FORMAT(' ONX: Ave TmTM = ',F15.2,', Min TmTM = ',F15.2,', Max TmTM = ',F15.2)
  ENDIF
  !
  CALL Delete(TmKxArr  )
!  CALL Delete(TmKTArr )
!  CALL Delete(NERIsArr)
  CALL Delete(TmMLArr )
  CALL Delete(TmREArr )
  CALL Delete(TmFOArr )
  CALL Delete(TmTMArr )
  !
#else
  !
  !WRITE(*,1001) TmKx
!  WRITE(*,1002) xTotNERIs
  !WRITE(*,1003) TmML
  !WRITE(*,1004) TmRE
  !WRITE(*,1005) TmFO
  !WRITE(*,1006) TmTM
  !
1001 FORMAT(' ONX: Tot TmK  = ',F15.2)
!1002 FORMAT(' ONX: Tot ERI  = ',F15.2)
1003 FORMAT(' ONX: Tot TmDO = ',F15.2)
1004 FORMAT(' ONX: Tot TmRE = ',F15.2)
1005 FORMAT(' ONX: Tot TmFO = ',F15.2)
1006 FORMAT(' ONX: Tot TmTM = ',F15.2)
  !
#endif !PARALLEL

!--------------------------------------------------------------------------------
! Clean up...
!--------------------------------------------------------------------------------
  !
  CALL Delete(Kx)
  !
  CALL Delete(BSc   )
  CALL Delete(GMc   )
  CALL Delete(BSp   )
  CALL Delete(GMp   )
  !
  CALL ShutDown(Prog)
  !
END PROGRAM ONX2


