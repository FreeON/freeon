!------------------------------------------------------------------------------
!    This code is part of the FreeON suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the FreeON group for review, and possible
!    dissemination in future releases.
!------------------------------------------------------------------------------

#include "MondoConfig.h"

PROGRAM ONX2

#ifndef PARALLEL
#undef ONX2_PARALLEL
#endif

  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg

  USE ONX2DataType
  USE ONX2List  , ONLY: AllocList,DeAllocList,MakeList,PrintList
  USE ONX2ComputK
  USE ONXParameters
  USE ONXInit   , ONLY: InitK
  USE ONXCtrSclg, ONLY: TrnMatBlk
#ifdef ONX2_PARALLEL
  USE MondoMPI
  USE FastMatrices
  USE ONXRng    , ONLY: RangeOfExchangeFASTMAT
  USE ONXFillOut, ONLY: FillOutBCSR!FillOutFASTMAT,
  USE ONXGet    , ONLY: Get_Essential_RowCol,GetOffArr,Reduce_FASTMAT
  USE PartDrv   , ONLY: PDrv_Initialize,PDrv_Finalize
#else
  USE ONXRng    , ONLY: RangeOfExchangeBCSR
  USE ONXFillOut, ONLY: FillOutBCSR
  USE ONXGet    , ONLY: GetOffArr
#endif

  IMPLICIT NONE

#ifdef ONX2_PARALLEL
  TYPE(FASTMAT),POINTER          :: KxFM
  TYPE(FASTMAT),POINTER          :: DFM
  TYPE(DBCSR)                    :: D
  TYPE(BCSR)                     :: Kx,T1,T2
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

#ifdef ONX2_PARALLEL
  TYPE(DBL_RNK2)                 :: GradTmp
  TYPE(INT_VECT)                 :: APt,BPt,CPt,DPt
#endif

#ifdef ONX2_PARALLEL
  TYPE(DBL_VECT)                 :: TmKxArr,TmMLArr,TmTMArr,TmALArr,TmDLArr,TmREArr,TmFOArr
!  REAL(DOUBLE),EXTERNAL          :: MondoTimer
  INTEGER                        :: CMin,CMax,DMin,DMax,iErr
  INTEGER                        :: ANbr,BNbr,CNbr,DNbr
  TYPE(INT_VECT)                 :: OPart,GPart
  INTEGER                        :: iONXPartExist,iGONXPartExist
#endif
  REAL(DOUBLE)                   :: TmAl,TmDl,TmKx,TmML,TmRE,TmFO,TmTM
  INTEGER                        :: OldFileID
  REAL(DOUBLE)                   :: Time1,Time2
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: InFile
  CHARACTER(LEN=*),PARAMETER     :: Prog='ONX2'

  TYPE(INT_RNK2) :: OffArrC,OffArrP
#ifdef ONX2_PARALLEL
  TYPE(CList), DIMENSION(:), POINTER :: ListC,ListD
#else
  TYPE(CList), DIMENSION(:), POINTER :: ListC
#endif

#ifdef ONX2_PARALLEL
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
#else
  CALL StartUp(Args,Prog)
#endif

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
     ! Check if a partition exist in Restart
#ifdef ONX2_PARALLEL
     CALL New(OPart,NPrc*4,0)
     CALL New(GPart,NPrc*4,0)
     iONXPartExist=-1000;iGONXPartExist=-1000;
     CALL Get(iONXPartExist ,'ONXPartExist' )
     IF(iONXPartExist .EQ.NPrc) CALL Get(OPart,'ONXPart')
     CALL Get(iGONXPartExist,'GONXPartExist')
     IF(iGONXPartExist.EQ.NPrc) CALL Get(GPart,'GONXPart')
     !write(*,*) 'iONXPartExist' ,iONXPartExist,MyID
     !write(*,*) 'iGONXPartExist',iGONXPartExist,MyID
#endif
     ! Close the old hdf up 
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(OldFileID)
     ! Reopen current group and HDF
     HDFFileID=OpenHDF(H5File)
     H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
     HDF_CurrentID=H5GroupID
     ! Put the partition in the new hdf if needed
#ifdef ONX2_PARALLEL
     IF(iONXPartExist.EQ.NPrc) THEN
        CALL Put(iONXPartExist ,'ONXPartExist' )
        CALL Put(OPart,'ONXPart')
        IF(MyID.EQ.0)WRITE(*,*) 'ONX2: Successfully load/save the ONX Restart partition!'
     ENDIF
     IF(iGONXPartExist.EQ.NPrc) THEN
        CALL Put(iGONXPartExist,'GONXPartExist')
        CALL Put(GPart,'GONXPart')
        IF(MyID.EQ.0)WRITE(*,*) 'ONX2: Successfully load/save the GONX Restart partition!'
     ENDIF
     CALL Delete(OPart)
     CALL Delete(GPart)
#endif
     CALL Delete(Stat)
  ELSE
     CALL Get(BSp,Tag_O=PrvBase)
     ! Get the current geometry here...
     CALL Get(GMp,Tag_O=CurGeom)
     CALL Get(BSiz ,'atsiz',Tag_O=CurBase)
     CALL Get(OffS ,'atoff',Tag_O=CurBase)
     CALL Get(NBasF,'nbasf',Tag_O=CurBase)
  ENDIF
  IF(SCFActn=='TD-SCF')THEN   
     NoSym=.TRUE.
  ENDIF
  !
  CALL Get(BSc,Tag_O=CurBase)
  CALL Get(GMc,Tag_O=CurGeom)
  !
  !------------------------------------------------
  ! Initialization and allocations.
#ifdef ONX2_PARALLEL
  NULLIFY(ListC,ListD,DFM,KxFM)
#else
  NULLIFY(ListC)
#endif
  !
  TmTM=Zero;TmML=Zero;TmAL=Zero;TmDL=Zero
  !
  CALL New(OffArrC,(/BSc%NCtrt,BSc%NKind/))
  CALL New(OffArrP,(/BSp%NCtrt,BSp%NKind/))
  !
  CALL GetOffArr(OffArrC,BSc)
  CALL GetOffArr(OffArrP,BSp)
  !
  CALL GetBufferSize(GMc,BSc,GMp,BSp)
  !
  !------------------------------------------------
  ! Get denstiy matrix.
  !
  SELECT CASE(SCFActn)
  CASE('StartResponse','FockPrimeBuild')
#ifndef PARALLEL
     CALL Get(D,TrixFile('DPrime'//Args%C%C(3),Args,0))
#else
     CALL PDrv_Initialize(DFM,TrixFile('DPrime'//Args%C%C(3),Args,0),'ONXPart',Args)
#endif
  CASE('TD-SCF')
#ifndef PARALLEL
     CALL Get(D,TrixFile(Args%C%C(3),Args,0))
#else
     CALL PDrv_Initialize(DFM,TrixFile(Args%C%C(3),Args,0),'ONXPart',Args)
#endif
  CASE('InkFok')
#ifndef PARALLEL
     CALL Get(D,TrixFile('DeltaD',Args,0))
#else
     CALL PDrv_Initialize(DFM,TrixFile('DeltaD',Args,0),'ONXPart',Args)
#endif
  CASE('BasisSetSwitch')
     ! We need to load the previous basis set setting.
     CALL Get(BSiz ,'atsiz',Tag_O=PrvBase)
     CALL Get(OffS ,'atoff',Tag_O=PrvBase)
     CALL Get(NBasF,'nbasf',Tag_O=PrvBase)
#ifndef PARALLEL
     CALL Get(D,TrixFile('D',Args,-1))
#else
     CALL PDrv_Initialize(DFM,TrixFile('D',Args,-1),'ONXPart',Args)
#endif
     ! We need to load the current basis set setting.
     CALL Get(BSiz ,'atsiz',Tag_O=CurBase)
     CALL Get(OffS ,'atoff',Tag_O=CurBase)
     CALL Get(NBasF,'nbasf',Tag_O=CurBase)
  CASE DEFAULT
#ifdef ONX2_PARALLEL
     CALL PDrv_Initialize(DFM,TrixFile('D',Args,0),'ONXPart',Args)
#else
     CALL Get(D,TrixFile('D',Args,0))
#endif
  END SELECT
  !
  !------------------------------------------------
  ! Normalize the density matrix
#ifdef ONX2_PARALLEL
  Time1 = MondoTimer()
  CALL TrnMatBlk(BSp,GMp,DFM)
  Time2 = MondoTimer()
#else
  CALL CPU_TIME(Time1)
  CALL TrnMatBlk(BSp,GMp,D)
  CALL CPU_TIME(Time2)
#endif
  TmTM = Time2-Time1
  !
  !------------------------------------------------
  ! Allocate the list(s) and get the buffer sizes.
  !WRITE(*,*) 'allocate List'
#ifdef ONX2_PARALLEL
  Time1 = MondoTimer()
  CALL Get_Essential_RowCol(DFM,CPt,CNbr,CMin,CMax,DPt,DNbr,DMin,DMax)
  !write(*,*) 'CMin',CMin,'CMax',CMax,'MyID',MyID
  !write(*,*) 'DMin',DMin,'DMax',DMax,'MyID',MyID
  !
  CALL AllocList(ListC,CMin,CMax)
  CALL AllocList(ListD,DMin,DMax)
  Time2 = MondoTimer()
#else
  CALL CPU_TIME(Time1)
  CALL AllocList(ListC,1,NAtoms)
  CALL CPU_TIME(Time2)
#endif
  TmAL = Time2-Time1
  !WRITE(*,*) 'allocate List: ok',Time2-Time1
  !
  !------------------------------------------------
  ! Make the distribution list(s).
  !WRITE(*,*) 'make List'
#ifdef ONX2_PARALLEL
  Time1 = MondoTimer()
  CALL MakeList(ListC,GMc,BSc,GMp,BSp,CS_OUT,CPt,CNbr,APt,ANbr)
  CALL MakeList(ListD,GMc,BSc,GMp,BSp,CS_OUT,DPt,DNbr,BPt,BNbr)
  Time2 = MondoTimer()
#else
  CALL CPU_TIME(Time1)
  CALL MakeList(ListC,GMc,BSc,GMp,BSp,CS_OUT)
  CALL CPU_TIME(Time2)
#endif
  TmML = Time2-Time1
  !WRITE(*,*) 'make List: ok',Time2-Time1
  !
  !------------------------------------------------
  ! Print list.
#ifdef ONX2_DBUG
  WRITE(*,*) 'Print List'
#ifdef ONX2_PARALLEL
  CALL PrintList(ListC)
  CALL PrintList(ListD)
#else
  CALL PrintList(ListC)
#endif
  WRITE(*,*) 'Print List:ok'
#endif
  !
  !------------------------------------------------
  ! Get range of K.
#ifdef ONX2_PARALLEL
  time1 = MondoTimer()
  CALL RangeOfExchangeFASTMAT(BSc,GMc,BSp,GMp,DFM)
  time2 = MondoTimer()
#else
  CALL CPU_TIME(time1)
  CALL RangeOfExchangeBCSR(BSc,GMc,BSp,GMp,D)
  CALL CPU_TIME(time2)
#endif
  TmRE = time2-time1
  !
  !------------------------------------------------
  ! Allocate the K matrix.
#ifdef ONX2_PARALLEL
  CALL New_FASTMAT(KxFM,0,(/0,0/),NSMat_O=DFM%NSMat)
#else
  CALL New(Kx,(/NRows+1,NCols,NElem/),NSMat_O=D%NSMat)
  CALL DBL_VECT_EQ_DBL_SCLR(NElem,Kx%MTrix%D(1),0.0D0)
  CALL InitK(BSc,GMc,Kx)
#endif
  !
  !------------------------------------------------
  ! Compute Exchange matrix.
  !WRITE(*,*) 'Compute Kx'
#ifdef ONX2_PARALLEL
  Time1 = MondoTimer()
  CALL ComputK(DFM,KxFM,ListC,ListD,OffArrC,OffArrP,GMc,BSc,GMp,BSp,CS_OUT)
  Time2 = MondoTimer()
#else
  CALL CPU_TIME(Time1)
  CALL ComputK(D,Kx,ListC,ListC,OffArrC,OffArrP,GMc,BSc,GMp,BSp,CS_OUT)
  CALL CPU_TIME(Time2)
#endif
  TmKx = Time2-Time1
  !WRITE(*,*) 'Compute Kx:ok',Time2-Time1
  !
  !------------------------------------------------
  ! Free up some space. Deallocate the list(s).
  !WRITE(*,*) 'deallocate List'
#ifdef ONX2_PARALLEL
  Time1 = MondoTimer()
  CALL DeAllocList(ListC)
  CALL DeAllocList(ListD)
  Time2 = MondoTimer()
#else
  CALL CPU_TIME(Time1)
  CALL DeAllocList(ListC)
  CALL CPU_TIME(Time2)
#endif
  !WRITE(*,*) 'deallocate List:ok',Time2-Time1
  !
  !------------------------------------------------
  ! 
#ifdef ONX2_PARALLEL
  !otime1 = MondoTimer()
  !oCALL FillOutFastMat(BSc,GMc,KxFM)
  !otime2 = MondoTimer()
#else
  IF(.NOT.NoSym)THEN
     CALL CPU_TIME(time1)
     CALL FillOutBCSR(BSc,GMc,Kx)
     CALL CPU_TIME(time2)
  ENDIF
#endif
  TmFO = time2-time1
  !
  !------------------------------------------------
  ! Normilization.
#ifdef ONX2_PARALLEL
  !otime1 = MondoTimer()
  !oCALL TrnMatBlk(BSc,GMc,KxFM)
  !otime2 = MondoTimer()
#else
  CALL CPU_TIME(time1)
  CALL TrnMatBlk(BSc,GMc,Kx)
  CALL CPU_TIME(time2)
#endif
  TmTM = time2-time1
  !
  !------------------------------------------------
  ! Redistribute partition informations.
#ifdef ONX2_PARALLEL
  CALL PDrv_Finalize(DFM,CollectInPar_O=.TRUE.)
  CALL Delete_FastMat1(DFM)
#else
  CALL Delete(D)
#endif
  !
  !------------------------------------------------
  !
#ifdef ONX2_PARALLEL
  !
  IF(SCFActn == 'InkFok') CALL Halt('InkFok in PARALLEL ONX is not supported.')
  ! Collect the data on the root.
  !=======
  CALL Reduce_FASTMAT(T1,KxFM)
  !=======
  !!!CALL Redistribute_FASTMAT(KxFM)
  !!!CALL Set_BCSR_EQ_DFASTMAT(T1,KxFM)
  !=======
  CALL Delete_FastMat1(KxFM)
  !
  time1 = MondoTimer()
  IF(MyID.EQ.ROOT) THEN
     ! The following needs the -Fac- variable 
     ! in ComputK to take into account for the 
     ! double counting of diag(K).
     CALL TrnMatBlk(BSc,GMc,T1)
     CALL XPose(T1,T2)
     CALL Add(T1,T2,Kx)
     CALL Delete(T1)
     CALL Delete(T2)
  ENDIF
  time2 = MondoTimer()
  TmTM = time2-time1
  !
#else
  ! Add in correction if incremental K build
  IF(SCFActn == 'InkFok')THEN
     CALL New(T1,NSMat_O=Kx%NSMat)
     CALL New(T2,NSMat_O=Kx%NSMat)
     CALL Get(T1,TrixFile('K',Args,-1))
     CALL Add(Kx,T1,T2)
     CALL Filter(Kx,T2)
     CALL Delete(T1)
     CALL Delete(T2)
  ELSE
     CALL New(T1,NSMat_O=Kx%NSMat)
     CALL Filter(T1,Kx)
!!     WRITE(*,*)SIZE(Kx%MTrix%D),Kx%NNon0
!!     WRITE(*,*)SIZE(T1%MTrix%D),T1%NNon0
     CALL SetEq(Kx,T1)
!!     WRITE(*,*)SIZE(Kx%MTrix%D),Kx%NNon0
  ENDIF
#endif
  !
  !------------------------------------------------
  ! Save on disc.
  IF(SCFActn=='StartResponse'.OR.SCFActn=='FockPrimeBuild')THEN
     CALL Put(Kx,TrixFile('KPrime'//TRIM(Args%C%C(3)),Args,0))
     CALL PChkSum(Kx,'Kx'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']',Prog)
     CALL PPrint( Kx,'Kx'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']')
     CALL Plot(   Kx,'Kx'//TRIM(Args%C%C(3))//'['//TRIM(SCFCycl)//']')
  ELSE
     CALL Put(Kx,TrixFile('K',Args,0))
     CALL PChkSum(Kx,'Kx['//TRIM(SCFCycl)//']',Prog)
     CALL PPrint( Kx,'Kx['//TRIM(SCFCycl)//']')
     CALL Plot(   Kx,'Kx['//TRIM(SCFCycl)//']')
  ENDIF
  !
  !------------------------------------------------
  ! Timing.
  !
#ifdef GONX2_INFO
#ifdef ONX2_PARALLEL
  !
  ! End Total Timing
  CALL New(TmKxArr,NPrc)
  CALL New(TmMLArr,NPrc)
  CALL New(TmTMArr,NPrc)
  CALL New(TmALArr,NPrc)
  CALL New(TmDLArr,NPrc)
  !
  CALL MPI_Gather(TmKx,1,MPI_DOUBLE_PRECISION,TmKxArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmML,1,MPI_DOUBLE_PRECISION,TmMLArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmTM,1,MPI_DOUBLE_PRECISION,TmTMArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmAL,1,MPI_DOUBLE_PRECISION,TmALArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL MPI_Gather(TmDL,1,MPI_DOUBLE_PRECISION,TmDLArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  !
  IF(MyID.EQ.ROOT) THEN
     !
     ! Imbalance stuff.
     CALL PImbalance(TmKxArr ,NPrc,Prog_O='ComputeK')
     !CALL PImbalance(TmKTArr,NPrc,Prog_O='GONX'     )
     !
     WRITE(*,1001) SUM(TmALArr%D )/DBLE(NPrc),MINVAL(TmALArr%D ),MAXVAL(TmALArr%D )
     WRITE(*,1002) SUM(TmMLArr%D )/DBLE(NPrc),MINVAL(TmMLArr%D ),MAXVAL(TmMLArr%D )
     WRITE(*,1003) SUM(TmTMArr%D )/DBLE(NPrc),MINVAL(TmTMArr%D ),MAXVAL(TmTMArr%D )
     WRITE(*,1004) SUM(TmKxArr%D )/DBLE(NPrc),MINVAL(TmKxArr%D ),MAXVAL(TmKxArr%D )
     WRITE(*,1005) SUM(TmDLArr%D )/DBLE(NPrc),MINVAL(TmDLArr%D ),MAXVAL(TmDLArr%D )
     !WRITE(*,1006) SUM(NERIsArr%D)           ,MINVAL(NERIsArr%D),MAXVAL(NERIsArr%D)
     !
1001 FORMAT(' ONX: Ave TmAL = ',F15.2,', Min TmAL = ',F15.2,', Max TmAL = ',F15.2)
1002 FORMAT(' ONX: Ave TmML = ',F15.2,', Min TmML = ',F15.2,', Max TmML = ',F15.2)
1003 FORMAT(' ONX: Ave TmTM = ',F15.2,', Min TmTM = ',F15.2,', Max TmTM = ',F15.2)
1004 FORMAT(' ONX: Ave TmKx = ',F15.2,', Min TmKx = ',F15.2,', Max TmKx = ',F15.2)
1005 FORMAT(' ONX: Ave TmDL = ',F15.2,', Min TmDL = ',F15.2,', Max TmDL = ',F15.2)
     !1006 FORMAT(' ONX: Tot ERI  = ',F15.2,', Min ERI  = ',F15.2,', Max ERI  = ',F15.2)
  ENDIF
  !
  CALL Delete(TmALArr)
  CALL Delete(TmMLArr)
  CALL Delete(TmTMArr)
  CALL Delete(TmKxArr)
  CALL Delete(TmDLArr)
  !
#else
  !
  WRITE(*,1001) TmAL
  WRITE(*,1002) TmML
  WRITE(*,1003) TmTM
  WRITE(*,1004) TmKx
  WRITE(*,1005) TmDL
  !WRITE(*,1006) ....
  !
1001 FORMAT(' ONX: Ave TmAL = ',F15.2)
1002 FORMAT(' ONX: Ave TmML = ',F15.2)
1003 FORMAT(' ONX: Ave TmTM = ',F15.2)
1004 FORMAT(' ONX: Ave TmKx = ',F15.2)
1005 FORMAT(' ONX: Ave TmDL = ',F15.2)
  !1006 FORMAT(' ONX: Tot ERI  = ',F15.2)
  !
#endif
#endif
  !
  !------------------------------------------------
  ! Deallocate Arrays
  !
  CALL Delete(Kx)
  CALL Delete(OffArrC)
  CALL Delete(OffArrP)
  !
  CALL Delete(BSc)
  CALL Delete(GMc)
  CALL Delete(BSp)
  CALL Delete(GMp)
  !
  CALL ShutDown(Prog)
  !
END PROGRAM ONX2
