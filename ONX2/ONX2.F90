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
  USE ONXRng    , ONLY: RangeOfExchangeBCSR
  USE ONXFillOut, ONLY: FillOutBCSR

USE ONXCtrSclg, ONLY: TrnMatBlk
  !
  !New
  USE ONX2ComputK
  !
  IMPLICIT NONE
  !
  TYPE(BCSR)                     :: D
  TYPE(BCSR)                     :: Kx,T1,T2
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
!New
  TYPE(CList ), DIMENSION(:), POINTER :: ListC
  TYPE(CList2), DIMENSION(:), POINTER :: ListC2

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

  CALL Get(D,TrixFile('D',Args,0))

  !
  CALL RangeOfExchangeBCSR(BSc,GMc,BSp,GMp,D)

  CALL TrnMatBlk(BSp,GMp,D)

  !
  CALL New(Kx,(/NRows+1,NCols,NElem/))
  CALL SetEq(Kx%MTrix,Zero)
  CALL InitK(BSc,GMc,Kx)
  !

  WRITE(*,*) 'WE ARE IN ONX2 WE ARE IN ONX2 WE ARE IN ONX2 WE ARE IN ONX2 WE ARE IN ONX2'

  WRITE(*,*) 'allocate List'
  CALL CPU_TIME(Time1)
  CALL AllocList(ListC2,NAtoms)
  CALL CPU_TIME(Time2)
  WRITE(*,*) 'allocate List: ok',Time2-Time1


  WRITE(*,*) 'make List'
  CALL CPU_TIME(Time1)
  CALL MakeList(ListC2,GMc,BSc,CS_OUT)
  CALL CPU_TIME(Time2)
  WRITE(*,*) 'make List: ok',Time2-Time1

  WRITE(*,*) 'Print List'
!  CALL PrintList2(ListC2)
!  write(*,*)'========================'
!  CALL PrintList(ListC)
  WRITE(*,*) 'Print List:ok'

  WRITE(*,*) 'Compute Kx'
  CALL CPU_TIME(Time1)
  CALL ComputK(D,Kx,ListC2,ListC2,GMc,BSc,CS_OUT)
!  CALL ComputK(D,K,ListC,ListC,GMc,BSc,CS_OUT)
  CALL CPU_TIME(Time2)
  WRITE(*,*) 'Compute Kx:ok',Time2-Time1

  WRITE(*,*) 'deallocate List'
  CALL CPU_TIME(Time1)
  CALL DeAllocList(ListC2)
!  CALL DeAllocList(ListC)
  CALL CPU_TIME(Time2)
  WRITE(*,*) 'deallocate List:ok',Time2-Time1

  CALL FillOutBCSR(BSc,GMc,Kx)

  !CALL Print_BCSR(Kx,'Kx',Unit_O=6)

  CALL TrnMatBlk(BSc,GMc,Kx)

  !
  ! Save on disc.
  CALL Put(Kx,TrixFile('K',Args,0))
  CALL PChkSum(Kx,'Kx['//TRIM(SCFCycl)//']',Prog)
  CALL PPrint( Kx,'Kx['//TRIM(SCFCycl)//']')
  CALL Plot(   Kx,'Kx['//TRIM(SCFCycl)//']')
  !
!--------------------------------------------------------------------------------
! Clean up...
!--------------------------------------------------------------------------------
  !
  CALL Delete(D)
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


