PROGRAM GONX2
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


  USE ONXCtrSclg, ONLY: TrnMatBlk


  !New
  USE ONX2DataType
  USE ONX2List, ONLY: AllocList,DeAllocList,MakeGList,PrintList
  USE GONX2ComputDK, ONLY: ComputDK
  !
  IMPLICIT NONE
  !
!#ifdef PARALLEL
!#else
  TYPE(BCSR)                     :: D
!#endif
  TYPE(BSET)                     :: BSc
  TYPE(CRDS)                     :: GMc
  TYPE(BSET)                     :: BSp
  TYPE(CRDS)                     :: GMp
  TYPE(ARGMT)                    :: Args
  TYPE(INT_VECT)                 :: Stat


  TYPE(DBL_RNK2)   :: GradX,GradAux

!--------------------------------------------------------------------------------
! Misc. variables and parameters...
!--------------------------------------------------------------------------------
  INTEGER                        :: OldFileID
  REAL(DOUBLE)                   :: time1,time2
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: InFile,RestartHDF
  CHARACTER(LEN=*),PARAMETER     :: Prog='GONX2'
!--------------------------------------------------------------------------------

  TYPE(CList2), DIMENSION(:), POINTER :: ListC,ListD



#ifdef PARALLEL
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
  WRITE(*,*) '-------- We are in GONX2 --------'


!!$  write(*,*) '1x',GMc%Carts%D(1,1)
!!$  write(*,*) '1y',GMc%Carts%D(2,1)
!!$  write(*,*) '1z',GMc%Carts%D(3,1)
!!$  write(*,*) '2x',GMc%Carts%D(1,2)
!!$  write(*,*) '2y',GMc%Carts%D(2,2)
!!$  write(*,*) '2z',GMc%Carts%D(3,2)
!!$  write(*,*) '3x',GMc%Carts%D(1,3)
!!$  write(*,*) '3y',GMc%Carts%D(2,3)
!!$  write(*,*) '3z',GMc%Carts%D(3,3)

!!$  GMc%Carts%D(1,1)
!!$  GMc%Carts%D(2,1)
!!$  GMc%Carts%D(3,1)
!!$  GMc%Carts%D(1,2)
!!$  GMc%Carts%D(2,2)
!  GMc%Carts%D(3,NAtoms)=GMc%Carts%D(3,NAtoms)+1.D-3

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
!!$#ifdef PARALLEL
!!$     CALL PDrv_Initialize(DFastMat,TrixFile('D',Args,0),'ONXPart',Args)
!!$#else
     CALL Get(D,TrixFile('D',Args,0))
!!$#endif
!!$  !ENDIF
!!$  END SELECT
  !
  !
!#ifdef PARALLEL
!  CALL TrnMatBlk(BSp,GMp,DFastMat)
!#else
  CALL TrnMatBlk(BSp,GMp,D       )
!#endif
  !
  !
  WRITE(*,*) 'allocate List'
!#ifdef PARALLEL
!  Time1 = MPI_WTIME()
!  CALL AllocList(ListC,NAtoms) !!!!!!!!!!!!!!!!!! Add atom1-atom2 or sthg like that !!!!!!!!!!!!!!!!!!
!  CALL AllocList(ListD,NAtoms) !!!!!!!!!!!!!!!!!! Add atom1-atom2 or sthg like that !!!!!!!!!!!!!!!!!!
!  Time2 = MPI_WTIME()
!#else
  CALL CPU_TIME(Time1)
  CALL AllocList(ListC,NAtoms)
  CALL CPU_TIME(Time2)
!#endif
  WRITE(*,*) 'allocate List: ok',Time2-Time1
  !
  !------------------------------------------------
  ! Make the distribution list(s).
  WRITE(*,*) 'make List'
!#ifdef PARALLEL
!  Time1 = MPI_WTIME()
!  CALL MakeList(ListC,GMc,BSc,CS_OUT) !!!!!!!!!!!!!!!!!! Add atom list !!!!!!!!!!!!!!!!!!
!  CALL MakeList(ListD,GMc,BSc,CS_OUT) !!!!!!!!!!!!!!!!!! Add atom list !!!!!!!!!!!!!!!!!!
!  Time2 = MPI_WTIME()
!#else
  CALL CPU_TIME(Time1)
  CALL MakeGList(ListC,GMc,BSc,CS_OUT) !Change that in MakeGList
  CALL CPU_TIME(Time2)
!#endif
  !TmML = Time2-Time1
  WRITE(*,*) 'make List: ok',Time2-Time1
  !
  !------------------------------------------------
  !
#ifdef ONX2_DBUG
  WRITE(*,*) 'Print List'
!#ifdef PARALLEL
!  CALL PrintList(ListC)
!  CALL PrintList(ListD)
!#else
  CALL PrintList(ListC)
!#endif
  WRITE(*,*) 'Print List:ok'
#endif
  !
  !
  !------------------------------------------------
  !
  CALL New(GradX  ,(/3,NAtoms/))
  CALL DBL_VECT_EQ_DBL_SCLR(3*NAtoms,GradX%D(1,1),0.0d0)
  !
  !
  !------------------------------------------------
  ! Compute Exchange Forces.
  WRITE(*,*) 'DKx'
!#ifdef PARALLEL
!  Time1 = MPI_WTIME()
!  CALL ComputK(DFastMat,KxFastMat,ListC,ListD,GMc,BSc,CS_OUT)
!  Time2 = MPI_WTIME()
!#else
  CALL CPU_TIME(Time1)
  CALL ComputDK(D,GradX,ListC,ListC,GMc,BSc,CS_OUT)
  CALL CPU_TIME(Time2)
!#endif
  !TmKx = Time2-Time1
  WRITE(*,*) 'DKx:ok',Time2-Time1

!!$  write(*,*) GradX%D
!!$  write(*,*) 'sum(GradX%D)=',sum(GradX%D)

  !
  !------------------------------------------------
  ! Free up some space. Deallocate the list(s).
  WRITE(*,*) 'deallocate List'
!#ifdef PARALLEL
!  Time1 = MPI_WTIME()
!  CALL DeAllocList(ListC)
!  CALL DeAllocList(ListD)
!  Time2 = MPI_WTIME()
!#else
  CALL CPU_TIME(Time1)
  CALL DeAllocList(ListC)
  CALL CPU_TIME(Time2)
!#endif
  WRITE(*,*) 'deallocate List:ok',Time2-Time1



  !
  !------------------------------------------------
  ! Redistribute partition informations.
!#ifdef PARALLEL
!  CALL PDrv_Finalize(DFastMat,CollectInPar_O=.TRUE.)
!  CALL Delete_FastMat1(DFastMat)
!#else
  CALL Delete(D)
!#endif
  !
  !



  !
  !------------------------------------------------
  !
  CALL New(GradAux,(/3,NAtoms/))

  !CALL Get(GradAux,'gradients',Tag_O=chGEO)
  CALL Get(GradAux,'gradients',Tag_O=CurGeom)
  GradAux%D=GradX%D+GradAux%D     
  !CALL Put(GradAux,'gradients',Tag_O=chGEO)
  CALL Put(GradAux,'gradients',Tag_O=CurGeom)

!!$  write(*,*) GradAux%D
!!$  write(*,*) 'sum(GradAux%D)=',sum(GradAux%D)

  CALL Delete(GradX  )
  CALL Delete(GradAux)


  !
  CALL Delete(BSc   )
  CALL Delete(GMc   )
  CALL Delete(BSp   )
  CALL Delete(GMp   )
  !
  CALL ShutDown(Prog)
  !
END PROGRAM GONX2

