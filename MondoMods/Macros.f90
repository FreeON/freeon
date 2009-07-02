!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
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
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------

MODULE Macros
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalObjects
  USE GlobalCharacters
  USE ParsingConstants
  USE PrettyPrint
  USE Parse
  USE Clock
  USE InOut
  USE Functionals
  USE AtomPairs
  USE Mechanics
#ifdef PARALLEL
  USE MondoMPI
#endif
#ifdef NAG
  USE F90_UNIX
#endif

  IMPLICIT NONE

  INTERFACE Init
    MODULE PROCEDURE Init_TIME, Init_DEBG, Init_MEMS
  END INTERFACE
  INTEGER HDFFileID,H5GroupID
  !-------------------------------------------------------------------------------
CONTAINS
  SUBROUTINE StartUp(Args,Prog,Serial_O)
    TYPE(ARGMT),INTENT(OUT)              :: Args
    CHARACTER(LEN=*),INTENT(IN)          :: Prog
    LOGICAL,OPTIONAL                     :: Serial_O
    INTEGER                              :: I
    TYPE(INT_VECT)                       :: SpaceTimeSplit

    REAL(DOUBLE)                         :: ETag
#ifdef PARALLEL
    LOGICAL                              :: Serial
    INTEGER                              :: ChkNPrc,iTAG
    CHARACTER(LEN=DCL)                   :: MONDO_HOST
    !-------------------------------------------------------------------------------
    CALL InitMPI()
#endif
    ! Get arguments
    CALL Get(Args)

    ! Get SCRATCH directory from last argument.
    !CALL GetEnv('MONDO_SCRATCH', MONDO_SCRATCH)
    !if(LEN(TRIM(MONDO_SCRATCH)) == 0) then
      !MONDO_SCRATCH = HAVE_MONDO_SCRATCH
    !endif
    !MONDO_SCRATCH=TRIM(MONDO_SCRATCH)//'/'
    MONDO_SCRATCH = Args%C%C(SIZE(Args%C%C, 1))

    ! The HDF5 file name
    H5File=TRIM(MONDO_SCRATCH)//TRIM(Args%C%C(1))//TRIM(InfF)
    !WRITE(*,*) "setting H5File to "//TRIM(H5File)
    InfFile=H5File
    ! Open the HDF file
    HDFFileID=OpenHDF(H5File)
    ! This is the global IO file ID
    HDF_CurrentID=HDFFileID
    ! Mark prog for failure
    CALL MarkFailure(Prog)
    ! Load global SCF status strings
    CALL LoadTopLevelGlobals(Args)
    ! Parse this programs debug level
    CALL Init(PrintFlags,Prog)
#ifdef PARALLEL
    ! Get the space-time parallel topology
    CALL New(SpaceTimeSplit,3)
    CALL Get(SpaceTimeSplit,'SpaceTime')
    CALL CloseHDF(HDFFileID)

    ! Create Cartesian topology for the clones
    MyClone=CartCommSplit(SpaceTimeSplit,Serial_O)
    IF(NClones>1)THEN
       CurClone=IntToChar(MyClone)
    ELSE
       CurClone=""
    ENDIF
    ! CALL AlignNodes(' After CartCommSplit ')
#else
    CALL Get(MyClone,'SpaceTime')
    CurClone = IntToChar(MyClone)
    CALL CloseHDF(HDFFileID)
#endif
    ! Each ROOT in each MONDO_COMM opens the HDF file --->FOR READ ONLY<---
    HDFFileID=OpenHDF(H5File)
    HDF_CurrentID=HDFFileID
    ! Operate at the top level of the archive
    ! Open a group that the ROOT of each clone accesses by default
    H5GroupID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(MyClone)))
    ! Default is now the cloned group id rather than the HDF file id
    HDF_CurrentID=H5GroupID
    !    WRITE(*,22)TRIM(Prog),MyId,MyClone
    ! 22 format(A10,' IN MACROS, GETTING TAG',I3,", ",I3)
    ! 23 format(A10,' IN MACROS, GOT     TAG',I3,", ",I3)
    !    CALL GET(ITAG,"TAG")
    !    WRITE(*,23)TRIM(Prog),MyId,MyClone
    !    CALL GET(ETAG,'E_NuclearTotal',StatsToChar(Current))
    !    CALL PUT(2*ITAG,"TAG")
    !    CALL PUT(2D0*ETAG,'E_NuclearTotal',StatsToChar(Current))
    !    IF(MyID==ROOT)WRITE(*,33),MyId,MyClone,iTAG,ETag,HDFFileID,HDF_CurrentID
    !    33 FORMAT(' MyId = ',I2,', MyClone = ',I2,', iTAG = ',I6,' ETAG = ',D12.6,', FileID = ',I10,' CURRENTID = ',I10)
    ! Load variables global at the group (clone) level
    CALL LoadGroupGlobals(Args)
    ! Initialize memory statistics
    CALL Init(MemStats)
    ! Start the clock ...
    CALL Elapsed_Time(PerfMon,'Init')
    ! Print time stamp
    IF(PrintFlags%Key>DEBUG_MEDIUM)THEN
#ifdef PARALLEL
      IF(MyId==ROOT) THEN
        CALL TimeStamp('Entering '//TRIM(Prog),Enter_O=.TRUE.)
      ENDIF
#else
      CALL TimeStamp('Entering '//TRIM(Prog))
#endif
    ENDIF
    !
    !   CALL ShutDown(Prog)
  END SUBROUTINE StartUp
#ifdef PARALLEL
  FUNCTION CartCommSplit(SpaceTime,Serial_O) RESULT(MyClone)
    INTEGER               :: IErr,MyClone,CART_COMM
    TYPE(INT_VECT)        :: SpaceTime
    LOGICAL,OPTIONAL      :: Serial_O
    CHARACTER(LEN=10)     :: Sub='SpaceNTime'
    LOGICAl               :: AllButRootMustDie
    INTEGER,DIMENSION(2)  :: Local
    !-----------------------------------------------------------------------!
    ! Determine parallel status
    IF(PRESENT(Serial_O))THEN
      IF(Serial_O)THEN
        InParallel=.FALSE.
      ELSE
        InParallel=.TRUE.
      ENDIF
    ELSE
      InParallel=.FALSE.
    ENDIF
    ! Create a SpaceTime%I(1) x SpaceTime%I(2) Cartesian communicator
    IF(InParallel)THEN
      CALL MPI_CART_CREATE(MPI_COMM_WORLD,2,(/SpaceTime%I(1),SpaceTime%I(2)/),  &
           (/.FALSE.,.FALSE./),.TRUE.,CART_COMM,IErr)
    ELSE
      CALL MPI_CART_CREATE(MPI_COMM_WORLD,2,(/1,SpaceTime%I(2)/),  &
           (/.FALSE.,.FALSE./),.TRUE.,CART_COMM,IErr)
    ENDIF
    CALL ErrChk(IErr,Sub)
    IF(CART_COMM/=MPI_COMM_NULL)THEN
      ! Find out which row (group) this PE belongs to
      CALL MPI_CART_COORDS(CART_COMM,MyId,2,Local,IErr)
      CALL ErrChk(IErr,Sub)
    ELSE
      ! This is a dead node
      CALL ShutDown('DeaddNode')
    ENDIF
    ! Offset the actuall clone
    MyClone=SpaceTime%I(3)+Local(2)
    !CALL AlignNodes('MyClone = '//TRIM(IntToChar(MyClone)))
    !
    ! Now split into SpaceTime%I(1) rows. Each row has SpaceTime%I(2) processors
    ! parallel in the spatial domain and using MONDO_COMM as their
    ! default communicator
    CALL MPI_CART_SUB(CART_COMM,(/.TRUE.,.FALSE./),MONDO_COMM,IErr)
    CALL ErrChk(IErr,Sub)
    ! Reload local rank and PE number for the new MONDO_COMM
    MyID=MRank()
    NPrc=MSize()
  END FUNCTION CartCommSplit
#endif
  SUBROUTINE ShutDown(Prog)
    CHARACTER(LEN=*),INTENT(IN) :: Prog

    !    IF(MyID==ROOT)WRITE(*,33)Prog,MyId,MyClone,H5GroupID,HDFFileID
    !    33 FORMAT(A10,'SHUTDOWN: MyId = ',I2,', MyClone = ',I2,', GroupID = ',I10,' FileID = ',I10)

    ! Following is a bit of fancy HDF footwork, nessesary for parallel-clones to work properly
    !
    ! Close the clone directories
    CALL CloseHDFGroup(H5GroupID)
#ifdef PARALLEL
    ! Close the global HDF file, possibly for multiple space-time root nodes
    CALL CloseHDF(HDFFileID)
    ! Revert back to global communicator, rank etc ...
    MONDO_COMM=MPI_COMM_WORLD
    MyID=MRank()
    ! and reopen the upper level HDF directory for just the world root node
    HDFFileID=OpenHDF(H5File)
#endif
    ! Here is the global file ID for the upper level HDF direcory, for the world root node
    HDF_CurrentID=HDFFileID
    IF(HasQM()) THEN
      CALL Delete(BSiz)
      CALL Delete(OffS)
    ENDIF
#ifdef PARALLEL
    IF(InParallel)THEN
      CALL Delete(Beg)
      CALL Delete(End)
    ENDIF
#endif

    IF(PrintFlags%Key>=DEBUG_MEDIUM)THEN
      CALL Elapsed_TIME(PerfMon,'Accum',Proc_O=Prog)
      CALL PPrint(PerfMon,Prog)
      IF(.FALSE.) THEN
        CALL Elapsed_TIME(PerfMon,'Accum',Proc_O=Prog)
        CALL PPrint(PerfMon,Prog,BareBones_O=.TRUE.)
      ENDIF
    ELSE
      IF(.FALSE.) THEN
        CALL Elapsed_TIME(PerfMon,'Accum',Proc_O=Prog)
        CALL PPrint(PerfMon,Prog,BareBones_O=.TRUE.)
      ENDIF
    ENDIF

    ! OMG,this is anoying!!
    ! Check MemStats.
!    IF(MemStats%Allocs /= MemStats%DeAllocs) THEN
!      CALL MondoLog(DEBUG_NONE, "ShutDown ("//TRIM(Prog)//")", "Possible memory leak. Allocs = " &
!        //TRIM(IntToChar(MemStats%Allocs))//", DeAllocs = " &
!        //TRIM(IntToChar(MemStats%DeAllocs)))
      !CALL Halt("[FIXME]")
!    ENDIF

    IF(PrintFlags%Key==DEBUG_MAXIMUM) THEN
      CALL PPrint(MemStats,Prog)
    ENDIF

    ! Now mark sucess of this program ...
    CALL Put(.FALSE.,'ProgramFailed')
    ! ... and close the HDF file ...
    CALL CloseHDF(HDFFileID)
    ! ... shutdown MPI and print a time stamp
#ifdef PARALLEL
    IF(PrintFlags%Key>DEBUG_MEDIUM.AND.MyId==ROOT)  &
         CALL TimeStamp('Exiting '//TRIM(Prog),Enter_O=.FALSE.)
    ! Shutdown MPI
    CALL FiniMPI()
#else
    IF(PrintFlags%Key>DEBUG_MEDIUM)  &
         CALL TimeStamp('Exiting '//TRIM(Prog),Enter_O=.FALSE.)
#endif
    STOP
  END SUBROUTINE ShutDown
  !==============================================================
  ! LOAD GLOBAL VARIABLES FROM THE TOP LEVEL OF THE HDF FILE
  !==============================================================
  SUBROUTINE LoadTopLevelGlobals(Args)
    TYPE(ARGMT) :: Args
    INTEGER :: I,ChkNPrc
    !-----------------------------------------------------------!
    IF(Args%NC>=2) THEN
      SCFActn=TRIM(Args%C%C(2))
    ENDIF
    Current=Args%I%I(1:3)
    Previous=Args%I%I(4:6)
    ! SCF Cycle
    SCFCycl=TRIM(IntToChar(Args%I%I(1)))
    IF(Args%I%I(3)-1<1)THEN
      PrvCycl=TRIM(IntToChar(Previous(1)))
    ELSE
      PrvCycl=TRIM(IntToChar(Args%I%I(1)-1))
    ENDIF
    NxtCycl=TRIM(IntToChar(Args%I%I(1)+1))
    ! Geometry
    CurGeom=TRIM(IntToChar(Args%I%I(3)))
    IF(Args%I%I(3)-1<1)THEN
      PrvGeom=TRIM(IntToChar(Previous(3)))
    ELSE
      PrvGeom=TRIM(IntToChar(Args%I%I(3)-1))
    ENDIF
    NxtGeom=TRIM(IntToChar(Args%I%I(3)+1))

    ! Basis
    CurBase = TRIM(IntToChar(Args%I%I(2)))
    PrvBase = TRIM(IntToChar(Previous(2)))
    ScrName = TRIM(MONDO_SCRATCH)//TRIM(Args%C%C(1))
    PWDName = TRIM("")

    ! Load global file names
    CALL Get(RecycleHDF, "RecycleHDF")
    CALL Get(logFile,'logfile')
    CALL Get(OutFile,'outputfile')
    CALL Get(InpFile,'inputfile')
    CALL Get(Restart,'OldInfo')
    CALL Get(MaxAtms,'maxatms',Tag_O=CurBase)
    CALL Get(MaxBlks,'maxblks',Tag_O=CurBase)
    CALL Get(MaxNon0,'maxnon0',Tag_O=CurBase)
#ifdef PARALLEL
    IF(InParallel)THEN
      CALL Get(MaxAtmsNode,'maxatmsnode', Tag_O=CurBase)
      CALL Get(MaxBlksNode,'maxblksnode', Tag_O=CurBase)
      CALL Get(MaxNon0Node,'maxnon0node', Tag_O=CurBase)
    ENDIF
#endif
    CALL Get(ModelChem,'ModelChemistry',Tag_O=CurBase)
    CALL Get(NClones,'clones')

  END SUBROUTINE LoadTopLevelGlobals
  !==============================================================
  ! LOAD GLOBAL VARIABLES FROM EACH CLONE/GROUP OF THE HDF FILE
  !==============================================================
  SUBROUTINE LoadGroupGlobals(Args)
    TYPE(ARGMT) :: Args
    INTEGER :: I,ChkNPrc
    !-----------------------------------------------------------!
#ifdef MMech
    IF(HasQM())THEN
#endif
      CALL Get(NEl   ,'nel'     ,Tag_O=CurGeom)
      CALL Get(NAlph ,'nelalpha',Tag_O=CurGeom)
      CALL Get(NBeta ,'nelbeta' ,Tag_O=CurGeom)
      CALL Get(TotCh ,'charge'  ,Tag_O=CurGeom)
      CALL Get(NAtoms,'natoms'  ,Tag_O=CurGeom)
      CALL New(BSiz,NAtoms)
      CALL New(OffS,NAtoms)
      CALL Get(NBasF,'nbasf',Tag_O=CurBase)
      CALL Get(BSiz,'atsiz',Tag_O=CurBase)
      CALL Get(OffS,'atoff',Tag_O=CurBase)
      ! Global value for max block size
      MaxBlkSize=0
      DO I=1,NAtoms
        MaxBlkSize=MAX(MaxBlkSize,BSiz%I(I))
      ENDDO
#ifdef PARALLEL
      IF(InParallel)THEN
        CALL New(OffSt,NPrc-1,0)
        CALL Get(OffSt,'dbcsroffsets',Tag_O=CurBase)
        !             CALL Get(ChkNPrc,'chknprc')
        !             IF(NPrc/=ChkNPrc) &
        !                  CALL Halt(' In StartUp() --- Inconsistency: NPrc = '  &
        !                  //TRIM(IntToChar(NPrc))//' ChkNPrc = ' &
        !                  //TRIM(IntToChar(ChkNPrc)))
        CALL New(Beg,NPrc-1,0)
        CALL New(End,NPrc-1,0)
        CALL Get(Beg,'beg',Tag_O=CurBase)
        CALL Get(End,'end',Tag_O=CurBase)
      ENDIF
#endif
#ifdef MMech
    ENDIF
#endif
    ! Load QM/MM switches
    CALL InitMMech()
    ! Load global thresholding values
    CALL SetThresholds(CurBase)
    !
    ! Note that these cell sets are duplicated in the geometry object CRDS
    ! They are maintained here ONLY for legacy purposes, and should not be
    ! used in future programming efforts.  They have been entirely removed
    ! in QCTC2, and MondoSCF2

    !    CALL Get(CS_IN ,'incells' ,Tag_O=CurBase)
    !    CALL Get(CS_OUT,'ovcells',Tag_O=CurBase)

    CALL Get(CS_IN ,'incells' ,Tag_O=CurGeom)
    CALL Get(CS_OUT,'ovcells',Tag_O=CurGeom)

  END SUBROUTINE LoadGroupGlobals
  !=========================================================
  ! MARK FAILURE OF PROG
  !=========================================================
  SUBROUTINE MarkFailure(Prog)
    CHARACTER(LEN=*) :: Prog
    CALL Put(.TRUE.,'ProgramFailed')
    CALL Put(Prog  ,'FailedProgram')
  END SUBROUTINE MarkFailure
  !=========================================================
  ! MARK SUCCESS OF PROG
  !=========================================================
  SUBROUTINE MarkSuccess()
    CALL Put(.FALSE.,'ProgramFailed')
  END SUBROUTINE MarkSuccess
  !-------------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------------
  SUBROUTINE Init_MEMS(A)
    TYPE(MEMS), INTENT(OUT) :: A
    A%Allocs=0
    A%DeAllocs=0
    A%MemTab=0
    A%MaxMem=0
    A%MaxAlloc=0
  END SUBROUTINE Init_MEMS
  !-------------------------------------------------------------------------------
  !
  !-------------------------------------------------------------------------------
  SUBROUTINE Init_TIME(A)
    TYPE(TIME),INTENT(OUT) :: A
#ifdef PARALLEL
    IF(InParallel)  &
         CALL AlignNodes()
#endif
    A%FLOP=Zero
    A%CPUS=CPUSec()
    A%Wall=WallSec()
  END SUBROUTINE Init_TIME
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE Init_DEBG(A,Prog)
    TYPE(DEBG), INTENT(OUT)      :: A
    CHARACTER(LEN=*), INTENT(IN) :: Prog
#ifdef PARALLEL
    IF(MyId==ROOT)THEN
#endif
      CALL OpenASCII(InpFile,Inp)
      IF(OptKeyQ(Inp,GLOBAL_DEBUG,DBG_PRT_INTS))THEN
        A%Int=DEBUG_INTEGRAL
      ELSE
        A%Int=DEBUG_NONE
      ENDIF
      IF(OptKeyQ(Inp,TRIM(Prog)  ,DBG_NONE).OR.         &
           OptKeyQ(Inp,GLOBAL_DEBUG,DBG_NONE) )THEN
        A%Key=DEBUG_NONE
      ELSEIF(OptKeyQ(Inp,TRIM(Prog)  ,DBG_MEDIUM).OR.   &
           OptKeyQ(Inp,GLOBAL_DEBUG,DBG_MEDIUM) )THEN
        A%Key=DEBUG_MEDIUM
      ELSEIF(OptKeyQ(Inp,TRIM(Prog),  DBG_MAXIMUM).OR.  &
           OptKeyQ(Inp,GLOBAL_DEBUG,DBG_MAXIMUM) )THEN
        A%Key=DEBUG_MAXIMUM
      ELSE
        A%Key=DEBUG_MINIMUM
      ENDIF

      IF(OptKeyQ(Inp,TRIM(Prog),  DBG_MATRICES).OR. &
           OptKeyQ(Inp,GLOBAL_DEBUG,DBG_MATRICES) )THEN
        A%Mat=DEBUG_MATRICES
      ELSEIF(OptKeyQ(Inp,TRIM(Prog),  PLT_MATRICES).OR. &
           OptKeyQ(Inp,GLOBAL_DEBUG,PLT_MATRICES) )THEN
        A%Mat=PLOT_MATRICES
      ELSE
        A%Mat=DEBUG_NONE
      ENDIF
      !
      IF(OptKeyQ(Inp,GLOBAL_DEBUG,DBG_PRT_MM) .OR. &
           OptKeyQ(Inp,TRIM(Prog)  ,DBG_PRT_MM)) THEN
        A%MM=DEBUG_MM
      ELSEIF(OptKeyQ(Inp,GLOBAL_DEBUG,DBG_PRT_FRC) .OR. &
           OptKeyQ(Inp,TRIM(Prog)  ,DBG_PRT_FRC)) THEN
        A%MM=DEBUG_FRC
      ELSE
        A%MM=DEBUG_NONE
      ENDIF
      !
      IF(OptKeyQ(Inp,TRIM(Prog),  DBG_CHKSUMS).OR. &
           OptKeyQ(Inp,GLOBAL_DEBUG,DBG_CHKSUMS))THEN
        A%Chk=DEBUG_CHKSUMS
      ELSE
        A%Chk=DEBUG_NONE
      ENDIF
      IF(OptKeyQ(Inp,TRIM(Prog),  DBG_MMA_STYLE).OR.     &
           OptKeyQ(Inp,GLOBAL_DEBUG,DBG_MMA_STYLE) )THEN
        A%Fmt=DEBUG_MMASTYLE
      ELSEIF(OptKeyQ(Inp,TRIM(Prog),  DBG_FLT_STYLE).OR. &
           OptKeyQ(Inp,GLOBAL_DEBUG,DBG_FLT_STYLE) )THEN
        A%Fmt=DEBUG_FLTSTYLE
      ELSE
        A%Fmt=DEBUG_DBLSTYLE
      ENDIF
#ifdef PARALLEL
    ENDIF
    IF(InParallel)CALL BCast(A)
#endif
    CLOSE(Inp)
  END SUBROUTINE Init_DEBG
END MODULE Macros
