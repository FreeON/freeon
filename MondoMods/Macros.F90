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
         CHARACTER(LEN=DCL)                   :: H5File
         LOGICAL,OPTIONAL                     :: Serial_O
         INTEGER                              :: I
#ifdef PARALLEL
         LOGICAL                              :: Serial
         INTEGER                              :: ChkNPrc,MyClone
         TYPE(INT_VECT)                       :: SpaceTimeSplit
	 CHARACTER(LEN=DCL)                   :: MONDO_HOST 
!-------------------------------------------------------------------------------
         IF(PRESENT(Serial_O))THEN
            IF(.NOT.Serial_O)CALL InitMPI()
         ENDIF


!	 CALL GetEnv('MONDO_HOST',MONDO_HOST)
!	 IF(InParallel)CALL AlignNodes(TRIM(MONDO_HOST))
#endif
!        Get arguments
         CALL Get(Args)
!        Get SCRATCH directoryfrom env
         CALL GetEnv('MONDO_SCRATCH',MONDO_SCRATCH)
         MONDO_SCRATCH=TRIM(MONDO_SCRATCH)//'/'
         IF(LEN(TRIM(MONDO_SCRATCH))==0)CALL Halt(' $(MONDO_SCRATCH) not set.')
!        The HDF5 file name 
         H5File=TRIM(MONDO_SCRATCH)//TRIM(Args%C%C(1))//TRIM(InfF)
         InfFile=H5File
         ! Open the HDF file
         HDFFileID=OpenHDF(H5File)
         ! This is the global IO file ID
         HDF_CurrentID=HDFFileID
#ifdef PARALLEL_CLONES
         ! Get the space-time parallel topology
         CALL Get(SpaceTimeSplit,'spacetimesplit')
         CALL CloseHDF(HDFFileID)
         ! Create Cartesian topology for the clones         
         MyClone=CartCommSplit(SpaceTimeSplit)
         ! Each ROOT in each MONDO_COMM opens the HDF file 
         HDFFileID=OpenHDF(H5File)
         ! Now we open a group that the ROOT of each clone accesses by default 
         H5GroupID=OpenHDFGroup(HDFFileID,"clone"//TRIM(IntToChar(MyClone)))
         ! Default is now the cloned group id rather than the HDF file id 
         HDF_CurrentID=H5GroupID
         ! Now back to busy as usual 
#endif
         ! Mark for failure
         CALL MarkFailure(Prog)
         ! Load global SCF status strings
         CALL Init_SCFStat(Args)
         ! Load QM/MM switches
         CALL InitMMech()
         ! Load global thresholding values
         CALL SetThresholds(CurBase)
         ! Initialize memory statistics
         CALL Init(MemStats)
         ! Load some more global variables
         CALL Init_Globals(Args)
         ! Initialize debug level
         CALL Init(PrintFlags,Prog)
         PrintFlags%Key=DEBUG_MAXIMUM
         ! Start the clock ...
         CALL Elapsed_Time(PerfMon,'Init')
!        Print time stamp
         IF(PrintFlags%Key>DEBUG_MEDIUM)THEN
#ifdef PARALLEL
            IF(MyId==ROOT) &
               CALL TimeStamp(' Entering '//TRIM(Prog),Enter_O=.TRUE.)
#else
            CALL TimeStamp('Entering '//TRIM(Prog))
#endif 
         ENDIF
       END SUBROUTINE StartUp
       !-------------------------------------------------------------------------------
       SUBROUTINE ShutDown(Prog)
         CHARACTER(LEN=*),INTENT(IN) :: Prog
         !-------------------------------------------------------------------------------
#ifdef MMech
         IF(HasQM()) THEN
#endif
         CALL Delete(BSiz)
         CALL Delete(OffS)
#ifdef MMech
         ENDIF
#endif
#ifdef PARALLEL
         IF(InParallel)THEN
            CALL Delete(Beg)
            CALL Delete(End)
         ENDIF
#endif
         IF(PrintFlags%Key>=DEBUG_MEDIUM)THEN
!        IF(PrintFlags%Key>=DEBUG_MEDIUM.OR.PrintFlags%Chk==DEBUG_CHKSUMS)THEN
            CALL Elapsed_TIME(PerfMon,'Accum',Proc_O=Prog)
            CALL PPrint(PerfMon,Prog)
         ENDIF
         IF(PrintFlags%Key==DEBUG_MAXIMUM) &
            CALL PPrint(MemStats,Prog)

#ifdef PARALLEL_CLONES
         CALL CloseHDFGroup(H5GroupID)
         HDF_CurrentID=HDFFileID
#endif
#ifdef PARALLEL
         IF(InParallel) &
            CALL FiniMPI()
         IF(PrintFlags%Key>DEBUG_MEDIUM.AND.MyId==ROOT)  &
            CALL TimeStamp('Exiting '//TRIM(Prog),Enter_O=.FALSE.)
#else
         IF(PrintFlags%Key>DEBUG_MEDIUM)  &
            CALL TimeStamp('Exiting '//TRIM(Prog),Enter_O=.FALSE.)

#endif 
         CALL MarkSuccess()
!        Print time stamp
         CALL CloseHDF(HDFFileID)
         STOP 
      END SUBROUTINE ShutDown     
      !=========================================================
      ! MARK FAILURE OF PROG
      !=========================================================
      SUBROUTINE MarkFailure(Prog)
         CHARACTER(LEN=*) :: Prog
         CALL Put(.TRUE.,'ProgramFailed')
         CALL Put(Prog  ,'FailedProgram')
      END SUBROUTINE         
      !=========================================================
      ! MARK SUCCESS OF PROG
      !=========================================================
      SUBROUTINE MarkSuccess()
         CALL Put(.FALSE.,'ProgramFailed')
      END SUBROUTINE         
      !=========================================================
      ! LOAD MISC GLOBAL VARIABLES
      !=========================================================
      SUBROUTINE Init_Globals(Args)
        TYPE(ARGMT) :: Args
        INTEGER :: I,ChkNPrc
        !--------------------------------------------------------!
         CALL Get(logFile,'logfile')
         CALL Get(OutFile,'outputfile')
         CALL Get(InpFile,'inputfile')
         ScrName=TRIM(MONDO_SCRATCH)//TRIM(Args%C%C(1))
!         CALL BCast(ScrName)
!        Load global scalars
#ifdef MMech
         IF(HasQM())THEN
#endif
            CALL Get(NEl,      'nel',           Tag_O=CurGeom)
            CALL Get(NAtoms,   'natoms',        Tag_O=CurGeom)
            CALL Get(MaxAtms,  'maxatms',       Tag_O=CurBase)
            CALL Get(MaxBlks,  'maxblks',       Tag_O=CurBase)
            CALL Get(MaxNon0,  'maxnon0',       Tag_O=CurBase)
            CALL Get(NBasF,    'nbasf',         Tag_O=CurBase)
            CALL Get(ModelChem,'ModelChemistry',Tag_O=CurBase)
            CALL New(BSiz,NAtoms)
            CALL New(OffS,NAtoms)
            CALL Get(BSiz,'atsiz',Tag_O=CurBase)
            CALL Get(OffS,'atoff',Tag_O=CurBase)
            ! Global value for max block size
            MaxBlkSize=0
            DO I=1,NAtoms 
               MaxBlkSize=MAX(MaxBlkSize,BSiz%I(I)) 
            ENDDO
#ifdef PARALLEL
            IF(InParallel)THEN
               CALL Get(MaxAtmsNode,'maxatmsnode', Tag_O=CurBase)
               CALL Get(MaxBlksNode,'maxblksnode', Tag_O=CurBase)
               CALL Get(MaxNon0Node,'maxnon0node', Tag_O=CurBase)
               CALL New(OffSt,NPrc-1,0)
               CALL Get(OffSt,'dbcsroffsets',Tag_O=CurBase)
               CALL Get(ChkNPrc,'chknprc')
               IF(NPrc/=ChkNPrc) &
                    CALL Halt(' In StartUp() --- Inconsistency: NPrc = '  &
                    //TRIM(IntToChar(NPrc))//' ChkNPrc = ' &
                    //TRIM(IntToChar(ChkNPrc)))
               CALL New(Beg,NPrc-1,0)
               CALL New(End,NPrc-1,0)
               CALL Get(Beg,'beg',Tag_O=CurBase)
               CALL Get(End,'end',Tag_O=CurBase)
            ENDIF
#endif
#ifdef MMech
         ENDIF
#endif
       END SUBROUTINE Init_Globals
      !=========================================================
      ! LOAD GLOBAL STRINGS FOR SCF STATE INFORMATION
      ! Clumsy, should ultimately be done away with....
      !=========================================================
      SUBROUTINE Init_SCFStat(Args)
         TYPE(ARGMT) :: Args
        !------------------------------------------------------!
        IF(Args%NC>=2) &
             SCFActn=TRIM(Args%C%C(2))
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
        CurBase=TRIM(IntToChar(Args%I%I(2)))
        PrvBase=TRIM(IntToChar(Previous(2)))
      END SUBROUTINE Init_SCFStat
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
