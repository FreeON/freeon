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
!-------------------------------------------------------------------------------
!  Global status variables
   INTEGER,DIMENSION(3) :: Current
   INTEGER,DIMENSION(3) :: Previous
   CHARACTER(LEN=3)     :: SCFCycl
   CHARACTER(LEN=3)     :: PrvCycl
   CHARACTER(LEN=3)     :: NxtCycl
   CHARACTER(LEN=20)    :: SCFActn
   CHARACTER(LEN=3)     :: CurBase
   CHARACTER(LEN=3)     :: PrvBase
   CHARACTER(LEN=3)     :: CurGeom
   CHARACTER(LEN=3)     :: PrvGeom
   CHARACTER(LEN=3)     :: NxtGeom
!-------------------------------------------------------------------------------
   CONTAINS
      SUBROUTINE StartUp(Args,Prog,Serial_O)
         TYPE(ARGMT),INTENT(OUT)              :: Args
         CHARACTER(LEN=*),INTENT(IN)          :: Prog
         LOGICAL,OPTIONAL                     :: Serial_O
         INTEGER                              :: I
#ifdef PARALLEL
         LOGICAL                              :: Serial
         INTEGER                              :: ChkNPrc
!-------------------------------------------------------------------------------
         IF(PRESENT(Serial_O))THEN
            Serial=Serial_O
         ELSE
            Serial=.TRUE.
         ENDIF
!        Fire up MPI
         IF(.NOT.Serial)CALL InitMPI()
#endif
!        Get arguments and open InfFile 
         CALL Get(Args)
         IF(Args%NC>=2) &
            SCFActn=TRIM(Args%C%C(2))
         Current=Args%I%I(1:3)
         Previous=Args%I%I(4:6)
!        SCF Cycle
         SCFCycl=TRIM(IntToChar(Args%I%I(1)))
         IF(Args%I%I(3)-1<1)THEN
            PrvCycl=TRIM(IntToChar(Previous(1)))        
         ELSE
            PrvCycl=TRIM(IntToChar(Args%I%I(1)-1))
         ENDIF
         NxtCycl=TRIM(IntToChar(Args%I%I(1)+1))
!        Geometry
         CurGeom=TRIM(IntToChar(Args%I%I(3)))
         IF(Args%I%I(3)-1<1)THEN
            PrvGeom=TRIM(IntToChar(Previous(3)))        
         ELSE
            PrvGeom=TRIM(IntToChar(Args%I%I(3)-1))
         ENDIF
         NxtGeom=TRIM(IntToChar(Args%I%I(3)+1))
!        Basis
         CurBase=TRIM(IntToChar(Args%I%I(2)))
         PrvBase=TRIM(IntToChar(Previous(2)))
!-----------------------------------------------------------------------------
!        Get PWD and SCRATCH directories from env
         CALL GetEnv('PWD',MONDO_PWD)   
         CALL GetEnv('MONDO_SCRATCH',MONDO_SCRATCH)
         IF(LEN(TRIM(MONDO_HOME))==0)CALL Halt(' $(MONDO_HOME) not set.')
         IF(LEN(TRIM(MONDO_SCRATCH))==0)CALL Halt(' $(MONDO_SCRATCH) not set.')
!        Set current working file names
         MONDO_PWD=TRIM(MONDO_PWD)//'/'
         MONDO_SCRATCH=TRIM(MONDO_SCRATCH)//'/'
         ScrName=TRIM(MONDO_SCRATCH)//TRIM(Args%C%C(1))
         PWDName=TRIM(MONDO_PWD)//TRIM(Args%C%C(1))
         InfFile=TRIM(ScrName)//TRIM(InfF)
!-----------------------------------------------------------------------------
!        Open HDF file and mark for failure
         CALL OpenHDF(TRIM(InfFile))
         CALL MarkFailure(Prog)
!-----------------------------------------------------------------------------
!        Load global characters
         CALL Get(logFile,'logfile')
         CALL Get(OutFile,'outputfile')
         CALL Get(InpFile,'inputfile')
!        Load global scalars
         CALL Get(NEl,      'nel')
         CALL Get(NAtoms,   'natoms')
         CALL Get(MaxAtms,  'maxatms',       Tag_O=CurBase)
         CALL Get(MaxBlks,  'maxblks',       Tag_O=CurBase)
         CALL Get(MaxNon0,  'maxnon0',       Tag_O=CurBase)
         CALL Get(NBasF,    'nbasf',         Tag_O=CurBase)
         CALL Get(ModelChem,'ModelChemistry',Tag_O=CurBase)
#ifdef PARALLEL
         CALL Get(MaxAtmsNode,'maxatmsnode', Tag_O=CurBase)
         CALL Get(MaxBlksNode,'maxblksnode', Tag_O=CurBase)
         CALL Get(MaxNon0Node,'maxnon0node', Tag_O=CurBase)
#endif
!        Initialize global objects
         CALL Init(MemStats)
         CALL Init(PrintFlags,Prog)
         CALL Elapsed_Time(PerfMon,'Init')
!        Load global objects
         CALL New(BSiz,NAtoms)
         CALL New(OffS,NAtoms)
         CALL Get(BSiz,'atsiz',Tag_O=CurBase)
         CALL Get(OffS,'atoff',Tag_O=CurBase)
!        Load global value for max block size
         MaxBlkSize=0
         DO I=1,NAtoms; MaxBlkSize=MAX(MaxBlkSize,BSiz%I(I)); ENDDO
!        Load global thresholding values
         CALL SetThresholds(CurBase)
#ifdef PARALLEL
         IF(InParallel)THEN        
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
         CALL Delete(BSiz)
         CALL Delete(OffS)
#ifdef PARALLEL
         IF(InParallel)THEN
            CALL Delete(Beg)
            CALL Delete(End)
         ENDIF
#endif
         IF(PrintFlags%Key>=DEBUG_MEDIUM)THEN
            CALL Elapsed_TIME(PerfMon,'Accum',Proc_O=Prog)
            CALL PPrint(PerfMon,Prog)
         ENDIF
         IF(PrintFlags%Key==DEBUG_MAXIMUM) &
            CALL PPrint(MemStats,Prog)
         CALL MarkSuccess()
         CALL CloseHDF()
!        Print time stamp
#ifdef PARALLEL
         IF(InParallel) &
            CALL FiniMPI()
         IF(PrintFlags%Key>DEBUG_MEDIUM.AND.MyId==ROOT)  &
            CALL TimeStamp('Exiting '//TRIM(Prog),Enter_O=.FALSE.)
#else
         IF(PrintFlags%Key>DEBUG_MEDIUM)  &
            CALL TimeStamp('Exiting '//TRIM(Prog),Enter_O=.FALSE.)
#endif 
         STOP 
      END SUBROUTINE ShutDown     

      SUBROUTINE MarkFailure(Prog)
         CHARACTER(LEN=*) :: Prog
         CALL Put(.TRUE.,'ProgramFailed')
         CALL Put(Prog  ,'FailedProgram')
      END SUBROUTINE         
 
      SUBROUTINE MarkSuccess()
         CALL Put(.FALSE.,'ProgramFailed')
      END SUBROUTINE         
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
