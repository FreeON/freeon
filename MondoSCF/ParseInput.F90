!    PARSING MODULE
!    Authors: Matt Challacombe, C.J. Tymczak
!----------------------------------------------------------------------------
MODULE ParseInput
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalObjects
   USE GlobalCharacters
   USE ProcessControl
   USE BasisSetParameters
   USE InOut
   USE MemMan
   USE Macros
   USE ParsingConstants
   USE ParsingKeys
   USE Parse
   USE Order
   USE SCFLocals
   USE PrettyPrint
   USE Functionals
   USE Overlay
   USE AtomPairs
#ifdef MMech
   USE Dynamo    
   USE IntCoo    
#endif
   IMPLICIT NONE
   INTEGER                    :: I,ILoc,NLoc,NOpts
   INTEGER,DIMENSION(MaxSets) :: Loc
   CHARACTER(LEN=3)           :: CSet
   CHARACTER(LEN=8)           :: CCrd
   CONTAINS
!============================================================================
!     Parse the input and create proto-HDF file
!============================================================================
      SUBROUTINE ParseInp(Ctrl)
         Implicit None
         CHARACTER (LEN=DEFAULT_CHR_LEN) :: OPLS_DATABASE,MM_COORDS,MM_SEQ
!
         TYPE(SCFControls), INTENT(INOUT) :: Ctrl
!
!        Read comand line, environement variables, create file names, init files etc
         CALL ParseCmndLine(Ctrl)
#ifdef MMech
         CALL ParseMech(Ctrl)
#endif
         CALL ParseGrad(Ctrl)
#ifdef MMech
         If(Ctrl%Mechanics(2)) Then
#endif
!        Read geometry and lattice variables, reorder, rescale etc  
         CALL ParseGeometry(Ctrl)  
!        Read in the basis sets
         CALL ParseBaseSets(Ctrl) 
!        Read in the SCF options
         CALL ParseMethods(Ctrl)
!   Call ParseAcc(Ctrl)
#ifdef MMech
         EndIf
#endif
!
#ifdef MMech
         If(Ctrl%Mechanics(1)) Then
           Call ParseMM(Ctrl)
         EndIf
#endif
!
! Parse accuracy
!
         Call ParseAcc(Ctrl)
!
!        Print Parsed options
         CALL ParsePrint(Ctrl)
!
      END SUBROUTINE ParseInp
!============================================================================
!     Parce The Command Lines
!============================================================================
      SUBROUTINE ParseCmndLine(Ctrl)
         TYPE(SCFControls)          :: Ctrl
         TYPE(ARGMT)                    :: Args
         INTEGER                        :: I,K,L,DotDex
         INTEGER,EXTERNAL               :: GetPID
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Line,MPILine,GenFile,OldInfo ,Mssg        
         CHARACTER(LEN=DEFAULT_CHR_LEN), &
                   PARAMETER            :: Soft='-s'
         CHARACTER(LEN=DEFAULT_CHR_LEN), &
                   DIMENSION(3)         :: SoftLinks
         CHARACTER(LEN=1)               :: S=CHAR(092)
         TYPE(CHR_VECT)                 :: TmpChars
!----------------------------------------------------------------------------
!        Get command line arguments and environment variables
!
         CALL Get(Args)
         IF(Args%NC==0)CALL MondoHalt(PRSE_ERROR,' No arguments to MondoSCF !')
!
         CALL GetEnv('PWD',MONDO_PWD)
         CALL GetEnv('MONDO_HOME',MONDO_HOME)
         CALL GetEnv('MONDO_SCRATCH',MONDO_SCRATCH)
         CALL GetEnv('MONDO_EXEC',MONDO_EXEC)
         CALL GetEnv('MONDO_HOST',MONDO_HOST)
         CALL GetEnv('MONDO_MACH',MONDO_MACH)
         CALL GetEnv('MONDO_SYST',MONDO_SYST)
         CALL GetEnv('MONDO_VRSN',MONDO_VRSN)
         CALL GetEnv('MONDO_PLAT',MONDO_PLAT)
         IF(LEN(TRIM(MONDO_HOME))==0)CALL MondoHalt(PRSE_ERROR,' $(MONDO_HOME) not set.')
         IF(LEN(TRIM(MONDO_SCRATCH))==0)CALL MondoHalt(PRSE_ERROR,' $(MONDO_SCRATCH) not set.')
!----------------------------------------------------------------------------
!        Set local path names etc
!
         MONDO_PWD=TRIM(MONDO_PWD)//'/'
         MONDO_HOME=TRIM(MONDO_HOME)//'/'
         MONDO_SCRATCH=TRIM(MONDO_SCRATCH)//'/'
         MONDO_EXEC=TRIM(MONDO_EXEC)//'/'
         PROCESS_ID=IntToChar(GetPID())
!----------------------------------------------------------------------------
!        Determine input and output names with full paths        
         InpFile=TRIM(MONDO_PWD)//TRIM(Args%C%C(1))
!        Determine if input file has a '.' in it.  If so, create SCF name 
!        from string up to the '.'
         DotDex=INDEX(Args%C%C(1),'.')
         IF(DotDex==0)THEN
            SCF_NAME=TRIM(Args%C%C(1))//'_'//TRIM(PROCESS_ID)
         ELSE
            SCF_NAME=Args%C%C(1)(1:DotDex-1)//'_'//TRIM(PROCESS_ID)
         ENDIF
         ScrName=TRIM(MONDO_SCRATCH)//TRIM(SCF_NAME)
         PWDName=TRIM(MONDO_PWD)//TRIM(SCF_NAME)
!        Create user defined or implicit file names
         IF(Args%NC==1)THEN
            OutFile=TRIM(PWDName)//OutF 
            LogFile=TRIM(PWDName)//LogF
            GeoFile=TRIM(PWDName)//GeoF
         ELSEIF(Args%NC==2)THEN
            OutFile=TRIM(MONDO_PWD)//TRIM(Args%C%C(2))
            LogFile=TRIM(PWDName)//LogF
            GeoFile=TRIM(PWDName)//GeoF
         ELSEIF(Args%NC==3)THEN
            OutFile=TRIM(MONDO_PWD)//TRIM(Args%C%C(2))
            LogFile=TRIM(MONDO_PWD)//TRIM(Args%C%C(3))
            GeoFile=TRIM(PWDName)//GeoF
         ELSEIF(Args%NC==4)THEN
            OutFile=TRIM(MONDO_PWD)//TRIM(Args%C%C(2))
            LogFile=TRIM(MONDO_PWD)//TRIM(Args%C%C(3))
            GeoFile=TRIM(MONDO_PWD)//TRIM(Args%C%C(4))
         ENDIF
!        Name the HDF5 info file
         InfFile=TRIM(ScrName)//InfF                
!        Set SCF and InfoFile names
         Ctrl%Info=InfFile
         Ctrl%Name=TRIM(SCF_NAME)
         CALL InitHDF(Ctrl%Info)
         CALL OpenHDF(Ctrl%Info)
!----------------------------------------------------------------------------
!        Check for guess, restart. Initialize HDF file
!
         CALL OpenASCII(InpFile,Inp,OldFileQ_O=.TRUE.)
         IF(OptKeyQ(Inp,GUESS_OPTION,GUESS_RESTART))THEN
            Ctrl%Rest=.TRUE.            
            Ctrl%SuperP=.FALSE.
            IF(.NOT.OptCharQ(Inp,RESTART_INFO,Ctrl%OldInfo))  &
               CALL Halt('Restart requested, but no hdf file specified.')
!           Check for absolute path 
            IF(INDEX(Ctrl%OldInfo,'/')==0)THEN
               Ctrl%OldInfo=TRIM(MONDO_PWD)//Ctrl%OldInfo
            ENDIF
            CALL Put(Ctrl%OldInfo,'oldinfo')
         ELSE
            Ctrl%SuperP=.TRUE.
            Ctrl%Rest=.FALSE.
         ENDIF
         CLOSE(UNIT=Inp,STATUS='KEEP')
!----------------------------------------------------------------------------
!        Initialize and open info file 
!
         CALL Put(InpFile,'inputfile')
         CALL Put(InfFile,'infofile')
         CALL Put(LogFile,'logfile')
         CALL Put(OutFile,'outputfile')
!        Initialize ASCII files
         CALL OpenASCII(OutFile,Out,NewFile_O=.TRUE.)
         CALL OpenASCII(LogFile,LgF,NewFile_O=.TRUE.)
         CLOSE(UNIT=Out,STATUS='KEEP')
         CLOSE(UNIT=LgF,STATUS='KEEP')
!----------------------------------------------------------------------------
!        Write banner and title to output file
!
         CALL OpenASCII(InpFile,Inp)
         CALL OpenASCII(OutFile,Out)
!         CALL PrintProtectL(Out)
         WRITE(Out,77)(Rtrn,I=1,15)
         WRITE(*,77)(Rtrn,I=1,15)
      77 FORMAT(A1,A1,                                                   &
         ' __    __                 _       ____________ ______ ',A1,    &
         '|  \  /  |               | |     /       /    |      |',A1,    & 
         "|   \/   | ___  _ __   __| | ___/   ____/   __|  |==='",A1,    &
         "|        |/ _ \| '_ \ / _  |/ _ \____  \   (__|  ____|",A1,    &
         '|  |\/|  | (_) | | | | (_| | (_) )     /\     |  |    ',A1,    &
         '|__|  |__|\___/|_| |_|\____|\___/_____/  \____|__|    ',A1,A1, &
         ' Version 1.0 alpha 4                                  ',A1,    &  
         ' A program suite for O(N) SCF theory and ab initio MD ',A1,    &
         ' Matt Challacombe, Eric Schwegler, C.J. Tymczak,      ',A1,    &
         ' Chee Kwan Gan, Karoly Nemeth and Anders Niklasson    ',A1,    &
         ' Los Alamos National Laboratory (LA-CC 01-2)          ',A1,    & 
         ' Copywrite 2001, University of California.            ',A1)
!        Write information on host, platform, etc
         Mssg='Compliled for '//TRIM(MONDO_PLAT)//', executing on '//TRIM(MONDO_HOST) &
           //Rtrn//' a '//TRIM(MONDO_MACH)//' machine'//' running '//TRIM(MONDO_SYST) &
           //' '//TRIM(MONDO_VRSN)       
         WRITE(*,*)TRIM(Mssg)
         WRITE(Out,*)TRIM(Mssg)
         WRITE(*,*)
         WRITE(Out,*)
!        Parse for title, put to output file
         CALL Align(BEGIN_TITLE,Inp)
         K=1
         DO I=1,1000
            READ(Inp,DEFAULT_CHR_FMT,END=1)Line
            IF(INDEX(Line,END_TITLE)/=0)GOTO 2
            WRITE(Out,DEFAULT_CHR_FMT)Line
         ENDDO
      1  CALL MondoHalt(PRSE_ERROR,' Found no <EndTitle> in input file '//TRIM(InpFile))
      2  CONTINUE
         WRITE(Out,*)' '
!         CALL PrintProtectR(Out)
!        Close input and output
         CLOSE(UNIT=Inp,STATUS='KEEP')
         CLOSE(UNIT=Out,STATUS='KEEP')
!        TimeStamp
         CALL TimeStamp('Starting MondoSCF')
#if defined(PARALLEL) && !defined(MPI2)
!----------------------------------------------------------------------------
!        Parse <OPTIONS> for mpirun invokation
 
         CALL OpenASCII(InpFile,Inp)
         IF(.NOT.OptCharQ(Inp,MPI_OPTION,MPILine))          & 
            CALL Halt(' mpi invokation command not found. ' &
                    //' Check input for option '            &
                    //TRIM(MPI_OPTION))
         CLOSE(UNIT=Inp,STATUS='KEEP')
         CALL LineToChars(MPILine,TmpChars)
         MPI_INVOKE=TmpChars%C(1)
         MPI_FLAGS=' '
         DO I=2,SIZE(TmpChars%C)
            MPI_FLAGS=TRIM(MPI_FLAGS)//Blnk//TRIM(TmpChars%C(I)) 
         ENDDO
         CALL Delete(TmpChars)           
#endif
!----------------------------------------------------------------------------
!        Parse <OPTIONS> file for print flags (uses Init in Macros...)
!
         CALL Init(PrintFlags,'MAIN')
         CALL OpenASCII(InpFile,Inp)
         IF(OptKeyQ(Inp,GLOBAL_DEBUG,DBG_PRT_INTS))THEN
            PrintFlags%Int=DEBUG_INTEGRAL
         ELSE
            PrintFlags%Int=DEBUG_NONE
         ENDIF                       
         IF(OptKeyQ(Inp,GLOBAL_DEBUG,DBG_PRT_SETS))THEN
            PrintFlags%Set=DEBUG_BASISSET
         ELSE
            PrintFlags%Set=DEBUG_NONE
         ENDIF                       
         CLOSE(UNIT=Inp,STATUS='KEEP')
!----------------------------------------------------------------------------
!        Tidy up 
!
         CALL CloseHDF()
         CALL Delete(Args)
!
      END SUBROUTINE ParseCmndLine
!
!============================================================================
!     Print Out the Parsed Information
!============================================================================
      SUBROUTINE ParsePrint(Ctrl)
        TYPE(SCFControls)                :: Ctrl
        INTEGER                          :: I,RestAccL,RestMeth,RestModl
        TYPE(INT_VECT)                   :: Stat
        TYPE(CRDS)                       :: GM
        CHARACTER(LEN=8)                 :: Cur
        CHARACTER(LEN=DEFAULT_CHR_LEN)   :: Method,Accuracy,Chemistry
        CHARACTER(LEN=BASESET_CHR_LEN)   :: BName   
        CHARACTER(LEN=2*DEFAULT_CHR_LEN) :: Mssg
!----------------------------------------------------------------------------
        ! Open input file
        CALL OpenASCII(OutFile,Out)
        CALL PrintProtectL(Out)
        WRITE(Out,*)'Current HDF file :: ',TRIM(Ctrl%Info)
        IF(Ctrl%Rest)THEN
           WRITE(Out,*)'Restart HDF file :: ',TRIM(Ctrl%OldInfo)
           CALL OpenHDF(Ctrl%OldInfo)
           CALL New(Stat,3)
           CALL Get(Stat,'current')
           Cur=IntToChar(Stat%I(2))
           CALL Get(RestAccL,'SCFAccuracy',Cur)
           CALL Get(RestMeth,'SCFMethod',Cur)
           CALL Get(RestModl,'ModelChemistry',Cur)
           CALL Get(BName,'bsetname',Cur)
           Cur=IntToChar(Stat%I(3)) 
           CALL Get(GM,Cur)
           CALL CloseHDF()
           Mssg='Restart using '//TRIM(BName)//'/'//TRIM(FunctionalName(RestModl)) &
                //' density from previous geometry #'//TRIM(Cur)
           WRITE(Out,*)TRIM(Mssg)
           CALL Delete(Stat)
           CALL Delete(GM)
        ENDIF
        WRITE(Out,*)
        ! Print the accuracy, method and model chemistry for each basis set
        WRITE(Out,*)'PROGRAM OF CALCULATIONS:',Rtrn
        DO I=1,Ctrl%NSet
#ifdef MMech
           IF(HasQM()) THEN
#endif
              IF(Ctrl%Method(I)==SDMM_R_SCF)THEN
                 Method='restricted simplified density matrix minimization using '
              ELSE
                 Method='restricted Roothaan-Hall solution of the SCF using '
              ENDIF
#ifdef MMech
           ENDIF
#endif
           IF(Ctrl%AccL(I)==1)THEN
              Accuracy='  a loose accuracy level, '
           ELSEIF(Ctrl%AccL(I)==2)THEN
              Accuracy='  a good accuracy level, '
           ELSEIF(Ctrl%AccL(I)==3)THEN
              Accuracy='  a tight accuracy level, '
           ELSEIF(Ctrl%AccL(I)==4)THEN
              Accuracy='  a very tight accuracy level, '      
           ENDIF
#ifdef MMech
           IF(HasQM()) THEN
#endif
              Chemistry=' and the '//TRIM(FunctionalName(Ctrl%Model(I)))//' model.'
              IF(I==1)THEN
                 Mssg=' A '//TRIM(Method)//Rtrn//TRIM(Accuracy)//' '//TRIM(Ctrl%BName(1))//TRIM(Chemistry)
              ELSE
                 Mssg=' Followed by a  '//TRIM(Method)//Rtrn//TRIM(Accuracy) &
                      //TRIM(Ctrl%BName(I))//TRIM(Chemistry)
              ENDIF
              WRITE(Out,*)TRIM(Mssg),Rtrn
#ifdef MMech
           ENDIF
#endif
        ENDDO
        CALL PrintProtectR(Out)
        CLOSE(UNIT=Out,STATUS='KEEP')
        ! Put SCF Method
#ifdef MMech
        IF(HasQM()) THEN
#endif
           CALL OpenHDF(InfFile)
           DO I=1,Ctrl%NSet
              CALL Put(Ctrl%AccL(I),'SCFAccuracy',Tag_O=IntToChar(I))
              CALL Put(Ctrl%Method(I),'SCFMethod',Tag_O=IntToChar(I))
           ENDDO
           CALL CloseHDF()
#ifdef MMech
        ENDIF
#endif
      END SUBROUTINE ParsePrint
!============================================================================
!     Parse The Methods
!============================================================================
      SUBROUTINE ParseMethods(Ctrl)
         TYPE(SCFControls)          :: Ctrl
         TYPE(TOLS)                 :: Thrsh ! Thresholds
!----------------------------------------------------------------------------
!        Open info file
         CALL OpenHDF(InfFile)
!        Open input file
         CALL OpenASCII(InpFile,Inp) 
!
!        Parse <OPTIONS.SCF>
!
         IF(OptKeyQ(Inp,INKFOCK_OPTION,INKFOCK_ON))THEN
            Ctrl%ShudInk=.TRUE.
         ELSE
            Ctrl%ShudInk=.FALSE.
         ENDIF       
!
         NOpts=0
         Ctrl%Method=RH_R_SCF ! default is restricted RH
!
         IF(OptKeyLocQ(Inp,SCF_OPTION,SCF_SDMM,MaxSets,NLoc,Loc))THEN
            NOpts=NOpts+NLoc
            DO ILoc=1,NLoc; Ctrl%Method(Loc(ILoc))=SDMM_R_SCF; ENDDO             
         ENDIF
!
         IF(OptKeyLocQ(Inp,SCF_OPTION,SCF_PM,MaxSets,NLoc,Loc))THEN
            NOpts=NOpts+NLoc
            DO ILoc=1,NLoc; Ctrl%Method(Loc(ILoc))=PM_R_SCF; ENDDO             
         ENDIF         
!
         IF(OptKeyLocQ(Inp,SCF_OPTION,SCF_SP2,MaxSets,NLoc,Loc))THEN
            NOpts=NOpts+NLoc
            DO ILoc=1,NLoc; Ctrl%Method(Loc(ILoc))=SP2_R_SCF; ENDDO             
         ENDIF
!
         IF(OptKeyLocQ(Inp,SCF_OPTION,SCF_SP4,MaxSets,NLoc,Loc))THEN
            NOpts=NOpts+NLoc
            DO ILoc=1,NLoc; Ctrl%Method(Loc(ILoc))=SP4_R_SCF; ENDDO             
         ENDIF
!
         IF(OptKeyLocQ(Inp,SCF_OPTION,SCF_NT4,MaxSets,NLoc,Loc))THEN
            NOpts=NOpts+NLoc
            DO ILoc=1,NLoc; Ctrl%Method(Loc(ILoc))=NT4_R_SCF; ENDDO             
         ENDIF
!
         IF(OptKeyLocQ(Inp,SCF_OPTION,SCF_RHHF,MaxSets,NLoc,Loc))THEN
            NOpts=NOpts+NLoc
            DO ILoc=1,NLoc; Ctrl%Method(Loc(ILoc))=RH_R_SCF; ENDDO
         ENDIF 
         IF(NOpts==1)THEN
            Ctrl%Method=Ctrl%Method(1)
         ELSEIF(NOpts>1.AND.NOpts/=Ctrl%NSet)THEN
            CALL MondoHalt(PRSE_ERROR,'Number of '//SCF_OPTION &
                           //' options does not match number of Basis sets.')
         ENDIF
!----------------------------------------------------------------------------
!        Parse <OPTIONS.MODEL> 
!
         NOpts=0
         Ctrl%Model=1 ! Default is HF
!        Exact exchage
         IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_ExactX,MaxSets,NLoc,Loc))THEN
            NOpts=NOpts+NLoc
            DO ILoc=1,NLoc; Ctrl%Model(Loc(ILoc))=EXACT_EXCHANGE; ENDDO
         ENDIF 
!        Slater-Dirac exchage
         IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_SD,MaxSets,NLoc,Loc))THEN
            NOpts=NOpts+NLoc
            DO ILoc=1,NLoc; Ctrl%Model(Loc(ILoc))=SD_EXCHANGE; ENDDO
         ENDIF 
!        X-alpha exchage
         IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_XA,MaxSets,NLoc,Loc))THEN
            NOpts=NOpts+NLoc
            DO ILoc=1,NLoc; Ctrl%Model(Loc(ILoc))=XA_EXCHANGE; ENDDO
         ENDIF 
!        PW91 exchange
         IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_PW91x,MaxSets,NLoc,Loc))THEN
            NOpts=NOpts+NLoc
            DO ILoc=1,NLoc; Ctrl%Model(Loc(ILoc))=PW91_EXCHANGE; ENDDO
         ENDIF       
!        PBE exchange
         IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_PBEx,MaxSets,NLoc,Loc))THEN
             NOpts=NOpts+NLoc
             DO ILoc=1,NLoc; Ctrl%Model(Loc(ILoc))=PBE_EXCHANGE; ENDDO
         ENDIF       
!        B88 exchange
         IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_B88x,MaxSets,NLoc,Loc))THEN
             NOpts=NOpts+NLoc
             DO ILoc=1,NLoc; Ctrl%Model(Loc(ILoc))=B88_EXCHANGE; ENDDO
         ENDIF       
!        Pure Slater exchange with VWN3 LSDA correlation 
         IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_VWN3xc,MaxSets,NLoc,Loc))THEN
             NOpts=NOpts+NLoc
             DO ILoc=1,NLoc; Ctrl%Model(Loc(ILoc))=PURE_VWN3_LSD; ENDDO
         ENDIF       
!        Pure Slater exchange with VWN5 LSDA correlation 
         IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_VWN5xc,MaxSets,NLoc,Loc))THEN
             NOpts=NOpts+NLoc
             DO ILoc=1,NLoc; Ctrl%Model(Loc(ILoc))=PURE_VWN5_LSD; ENDDO
         ENDIF       
!        Pure Slater exchange with PW91 LSDA correlation 
         IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_PW91xc,MaxSets,NLoc,Loc))THEN
             NOpts=NOpts+NLoc
             DO ILoc=1,NLoc; Ctrl%Model(Loc(ILoc))=PURE_PW91_LSD; ENDDO
         ENDIF       
!        Pure PBE GGA exchange-correlation 
         IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_PBExc,MaxSets,NLoc,Loc))THEN
             NOpts=NOpts+NLoc
             DO ILoc=1,NLoc; Ctrl%Model(Loc(ILoc))=PURE_PBE_GGA; ENDDO
         ENDIF
!        Pure BLYP GGA exchange-correlation 
         IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_BLYPxc,MaxSets,NLoc,Loc))THEN
             NOpts=NOpts+NLoc
             DO ILoc=1,NLoc; Ctrl%Model(Loc(ILoc))=PURE_BLYP_GGA; ENDDO
         ENDIF       
!        Hybrid B3LYP/VWN3 exchange-correlation 
         IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_B3LYP_VWN3xc,MaxSets,NLoc,Loc))THEN
             NOpts=NOpts+NLoc
             DO ILoc=1,NLoc; Ctrl%Model(Loc(ILoc))=HYBRID_B3LYP_VWN3; ENDDO
         ENDIF       
!        Hybrid B3LYP/VWN5 exchange-correlation 
         IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_B3LYP_VWN5xc,MaxSets,NLoc,Loc))THEN
             NOpts=NOpts+NLoc
             DO ILoc=1,NLoc; Ctrl%Model(Loc(ILoc))=HYBRID_B3LYP_VWN5; ENDDO
         ENDIF       
!        Hybrid B3LYP/PW91 exchange-correlation 
         IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_B3LYP_PW91xc,MaxSets,NLoc,Loc))THEN
             NOpts=NOpts+NLoc
             DO ILoc=1,NLoc; Ctrl%Model(Loc(ILoc))=HYBRID_B3LYP_PW91; ENDDO
         ENDIF       
!        Hybrid PBE0 exchange-correlation 
         IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_PBE0xc,MaxSets,NLoc,Loc))THEN
             NOpts=NOpts+NLoc
             DO ILoc=1,NLoc; Ctrl%Model(Loc(ILoc))=HYBRID_PBE0; ENDDO
         ENDIF       
         IF(NOpts==1)THEN
            Ctrl%Model=Ctrl%Model(1)
         ELSEIF(NOpts>1.AND.NOpts/=Ctrl%NSet)THEN
            CALL MondoHalt(PRSE_ERROR,'Number of '//MODEL_OPTION &
                           //' options does not match number of Basis sets.')
         ENDIF
!        Put model chemistry
         DO I=1,Ctrl%NSet
            CALL Put(Ctrl%Model(I),'ModelChemistry',Tag_O=IntToChar(I))
         ENDDO
!!----------------------------------------------------------------------------
!!        Parse <OPTIONS.ACCURACY> 
!!
!         Ctrl%AccL=2 ! Default is "good"    
!         NOpts=0
!!
!         IF(OptKeyLocQ(Inp,ACCURACY_OPTION,ACCURACY_CHEEZY,MaxSets,NLoc,Loc))THEN
!            NOpts=NOpts+NLoc
!            DO ILoc=1,NLoc; Ctrl%AccL(Loc(ILoc))=1; ENDDO
!         ENDIF
!         IF(OptKeyLocQ(Inp,ACCURACY_OPTION,ACCURACY_GOOD,MaxSets,NLoc,Loc))THEN
!            NOpts=NOpts+NLoc
!            DO ILoc=1,NLoc; Ctrl%AccL(Loc(ILoc))=2; ENDDO
!         ENDIF
!         IF(OptKeyLocQ(Inp,ACCURACY_OPTION,ACCURACY_TIGHT,MaxSets,NLoc,Loc))THEN
!            NOpts=NOpts+NLoc
!            DO ILoc=1,NLoc; Ctrl%AccL(Loc(ILoc))=3; ENDDO
!         ENDIF
!         IF(OptKeyLocQ(Inp,ACCURACY_OPTION,ACCURACY_RETENTIVE,MaxSets,NLoc,Loc))THEN
!            NOpts=NOpts+NLoc
!            DO ILoc=1,NLoc; Ctrl%AccL(Loc(ILoc))=4; ENDDO
!         ENDIF
!         IF(NOpts>1.AND.NOpts/=Ctrl%NSet) &
!            CALL MondoHalt(PRSE_ERROR,'Number of '//ACCURACY_OPTION &
!                           //' options does not match number of Basis sets.')
!!        Set thresholds
!         IF(NOpts==Ctrl%NSet)THEN
!!           All thresholds are user determined
!            write(*,*) 'thrsh 1'
!            DO I=1,Ctrl%NSet
!               Thrsh%Cube=CubeNeglect(Ctrl%AccL(I))
!               Thrsh%Trix=TrixNeglect(Ctrl%AccL(I))
!               Thrsh%Dist=DistNeglect(Ctrl%AccL(I))
!               Thrsh%TwoE=TwoENeglect(Ctrl%AccL(I))
!               Thrsh%ETol=ETol(Ctrl%AccL(I))
!               Thrsh%DTol=DTol(Ctrl%AccL(I))
!               CALL Put(Thrsh,Tag_O=IntToChar(I))
!            ENDDO
!         ELSE ! Default, cheezy for all sets except last, which is good or user defined.
!            write(*,*) 'thrsh 2'
!            DO I=1,Ctrl%NSet-1
!               Thrsh%Cube=CubeNeglect(1)
!               Thrsh%Trix=TrixNeglect(1)
!               Thrsh%Dist=DistNeglect(1)
!               Thrsh%TwoE=TwoENeglect(1)
!               Thrsh%ETol=ETol(1)
!               Thrsh%DTol=DTol(1)
!               CALL Put(Thrsh,Tag_O=IntToChar(I))
!            ENDDO
!            write(*,*) 'thrsh 3'
!            Thrsh%Cube=CubeNeglect(Ctrl%AccL(1))
!            Thrsh%Trix=TrixNeglect(Ctrl%AccL(1))
!            Thrsh%Dist=DistNeglect(Ctrl%AccL(1))
!            Thrsh%TwoE=TwoENeglect(Ctrl%AccL(1))
!            Thrsh%ETol=ETol(Ctrl%AccL(1))
!            Thrsh%DTol=DTol(Ctrl%AccL(1))
!            CALL Put(Thrsh,Tag_O=IntToChar(Ctrl%NSet))
!         ENDIF
!
!        Close files
         CLOSE(UNIT=Inp,STATUS='KEEP')
         CALL CloseHDF()
      END SUBROUTINE ParseMethods
!============================================================================
!     Parse the Geometry Variable
!============================================================================
      SUBROUTINE ParseGeometry(Ctrl)
         TYPE(SCFControls),INTENT(INOUT) :: Ctrl
         TYPE(CRDS)                      :: GM
         INTEGER                         :: I,J,K,NKind,NUnPEl,QMCHARGE
         LOGICAL                         :: ReOrder,HilbertOrder
         CHARACTER(LEN=2)                :: At
         CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Line
!----------------------------------------------------------------------------
!        Open infofile
!
         CALL OpenHDF(InfFile)
!----------------------------------------------------------------------------
!        Open input file for parsing
!
         CALL OpenASCII(InpFile,Inp)
!----------------------------------------------------------------------------
!        Parse <OPTIONS> for <GEOMETRY> keys 
!
         GM%InAU=OptKeyQ(Inp,GEOMETRY,IN_AU)
         IF(OptKeyQ(Inp,GEOMETRY,NO_ORDER))THEN
            GM%Ordrd=SFC_NONE
         ELSEIF(OptKeyQ(Inp,GEOMETRY,Z_ORDER))THEN
            GM%Ordrd=SFC_PEANO
         ELSEIF(OptKeyQ(Inp,GEOMETRY,RANDOM_ORDER))THEN
            GM%Ordrd=SFC_RANDOM
         ELSEIF(OptKeyQ(Inp,GEOMETRY,H_Order))THEN
            GM%Ordrd=SFC_HILBERT 
         ELSEIF(OptKeyQ(Inp,GEOMETRY,TRAVEL_Order))THEN
            GM%Ordrd=SFC_TRAVEL 
         ELSEIF(OptKeyQ(Inp,GEOMETRY,TABLETRAV_ORDER))THEN
         GM%Ordrd=SFC_TABLETRAV
         ELSE
            CALL MondoHalt(PRSE_ERROR,'Unrecognized ordering.')
         ENDIF
!---------------------------------------------------------------------------- 
!        Parse <OPTIONS> for <TOT_CHARGE> and <MULTIPLICITY>  
!
         IF(.NOT.OptIntQ(Inp,TOTAL_CHARGE,QMCHARGE)) THEN
            CALL MondoHalt(PRSE_ERROR,TOTAL_CHARGE//' Not found in input.')
         ELSE
            GM%TotCh=DBLE(QMCHARGE)
         ENDIF
         IF(.NOT.OptIntQ(Inp,MULTIPLICITY,GM%Multp)) &
            CALL MondoHalt(PRSE_ERROR,MULTIPLICITY//' Not found in input.')
!
!---------------------------------------------------------------------------- 
!        Parse <OPTIONS> for <Extrap=>
!     
         IF(OptKeyQ(Inp,EXTRAPOLATE,EXTRAP_NO_EXTRAP))THEN
            Ctrl%Extrap=EXTRAP_GEOM_RSTRT
         ELSEIF(OptKeyQ(Inp,EXTRAPOLATE,EXTRAP_PROJECT))THEN
            Ctrl%Extrap=EXTRAP_GEOM_PRJCT
         ELSE
            Ctrl%Extrap=EXTRAP_GEOM_INTRP
         ENDIF
!---------------------------------------------------------------------------- 
!        Parse <OPTIONS> for <Vis=>
         IF(OptKeyQ(Inp,VISUALIZE,VIS_RHOPOT))THEN
            Ctrl%Vis=VIS_DX_RHOPOT
         ELSE
            Ctrl%Vis=VIS_DX_NO_VIS
         ENDIF
#ifdef PERIODIC
!----------------------------------------------------------------------------
!        Parse <OPTIONS> for <PERIODIC> keys
!        Which Directions are Periodic
!
         IF(OptKeyQ(Inp,PBOUNDRY,PBC_OFF)) THEN
            GM%PBC%AutoW(1)=.FALSE.
            GM%PBC%AutoW(2)=.FALSE.
            GM%PBC%AutoW(3)=.FALSE.
         ELSEIF(OptKeyQ(Inp,PBOUNDRY,PBC_X)) THEN
            GM%PBC%AutoW(1)=.TRUE.
            GM%PBC%AutoW(2)=.FALSE.
            GM%PBC%AutoW(3)=.FALSE.
         ELSEIF(OptKeyQ(Inp,PBOUNDRY,PBC_Y)) THEN
            GM%PBC%AutoW(1)=.FALSE.
            GM%PBC%AutoW(2)=.TRUE.
            GM%PBC%AutoW(3)=.FALSE.
         ELSEIF(OptKeyQ(Inp,PBOUNDRY,PBC_Z)) THEN
            GM%PBC%AutoW(1)=.FALSE.
            GM%PBC%AutoW(2)=.FALSE.
            GM%PBC%AutoW(3)=.TRUE.
         ELSEIF(OptKeyQ(Inp,PBOUNDRY,PBC_XY)) THEN
            GM%PBC%AutoW(1)=.TRUE.
            GM%PBC%AutoW(2)=.TRUE.
            GM%PBC%AutoW(3)=.FALSE.
         ELSEIF(OptKeyQ(Inp,PBOUNDRY,PBC_XZ)) THEN
            GM%PBC%AutoW(1)=.TRUE.
            GM%PBC%AutoW(2)=.FALSE.
            GM%PBC%AutoW(3)=.TRUE.
         ELSEIF(OptKeyQ(Inp,PBOUNDRY,PBC_YZ)) THEN
            GM%PBC%AutoW(1)=.FALSE.
            GM%PBC%AutoW(2)=.TRUE.
            GM%PBC%AutoW(3)=.TRUE.
         ELSEIF(OptKeyQ(Inp,PBOUNDRY,PBC_XYZ)) THEN
            GM%PBC%AutoW(1)=.TRUE.
            GM%PBC%AutoW(2)=.TRUE.
            GM%PBC%AutoW(3)=.TRUE.
         ELSE
            CALL OpenASCII(OutFile,Out)
            WRITE(Out,*) '** Auto-Wraps at default value => (Off) **'
            CLOSE(UNIT=Out,STATUS='KEEP')
            GM%PBC%AutoW(1)=.FALSE.
            GM%PBC%AutoW(2)=.FALSE.
            GM%PBC%AutoW(3)=.FALSE.
         ENDIF
!----------------------------------------------------------------------------
!        Calculate the Dimension
!         
         GM%PBC%Dimen=0
         DO I=1,3
            IF(GM%PBC%AutoW(I)) GM%PBC%Dimen=GM%PBC%Dimen+1
         ENDDO
!----------------------------------------------------------------------------
!        To wrap or not to wrap atoms into the box
!
         IF(OptKeyQ(Inp,PBOUNDRY,ATOMW_ON)) THEN
            GM%PBC%AtomW=.TRUE.
         ELSEIF(OptKeyQ(Inp,PBOUNDRY,ATOMW_OFF)) THEN
            GM%PBC%AtomW=.FALSE.
         ELSE
            CALL OpenASCII(OutFile,Out)
            WRITE(Out,*) '** Atom-Wrap at default value => (AtomWrap-Off) **'
            CLOSE(UNIT=Out,STATUS='KEEP')
            GM%PBC%AtomW=.FALSE.
         ENDIF
!----------------------------------------------------------------------------
!        Input Type on the BoxShape Vectors
!
         IF(OptKeyQ(Inp,PBOUNDRY,LVF_VEC)) THEN
            GM%PBC%InVecForm=.TRUE.
         ELSEIF(OptKeyQ(Inp,PBOUNDRY,LVF_ANG)) THEN
            GM%PBC%InVecForm=.FALSE.
         ELSE
            CALL OpenASCII(OutFile,Out)
            WRITE(Out,*) '** Lattice Vector Format at default value => (Vector Format) **'
            CLOSE(UNIT=Out,STATUS='KEEP')
            GM%PBC%InVecForm=.TRUE.
         ENDIF

!----------------------------------------------------------------------------
!        Input Type on the Coordinates, Atomic or Fractional
!
         IF(OptKeyQ(Inp,PBOUNDRY,CRT_ATOM)) THEN
            GM%PBC%InAtomCrd=.TRUE.
         ELSEIF (OptKeyQ(Inp,PBOUNDRY,CRT_FRAC)) THEN
            GM%PBC%InAtomCrd=.FALSE.
         ELSE
            CALL OpenASCII(OutFile,Out)
            WRITE(Out,*) '** Coodinate Format at default value => (Atomic Coord) **'
            CLOSE(UNIT=Out,STATUS='KEEP')
            GM%PBC%InAtomCrd=.TRUE.
         ENDIF
!-------------------------------------------------------------
!        Intput Type of Translation
!
         IF(OptKeyQ(Inp,PBOUNDRY,TRAN_COM)) THEN
            GM%PBC%Translate = .TRUE.
            GM%PBC%Trans_COM = .TRUE.
         ELSE
            GM%PBC%Translate = .FALSE.
            GM%PBC%Trans_COM = .FALSE.
         ENDIF
!-------------------------------------------------------------
!        Intput permability
!
         IF(.NOT. OptDblQ(Inp,EPSILON,GM%PBC%Epsilon)) THEN
            GM%PBC%Epsilon = 1.D32
         ENDIF
#endif
!----------------------------------------------------------------------------
!        Parse <OPTIONS> for <GEOMETRY> format
!
         IF(OptKeyQ(Inp,GEOMETRY,MSI_FORMAT))THEN
            CALL ParseCoordinates_MSI(Ctrl,GM)
         ELSEIF(OptKeyQ(Inp,GEOMETRY,XMOL_FORMAT))THEN
            CALL ParseCoordinates_XMOL(Ctrl,GM)
         ELSE
!          CALL New(GM)
           CALL ParseCoordinates_MONDO(Ctrl,GM)
         ENDIF
!
         CLOSE(UNIT=Inp,STATUS='KEEP')
         CALL CloseHDF()
!  
      END SUBROUTINE ParseGeometry
!---------------------------------------------------------------------------- 
!
!---------------------------------------------------------------------------- 
      SUBROUTINE ParseCoordinates_MONDO(Ctrl,GM)
         TYPE(SCFControls),INTENT(INOUT) :: Ctrl
         TYPE(CRDS)                      :: GM
         TYPE(INT_VECT)                  :: Kinds
         REAL(DOUBLE)                    :: R,Rx,Ry,Rz
         REAL(DOUBLE),DIMENSION(6)       :: Carts !,Vects
         INTEGER                         :: I,J,K,NumGeom
         CHARACTER(LEN=2)                :: At
         CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Line, LineLowCase
         LOGICAL                         :: LastConfig
!---------------------------------------------------------------------------- 
!        Find number of atoms and atom types
!
         CALL AlignLowCase('begingeometry',Inp)
         NAtoms=0
         DO 
           READ(Inp,DEFAULT_CHR_FMT,END=1)Line
           LineLowCase = Line
           Call LowCase(LineLowCase)
           IF(INDEX(LineLowCase,'endgeometry')/=0.OR. &
              INDEX(LineLowCase,'next_geometry_mondo_default')/=0)EXIT
           NAtoms=NAtoms+1
         ENDDO
         GM%NAtms=NAtoms
!        Allocate the geometry
         CALL New(GM)
!
#ifdef PERIODIC
!---------------------------------------------------------------------------- 
!        Parse <PERIODIC> for Lattice Vectors
!
         CALL ParsePeriodic_MONDO(Ctrl,GM)
#endif
!---------------------------------------------------------------------------- 
!        Parse <GEOMETRY> for coordinates
!
         CALL AlignLowCase('begingeometry',Inp)
         NumGeom   = 0
         LastConfig=.FALSE.
         DO
            NumGeom = NumGeom+1
            GM%Confg= NumGeom
            NAtoms=0
            DO 
               READ(Inp,DEFAULT_CHR_FMT,END=1)Line
               LineLowCase = Line
               Call LowCase(LineLowCase)
               IF(INDEX(LineLowCase,'next_geometry_mondo_default')/=0)EXIT
               IF(INDEX(LineLowCase,'endgeometry')/=0)THEN
                  LastConfig=.TRUE.
                  EXIT
               ELSE
                  LastConfig=.FALSE.
               ENDIF
               NAtoms=NAtoms+1
               CALL LineToGeom(Line,At,Carts)
               GM%Carts%D(:,NAtoms)=Carts(1:3) 
               GM%Vects%D(:,NAtoms)=Carts(4:6)
               DO J=1,104
                  IF(At==Ats(J))THEN
                     GM%AtNum%D(NAtoms)=J
                     GM%AtMss%D(NAtoms)=AtsMss(J)
                     EXIT
                  ENDIF
               ENDDO
            ENDDO
            IF(NAtoms/=GM%NAtms) &
                 CALL MondoHalt(PRSE_ERROR,'Atom number mismatch in ParseCoordinates_MONDO')
!           Reorder with SFC
            CALL ReorderCoordinates(GM)
!           Determine the number of atom types (kinds)
!           Do this AFTER possible reodering or it will mess up the basis set.
!           NOTE: Should recompute basis for EACH geometry.  Things are quite
!           messed up at present.  Need to really clean house on front end.
            CALL FindKind(GM)
#ifdef PERIODIC
!           Convert to AU and Compute Fractioan and Atomic Coordinates  
            IF(GM%PBC%InAtomCrd) THEN
               IF(.NOT.GM%InAU) THEN 
                  GM%Carts%D    = GM%Carts%D*AngstromsToAU 
               ENDIF
               CALL CalFracCarts(GM)
            ELSE
               GM%BoxCarts%D=GM%Carts%D
               GM%BoxVects%D=GM%Vects%D
               CALL CalAtomCarts(GM)
            ENDIF
!
            IF(GM%PBC%Trans_COM) THEN
               CALL CalTransVec(GM)
            ENDIF
            CALL Translate(GM,GM%PBC%TransVec)
            CALL WrapAtoms(GM)
#else
!           Convert to AU
            IF(.NOT.GM%InAU) THEN
               GM%Carts%D=GM%Carts%D*AngstromsToAU
            ENDIF
#endif
!
!           Compute spin coordinates
            CALL SpinCoords(GM) 
!           Determine a bounding box for the system
            GM%BndBox%D=SetBox(GM%Carts)
#ifdef PERIODIC
!           ReSet the Cell Center
            DO I=1,3
               IF(.NOT. GM%PBC%AutoW(I)) THEN
                  GM%PBC%CellCenter(I) = Half*(GM%BndBox%D(I,2)+GM%BndBox%D(I,1))
               ENDIF
            ENDDO
#endif
!           Output the coordinates
            CALL Put(GM,Tag_O=TRIM(IntToChar(NumGeom)))
!           Print the coordinates
            IF(PrintFlags%Key>DEBUG_NONE) CALL PPrint(GM)
!            CALL PPrint(GM,'Graphite_98.xyz',6,'XYZ')
!            IF(.TRUE.) STOP
!           Exit 
            IF(LastConfig)EXIT
         ENDDO
!
         IF(Ctrl%Grad == GRAD_ONE_FORCE .OR. Ctrl%Grad == GRAD_NO_GRAD) THEN
            Ctrl%NGeom = NumGeom
         ELSE
            IF(NumGeom .GT. 1) THEN
               CALL MondoHalt(PRSE_ERROR,'Only the initial Geometry should be supplied for MD or Opt')
            ENDIF
         ENDIF
         CALL Put(NumGeom,'NumberOfGeometries')
         NAtoms=GM%NAtms

!---------------------------------------------------------------------------- 
!        Finish up
         CALL Delete(GM)
!        CLOSE(UNIT=Inp,STATUS='KEEP') !!!! it's done above, was done even before
         CALL Put(Ctrl%NGeom,'nconfig') !!! currently only one MM geom
!        CALL CloseHDF()
         CALL PPrint(MemStats,'ParseGeometry')
         RETURN
!
1        CALL Halt('While parsing '//TRIM(InpFile)//', failed to find '     &
                   //TRIM(END_GEOMETRY)//'. You may be missing blank '  &
                   //'line at the end of the input file.')
      END SUBROUTINE ParseCoordinates_MONDO
!---------------------------------------------------------------------------- 
!
!---------------------------------------------------------------------------- 
      SUBROUTINE ParseCoordinates_XMOL(Ctrl,GM)
         TYPE(SCFControls),INTENT(INOUT) :: Ctrl
         TYPE(CRDS)                      :: GM
         TYPE(INT_VECT)                  :: Kinds
         REAL(DOUBLE)                    :: R,Rx,Ry,Rz
         REAL(DOUBLE),DIMENSION(3)       :: Carts(3)
         INTEGER                         :: I,J,K
         CHARACTER(LEN=2)                :: At
         CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Line
         LOGICAL                         :: LastConfig
         TYPE(CHR_VECT)                  :: Chars
!---------------------------------------------------------------------------- 
!        Find number of atoms and atom types
!
         CALL Align(BEGIN_GEOMETRY,Inp)
         READ(Inp,DEFAULT_CHR_FMT,END=1)Line
         READ(Inp,DEFAULT_CHR_FMT,END=1)Line
         NAtoms=0
         DO 
           READ(Inp,DEFAULT_CHR_FMT,END=1)Line
           CALL LineToChars(Line,Chars)
           K=SIZE(Chars%C)
           CALL Delete(Chars)
           IF(INDEX(Line,END_GEOMETRY)/=0.OR.K==1)EXIT
           NAtoms=NAtoms+1
         ENDDO
         GM%NAtms=NAtoms
!        Allocate the geometry
         CALL New(GM)
#ifdef PERIODIC
!---------------------------------------------------------------------------- 
!        Parse <PERIODIC> for Lattice Vectors
!
            CALL ParsePeriodic_XMOL(Ctrl,GM)
#endif
!---------------------------------------------------------------------------- 
!        Parse <GEOMETRY> for coordinates
!
         CALL Align(BEGIN_GEOMETRY,Inp)
         READ(Inp,DEFAULT_CHR_FMT,END=1)Line
         READ(Inp,DEFAULT_CHR_FMT,END=1)Line
!
         Ctrl%NGeom=0
         LastConfig=.FALSE.
         DO
            Ctrl%NGeom=Ctrl%NGeom+1
            GM%Confg=Ctrl%NGeom
            NAtoms=0
            DO 
               READ(Inp,DEFAULT_CHR_FMT,END=1)Line
               CALL LineToChars(Line,Chars)
               K=SIZE(Chars%C)
               CALL Delete(Chars)
               IF(INDEX(Line,END_GEOMETRY)/=0)THEN
                  LastConfig=.TRUE.
                  EXIT
               ELSE
                  LastConfig=.FALSE.
               ENDIF
               IF(K==1)THEN
                  READ(Inp,DEFAULT_CHR_FMT,END=1)Line
                  EXIT 
               ENDIF                 
               NAtoms=NAtoms+1
               CALL LineToGeom(Line,At,Carts)
               GM%Carts%D(:,NAtoms)=Carts(:)
               DO J=1,104
                  IF(At==Ats(J))THEN
                     GM%AtNum%D(NAtoms)=J
                     GM%AtMss%D(NAtoms)=AtsMss(J)
                     EXIT
                  ENDIF
               ENDDO
            ENDDO
            IF(NAtoms/=GM%NAtms) &
               CALL MondoHalt(PRSE_ERROR,'Atom number mismatch in ParseCoordinates_XMOL')
!           Reorder with SFC
            CALL ReorderCoordinates(GM)
!           Determine the number of atom types (kinds)
            CALL FindKind(GM)
!           Convert to AU
            IF(.NOT.GM%InAU) &
               GM%Carts%D=GM%Carts%D*AngstromsToAU
!           Compute spin coordinates
            CALL SpinCoords(GM) 
!           Determine a bounding box for the system
            GM%BndBox%D=SetBox(GM%Carts)
!           Output the coordinates
            CALL Put(GM,Tag_O=IntToChar(Ctrl%NGeom))
            CALL Put(Ctrl%NGeom,'NumberOfGeometries')
            IF(LastConfig)EXIT
         ENDDO
         NAtoms=GM%NAtms
!---------------------------------------------------------------------------- 
!        Finish up
         CALL Delete(GM)
         CLOSE(UNIT=Inp,STATUS='KEEP')
         CALL Put(Ctrl%NGeom,'nconfig')
         CALL CloseHDF()
         CALL PPrint(MemStats,'ParseGeometry')
         RETURN
       1 CALL Halt('While parsing '//TRIM(InpFile)//', failed to find '     &
                   //TRIM(END_GEOMETRY)//'. \n  You may be missing blank '  &
                   //'line at the end of the input file.')
!
      END SUBROUTINE ParseCoordinates_XMOL
!
      SUBROUTINE ParseCoordinates_MSI(Ctrl,GM)
         TYPE(SCFControls),INTENT(INOUT) :: Ctrl
         TYPE(CRDS)                      :: GM
         TYPE(INT_VECT)                  :: Kinds
         REAL(DOUBLE)                    :: R,Rx,Ry,Rz
         REAL(DOUBLE),DIMENSION(3)       :: Carts(3)
         TYPE(CHR_VECT)                  :: Chars
         INTEGER                         :: I,J,K
         CHARACTER(LEN=2)                :: At
         CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Line
         LOGICAL                         :: LastConfig
!---------------------------------------------------------------------------- 
!        Find number of atoms and atom types
!
         CALL Align(BEGIN_GEOMETRY,Inp)
         NAtoms=0
         DO 
           READ(Inp,DEFAULT_CHR_FMT,END=1)Line
           IF(INDEX(Line,BEGIN_NEXT_GEOMETRY_MSI_ARCHIVE)/=0)THEN
               READ(Inp,DEFAULT_CHR_FMT,END=1)Line
               READ(Inp,DEFAULT_CHR_FMT,END=1)Line
           ENDIF
           IF(INDEX(Line,END_NEXT_GEOMETRY_MSI_ARCHIVE)/=0)EXIT
           NAtoms=NAtoms+1
         ENDDO
         GM%NAtms=NAtoms
!        Allocate the geometry
         CALL New(GM)
#ifdef PERIODIC
!---------------------------------------------------------------------------- 
!        Parse <PERIODIC> for Lattice Vectors
!
            CALL ParsePeriodic_MSI(Ctrl,GM)
#endif
!---------------------------------------------------------------------------- 
!        Parse <GEOMETRY> for coordinates
!
         CALL Align(BEGIN_GEOMETRY,Inp)
         Ctrl%NGeom=0
         LastConfig=.FALSE.
         DO
            Ctrl%NGeom=Ctrl%NGeom+1
            GM%Confg=Ctrl%NGeom
            NAtoms=0
            DO 
               READ(Inp,DEFAULT_CHR_FMT,END=1)Line
               IF(INDEX(Line,BEGIN_NEXT_GEOMETRY_MSI_ARCHIVE)/=0)THEN
                  READ(Inp,DEFAULT_CHR_FMT,END=1)Line
                  READ(Inp,DEFAULT_CHR_FMT,END=1)Line
               ENDIF
               IF(INDEX(Line,END_NEXT_GEOMETRY_MSI_ARCHIVE)/=0)THEN
                  READ(Inp,DEFAULT_CHR_FMT,END=1)Line
                  EXIT
               ENDIF
               IF(INDEX(Line,END_GEOMETRY)/=0)THEN
                  LastConfig=.TRUE.
                  EXIT
               ELSE
                  LastConfig=.FALSE.
               ENDIF
               NAtoms=NAtoms+1
               CALL LineToDbls(Line,3,Carts)
               GM%Carts%D(:,NAtoms)=Carts(:)
               CALL LineToChars(Line,Chars)
               At=TRIM(ADJUSTL(Chars%C(8)))
               CALL LowCase(At)
               DO J=1,104
                  IF(TRIM(At)==TRIM(Ats(J)))THEN
                     GM%AtNum%D(NAtoms)=J
                     GM%AtMss%D(NAtoms)=AtsMss(J)
                     EXIT
                  ENDIF
               ENDDO
            ENDDO
            IF(LastConfig)EXIT
            IF(NAtoms/=GM%NAtms) &
               CALL MondoHalt(PRSE_ERROR,'Atom number mismatch in ParseCoordinates_MSI')
!           Reorder with SFC
            CALL ReorderCoordinates(GM)
!           Determine the number of atom types (kinds)
            CALL FindKind(GM)
!           Convert to AU
            IF(.NOT.GM%InAU) &
               GM%Carts%D=GM%Carts%D*AngstromsToAU
!           Compute spin coordinates
            CALL SpinCoords(GM) 
!           Determine a bounding box for the system
            GM%BndBox%D=SetBox(GM%Carts)
!           Compute nuclear-nuclear repulsion energy
!            GM%ENucN=ENukeNuke(GM)
!           Output the coordinates
!            CALL OpenASCII(OutFile,Out)
!            CALL PPrint(GM)            
!            CLOSE(Out)
            CALL Put(GM,Tag_O=IntToChar(Ctrl%NGeom))
            CALL Put(Ctrl%NGeom,'NumberOfGeometries')
         ENDDO
         NAtoms=GM%NAtms
!---------------------------------------------------------------------------- 
!        Finish up
!
         CALL Delete(GM)
         CALL Delete(Chars)
         CLOSE(UNIT=Inp,STATUS='KEEP')
         CALL Put(Ctrl%NGeom,'nconfig')
         CALL CloseHDF()
         CALL PPrint(MemStats,'ParseGeometry')
         RETURN
1        CALL Halt('While parsing '//TRIM(InpFile)//', failed to find '     &
                   //TRIM(END_GEOMETRY)//'. You may be missing blank '  &
                   //'line at the end of the input file.')
!
      END SUBROUTINE ParseCoordinates_MSI
#ifdef PERIODIC
!============================================================================
!
!============================================================================
      SUBROUTINE ParsePeriodic_MONDO(Ctrl,GM)
        TYPE(SCFControls),INTENT(INOUT) :: Ctrl
        TYPE(CRDS)                      :: GM
        INTEGER                         :: NLvec,NTvec,I,J,K
!
        NLvec = 0
        NTvec = 0
        GM%PBC%TransVec = Zero
        GM%PBC%BoxShape = Zero
        GM%PBC%InvBoxSh = Zero
!
        DO I=1,3
           GM%PBC%BoxShape(I,I) = One
           GM%PBC%InvBoxSh(I,I) = One
        ENDDO
!---------------------------------------------------------------------------- 
!       Get the Lattice and Translate Vectors
!
        IF(FindKey(BEGIN_PERIODIC,Inp)) THEN
           CALL Align(BEGIN_PERIODIC,Inp)
           IF(GM%PBC%InVecForm) THEN
              CALL GetLatVec(GM,NLvec,NTvec)
           ELSE
              CALL GetLatAng(GM,NLvec,NTvec)
           ENDIF
        ELSE
           IF(GM%PBC%Dimen .NE. 0) THEN
              CALL MondoHalt(PRSE_ERROR,'No Lattice Vectors where supplied')
           ENDIF
        ENDIF
!--------------------------------------------------------------
!       Logic for the Translate and Lattice Vector
!
        IF(.NOT. GM%PBC%Trans_COM) THEN
           IF(NTvec == 0) THEN
              GM%PBC%Translate = .FALSE.
           ELSEIF(NTvec == 1) THEN
              GM%PBC%Translate = .TRUE.
           ELSE
              CALL MondoHalt(PRSE_ERROR,'Number of Translate Vectors is Incorrect')
           ENDIF
        ENDIF
!--------------------------------------------------------------
        IF(NLvec .LT. GM%PBC%Dimen) THEN
           CALL MondoHalt(PRSE_ERROR,'Number of Lattice Vectors is Incorrect')      
        ENDIF

!---------------------------------------------------------------------------- 
!       Convert the lattice and translate vectors to AU 
!
        IF(.NOT.GM%InAU) THEN 
           GM%PBC%BoxShape = GM%PBC%BoxShape*AngstromsToAU
           GM%PBC%TransVec = GM%PBC%TransVec*AngstromsToAU
        ENDIF
!----------------------------------------------------------------------------
!       Calculate the Box Volume
!
        GM%PBC%CellVolume = One
        DO I=1,3
           IF(GM%PBC%AutoW(I)) THEN
              GM%PBC%CellVolume = GM%PBC%CellVolume*GM%PBC%BoxShape(I,I)
           ENDIF
        ENDDO
!----------------------------------------------------------------------------
!       Calculate the Dipole and Quadripole Factors
!
        IF(GM%PBC%Dimen < 2) THEN
           GM%PBC%DipoleFAC = Zero
           GM%PBC%QupoleFAC = Zero
        ELSEIF(GM%PBC%Dimen ==2) THEN
           GM%PBC%DipoleFAC = (Four*Pi/GM%PBC%CellVolume)*(One/(GM%PBC%Epsilon+One)) 
           GM%PBC%QupoleFAC =  Zero
           IF(ABS(GM%PBC%DipoleFAC) .LT. 1.D-14) GM%PBC%DipoleFAC = Zero
        ELSEIF(GM%PBC%Dimen ==3) THEN
           GM%PBC%DipoleFAC = -(Four*Pi/GM%PBC%CellVolume)*(One/Three - One/(Two*GM%PBC%Epsilon+One)) 
           GM%PBC%QupoleFAC =  (Two*Pi/GM%PBC%CellVolume)*(One/Three - One/(Two*GM%PBC%Epsilon+One)) 
           IF(ABS(GM%PBC%DipoleFAC) .LT. 1.D-14) GM%PBC%DipoleFAC = Zero
           IF(ABS(GM%PBC%QupoleFAC) .LT. 1.D-14) GM%PBC%QupoleFAC = Zero
        ENDIF
!----------------------------------------------------------------------------
!       Calculate the Cell Center
!
        DO I=1,3
           GM%PBC%CellCenter(I) = Zero
           IF(GM%PBC%AutoW(I)) THEN
              DO J=1,3
                 IF(GM%PBC%AutoW(J)) THEN
                    GM%PBC%CellCenter(I) =  GM%PBC%CellCenter(I) + Half*GM%PBC%BoxShape(I,J)
                 ENDIF
              ENDDO
           ENDIF
        ENDDO
!----------------------------------------------------------------------------
!       Calculate The Inverse of BoxShape  InvBoxSh = [BoxShape]^(-1)
!
        GM%PBC%InvBoxSh(1,1)=One/GM%PBC%BoxShape(1,1)
        GM%PBC%InvBoxSh(2,1)=Zero
        GM%PBC%InvBoxSh(3,1)=Zero
!
        GM%PBC%InvBoxSh(1,2)=-GM%PBC%BoxShape(1,2)/(GM%PBC%BoxShape(1,1)*GM%PBC%BoxShape(2,2))
        GM%PBC%InvBoxSh(2,2)=One/GM%PBC%BoxShape(2,2)
        GM%PBC%InvBoxSh(3,2)=Zero
! 
        GM%PBC%InvBoxSh(1,3)=(GM%PBC%BoxShape(1,2)*GM%PBC%BoxShape(2,3)  &
                            - GM%PBC%BoxShape(2,2)*GM%PBC%BoxShape(1,3)) &
                            /(GM%PBC%BoxShape(1,1)*GM%PBC%BoxShape(2,2)*GM%PBC%BoxShape(3,3))
        GM%PBC%InvBoxSh(2,3)=-GM%PBC%BoxShape(2,3)/(GM%PBC%BoxShape(2,2)*GM%PBC%BoxShape(3,3))
        GM%PBC%InvBoxSh(3,3)=One/GM%PBC%BoxShape(3,3)
!
      END SUBROUTINE ParsePeriodic_MONDO
!============================================================================
!
!============================================================================
      SUBROUTINE ParsePeriodic_MSI(Ctrl,GM)
        TYPE(SCFControls),INTENT(INOUT) :: Ctrl
        TYPE(CRDS)                      :: GM
        CALL MondoHalt(PRSE_ERROR,'MSI Format with Periodic not supported')   
      END SUBROUTINE ParsePeriodic_MSI
!
!
      SUBROUTINE ParsePeriodic_XMOL(Ctrl,GM)
        TYPE(SCFControls),INTENT(INOUT) :: Ctrl
        TYPE(CRDS)                      :: GM
        CALL MondoHalt(PRSE_ERROR,'XMOL Format with Periodic not supported')   
      END SUBROUTINE ParsePeriodic_XMOL
!============================================================================
!
!============================================================================
      SUBROUTINE GetLatVec(GM,NLvec,NTvec)
        TYPE(CRDS)                      :: GM
        INTEGER                         :: NLvec,NTvec
        CHARACTER(LEN=2)                :: At
        CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Line
        REAL(DOUBLE),DIMENSION(3)       :: Vec
!
        DO 
           READ(Inp,DEFAULT_CHR_FMT,END=1) Line
           IF(INDEX(Line,END_PERIODIC)==0)THEN
              CALL LineToGeom(Line,At,Vec)
              IF(At==TRAN_VEC) THEN
                 NTvec = NTvec + 1
                 GM%PBC%TransVec(1:3)  =Vec(1:3)
              ENDIF
              IF(At==ALAT_VEC) THEN
                 NLvec = NLvec+1
                 GM%PBC%BoxShape(1:3,1)=Vec(1:3)
              ENDIF
              IF(At==BLAT_VEC) THEN
                 NLvec = NLvec+1
                 GM%PBC%BoxShape(1:3,2)=Vec(1:3)
              ENDIF
              IF(At==CLAT_VEC) THEN
                 NLvec = NLvec+1
                 GM%PBC%BoxShape(1:3,3)=Vec(1:3)
              ENDIF
           ELSE
              EXIT
           ENDIF
        ENDDO
        RETURN
1       CALL Halt('While parsing '//TRIM(InpFile)//', failed to find '     &
                   //TRIM(END_PERIODIC)//'. You may be missing blank '  &
                   //'line at the end of the input file.')
!
      END SUBROUTINE GetLatVec
!============================================================================
!
!============================================================================
      SUBROUTINE GetLatAng(GM,NLvec,NTvec)
        TYPE(CRDS)                      :: GM
        INTEGER                         :: NLvec,NTvec,Dimen,I,J
        CHARACTER(LEN=2)                :: At
        CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Line
        REAL(DOUBLE),PARAMETER          :: DegToRad =  1.745329251994329576923D-2
        REAL(DOUBLE),DIMENSION(6)       :: Vec
        REAL(DOUBLE)                    :: AngAB,AngAC,AngBC,Error
!
        DO 
           READ(Inp,DEFAULT_CHR_FMT,END=1) Line
           IF(INDEX(Line,END_PERIODIC)==0) THEN
              CALL LineToGeom(Line,At,Vec)
              IF(At==TRAN_VEC) THEN
                 NTvec = NTvec + 1
                 GM%PBC%TransVec(1:3)=Vec(1:3)
              ELSE
                 J = 0
                 DO I=1,6
                    IF(Vec(I) == Zero) EXIT
                    J = J + 1
                 ENDDO
                 IF(J==1) THEN
                    NLvec = 1
                    IF(GM%PBC%AutoW(1)) THEN
                       GM%PBC%BoxShape(1,1)=Vec(1)            
                    ELSEIF(GM%PBC%AutoW(2)) THEN
                       GM%PBC%BoxShape(2,2)=Vec(1)
                    ELSEIF(GM%PBC%AutoW(3)) THEN                             
                       GM%PBC%BoxShape(3,3)=Vec(1)
                    ENDIF
                 ELSEIF(J==3) THEN                 
                    NLvec = 2
                    IF(GM%PBC%AutoW(1) .AND. GM%PBC%AutoW(2)) THEN
                       GM%PBC%BoxShape(1,1)=Vec(1)
                       GM%PBC%BoxShape(1,2)=Vec(2)*COS(DegToRad*Vec(3))
                       GM%PBC%BoxShape(2,2)=Vec(2)*SIN(DegToRad*Vec(3)) 
                    ELSEIF(GM%PBC%AutoW(1) .AND. GM%PBC%AutoW(3)) THEN
                       GM%PBC%BoxShape(1,1)=Vec(1)
                       GM%PBC%BoxShape(1,3)=Vec(2)*COS(DegToRad*Vec(3))
                       GM%PBC%BoxShape(3,3)=Vec(2)*SIN(DegToRad*Vec(3))
                    ELSEIF(GM%PBC%AutoW(2) .AND. GM%PBC%AutoW(3)) THEN
                       GM%PBC%BoxShape(2,2)=Vec(1)
                       GM%PBC%BoxShape(2,3)=Vec(2)*COS(DegToRad*Vec(3))
                       GM%PBC%BoxShape(3,3)=Vec(2)*SIN(DegToRad*Vec(3))
                    ENDIF
                 ELSEIF(J==6) THEN
                    NLvec = 3
                    GM%PBC%BoxShape(1,1)=Vec(1)
!
                    GM%PBC%BoxShape(1,2)=Vec(2)*COS(DegToRad*Vec(6))
                    GM%PBC%BoxShape(2,2)=Vec(2)*SIN(DegToRad*Vec(6))
!
                    GM%PBC%BoxShape(1,3)=Vec(3)*COS(DegToRad*Vec(5))
                    GM%PBC%BoxShape(2,3)=(Vec(2)*Vec(3)*COS(DegToRad*Vec(4)) &
                                       -GM%PBC%BoxShape(1,2)*GM%PBC%BoxShape(1,3))/ GM%PBC%BoxShape(2,2)
                    GM%PBC%BoxShape(3,3)=SQRT(Vec(3)**2-GM%PBC%BoxShape(1,3)**2-GM%PBC%BoxShape(2,3)**2)
!
                    AngAB = ACOS(GM%PBC%BoxShape(1,1)*GM%PBC%BoxShape(1,2)/(Vec(1)*Vec(2)))/DegToRad
                    AngAC = ACOS(GM%PBC%BoxShape(1,1)*GM%PBC%BoxShape(1,3)/(Vec(1)*Vec(3)))/DegToRad
                    AngBC = GM%PBC%BoxShape(1,2)*GM%PBC%BoxShape(1,3)+GM%PBC%BoxShape(2,2)*GM%PBC%BoxShape(2,3)
                    AngBC = ACOS(AngBC/(Vec(2)*Vec(3)))/DegToRad
                    Error = ABS(AngAB-Vec(6))+ABS(AngAC-Vec(5))+ ABS(AngBC-Vec(4))
!
                    IF(Error .GT. 1.D-6 ) THEN
                       CALL MondoHalt(PRSE_ERROR,'Angles Are Inccorect')
                    ENDIF
                 ELSE
                    CALL MondoHalt(PRSE_ERROR,'Number of Magnitudes and Angles supplied is incorrect')
                 ENDIF
              ENDIF
           ELSE
              EXIT
           ENDIF
        ENDDO
!
        DO I=1,3
           DO J=1,3
              IF(ABS(GM%PBC%BoxShape(I,J)).LT. 1.D-12) GM%PBC%BoxShape(I,J)=Zero
           ENDDO
        ENDDO
!
        RETURN
1       CALL Halt('While parsing '//TRIM(InpFile)//', failed to find '     &
                   //TRIM(END_PERIODIC)//'. You may be missing blank '  &
                   //'line at the end of the input file.')
!
      END SUBROUTINE GetLatAng
!============================================================================
!
!============================================================================
      SUBROUTINE  CalTransVec(GM)
        TYPE(CRDS)                  :: GM
        REAL(DOUBLE),DIMENSION(1:3) :: CMVec
        INTEGER                     :: I
!
        CMVec(:)=Zero
        DO I=1,GM%NAtms
           CMVec(:) = CMVec(:)+GM%BoxCarts%D(:,I)
        ENDDO
        CMVec(:)           = Half-CMVec(:)/DBLE(GM%NAtms)
        GM%PBC%TransVec(:) = FracToAtom(GM,CMVec(:))
!
        DO I=1,3
           IF(.NOT. GM%PBC%AutoW(I))  GM%PBC%TransVec(I) = Zero
        ENDDO
!
      END SUBROUTINE CalTransVec
#endif
!============================================================================
!
!============================================================================
      SUBROUTINE SpinCoords(GM) 
         TYPE(CRDS)     :: GM
         INTEGER        :: I,J,NUnPEl
!----------------------------------------------------------------------------
!        Calculate the electronic coordinates
         GM%NElec=0
         DO I=1,GM%NAtms
            GM%NElec=GM%NElec+GM%AtNum%D(I)
         ENDDO
         GM%NElec=GM%NElec-GM%TotCh

         NUnPEl=GM%Multp-1
         IF(NUnPEl.NE.0) &
            CALL MondoHalt(PRSE_ERROR,'Open shell not supported yet.'   &
                           //' NElectrons = '//TRIM(IntToChar(GM%NElec)) &
                           //' NUnPairedE = '//TRIM(IntToChar(NUnPEl)))
         GM%NAlph=DBLE(GM%NElec+NUnPEl)*Half
         GM%NBeta=DBLE(GM%NElec-NUnPEl)*Half
      END SUBROUTINE SpinCoords
!============================================================================
!
!============================================================================
      SUBROUTINE FindKind(GM) 
         TYPE(CRDS)     :: GM
         TYPE(INT_VECT) :: Kinds
         INTEGER        :: I,J,K
!----------------------------------------------------------------------------
         CALL New(Kinds,GM%NAtms)
         GM%NKind=1
         Kinds%I(1)=GM%AtNum%D(1)
         DO I=2,GM%NAtms
            DO J=1,GM%NKind
               IF(Kinds%I(J)==GM%AtNum%D(I))GOTO 3
            ENDDO
            GM%NKind=GM%NKind+1
            Kinds%I(GM%NKind)=GM%AtNum%D(I)
         3 CONTINUE
         ENDDO
         DO K=1,GM%NKind
            DO I=1,GM%NAtms
               IF(Kinds%I(K)==GM%AtNum%D(I))GM%AtTyp%I(I)=K
           ENDDO
         ENDDO
         CALL Delete(Kinds)
      END SUBROUTINE FindKind
!

!============================================================================
!
!============================================================================
!      FUNCTION ENukeNuke(GM) RESULT(ENucN)
!         TYPE(CRDS)   :: GM
!         REAL(DOUBLE) :: ENucN,R
!         INTEGER      :: I,J
!----------------------------------------------------------------------------
!         ENucN=0.D0
!         DO I=1,NAtoms
!            DO J=1,I-1
!               R=DSQRT((GM%Carts%D(1,i)-GM%Carts%D(1,j))**2 & 
!                      +(GM%Carts%D(2,i)-GM%Carts%D(2,j))**2 & 
!                      +(GM%Carts%D(3,i)-GM%Carts%D(3,j))**2 )
!               ENucN=ENucN+DBLE(GM%AtNum%D(I)*GM%AtNum%D(J))/R
!            ENDDO
!         ENDDO
!      END FUNCTION ENukeNuke
!============================================================================
!
!============================================================================
      FUNCTION SetBox(Carts,Carts2) RESULT(Box)
         TYPE(DBL_RNK2)               :: Carts
         TYPE(DBL_RNK2),OPTIONAL,INTENT(IN) :: Carts2
         REAL(DOUBLE), DIMENSION(3,2) :: Box
         INTEGER                      :: J
!----------------------------------------------------------------------------
         Box(:,1)=+1.D8
         Box(:,2)=-1.D8
         DO J=1,SIZE(Carts%D,2)
            Box(:,1)=MIN(Box(:,1),Carts%D(:,J))
            Box(:,2)=MAX(Box(:,2),Carts%D(:,J))
         ENDDO
       IF(PRESENT(Carts2)) THEN
         DO J=1,SIZE(Carts2%D,2)
            Box(:,1)=MIN(Box(:,1),Carts2%D(:,J))
            Box(:,2)=MAX(Box(:,2),Carts2%D(:,J))
         ENDDO
       ENDIF
      END FUNCTION SetBox
!============================================================================
!
!============================================================================
      SUBROUTINE ReorderCoordinates(GM)
         TYPE(CRDS)     :: GM
         TYPE(DBL_VECT) :: DTemp
         TYPE(INT_VECT) :: ITemp,Kinds,Point
         INTEGER        :: J
!----------------------------------------------------------------------------
         IF(GM%Ordrd==SFC_NONE)RETURN
!
         CALL New(Point,NAtoms)
         CALL New(DTemp,NAtoms)
         CALL New(ITemp,NAtoms)
         CALL SFCOrder(NAtoms,GM%Carts,Point,GM%Ordrd)
!        Reorder Coordinates
         DO I=1,3
            DO J=1,NAtoms
               DTemp%D(J)=GM%Carts%D(I,J)
            ENDDO
            DO J=1,NAtoms
               GM%Carts%D(I,J)=DTemp%D(Point%I(J))
            ENDDO
         ENDDO
!
#ifdef PERIODIC
!        Reorder Fractional Coordinates
         DO I=1,3
            DO J=1,NAtoms
               DTemp%D(J)=GM%BoxCarts%D(I,J)
            ENDDO
            DO J=1,NAtoms
               GM%BoxCarts%D(I,J)=DTemp%D(Point%I(J))
            ENDDO
         ENDDO
#endif

!        Reorder Atomic Number and Mass
         DO J=1,NAtoms
            ITemp%I(J)=GM%AtNum%D(J)
         ENDDO
         DO J=1,NAtoms
            GM%AtNum%D(J)=ITemp%I(Point%I(J))
         ENDDO
!
         DO J=1,NAtoms
            DTemp%D(J)=GM%AtMss%D(J)
         ENDDO
         DO J=1,NAtoms
            GM%AtMss%D(J)=DTemp%D(Point%I(J))
         ENDDO
!
         CALL Delete(Point)
         CALL Delete(DTemp)
         CALL Delete(ITemp)
     END SUBROUTINE ReorderCoordinates
!============================================================================
!     Parse the Basis Sets
!============================================================================
      SUBROUTINE ParseBaseSets(Ctrl)
         USE BasisSetParameters
         IMPLICIT NONE
!----------------------------------------------------------------------------
         TYPE(SCFControls)                     :: Ctrl
         TYPE(BSET)                            :: BS
         TYPE(BSET),DIMENSION(MaxSets)         :: Base
         TYPE(DBL_VECT)                        :: ZAtNum
         TYPE(INT_VECT)                        :: Types,BSiz_2,OffS_2
!        TYPE(INT_VECT)                        :: ZAtNum,Types,BSiz_2,OffS_2
         REAL(DOUBLE), DIMENSION(1:MaxASymt+2) :: Dbls
         CHARACTER(LEN=DEFAULT_CHR_LEN)        :: Line
         INTEGER                               :: I,J,K,L,N,NS,NP,NC,NK, &
                                                  MinL,MaxL,ISet,KFound
!----------------------------------------------------------------------------
!        Open infofile
         CALL OpenHDF(InfFile)
!        Parse basis sets from input file
         CALL OpenASCII(InpFile,Inp) 
!        Parse <OPTIONS.BASIS_SETS>
         Ctrl%NSet=0
         Base(:)%BType=0
         DO I=1,NSupSets
            IF(OptKeyLocQ(Inp,BASIS_SETS,TRIM(CSets(1,I)),MaxSets,NLoc,Loc))THEN
               Ctrl%NSet=Ctrl%NSet+NLoc
               DO ILoc=1,NLoc 
                  Base(Loc(ILoc))%BName=ADJUSTL(CSets(2,I))
                  Base(Loc(ILoc))%BType=I
               ENDDO
            ENDIF       
         ENDDO
         CLOSE(UNIT=Inp,STATUS='KEEP')
!
         IF(Ctrl%NSet==0)CALL Halt('<'//TRIM(BASIS_SETS)  &
           //'> keyword not found in '//TRIM(InpFile))  
!
!        If using a guess of AO density superposition, check to see
!        if the first basis set is minimial. If not, insert one. 
!
#ifdef LEGACY
         IF(INDEX(Base(1)%BName,'STO')==0.AND.Ctrl%SuperP)THEN
            CALL New(Types,Ctrl%NSet)
            Types%I(1:Ctrl%NSet)=Base(1:Ctrl%NSet)%BType
            Base(1)%BType=1 
            Base(2:Ctrl%NSet+1)%BType=Types%I(1:Ctrl%NSet)
            Ctrl%NSet=Ctrl%NSet+1
            DO I=1,Ctrl%NSet
               Base(I)%BName=ADJUSTL(CSets(2,Base(I)%BType))
               Ctrl%BName(I)=Base(I)%BName
            ENDDO               
            CALL Delete(Types)
         ELSE
            DO I=1,Ctrl%NSet
               Ctrl%BName(I)=Base(I)%BName
            ENDDO               
         ENDIF
#else
            DO I=1,Ctrl%NSet
               Ctrl%BName(I)=Base(I)%BName
            ENDDO               
#endif
!----------------------------------------------------------------------------
!         Allocate a temp basis set
!
          BS%LMNLen=MaxAsymt
          BS%NCtrt =MaxCntrx
          BS%NPrim =MaxPrmtv
          CALL Get(BS%NAtms,'natoms',Tag_O='1')
          CALL Get(BS%NKind,'nkind',Tag_O='1')
          CALL New(BS)
          CALL New(ZAtNum,NAtoms)
          CALL Get(ZAtNum,'atomicnumbers',Tag_O='1')
!----------------------------------------------------------------------------
!         Count kinds and load kind array
!
          BS%NKind=1
          BS%Kinds%I(1)=ZAtNum%D(1)
          DO I=2,NAtoms 
             DO J=1,BS%NKind
                IF(BS%Kinds%I(J)==ZAtNum%D(I))GOTO 10
             ENDDO
             BS%NKind=BS%NKind+1
             BS%Kinds%I(BS%NKind)=ZAtNum%D(I)
          10 CONTINUE
          ENDDO
!----------------------------------------------------------------------------
!        Go over each selected basis set
!
         DO ISet=1,Ctrl%NSet    
            CSet=IntToChar(ISet)
!           Open basis set file 
            IF(Base(ISet)%BType==0) &
               CALL MondoHalt(-9999,' Basis set # '//TRIM(IntToChar(ISet)) &
                                     //' Undefined in input!')
            BasFile=TRIM(MONDO_HOME)//'BasisSets/' &
                  //TRIM(CSets(2,Base(ISet)%BType))//BasF
            CALL OpenASCII(BasFile,Bas,Rewind_O=.TRUE.)
            IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
               CALL PPrint('Base('//TRIM(IntToChar(ISet))       &
                        //')%BName   = '//TRIM(Base(ISet)%BName))
               CALL PPrint('Opening BasFile = '//TRIM(BasFile))
            ENDIF
!           zero counters
            BS%NASym=0
            DO I=1,BS%NKind
               BS%NCFnc%I(I)=0
               DO J=1,MaxCntrx
                  BS%NPFnc%I(J,I)=0
               ENDDO
            ENDDO
!           Parse basis set 
            BS%NCtrt=0
            BS%NPrim=0
            KFound=0
            DO 
               READ(Bas,DEFAULT_CHR_FMT,END=99)Line                 
               DO NK=1,BS%NKind            
                  IF(KeyQ(Line,Ats(BS%Kinds%I(NK))).AND.KeyQ(Line,'0'))THEN                    
                     NC=0
                     KFound=KFound+1
                     DO 
                        READ(Bas,DEFAULT_CHR_FMT,END=99)Line         
                        IF(KeyQ(Line,Stars))GOTO 100
                        NC=NC+1                 
                        DO K=1,MaxLTyps
                           IF(KeyQ(Line,CLTyps(K)))THEN
                              BS%ASymm%I(1,NC,NK)=LTyps(1,K)                     
                              BS%ASymm%I(2,NC,NK)=LTyps(2,K)
                              DO L=1,MaxPrmtv
                                 IF(KeyQ(Line,TRIM(IntToChar(L))))THEN
                                    NP=L
                                    GOTO 101
                                 ENDIF
                              ENDDO
                           ENDIF
                        ENDDO                                    
                        CALL Halt('Logic failure in ParseBaseSets for ' &
                                 //TRIM(BasFile)//Rtrn//'Line = '//TRIM(Line))
                    101 CONTINUE
                        BS%NCFnc%I(NK)=NC
                        BS%NPFnc%I(NC,NK)=NP
                        BS%NCtrt=MAX(BS%NCtrt,NC)
                        BS%NPrim=MAX(BS%NPrim,NP)
                        MinL=BS%ASymm%I(1,NC,NK)
                        MaxL=BS%ASymm%I(2,NC,NK)
                        BS%NAsym=MAX(BS%NAsym,MaxL)
                         DO NP=1,BS%NPFnc%I(NC,NK)
                           READ(Bas,DEFAULT_CHR_FMT,END=99)Line
                           N=MaxL-MinL+2
                           CALL LineToDbls(Line,N,Dbls)
                           BS%Expnt%D(NP,NC,NK)=Dbls(1)
                           K=1
                           BS%CCoef%D(1:,NP,NC,NK)=Zero
                           DO NS=MinL,MaxL
                              K=K+1
                              BS%CCoef%D(NS+1,NP,NC,NK)=Dbls(K) 
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
           100 CONTINUE
            ENDDO            
         99 CONTINUE
            CLOSE(UNIT=Bas,STATUS='KEEP')
            IF(KFound/=BS%NKind)            & 
               CALL MondoHalt(PRSE_ERROR,   &
               'Unable to find all the atom types in ' &
               //TRIM(BasFile))
            CALL NormalizeBaseSets(BS,Base(ISet))
            CALL CalcDimensions(Base(ISet),BSiz_2,OffS_2)
            CALL Put(Base(ISet),Tag_O=CSet)
            CALL Put(BSiz_2,'atsiz',Tag_O=CSet)
            CALL Put(OffS_2,'atoff',Tag_O=CSet)
            CALL Delete(BSiz_2)
            CALL Delete(OffS_2)
         ENDDO
!----------------------------------------------------------------------------
!        Tidy up
!
         CALL Delete(BS)
         DO I=1,Ctrl%NSet
            CALL Delete(Base(I))
         ENDDO
         CALL Put(Ctrl%NSet,'NumberOfSets')
         IF(PrintFlags%Key==DEBUG_MAXIMUM)  &
            CALL PPrint(MemStats,'ParseBaseSets')
         CALL CloseHDF()
      END SUBROUTINE ParseBaseSets
!============================================================================
      SUBROUTINE NormalizeBaseSets(A,B)
         IMPLICIT NONE
         TYPE(BSET), INTENT(INOUT)   :: A
         TYPE(BSET), INTENT(OUT)     :: B
         REAL(DOUBLE)                :: Z,C,Expnt,RNorm,ZA,ZB,CA,CB,SQNrm
         INTEGER                     :: K,L,M,N,NC,NK,NP,LMN,MinL,MaxL,PFA,PFB
         REAL(DOUBLE),PARAMETER, &    ! Fact=Sqrt[Pi](2*L-1)!! 2^(-L)
               DIMENSION(0:10)       :: Fact=(/0.17724538509055160273D1, &
                                               0.8862269254527580136D0,  &
                                               0.13293403881791370205D1, &
                                               0.3323350970447842551D1,  &
                                               0.11631728396567448929D2, &
                                               0.5234277778455352018D2,  &
                                               0.28788527781504436100D3, &
                                               0.1871254305797788347D4,  &
                                               0.14034407293483412599D5, &
                                               0.1192924619946090071D6,  &
                                               0.11332783889487855673D7/) 
         B%NKind=A%NKind
         B%NAtms=A%NAtms
         B%NPrim=A%NPrim
         B%NCtrt=A%NCtrt
         B%NAsym=A%NAsym
         B%LMNLen=LHGTF(B%NAsym)
         CALL New(B)
         DO NK=1,B%NKind
            B%Kinds%I(NK)=A%Kinds%I(NK)
            B%NCFnc%I(NK)=A%NCFnc%I(NK)
            DO NC=1,B%NCFnc%I(NK)
               B%NPFnc%I(NC,NK)=A%NPFnc%I(NC,NK)
               B%ASymm%I(1:2,NC,NK)=A%ASymm%I(1:2,NC,NK)
               DO NP=1,B%NPFnc%I(NC,NK)
                  B%Expnt%D(NP,NC,NK)=A%Expnt%D(NP,NC,NK)
               ENDDO
            ENDDO
         ENDDO
         CALL BSetIndx(B)
         DO K=1,B%NKind
            DO NC=1,B%NCFnc%I(K)
               MinL=B%ASymm%I(1,NC,K)
               MaxL=B%ASymm%I(2,NC,K) 
               DO NP=1,B%NPFnc%I(NC,K)
                  Z=A%Expnt%D(NP,NC,K)
                  DO L=MinL,MaxL
                     C=A%CCoef%D(L+1,NP,NC,K)
                     A%CCoef%D(L+1,NP,NC,K)=C*Z**(Half*DBLE(L)+0.75D0)
                  ENDDO
               ENDDO
               DO LMN=B%LStrt%I(NC,K),B%LStop%I(NC,K)            
                  L=B%LxDex%I(LMN)
                  M=B%LyDex%I(LMN)
                  N=B%LzDex%I(LMN)
                  Expnt=-(DBLE(L+M+N)+1.5D0)
                  RNorm=0.0D0
                  DO PFA=1,B%NPFnc%I(NC,K)
                     ZA=B%Expnt%D(PFA,NC,K)
                     CA=A%CCoef%D(L+M+N+1,PFA,NC,K)
                     DO PFB=1,B%NPFnc%I(NC,K)
                        ZB=B%Expnt%D(PFB,NC,K)
                        CB=A%CCoef%D(L+M+N+1,PFB,NC,K)
                        RNorm=RNorm+CA*CB*(ZA+ZB)**Expnt
                     ENDDO
                  ENDDO
                  IF(RNorm==Zero) &
                     CALL MondoHalt(-999,' Zero integrals durring basis set normilizaiton! ')
                  RNorm=RNorm*Fact(L)*Fact(M)*Fact(N)
                  SqNrm=1.0D0/DSQRT(RNorm)
                  DO NP=1,B%NPFnc%I(NC,K)
                     B%CCoef%D(LMN,NP,NC,K)=A%CCoef%D(L+M+N+1,NP,NC,K)*SqNrm
                  ENDDO
              ENDDO 
            ENDDO
         ENDDO 
      END SUBROUTINE NormalizeBaseSets

      SUBROUTINE CalcDimensions(BS,BSiz_2,OffS_2)
         TYPE(BSET),     INTENT(INOUT):: BS
         TYPE(INT_VECT), INTENT(OUT)  :: BSiz_2,OffS_2
         TYPE(CRDS)                   :: GM
         INTEGER                      :: NA,NK,NC,Stride
!----------------------------------------------------------------------------
!        Get the geometry and compute matrix dimensions and bloking indecies
!
         CALL Get(GM,Tag_O='1')
         CALL New(BSiz_2,GM%NAtms)
         CALL New(OffS_2,GM%NAtms)
         BS%NBasF=0
         OffS_2%I(1)=1
         DO NA=1,GM%NAtms
            BSiz_2%I(NA)=0
            NK=GM%AtTyp%I(NA)
            IF(BS%NCFnc%I(NK)==0)  &
               CALL MondoHalt(PRSE_ERROR,'BS%NCFnc%I(NK)=0 in CalcDimensions.')
            DO NC=1,BS%NCFnc%I(NK)
               Stride=BS%LStop%I(NC,NK)-BS%LStrt%I(NC,NK)+1
               BS%NBasF=BS%NBasF+Stride
               BSiz_2%I(NA)=BSiz_2%I(NA)+Stride
            ENDDO
            IF(NA.GE.2)OffS_2%I(NA)=OffS_2%I(NA-1)+BSiz_2%I(NA-1)         
         ENDDO
         CALL Delete(GM)
       END SUBROUTINE CalcDimensions
!
!----------------------------------------------------------------------------
!
#ifdef MMech
      SUBROUTINE ParseMM(Ctrl)
        Implicit None
        TYPE(SCFControls)          :: Ctrl
        TYPE(CRDS)                 :: GM_MM,GM
        REAL(DOUBLE)               :: sum,dx,dy,dz,dd
        CHARACTER(LEN=DEFAULT_CHR_LEN) :: OPLS_DATABASE,MM_COORDS,MM_SEQ
        CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Line, LineLowCase
        CHARACTER(LEN=DEFAULT_CHR_LEN)  :: QMMM_DATA
        TYPE(INT_VECT) :: ATMMARK,ACTIVE_BOND,ACTIVE_ANGLE
        TYPE(INT_VECT) :: ACTIVE_DIHEDRAL,ACTIVE_IMPROPER
        INTEGER :: I,J,K,L,N1,N2,N3,N4
        TYPE(DBL_VECT) :: CHARGE14,LJEPS14,LJEPS,LJRAD
!
         CALL OpenASCII(InpFile,Inp)
         CALL OpenASCII(OutFile,Out)
!  
! Find names of datafiles for MM calculation
!
         CALL AlignLowCase('begin_mm_filenames',Inp)
           READ(Inp,DEFAULT_CHR_FMT) OPLS_DATABASE
           READ(Inp,DEFAULT_CHR_FMT) MM_SEQ
           READ(Inp,DEFAULT_CHR_FMT) MM_COORDS
           IF(Ctrl%Mechanics(1).AND.Ctrl%Mechanics(2)) &
              READ(Inp,DEFAULT_CHR_FMT) QMMM_DATA
         DO 
           READ(Inp,DEFAULT_CHR_FMT,END=1)Line
           LineLowCase = Line
           Call LowCase(LineLowCase)
           IF(INDEX(LineLowCase,'end_mm_filenames')/=0) GO TO 2
         ENDDO
      1  CALL MondoHalt(PRSE_ERROR,' Found no <end_mm_filenames> in input file '//TRIM(InpFile))
      2  CONTINUE
!
! Check whether filenames are given 
!
         N1 = LEN_TRIM(OPLS_DATABASE) 
         N2 = LEN_TRIM(MM_SEQ)
         N3 = LEN_TRIM(MM_COORDS)
         N4 = LEN_TRIM(QMMM_DATA)
         IF(N1 == 0) CALL MondoHalt(PRSE_ERROR,' Found no <opls_database> in input file '//TRIM(InpFile))
         IF(N2 == 0) CALL MondoHalt(PRSE_ERROR,' Found no <mm_seq> in input file '//TRIM(InpFile))
         IF(N3 == 0) CALL MondoHalt(PRSE_ERROR,' Found no <mm_coords> in input file '//TRIM(InpFile))
         IF(Ctrl%Mechanics(1).AND.Ctrl%Mechanics(2)) THEN
         IF(N4 == 0) CALL MondoHalt(PRSE_ERROR,' Found no <qmmmdata> in input file '//TRIM(InpFile))
         ENDIF
!
! Now process MM data files
!
      CALL MM_FILE_PROCESS ( "tmp_opls_bin", OPLS_DATABASE(1:N1))
      CALL MM_SYSTEM_CONSTRUCT ( "tmp_opls_bin", MM_SEQ(1:N2) )
      CALL COORDINATES_READ ( MM_COORDS(1:N3) )
!
!     write(*,*) 'bef sort'
!     CALL SORT_INTO_BOX(7.D0,ATMCRD,MM_NATOMS) !!!
!     stop
      CALL NEW(ATMMARK,MM_NATOMS)
      CALL NEW(ACTIVE_BOND,NBONDS)
      CALL NEW(ACTIVE_ANGLE,NANGLES)
      CALL NEW(ACTIVE_DIHEDRAL,NDIHEDRALS)
      CALL NEW(ACTIVE_IMPROPER,NIMPROPERS)
!
      IF(Ctrl%Mechanics(1).AND.Ctrl%Mechanics(2)) THEN
!
! Calling QMMMDATA_READ also modifies ATMCHG (MM charges)
!
        CALL QMMMDATA_READ(QMMM_DATA,ATMMARK,ACTIVE_BOND,ACTIVE_ANGLE, &
        ACTIVE_DIHEDRAL,ACTIVE_IMPROPER)
!
      ELSE
!
        ATMMARK%I(1:MM_NATOMS)=1
        ACTIVE_BOND%I(1:NBONDS)=1
        ACTIVE_ANGLE%I(1:NANGLES)=1
        ACTIVE_DIHEDRAL%I(1:NDIHEDRALS)=1
        ACTIVE_IMPROPER%I(1:NIMPROPERS)=1
!
      ENDIF
!
! SAVE 14 charges and Lennard-Jones parameters
! for the calculation of exclusion energy terms
!
      CALL NEW(CHARGE14,MM_NATOMS)
      CALL NEW(LJEPS14,MM_NATOMS)
      CALL NEW(LJEPS,MM_NATOMS)
      CALL NEW(LJRAD,MM_NATOMS)
        CHARGE14%D(:)=ATMCHG14(:)
        LJEPS14%D(:)=ATMEPS14(:)
        LJEPS%D(:)=ATMEPS(:)
        LJRAD%D(:)=ATMSIG(:)
      CALL OpenHDF(Ctrl%Info)
        CALL PUT(CHARGE14,'CHARGE14')
        CALL PUT(LJEPS14,'LJEPS14')
        CALL PUT(LJEPS,'LJEPS')
        CALL PUT(LJRAD,'LJRAD')
      CALL CloseHDF
      CALL DELETE(CHARGE14)
      CALL DELETE(LJEPS14)
      CALL DELETE(LJEPS)
      CALL DELETE(LJRAD)
!
! Calculate topology information (currently only 
! for MM and QMMM calculations available)
!
      IF(Ctrl%Mechanics(1)) THEN
        CALL TOPOLOGIES_MM(MM_NATOMS,NBONDS, &
             BONDS(:)%I,BONDS(:)%J,Ctrl%Info)
      ENDIF
!
! Calculate total MM charge
!
        SUM=Zero
      DO I=1,MM_NATOMS
        SUM=SUM+ATMCHG(I)
      ENDDO
        GM_MM%TotCh=SUM
        write(*,*) 'mm total charge ',GM_MM%TotCh
!
      GM_MM%NAtms = MM_Natoms
       CALL New(GM_MM)
      GM_MM%InAu = .true.
      GM_MM%Nkind = NTYPES
      GM_MM%Carts%D(:,:) = ATMCRD(:,:)*AngstromsToAU
      GM_MM%AtMss%D(:) = ATMMAS(:)
      GM_MM%AtNum%D(:) = ATMCHG(:)
      IF(HasQm()) Then
        CALL Get(GM,Tag_O=IntToChar(Ctrl%NGeom))
        GM_MM%BndBox%D=SetBox(GM_MM%Carts,GM%Carts)
      ELSE
        GM_MM%BndBox%D=SetBox(GM_MM%Carts)
      ENDIF
!
! Print out MM coorinates into output file
!
        WRITE(Out,*) 
        WRITE(Out,*) 'MM system coordinates:'
        WRITE(Out,*) 
      DO I=1,MM_Natoms
        WRITE(Out,110) I,ATMNAM(I),ATMCRD(1:3,I)
      ENDDO
110   FORMAT(I7,2X,A8,3F12.6)
!
      CALL OpenHDF(Ctrl%Info)
        CALL Put(GM_MM,Tag_O='GM_MM'//'1')
!       IF(Ctrl%Mechanics(1).AND.Ctrl%Mechanics(2)) THEN
          IF(MM_NATOMS/=0) THEN
            CALL Put(ATMMARK,'ATMMARK',N_O=MM_NATOMS)
          ENDIF
          IF(NBONDS/=0) CALL Put(ACTIVE_BOND,'ACTIVE_BOND',N_O=NBONDS)
          IF(NANGLES/=0) CALL Put(ACTIVE_ANGLE,'ACTIVE_ANGLE',N_O=NANGLES)
          IF(NDIHEDRALS/=0) CALL Put(ACTIVE_DIHEDRAL,'ACTIVE_DIHEDRAL',N_O=NDIHEDRALS)
          IF(NIMPROPERS/=0) CALL Put(ACTIVE_IMPROPER,'ACTIVE_IMPROPER',N_O=NIMPROPERS)
          CALL DELETE(ATMMARK)
          CALL DELETE(ACTIVE_BOND)
          CALL DELETE(ACTIVE_ANGLE)
          CALL DELETE(ACTIVE_DIHEDRAL)
          CALL DELETE(ACTIVE_IMPROPER)
!       ENDIF
      CALL CloseHDF()
!
!! test coulomb energy
!!
!        SUM=zero
!      Do i=1,MM_NATOMS
!      Do j=i+1,MM_NATOMS
!        DX=GM_MM%Carts%D(1,i)-GM_MM%Carts%D(1,j)
!        DY=GM_MM%Carts%D(2,i)-GM_MM%Carts%D(2,j)
!        DZ=GM_MM%Carts%D(3,i)-GM_MM%Carts%D(3,j)
!        DD=SQRT(DX**2+DY**2+DZ**2) !!!!*AngstromsToAU  
!        SUM=SUM+GM_MM%AtNum%D(i)*GM_MM%AtNum%D(j)/DD
!      enddo
!      enddo
!        write(*,*) 'coulomb energy= ',sum 
!        write(out,*) 'coulomb energy= ',sum 
!
      CALL DELETE(GM_MM)
!
      CLOSE(UNIT=Out,STATUS='KEEP')
!
      END SUBROUTINE ParseMM
#endif
!
!------------------------------------------------------------------------
!
#ifdef MMech
!
      SUBROUTINE ParseMech(Ctrl)
         TYPE(SCFControls)          :: Ctrl
!----------------------------------------------------------------------------
!        Open info file
         CALL OpenHDF(InfFile)
!        Open input file
         CALL OpenASCII(InpFile,Inp) 
!
!        Parse <OPTIONS> for QM_MM
!
         CALL OpenASCII(InpFile,Inp)
!
         Ctrl%Mechanics(1) = .false.
         Ctrl%Mechanics(2) = .true.  !!!! default is QM 
         IF(    OptKeyQ(Inp,mechanics_option,pureMM))THEN
            Ctrl%Mechanics(1) = .true.
            Ctrl%Mechanics(2) = .false.
         ELSEIF(OptKeyQ(Inp,mechanics_option,pureQM))THEN
            Ctrl%Mechanics(1) = .false.
            Ctrl%Mechanics(2) = .true.
         ELSEIF(OptKeyQ(Inp,mechanics_option,QMandMM))THEN
            Ctrl%Mechanics(1) = .true.
            Ctrl%Mechanics(2) = .true.
         ENDIF
!
        CALL OpenHDF(InfFile)
        CALL Put(Ctrl%Mechanics(1),'Ctrl_Mechanics1')
        CALL Put(Ctrl%Mechanics(2),'Ctrl_Mechanics2')
        CALL CLOSEHDF()
!
        CLOSE(UNIT=Inp,STATUS='KEEP')
!
      END SUBROUTINE ParseMech
#endif
!---------------------------------------------------------------------------- 
      SUBROUTINE ParseGrad(Ctrl)
         TYPE(SCFControls)          :: Ctrl
         TYPE(TOLS)                 :: Thrsh ! Thresholds
!----------------------------------------------------------------------------
!        Open info file
         CALL OpenHDF(InfFile)
!        Open input file
         CALL OpenASCII(InpFile,Inp) 
!        Parse <OPTIONS> for <Grad=>
!     
         IF(OptKeyQ(Inp,GRADIENTS,FORCE))THEN
            Ctrl%Grad=GRAD_ONE_FORCE
            Ctrl%NGeom=1
         ELSEIF(OptKeyQ(Inp,DYNAMICS,MD_VERLET))THEN
            Ctrl%Grad=GRAD_MD
            Ctrl%MDC%MD_Algor=1
            IF(.NOT. OptIntQ(Inp,MAX_STEPS,Ctrl%NGeom)) THEN
               Ctrl%NGeom=1
            ENDIF
            IF(.NOT. OptDblQ(Inp,MD_TIME_STEP,Ctrl%MDC%TimeStep)) THEN
               Ctrl%MDC%TimeStep    = 1.0D-2
            ENDIF
            IF(.NOT. OptDblQ(Inp,MD_VEL_SCALE,Ctrl%MDC%VelScaling)) THEN
               Ctrl%MDC%VelScaling  = One
            ENDIF
            IF(.NOT. OptDblQ(Inp,MD_TMP_SCALE,Ctrl%MDC%TempScaling)) THEN
               Ctrl%MDC%TempScaling = One
            ENDIF
         ELSEIF(OptKeyQ(Inp,DYNAMICS,MD_PRECOR))THEN
            Ctrl%Grad=GRAD_MD
            Ctrl%MDC%MD_Algor=2
            CALL MondoHalt(PRSE_ERROR,'Predictor-Corrector Algorithmn not implimented')
         ELSEIF(OptKeyQ(Inp,OPTIMIZATION,OPT_QUNEW))THEN
            IF(OptKeyQ(Inp,OPTIMIZATION,OPT_ONE_BASE))THEN
               Ctrl%Grad=GRAD_QNEW_ONE_OPT
            ELSE
               Ctrl%Grad=GRAD_QNEW_OPT
            ENDIF
            IF(.NOT.OptIntQ(Inp,MAX_STEPS,Ctrl%NGeom))THEN
               Ctrl%NGeom=500
            ENDIF
         ELSEIF(OptKeyQ(Inp,OPTIMIZATION,OPT_TSTATE))THEN
            Ctrl%Grad=GRAD_TS_SEARCH
            IF(.NOT.OptIntQ(Inp,MAX_STEPS,Ctrl%NGeom))THEN
               Ctrl%NGeom=500
            ENDIF
         ELSE
            Ctrl%Grad=GRAD_NO_GRAD
         ENDIF
         IF(Ctrl%Grad>=GRAD_QNEW_OPT)THEN
            IF(.NOT.OptIntQ(Inp,MAX_STEPS,Ctrl%NGeom))THEN
               Ctrl%NGeom=500
            ENDIF
         ENDIF
        CALL CLOSEHDF()
        CLOSE(UNIT=Inp,STATUS='KEEP')
!
     END SUBROUTINE PARSEGRAD
!
!----------------------------------------------------------------------------
#ifdef MMech
!
      SUBROUTINE QMMMDATA_READ(FILE,ATMMARK,ACTIVE_BOND,ACTIVE_ANGLE, &
        ACTIVE_DIHEDRAL,ACTIVE_IMPROPER)
      IMPLICIT NONE
      INTEGER III,K,L,M,N,N4,ISUB,IRES,II,J
      CHARACTER(LEN=DEFAULT_CHR_LEN)  :: FILE
      TYPE(INT_VECT) :: ATMMARK,ACTIVE_BOND,ACTIVE_ANGLE
      TYPE(INT_VECT) :: ACTIVE_DIHEDRAL,ACTIVE_IMPROPER
!
      N4 = LEN_TRIM(FILE)
!
! Read QM-MM distinquishing list and modified charges
! which must be zero on all QM atoms and 'scaled' for the MM ones.
!
      OPEN(UNIT=60,FILE=FILE(1:N4),STATUS='OLD')
      DO ISUB=1,NSUBSYS
        DO IRES=SUBIND(ISUB)+1,SUBIND(ISUB+1)
            READ(60,750) 
          DO I=RESIND(IRES)+1,RESIND(IRES+1)
            READ(60,700) ATMNAM(I),ATMMARK%I(I),ATMCHG(I)
          ENDDO
        ENDDO
      ENDDO
!
600   FORMAT(I10,3F15.7)
700   FORMAT(10X,3X,A8,I4,2F20.10)
750   FORMAT(I10,A8,F20.10)
!
! bonds
!
         READ(60,*) 
      DO I=1,NBONDS
       READ(60,61) II,K,L,ACTIVE_BOND%I(I)
      ENDDO
61    Format(I10,3X,2I5,3X,I5)
!
! angles
!
         READ(60,*) 
      DO II=1,NANGLES
        READ(60,62) III,I,J,K,ACTIVE_ANGLE%I(II)
      ENDDO
62    Format(I10,3X,3I5,3X,I5)
!
! dihedrals
!
         READ(60,*) 
      DO II=1,NDIHEDRALS
        READ(60,81) III,I,J,K,L,ACTIVE_DIHEDRAL%I(II)
      ENDDO
81  FORMAT(I10,3X,4I5,3X,I5,F12.4)
!
! impropers
!
        READ(60,*) 
      DO II=1,NIMPROPERS
        READ(60,81) III,I,J,K,L,ACTIVE_IMPROPER%I(II)
      ENDDO
!
      CLOSE(60)
!
      END SUBROUTINE QMMMDATA_READ
#endif
!============================================================================
!     Parse Accuracy for MM
!============================================================================
      SUBROUTINE ParseAcc(Ctrl)
         TYPE(SCFControls)          :: Ctrl
         TYPE(TOLS)                 :: Thrsh ! Thresholds
!----------------------------------------------------------------------------
!        Open info file
         CALL OpenHDF(InfFile)
!        Open input file
         CALL OpenASCII(InpFile,Inp) 
!
!        Parse <OPTIONS.ACCURACY> 
!
         Ctrl%AccL=2 ! Default is "good"    
         NOpts=0
!
         IF(OptKeyLocQ(Inp,ACCURACY_OPTION,ACCURACY_CHEEZY,MaxSets,NLoc,Loc))THEN
            NOpts=NOpts+NLoc
            DO ILoc=1,NLoc; Ctrl%AccL(Loc(ILoc))=1; ENDDO
         ENDIF
         IF(OptKeyLocQ(Inp,ACCURACY_OPTION,ACCURACY_GOOD,MaxSets,NLoc,Loc))THEN
            NOpts=NOpts+NLoc
            DO ILoc=1,NLoc; Ctrl%AccL(Loc(ILoc))=2; ENDDO
         ENDIF
         IF(OptKeyLocQ(Inp,ACCURACY_OPTION,ACCURACY_TIGHT,MaxSets,NLoc,Loc))THEN
            NOpts=NOpts+NLoc
            DO ILoc=1,NLoc; Ctrl%AccL(Loc(ILoc))=3; ENDDO
         ENDIF
         IF(OptKeyLocQ(Inp,ACCURACY_OPTION,ACCURACY_RETENTIVE,MaxSets,NLoc,Loc))THEN
            NOpts=NOpts+NLoc
            DO ILoc=1,NLoc; Ctrl%AccL(Loc(ILoc))=4; ENDDO
         ENDIF
#ifdef MMech
         IF(Ctrl%Mechanics(2)) THEN
#endif        
           IF(NOpts>1.AND.NOpts/=Ctrl%NSet) &
              CALL MondoHalt(PRSE_ERROR,'Number of '//ACCURACY_OPTION &
                           //' options does not match number of Basis sets.')
#ifdef MMech
         ENDIF
#endif        
!        Set thresholds
!
#ifdef MMech
!
!!! MMonly
!
         IF(Ctrl%Mechanics(1).AND..NOT.Ctrl%Mechanics(2)) THEN
           Ctrl%NSet=1 !!!!!!!!!!!!!TEMPORARY!!!!
           Thrsh%Cube=CubeNeglect(Ctrl%AccL(1))
           Thrsh%Trix=TrixNeglect(Ctrl%AccL(1))
           Thrsh%Dist=DistNeglect(Ctrl%AccL(1))
           Thrsh%TwoE=TwoENeglect(Ctrl%AccL(1))
           Thrsh%ETol=ETol(Ctrl%AccL(1))
           Thrsh%DTol=DTol(Ctrl%AccL(1))
           CALL Put(Thrsh,Tag_O=IntToChar(1))
         ELSE
#endif
!    
!!! not MMonly
         IF(NOpts==Ctrl%NSet)THEN
!           All thresholds are user determined
            DO I=1,Ctrl%NSet
               Thrsh%Cube=CubeNeglect(Ctrl%AccL(I))
               Thrsh%Trix=TrixNeglect(Ctrl%AccL(I))
               Thrsh%Dist=DistNeglect(Ctrl%AccL(I))
               Thrsh%TwoE=TwoENeglect(Ctrl%AccL(I))
               Thrsh%ETol=ETol(Ctrl%AccL(I))
               Thrsh%DTol=DTol(Ctrl%AccL(I))
               CALL Put(Thrsh,Tag_O=IntToChar(I))
            ENDDO
         ELSE ! Default, cheezy for all sets except last, which is good or user defined.
            DO I=1,Ctrl%NSet-1
               Thrsh%Cube=CubeNeglect(1)
               Thrsh%Trix=TrixNeglect(1)
               Thrsh%Dist=DistNeglect(1)
               Thrsh%TwoE=TwoENeglect(1)
               Thrsh%ETol=ETol(1)
               Thrsh%DTol=DTol(1)
               CALL Put(Thrsh,Tag_O=IntToChar(I))
            ENDDO
            Thrsh%Cube=CubeNeglect(Ctrl%AccL(1))
            Thrsh%Trix=TrixNeglect(Ctrl%AccL(1))
            Thrsh%Dist=DistNeglect(Ctrl%AccL(1))
            Thrsh%TwoE=TwoENeglect(Ctrl%AccL(1))
            Thrsh%ETol=ETol(Ctrl%AccL(1))
            Thrsh%DTol=DTol(Ctrl%AccL(1))
            CALL Put(Thrsh,Tag_O=IntToChar(Ctrl%NSet))
         ENDIF
!
#ifdef MMech
         ENDIF !!! MMonly
#endif
!
!        Close files
         CLOSE(UNIT=Inp,STATUS='KEEP')
         CALL CloseHDF()
      END SUBROUTINE ParseAcc
!
END MODULE
