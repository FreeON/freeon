MODULE DrvSCFs
  USE DerivedTypes
  USE GlobalScalars
#ifdef NAG
  USE F90_UNIX_PROC
#endif
  USE PrettyPrint
  USE SCFLocals
  USE Overlay
  USE ParsingKeys
  USE Macros
  USE Functionals
  IMPLICIT NONE 
!========================================================================================
!
!========================================================================================
  INTEGER,               PARAMETER :: DCL=DEFAULT_CHR_LEN
  CHARACTER(LEN=DCL),    PARAMETER :: False  = '.FALSE.'   ! False flag
  CHARACTER(LEN=DCL),    PARAMETER :: True   = '.TRUE.'    ! True flag
  CHARACTER(LEN=DCL),    PARAMETER :: Direct = 'Direct'    ! Do direct SCF
  CHARACTER(LEN=DCL),    PARAMETER :: Restart= 'Restart'   ! Restart from InfoFile
  CHARACTER(LEN=DCL),    PARAMETER :: Switch = 'Switch'    ! Do a basis set switch
  CHARACTER(LEN=DCL),    PARAMETER :: InkFok = 'InkFok'    ! Do incremental Fock build
  CHARACTER(LEN=DCL),    PARAMETER :: Core   = 'Core'      ! Use core Hamiltontian guess
  CHARACTER(LEN=DCL),    PARAMETER :: Minimal= 'Minimal'   ! Use superposition of AO DMs guess
  CHARACTER(LEN=DCL),    PARAMETER :: Blank  =' '
  CHARACTER(LEN=DCL),DIMENSION(10) :: CtrlVect
  CONTAINS
!========================================================================================
!
!========================================================================================
     SUBROUTINE SCFCycle(Ctrl)
        TYPE(SCFControls)  :: Ctrl
        INTEGER            :: ICyc,IGeo
!----------------------------------------------------------------------------------------
        ICyc=Ctrl%Current%I(1)
        IGeo=Ctrl%Current%I(2)
!       Hack ...
        Ctrl%DIIS=.TRUE.
        Ctrl%InkFok=.FALSE.
!       Do an SCF cycle 
        IF(ICyc>0)THEN
!          Restricted direct SCF cycle
           CALL DirectRSCF(Ctrl)
        ELSEIF(ICyc==0.AND.Ctrl%Rest.OR.(IGeo==1.AND.Ctrl%SuperP))THEN
!          Restricted SCF start up from guess or previous density 
           CALL SuperRSCF(Ctrl)
        ELSEIF(ICyc==0.AND.IGeo==1)THEN
!          Resctricted direct SCF from core Hamiltonian
           CALL CoreRSCF(Ctrl)
        ELSEIF(ICyc==0.AND.IGeo>1)THEN
!          Restricted SCF with basis set switch 
           CALL SwitchRSCF(Ctrl)
        ELSE
           WRITE(*,*)' Previous state = ',Ctrl%Previous%I
           WRITE(*,*)' Current  state = ',Ctrl%Current%I
           CALL MondoHalt(DRIV_ERROR,' bad option in SCFCycle ')
        ENDIF
   END SUBROUTINE
!========================================================================================
!
!========================================================================================
     FUNCTION SetCtrlVect(Ctrl,Actn1_O,Actn2_O) RESULT(CVect)
        TYPE(SCFControls),        INTENT(IN) :: Ctrl
        CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: Actn1_O,Actn2_O
        CHARACTER(LEN=DCL),DIMENSION(10)     :: CVect
        CVect(1)=Ctrl%Name
        CVect(2)=" "
        CVect(3)=" "
        IF(PRESENT(Actn1_O))CVect(2)=Actn1_O
        IF(PRESENT(Actn2_O))CVect(3)=Actn2_O
        CVect(4:10)=(/(IntToChar(Ctrl%Current%I(1))), &
                      (IntToChar(Ctrl%Current%I(2))), &
                      (IntToChar(Ctrl%Current%I(3))), &
                      (IntToChar(Ctrl%Previous%I(1))),&
                      (IntToChar(Ctrl%Previous%I(2))),&
                      (IntToChar(Ctrl%Previous%I(3))),&
                      (IntToChar(Ctrl%NCGC))/)              
     END FUNCTION SetCtrlVect
!========================================================================================
!    Perform a direct resctricted SCF cycle
!========================================================================================
     SUBROUTINE DirectRSCF(Ctrl)
        TYPE(SCFControls)  :: Ctrl
!---------------------------------------------------------------------------------------
        CALL LogSCF(Ctrl%Current,'Restricted SCF with DIIS.')
        CtrlVect=SetCtrlVect(Ctrl,Direct)
        CALL Invoke('MakeRho',     CtrlVect)
        CALL Invoke('QCTC',      CtrlVect)
        IF(HasDFT(Ctrl%Model(Ctrl%Current%I(2)))) &
           CALL Invoke('HiCu',     CtrlVect)
        IF(HasHF(Ctrl%Model(Ctrl%Current%I(2))))&
           CALL Invoke('ONX',     CtrlVect)
        CALL Invoke('FBuild',     CtrlVect,MPIRun_O=.TRUE.)
        IF(Ctrl%DIIS) &
           CALL Invoke('DIIS',     CtrlVect,MPIRun_O=.TRUE.)
        IF(Ctrl%Method(Ctrl%Current%I(2))==RH_R_SCF)THEN
           CALL Invoke('RHeqs',    CtrlVect)
        ELSEIF(Ctrl%Method(Ctrl%Current%I(2))==SDMM_R_SCF)THEN
           CALL Invoke('SDMM',     CtrlVect,MPIRun_O=.TRUE.)                                     
        ENDIF
        CALL Invoke('SCFstats',   CtrlVect,MPIRun_O=.TRUE.)                                     
     END SUBROUTINE DirectRSCF
!========================================================================================
!    Perform a direct resctricted SCF cycle
!========================================================================================
     SUBROUTINE CoreRSCF(Ctrl)
        TYPE(SCFControls)  :: Ctrl
!---------------------------------------------------------------------------------------
        CALL LogSCF(Ctrl%Current,'Guess is Core Hamiltonian')
        CtrlVect=SetCtrlVect(Ctrl,Core)
        CALL Invoke('MakeRho',     CtrlVect)
        CALL Invoke('QCTC',      CtrlVect)
        CALL Invoke('FBuild',     CtrlVect,MPIRun_O=.TRUE.)
!
!        CALL MondoHalt(-100,'STOPPING FOR Fock in Core')
!
        IF(Ctrl%Method(Ctrl%Current%I(2))==RH_R_SCF)THEN
           CALL Invoke('RHeqs',    CtrlVect)
        ELSEIF(Ctrl%Method(Ctrl%Current%I(2))==SDMM_R_SCF)THEN
           CALL Invoke('SDMM',     CtrlVect,MPIRun_O=.TRUE.)
        ENDIF
     END SUBROUTINE CoreRSCF
!========================================================================================
!    Perform a direct resctricted SCF cycle starting from a density that is an 
!    atomic superposition of guess density 
!========================================================================================
     SUBROUTINE SuperRSCF(Ctrl)
        TYPE(SCFControls)  :: Ctrl
!---------------------------------------------------------------------------------------
        CALL LogSCF(Ctrl%Current,'Guess is superposition of atomic density matrices')
        IF(Ctrl%Rest)THEN
           CtrlVect=SetCtrlVect(Ctrl,Restart)
           Ctrl%Rest=.FALSE.
        ELSE
           CtrlVect=SetCtrlVect(Ctrl,Direct)
        ENDIF
        CALL Invoke('P2Use',       CtrlVect,MPIRun_O=.TRUE.)
        CtrlVect=SetCtrlVect(Ctrl,Direct)
        CALL Invoke('MakeRho',     CtrlVect)
        CALL Invoke('QCTC',      CtrlVect)
        IF(HasDFT(Ctrl%Model(Ctrl%Current%I(2)))) &
           CALL Invoke('HiCu',     CtrlVect)
        IF(HasHF(Ctrl%Model(Ctrl%Current%I(2))))  &
           CALL Invoke('ONX',     CtrlVect)
        CALL Invoke('FBuild',     CtrlVect,MPIRun_O=.TRUE.)
!
!        CALL MondoHalt(-100,'STOPPING FOR Fock in Super:Old')
!
        IF(Ctrl%Method(Ctrl%Current%I(2))==RH_R_SCF)THEN
           CALL Invoke('RHeqs',    CtrlVect)
        ELSEIF(Ctrl%Method(Ctrl%Current%I(2))==SDMM_R_SCF)THEN
           CALL Invoke('SDMM',     CtrlVect,MPIRun_O=.TRUE.)                                     
        ENDIF
        CALL Invoke('SCFstats',   CtrlVect,MPIRun_O=.TRUE.)                                     
     END SUBROUTINE SuperRSCF
!========================================================================================
!    Perform a resctricted SCF cycle with basis set switch
!========================================================================================
     SUBROUTINE SwitchRSCF(Ctrl)
        TYPE(SCFControls)  :: Ctrl
        CHARACTER(LEN=DCL),DIMENSION(2) :: ShortCtrlVect        
        INTEGER,DIMENSION(3)            :: NextStats
!---------------------------------------------------------------------------------------
        CALL LogSCF(Ctrl%Current,'Restricted SCF start via basis set switch')
!       Build new rho from last density matrix and rename it, unless same basis different
!       geometry
        IF(Ctrl%Current%I(3)==Ctrl%Previous%I(3))THEN
           CtrlVect=SetCtrlVect(Ctrl,Switch)
           CtrlVect(3)=TRIM(IntToChar(Ctrl%Previous%I(1)+1))
           CtrlVect(4)=TRIM(IntToChar(Ctrl%Previous%I(2)))
           CtrlVect(5)=TRIM(IntToChar(Ctrl%Previous%I(3)))
           CALL Invoke('MakeRho',CtrlVect)
           ShortCtrlVect(1)=TrixFile('Rho',OffSet_O=1,Name_O=Ctrl%Name,Stats_O=Ctrl%Previous%I)
           ShortCtrlVect(2)=TrixFile('Rho',OffSet_O=0,Name_O=Ctrl%Name,Stats_O=Ctrl%Current%I )
           CALL Invoke('/bin/cp',ShortCtrlVect,AbsPath_O=.TRUE.)
        ENDIF
!
        CtrlVect=SetCtrlVect(Ctrl,Switch)
        ShortCtrlVect(1)=TrixFile('D',OffSet_O=1,Name_O=Ctrl%Name,Stats_O=Ctrl%Previous%I)
        ShortCtrlVect(2)=TrixFile('D',OffSet_O=0,Name_O=Ctrl%Name,Stats_O=Ctrl%Current%I )
        CALL Invoke('/bin/cp',ShortCtrlVect,AbsPath_O=.TRUE.)
        ShortCtrlVect(1)=TrixFile('OrthoD',OffSet_O=1,Name_O=Ctrl%Name,Stats_O=Ctrl%Previous%I)
        ShortCtrlVect(2)=TrixFile('OrthoD',OffSet_O=0,Name_O=Ctrl%Name,Stats_O=Ctrl%Current%I )
        CALL Invoke('/bin/cp',ShortCtrlVect,AbsPath_O=.TRUE.)

!       If same basis, different geometry do DM extrapolation 

        IF(Ctrl%Current%I(3)/=Ctrl%Previous%I(3))THEN
           CALL Invoke('P2Use',    CtrlVect,MPIRun_O=.TRUE.)
           CALL Invoke('MakeRho',     CtrlVect)
        ENDIF
        CALL Invoke('QCTC' ,     CtrlVect)
        IF(HasDFT(Ctrl%Model(Ctrl%Current%I(2))))  &
           CALL Invoke('HiCu',     CtrlVect)
        IF(HasHF(Ctrl%Model(Ctrl%Current%I(2))))   &
           CALL Invoke('ONX',     CtrlVect)
        CALL Invoke('FBuild',     CtrlVect,MPIRun_O=.TRUE.)
        IF(Ctrl%Method(Ctrl%Current%I(2))==RH_R_SCF)THEN
           CALL Invoke('RHeqs',    CtrlVect)
        ELSEIF(Ctrl%Method(Ctrl%Current%I(2))==SDMM_R_SCF)THEN
           CALL Invoke('SDMM',     CtrlVect,MPIRun_O=.TRUE.)
        ENDIF
        CALL Invoke('SCFstats',   CtrlVect,MPIRun_O=.TRUE.)                                     
     END SUBROUTINE SwitchRSCF
!========================================================================================
!
!========================================================================================
     SUBROUTINE OneEMats(Ctrl)
        TYPE(SCFControls)             :: Ctrl
        INTEGER                       :: ISet,I
        CtrlVect=SetCtrlVect(Ctrl,Core)
!        IF(PrintFlags%Int==DEBUG_INTEGRAL) &
!           CALL Invoke('PrintInt', CtrlVect)
        CALL Invoke('MakeS',        CtrlVect)
        IF(Ctrl%Method(Ctrl%Current%I(2))==RH_R_SCF)THEN
           CALL Invoke('LowdinO',  CtrlVect)
        ELSEIF(Ctrl%Method(Ctrl%Current%I(2))==SDMM_R_SCF)THEN
           CALL Invoke('AInv',CtrlVect)
        ELSE
           CALL MondoHalt(MISC_ERROR,' Neither FactoredS or LOrthog invoked in OneEMats ')
        ENDIF 
        CALL Invoke('MakeT' ,      CtrlVect)
!
#ifdef PERIODIC
        IF(CtrlVect(2) == Core)THEN
           CtrlVect=SetCtrlVect(Ctrl,Direct)
           CALL Invoke('P2Use',       CtrlVect,MPIRun_O=.TRUE.)
           CALL Invoke('MakeRho',     CtrlVect)
           CALL Invoke('QCTC',      CtrlVect)
        ELSE
           CALL Invoke('MakeRho',     CtrlVect)
           CALL Invoke('QCTC',     CtrlVect)
        ENDIF
#else
        CALL Invoke('MakeRho',     CtrlVect)
        CALL Invoke('QCTC' ,     CtrlVect)
#endif
!
     END SUBROUTINE OneEMats
!========================================================================================
!
!========================================================================================
   SUBROUTINE LogSCF(Stats,Mssg)
      TYPE(INT_VECT),   INTENT(IN) :: Stats
      CHARACTER(LEN=*), INTENT(IN) :: Mssg
      CALL Logger('------------------------------------------------------------')
      CALL Logger('  Cycl = '//TRIM(IntToChar(Stats%I(1))) &
                //', Base = '//TRIM(IntToChar(Stats%I(2))) &
                //', Geom = '//TRIM(IntToChar(Stats%I(3))) )
      CALL Logger(Mssg)
      CALL Logger('------------------------------------------------------------')
   END SUBROUTINE LogSCF
#ifdef NEEDS_LOTS_OF_WORK
!------------------------------------------------------------------------------------
!       Delete files from previous SCF cycle that are no longer needed 
!       (Needs work to figure out why /bin/rm sometimes does not work when Invoked)
!
        IF(TidyFiles.AND.Ctrl%Current%I(1)>1)THEN
           RmFile=TRIM(Ctrl%Name(Ctrl%Current%I(2)))//'_Cyc'//TRIM(ChCt%CCyc_M1)//'.ExtrapF*'
!           CALL Invoke('rm',(/RmFile/),ExeDir_O='/bin/')
           CALL Logger('/bin/rm '//TRIM(RmFile))
!           CALL SYSTEM('/bin/rm '//TRIM(RmFile))
           RmFile=TRIM(Ctrl%Name(Ctrl%Current%I(2)))//'_Cyc'//TRIM(ChCt%CCyc_M1)//'.J'
!           CALL Invoke('rm',(/RmFile/),ExeDir_O='/bin/')
           CALL Logger('/bin/rm '//TRIM(RmFile))
!           CALL SYSTEM('/bin/rm '//TRIM(RmFile))
           RmFile=TRIM(Ctrl%Name(Ctrl%Current%I(2)))//'_Cyc'//TRIM(ChCt%CCyc_M1)//'.K'
!           CALL Invoke('rm',(/RmFile/),ExeDir_O='/bin/')
           CALL Logger('/bin/rm '//TRIM(RmFile))
!           CALL SYSTEM('/bin/rm '//TRIM(RmFile))
           RmFile=TRIM(Ctrl%Name(Ctrl%Current%I(2)))//'_Cyc'//TRIM(ChCt%CCyc_M1)//'.Rho'
!           CALL Invoke('rm',(/RmFile/),ExeDir_O='/bin/')
           CALL Logger('/bin/rm '//TRIM(RmFile))
           CALL SYSTEM('/bin/rm '//TRIM(RmFile))
!           CALL SYSTEM('/bin/ls -l '//TRIM(MONDO_SCRATCH)//  &
!                       TRIM(SCFName)//'*.* >> '//TRIM(LogFile))
        ELSEIF(TidyFiles.AND.Ctrl%Current%I(1)==0)THEN
           RmFile=TRIM(Ctrl%Name(Ctrl%Current%I(2)))//'_Cyc0.F'
!           CALL Invoke('rm',(/RmFile/),ExeDir_O='/bin/')
           CALL Logger('/bin/rm '//TRIM(RmFile))
           CALL SYSTEM('/bin/rm '//TRIM(RmFile))
           RmFile=TRIM(Ctrl%Name(Ctrl%Current%I(2)))//'_Cyc0.J'
!           CALL Invoke('rm',(/RmFile/),ExeDir_O='/bin/')
           CALL Logger('/bin/rm '//TRIM(RmFile))
           CALL SYSTEM('/bin/rm '//TRIM(RmFile))
           RmFile=TRIM(Ctrl%Name(Ctrl%Current%I(2)))//'_Cyc0.K'
!           CALL Invoke('rm',(/RmFile/),ExeDir_O='/bin/')
           CALL Logger('/bin/rm '//TRIM(RmFile))
           CALL SYSTEM('/bin/rm '//TRIM(RmFile))
           RmFile=TRIM(Ctrl%Name(Ctrl%Current%I(2)))//'_Cyc0.Rho'
!           CALL Invoke('rm',(/RmFile/),ExeDir_O='/bin/')
           CALL Logger('/bin/rm '//TRIM(RmFile))
           CALL SYSTEM('/bin/rm '//TRIM(RmFile))
!           CALL SYSTEM('/bin/ls -l '//TRIM(MONDO_SCRATCH)//  &
!                       TRIM(SCFName)//'*.* >> '//TRIM(LogFile))
        ENDIF
!
!        IF(TidyFiles)THEN
!              RmFile=TRIM(Ctrl%Name(Ctrl%Current%I(2)-1))//'*'
!              CALL Logger('/bin/rm '//TRIM(RmFile))
!              CALL SYSTEM('/bin/rm '//TRIM(RmFile))
!              CALL SYSTEM('/bin/ls -l '//TRIM(MONDO_SCRATCH)//  &
!                          TRIM(SCFName)//'*.* >> '//TRIM(LogFile))
!           ENDIF
#endif
!========================================================================================
!
!========================================================================================
END MODULE




