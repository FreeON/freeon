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
  USE Functionals
  IMPLICIT NONE 
  CONTAINS
!========================================================================================
!
!========================================================================================
     SUBROUTINE SCFCycle(Ctrl)
        TYPE(SCFControls)  :: Ctrl
        INTEGER            :: ICyc,IGeo
!----------------------------------------------------------------------------------------
        ICyc=Ctrl%Current(1)
        IGeo=Ctrl%Current(3)
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
           WRITE(*,*)' Previous state = ',Ctrl%Previous
           WRITE(*,*)' Current  state = ',Ctrl%Current
           CALL MondoHalt(DRIV_ERROR,' bad option in SCFCycle ')
        ENDIF
   END SUBROUTINE
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
        IF(HasDFT(Ctrl%Model(Ctrl%Current(2)))) &
           CALL Invoke('HiCu',     CtrlVect)
        IF(HasHF(Ctrl%Model(Ctrl%Current(2))))&
           CALL Invoke('ONX',     CtrlVect)
        CALL Invoke('FBuild',     CtrlVect,MPIRun_O=.TRUE.)
        IF(Ctrl%DIIS) &
           CALL Invoke('DIIS',     CtrlVect,MPIRun_O=.TRUE.)
        IF(Ctrl%Method(Ctrl%Current(2))==RH_R_SCF)THEN
           CALL Invoke('RHeqs',    CtrlVect)
        ELSEIF(Ctrl%Method(Ctrl%Current(2))==SDMM_R_SCF)THEN
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
        IF(Ctrl%Method(Ctrl%Current(2))==RH_R_SCF)THEN
           CALL Invoke('RHeqs',    CtrlVect)
        ELSEIF(Ctrl%Method(Ctrl%Current(2))==SDMM_R_SCF)THEN
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
        IF(HasDFT(Ctrl%Model(Ctrl%Current(2)))) &
           CALL Invoke('HiCu',     CtrlVect)
        IF(HasHF(Ctrl%Model(Ctrl%Current(2))))  &
           CALL Invoke('ONX',     CtrlVect)
        CALL Invoke('FBuild',     CtrlVect,MPIRun_O=.TRUE.)
        IF(Ctrl%Method(Ctrl%Current(2))==RH_R_SCF)THEN
           CALL Invoke('RHeqs',    CtrlVect)
        ELSEIF(Ctrl%Method(Ctrl%Current(2))==SDMM_R_SCF)THEN
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
        IF(Ctrl%Current(3)==Ctrl%Previous(3))THEN
           CtrlVect=SetCtrlVect(Ctrl,Switch)
           CtrlVect(3)=TRIM(IntToChar(Ctrl%Previous(1)+1))
           CtrlVect(4)=TRIM(IntToChar(Ctrl%Previous(2)))
           CtrlVect(5)=TRIM(IntToChar(Ctrl%Previous(3)))
           CALL Invoke('MakeRho',CtrlVect)
           ShortCtrlVect(1)=TrixFile('Rho',OffSet_O=1,Name_O=Ctrl%Name,Stats_O=Ctrl%Previous)
           ShortCtrlVect(2)=TrixFile('Rho',OffSet_O=0,Name_O=Ctrl%Name,Stats_O=Ctrl%Current )
           CALL Invoke('/bin/cp',ShortCtrlVect,AbsPath_O=.TRUE.)
        ENDIF
!
        CtrlVect=SetCtrlVect(Ctrl,Switch)
        ShortCtrlVect(1)=TrixFile('D',OffSet_O=1,Name_O=Ctrl%Name,Stats_O=Ctrl%Previous)
        ShortCtrlVect(2)=TrixFile('D',OffSet_O=0,Name_O=Ctrl%Name,Stats_O=Ctrl%Current )
        CALL Invoke('/bin/cp',ShortCtrlVect,AbsPath_O=.TRUE.)
        ShortCtrlVect(1)=TrixFile('OrthoD',OffSet_O=1,Name_O=Ctrl%Name,Stats_O=Ctrl%Previous)
        ShortCtrlVect(2)=TrixFile('OrthoD',OffSet_O=0,Name_O=Ctrl%Name,Stats_O=Ctrl%Current )
        CALL Invoke('/bin/cp',ShortCtrlVect,AbsPath_O=.TRUE.)

!       If same basis, different geometry do DM extrapolation 

        IF(Ctrl%Current(3)/=Ctrl%Previous(3))THEN
           CALL Invoke('P2Use',    CtrlVect,MPIRun_O=.TRUE.)
           CALL Invoke('MakeRho',     CtrlVect)
        ENDIF
        CALL Invoke('QCTC' ,     CtrlVect)
        IF(HasDFT(Ctrl%Model(Ctrl%Current(2))))  &
           CALL Invoke('HiCu',     CtrlVect)
        IF(HasHF(Ctrl%Model(Ctrl%Current(2))))   &
           CALL Invoke('ONX',     CtrlVect)
        CALL Invoke('FBuild',     CtrlVect,MPIRun_O=.TRUE.)
        IF(Ctrl%Method(Ctrl%Current(2))==RH_R_SCF)THEN
           CALL Invoke('RHeqs',    CtrlVect)
        ELSEIF(Ctrl%Method(Ctrl%Current(2))==SDMM_R_SCF)THEN
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
        IF(Ctrl%Method(Ctrl%Current(2))==RH_R_SCF)THEN
           CALL Invoke('LowdinO',  CtrlVect)
        ELSEIF(Ctrl%Method(Ctrl%Current(2))==SDMM_R_SCF)THEN
           CALL Invoke('AInv',CtrlVect)
        ELSE
           CALL MondoHalt(MISC_ERROR,' Neither FactoredS or LOrthog invoked in OneEMats ')
        ENDIF 
        CALL Invoke('MakeT' ,      CtrlVect)
!
#ifdef PERIODIC
!        IF(CtrlVect(2) == Core)THEN
!           CtrlVect=SetCtrlVect(Ctrl,Direct)
!           CALL Invoke('P2Use',       CtrlVect,MPIRun_O=.TRUE.)
!           CALL Invoke('MakeRho',     CtrlVect)
!           CALL Invoke('QCTC',      CtrlVect)
!        ELSE
!           CALL Invoke('MakeRho',     CtrlVect)
!           CALL Invoke('QCTC',     CtrlVect)
!        ENDIF
#else
!        CALL Invoke('MakeRho',     CtrlVect)
!        CALL Invoke('QCTC' ,     CtrlVect)
#endif
!
     END SUBROUTINE OneEMats
!========================================================================================
!
!========================================================================================
   SUBROUTINE LogSCF(Stats,Mssg)
      INTEGER,DIMENSION(3),INTENT(IN) :: Stats
      CHARACTER(LEN=*), INTENT(IN) :: Mssg
      CALL Logger('------------------------------------------------------------')
      CALL Logger('  Cycl = '//TRIM(IntToChar(Stats(1))) &
                //', Base = '//TRIM(IntToChar(Stats(2))) &
                //', Geom = '//TRIM(IntToChar(Stats(3))) )
      CALL Logger(Mssg)
      CALL Logger('------------------------------------------------------------')
   END SUBROUTINE LogSCF
#ifdef NEEDS_LOTS_OF_WORK
!------------------------------------------------------------------------------------
!       Delete files from previous SCF cycle that are no longer needed 
!       (Needs work to figure out why /bin/rm sometimes does not work when Invoked)
!
        IF(TidyFiles.AND.Ctrl%Current(1)>1)THEN
           RmFile=TRIM(Ctrl%Name(Ctrl%Current(2)))//'_Cyc'//TRIM(ChCt%CCyc_M1)//'.ExtrapF*'
!           CALL Invoke('rm',(/RmFile/),ExeDir_O='/bin/')
           CALL Logger('/bin/rm '//TRIM(RmFile))
!           CALL SYSTEM('/bin/rm '//TRIM(RmFile))
           RmFile=TRIM(Ctrl%Name(Ctrl%Current(2)))//'_Cyc'//TRIM(ChCt%CCyc_M1)//'.J'
!           CALL Invoke('rm',(/RmFile/),ExeDir_O='/bin/')
           CALL Logger('/bin/rm '//TRIM(RmFile))
!           CALL SYSTEM('/bin/rm '//TRIM(RmFile))
           RmFile=TRIM(Ctrl%Name(Ctrl%Current(2)))//'_Cyc'//TRIM(ChCt%CCyc_M1)//'.K'
!           CALL Invoke('rm',(/RmFile/),ExeDir_O='/bin/')
           CALL Logger('/bin/rm '//TRIM(RmFile))
!           CALL SYSTEM('/bin/rm '//TRIM(RmFile))
           RmFile=TRIM(Ctrl%Name(Ctrl%Current(2)))//'_Cyc'//TRIM(ChCt%CCyc_M1)//'.Rho'
!           CALL Invoke('rm',(/RmFile/),ExeDir_O='/bin/')
           CALL Logger('/bin/rm '//TRIM(RmFile))
           CALL SYSTEM('/bin/rm '//TRIM(RmFile))
!           CALL SYSTEM('/bin/ls -l '//TRIM(MONDO_SCRATCH)//  &
!                       TRIM(SCFName)//'*.* >> '//TRIM(LogFile))
        ELSEIF(TidyFiles.AND.Ctrl%Current(1)==0)THEN
           RmFile=TRIM(Ctrl%Name(Ctrl%Current(2)))//'_Cyc0.F'
!           CALL Invoke('rm',(/RmFile/),ExeDir_O='/bin/')
           CALL Logger('/bin/rm '//TRIM(RmFile))
           CALL SYSTEM('/bin/rm '//TRIM(RmFile))
           RmFile=TRIM(Ctrl%Name(Ctrl%Current(2)))//'_Cyc0.J'
!           CALL Invoke('rm',(/RmFile/),ExeDir_O='/bin/')
           CALL Logger('/bin/rm '//TRIM(RmFile))
           CALL SYSTEM('/bin/rm '//TRIM(RmFile))
           RmFile=TRIM(Ctrl%Name(Ctrl%Current(2)))//'_Cyc0.K'
!           CALL Invoke('rm',(/RmFile/),ExeDir_O='/bin/')
           CALL Logger('/bin/rm '//TRIM(RmFile))
           CALL SYSTEM('/bin/rm '//TRIM(RmFile))
           RmFile=TRIM(Ctrl%Name(Ctrl%Current(2)))//'_Cyc0.Rho'
!           CALL Invoke('rm',(/RmFile/),ExeDir_O='/bin/')
           CALL Logger('/bin/rm '//TRIM(RmFile))
           CALL SYSTEM('/bin/rm '//TRIM(RmFile))
!           CALL SYSTEM('/bin/ls -l '//TRIM(MONDO_SCRATCH)//  &
!                       TRIM(SCFName)//'*.* >> '//TRIM(LogFile))
        ENDIF
!
!        IF(TidyFiles)THEN
!              RmFile=TRIM(Ctrl%Name(Ctrl%Current(2)-1))//'*'
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




