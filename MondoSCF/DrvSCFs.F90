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
    SUBROUTINE SCFCycle(Ctrl,Begin)
      TYPE(SCFControls)     :: Ctrl
      INTEGER               :: ICyc,IBas,IGeo
      INTEGER, DIMENSION(3) :: Begin
!----------------------------------------------------------------------------------------
      ICyc=Ctrl%Current(1)
      IBas=Ctrl%Current(2)
      IGeo=Ctrl%Current(3)
!     Hack ...
      Ctrl%DIIS=.TRUE.
      Ctrl%InkFok=.FALSE.
!     Do an SCF cycle 
      IF(Icyc==Begin(1)) THEN
         IF(IBas==Begin(2)) THEN
            IF(IGeo==Begin(3)) THEN
               IF(Ctrl%SuperP .OR. Ctrl%Rest) THEN
!                 Restricted direct SCF from Super position of Atomic orbitals
                  CALL SuperRSCF(Ctrl)
               ELSE
!                 Restricted direct SCF from core Hamiltonian
                  CALL CoreRSCF(Ctrl)
               ENDIF
            ELSE
!              Restricted SCF with geometry switch 
               CALL SwitchGeometryRSCF(Ctrl)
            ENDIF
         ELSE
!            Restricted SCF with basis-set switch 
            CALL SwitchBasisRSCF(Ctrl)
         ENDIF
      ELSE
!        Restricted direct SCF cycle
         CALL DirectRSCF(Ctrl)
      ENDIF
    END SUBROUTINE SCFCycle
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
!
      CALL CleanScratch(Ctrl,'Cycl')              
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
    SUBROUTINE SwitchBasisRSCF(Ctrl)
      TYPE(SCFControls)               :: Ctrl
      CHARACTER(LEN=DCL),DIMENSION(2) :: ShortCtrlVect        
      INTEGER,DIMENSION(3)            :: NextStats
!---------------------------------------------------------------------------------------
      CALL LogSCF(Ctrl%Current,'Restricted SCF start via basis set switch')
!     Build new rho from last density matrix and rename it
      CtrlVect=SetCtrlVect(Ctrl,Switch)
      CtrlVect(3)=TRIM(IntToChar(Ctrl%Previous(1)+1))
      CtrlVect(4)=TRIM(IntToChar(Ctrl%Previous(2)))
      CtrlVect(5)=TRIM(IntToChar(Ctrl%Previous(3)))
      CALL Invoke('MakeRho',CtrlVect)
      ShortCtrlVect(1)=TrixFile('Rho',OffSet_O=1,Name_O=Ctrl%Name,Stats_O=Ctrl%Previous)
      ShortCtrlVect(2)=TrixFile('Rho',OffSet_O=0,Name_O=Ctrl%Name,Stats_O=Ctrl%Current )
      CALL Invoke('/bin/cp',ShortCtrlVect,AbsPath_O=.TRUE.)
!     Copy the old density matrix to the new density matrix   
      CtrlVect=SetCtrlVect(Ctrl,Switch)
      ShortCtrlVect(1)=TrixFile('D',OffSet_O=1,Name_O=Ctrl%Name,Stats_O=Ctrl%Previous)
      ShortCtrlVect(2)=TrixFile('D',OffSet_O=0,Name_O=Ctrl%Name,Stats_O=Ctrl%Current )
      CALL Invoke('/bin/cp',ShortCtrlVect,AbsPath_O=.TRUE.)
      ShortCtrlVect(1)=TrixFile('OrthoD',OffSet_O=1,Name_O=Ctrl%Name,Stats_O=Ctrl%Previous)
      ShortCtrlVect(2)=TrixFile('OrthoD',OffSet_O=0,Name_O=Ctrl%Name,Stats_O=Ctrl%Current )
      CALL Invoke('/bin/cp',ShortCtrlVect,AbsPath_O=.TRUE.)
!
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
!
      CALL CleanScratch(Ctrl,'Base')  
    END SUBROUTINE SwitchBasisRSCF
!========================================================================================
!    Perform a resctricted SCF cycle with a geometry switch
!    For Now-Start from begining
!========================================================================================
    SUBROUTINE SwitchGeometryRSCF(Ctrl)
      TYPE(SCFControls)               :: Ctrl
      CHARACTER(LEN=DCL),DIMENSION(2) :: ShortCtrlVect        
      INTEGER,DIMENSION(3)            :: NextStats
!---------------------------------------------------------------------------------------
      CALL LogSCF(Ctrl%Current,'Restricted SCF start via geometry switch')
      CtrlVect=SetCtrlVect(Ctrl,Direct)
      CALL Invoke('P2Use',    CtrlVect,MPIRun_O=.TRUE.)
      CALL Invoke('MakeRho',     CtrlVect)
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
!
      CALL CleanScratch(Ctrl,'Geom')                              
    END SUBROUTINE SwitchGeometryRSCF
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
    END SUBROUTINE OneEMats
!========================================================================================
!
!========================================================================================
    SUBROUTINE LogSCF(Stats,Mssg)
      INTEGER,DIMENSION(3),INTENT(IN) :: Stats
      CHARACTER(LEN=*), INTENT(IN)    :: Mssg
      CALL Logger('------------------------------------------------------------')
      CALL Logger('  Cycl = '//TRIM(IntToChar(Stats(1))) &
                //', Base = '//TRIM(IntToChar(Stats(2))) &
                //', Geom = '//TRIM(IntToChar(Stats(3))) )
      CALL Logger(Mssg)
      CALL Logger('------------------------------------------------------------')
    END SUBROUTINE LogSCF
!========================================================================================
!   Clean the Scratch Directory
!========================================================================================
    SUBROUTINE CleanScratch(Ctrl,Action)
      TYPE(SCFControls)      :: Ctrl
      CHARACTER(LEN=4)       :: Action
      CHARACTER(LEN=120)     :: RemoveFile
      INTEGER                :: IGeo,IBas,ICyc
!      
      ICyc=Ctrl%Current(1)
      IBas=Ctrl%Current(2)
      IGeo=Ctrl%Current(3)
      IF(Action=='Cycl') THEN
         RemoveFile=TRIM(Ctrl%Name) //'_Geom#'//TRIM(INTTOCHAR(IGeo))  & 
                                    //'_Base#'//TRIM(INTTOCHAR(IBas))  &
                                    //'_Cycl#'//TRIM(INTTOCHAR(ICyc-1)) //'.Rho'                                
         CALL SYSTEM('/bin/rm '//RemoveFile)
         RemoveFile=TRIM(Ctrl%Name) //'_Geom#'//TRIM(INTTOCHAR(IGeo))  & 
                                    //'_Base#'//TRIM(INTTOCHAR(IBas))  &
                                    //'_Cycl#'//TRIM(INTTOCHAR(ICyc-1)) //'.J'                                
         CALL SYSTEM('/bin/rm '//RemoveFile)
         RemoveFile=TRIM(Ctrl%Name) //'_Geom#'//TRIM(INTTOCHAR(IGeo))  & 
                                    //'_Base#'//TRIM(INTTOCHAR(IBas))  &
                                    //'_Cycl#'//TRIM(INTTOCHAR(ICyc-1)) //'.K*'                                
         CALL SYSTEM('/bin/rm '//RemoveFile)
         RemoveFile=TRIM(Ctrl%Name) //'_Geom#'//TRIM(INTTOCHAR(IGeo))  & 
                                    //'_Base#'//TRIM(INTTOCHAR(IBas))  &
                                    //'_Cycl#'//TRIM(INTTOCHAR(ICyc-1)) //'.D'                                
         CALL SYSTEM('/bin/rm '//RemoveFile)
         RemoveFile=TRIM(Ctrl%Name) //'_Geom#'//TRIM(INTTOCHAR(IGeo))  & 
                                    //'_Base#'//TRIM(INTTOCHAR(IBas))  &
                                    //'_Cycl#'//TRIM(INTTOCHAR(ICyc-1)) //'.F'                                
         CALL SYSTEM('/bin/rm '//RemoveFile)
      ELSEIF(Action=='Base') THEN
         RemoveFile=TRIM(Ctrl%Name) //'_Geom#'//TRIM(INTTOCHAR(IGeo))    & 
                                    //'_Base#'//TRIM(INTTOCHAR(IBas-1)) //'*' 
         CALL SYSTEM('/bin/rm '//RemoveFile)  
      ELSEIF(Action=='Geom') THEN
         RemoveFile=TRIM(Ctrl%Name) //'_Geom#'//TRIM(INTTOCHAR(IGeo-1))//'*' 
         CALL SYSTEM('/bin/rm '//RemoveFile) 
      ENDIF
    END SUBROUTINE CleanScratch
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




