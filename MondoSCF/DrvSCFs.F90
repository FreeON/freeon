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
#ifdef MMech
  USE Dynamo
  USE Mechanics
#endif
  USE LinAlg
  IMPLICIT NONE 
  CONTAINS
!========================================================================================
!
!========================================================================================
     SUBROUTINE OneSCF(Ctrl,Sum_O)
        TYPE(SCFControls)  :: Ctrl
        INTEGER            :: J
        INTEGER            :: ISCF,ICyc
        LOGICAL,OPTIONAL   :: Sum_O
        LOGICAL            :: Summry
!
        IF(PRESENT(Sum_O))THEN
           Summry=Sum_O
        ELSE
           Summry=.TRUE.
        ENDIF
!
!       Do QM Calc on the SCF/DFT level
!
        DO ICyc=0,MaxSCFs
           Ctrl%Current(1)=ICyc
           CALL SetGlobalCtrlIndecies(Ctrl)
           IF(ICyc==0)CALL OneEMats(Ctrl)
!          Do an SCF cycle
           CALL SCFCycle(Ctrl)          
           IF(ConvergedQ(Ctrl))THEN
              IF(Summry) CALL SCFSummry(Ctrl)
              CALL VisDX(Ctrl)
              CALL CleanScratch(Ctrl,'CleanOnConverge')
              Ctrl%Previous=Ctrl%Current
              RETURN
           ENDIF
           Ctrl%Previous=Ctrl%Current
        ENDDO
        CALL MondoHalt(DRIV_ERROR,'Failed to converge SCF in ' &
                                  //TRIM(IntToChar(MaxSCFs))   &
                                  //' SCF iterations.')
    END SUBROUTINE OneSCF
!========================================================================================
!
!========================================================================================
    SUBROUTINE SCFCycle(Ctrl)
      TYPE(SCFControls)     :: Ctrl
!------------------------------------------------------
      CALL DensityBuild(Ctrl)
      CALL FockBuild(Ctrl)
      CALL SolveSCF(Ctrl)
      CALL CleanScratch(Ctrl)
    END SUBROUTINE SCFCycle
!==========================================================================
!   Build a density by hook or by crook
!==========================================================================
    SUBROUTINE DensityBuild(Ctrl)
      TYPE(SCFControls)    :: Ctrl
      LOGICAL,PARAMETER    :: DensityProject=.FALSE.
!---------------------------------------------------------------------------------------
!
      IF(CCyc==0.AND.CBas==PBas.AND.CGeo/=1)THEN
!	  .AND.Ctrl%Extrap>EXTRAP_GEOM_RSTRT)THEN
         ! Density projection is now the default guess for new geometries
!         IF(Ctrl%Extrap==EXTRAP_GEOM_PRJCT)THEN
!           Projection of density matrix between geometries
            CALL LogSCF(Ctrl%Current,'Geometry projection from configuration #' &
                                   //TRIM(PrvGeom)//' to configuration# ' &
                                   //TRIM(CurGeom)//'.')
!           Create density from last SCF 
            CtrlVect=SetCtrlVect(Ctrl,'Project')
!         ELSEIF(Ctrl%Extrap==EXTRAP_GEOM_INTRP)THEN
!           Extrapolation of density matrix between geometries
!            CALL LogSCF(Ctrl%Current,'Geometry extrapolation from configuration #' &
!                                   //TRIM(PrvGeom)//' to configuration# ' &
!                                   //TRIM(CurGeom)//'.')
!           Create density from last SCFs DM (ICyc+1)
!            Ctrl%Previous(1)=Ctrl%Previous(1)+1
!            CtrlVect=SetCtrlVect(Ctrl,'Extrapolate')
!         ENDIF
         CALL Invoke('P2Use',CtrlVect,MPIRun_O=.TRUE.)
         CALL Invoke('MakeRho',CtrlVect)
         CALL CleanScratch(Ctrl,'CleanLastGeom')
      ELSEIF(Ctrl%Rest)THEN
!        Restart from a previous density
         CALL LogSCF(Ctrl%Current,'Restarting from '//TRIM(Ctrl%OldInfo),.TRUE.)
         ! Make sure previous and current are set correctly.
         Ctrl%Previous=Ctrl%Current
         CtrlVect=SetCtrlVect(Ctrl,'Restart')
         CALL Invoke('P2Use',CtrlVect,MPIRun_O=.TRUE.)
         CALL Invoke('MakeRho',CtrlVect)
#ifdef PARALLEL
         CALL Invoke('ParaMakeRho',CtrlVect,MPIRun_O=.TRUE.)
#endif
         Ctrl%Rest=.FALSE.
      ELSEIF(CCyc==0.AND.CBas/=PBas)THEN
!        Basis set switch
         CALL LogSCF(Ctrl%Current,'Switching basis sets from ' &
                                //TRIM(Ctrl%BName(PBas))//' to ' &
                                //TRIM(Ctrl%BName(CBas))//'.')
!        Create density from last SCFs DM (ICyc+1)
         Ctrl%Previous(1)=Ctrl%Previous(1)+1
         CtrlVect=SetCtrlVect(Ctrl,'BasisSetSwitch')
         CALL Invoke('MakeRho',CtrlVect)
     ELSEIF(CCyc==0.AND.CBas==1.AND.CGeo==1)THEN
!        Start up density from atomic superposition in diagonal minimial basis 
         CALL LogSCF(Ctrl%Current,' Density guess from superposition of AOs.')
!        Create density from guess
         CtrlVect=SetCtrlVect(Ctrl,'DensitySuperposition')
         CALL Invoke('P2Use',CtrlVect,MPIRun_O=.TRUE.)
         CALL Invoke('MakeRho',CtrlVect)
#ifdef PARALLEL
         write(*,*) 'calling ParaMakeRho...'
         CALL Invoke('ParaMakeRho',CtrlVect,MPIRun_O=.TRUE.)
#endif
      ELSEIF(Ctrl%InkFok)THEN
!        First do a full density build 
         CALL LogSCF(Ctrl%Current,'Full density build followed by delta density build.',.TRUE.)
         CtrlVect=SetCtrlVect(Ctrl,'Direct')
         CALL Invoke('MakeRho',CtrlVect)
         CtrlVect=SetCtrlVect(Ctrl,'InkFok')
         CALL Invoke('MakeRho',CtrlVect)
      ELSE
!        Regular old density build
         CALL LogSCF(Ctrl%Current,'Standard density build',.TRUE.)
         CtrlVect=SetCtrlVect(Ctrl,'Direct')
         CALL Invoke('MakeRho',CtrlVect)
#ifdef PARALLEL
         write(*,*) 'calling ParaMakeRho...'
         CALL Invoke('ParaMakeRho',CtrlVect,MPIRun_O=.TRUE.)
#endif
      ENDIF
    END SUBROUTINE DensityBuild
!==========================================================================
!   Build a Fock matrix with extrapolation of the 
!   orthongoal Fock matrix using Pulay DIIS.
!==========================================================================
    SUBROUTINE FockBuild(Ctrl)
      TYPE(SCFControls)  :: Ctrl
      INTEGER            :: Modl
      LOGICAL            :: DoDIIs
!-----------------------------------------------------------
      IF(Ctrl%InkFok)THEN
         CALL LogSCF(Current,'Building an incremental Fock matrix')
      ELSE
         CALL LogSCF(Current,'Building the Fock matrix')
      ENDIF
      DoDIIS=CCyc>0
      Modl=Ctrl%Model(CBas)
      CALL Invoke('QCTC',CtrlVect)
      IF(HasDFT(Modl))CALL Invoke('HiCu',CtrlVect,MPIRun_O=.TRUE.)
      IF(HasHF(Modl)) CALL Invoke('ONX',CtrlVect)
      CALL Invoke('FBuild',CtrlVect,MPIRun_O=.TRUE.)
      IF(DoDIIS) &
      CALL Invoke('DIIS',CtrlVect,MPIRun_O=.TRUE.)
      IF(CtrlVect(2)=='BasisSetSwitch') &  
         CALL CleanScratch(Ctrl,'CleanLastBase')
    END SUBROUTINE FockBuild
!==========================================================================
!   Solve the SCF equations 
!==========================================================================
    SUBROUTINE SolveSCF(Ctrl)
      TYPE(SCFControls)  :: Ctrl
!-----------------------------------------------------------
      IF(Ctrl%Method(CBas)==RH_R_SCF)THEN
         CALL LogSCF(Current,'Solving SCF equations with Roothann-Hall')
         CALL Invoke('RHeqs',CtrlVect)
      ELSEIF(Ctrl%Method(CBas)==SDMM_R_SCF) THEN
         CALL LogSCF(Current,'Solving SCF equations with SDMM')
         CALL Invoke('SDMM',CtrlVect,MPIRun_O=.TRUE.) 
      ELSEIF(Ctrl%Method(CBas)==PM_R_SCF) THEN
         CALL LogSCF(Current,'Solving SCF equations with PM')
         CALL Invoke('PM',CtrlVect,MPIRun_O=.TRUE.) 
      ELSEIF(Ctrl%Method(CBas)==SP2_R_SCF) THEN
         CALL LogSCF(Current,'Solving SCF equations with SP2')
         CALL Invoke('SP2',CtrlVect,MPIRun_O=.TRUE.) 
      ELSEIF(Ctrl%Method(CBas)==SP4_R_SCF) THEN
         CALL LogSCF(Current,'Solving SCF equations with SP4')
         CALL Invoke('SP4',CtrlVect,MPIRun_O=.TRUE.) 
      ELSEIF(Ctrl%Method(CBas)==TS4_R_SCF) THEN
         CALL LogSCF(Current,'Solving SCF equations with TS4')
         CALL Invoke('TS4',CtrlVect,MPIRun_O=.TRUE.) 
      ELSE
         CALL MondoHalt(99,'No Method Chosen')
      ENDIF
!
      CALL Invoke('SCFstats',CtrlVect,MPIRun_O=.TRUE.)                                     
    END SUBROUTINE SolveSCF
!========================================================================================
!
!========================================================================================
    SUBROUTINE OneEMats(Ctrl)
      TYPE(SCFControls)             :: Ctrl
      INTEGER                       :: ISet,I
!
#ifdef PERIODIC
      CALL LogSCF(Current,'Peridic Far-Field Tensor',.TRUE.)
      CtrlVect=SetCtrlVect(Ctrl,'MakingPFFT')
      CALL Invoke('MakePFFT',CtrlVect,MPIRun_O=.TRUE.)
#endif
      CALL LogSCF(Current,'One-electron matrices.',.TRUE.)
      CtrlVect=SetCtrlVect(Ctrl,'OneElectron')
      CALL Invoke('MakeS',CtrlVect,MPIRun_O = .TRUE.)
      IF(Ctrl%Method(Ctrl%Current(2))==RH_R_SCF)THEN
         CALL Invoke('LowdinO',  CtrlVect)
      ELSE
         CALL Invoke('AInv',CtrlVect)
      ENDIF
      CALL Invoke('MakeT',CtrlVect,MPIRun_O=.TRUE.)
!
    END SUBROUTINE OneEMats
!========================================================================================
!
!========================================================================================
    SUBROUTINE VisDX(Ctrl)
      TYPE(SCFControls)             :: Ctrl
      IF(Ctrl%Vis/=VIS_DX_RHOPOT)RETURN
      Ctrl%Current(1)=Ctrl%Current(1)+1
      CtrlVect=SetCtrlVect(Ctrl,'Visualization')
      Ctrl%Current(1)=Ctrl%Current(1)-1
      CALL Invoke('MakeRho',CtrlVect)
      CALL Invoke('VisDX',CtrlVect)
    END SUBROUTINE VisDX
!========================================================================================
!
!========================================================================================
    SUBROUTINE LogSCF(Stats,Mssg,Banner_O)
      INTEGER,DIMENSION(3),INTENT(IN) :: Stats
      CHARACTER(LEN=*), INTENT(IN)    :: Mssg
      CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Mssg2
      LOGICAL,OPTIONAL                :: Banner_O
      LOGICAL                         :: Banner
      IF(PRESENT(Banner_O))THEN
         Banner=Banner_O
      ELSE
         Banner=.FALSE.
      ENDIF
      IF(Banner)THEN
         CALL Logger('!------------------------------------------------------------')
         CALL Logger('  Cycl = '//TRIM(IntToChar(Stats(1))) &
                   //', Base = '//TRIM(IntToChar(Stats(2))) &
                   //', Geom = '//TRIM(IntToChar(Stats(3))) )
      ENDIF
      Mssg2='  '//TRIM(Mssg)
      CALL Logger(Mssg2)
    END SUBROUTINE LogSCF
!========================================================================================
!   Clean the Scratch Directory
!========================================================================================
    SUBROUTINE CleanScratch(Ctrl,Action_O)
      TYPE(SCFControls)              :: Ctrl
      CHARACTER(LEN=*), OPTIONAL     :: Action_O
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: Action,RemCur,RemPrv,RemoveFile
      INTEGER                        :: Modl
!---------------------------------------------------------------------------------------      
      IF(PRESENT(Action_O))THEN
         Action=Action_O
      ELSE
         Action='CleanLastCycl'
      ENDIF
      IF(Action=='CleanLastCycl') THEN
         Modl=Ctrl%Model(CBas)
         RemCur=TRIM(ScrName)//'_Geom#' &
              //TRIM(CurGeom)//'_Base#' &
              //TRIM(CurBase)//'_Cycl#' &
              //TRIM(CurCycl)
         IF(PrvGeom/=CurGeom) THEN
           RemPrv=TRIM(ScrName)//'_Geom#' &
                //TRIM(PrvGeom)//'_Base#' &
                //TRIM(PrvBase)//'_Cycl#' &
                //TRIM(PrvCycl)
         ELSE IF(PrvBase/=CurBase) THEN
           RemPrv=TRIM(ScrName)//'_Geom#' &
                //TRIM(CurGeom)//'_Base#' &
                //TRIM(PrvBase)//'_Cycl#' &
                //TRIM(PrvCycl)
         ELSE
           RemPrv=TRIM(ScrName)//'_Geom#' &
                //TRIM(CurGeom)//'_Base#' &
                //TRIM(CurBase)//'_Cycl#' &
                //TRIM(PrvCycl)
         ENDIF
         RemoveFile=TRIM(RemPrv)//'.Rho'                                
!        RemoveFile=TRIM(RemCur)//'.Rho'                                
         CALL SYSTEM('/bin/rm -f  '//RemoveFile)
         RemoveFile=TRIM(RemPrv)//'.J'                                
         CALL SYSTEM('/bin/rm -f  '//RemoveFile)
         IF(HasDFT(Modl))THEN
            RemoveFile=TRIM(RemCur)//'.Kxc'                                
            CALL SYSTEM('/bin/rm -f  '//RemoveFile)
         ENDIF
         IF(HasHF(Modl))THEN
            RemoveFile=TRIM(RemPrv)//'.K'                                
            CALL SYSTEM('/bin/rm -f  '//RemoveFile)
         ENDIF
         IF(CCyc>0)THEN
            RemoveFile=TRIM(RemPrv)//'.OrthoD'                                
            CALL SYSTEM('/bin/rm -f  '//RemoveFile)
!            RemoveFile=TRIM(RemCur)//'.F_DIIS'                                
!            CALL SYSTEM('/bin/rm -f  '//RemoveFile)
            RemoveFile=TRIM(RemPrv)//'.F'                                
            CALL SYSTEM('/bin/rm -f  '//RemoveFile)
         ENDIF
      ELSEIF(Action=='CleanOnConverge')THEN
         RemoveFile=TRIM(ScrName)//'*.E' 
         CALL SYSTEM('/bin/rm -f  '//RemoveFile)  
         RemoveFile=TRIM(ScrName)//'*.OrthoF' 
         CALL SYSTEM('/bin/rm -f  '//RemoveFile)  
         RemoveFile=TRIM(RemCur)//'.D'                                
         CALL SYSTEM('/bin/rm -f  '//RemoveFile)
         IF(Ctrl%Method(Ctrl%Current(2))==SDMM_R_SCF)THEN
            RemoveFile=TRIM(ScrName)//'*.Z' 
            CALL SYSTEM('/bin/rm -f  '//RemoveFile)  
            RemoveFile=TRIM(ScrName)//'*.ZT' 
            CALL SYSTEM('/bin/rm -f  '//RemoveFile)  
         ELSE
            RemoveFile=TRIM(ScrName)//'*.X' 
            CALL SYSTEM('/bin/rm -f  '//RemoveFile)  
         ENDIF
         RemoveFile=TRIM(ScrName)//'*.T' 
         CALL SYSTEM('/bin/rm -f  '//RemoveFile)  
      ELSEIF(Action=='CleanLastGeom')THEN
         RemoveFile=TRIM(ScrName) //'_Geom#'//TRIM(PrvGeom)//'*.S' 
         CALL SYSTEM('/bin/rm -f  '//RemoveFile) 
         RemoveFile=TRIM(ScrName) //'_Geom#'//TRIM(PrvGeom)//'*.F' 
         CALL SYSTEM('/bin/rm -f  '//RemoveFile) 
         RemoveFile=TRIM(ScrName) //'_Geom#'//TRIM(PrvGeom)//'*.D' 
!         CALL SYSTEM('/bin/rm -f  '//RemoveFile) 
!         RemoveFile=TRIM(ScrName) //'_Geom#'//TRIM(PrvGeom)//'*.OrthoD' 
         CALL SYSTEM('/bin/rm -f  '//RemoveFile) 
!        RemoveFile=TRIM(ScrName) //'_Geom#'//TRIM(PrvGeom)//'*.Rho' 
         CALL SYSTEM('/bin/rm -f  '//RemoveFile) 
      ELSEIF(Action=='CleanLastBase')THEN
         RemoveFile=TRIM(ScrName) //'*_Base#'//TRIM(PrvBase)//'*.S' 
         CALL SYSTEM('/bin/rm -f  '//RemoveFile) 
         RemoveFile=TRIM(ScrName) //'*_Base#'//TRIM(PrvBase)//'*.F' 
         CALL SYSTEM('/bin/rm -f  '//RemoveFile) 
         RemoveFile=TRIM(ScrName) //'*_Base#'//TRIM(PrvBase)//'*.D' 
         CALL SYSTEM('/bin/rm -f  '//RemoveFile) 
         RemoveFile=TRIM(ScrName) //'*_Base#'//TRIM(PrvBase)//'*.OrthoD' 
         CALL SYSTEM('/bin/rm -f  '//RemoveFile) 
!        RemoveFile=TRIM(ScrName) //'*_Base#'//TRIM(PrvGeom)//'*.Rho' 
         CALL SYSTEM('/bin/rm -f  '//RemoveFile) 
      ENDIF
    END SUBROUTINE CleanScratch
!========================================================================================
!
!========================================================================================
      FUNCTION ConvergedQ(Ctrl)
         TYPE(SCFControls) :: Ctrl
         LOGICAL           :: ConvergedQ
         REAL(DOUBLE)      :: ETotA,ETotB,dETot,ETotQ, &
                              DMaxA,DMaxB,dDMax,DMaxQ, &
                              DIISA,DIISB,dDIIS,DIISQ, &
                              ETest,DTest
         TYPE(INT_VECT)    :: Stat
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
!-----------------------------------------------------------------------
         ConvergedQ=.FALSE.
         IF(CCyc==0)RETURN
!-----------------------------------------------------------------------
!        Mark the current status
         CALL New(Stat,3)
         Stat%I=Ctrl%Current
         CALL Put(Stat,'current')
         CALL Delete(Stat)
!        Gather convergence parameters
         CALL Get(EtotA,'Etot',StatsToChar(Ctrl%Previous))
         CALL Get(EtotB,'Etot',StatsToChar(Ctrl%Current))
         CALL Get(DMaxA,'DMax',StatsToChar(Ctrl%Previous))
         CALL Get(DMaxB,'DMax',StatsToChar(Ctrl%Current))
         IF(CCyc==0)THEN
            DIISA=1.D0
            DIISB=1.D0
         ELSEIF(CCyc==1)THEN
            CALL Get(DIISB,'diiserr',StatsToChar(Ctrl%Current))
            DIISB=DIISA
         ELSE
            CALL Get(DIISB,'diiserr',StatsToChar(Ctrl%Current))
            CALL Get(DIISA,'diiserr',StatsToChar(Ctrl%Previous))
         ENDIF
!        Close the InfFile
!        Absolute numbers
         dETot=ABS(ETotA-ETotB)
         dDMax=ABS(DMaxA-DMaxB)
         dDIIS=ABS(DIISA-DIISB)
!        Relative numbers (Quotients)
         ETotQ=dETot/ABS(ETotB)
         DMaxQ=dDMax/ABS(DMaxB+1.D-50)
         DIISQ=dDIIS/ABS(DIISB+1.D-50)
!----------------------------------------------------------------------------
!        Check for convergence
         ETest=ETol(Ctrl%AccL(CBas))
         DTest=DTol(Ctrl%AccL(CBas))
!        Check for absolute convergence below thresholds
!        and approach from correct direction.
         CALL OpenASCII(OutFile,Out)
         IF(((DMaxB<dTest.AND.ETotQ<ETest).OR.DMaxB<5D-1*dTest).AND.ETotB<ETotA)THEN
            Mssg='Normal SCF convergence.'
            ConvergedQ=.TRUE.
         ENDIF
!        Accept convergence from wrong side if thresholds are tightend.
         IF(DMaxB<dTest*75D-2.AND.ETotQ<ETest*3D-1)THEN
            Mssg='Normal SCF convergence.'
            ConvergedQ=.TRUE.
         ENDIF
         !WRITE(Out,*)'ETest = ',ETest
         !WRITE(Out,*)'DTest = ',DTest
         !WRITE(Out,*)'ETotQ = ',ETotQ
         !WRITE(Out,*)'ETotA = ',ETotA
         !WRITE(Out,*)'ETotB = ',ETotB
         !WRITE(Out,*)'DIISQ = ',DIISQ
         !WRITE(Out,*)'DMaxQ = ',DMaxQ
         !WRITE(Out,*)'DIISA = ',DIISA
         !WRITE(Out,*)'DIISB = ',DIISB
         !WRITE(Out,*)'DMaxA = ',DMaxA
         !WRITE(Out,*)'DMaxB = ',DMaxB
         CLOSE(Out)
!        Look for non-decreasing errors due to incomplete numerics
         IF(DIISQ<1.D-1.AND.DMaxQ<1.D-1.AND.CCyc>2)THEN
            IF(DIISB>DIISA.AND.DMaxB>DMaxA)THEN
               Mssg='SCF hit DIIS&DMax increase.'
               ConvergedQ=.TRUE.
            ENDIF
         ELSEIF(DIISQ<1.D-2.AND.DMaxQ<1.D-2.AND.CCyc>2)THEN
            IF(DIISB>DIISA)THEN
               Mssg='SCF hit DIIS increase.'
               ConvergedQ=.TRUE.
            ELSEIF(DMaxQ<1D-1.AND.DMaxB>DMaxA)THEN
               Mssg='SCF hit DMAX increase.'
               ConvergedQ=.TRUE.
            ENDIF
         ELSEIF((DIISQ<1D-3.OR.DMaxQ<1D-3).AND.CCyc>2)THEN
            Mssg='SCF convergence due to DIIS stagnation.'
            ConvergedQ=.TRUE.
         ENDIF
!        Logic for incremental Fock builds
         IF(.NOT.ConvergedQ.AND.DMaxB<1D-2.AND. &
            Ctrl%ShudInk.AND..NOT.Ctrl%BeenInkn)THEN
!           Turn on incremental Fock builds if we should, but havnt yet
            Ctrl%InkFok=.TRUE.
            Ctrl%BeenInkn=.TRUE.
         ELSEIF(ETotB>ETotA.AND.Ctrl%InkFok)THEN
!           If approaching convergence from the wrong direction, 
!           turn off IncFok.
            Ctrl%InkFok=.FALSE.
         ELSEIF(ConvergedQ.AND.Ctrl%InkFok)THEN
!           Turn full builds back on if convergence reached
!           with incremental methods.  
            ConvergedQ=.FALSE.
            Ctrl%InkFok=.FALSE.
            Ctrl%BeenInkn=.TRUE.
         ELSEIF(ConvergedQ.AND.(.NOT.Ctrl%InkFok))THEN
!           If truely converged after full Fock builds, then 
!           we are done with benn inkn, so reset the flag.
            Ctrl%BeenInkn=.FALSE.
         ENDIF            
!        Convergence announcement
         IF(ConvergedQ.AND.PrintFlags%Key>DEBUG_NONE)THEN
            WRITE(*,*)TRIM(Mssg)             
            CALL OpenASCII(OutFile,Out)
            WRITE(Out,*)TRIM(Mssg)             
            CLOSE(Out)
         ENDIF
!        Load statistics 
         Ctrl%NCyc(CBas)=CCyc
         Ctrl%Stats(CCyc,:)=(/ETotB,ETotQ,dDMax,DIISB,dETot/)
      END FUNCTION ConvergedQ
!
      SUBROUTINE SCFSummry(Ctrl)
         TYPE(SCFControls) :: Ctrl
         TYPE(CRDS)        :: GM 
         INTEGER           :: I,K,IBas,ICyc,IGeo,NCyc,Mthd,BTyp
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
!------------------------------------------------------------------
         IF(PrintFlags%Key/=DEBUG_MAXIMUM)RETURN

!        Simplify notation
         ICyc=Ctrl%Current(1)
         IBas=Ctrl%Current(2)
         IGeo=Ctrl%Current(3)
         Mthd=Ctrl%Method(IBas)
         NCyc=Ctrl%NCyc(IBas)
         IF(PrintFlags%Key<=DEBUG_MINIMUM)THEN
            CALL OpenASCII(OutFile,Out)         
            CALL PrintProtectL(Out)
            IF(IGeo==1.AND.IBas==1)THEN
               WRITE(Out,11)
               WRITE(Out,21)
               WRITE(Out,22)
               WRITE(Out,11)
            ENDIF
            WRITE(Out,23)IGeo,Ctrl%Stats(NCyc,1),NCyc+1,Ctrl%Stats(NCyc,3)
            CALL PrintProtectR(Out)
            CLOSE(UNIT=Out,STATUS='KEEP')            
         ELSEIF(PrintFlags%Key==DEBUG_MEDIUM)THEN
 !           CALL Get(GM,Tag_O=TRIM(IntToChar(IGeo)))
 !           CALL PPrint(GM)
 !           CALL Delete(GM)
         ELSEIF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
!            CALL Get(GM,Tag_O=TRIM(IntToChar(IGeo)))
!            CALL PPrint(GM)
!            CALL Delete(GM)
            CALL OpenASCII(OutFile,Out)         
            CALL PrintProtectL(Out)
            Mssg='Convergence statistics for R'         &
               //TRIM(FunctionalName(Ctrl%Model(IBas))) &
               //'/'//TRIM(Ctrl%BName(IBas))//':'       
            WRITE(Out,*)TRIM(Mssg)
            WRITE(Out,11)
            WRITE(Out,40)
            WRITE(Out,41)
            DO I=0,NCyc
               WRITE(Out,42)I,Ctrl%InkFok,(Ctrl%Stats(I,K),K=1,4)
            ENDDO
            WRITE(Out,11)
            CALL PrintProtectR(Out)
            CLOSE(UNIT=Out,STATUS='KEEP')            
         ENDIF
      10 FORMAT(72('-'))
      11 FORMAT(72('='))
!
      21 FORMAT(' Ab initio ')
      22 FORMAT('  Config         <E_SCF>       Cycles      DeltaD ')
      23 FORMAT(1x,I5,2x,F20.8,5x,I2,7x,D10.4)
!
      40 FORMAT('SCF  ',2x,'Ink',7x,'Total ',12x,'Delta',7x,'Max    ',5x,'DIIS   ')
      41 FORMAT('Cycle',2x,'Fok',7x,'Energy',12x,'E_Tot',7x,'Delta P',5x,'Error  ')
      42 FORMAT(1x,I2,',   ',L2,', ',F20.8,4(', ',D10.4))
!
      END SUBROUTINE SCFSummry
!
!-----------------------------------------------------------------------------
#ifdef MMech
!
      SUBROUTINE   MM_CoulombEnergy(Ctrl)
      IMPLICIT NONE
      TYPE(SCFControls)  :: Ctrl
      TYPE(INT_VECT)     :: Stat
      LOGICAL            :: CalcMMForce
!
      IF(Ctrl%Grad==GRAD_NO_GRAD) THEN
        CalcMMForce=.FALSE.
      ELSE
        CalcMMForce=.TRUE.
      ENDIF
!
!compute the nuclear energy over the charge distribution given by
! MM charges        
!
       CALL New(Stat,3)
       Stat%I=Ctrl%Current
       CALL Put(Stat,'current')
       Call Delete(Stat)
!
#ifdef PERIODIC
       CtrlVect=SetCtrlVect(Ctrl,'MakingPFFT')
       CALL Invoke('MakePFFT',CtrlVect,MPIRun_O=.TRUE.)
#endif
!
       IF(CalcMMForce) THEN
         CtrlVect=SetCtrlVect(Ctrl,'ForceEvaluation')
         CALL Invoke('MakeRho',CtrlVect)
         CALL Invoke('QCTC',CtrlVect)
         CALL Invoke('JForce',CtrlVect)
       ELSE
         CtrlVect=SetCtrlVect(Ctrl,'PureMM')
         CALL Invoke('MakeRho',CtrlVect)
         CALL Invoke('QCTC',CtrlVect)
       ENDIF
!
      END SUBROUTINE MM_CoulombEnergy
#endif
!-------------------------------------------------------------
#ifdef MMech

   SUBROUTINE MM_ENERG(Ctrl)
   IMPLICIT NONE
   TYPE(ScfControls) :: Ctrl
   TYPE(CRDS) :: GMLoc
   REAL(DOUBLE)       :: MM_COUL,CONVF,E_C_EXCL,E_LJ_EXCL
   REAL(DOUBLE)       :: EELECT,ELJ,CONVF2,ETOTMM
   REAL(DOUBLE)       :: EBond,EAngle,ETorsion,EOutOfPlane
   REAL(DOUBLE)       :: LJCutOff
   TYPE(DBL_RNK2)     :: GrdMM
   TYPE(DBL_VECT)     :: GrdToT,GrdAux,GrdToTInt
   INTEGER :: I,J,I1,I2,MMNatms
   LOGICAL            :: CalcMMForce
!
      Call InitMMech()
!
      CONVF=1000.D0*JtoHartree/C_Avogadro
      CONVF2=1000.D0*JtoHartree/C_Avogadro/AngstromsToAU
      LJCutOff=10.D0 !!! in Angstroms
!     CurG=IntToChar(Ctrl%Current(3)) !!! Current Geom
!
      IF(Ctrl%Grad==GRAD_NO_GRAD) THEN
        CalcMMForce=.FALSE.
      ELSE
        CalcMMForce=.TRUE.
      ENDIF
!
! Load MM data
!
      CALL GET(GMLoc,'GM_MM'//CurGeom)
      MMNatms=GMLoc%Natms
!
! Initialize and Zero GrdMM for the case MMOnly
! (this initial gradient gets read by jforce later)
!
      IF(MMonly()) THEN
        CALL New(GrdAux,3*MMNatms)
        GrdAux%D(:)=Zero
          CALL Put(GrdAux,'GradE',Tag_O=CurGeom)
        CALL Delete(GrdAux)
      ENDIF
!
! GM_MM must be on HDF before CoulombEnergy is called
!
! MM CoulombEnergy is calculated here only for case MMOnly 
! Otherwise it is calculated in QCTC and put into HDF later
!
      IF(MMOnly()) Then
        CALL MM_CoulombEnergy(Ctrl)
        CALL GET(MM_COUL,'MM_COUL',Tag_O=CurGeom)
        MM_COUL=MM_COUL/CONVF
      ENDIF
!
! Do MM Covalent terms
!
      EBond=Zero
      EAngle=Zero
      ETorsion=Zero
      EOutOfPlane=Zero
      ELJ=Zero
      E_LJ_EXCL=Zero
      E_C_EXCL=Zero
!
! Convert coordinates into Angstroms
!
      GMLoc%Carts%D=GMLoc%Carts%D/AngstromsToAU
#ifdef PERIODIC
      GMLoc%PBC%BoxShape=GMLoc%PBC%BoxShape/AngstromsToAU
      GMLoc%AbCarts%D=GMLoc%AbCarts%D/AngstromsToAU
      GMLoc%AbBoxCarts%D=GMLoc%AbBoxCarts%D/AngstromsToAU
#endif
!
      IF(CalcMMForce) THEN
        CALL New(GrdMM,(/3,MMNatms/))
        GrdMM%D(:,:)=Zero
!   
          CALL Bond_Energy(EBond,GMLoc%Carts%D,GradLoc=GrdMM)
          CALL Angle_Energy(EAngle,GMLoc%Carts%D,GradLoc=GrdMM)
          CALL Torsion_Energy(ETorsion,GMLoc%Carts%D,GradLoc=GrdMM)
          CALL OutOfPlane_Energy(EOutOfPlane,GMLoc%Carts%D,GradLoc=GrdMM)
!   
          CALL ENERGY_LENNARD_JONES(GMLoc,ELJ,LJCutOff,GrdMM)
          CALL EXCL(GMLoc,E_LJ_EXCL,E_C_EXCL,GrdMM)
!   
      ELSE
!   
           CALL Bond_Energy(EBond,GMLoc%Carts%D)
           CALL Angle_Energy(EAngle,GMLoc%Carts%D)
           CALL Torsion_Energy(ETorsion,GMLoc%Carts%D)
           CALL OutOfPlane_Energy(EOutOfPlane,GMLoc%Carts%D)
!   
           CALL ENERGY_LENNARD_JONES(GMLoc,ELJ,LJCutOff)
           CALL EXCL(GMLoc,E_LJ_EXCL,E_C_EXCL)
!   
      ENDIF
!   
      GMLoc%Carts%D=AngstromsToAU*GMLoc%Carts%D
#ifdef PERIODIC
      GMLoc%PBC%BoxShape=AngstromsToAU*GMLoc%PBC%BoxShape
      GMLoc%AbCarts%D=AngstromsToAU*GMLoc%AbCarts%D
      GMLoc%AbBoxCarts%D=AngstromsToAU*GMLoc%AbBoxCarts%D
#endif
      IF(CalcMMForce) THEN
!
! Convert MM _gradients_ into atomic units and add the rest of forces
!
      CALL New(GrdTot,3*MMNatms)
      CALL Get(GrdTot,'GradE',Tag_O=CurGeom)
! 
      DO I=1,MMNatms
        I1=3*(I-1)+1
        I2=I1+2   
        GrdTot%D(I1:I2)=GrdTot%D(I1:I2)+KJPerMolPerAngstToHPerBohr*GrdMM%D(1:3,I)
      ENDDO
!
      CALL Put(GrdTot,'GradE',Tag_O=CurGeom)
!
! print forces in KJ/mol/A or H/Bohr
!
      CALL Print_Force(GMLoc,GrdTot,'GrdTot in au ')
      GrdTot%D(:)=GrdTot%D(:)/KJPerMolPerAngstToHPerBohr
      CALL Print_Force(GMLoc,GrdTot,'GrdTot in KJ/mol/A')
      GrdTot%D(:)=GrdTot%D(:)*KJPerMolPerAngstToHPerBohr
!
      CALL Delete(GMLoc)
      CALL Delete(GrdTot)
      CALL Delete(GrdMM)
!
    ENDIF
!
    CALL OpenASCII(OutFile,Out)
!
      WRITE(OUT,*) 'MM Energies in KJ/mol:'
      WRITE(OUT,*) 'EBond= ',EBond
      WRITE(OUT,*) 'EAngle= ',EAngle
      WRITE(OUT,*) 'ETorsion= ',ETorsion
      WRITE(OUT,*) 'EOutOfPlane= ',EOutOfPlane
      WRITE(OUT,*) 'E_Lennard_Jones    TOTAL= ',ELJ
      WRITE(OUT,*) 'E_Lennard_Jones EXCLUDED= ',E_LJ_EXCL
      WRITE(OUT,*) 'E_Lennard_Jones         = ',ELJ-E_LJ_EXCL
    IF(MMOnly()) Then
      WRITE(OUT,*) 'E_MM_Coulomb    TOTAL= ',MM_COUL
      WRITE(OUT,*) 'E_MM_Coulomb EXCLUDED= ',E_C_EXCL
      WRITE(OUT,*) 'E_MM_Coulomb         = ',MM_COUL-E_C_EXCL
      ETOTMM=MM_COUL+EBond+EAngle+ETorsion+EOutOfPlane+ELJ-E_LJ_EXCL-E_C_EXCL
      WRITE(OUT,*) 'E_Total= ',ETOTMM
    ENDIF
!
! convert covalent energies into atomic unit
!
    EBond=EBond*CONVF
    EAngle=EAngle*CONVF
    ETorsion=ETorsion*CONVF
    EOutOfPlane=EOutOfPlane*CONVF
    ELJ=ELJ*CONVF
    MM_COUL=MM_COUL*CONVF
    E_LJ_EXCL=E_LJ_EXCL*CONVF
    E_C_EXCL=E_C_EXCL*CONVF
!
      WRITE(OUT,*) 
      WRITE(OUT,*) 'MM Energies in atomic unit:'
      WRITE(OUT,*) 'EBond= ',EBond
      WRITE(OUT,*) 'EAngle= ',EAngle
      WRITE(OUT,*) 'ETorsion= ',ETorsion
      WRITE(OUT,*) 'EOutOfPlane= ',EOutOfPlane
      WRITE(OUT,*) 'E_Lennard_Jones    TOTAL= ',ELJ
      WRITE(OUT,*) 'E_Lennard_Jones EXCLUDED= ',E_LJ_EXCL
      WRITE(OUT,*) 'E_Lennard_Jones         = ',ELJ-E_LJ_EXCL
    IF(MMOnly()) Then
      WRITE(OUT,*) 'E_MM_Coulomb    TOTAL= ',MM_COUL
      WRITE(OUT,*) 'E_MM_Coulomb EXCLUDED= ',E_C_EXCL
      WRITE(OUT,*) 'E_MM_Coulomb         = ',MM_COUL-E_C_EXCL
      ETOTMM=MM_COUL+EBond+EAngle+ETorsion+EOutOfPlane+ELJ-E_LJ_EXCL-E_C_EXCL
      WRITE(OUT,*) 'E_Total= ',ETOTMM
      CALL Put(ETOTMM,'ETot',StatsToChar(Ctrl%Current))
    ENDIF
!
    CLOSE(UNIT=Out,STATUS='KEEP')
!
      CALL PUT(EBond,'MM_EBond',Tag_O=CurGeom)
      CALL PUT(EAngle,'MM_EAngle',Tag_O=CurGeom)
      CALL PUT(ETorsion,'MM_ETorsion',Tag_O=CurGeom)
      CALL PUT(EOutOfPlane,'MM_EOutOfPlane',Tag_O=CurGeom)
      CALL PUT(ELJ,'MM_ELJ',Tag_O=CurGeom)
      CALL PUT(E_LJ_EXCL,'E_LJ_EXCL',Tag_O=CurGeom)
      CALL PUT(E_C_EXCL,'E_C_EXCL',Tag_O=CurGeom)
!
END SUBROUTINE MM_ENERG
#endif
!-----------------------------------------------------
!
   SUBROUTINE Mulliken_Analysis(Ctrl)
     IMPLICIT NONE
     TYPE(Scfcontrols) :: Ctrl
     TYPE(BCSR)        :: Pmat,Smat,Tmat
     TYPE(BSET)        :: BS
     INTEGER           :: I,J,P,Q,S,II,JJ,SS,N,OUT
     TYPE(DBL_VECT)    :: Population,NuclCharge
     REAL(DOUBLE)      :: SUMM
     TYPE(ARGMT)       :: Args
!
     CALL Get(Args)
     CALL SetGlobalCtrlIndecies(Ctrl)
!
! Get basis set info
!
     CALL Get(BS,Tag_O=CurBase)
!
! Get BSiz and OffS
!
     NBasF=BS%NBasF
     PrintFlags%Mat=DEBUG_MATRICES
     CALL New(BSiz,natoms)
     CALL New(OffS,natoms)
     CALL Get(BSiz,'atsiz',Tag_O=CurBase)
     CALL Get(OffS,'atoff',Tag_O=CurBase)
!
! Get P and S matrices from HDF
!
     CALL Get(Smat,TrixFile('S',Args,Stats_O=Current))
     CALL Get(Pmat,TrixFile('D',Args,1,Stats_O=Current))
!
! Calc. population
!
     CALL New(Tmat)
     CALL Multiply(Pmat,Smat,Tmat)
!
     CALL New(Population,Tmat%Natms)
     Population%D(:)=Zero
!
     DO I=1,Tmat%Natms
	   S=BSiz%I(I)
       DO J=Tmat%RowPt%I(I),Tmat%RowPt%I(I+1)-1
	 P=Tmat%ColPt%I(J) ! col. index
	 IF(I==P) THEN
	   Q=Tmat%BlkPt%I(J) ! real array of block starts here
	      SS=0
	      SUMM=Zero
	   DO JJ=1,S
	      DO II=1,S
	      SS=SS+1
		IF(JJ==II) SUMM=SUMM+Tmat%MTrix%D(Q-1+SS)
	      ENDDO
	   ENDDO
	      Population%D(I)=SUMM
	 ENDIF
       ENDDO 
     ENDDO 
!
! For RHF, Population must be multiplied by two
!
     Population%D=Two*Population%D
!
! Get nuclear charges
!
     CALL New(NuclCharge,Natoms)
     CALL Get(NuclCharge,'atomicnumbers',Tag_O=CurGeom)
!
! Write Mulliken charges into output file and 
! save them into HDF
!
     CALL OpenASCII(OutFile,Out)
!
       WRITE(OUT,*)
       WRITE(OUT,*) 'Mulliken charges: '
       WRITE(OUT,100) SUM(NuclCharge%D-Population%D)
100  FORMAT('Total charge= ',F10.6)
       WRITE(OUT,*) ' Atom ','   Charge= '
       DO I=1,Natoms     
	 WRITE(OUT,200) I,NuclCharge%D(I)-Population%D(I)
200  FORMAT(I6,F10.3)
       ENDDO
!
     CALL Put(Population,'Population',Tag_O='#Geom_'//TRIM(CurGeom)//'_#Base_'//CurBase//'_#Cycle'//CurCycl)
!
     CLOSE(Out)
!
     CALL Delete(NuclCharge)
     CALL Delete(OffS)
     CALL Delete(BSiz)
     CALL Delete(Population)
     CALL Delete(Tmat)
     CALL Delete(Smat)
     CALL Delete(Pmat)
     CALL Delete(BS)
!
   END SUBROUTINE Mulliken_Analysis
!
!-----------------------------------------------------------------------------
END MODULE




