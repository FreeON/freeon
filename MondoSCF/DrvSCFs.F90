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
#endif
  IMPLICIT NONE 
  CONTAINS
!========================================================================================
!
!========================================================================================
     SUBROUTINE OneSCF(Ctrl,Sum_O)
        TYPE(SCFControls)  :: Ctrl
!       TYPE(CRDS)         :: GM_MM
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
! Do QM Calc on the SCF/DFT level
!
        DO ICyc=0,MaxSCFs
           Ctrl%Current(1)=ICyc
           CALL SetGlobalCtrlIndecies(Ctrl)
           IF(ICyc==0)CALL OneEMats(Ctrl)
!          Do an SCF cycle
           CALL SCFCycle(Ctrl)          
           IF(ConvergedQ(Ctrl))THEN
              IF(Summry)CALL SCFSummry(Ctrl)
              CALL VisDX(Ctrl)
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
      IF(CCyc==0.AND.CBas==PBas.AND.CGeo/=1.AND. &
         Ctrl%Extrap>EXTRAP_GEOM_RSTRT)THEN
         IF(Ctrl%Extrap==EXTRAP_GEOM_PRJCT)THEN
!           Projection of density matrix between geometries
            CALL LogSCF(Ctrl%Current,'Geometry projection from configuration #' &
                                   //TRIM(PrvGeom)//' to configuration# ' &
                                   //TRIM(CurGeom)//'.')
!           Create density from last SCF 
            CtrlVect=SetCtrlVect(Ctrl,'Project')
         ELSEIF(Ctrl%Extrap==EXTRAP_GEOM_INTRP)THEN
!           Extrapolation of density matrix between geometries
            CALL LogSCF(Ctrl%Current,'Geometry extrapolation from configuration #' &
                                   //TRIM(PrvGeom)//' to configuration# ' &
                                   //TRIM(CurGeom)//'.')
!           Create density from last SCFs DM (ICyc+1)

            Ctrl%Previous(1)=Ctrl%Previous(1)+1
            CtrlVect=SetCtrlVect(Ctrl,'Extrapolate')
         ENDIF
         CALL Invoke('P2Use',CtrlVect,MPIRun_O=.TRUE.)
         CALL Invoke('MakeRho',CtrlVect)
         CALL CleanScratch(Ctrl,'CleanLastGeom')
      ELSEIF(Ctrl%Rest)THEN
!        Restart from a previous density
         CALL LogSCF(Ctrl%Current,'Restarting from '//TRIM(Ctrl%OldInfo),.TRUE.)
         CtrlVect=SetCtrlVect(Ctrl,'Restart')
         CALL Invoke('P2Use',CtrlVect,MPIRun_O=.TRUE.)
         CALL Invoke('MakeRho',CtrlVect)
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
      IF(HasDFT(Modl))CALL Invoke('HiCu',CtrlVect)
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
      ELSE
         CALL LogSCF(Current,'Solving SCF equations with SDMM')
         CALL Invoke('SDMM',CtrlVect,MPIRun_O=.TRUE.) 
      ENDIF
      CALL Invoke('SCFstats',CtrlVect,MPIRun_O=.TRUE.)                                     
    END SUBROUTINE SolveSCF
!========================================================================================
!
!========================================================================================
    SUBROUTINE OneEMats(Ctrl)
      TYPE(SCFControls)             :: Ctrl
      INTEGER                       :: ISet,I
      CALL LogSCF(Current,'One-electron matrices.',.TRUE.)
      CtrlVect=SetCtrlVect(Ctrl,'OneElectron')
      CALL Invoke('MakeS',CtrlVect,MPIRun_O = .TRUE.)
      IF(Ctrl%Method(Ctrl%Current(2))==RH_R_SCF)THEN
         CALL Invoke('LowdinO',  CtrlVect)
      ELSEIF(Ctrl%Method(Ctrl%Current(2))==SDMM_R_SCF)THEN
         CALL Invoke('AInv',CtrlVect)
      ELSE
         CALL MondoHalt(MISC_ERROR,' Neither FactoredS or LOrthog invoked in OneEMats ')
      ENDIF
      CALL Invoke('MakeT',CtrlVect,MPIRun_O=.TRUE.)
    END SUBROUTINE OneEMats

    SUBROUTINE VisDX(Ctrl)
      TYPE(SCFControls)             :: Ctrl
      IF(Ctrl%Vis/=VIS_DX_RHOPOT)RETURN
      Ctrl%Current(1)=Ctrl%Current(1)+1
      CtrlVect=SetCtrlVect(Ctrl,'Visualization')
      Ctrl%Current(1)=Ctrl%Current(1)-1
      CALL Invoke('MakeRho',CtrlVect)
      CALL Invoke('PotMapRho',CtrlVect)
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
         RemPrv=TRIM(ScrName)//'_Geom#' &
              //TRIM(CurGeom)//'_Base#' &
              //TRIM(CurBase)//'_Cycl#' &
              //TRIM(PrvCycl)
         RemoveFile=TRIM(RemCur)//'.Rho'                                
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
         RemoveFile=TRIM(ScrName) //'_Geom#'//TRIM(PrvGeom)//'*.Rho' 
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
         RemoveFile=TRIM(ScrName) //'*_Base#'//TRIM(PrvGeom)//'*.Rho' 
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
         CALL OpenHDF(Ctrl%Info)
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
         CALL CloseHDF()
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
         IF(DMaxB<dTest.AND.ETotQ<ETest.AND.ETotB<ETotA)THEN
            Mssg='Normal SCF convergence.'
            ConvergedQ=.TRUE.
         ENDIF
!        Accept convergence from wrong side if thresholds are tightend.
         IF(DMaxB<dTest*6D-1.AND.ETotQ<ETest*1D-1)THEN
            Mssg='Normal SCF convergence.'
            ConvergedQ=.TRUE.
         ENDIF
!        Check to see if convergence is in an asymptotic regime
         IF(DIISB<1.D-2.AND.DMaxB<5.D-1)THEN
!           Look for non-decreasing error due to incomplete numerics
            IF(DIISQ<1.D-1.AND.DMaxQ<1.D-1.AND.CCyc>2)THEN                
               IF(DIISB>DIISA.AND.DMaxB>DMaxA)THEN
                  Mssg='SCF hit DIIS/DMax increase.'
                  ConvergedQ=.TRUE.
               ENDIF
            ENDIF
!           Look for convergence stall-outs 
            IF(DIISQ<9.D-2.AND.DMaxQ<1.D-2)THEN
               Mssg='SCF convergence stalled.'
               ConvergedQ=.TRUE.
            ENDIF 
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
            IF(PrintFlags%Key>DEBUG_MINIMUM)THEN
               CALL OpenASCII(OutFile,Out)
               WRITE(Out,*)TRIM(Mssg)             
               CLOSE(Out)
            ENDIF
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
 !           CALL OpenHDF(InfFile)
 !           CALL Get(GM,Tag_O=TRIM(IntToChar(IGeo)))
 !           CALL CloseHDF()
 !           CALL PPrint(GM)
 !           CALL Delete(GM)
         ELSEIF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
!            CALL OpenHDF(InfFile)
!            CALL Get(GM,Tag_O=TRIM(IntToChar(IGeo)))
!            CALL CloseHDF()
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
      SUBROUTINE   MM_COULOMBENERGY(Ctrl)
      IMPLICIT NONE
      TYPE(SCFControls)  :: Ctrl
      TYPE(INT_VECT)     :: Stat

!compute the nuclear energy over the charge distribution given by
! MM charges        
!
       CALL OpenHDF(InfFile)
       CALL New(Stat,3)
       Stat%I=Ctrl%Current
       CALL Put(Stat,'current')
       Call Delete(Stat)
       CALL CLOSEHDF()
!
       CtrlVect=SetCtrlVect(Ctrl,'PureMM')
       CALL Invoke('MakeRho',CtrlVect)
       CALL Invoke('QCTC',CtrlVect)
!
      END SUBROUTINE MM_COULOMBENERGY
#endif
!-------------------------------------------------------------
#ifdef MMech

   SUBROUTINE MM_ENERG(Ctrl)
   IMPLICIT NONE
   TYPE(Scfcontrols) Ctrl
   REAL(DOUBLE)       :: MM_COUL,CONVF,E_C_EXCL,E_LJ_EXCL
   REAL(DOUBLE)       :: EBOND,EELECT,ELJ  
   CHARACTER(LEN=3)   :: Cur
   TYPE(INT_VECT)     :: Stat
!  TYPE(DBL_RNK2)     :: Grd1,Grd2
   INTEGER :: I
!
       CONVF=1000.D0*JtoHartree/C_Avogadro
       CALL New(Stat,3)
       Stat%I=Ctrl%Current
       Cur=IntToChar(Stat%I(2))
       Call Delete(Stat)
!
! GM_MM must be on HDF before COULOMBENERGY is called
!
! MM coulombenergy is calculated here only for case MMOnly 
! Otherwise it is calculated in QCTC and put into HDF later
!
       IF(Ctrl%Mechanics(1).AND..NOT.Ctrl%Mechanics(2)) Then
         CALL MM_COULOMBENERGY(Ctrl)
         CALL OpenHDF(Ctrl%Info)
         CALL GET(MM_COUL,'MM_COUL',Tag_O=Cur)
         CALL CloseHDF()
         MM_COUL=MM_COUL/CONVF
       ENDIF
!
! Do MM Covalent terms
!
       CALL ENERGY_BOND ( EBOND, VIRIAL,InfFile=InfFile )
       CALL ENERGY_ANGLE ( EANGLE ,InfFile=InfFile)
       CALL ENERGY_DIHEDRAL ( EDIHEDRAL ,InfFile=InfFile)
       CALL ENERGY_IMPROPER ( EIMPROPER ,InfFile=InfFile)
       CALL ENERGY_LENNARD_JONES(ELJ,Ctrl%Current(2),7.D0)
!
! calculate exclusion energies
!
       CALL EXCL(MM_Natoms,Cur,InfFile,E_LJ_EXCL,E_C_EXCL)
!
!        call new(grd1,(/3,MM_Natoms/))
!        call new(grd2,(/3,MM_Natoms/))
!        grd1%d(:,:)=zero
!        grd2%d(:,:)=zero
!      CALL ENERGY_NON_BONDING_CALCULATE( EELECT, ELJ, VIRIAL,grd1%D(:,:))
!      write(*,*) 'elj aft traditional= ',elj,Ctrl%Current(2)
!
!   do i=1,mm_Natoms
!     write(*,*) 'grd= ',i,grd1%D(1,i),grd2%D(1,i),grd1%D(1,i)-grd2%D(1,i)
!     write(*,*) 'grd= ',i,grd1%D(2,i),grd2%D(2,i),grd1%D(2,i)-grd2%D(2,i)
!     write(*,*) 'grd= ',i,grd1%D(3,i),grd2%D(3,i),grd1%D(3,i)-grd2%D(3,i)
!   enddo
!        call delete(grd1)
!        call delete(grd2)
!
       CALL OpenASCII(OutFile,Out)
!
       write(out,*) 'MM Energies in KJ/mol:'
       write(out,*) 'Ebond= ',ebond
       write(out,*) 'Eangle= ',eangle
       write(out,*) 'Edihedral= ',edihedral
       write(out,*) 'Eimproper= ',eimproper
       write(out,*) 'E_Lennard_Jones    TOTAL= ',ELJ
       write(out,*) 'E_Lennard_Jones EXCLUDED= ',E_LJ_EXCL
       write(out,*) 'E_Lennard_Jones         = ',ELJ-E_LJ_EXCL
       IF(Ctrl%Mechanics(1).AND..NOT.Ctrl%Mechanics(2)) Then
       write(out,*) 'E_MM_Coulomb    TOTAL= ',MM_COUL
       write(out,*) 'E_MM_Coulomb EXCLUDED= ',E_C_EXCL
       write(out,*) 'E_MM_Coulomb         = ',MM_COUL-E_C_EXCL
       write(out,*) 'E_Total= ',MM_COUL+ebond+eangle+edihedral+eimproper+ELJ-E_LJ_EXCL-E_C_EXCL
       ENDIF
!
! convert covalent energies into atomic unit
!
       EBOND=EBOND*CONVF
       EANGLE=EANGLE*CONVF
       EDIHEDRAL=EDIHEDRAL*CONVF
       EIMPROPER=EIMPROPER*CONVF
       ELJ=ELJ*CONVF
       MM_COUL=MM_COUL*CONVF
       E_LJ_EXCL=E_LJ_EXCL*CONVF
       E_C_EXCL=E_C_EXCL*CONVF
!
       write(out,*) 
       write(out,*) 'MM Energies in atomic unit:'
       write(out,*) 'Ebond= ',ebond
       write(out,*) 'Eangle= ',eangle
       write(out,*) 'Edihedral= ',edihedral
       write(out,*) 'Eimproper= ',eimproper
       write(out,*) 'E_Lennard_Jones    TOTAL= ',ELJ
       write(out,*) 'E_Lennard_Jones EXCLUDED= ',E_LJ_EXCL
       write(out,*) 'E_Lennard_Jones         = ',ELJ-E_LJ_EXCL
       IF(Ctrl%Mechanics(1).AND..NOT.Ctrl%Mechanics(2)) Then
       write(out,*) 'E_MM_Coulomb    TOTAL= ',MM_COUL
       write(out,*) 'E_MM_Coulomb EXCLUDED= ',E_C_EXCL
       write(out,*) 'E_MM_Coulomb         = ',MM_COUL-E_C_EXCL
       write(out,*) 'E_Total= ',MM_COUL+ebond+eangle+edihedral+eimproper+ELJ-E_LJ_EXCL-E_C_EXCL
       ENDIF
!
       GM_MM%NAtms = MM_Natoms
!
        CALL OpenHDF(Ctrl%Info)
        CALL PUT(EBOND,'MM_EBOND',Tag_O=Cur)
        CALL PUT(EANGLE,'MM_EANGLE',Tag_O=Cur)
        CALL PUT(EDIHEDRAL,'MM_EDIHEDRAL',Tag_O=Cur)
        CALL PUT(EIMPROPER,'MM_EIMPROPER',Tag_O=Cur)
        CALL PUT(VIRIAL,'MM_VIRIAL',Tag_O=Cur)
        CALL PUT(ELJ,'MM_ELJ',Tag_O=Cur)
        CALL PUT(E_LJ_EXCL,'E_LJ_EXCL',Tag_O=Cur)
        CALL PUT(E_C_EXCL,'E_C_EXCL',Tag_O=Cur)
        CALL CloseHDF()
!
        CLOSE(UNIT=Out,STATUS='KEEP')
!
   END SUBROUTINE MM_ENERG
#endif

!-----------------------------------------------------------------------------
END MODULE




