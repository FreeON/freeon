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
     SUBROUTINE OneSCF(Ctrl,Sum_O)
        TYPE(SCFControls)  :: Ctrl
        INTEGER            :: ISCF,ICyc
        LOGICAL,OPTIONAL   :: Sum_O
        LOGICAL            :: Summry
        IF(PRESENT(Sum_O))THEN
           Summry=Sum_O
        ELSE
           Summry=.TRUE.
        ENDIF
        ISCF=Ctrl%Current(1)       
        DO ICyc=ISCF,MaxSCFs
           Ctrl%Current(1)=ICyc
           CALL SetGlobalCtrlIndecies(Ctrl)
           IF(ICyc==0)CALL OneEMats(Ctrl)
           CALL SCFCycle(Ctrl)          
           IF(ConvergedQ(Ctrl))THEN
              IF(Summry)CALL SCFSummry(Ctrl)
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
!      CALL CleanScratch(Ctrl)
    END SUBROUTINE SCFCycle
!==========================================================================
!   Build a density by hook or by crook
!==========================================================================
    SUBROUTINE DensityBuild(Ctrl)
      TYPE(SCFControls)    :: Ctrl
      INTEGER,DIMENSION(3) :: Tmp
!---------------------------------------------------------------------------------------
      IF(CCyc==0.AND.CGeo==1.AND.CBas>1)THEN
!        Basis set switch
         CALL LogSCF(Ctrl%Current,' Switching basis sets from ' &
                                //TRIM(Ctrl%BName(PBas))//' to ' &
                                //TRIM(Ctrl%BName(CBas))//'.')
!        Create density from last SCF, fake out MakeRho....
         Tmp=Ctrl%Current  
         Ctrl%Current=Ctrl%Previous       
         CtrlVect=SetCtrlVect(Ctrl,'BasisSetSwitch')
         Ctrl%Current=Tmp
         CALL Invoke('MakeRho',CtrlVect)
      ELSEIF(CCyc==0.AND.(CGeo/=PGeo))THEN
!-------------------------------------------------------------------------
!        Projection of density matrix between geometries
!         CALL LogSCF(Ctrl%Current,' Geometry projection from configuration #' &
!                                //TRIM(PrvGeom)//' to configuration# ' &
!                                //TRIM(CurGeom)//'.')
!        Create density from last SCF 
!         CtrlVect=SetCtrlVect(Ctrl,'Project')
!-------------------------------------------------------------------------
!        Extrapolation of density matrix between geometries
         CALL LogSCF(Ctrl%Current,' Geometry extrapolation from configuration #' &
                                //TRIM(PrvGeom)//' to configuration# ' &
                                //TRIM(CurGeom)//'.')
!        Create density from last SCF 
         CtrlVect=SetCtrlVect(Ctrl,'Extrapolate')
         CALL Invoke('P2Use',CtrlVect)
         CALL Invoke('MakeRho',CtrlVect)
     ELSEIF(CCyc==0.AND.CBas==1.AND.CGeo==1)THEN
!        Start up density from atomic superposition in diagonal minimial basis 
         CALL LogSCF(Ctrl%Current,' Density guess from superposition of AOs ')
!        Create density from guess
         CtrlVect=SetCtrlVect(Ctrl,'DensitySuperposition')
         CALL Invoke('P2Use',CtrlVect)
         CALL Invoke('MakeRho',CtrlVect)
      ELSE
!        Regular old density build
         CALL LogSCF(Ctrl%Current,'Direct SCF')
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
      LOGICAL            :: DoDIIS
!-----------------------------------------------------------
      DoDIIS=CCyc>0
!      IF(DoDIIS)THEN
!         CALL LogSCF(Current,'Fock build with Pulay DIIS')
!      ELSE
!         CALL LogSCF(Current,'Fock build with Pulay DIIS')
!      ENDIF
      Modl=Ctrl%Model(CBas)
      CALL Invoke('QCTC',CtrlVect)
      IF(HasDFT(Modl))CALL Invoke('HiCu',CtrlVect)
      IF(HasHF(Modl)) CALL Invoke('ONX',CtrlVect)
      CALL Invoke('FBuild',CtrlVect,MPIRun_O=.TRUE.)
      IF(DoDIIS) &
      CALL Invoke('DIIS',CtrlVect,MPIRun_O=.TRUE.)
    END SUBROUTINE FockBuild
!==========================================================================
!   Solve the SCF equations 
!==========================================================================
    SUBROUTINE SolveSCF(Ctrl)
      TYPE(SCFControls)  :: Ctrl
!-----------------------------------------------------------
!      CtrlVect=SetCtrlVect(Ctrl,'SCFSolve')
      IF(Ctrl%Method(CBas)==RH_R_SCF)THEN
!         CALL LogSCF(Current,'Solving SCF equations with Roothaan-Hall')
         CALL Invoke('RHeqs',CtrlVect)
      ELSE
!         CALL LogSCF(Current,'Solving SCF equations with SDMM')
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
      CtrlVect=SetCtrlVect(Ctrl,'OneElectron')
      CALL Invoke('MakeS',CtrlVect,MPIRun_O = .TRUE.)

      IF(Ctrl%Method(Ctrl%Current(2))==RH_R_SCF)THEN
         CALL Invoke('LowdinO',  CtrlVect)
      ELSEIF(Ctrl%Method(Ctrl%Current(2))==SDMM_R_SCF)THEN
         CALL Invoke('AInv',CtrlVect)
      ELSE
         CALL MondoHalt(MISC_ERROR,' Neither FactoredS or LOrthog invoked in OneEMats ')
      ENDIF
      CALL Invoke('MakeT',CtrlVect)
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
      TYPE(SCFControls)              :: Ctrl
      CHARACTER(LEN=4)               :: Action
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: RemTag,RemoveFile
      INTEGER                        :: IGeo,IBas,ICyc
!---------------------------------------------------------------------------------------      
      ICyc=Ctrl%Current(1)
      IBas=Ctrl%Current(2)
      IGeo=Ctrl%Current(3)
      IF(Action=='Cycl') THEN
         RemTag=TRIM(ScrName)//'_Geom#'//TRIM(IntToChar(IGeo))  & 
                             //'_Base#'//TRIM(IntToChar(IBas))  &
                             //'_Cycl#'//TRIM(IntToChar(ICyc-1)) 
         RemoveFile=TRIM(RemTag)//'.Rho'                                
         CALL SYSTEM('/bin/rm '//RemoveFile)
         RemoveFile=TRIM(RemTag)//'.J'                                
         CALL SYSTEM('/bin/rm '//RemoveFile)
         RemoveFile=TRIM(RemTag)//'.K'                                
         CALL SYSTEM('/bin/rm '//RemoveFile)
         RemoveFile=TRIM(RemTag)//'.D'                                
         CALL SYSTEM('/bin/rm '//RemoveFile)
         RemoveFile=TRIM(RemTag)//'.F'                                
         CALL SYSTEM('/bin/rm '//RemoveFile)
      ELSEIF(Action=='Base') THEN
         RemoveFile=TRIM(ScrName) //'_Geom#'//TRIM(IntToChar(IGeo))    & 
                                    //'_Base#'//TRIM(IntToChar(IBas-1)) //'*' 
         CALL SYSTEM('/bin/rm '//RemoveFile)  
      ELSEIF(Action=='Geom') THEN
         RemoveFile=TRIM(ScrName) //'_Geom#'//TRIM(IntToChar(IGeo-1))//'*' 
         CALL SYSTEM('/bin/rm '//RemoveFile) 
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
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
!-----------------------------------------------------------------------
         ConvergedQ=.FALSE.
         IF(CCyc==0)RETURN
!-------------------------------------
!        Gather convergence parameters
         CALL OpenHDF(Ctrl%Info)
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
         DMaxQ=dDMax/ABS(DMaxB)
         DIISQ=dDIIS/ABS(DIISB)
!------------------------------
!        Check for convergence
         ETest=ETol(Ctrl%AccL(CBas))
         DTest=DTol(Ctrl%AccL(CBas))
!        Check for absolute convergence below thresholds
         IF(dDMax<dTest.AND.ETotQ<ETest)THEN
            Mssg='Normal SCF convergence.'
            ConvergedQ=.TRUE.
         ENDIF
!        Check to see if convergence is in an asymptotic regime
         IF(DIISB<1.D-2.AND.DMAXB<5.D-1)THEN
!           Look for non-decreasing error due to incomplete numerics
            IF(DIISQ<1.D-1.AND.DMaxQ<1.D-1.AND.CCyc>2)THEN                
               IF(DIISB>DIISA.AND.DMAXB>DMAXA)THEN
                  Mssg='SCF hit DIIS/DMax increase.'
                  ConvergedQ=.TRUE.
               ENDIF
            ENDIF
!           Look for convergence stall-outs 
            IF(DIISQ<9.D-2.AND.DMAXQ<1.D-2)THEN
               WRITE(*,*)' DIISQ = ',DIISQ
               WRITE(*,*)' DMAXQ = ',DMAXQ
               Mssg='SCF convergence stalled.'
               ConvergedQ=.TRUE.
            ENDIF 
         ENDIF
         IF(ConvergedQ.AND.PrintFlags%Key>DEBUG_NONE)THEN
            WRITE(*,*)TRIM(Mssg)             
            IF(PrintFlags%Key>DEBUG_MINIMUM)THEN
               CALL OpenASCII(OutFile,Out)
               WRITE(Out,*)TRIM(Mssg)             
               CLOSE(Out)
            ENDIF
         ENDIF
!--------------------------------------------------------
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
END MODULE




