MODULE ChkSCFs
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalCharacters
   USE ProcessControl
   USE InOut
   USE Parse
   USE PrettyPrint
   USE SCFLocals
   USE ParsingKeys
   USE Functionals
   IMPLICIT NONE
   CONTAINS
!==============================================================================
!     Simple backtracking algorithm for quasi-Newton minimization from NR 9.7
!==============================================================================
      FUNCTION ChkStep(Ctrl,GM,Bak)
         TYPE(SCFControls)              :: Ctrl
         TYPE(CRDS)                     :: GM
         INTEGER                        :: ChkStep,Bak
         REAL(DOUBLE),DIMENSION(20),SAVE :: StepSz
         REAL(DOUBLE),DIMENSION(20),SAVE :: G
         REAL(DOUBLE)                   :: E0,E1,GPrime,GradEDotDeltaX,EUncert
!-----------------------------------------------------------------------
         WRITE(*,*)' ChkStep ChkStep ChkStep ChkStep ChkStep ChkStep ChkStep '
         ChkStep=0
!        Open the InfFile
         CALL OpenHDF(Ctrl%Info)
!        Get the current geometry
         CALL Get(GM,Tag_O=IntToChar(Ctrl%Current(3)))
!        Get the previous and current total energy
         CALL Get(E0,'Etot',Tag_O=StatsToChar(Ctrl%Previous))
         CALL Get(E1,'Etot',Tag_O=StatsToChar(Ctrl%Current))
         GM%ETotal=E1
         CALL Get(GM%GradRMS,'RMSGrad',Tag_O=IntToChar(Ctrl%Current(3)))
         CALL Get(GM%GradMax,'MaxGrad',Tag_O=IntToChar(Ctrl%Current(3)))
         GM%Confg=Ctrl%Current(3)
!        Check for non-decreasing stepsize 
         CALL Get(GradEDotDeltaX,'GradEDotDeltaX')
         WRITE(*,*)'E'//TRIM(StatsToChar(Ctrl%Previous))//' = '//TRIM(DblToMedmChar(E0))
         WRITE(*,*)'E'//TRIM(StatsToChar(Ctrl%Current))//' = '//TRIM(DblToMedmChar(E1))
         EUncert=1.D1*ETol(Ctrl%AccL(Ctrl%Current(2)))
         IF(E1<E0+EUncert)THEN
            Bak=1
            ChkStep=1
            StepSz(Bak)=One
            CALL PPrint(GM,GeoFile,Geo,'XYZ')
         ELSE
            ChkStep=0
            Bak=Bak+1
            G(1)=E0         
            G(Bak)=E1
            CALL Get(GPrime,'GradEDotNewStep')
            IF(Bak==2)THEN
!              Improved quadratic estimate
               StepSz(2)=MAX(1.D-1,-Half*GPrime/(G(2)-G(1)-GPrime))
            ELSEIF(Bak==3)THEN
               StepSz(Bak)=5.D-2
            ELSEIF(Bak==4)THEN
               StepSz(Bak)=1.D-2
            ELSE
!              Accept last step, cross fingers and get on with it
               Bak=1
               ChkStep=1
               StepSz(Bak)=One
               CALL PPrint(GM,GeoFile,Geo,'XYZ')
            ENDIF
         ENDIF
         WRITE(*,*)' StepSize'//TRIM(IntToChar(Bak))//' = ',DblToShrtChar(StepSz(Bak))
         CALL Put(StepSz(Bak),'StepSize')
         CALL CloseHDF()
      END FUNCTION ChkStep
!==============================================================================
!
!==============================================================================
      FUNCTION ConvergedQ(Ctrl)
         TYPE(SCFControls) :: Ctrl
         INTEGER           :: IBas,ICyc,IGeo
         INTEGER,SAVE      :: DoubleChk
         LOGICAL           :: ConvergedQ
         REAL(DOUBLE)      :: ETot,ETotA,ETotB,DMax,      &
                              DeltaE,DeltaD,DIISA,DIISB,  &
                              DMaxA,DMaxB,Delta_DMax, Delta_ETot
         REAL(DOUBLE)      :: ConvQ_DIIS,ConvQ_ETot,ConvQ_DMax
         TYPE(INT_VECT)    :: Stat
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,CrntTag,PrevTag
!-----------------------------------------------------------------------
         ConvergedQ=.FALSE.
!        Open the InfFile
         CALL OpenHDF(Ctrl%Info)
!        Put the current status
         CALL New(Stat,3)
         Stat%I=Ctrl%Previous
         CALL Put(Stat,'PreviousStatus')
         Stat%I=Ctrl%Current
         CALL Put(Stat,'CurrentStatus')
         CALL Delete(Stat)
!        Notation ...
         ICyc=Ctrl%Current(1)
         IBas=Ctrl%Current(2)
         IGeo=Ctrl%Current(3)
!
         IF(ICyc==1)THEN
            DoubleChk=0
         ENDIF
!---------------------------------------------------------------------
!        Gather convergence parameters
!
         DIISA=1.0D4
         DIISB=1.0D1
         ETotA=1.0D4
         ETotB=1.0D1
         DMaxA=1.0D4
         DMaxB=1.0D1         
         IF(.NOT.(ICyc==0.AND.IBas==1))THEN
            CrntTag=StatsToChar(Ctrl%Current)
            PrevTag=StatsToChar(Ctrl%Previous)
            IF(ICyc>1)THEN
               CALL Get(DMaxA,'DMax',Tag_O=TRIM(PrevTag))
!               WRITE(*,*)' DMAX',TRIM(PrevTag),' = ',DMaxA
               CALL Get(EtotA,'Etot',Tag_O=TRIM(PrevTag))
            ENDIF
            CALL Get(DMaxB,'DMax',Tag_O=TRIM(CrntTag))
            CALL Get(EtotB,'Etot',Tag_O=TRIM(CrntTag))
            IF(ICyc>0) &
            CALL Get(DIISB,'diiserr',Tag_O=TRIM(CrntTag))
            IF(ICyc>1) &
            CALL Get(DIISA,'diiserr',Tag_O=TRIM(PrevTag))
         ENDIF
!--------------------------------------------------------------------
!        IO 
!
         Mssg=' DIIS Err = '//TRIM(DblToShrtChar(DIISB)) &
            //', MAX/P = '//TRIM(DblToShrtChar(DMAXB))   &
           //', <SCF> = '//TRIM(DblToMedmChar(ETotB))
         WRITE(*,*)TRIM(Mssg)
!---------------------------------------------------------------------
!        Check for convergence
!
         Delta_DMax=ABS(DMaxA-DMaxB)
         Delta_ETot=ABS(ETotA-ETotB)
         ConvQ_ETot=ABS((ETotA-ETotB)/ETotB)
         ConvQ_DMax=ABS((DMaxA-DMaxB)/DMaxB)
        ConvQ_DIIS=ABS((DIISA-DIISB)/DIISB)
!         WRITE(*,*)' ConvQ_DIIS = ',ConvQ_DIIS,' ConvQ_DMAX = ',ConvQ_DMAX
!        Could happen ...
         IF(ConvQ_ETot<1.D-14)THEN
            ConvergedQ=.TRUE.
            WRITE(*,*)' Met Convergence criteria A '
         ENDIF
!        Check to see if convergence is in an asymptotic regime
!         WRITE(*,*)' DIISB = ',DIISB
!         WRITE(*,*)' DMAXB = ',DMAXB
         IF(DIISB<1.D-2.AND.DMAXB<5.D-1)THEN
!           Look for non-decreasing error due to incomplete numerics
            IF(ConvQ_DIIS<7.D-1.AND.ConvQ_DMax<7.D-1.AND.ICyc>2)THEN                

!                     WRITE(*,*)' ConvQ_DIIS = ',ConvQ_DIIS
!                     WRITE(*,*)' ConvQ_DMAX = ',ConvQ_DMAX
!                     WRITE(*,*)' DIIS I     = ',DIISA
!                     WRITE(*,*)' DIIS I+1   = ',DIISB
!                     WRITE(*,*)' DMAX I     = ',DMAXA
!                     WRITE(*,*)' DMAX I+1   = ',DMAXB


               IF(DIISB>DIISA.AND.DMAXB>DMAXA)THEN
                  DoubleChk=DoubleChk+1              
!                  WRITE(*,*)' DOUBLE CHECK = ',DoubleChk
                  IF(DoubleChk==1)THEN
                     ConvergedQ=.TRUE.
                     WRITE(*,*)' Met Convergence criteria B: DIIS/DMAX increase.'
                  ENDIF
               ENDIF
            ENDIF
!           Look for convergence stall-outs 
            IF(ConvQ_DIIS<9.D-2.AND.ConvQ_DMAX<9.D-2)THEN
               DoubleChk=DoubleChk+1              
               IF(DoubleChk==1)THEN
                  WRITE(*,*)' Met Convergence criteria C: DIIS/DMAX stall out.'
                  ConvergedQ=.TRUE.
               ENDIF
            ENDIF 
!           Check for absolute convergence below thresholds
            IF(ICyc>0.AND.Delta_DMax<DTol(Ctrl%AccL(IBas)) &
                     .AND.ConvQ_ETot<ETol(Ctrl%AccL(IBas)))THEN 
               WRITE(*,*)' Met specified convergence criteria.'
!               WRITE(*,*)' DeltaDDMAX < ',DTol(Ctrl%AccL(IBas))
!               WRITE(*,*)' RelErrETOT < ',ETol(Ctrl%AccL(IBas))
               ConvergedQ=.TRUE.
            ENDIF
         ENDIF             
!--------------------------------------------------------
!        Load statistics 
!
         Ctrl%EErr(IBas)=ConvQ_ETot
         Ctrl%Stats(ICyc,1)=ETotB
         Ctrl%Stats(ICyc,2)=ConvQ_ETot
         Ctrl%Stats(ICyc,3)=Delta_DMax
         Ctrl%Stats(ICyc,4)=DIISB
         Ctrl%Stats(ICyc,5)=Zero
!        Close the InfFile
         CALL CloseHDF()
!
         Ctrl%NCyc(IBas)=ICyc
         Ctrl%Previous=Ctrl%Current

!if(icyc<8)convergedq=.false.

         RETURN
!
      END FUNCTION ConvergedQ
!
      FUNCTION StatsToChar(Stats) RESULT(StatString)
         INTEGER,DIMENSION(3) :: Stats
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: StatString
         StatString='_'//TRIM(IntToChar(Stats(3)))  &
                   //'_'//TRIM(IntToChar(Stats(2))) &
                   //'_'//TRIM(IntToChar(Stats(1)))
      END FUNCTION StatsToChar
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
!      
END MODULE
