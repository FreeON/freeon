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
      FUNCTION ConvergedQ(Ctrl)
         TYPE(SCFControls) :: Ctrl
         INTEGER           :: IBas,ICyc,IGeo
         LOGICAL           :: ConvergedQ
         REAL(DOUBLE)      :: ETot,ETotA,ETotB,DMax,      &
                              DeltaE,DeltaD,DIISA,DIISB,  &
                              DMaxA,DMaxB,Delta_DMax, Delta_ETot
         REAL(DOUBLE)      :: ConvQ_DIIS,ConvQ_ETot,ConvQ_DMax
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg,CrntTag,PrevTag
!-----------------------------------------------------------------------
         ConvergedQ=.FALSE.
!        Open the InfFile
         CALL OpenHDF(Ctrl%Info)
!        Put the current status
         CALL Put(Ctrl%Previous,'PreviousStatus')
         CALL Put(Ctrl%Current,'CurrentStatus')
!        Notation ...
         ICyc=Ctrl%Current%I(1)
         IBas=Ctrl%Current%I(2)
         IGeo=Ctrl%Current%I(3)
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
            IF((Ctrl%SuperP).OR.ICyc>1)THEN
               CALL Get(DMaxA,'DMax',Tag_O=TRIM(PrevTag))
               CALL Get(EtotA,'Etot',Tag_O=TRIM(PrevTag))
            ENDIF
            CALL Get(DMaxB,'DMax',Tag_O=TRIM(CrntTag))
            CALL Get(EtotB,'Etot',Tag_O=TRIM(CrntTag))
            IF(ICyc>0) &
            CALL Get(DIISB,'diiserr',Tag_O=TRIM(CrntTag))
            IF(ICyc>1) &
            CALL Get(DIISA,'diiserr',Tag_O=TRIM(PrevTag))
         ENDIF
!---------------------------------------------------------------------
!        Check for convergence
!
         Delta_DMax=ABS(DMaxA-DMaxB)
         Delta_ETot=ABS(ETotA-ETotB)
         ConvQ_ETot=ABS((ETotA-ETotB)/ETotB)
         ConvQ_DMax=ABS((DMaxA-DMaxB)/DMaxB)
         ConvQ_DIIS=ABS((DIISA-DIISB)/DIISB)
         WRITE(*,*)' ConvQ_DIIS = ',ConvQ_DIIS,' ConvQ_DMAX = ',ConvQ_DMAX
!        Could happen ...
         IF(ConvQ_ETot<1.D-14)THEN
            ConvergedQ=.TRUE.
            WRITE(*,*)' Met Convergence criteria A '
         ENDIF
!        Check to see if convergence is in an asymptotic regime
         IF(ConvQ_DIIS<1.D0.AND.ConvQ_DMax<1.D0.AND.ICyc>2)THEN                
!           Look for non-decreasing error stagnation due to incomplete numerics
            IF(DIISB>DIISA.OR.DMAXB>DMAXA)THEN
               ConvergedQ=.TRUE.
               WRITE(*,*)' Met Convergence criteria B: DIIS stalled.'
            ENDIF
         ENDIF
!        Look for convergence stall-outs 
         IF(ConvQ_DIIS<0.4D0.AND.ConvQ_DMAX<0.4D0)THEN
            WRITE(*,*)' Met Convergence criteria C: DIIS stalled.'
            ConvergedQ=.TRUE.
         ENDIF              
!        Check for absolute convergence below thresholds
         IF(ICyc>0.AND.Delta_DMax<DTol(Ctrl%AccL(IBas)).AND.ConvQ_ETot<ETol(Ctrl%AccL(IBas)))THEN 
            ConvergedQ=.TRUE.
            Ctrl%Fail=SCF_CONVERGED
            WRITE(*,*)' Met Convergence criteria D:  '
         ENDIF
!--------------------------------------------------------------------
!        IO 
!
         Mssg=' DIIS Err = '//TRIM(DblToMedmChar(DIISB)) &
           //', MAX/P = '//TRIM(DblToMedmChar(DMAXB))
         WRITE(*,*)TRIM(Mssg)
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
         Ctrl%Previous%I=Ctrl%Current%I
!
         RETURN
!
      END FUNCTION ConvergedQ
!
      FUNCTION StatsToChar(Stats) RESULT(StatString)
         TYPE(INT_VECT)  :: Stats
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: StatString
         StatString='_'//TRIM(IntToChar(Stats%I(3)))  &
                   //'_'//TRIM(IntToChar(Stats%I(2))) &
                   //'_'//TRIM(IntToChar(Stats%I(1)))
      END FUNCTION StatsToChar
!
      SUBROUTINE SCFSummry(Ctrl)
         TYPE(SCFControls) :: Ctrl
         TYPE(CRDS)        :: GM 
         INTEGER           :: I,K,IBas,ICyc,IGeo,NCyc,Mthd,BTyp
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
!------------------------------------------------------------------
!        Simplify notation
         ICyc=Ctrl%Current%I(1)
         IBas=Ctrl%Current%I(2)
         IGeo=Ctrl%Current%I(3)
         Mthd=Ctrl%Method(IBas)
         NCyc=Ctrl%NCyc(IBas)
         CALL OpenASCII(OutFile,Out)
         CALL PrintProtectL(Out)
         IF(PrintFlags%Key<=DEBUG_MINIMUM)THEN
            IF(IGeo==1.AND.IBas==1)THEN
               WRITE(Out,11)
               WRITE(Out,21)
               WRITE(Out,22)
               WRITE(Out,11)
            ENDIF
            WRITE(Out,23)IGeo,Ctrl%Stats(NCyc,1),NCyc+1,Ctrl%Stats(NCyc,3)
         ELSEIF(PrintFlags%Key==DEBUG_MEDIUM)THEN
            CALL OpenHDF(InfFile)
            CALL Get(GM,Tag_O=TRIM(IntToChar(IGeo)))
            CALL CloseHDF()
            CALL PPrint(GM)
            CALL Delete(GM)
         ELSEIF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
            WRITE(Out,11)
            CALL OpenHDF(InfFile)
            CALL Get(GM,Tag_O=TRIM(IntToChar(IGeo)))
            CALL CloseHDF()
            CALL PPrint(GM)
            CALL Delete(GM)
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
         ENDIF
         CALL PrintProtectR(Out)
         CLOSE(UNIT=Out,STATUS='KEEP')            
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
