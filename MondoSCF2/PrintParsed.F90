MODULE PrintParsed
  USE InOut
#ifdef NAG
  USE F90_UNIX
#endif
  USE PrettyPrint
  USE ControlStructures
  IMPLICIT NONE
CONTAINS
  !===============================================================================================
  ! PRINT THE MONDOSCF BANNER TO THE OUTPUT FILE
  !===============================================================================================
  SUBROUTINE PrintsStartUp(N)
    TYPE(FileNames) :: N
    CHARACTER(LEN=DCL)  :: OFile,M_HOST,M_MACH,M_SYST,M_VRSN,M_PLAT
    INTEGER             :: I
    !-----------------------------------------------------------------------------------------------         
    CALL GetEnv('MONDO_HOST',M_HOST)
    CALL GetEnv('MONDO_MACH',M_MACH)
    CALL GetEnv('MONDO_SYST',M_SYST)
    CALL GetEnv('MONDO_VRSN',M_VRSN)
    CALL GetEnv('MONDO_PLAT',M_PLAT)   
    WRITE(*,*)' OUTFILE = ',N%OFile
    CALL OpenASCII(N%OFile,Out)
    CALL PrintProtectL(Out)
    ! Print MondoSCF banner and authorship
    WRITE(Out,77)(Rtrn,I=1,15)
77  FORMAT(A1,A1,                                                        &
         ' __    __                 _       ____________ ______ ',A1,    &
         '|  \  /  |               | |     /       /    |      |',A1,    & 
         "|   \/   | ___  _ __   __| | ___/   ____/   __|  |==='",A1,    &
         "|        |/ _ \| '_ \ / _  |/ _ \____  \   (__|  ____|",A1,    &
         '|  |\/|  | (_) | | | | (_| | (_) )     /\     |  |    ',A1,    &
         '|__|  |__|\___/|_| |_|\____|\___/_____/  \____|__|    ',A1,A1, &
         ' Version 1.0 alpha 8                                  ',A1,    &  
         ' A program suite for O(N) SCF theory and ab initio MD ',A1,    &
         ' Matt Challacombe, Eric Schwegler, C.J. Tymczak,      ',A1,    &
         ' Chee Kwan Gan, Karoly Nemeth, Anders Niklasson,      ',A1,    &
         ' Graeme Henkelman, and Valery Weber                   ',A1,    &
         ' Los Alamos National Laboratory (LA-CC 01-2)          ',A1,    & 
         ' Copyright 2001, University of California.            ',A1)
    ! Write information on host, platform, etc
    Mssg='Compiled for '//TRIM(M_PLAT)//', executing on '//TRIM(M_HOST)     &
         //Rtrn//' a '//TRIM(M_MACH)//' machine'//' running '//TRIM(M_SYST) &
         //' '//TRIM(M_VRSN)       
    WRITE(Out,*)TRIM(Mssg)
    WRITE(Out,*)    
    WRITE(Out,*)'Current HDF :: ',TRIM(N%HFile)
    CALL PrintProtectR(Out)
    CLOSE(Out)
    CALL TimeStamp('Starting MondoSCF')
  END SUBROUTINE PrintsStartUp
END MODULE PrintParsed

#ifdef DLKFJLSDJFLSDJFLSFDJ





    ! Print the MondoSCF banner
    CALL PrintMondoSCF(OutFile)





    ! Debug if asked
    IF(PrintFlags%Key==DEBUG_MAXIMUM)THEN
       CALL OpenASCII(OutFile,Out)           
       CALL PrintProtectL(Out)

#if defined(PARALLEL) && !defined(MPI2)
     Mssg='MPI Invokation : '//TRIM(MPI_INVOKE)
       WRITE(Out,*)TRIM(Mssg)
#endif
       CALL New(CBeg,NPrc-1,0)
       CALL New(CEnd,NPrc-1,0)
       DO I=0,NPrc-1
          CBeg%C(I)=IntToChar(Beg%I(I))
          CEnd%C(I)=IntToChar(End%I(I))
       ENDDO
       WRITE(Out,*)'Atomic partitioning = ',                                 &
            ('['//TRIM(CBeg%C(K))//'-'//TRIM(CEnd%C(K))//'], ',K=0,NPrc-2),  &
            '['//TRIM(CBeg%C(NPrc-1))//'-'//TRIM(CEnd%C(NPrc-1)),']'

       DO K=0,NPrc-1
          CBeg%C(K)=IntToChar(SUM(Bsiz%I(Beg%I(K):End%I(K))))
       ENDDO
       WRITE(Out,*)'Basis Functions per node = ',                &
            (TRIM(CBeg%C(K))//', ',K=0,NPrc-2),TRIM(CBeg%C(NPrc-1))
       CALL PrintProtectR(Out)
       CLOSE(UNIT=Out,STATUS='KEEP')
       CALL Delete(CBeg)
       CALL Delete(CEnd)
    ENDIF
    !-------------------------------------------------------
    !        Tidy up

!============================================================================
!     Print Out the Parsed Information
!============================================================================
      SUBROUTINE ParsePrint(Ctrl)
        TYPE(SCFControls)                :: Ctrl
        INTEGER                          :: I,RestAccL,RestMeth,RestModl
        CHARACTER(LEN=8)                 :: Cur
        CHARACTER(LEN=DEFAULT_CHR_LEN)   :: Method,Accuracy,Chemistry
        CHARACTER(LEN=BASESET_CHR_LEN)   :: BName   
        CHARACTER(LEN=2*DEFAULT_CHR_LEN) :: Mssg
!----------------------------------------------------------------------------
        CALL PrintProtectL(Out)
        CALL OpenASCII(OutFile,Out)
        WRITE(Out,*)'Current HDF file :: ',TRIM(Ctrl%Info)
        IF(Ctrl%Rest)THEN
           WRITE(Out,*)'Restart HDF file :: ',TRIM(Ctrl%OldInfo)
           CALL OpenHDF(Ctrl%OldInfo)
           Cur=IntToChar(Ctrl%Previous(2))
           CALL Get(RestAccL,'SCFAccuracy',Cur)
           CALL Get(RestMeth,'SCFMethod',Cur)
           CALL Get(RestModl,'ModelChemistry',Cur)
           CALL Get(BName,'bsetname',Cur)
           CALL CloseHDF()
           CALL OpenHDF(InfFile)
           Cur=IntToChar(Ctrl%Previous(3)) 
           Mssg='Restart using '//TRIM(BName)//'/'//TRIM(FunctionalName(RestModl)) &
                //' density and coordinates from previous geometry #'//TRIM(Cur)
           WRITE(Out,*)TRIM(Mssg)
        ENDIF
        WRITE(Out,*)
        ! Print the accuracy, method and model chemistry for each basis set
        WRITE(Out,*)'PROGRAM OF CALCULATIONS:',Rtrn
        DO I=1,Ctrl%NSet
!
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
              IF(Ctrl%Method(I)==SDMM_R_SCF)THEN
                 Method='restricted simplified density matrix minimization using '
              ELSEIF(Ctrl%Method(I)==PM_R_SCF)THEN
                 Method='restricted Palser-Manolopoulos (PM) purification using '
              ELSEIF(Ctrl%Method(I)==SP2_R_SCF)THEN
                 Method='restricted Quadratic Spectral Projection (SP2) purification using '
              ELSEIF(Ctrl%Method(I)==SP4_R_SCF)THEN
                 Method='restricted Quartic Spectral Projection (SP4) purification using '
              ELSEIF(Ctrl%Method(I)==TS4_R_SCF)THEN
                 Method='restricted Quartic Trace Setting (TS4) purification using '
              ELSE
                 Method='restricted Roothaan-Hall solution of the SCF using '
              ENDIF
!
              Chemistry=' and the '//TRIM(FunctionalName(Ctrl%Model(I)))//' model.'
              IF(I==1)THEN
                 Mssg=' A '//TRIM(Method)//Rtrn//TRIM(Accuracy)//' '//TRIM(Ctrl%BName(1))//TRIM(Chemistry)
              ELSE
                 Mssg=' Followed by a  '//TRIM(Method)//Rtrn//TRIM(Accuracy) &
                      //TRIM(Ctrl%BName(I))//TRIM(Chemistry)
              ENDIF
              WRITE(Out,*)TRIM(Mssg),Rtrn
!
! Put SCF Method
!
              CALL Put(Ctrl%AccL(I),'SCFAccuracy',Tag_O=IntToChar(I))
              CALL Put(Ctrl%Method(I),'SCFMethod',Tag_O=IntToChar(I))
!
#ifdef MMech
        ENDIF
#endif
        ENDDO
        CALL PrintProtectR(Out)
        CLOSE(Out,STATUS='KEEP')
!
      END SUBROUTINE ParsePrint


  !===============================================================================

  !===============================================================================
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

  SUBROUTINE SCFSummry(Ctrl)
    TYPE(SCFControls) :: Ctrl
    TYPE(CRDS)        :: GM 
    INTEGER           :: I,K,IBas,ICyc,IGeo,NCyc,Mthd,BTyp
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
    !-------------------------------------------------------------------------------
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
10  FORMAT(72('-'))
11  FORMAT(72('='))

21  FORMAT(' Ab initio ')
22  FORMAT('  Config         <E_SCF>       Cycles      DeltaD ')
23  FORMAT(1x,I5,2x,F20.8,5x,I2,7x,D10.4)

40  FORMAT('SCF  ',2x,'Ink',7x,'Total ',12x,'Delta',7x,'Max    ',5x,'DIIS   ')
41  FORMAT('Cycle',2x,'Fok',7x,'Energy',12x,'E_Tot',7x,'Delta P',5x,'Error  ')
42  FORMAT(1x,I2,',   ',L2,', ',F20.8,4(', ',D10.4))

  END SUBROUTINE SCFSummry
#endif
