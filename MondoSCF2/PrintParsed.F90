    ! Print the MondoSCF banner
    CALL PrintMondoSCF(OutFile)



!===============================================================================================
! PRINT THE MONDOSCF BANNER TO THE OUTPUT FILE
!===============================================================================================
  SUBROUTINE PrintMondoSCF(OFile)
    CHARACTER(LEN=DCL)  :: OFile,M_HOST,M_MACH,M_SYST,M_VRSN,M_PLAT
    INTEGER             :: I
!-----------------------------------------------------------------------------------------------         
    CALL TimeStamp('Starting MondoSCF')
    CALL GetEnv('MONDO_HOST',M_HOST)
    CALL GetEnv('MONDO_MACH',M_MACH)
    CALL GetEnv('MONDO_SYST',M_SYST)
    CALL GetEnv('MONDO_VRSN',M_VRSN)
    CALL GetEnv('MONDO_PLAT',M_PLAT)   
    CALL OpenASCII(OFile,Out)
    CALL PrintProtectL(Out)
    ! Print MondoSCF banner and authorship
    WRITE(Out,77)(Rtrn,I=1,15)
    WRITE(*,77)(Rtrn,I=1,15)
77  FORMAT(A1,A1,                                                   &
         ' __    __                 _       ____________ ______ ',A1,    &
         '|  \  /  |               | |     /       /    |      |',A1,    & 
         "|   \/   | ___  _ __   __| | ___/   ____/   __|  |==='",A1,    &
         "|        |/ _ \| '_ \ / _  |/ _ \____  \   (__|  ____|",A1,    &
         '|  |\/|  | (_) | | | | (_| | (_) )     /\     |  |    ',A1,    &
         '|__|  |__|\___/|_| |_|\____|\___/_____/  \____|__|    ',A1,A1, &
         ' Version 1.0 alpha 5                                  ',A1,    &  
         ' A program suite for O(N) SCF theory and ab initio MD ',A1,    &
         ' Matt Challacombe, Eric Schwegler, C.J. Tymczak,      ',A1,    &
         ' Chee Kwan Gan, Karoly Nemeth, Anders Niklasson,      ',A1,    &
         ' and Hugh Nymeyer                                     ',A1,    &
         ' Los Alamos National Laboratory (LA-CC 01-2)          ',A1,    & 
         ' Copyright 2001, University of California.            ',A1)
    ! Write information on host, platform, etc
    Mssg='Compiled for '//TRIM(M_PLAT)//', executing on '//TRIM(M_HOST)     &
         //Rtrn//' a '//TRIM(M_MACH)//' machine'//' running '//TRIM(M_SYST) &
         //' '//TRIM(M_VRSN)       
    WRITE(*,*)TRIM(Mssg)
    WRITE(Out,*)TRIM(Mssg)
    WRITE(*,*)
    WRITE(Out,*)
    CALL PrintProtectR(Out)
    CLOSE(Out)
  END SUBROUTINE PrintMondoSCF


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


