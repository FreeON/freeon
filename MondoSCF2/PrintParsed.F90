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
         ' Version 1.0 alpha 10 "Gyokuro"                       ',A1,    &  
         ' A program suite for O(N) SCF theory and ab initio MD ',A1,    &
         ' Matt Challacombe, C.J. Tymczak,                      ',A1,    &
         ' Chee Kwan Gan, Karoly Nemeth, Valery Weber,          ',A1,    &  
         ' Eric Schwegler and Anders Niklasson                  ',A1,    &
         ' Los Alamos National Laboratory                       ',A1,    & 
         ' LA-CC-04-086 (formerly 01-2)                         ',A1,    & 
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
