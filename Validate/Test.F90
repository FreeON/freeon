PROGRAM TestMondoSCFOutput
   USE DerivedTypes
   USE GlobalScalars
   USE GlobalObjects
   USE GlobalCharacters
   USE ProcessControl
   USE InOut
   USE MemMan
   USE ParsingConstants
   USE Parse
   USE PrettyPrint
!
   TYPE(ARGMT) :: Arg
   INTEGER                        :: I,I1,I2
   INTEGER,PARAMETER              :: Chk=22
   REAL(DOUBLE)                   :: D1,D2
   CHARACTER(LEN=DEFAULT_CHR_LEN) :: PWD,FileToCheck,RefFile, &
                                     Line1,Line2
!------------------------------------------------------------
   CALL Get(Arg)
   CALL GetEnv('PWD',PWD)   
!
   FileToCheck=TRIM(PWD)//'/'//TRIM(Arg%C%C(1))
   RefFile=TRIM(PWD)//'/'//TRIM(Arg%C%C(2))
   LogFile=TRIM(PWD)//'/Validation.log'
!
   CALL OpenASCII(FileToCheck,Inp,ReWind_O=.TRUE.)
   CALL OpenASCII(RefFile,Chk,ReWind_O=.TRUE.)
!
   DO I=1,100000
      READ(Inp,DEFAULT_CHR_FMT,END=1)Line1
      READ(Chk,DEFAULT_CHR_FMT,END=1)Line2
!---------------------------------------------------------------
!     Test total energies
!
      I1=0;I2=0
      IF(INDEX(Line1,'<SCF>')/=0)I1=INDEX(Line1,'=')
      IF(INDEX(Line2,'<SCF>')/=0)I2=INDEX(Line2,'=')
      IF(I1/=0.AND.I2/=0)THEN
         D1=CharToDbl(TRIM(Line1(I1+1:)))
         D2=CharToDbl(TRIM(Line2(I2+1:)))
         RD=ABS((D1-D2)/D2)
         IF(RD>1.D-10)THEN
            WRITE(*,*)' VALIDATION TEST FAILED! '
            CALL OpenASCII(LogFile,Log)
            WRITE(Log,*)'Test failed on comparison of total energies: '
            WRITE(Log,*)TRIM(Line1) 
            WRITE(Log,*)TRIM(Line2)
            CLOSE(Log)
            STOP
         ENDIF
      ELSEIF(I1/=0.AND.I2==0.OR.I1==0.AND.I2/=0)THEN
        WRITE(*,*)' VALIDATION TEST FAILED! '
        CALL OpenASCII(LogFile,Log)
        WRITE(Log,*)'Test failed by lines out of sinc '
        WRITE(Log,*)TRIM(Line1) 
        WRITE(Log,*)TRIM(Line2)
        CLOSE(Log)
        STOP
      ENDIF
!---------------------------------------------------------------
!     Test check sums
!
      I1=0;I2=0
      IF(INDEX(Line1,'CheckSum')/=0)I1=INDEX(Line1,'=')
      IF(INDEX(Line2,'CheckSum')/=0)I2=INDEX(Line2,'=')
      IF(I1/=0.AND.I2/=0)THEN
         D1=CharToDbl(TRIM(Line1(I1+1:)))
         D2=CharToDbl(TRIM(Line2(I2+1:)))
         RD=ABS((D1-D2)/D2)
         IF(RD>1.D-8)THEN
            WRITE(*,*)' VALIDATION TEST FAILED! '
            CALL OpenASCII(LogFile,Log)
            WRITE(Log,*)'Test failed on comparison of check sums: '
            WRITE(Log,*)TRIM(Line1) 
            WRITE(Log,*)TRIM(Line2)
            CLOSE(Log)
            STOP
         ENDIF
      ELSEIF(I1/=0.AND.I2==0.OR.I1==0.AND.I2/=0)THEN
        WRITE(*,*)' VALIDATION TEST FAILED! '
        CALL OpenASCII(LogFile,Log)
        WRITE(Log,*)'Test failed by lines out of sinc '
        WRITE(Log,*)TRIM(Line1) 
        WRITE(Log,*)TRIM(Line2)
        CLOSE(Log)
        STOP
      ENDIF
   ENDDO
 1 CALL OpenASCII(LogFile,Log)
   WRITE(Log,*)'Test of '//TRIM(Arg%C%C(1))//' succesful.'
   CLOSE(Log)
END PROGRAM TestMondoSCFOutput
