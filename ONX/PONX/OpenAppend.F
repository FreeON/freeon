C
C     Matt Challacombe
C     Los Alamos National Laboratory
C     Copywrite 1997, The University of California
C
      SUBROUTINE OpenAppend(iUnit,FileName)
C
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      CHARACTER*(*) FileName
      CHARACTER*128 File2
      CHARACTER*72  Dummy
      CHARACTER*4   IOErr
C    
      LOGICAL Opened,Present
C
      File2=FileName
      jf=LineEnd(File2)
      INQUIRE(FILE=File2(1:jf),OPENED=Opened,EXIST=Present,
     >        ERR=101,IOSTAT=IOS)
C
      IF(Opened)THEN
         WRITE(6,*)' Already opened:',FileName
         STOP
      ELSEIF(Present)THEN
         OPEN(UNIT=iUnit,FILE=FileName,ACCESS='SEQUENTIAL',
     >        FORM='FORMATTED',ERR=101,IOSTAT=IOS,STATUS='OLD')
      ELSEIF(.NOT.Present)THEN
         OPEN(UNIT=iUnit,FILE=FileName,ACCESS='SEQUENTIAL',
     >        FORM='FORMATTED',ERR=101,IOSTAT=IOS,STATUS='NEW')         
      ENDIF

      DO WHILE(2+2.EQ.4)
         READ(iUnit,100,END=99)Dummy
      ENDDO

 99   BACKSPACE(UNIT=iUnit)
C
      RETURN
C
 100  FORMAT(A72)
 101  WRITE(IOErr,'(I4)')IOS
      CALL MondoHalt('OpenAppend: '//File2(1:jf)//' IOErr = '//IOErr)
      END
