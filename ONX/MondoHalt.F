C
C     Matt Challacombe
C     Los Alamos National Laboratory
C     Copywrite 1997, The University of California
C
      SUBROUTINE MondoHalt(Message)
C
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
C
      CHARACTER*(*) Message
C
      INCLUDE 'Units.Inc'
      INCLUDE 'Characters.Inc'
C  
      LOGICAL Op
C
      WRITE(iErr,*)Message
C
      CALL OpenAppend(iOut,OutFile)
      WRITE(iOut,*)Message
      CLOSE(UNIT=iOut,STATUS='KEEP')
C
      CALL OpenAppend(iLog,LogFile)
      WRITE(iLog,*)'DEAD'
      CLOSE(UNIT=iLog,STATUS='KEEP')
C
      CALL Trap()
C
      STOP
      END
C
