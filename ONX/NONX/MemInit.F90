SUBROUTINE MemInit
  USE InOut
  USE ONXMemory
  IMPLICIT NONE

  IF (ErrorCode==eInit) THEN
    MXINT = 1000
    MAXD  = 100
    MAXC  = 10
    MAXT  = 10
    MAXP  = 8
  ELSEIF (ErrorCode==eMXINT) THEN
    MXINT = MXINT * 5
  ELSEIF (ErrorCode==eMAXD) THEN
    MAXD = MAXD * 2
  ELSEIF (ErrorCode==eMAXC) THEN
    MAXC = MAXC + 1
  ELSEIF (ErrorCode==eMAXT) THEN
    MAXT = MAXT + 1
  ELSE
    write(*,*) "Error code = ",ErrorCode
    CALL Halt(' Unknown error code in ONX')
  END IF

END SUBROUTINE MemInit
