      SUBROUTINE LSearch(N,Lng,Target,Buff)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      REAL*8 Buff(*),Target
      INTEGER N,Lng,I
      IF(Target.LT.Buff(Lng)) THEN
         N=Lng
         RETURN
      ELSEIF(Target.GT.Buff(1)) THEN
         N=0
         RETURN
      ELSE
         N=Lng+1 
 1000    CONTINUE
         N=N-1
         IF(Target.GT.Buff(N)) GOTO 1000
      ENDIF
      RETURN
      END
