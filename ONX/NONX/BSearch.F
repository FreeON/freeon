      SUBROUTINE BSearch(N,Lng,Target,Buff)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      REAL*8 Buff(*),Target
      INTEGER K,L,MID,N,Lng
      IF(Target.GT.Buff(1)) THEN
        N=0
      ELSEIF(Target.LT.Buff(Lng)) THEN
        N=Lng
      ELSE
        K=1
        L=Lng
        MID=Lng/2
 1000   CONTINUE
        IF(Target.LE.Buff(MID)) THEN
          K=MID
        ELSE
          L=MID
        ENDIF
        MID=(K+L)/2
        IF((K.LT.MID).AND.(MID.LT.L)) GOTO 1000
        N=MID
      ENDIF
      RETURN
      END
