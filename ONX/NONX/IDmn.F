      FUNCTION IDmn(L)
c
c Number of locations needed to contract the VRRs
c
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER DimData(0:4)
      DATA DimData /1,4,17,24,36/
      IF(L.LE.4) THEN
         IDmn=DimData(L)
      ELSE
         CALL MondoHalt('L too big in IDmn')
      ENDIF
      RETURN
      END
