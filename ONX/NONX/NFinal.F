      FUNCTION NFinal(IType)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INTEGER NFinal,IType
      INTEGER T2Size(28)
      SAVE T2Size
      DATA T2Size /  1,
     >               4, 3,
     >              10, 9, 6,
     >              20,19,16,10,
     >              35,34,31,25,15,
     >              56,55,52,46,36,21,
     >              84,83,80,74,64,49,28 /
      IF(IType.LE.28) THEN
         NFinal=T2Size(IType)
      ELSE
         WRITE(*,*) "IType=",IType
         CALL MondoHalt('Illegal IType in NFinal')
      ENDIF
      RETURN
      END
