      SUBROUTINE SkipOut(LenA,LenB,SchA,SchB,IndA,Thresh)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      REAL*8 Thresh,SchA(*),SchB(*)
      INTEGER LenA,LenB,IndA(0:500),I
      CALL BSearch(IndA(0),LenA,Thresh/SchB(1),SchA)
      CALL BSearch(IndA(1),LenB,Thresh/SchA(1),SchB)
      DO I=2,IndA(0)
         CALL LSearch(IndA(I),IndA(I-1),Thresh/SchA(I),SchB)
      ENDDO
      RETURN
      END
