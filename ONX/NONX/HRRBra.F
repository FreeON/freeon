      SUBROUTINE HRRBra(W1,W2,x,y,z,N,LenBi,LenBf,NumK,IType)
c----------------------------------------------------------------------
c  Input:
c     W1     Contracted (p|Ket) integrals (from HRRKet step)
c     AC     Buffer of distances between bra basis functions
c     N      Number of integrals to compute
c     LenBi  Size of (Bra| space used in the VRR step
c     LenBf  Size of (Bra| space returned by the HRR step
c     NumK   Number of |Ket) locations to loop over
c     IType  Distribution type on the (Bra|
c
c  Output:
c     W2     (Bra|Ket) transformed integrals, the final 2-e ints
c
c  IType = TKet = KType*100 + LType
c
c----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER LenBi,LenBf,NumK
      INTEGER IType,N,I,J
      REAL*8 W1(N,LenBi,NumK)
      REAL*8 W2(N,LenBf,NumK)
      REAL*8 x,y,z

      IF(IType.EQ.0201) THEN
        DO 1000 J=1,NumK
        DO 1000 I=1,N
          W2(I,1,J)=W1(I,1,J)
          W2(I,2,J)=W1(I,2,J)
          W2(I,3,J)=W1(I,3,J)
          W2(I,4,J)=W1(I,4,J)
 1000   CONTINUE
        RETURN
      ELSEIF(IType.EQ.0202) THEN
        DO 1010 J=1,NumK
        DO 1010 I=1,N
          W2(I,1,J)=W1(I,1,J)
          W2(I,2,J)=W1(I,2,J)+x*W1(I,7,J)
          W2(I,3,J)=W1(I,3,J)+y*W1(I,7,J)
          W2(I,4,J)=W1(I,4,J)+z*W1(I,7,J)
          W2(I,5,J)=W1(I,5,J)
          W2(I,6,J)=W1(I,6,J)+x*W1(I,8,J)
          W2(I,7,J)=W1(I,10,J)+y*W1(I,8,J)
          W2(I,8,J)=W1(I,14,J)+z*W1(I,8,J)
          W2(I,9,J)=W1(I,9,J)
          W2(I,10,J)=W1(I,10,J)+x*W1(I,12,J)
          W2(I,11,J)=W1(I,11,J)+y*W1(I,12,J)
          W2(I,12,J)=W1(I,15,J)+z*W1(I,12,J)
          W2(I,13,J)=W1(I,13,J)
          W2(I,14,J)=W1(I,14,J)+x*W1(I,17,J)
          W2(I,15,J)=W1(I,15,J)+y*W1(I,17,J)
          W2(I,16,J)=W1(I,16,J)+z*W1(I,17,J)
 1010   CONTINUE
        RETURN
      ELSEIF(IType.EQ.0301) THEN
        DO 1100 J=1,NumK
        DO 1100 I=1,N
          W2(I,1,J)=W1(I,1,J)
          W2(I,2,J)=W1(I,2,J)
          W2(I,3,J)=W1(I,3,J)
 1100   CONTINUE
        RETURN
      ELSEIF(IType.EQ.0302) THEN
        DO 1110 J=1,NumK
        DO 1110 I=1,N
          W2(I,1,J)=W1(I,1,J)
          W2(I,2,J)=W1(I,6,J)+x*W1(I, 2,J)
          W2(I,3,J)=W1(I,7,J)+y*W1(I, 2,J)
          W2(I,4,J)=W1(I,10,J)+z*W1(I, 2,J)
          W2(I,5,J)=W1(I,5,J)
          W2(I,6,J)=W1(I,7,J)+x*W1(I, 3,J)
          W2(I,7,J)=W1(I,8,J)+y*W1(I, 3,J)
          W2(I,8,J)=W1(I,11,J)+z*W1(I, 3,J)
          W2(I,9,J)=W1(I,9,J)
          W2(I,10,J)=W1(I,10,J)+x*W1(I, 4,J)
          W2(I,11,J)=W1(I,11,J)+y*W1(I, 4,J)
          W2(I,12,J)=W1(I,12,J)+z*W1(I, 4,J)
 1110   CONTINUE
        RETURN
      ELSEIF(IType.EQ.0303) THEN
        DO 1111 J=1,NumK
        DO 1111 I=1,N
          W2(I,1,J)=W1(I,4,J)+x*W1(I,1,J)
          W2(I,2,J)=W1(I,5,J)+y*W1(I,1,J)
          W2(I,3,J)=W1(I,7,J)+z*W1(I,1,J)
          W2(I,4,J)=W1(I,5,J)+x*W1(I,2,J)
          W2(I,5,J)=W1(I,6,J)+y*W1(I,2,J)
          W2(I,6,J)=W1(I,8,J)+z*W1(I,2,J)
          W2(I,7,J)=W1(I,7,J)+x*W1(I,3,J)
          W2(I,8,J)=W1(I,8,J)+y*W1(I,3,J)
          W2(I,9,J)=W1(I,9,J)+z*W1(I,3,J)
 1111   CONTINUE
        RETURN
      ELSEIF(IType.EQ.0601) THEN
        DO 2200 J=1,NumK
        DO 2200 I=1,N
          W2(I,1,J)=W1(I,1,J)
          W2(I,2,J)=W1(I,2,J)
          W2(I,3,J)=W1(I,3,J)
          W2(I,4,J)=W1(I,4,J)
          W2(I,5,J)=W1(I,5,J)
          W2(I,6,J)=W1(I,6,J)
 2200   CONTINUE
        RETURN
      ELSEIF(IType.EQ.0602) THEN
        DO 2210 J=1,NumK
        DO 2210 I=1,N
          W2(I,1,J)=W1(I,1,J)
          W2(I,2,J)=W1(I,10,J)+x*W1(I,2,J)
          W2(I,3,J)=W1(I,11,J)+y*W1(I,2,J)
          W2(I,4,J)=W1(I,15,J)+z*W1(I,2,J)
          W2(I,5,J)=W1(I,5,J)
          W2(I,6,J)=W1(I,11,J)+x*W1(I,3,J)
          W2(I,7,J)=W1(I,12,J)+y*W1(I,3,J)
          W2(I,8,J)=W1(I,16,J)+z*W1(I,3,J)
          W2(I,9,J)=W1(I,9,J)
          W2(I,10,J)=W1(I,12,J)+x*W1(I,4,J)
          W2(I,11,J)=W1(I,14,J)+y*W1(I,4,J)
          W2(I,12,J)=W1(I,18,J)+z*W1(I,4,J)
          W2(I,13,J)=W1(I,13,J)
          W2(I,14,J)=W1(I,15,J)+x*W1(I,6,J)
          W2(I,15,J)=W1(I,16,J)+y*W1(I,6,J)
          W2(I,16,J)=W1(I,19,J)+z*W1(I,6,J)
          W2(I,17,J)=W1(I,17,J)
          W2(I,18,J)=W1(I,16,J)+x*W1(I,7,J)
          W2(I,19,J)=W1(I,18,J)+y*W1(I,7,J)
          W2(I,20,J)=W1(I,20,J)+z*W1(I,7,J)
          W2(I,21,J)=W1(I,21,J)
          W2(I,22,J)=W1(I,19,J)+x*W1(I,8,J)
          W2(I,23,J)=W1(I,20,J)+y*W1(I,8,J)
          W2(I,24,J)=W1(I,22,J)+z*W1(I,8,J)
 2210   CONTINUE
        RETURN
      ELSEIF(IType.EQ.0603) THEN
        DO 2211 J=1,NumK
        DO 2211 I=1,N
          W2(I,1,J)=W1(I,7,J)+x*W1(I,1,J)
          W2(I,2,J)=W1(I,8,J)+y*W1(I,1,J)
          W2(I,3,J)=W1(I,11,J)+z*W1(I,1,J)
          W2(I,4,J)=W1(I,8,J)+x*W1(I,2,J)
          W2(I,5,J)=W1(I,9,J)+y*W1(I,2,J)
          W2(I,6,J)=W1(I,12,J)+z*W1(I,2,J)
          W2(I,7,J)=W1(I,9,J)+x*W1(I,3,J)
          W2(I,8,J)=W1(I,10,J)+y*W1(I,3,J)
          W2(I,9,J)=W1(I,13,J)+z*W1(I,3,J)
          W2(I,10,J)=W1(I,11,J)+x*W1(I,4,J)
          W2(I,11,J)=W1(I,12,J)+y*W1(I,4,J)
          W2(I,12,J)=W1(I,14,J)+z*W1(I,4,J)
          W2(I,13,J)=W1(I,12,J)+x*W1(I,5,J)
          W2(I,14,J)=W1(I,13,J)+y*W1(I,5,J)
          W2(I,15,J)=W1(I,15,J)+z*W1(I,5,J)
          W2(I,16,J)=W1(I,14,J)+x*W1(I,6,J)
          W2(I,17,J)=W1(I,15,J)+y*W1(I,6,J)
          W2(I,18,J)=W1(I,16,J)+z*W1(I,6,J)
 2211   CONTINUE
        RETURN
      ELSEIF(IType.EQ.0606) THEN
        DO 2222 J=1,NumK
        DO 2222 I=1,N
        W2(I,1,J)=W1(I,17,J)+x*W1(I,7,J)+x*(W1(I,7,J)+x*W1(I,1,J))
        W2(I,2,J)=W1(I,18,J)+x*W1(I,8,J)+y*(W1(I,7,J)+x*W1(I,1,J))
        W2(I,3,J)=W1(I,19,J)+y*W1(I,8,J)+y*(W1(I,8,J)+y*W1(I,1,J))
        W2(I,4,J)=W1(I,22,J)+x*W1(I,11,J)+z*(W1(I,7,J)+x*W1(I,1,J))
        W2(I,5,J)=W1(I,23,J)+z*W1(I,8,J)+y*(W1(I,11,J)+z*W1(I,1,J))
        W2(I,6,J)=W1(I,26,J)+z*W1(I,11,J)+z*(W1(I,11,J)+z*W1(I,1,J))
        W2(I,7,J)=W1(I,18,J)+x*W1(I,8,J)+x*(W1(I,8,J)+x*W1(I,2,J))
        W2(I,8,J)=W1(I,19,J)+y*W1(I,8,J)+x*(W1(I,9,J)+y*W1(I,2,J))
        W2(I,9,J)=W1(I,20,J)+y*W1(I,9,J)+y*(W1(I,9,J)+y*W1(I,2,J))
        W2(I,10,J)=W1(I,23,J)+x*W1(I,12,J)+z*(W1(I,8,J)+x*W1(I,2,J))
        W2(I,11,J)=W1(I,24,J)+y*W1(I,12,J)+z*(W1(I,9,J)+y*W1(I,2,J))
        W2(I,12,J)=W1(I,27,J)+z*W1(I,12,J)+z*(W1(I,12,J)+z*W1(I,2,J))
        W2(I,13,J)=W1(I,19,J)+x*W1(I,9,J)+x*(W1(I,9,J)+x*W1(I,3,J))
        W2(I,14,J)=W1(I,20,J)+y*W1(I,9,J)+x*(W1(I,10,J)+y*W1(I,3,J))
        W2(I,15,J)=W1(I,21,J)+y*W1(I,10,J)+y*(W1(I,10,J)+y*W1(I,3,J))
        W2(I,16,J)=W1(I,24,J)+z*W1(I,9,J)+x*(W1(I,13,J)+z*W1(I,3,J))
        W2(I,17,J)=W1(I,25,J)+y*W1(I,13,J)+z*(W1(I,10,J)+y*W1(I,3,J))
        W2(I,18,J)=W1(I,28,J)+z*W1(I,13,J)+z*(W1(I,13,J)+z*W1(I,3,J))
        W2(I,19,J)=W1(I,22,J)+x*W1(I,11,J)+x*(W1(I,11,J)+x*W1(I,4,J))
        W2(I,20,J)=W1(I,23,J)+x*W1(I,12,J)+y*(W1(I,11,J)+x*W1(I,4,J))
        W2(I,21,J)=W1(I,24,J)+y*W1(I,12,J)+y*(W1(I,12,J)+y*W1(I,4,J))
        W2(I,22,J)=W1(I,26,J)+z*W1(I,11,J)+x*(W1(I,14,J)+z*W1(I,4,J))
        W2(I,23,J)=W1(I,27,J)+z*W1(I,12,J)+y*(W1(I,14,J)+z*W1(I,4,J))
        W2(I,24,J)=W1(I,29,J)+z*W1(I,14,J)+z*(W1(I,14,J)+z*W1(I,4,J))
        W2(I,25,J)=W1(I,23,J)+x*W1(I,12,J)+x*(W1(I,12,J)+x*W1(I,5,J))
        W2(I,26,J)=W1(I,24,J)+y*W1(I,12,J)+x*(W1(I,13,J)+y*W1(I,5,J))
        W2(I,27,J)=W1(I,25,J)+y*W1(I,13,J)+y*(W1(I,13,J)+y*W1(I,5,J))
        W2(I,28,J)=W1(I,27,J)+z*W1(I,12,J)+x*(W1(I,15,J)+z*W1(I,5,J))
        W2(I,29,J)=W1(I,28,J)+z*W1(I,13,J)+y*(W1(I,15,J)+z*W1(I,5,J))
        W2(I,30,J)=W1(I,30,J)+z*W1(I,15,J)+z*(W1(I,15,J)+z*W1(I,5,J))
        W2(I,31,J)=W1(I,26,J)+x*W1(I,14,J)+x*(W1(I,14,J)+x*W1(I,6,J))
        W2(I,32,J)=W1(I,27,J)+y*W1(I,14,J)+x*(W1(I,15,J)+y*W1(I,6,J))
        W2(I,33,J)=W1(I,28,J)+y*W1(I,15,J)+y*(W1(I,15,J)+y*W1(I,6,J))
        W2(I,34,J)=W1(I,29,J)+z*W1(I,14,J)+x*(W1(I,16,J)+z*W1(I,6,J))
        W2(I,35,J)=W1(I,30,J)+z*W1(I,15,J)+y*(W1(I,16,J)+z*W1(I,6,J))
        W2(I,36,J)=W1(I,31,J)+z*W1(I,16,J)+z*(W1(I,16,J)+z*W1(I,6,J))
 2222   CONTINUE
        RETURN
c
c Codes for f-functions:
c
c     ELSEIF(IType.EQ.1001) THEN
c     ELSEIF(IType.EQ.1002) THEN
c     ELSEIF(IType.EQ.1003) THEN
c     ELSEIF(IType.EQ.1006) THEN
c     ELSEIF(IType.EQ.1010) THEN
      ELSE
         WRITE(*,*) "IType = ",IType
         STOP 'Illegal IType in hrrbra_old'
      ENDIF
      RETURN
      END
