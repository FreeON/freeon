      SUBROUTINE HRRBra(W1,W2,AC,N,Ind,LenBi,LenBf,NumK,IType)
c---------------------------------------------------------------------- 
c  Input:
c     W1     Contracted (p|q) integrals (from VRR+contraction steps)
c     AC     Buffer of distances between bra basis functions
c     N      Number of integrals to compute
c     Ind    Buffer of indexes to get the right AC info
c     LenBi  Size of (Bra| space used in the VRR step
c     LenBf  Size of (Bra| space returned by the HRR step
c     NumK   Number of |Ket) locations to loop over
c     IType  Distribution type on the (Bra|
c
c  Output:
c     W2     transformed integrals
c
c  IType = LAmax*1000+LAmin*100+LCmax*10+LCmin
c 
c---------------------------------------------------------------------- 
      IMPLICIT NONE
      INTEGER LenBi,LenBf,NumK
      INTEGER IType,N,I,J
      INTEGER Ind(2,N)
      REAL*8 W1(N,LenBi,NumK)
      REAL*8 W2(N,LenBf,NumK)
      REAL*8 AC(3,*),x,y,z

      IF(IType.EQ.1000) THEN
        DO 1000 J=1,NumK
        DO 1000 I=1,N
          W2(I,1,J)=W1(I,1,J)
          W2(I,2,J)=W1(I,2,J)
          W2(I,3,J)=W1(I,3,J)
          W2(I,4,J)=W1(I,4,J)
 1000   CONTINUE
        RETURN
      ELSEIF(IType.EQ.1010) THEN
        DO 1010 J=1,NumK
        DO 1010 I=1,N
          x=AC(1,Ind(1,I))
          y=AC(2,Ind(1,I))
          z=AC(3,Ind(1,I))
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
      ELSEIF(IType.EQ.1100) THEN
        DO 1100 J=1,NumK
        DO 1100 I=1,N
          W2(I,1,J)=W1(I,2,J)
          W2(I,2,J)=W1(I,3,J)
          W2(I,3,J)=W1(I,4,J)
 1100   CONTINUE
        RETURN
      ELSEIF(IType.EQ.1110) THEN
        DO 1110 J=1,NumK
        DO 1110 I=1,N
          x=AC(1,Ind(1,I))
          y=AC(2,Ind(1,I))
          z=AC(3,Ind(1,I))
          W2(I,1,J)=W1(I,5,J)
          W2(I,2,J)=W1(I,6,J)+x*W1(I,8,J)
          W2(I,3,J)=W1(I,10,J)+y*W1(I,8,J)
          W2(I,4,J)=W1(I,14,J)+z*W1(I,8,J)
          W2(I,5,J)=W1(I,9,J)
          W2(I,6,J)=W1(I,10,J)+x*W1(I,12,J)
          W2(I,7,J)=W1(I,11,J)+y*W1(I,12,J)
          W2(I,8,J)=W1(I,15,J)+z*W1(I,12,J)
          W2(I,9,J)=W1(I,13,J)
          W2(I,10,J)=W1(I,14,J)+x*W1(I,17,J)
          W2(I,11,J)=W1(I,15,J)+y*W1(I,17,J)
          W2(I,12,J)=W1(I,16,J)+z*W1(I,17,J)
 1110   CONTINUE
        RETURN
      ELSEIF(IType.EQ.1111) THEN
        DO 1111 J=1,NumK
        DO 1111 I=1,N
          x=AC(1,Ind(1,I))
          y=AC(2,Ind(1,I))
          z=AC(3,Ind(1,I))
          W2(I,1,J)=W1(I,6,J)+x*W1(I,8,J)
          W2(I,2,J)=W1(I,10,J)+y*W1(I,8,J)
          W2(I,3,J)=W1(I,14,J)+z*W1(I,8,J)
          W2(I,4,J)=W1(I,10,J)+x*W1(I,12,J)
          W2(I,5,J)=W1(I,11,J)+y*W1(I,12,J)
          W2(I,6,J)=W1(I,15,J)+z*W1(I,12,J)
          W2(I,7,J)=W1(I,14,J)+x*W1(I,17,J)
          W2(I,8,J)=W1(I,15,J)+y*W1(I,17,J)
          W2(I,9,J)=W1(I,16,J)+z*W1(I,17,J)
 1111   CONTINUE
        RETURN
      ELSEIF(IType.EQ.2200) THEN
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
      ELSEIF(IType.EQ.2210) THEN
        DO 2210 J=1,NumK
        DO 2210 I=1,N
          x=AC(1,Ind(1,I))
          y=AC(2,Ind(1,I))
          z=AC(3,Ind(1,I))
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
      ELSEIF(IType.EQ.2211) THEN
        DO 2211 J=1,NumK
        DO 2211 I=1,N
          x=AC(1,Ind(1,I))
          y=AC(2,Ind(1,I))
          z=AC(3,Ind(1,I))
          W2(I,1,J)=W1(I,10,J)+x*W1(I,2,J)
          W2(I,2,J)=W1(I,11,J)+y*W1(I,2,J)
          W2(I,3,J)=W1(I,15,J)+z*W1(I,2,J)
          W2(I,4,J)=W1(I,11,J)+x*W1(I,3,J)
          W2(I,5,J)=W1(I,12,J)+y*W1(I,3,J)
          W2(I,6,J)=W1(I,16,J)+z*W1(I,3,J)
          W2(I,7,J)=W1(I,12,J)+x*W1(I,4,J)
          W2(I,8,J)=W1(I,14,J)+y*W1(I,4,J)
          W2(I,9,J)=W1(I,18,J)+z*W1(I,4,J)
          W2(I,10,J)=W1(I,15,J)+x*W1(I,6,J)
          W2(I,11,J)=W1(I,16,J)+y*W1(I,6,J)
          W2(I,12,J)=W1(I,19,J)+z*W1(I,6,J)
          W2(I,13,J)=W1(I,16,J)+x*W1(I,7,J)
          W2(I,14,J)=W1(I,18,J)+y*W1(I,7,J)
          W2(I,15,J)=W1(I,20,J)+z*W1(I,7,J)
          W2(I,16,J)=W1(I,19,J)+x*W1(I,8,J)
          W2(I,17,J)=W1(I,20,J)+y*W1(I,8,J)
          W2(I,18,J)=W1(I,22,J)+z*W1(I,8,J)
 2211   CONTINUE
        RETURN
      ELSEIF(IType.EQ.2222) THEN
        DO 2222 J=1,NumK
        DO 2222 I=1,N
        x=AC(1,Ind(1,I))
        y=AC(2,Ind(1,I))
        z=AC(3,Ind(1,I))
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
        W2(I,22,J)=W1(I,26,J)+z*W1(I,11,J)+x*(W1(I,12,J)+z*W1(I,4,J))
        W2(I,23,J)=W1(I,27,J)+z*W1(I,12,J)+y*(W1(I,12,J)+z*W1(I,4,J))
        W2(I,24,J)=W1(I,29,J)+z*W1(I,14,J)+z*(W1(I,12,J)+z*W1(I,4,J))
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
c     ELSEIF(IType.EQ.3300) THEN
c     ELSEIF(IType.EQ.3310) THEN
c     ELSEIF(IType.EQ.3311) THEN
c     ELSEIF(IType.EQ.3322) THEN
c     ELSEIF(IType.EQ.3333) THEN
      ELSE
         WRITE(*,*) "IType = ",IType
         CALL MondoHalt('Illegal IType in BraHRR')
      ENDIF
      RETURN
      END
