      SUBROUTINE HRRKetOld(W,BD,N,Ind,LenB,NumB,IType)
c----------------------------------------------------------------------
c  Input:
c     W      Contracted (Bra|q) integrals (from VRR+contraction steps)
c     BD     Buffer of distances between ket basis functions
c     N      Number of integrals to compute
c     Ind    Buffer of indexes to get the right BD info
c     LenB   Size of the bra space present
c     NumB   Number of (bra| locations to loop over
c     IType  Distribution type on the |Ket)
c
c  Output:
c     W      (Bra|Ket) transformed integrals
c
c  IType = LBmax*1000+LBmin*100+LDmax*10+LDmin
c
c----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N,LenB,NumB,IType,NN
      INTEGER I,I0
      INTEGER Ind(2,N)
      REAL*8 W(N*LenB,*)
      REAL*8 T(36)
      REAL*8 BD(3,*),x,y,z

      NN=N*NumB
      IF(    IType.EQ.1000) THEN
         RETURN
      ELSEIF(IType.EQ.1010) THEN
        DO 1010 I=1,NN
          I0=MOD(I-1,N)+1
          x=BD(1,Ind(2,I0))
          y=BD(2,Ind(2,I0))
          z=BD(3,Ind(2,I0))
          T(1)=W(I,1)
          T(2)=W(I,2)+x*W(I,7)
          T(3)=W(I,3)+y*W(I,7)
          T(4)=W(I,4)+z*W(I,7)
          T(5)=W(I,5)
          T(6)=W(I,6)+x*W(I, 8)
          T(7)=W(I,10)+y*W(I,8)
          T(8)=W(I,14)+z*W(I,8)
          T(9)=W(I,9)
          T(10)=W(I,10)+x*W(I,12)
          T(11)=W(I,11)+y*W(I,12)
          T(12)=W(I,15)+z*W(I,12)
          T(13)=W(I,13)
          T(14)=W(I,14)+x*W(I,17)
          T(15)=W(I,15)+y*W(I,17)
          T(16)=W(I,16)+z*W(I,17)
          W(I, 1) = T(1)
          W(I, 2) = T(2)
          W(I, 3) = T(3)
          W(I, 4) = T(4)
          W(I, 5) = T(5)
          W(I, 6) = T(6)
          W(I, 7) = T(7)
          W(I, 8) = T(8)
          W(I, 9) = T(9)
          W(I,10) = T(10)
          W(I,11) = T(11)
          W(I,12) = T(12)
          W(I,13) = T(13)
          W(I,14) = T(14)
          W(I,15) = T(15)
          W(I,16) = T(16)
 1010   CONTINUE
        RETURN
      ELSEIF(IType.EQ.1100) THEN
         RETURN
      ELSEIF(IType.EQ.1110) THEN
        DO 1110 I=1,NN
          I0=MOD(I-1,N)+1
          x=BD(1,Ind(2,I0))
          y=BD(2,Ind(2,I0))
          z=BD(3,Ind(2,I0))
          T(2)=W(I,6)+x*W(I,2)
          T(3)=W(I,7)+y*W(I,2)
          T(4)=W(I,10)+z*W(I,2)
          T(6)=W(I,7)+x*W(I,3)
          T(7)=W(I,8)+y*W(I,3)
          T(8)=W(I,11)+z*W(I,3)
          T(10)=W(I,10)+x*W(I,4)
          T(11)=W(I,11)+y*W(I,4)
          T(12)=W(I,12)+z*W(I,4)
          W(I,2)=T(2)
          W(I,3)=T(3)
          W(I,4)=T(4)
          W(I,6)=T(6)
          W(I,7)=T(7)
          W(I,8)=T(8)
          W(I,10)=T(10)
          W(I,11)=T(11)
          W(I,12)=T(12)
 1110   CONTINUE
        RETURN
      ELSEIF(IType.EQ.1111) THEN
        DO 1111 I=1,NN
          I0 = MOD(I-1,N)+1
          x=BD(1,Ind(2,I0))
          y=BD(2,Ind(2,I0))
          z=BD(3,Ind(2,I0))
          T(1)=W(I,4)+x*W(I,1)
          T(2)=W(I,5)+y*W(I,1)
          T(3)=W(I,7)+z*W(I,1)
          T(4)=W(I,5)+x*W(I,2)
          T(5)=W(I,6)+y*W(I,2)
          T(6)=W(I,8)+z*W(I,2)
          T(7)=W(I,7)+x*W(I,3)
          T(8)=W(I,8)+y*W(I,3)
          T(9)=W(I,9)+z*W(I,3)
          W(I,1)=T(1)
          W(I,2)=T(2)
          W(I,3)=T(3)
          W(I,4)=T(4)
          W(I,5)=T(5)
          W(I,6)=T(6)
          W(I,7)=T(7)
          W(I,8)=T(8)
          W(I,9)=T(9)
 1111   CONTINUE
        RETURN
      ELSEIF(IType.EQ.2200) THEN
        RETURN
      ELSEIF(IType.EQ.2210) THEN
        DO 2210 I=1,NN
          I0 = MOD(I-1,N)+1
          x=BD(1,Ind(2,I0))
          y=BD(2,Ind(2,I0))
          z=BD(3,Ind(2,I0))
          T(1)=W(I,1)
          T(2)=W(I,10)+x*W(I,2)
          T(3)=W(I,11)+y*W(I,2)
          T(4)=W(I,15)+z*W(I,2)
          T(5)=W(I,5)
          T(6)=W(I,11)+x*W(I,3)
          T(7)=W(I,12)+y*W(I,3)
          T(8)=W(I,16)+z*W(I,3)
          T(9)=W(I,9)
          T(10)=W(I,12)+x*W(I,4)
          T(11)=W(I,14)+y*W(I,4)
          T(12)=W(I,18)+z*W(I,4)
          T(13)=W(I,13)
          T(14)=W(I,15)+x*W(I,6)
          T(15)=W(I,16)+y*W(I,6)
          T(16)=W(I,19)+z*W(I,6)
          T(17)=W(I,17)
          T(18)=W(I,16)+x*W(I,7)
          T(19)=W(I,18)+y*W(I,7)
          T(20)=W(I,20)+z*W(I,7)
          T(21)=W(I,21)
          T(22)=W(I,19)+x*W(I,8)
          T(23)=W(I,20)+y*W(I,8)
          T(24)=W(I,22)+z*W(I,8)
          W(I,1) = T(1)
          W(I,2) = T(2)
          W(I,3) = T(3)
          W(I,4) = T(4)
          W(I,5) = T(5)
          W(I,6) = T(6)
          W(I,7) = T(7)
          W(I,8) = T(8)
          W(I,9) = T(9)
          W(I,10) = T(10)
          W(I,11) = T(11)
          W(I,12) = T(12)
          W(I,13) = T(13)
          W(I,14) = T(14)
          W(I,15) = T(15)
          W(I,16) = T(16)
          W(I,17) = T(17)
          W(I,18) = T(18)
          W(I,19) = T(19)
          W(I,20) = T(20)
          W(I,21) = T(21)
          W(I,22) = T(22)
          W(I,23) = T(23)
          W(I,24) = T(24)
 2210   CONTINUE
        RETURN
      ELSEIF(IType.EQ.2211) THEN
        DO 2211 I=1,NN
          I0=MOD(I-1,N)+1
          x=BD(1,Ind(2,I0))
          y=BD(2,Ind(2,I0))
          z=BD(3,Ind(2,I0))
          T(1)=W(I,7)+x*W(I,1)
          T(2)=W(I,8)+y*W(I,1)
          T(3)=W(I,11)+z*W(I,1)
          T(4)=W(I,8)+x*W(I,2)
          T(5)=W(I,9)+y*W(I,2)
          T(6)=W(I,12)+z*W(I,2)
          T(7)=W(I,9)+x*W(I,3)
          T(8)=W(I,10)+y*W(I,3)
          T(9)=W(I,13)+z*W(I,3)
          T(10)=W(I,11)+x*W(I,4)
          T(11)=W(I,12)+y*W(I,4)
          T(12)=W(I,14)+z*W(I,4)
          T(13)=W(I,12)+x*W(I,5)
          T(14)=W(I,13)+y*W(I,5)
          T(15)=W(I,15)+z*W(I,5)
          T(16)=W(I,14)+x*W(I,6)
          T(17)=W(I,15)+y*W(I,6)
          T(18)=W(I,16)+z*W(I,6)
          W(I,1) = T(1)
          W(I,2) = T(2)
          W(I,3) = T(3)
          W(I,4) = T(4)
          W(I,5) = T(5)
          W(I,6) = T(6)
          W(I,7) = T(7)
          W(I,8) = T(8)
          W(I,9) = T(9)
          W(I,10) = T(10)
          W(I,11) = T(11)
          W(I,12) = T(12)
          W(I,13) = T(13)
          W(I,14) = T(14)
          W(I,15) = T(15)
          W(I,16) = T(16)
          W(I,17) = T(17)
          W(I,18) = T(18)
 2211   CONTINUE
        RETURN
      ELSEIF(IType.EQ.2222) THEN
        DO 2222 I=1,NN
          I0 = MOD(I-1,N)+1
          x=BD(1,Ind(2,I0))
          y=BD(2,Ind(2,I0))
          z=BD(3,Ind(2,I0))
          T(1)=W(I,17)+x*W(I,7)+x*(W(I,7)+x*W(I,1))
          T(2)=W(I,18)+x*W(I,8)+y*(W(I,7)+x*W(I,1))
          T(3)=W(I,19)+y*W(I,8)+y*(W(I,8)+y*W(I,1))
          T(4)=W(I,22)+x*W(I,11)+z*(W(I,7)+x*W(I,1))
          T(5)=W(I,23)+z*W(I,8)+y*(W(I,11)+z*W(I,1))
          T(6)=W(I,26)+z*W(I,11)+z*(W(I,11)+z*W(I,1))
          T(7)=W(I,18)+x*W(I,8)+x*(W(I,8)+x*W(I,2))
          T(8)=W(I,19)+y*W(I,8)+x*(W(I,9)+y*W(I,2))
          T(9)=W(I,20)+y*W(I,9)+y*(W(I,9)+y*W(I,2))
          T(10)=W(I,23)+x*W(I,12)+z*(W(I,8)+x*W(I,2))
          T(11)=W(I,24)+y*W(I,12)+z*(W(I,9)+y*W(I,2))
          T(12)=W(I,27)+z*W(I,12)+z*(W(I,12)+z*W(I,2))
          T(13)=W(I,19)+x*W(I,9)+x*(W(I,9)+x*W(I,3))
          T(14)=W(I,20)+y*W(I,9)+x*(W(I,10)+y*W(I,3))
          T(15)=W(I,21)+y*W(I,10)+y*(W(I,10)+y*W(I,3))
          T(16)=W(I,24)+z*W(I,9)+x*(W(I,13)+z*W(I,3))
          T(17)=W(I,25)+y*W(I,13)+z*(W(I,10)+y*W(I,3))
          T(18)=W(I,28)+z*W(I,13)+z*(W(I,13)+z*W(I,3))
          T(19)=W(I,22)+x*W(I,11)+x*(W(I,11)+x*W(I,4))
          T(20)=W(I,23)+x*W(I,12)+y*(W(I,11)+x*W(I,4))
          T(21)=W(I,24)+y*W(I,12)+y*(W(I,12)+y*W(I,4))
          T(22)=W(I,26)+z*W(I,11)+x*(W(I,14)+z*W(I,4))
          T(23)=W(I,27)+z*W(I,12)+y*(W(I,14)+z*W(I,4))
          T(24)=W(I,29)+z*W(I,14)+z*(W(I,14)+z*W(I,4))
          T(25)=W(I,23)+x*W(I,12)+x*(W(I,12)+x*W(I,5))
          T(26)=W(I,24)+y*W(I,12)+x*(W(I,13)+y*W(I,5))
          T(27)=W(I,25)+y*W(I,13)+y*(W(I,13)+y*W(I,5))
          T(28)=W(I,27)+z*W(I,12)+x*(W(I,15)+z*W(I,5))
          T(29)=W(I,28)+z*W(I,13)+y*(W(I,15)+z*W(I,5))
          T(30)=W(I,30)+z*W(I,15)+z*(W(I,15)+z*W(I,5))
          T(31)=W(I,26)+x*W(I,14)+x*(W(I,14)+x*W(I,6))
          T(32)=W(I,27)+y*W(I,14)+x*(W(I,15)+y*W(I,6))
          T(33)=W(I,28)+y*W(I,15)+y*(W(I,15)+y*W(I,6))
          T(34)=W(I,29)+z*W(I,14)+x*(W(I,16)+z*W(I,6))
          T(35)=W(I,30)+z*W(I,15)+y*(W(I,16)+z*W(I,6))
          T(36)=W(I,31)+z*W(I,16)+z*(W(I,16)+z*W(I,6))
          W(I,1) = T(1)
          W(I,2) = T(2)
          W(I,3) = T(3)
          W(I,4) = T(4)
          W(I,5) = T(5)
          W(I,6) = T(6)
          W(I,7) = T(7)
          W(I,8) = T(8)
          W(I,9) = T(9)
          W(I,10) = T(10)
          W(I,11) = T(11)
          W(I,12) = T(12)
          W(I,13) = T(13)
          W(I,14) = T(14)
          W(I,15) = T(15)
          W(I,16) = T(16)
          W(I,17) = T(17)
          W(I,18) = T(18)
          W(I,19) = T(19)
          W(I,20) = T(20)
          W(I,21) = T(21)
          W(I,22) = T(22)
          W(I,23) = T(23)
          W(I,24) = T(24)
          W(I,25) = T(25)
          W(I,26) = T(26)
          W(I,27) = T(27)
          W(I,28) = T(28)
          W(I,29) = T(29)
          W(I,30) = T(30)
          W(I,31) = T(31)
          W(I,32) = T(32)
          W(I,33) = T(33)
          W(I,34) = T(34)
          W(I,35) = T(35)
          W(I,36) = T(36)
 2222   CONTINUE
        RETURN
c     ELSEIF(IType.EQ.3300) THEN
c     ELSEIF(IType.EQ.3310) THEN
c     ELSEIF(IType.EQ.3311) THEN
c     ELSEIF(IType.EQ.3322) THEN
c     ELSEIF(IType.EQ.3333) THEN
      ELSE
         write(*,*) "IType = ",IType
         CALL MondoHalt('Illegal IType in KetHRR')
      ENDIF
      RETURN
      END
