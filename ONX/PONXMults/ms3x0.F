      SUBROUTINE MS3x0(N,iBgn,iEnd,IndA,KBra,KKet,Mults,LenBf,LenKf,
     >                 ESw,ThreshM,ThreshP,Bfn,BfnBufB,BfnBufK,
     >                 PrmBufB,PrmBufK,CBufB,SBufB,CBufK,SBufK,C)
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "pmonx.inc"
      REAL*8 ESw,ThreshM,ThreshP
      REAL*8 PrmBufB(MXB,KBra+MInfo,*)
      REAL*8 PrmBufK(MXB,KKet+MInfo,*)
      REAL*8 CBufB(*),SBufB(*)
      REAL*8 CBufK(*),SBufK(*)
      REAL*8 C(*)
      INTEGER Bfn(2,*)
      INTEGER BfnBufB(*),BfnBufK(*)
      INTEGER iBgn,iEnd
      INTEGER IndA(0:*)
      LOGICAL Mults
      INCLUDE "bkonx.inc"

c      return
 

      DO I=iBgn,iEnd
         J  = IndA(I)
         Px = PrmBufB(1,KBra+3,I)
         Py = PrmBufB(2,KBra+3,I)
         Pz = PrmBufB(3,KBra+3,I)
         Sx = PrmBufK(4,KKet+3,J)
         Sy = PrmBufK(5,KKet+3,J)
         Sz = PrmBufK(6,KKet+3,J)
         PS = DSQRT((Px-Sx)*(Px-Sx)+
     &              (Py-Sy)*(Py-Sy)+
     &              (Pz-Sz)*(Pz-Sz))
         PQ = PS-PrmBufK(7,KKet+3,J)
c   MAC 1
         IF(PQ.GT.1.0D-2) THEN
            r1xPQ  = 1.0D0/PQ
            r1xPQ2 = r1xPQ*r1xPQ
            r1xPQ3 = r1xPQ*r1xPQ2
            r1xPQ4 = r1xPQ*r1xPQ3
            r1xPQ5 = r1xPQ*r1xPQ4

            rep = PrmBufB(1,KBra+1,I)
            ep0 = PrmBufB(2,KBra+1,I)
            ep1 = PrmBufB(3,KBra+1,I)
            ep2 = PrmBufB(4,KBra+1,I)
            ep3 = PrmBufB(5,KBra+1,I)
            req = PrmBufK(1,KKet+2,J)
            eq0 = PrmBufK(2,KKet+2,J)

            rM2 = r1xPQ2*req*ep0 +
     >            r1xPQ3*req*ep1*2.0D0 +
     >            r1xPQ4*req*ep2*6.0D0 +
     >            r1xPQ5*(req*ep3*2.4D1+rep*eq0)

c   MAC 2
            IF(rM2.GT.ThreshM) GOTO 2000
            Zp = PrmBufB(8,KBra+3,I)
            Zq = PrmBufK(9,KKet+3,J)
            T  = Zp*Zq/(Zp+Zq)*PQ*PQ
c   PAC 1
            IF(T.GT.ESw) GOTO 1000
            ST=DSQRT(T)
            epm=PrmBufB(MXB,KBra+1,I)
            rP2=epm*eq0*r1xPQ*DERFC(ST)
c   PAC 2
            IF(rP2.GT.ThreshP) GOTO 2000
 1000       CONTINUE
            IndA(I)=-IndA(I)
            Mults=.TRUE.
         ENDIF
 2000    CONTINUE
      ENDDO

      IF(Mults) THEN
         N=0
         DO I=iBgn,iEnd
            Ixa = BfnBufB(I)
            DO J=1,-IndA(I)
               Ixb = BfnBufK(J)
               IF(Ixa.GE.Ixb) THEN
                  N=N+1
                  Bfn(1,N)=Ixa
                  Bfn(2,N)=Ixb
               ENDIF
            ENDDO
         ENDDO
         IF(N.EQ.0) THEN
            Mults=.FALSE.
         ELSE
            CALL MD3x0(N,iBgn,iEnd,IndA,KBra,KKet,LenBf,LenKf,
     &                 BfnBufB,BfnBufK,
     &                 PrmBufB,CBufB,SBufB,
     &                 PrmBufK,CBufK,SBufK,C)
         ENDIF
      ENDIF

      RETURN
      END
