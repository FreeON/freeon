      SUBROUTINE ONX1111(N,iBgn,iEnd,Indx,IndA,Bfn,BfnBufB,BfnBufK,
     >                   KBra,KKet,
     >                   PrmBufB,PrmBufK,C)
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "pmonx.inc"
      REAL*8 C(N)
      REAL*8 F0_0(0:486),F0_1(0:486),F0_2(0:486),F0_3(0:486)
      REAL*8 PrmBufB(MXB,KBra+MInfo,*)
      REAL*8 PrmBufK(MXB,KKet+MInfo,*)
      integer bfn(2,N)
      INTEGER Indx(2,MXBATCH),IndA(0:*)
      INTEGER BfnBufB(*),BfnBufK(*)
      INCLUDE "bkonx.inc"
      INCLUDE "Gamma_0.inc"
      INCLUDE "Gamma_A.inc"
    
      Tsh0=Tsh/Dcd
      DO I0=1,N
         C(I0)=0.0D0
         I=Indx(1,I0)
         J=Indx(2,I0)
         DO 2000 IKK=1,KKet
            Tsh2 = Tsh0/PrmBufK(9,IKK,J)
            Eta  = PrmBufK(1,IKK,J)
            Qx   = PrmBufK(2,IKK,J)
            Qy   = PrmBufK(3,IKK,J)
            Qz   = PrmBufK(4,IKK,J)
            Uq   = PrmBufK(5,IKK,J)
            rtmp = 0.0D0
            DO 1000 IKB=1,KBra
               IF(PrmBufB(9,IKB,I).LE.Tsh2) GOTO 1010
               Zeta   = PrmBufB(1,IKB,I)
               PQx    = PrmBufB(2,IKB,I) - Qx
               PQy    = PrmBufB(3,IKB,I) - Qy
               PQz    = PrmBufB(4,IKB,I) - Qz
               Up     = PrmBufB(5,IKB,I)
               r1xZpE = 1.0D0/(Zeta+Eta)
               T=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
               IF(T.lt.Switch) THEN
                  m=INT(T*Grid)
                  T2=T*T
                  T3=T2*T
                  R1=F0_0(m)+T*F0_1(m)+T2*F0_2(m)+T3*F0_3(m)
               ELSE
                  R1=F0Asymp/DSQRT(T)
               ENDIF
               rtmp=rtmp+R1*Up*DSQRT(r1xZpE)
 1000       CONTINUE
 1010       CONTINUE
            C(I0)=C(I0)+rtmp*Uq
 2000    CONTINUE
 2010    CONTINUE
      ENDDO
      RETURN
      END
