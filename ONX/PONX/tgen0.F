      SUBROUTINE TGen0(Ltot,KBra,KKet,N,Ind,
     >                 PrmBufB,PrmBufK,R)
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "pmonx.inc"
      REAL*8 PrmBufB(MXB,KBra+MInfo,*)
      REAL*8 PrmBufK(MXB,KKet+MInfo,*)
      REAL*8 R(*)
      REAL*8 F0_0(0:486),F0_1(0:486),F0_2(0:486),F0_3(0:486)
      INTEGER Ind(2,*)
      INCLUDE "Gamma_A.inc"
      INCLUDE "Gamma_0.inc"

      I1=0
      DO IKK=1,KKet
         DO IKB=1,KBra
            DO I0=1,N

               I=Ind(1,I0)
               J=Ind(2,I0)
               I1=I1+1

               Eta = PrmBufK(1,IKK,J)
               Zeta = PrmBufB(1,IKB,I)
               r1xZpE=1.0D0/(Zeta+Eta)

               PQx = PrmBufB(2,IKB,I)-PrmBufK(2,IKK,J)
               PQy = PrmBufB(3,IKB,I)-PrmBufK(3,IKK,J)
               PQz = PrmBufB(4,IKB,I)-PrmBufK(4,IKK,J)

               T=(PQx*PQx+PQy*PQy+PQz*PQz)*
     &            Eta*Zeta*r1xZpE

               Rkk = PrmBufB(5,IKB,I)*PrmBufK(5,IKK,J)
 
               IF(T.LT.Switch) THEN
                  m=INT(T*Grid)
                  T2=T*T
                  T3=T2*T
                  G0=F0_0(m)+T*F0_1(m)+T2*F0_2(m)+T3*F0_3(m)
                  R(I1)=Rkk*G0*DSQRT(r1xZpE)
               ELSE
                  R(I1)=Rkk*F0Asymp*DSQRT(r1xZpE/T)
               ENDIF

            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
