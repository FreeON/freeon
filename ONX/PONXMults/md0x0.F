      SUBROUTINE MD0x0(iBgn,iEnd,IndA,KBra,KKet,BfnBufB,BfnBufK,
     >                 PrmBufB,CBufB,SBufB,PrmBufK,CBufK,SBufK,C)
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "pmonx.inc"
      REAL*8 PrmBufB(MXB,KBra+MInfo,*)
      REAL*8 PrmBufK(MXB,KKet+MInfo,*)
      REAL*8 CBufB(*),SBufB(*)
      REAL*8 CBufK(*),SBufK(*)
      REAL*8 C(*)
      INTEGER IndA(0:*)
      INTEGER BfnBufB(*),BfnBufK(*)
      INCLUDE "bkonx.inc"
  
      KBra3=KBra+3
      KKet3=KKet+3
      I0=0
      DO I=iBgn,iEnd
         Ixa   = BfnBufB(I)
         Px    = PrmBufB(1,KBra3,I)
         Py    = PrmBufB(2,KBra3,I)
         Pz    = PrmBufB(3,KBra3,I)
         rMomP = CBufB(I)
         DO J=1,-IndA(I)
            Ixb = BfnBufK(J)
            IF(Ixa.GE.Ixb) THEN
               I0    = I0+1
               Qx    = PrmBufK(1,KKet3,J)
               Qy    = PrmBufK(2,KKet3,J)
               Qz    = PrmBufK(3,KKet3,J)
               rMomQ = CBufK(J)
               PQ    = DSQRT((Px-Qx)*(Px-Qx)+
     &                       (Py-Qy)*(Py-Qy)+
     &                       (Pz-Qz)*(Pz-Qz))
               C(I0) = rMomP*rMomQ/PQ
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
