      SUBROUTINE MD0x4(N,iBgn,iEnd,IndA,KBra,KKet,LenBf,LenKf,
     >                 BfnBufB,BfnBufK,
     >                 PrmBufB,CBufB,SBufB,PrmBufK,CBufK,SBufK,C)
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "pmonx.inc"
      REAL*8 PrmBufB(MXB,KBra+MInfo,*)
      REAL*8 PrmBufK(MXB,KKet+MInfo,*)
      REAL*8 CBufB(0: 0,LenBf,*),SBufB(0: 0,LenBf,*)
      REAL*8 CBufK(0:14,LenKf,*),SBufK(0:14,LenKf,*)
      REAL*8 C(N,*)
      INTEGER IndA(0:*)
      INTEGER BfnBufB(*),BfnBufK(*)
      INCLUDE "bkonx.inc"
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      PARAMETER (MaxEll=24)
      PARAMETER (LenEll=(MaxEll+3)*MaxEll/2)
C
C     Work arrays:
C
      REAL*8 CKet(0:LenEll)
      REAL*8 SKet(0:LenEll)
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      KBra3=KBra+3
      KKet3=KKet+3
      I0=0
      DO I=iBgn,iEnd
         Ixa   = BfnBufB(I)
         Px    = PrmBufB(1,KBra3,I)
         Py    = PrmBufB(2,KBra3,I)
         Pz    = PrmBufB(3,KBra3,I)
         DO J=1,-IndA(I)
            Ixb = BfnBufK(J)
            IF(Ixa.GE.Ixb) THEN
               I0    = I0+1
               PQx   = Px - PrmBufK(1,KKet3,J)
               PQy   = Py - PrmBufK(2,KKet3,J)
               PQz   = Pz - PrmBufK(3,KKet3,J)

               INCLUDE "Interact0x4.inc"

               DO KK=1,LenKf

                  INCLUDE "ContKet0x4.inc"

                  C(I0,KK)=CKet(0)*CBufB(0,1,I)

               ENDDO
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
