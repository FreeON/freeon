      SUBROUTINE MD3x1(N,iBgn,iEnd,IndA,KBra,KKet,LenBf,LenKf,
     >                 BfnBufB,BfnBufK,
     >                 PrmBufB,CBufB,SBufB,PrmBufK,CBufK,SBufK,C)
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "pmonx.inc"
      REAL*8 PrmBufB(MXB,KBra+MInfo,*)
      REAL*8 PrmBufK(MXB,KKet+MInfo,*)
      REAL*8 CBufB(0:9,LenBf,*),SBufB(0:9,LenBf,*)
      REAL*8 CBufK(0:2,LenKf,*),SBufK(0:2,LenKf,*)
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

               INCLUDE "Interact3x1.inc"

               DO KK=1,LenKf
 
                  INCLUDE "ContKet3x1.inc"

                  DO KB=1,LenBf
                     IS=KB+(KK-1)*LenBf

      C(I0,IS)=         CKet(0)*CBufB(0,KB,I)+SKet(0)*SBufB(0,KB,I)
      C(I0,IS)=C(I0,IS)+CKet(1)*CBufB(1,KB,I)+SKet(1)*SBufB(1,KB,I)
      C(I0,IS)=C(I0,IS)+CKet(2)*CBufB(2,KB,I)+SKet(2)*SBufB(2,KB,I)
      C(I0,IS)=C(I0,IS)+CKet(3)*CBufB(3,KB,I)+SKet(3)*SBufB(3,KB,I)
      C(I0,IS)=C(I0,IS)+CKet(4)*CBufB(4,KB,I)+SKet(4)*SBufB(4,KB,I)
      C(I0,IS)=C(I0,IS)+CKet(5)*CBufB(5,KB,I)+SKet(5)*SBufB(5,KB,I)
      C(I0,IS)=C(I0,IS)+CKet(6)*CBufB(6,KB,I)+SKet(6)*SBufB(6,KB,I)
      C(I0,IS)=C(I0,IS)+CKet(7)*CBufB(7,KB,I)+SKet(7)*SBufB(7,KB,I)
      C(I0,IS)=C(I0,IS)+CKet(8)*CBufB(8,KB,I)+SKet(8)*SBufB(8,KB,I)
      C(I0,IS)=C(I0,IS)+CKet(9)*CBufB(9,KB,I)+SKet(9)*SBufB(9,KB,I)


                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
