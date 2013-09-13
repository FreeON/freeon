      SUBROUTINE XLateV(NcD,Nc,EllP,EllQ,Len,Sgn,Idx,
     >                  FactOlm0,FactOlm2,P,Q,Cq,Sq,Cp,Sp)
C
      IMPLICIT NONE
C
      INTEGER FEll,NcD,Nc,EllP,EllQ
      INTEGER Ell,Len,I,L
C
      REAL*8 FactOlm0(0:*)
      REAL*8 FactOlm2(0:*)
C
      REAL*8 Q(3,*),P(3)
      REAL*8 Sgn(4,*)
      INTEGER Idx(3,*)
      REAL*8 Cp(0:*),Sp(0:*)
      REAL*8 Cq(NcD,0:*),Sq(NcD,0:*)
C
      INTEGER iClust,MxEll,MxLen
      PARAMETER (iClust=2048)
      PARAMETER (MxEll=20)
      PARAMETER (MxLen=MxEll*(MxEll+3)/2+1)
C
      REAL*8 X(iClust),Y(iClust),Z(iClust),CoTan(iClust)
      REAL*8 LegendreP(iClust*(MxEll+1)**2)
      REAL*8 Sine(iClust*(MxEll+1))
      REAL*8 Cosine(iClust*(MxEll+1))
      REAL*8 RToTh(iClust*(MxEll+1))
      REAL*8 Cpq(iClust*(MxLen+1)),Spq(iClust*(MxLen+1))
C
      INTEGER LSP,LenSP
      LSP(L)=L*(L+3)/2
C
      WRITE(*,*)' Len = ',Len
      STOP


      Ell=MAX(EllP,EllQ)
      LenSP=LSP(L)
C
      IF(Nc.GT.iClust)STOP '1 in XLateV '
      IF(NcD.GT.iClust)STOP '1.5 in XLateV '
      IF(Ell.GT.MxEll)STOP '2 in XLateV '
      IF(LenSP.GT.MxLen)STOP '3 in XLateV '
C
      DO I=1,Nc
         X(I)=(Q(1,I)-P(1))
         Y(I)=(Q(2,I)-P(2))
         Z(I)=(Q(3,I)-P(3))
      ENDDO
C
      CALL RegularV(Nc,Ell,X,Y,Z,CoTan,
     >              LegendreP,Cosine,Sine,RToTh,
     >              FactOlm0,FactOlm2,Cpq,Spq)
C
      CALL TranslateV(NcD,Nc,Len,Sgn,Idx,Cq,Sq,Cpq,Spq,Cp,Sp)
C
      RETURN
      END
C
      SUBROUTINE TranslateV(NcD,Nc,Len,Sgn,Idx,Cq,Sq,Cpq,Spq,Cp,Sp)
C
      INTEGER NcD,Nc,Len,I,J
C
      REAL*8 Sgn(4,*),Cq(NcD,0:*),Sq(NcD,0:*)
      INTEGER Idx(3,*)
      REAL*8 Cpq(Nc,0:*),Spq(Nc,0:*),Cp(0:*),Sp(0:*)
C
      REAL*8 Cdum,Sdum
      INTEGER idx1,idx2,idx3
      REAL*8 Flt1,Flt2,Flt3,Flt4
C
      DO J=1,Len
         idx1=Idx(1,J)
         idx2=Idx(2,J)
         idx3=Idx(3,J)
         Flt1=Sgn(1,J)
         Flt2=Sgn(2,J)
         Flt3=Sgn(3,J)
         Flt4=Sgn(4,J)
         Cdum=0D0
         Sdum=0D0

c         WRITE(*,33)idx1,idx2,idx3,flt1,flt2
c 33      FORMAT('B ',3(I3,', '),2(d12.6))

         DO I=1,Nc

          WRITE(*,22)I,J,idx1,idx2,idx3,Cq(I,idx2),
     >              Cpq(I,idx3),Sq(I,idx2),Spq(I,idx3)
 22       FORMAT(5(I4,", "),4(D12.6,", "))



            Cdum=Cdum+Flt1*Cq(I,idx2)*Cpq(I,idx3)
     >               +Flt2*Sq(I,idx2)*Spq(I,idx3)
            Sdum=Sdum+Flt3*Cq(I,idx2)*Spq(I,idx3)
     >               +Flt4*Sq(I,idx2)*Cpq(I,idx3)
         ENDDO
         Cp(idx1)=Cp(idx1)+Cdum
         Sp(idx1)=Sp(idx1)+Sdum
      ENDDO
C
      RETURN
      END
C
      SUBROUTINE RegularV(Nc,Ell,x,y,z,CoTan,
     $                    LegendreP,Cosine,Sine,RToTh,
     $                    FactOlm0,FactOlm2,C,S)
C
      IMPLICIT NONE
C
      INTEGER FEll,Nc,Ell,Len
C
      REAL*8 X(*),Y(*),Z(*),CoTan(*)
      REAL*8 C(Nc,0:*),S(Nc,0:*)
C
      REAL*8 LegendreP(Nc,0:Ell,0:*)
      REAL*8 Cosine(Nc,0:*)
      REAL*8 Sine(Nc,0:*)
      REAL*8 RToTh(Nc,0:*)
      REAL*8 FactOlm0(0:*)
      REAL*8 FactOlm2(0:*)
C
      INTEGER I,L,M,M1,M2,LDX,LM,LDex0
      REAL*8  X2,Y2,R,Rxy,OneOvRxy,TwoC,RL,RS,Sq,Flt
C
      INTEGER LSP,LTD
      LSP(L)=L*(L+3)/2
      LTD(L)=L*(L+1)/2
C
C      REAL*8 RSum(100)
C
C      DO I=1,Nc
C         RSum(I)=0D0
C      ENDDO

      DO I=1,Nc
C
         Sine(I,0)=0D0
         Cosine(I,0)=1D0
C
         x2=x(I)*x(I)
         y2=y(I)*y(I)
         R=DSQRT(x2+y2+z(I)*z(I))
         IF(R==0D0)THEN
            Sine(I,1)=0D0
            Cosine(I,1)=0D0
            CoTan(I)=0D0
         ELSE
            CoTan(I)=z(I)/R
            Rxy=DSQRT(x2+y2)
            IF(Rxy.NE.0D0)THEN
               OneOvRxy=1D0/DSQRT(x2+y2)
               Sine(I,1)=x(I)*OneOvRxy
               Cosine(I,1)=y(I)*OneOvRxy
            ELSE
               Sine(I,1)=0.70710678118654752D0
               Cosine(I,1)=0.70710678118654752D0
            ENDIF
         ENDIF
C
         TwoC=2D0*Cosine(I,1)
         DO m=2,Ell
            m1=m-1
            m2=m-2
            Sine(I,m)=TwoC*Sine(I,m1)-Sine(I,m2)
            Cosine(I,m)=TwoC*Cosine(I,m1)-Cosine(I,m2)

C            RSum(I)=RSum(I)+Sine(I,m)**2+Cosine(I,m)**2

         ENDDO
C
         Rl=1.0D0
         RS=1.0D0
         Sq=DSQRT(1.0D0-CoTan(I)*CoTan(I))
         DO m=0,Ell
            LegendreP(I,m,m)=FactOlm0(m)*RS
            RToTh(I,m)=Rl
            Rl=Rl*R
            RS=RS*Sq
         ENDDO
C
         DO m=0,Ell-1
            LegendreP(I,m,m+1)=CoTan(I)*LegendreP(I,m,m)
         ENDDO
C
      ENDDO
C
      DO l=2,Ell
         Flt=DBLE(2*l-1)
         LDex0=LTD(L)
         DO m=0,l-2
            DO I=1,Nc
               LegendreP(I,m,l)=(LegendreP(I,m,l-1)*CoTan(I)*Flt
     $                         - LegendreP(I,m,l-2))*FactOlm2(LDex0+M)
            ENDDO
         ENDDO
      ENDDO
C
      DO l=0,ell
         ldx=(l+1)*l/2
         DO m=0,l
            lm=ldx+m
            DO I=1,Nc
               S(I,lm)=RToTh(I,l)*LegendreP(I,m,l)*Sine(I,m)
               C(I,lm)=RToTh(I,l)*LegendreP(I,m,l)*Cosine(I,m)
C               RSUM(I)=RSUM(I)+C(I,lm)**2+S(I,lm)**2
            ENDDO
         ENDDO
      ENDDO
C
C      DO I=1,Nc
C         WRITE(*,*)X(I),Y(I),Z(I)
C         WRITE(*,22)(C(I,M),M=0,LSP(Ell))
C 22      FORMAT('B ',100(D12.6,", "))
C         WRITE(*,*)'B',I,Ell,RSum(I)
C      ENDDO
C
      RETURN
      END
