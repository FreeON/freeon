      SUBROUTINE ONXSch(rBuf,NASym,Max2L,Max4L,Len2L,
     $                  AuxR,LAmin,LAmax,LBmin,LBmax,
     $                  Ax,Ay,Az,Bx,By,Bz,AB2,
     $                  IType,JType,NICase,NKCase,LenC,
     $                  LMNdx,IndL,IndM,IndN,
     $                  EBra,Cr,Sr,
     $                  Ex,Ey,Ez,R,rInt,Thresh)
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "pmonx.inc"
      REAL*8 rBuf(MXB,*)
      REAL*8 rBufT(MXB,MXCNT*MXCNT+MInfo)
      REAL*8 Ax,Ay,Az,Bx,By,Bz,AB2
      REAL*8 EBra(len1,*)
      REAL*8 EBraT(len1)
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C     accumulated moments:
C
c      LenSP=NASym*(NASym+3)/2  ;  TOTAL LENGTH=LenSP+1
c      LenSP2=NASym*(2*NASym+3)
C
      REAL*8 Cr(0:LenS,*)
      REAL*8 Sr(0:LenS,*)
C
      PARAMETER (MaxEll=12)
      PARAMETER (LenEll=(MaxEll+3)*MaxEll/2)
C
C     Work arrays:
C
      REAL*8 Factorial(0:2*MaxEll)
      REAL*8 FactOlm0(0:MaxEll)
      REAL*8 FactOlm1(0:MaxEll)
      REAL*8 FactOlm2(0:MaxEll,0:MaxEll)
      REAL*8 FactMlm0(0:MaxEll)
      REAL*8 FactMlm1(0:MaxEll)
      REAL*8 FactMlm2(0:MaxEll,0:MaxEll)
C
      REAL*8 Cosine(  0:MaxEll)
      REAL*8 Sine(    0:MaxEll)
      REAL*8 RToTh(   0:MaxEll)
      REAL*8 LegendreP(0:MaxEll,0:MaxEll)
C
      REAL*8 Cpr(0:LenEll)
      REAL*8 Spr(0:LenEll)
C
C     temporary distribution moments:
C
      REAL*8 Cp(0:LenEll)
      REAL*8 Sp(0:LenEll)
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      Real*8 Etmp1(0:8)
      REAL*8 Etmp2(0:8)
      REAL*8 Ex(0:NASym,0:NASym,0:2*NASym)
      REAL*8 Ey(0:NASym,0:NASym,0:2*NASym)
      REAL*8 Ez(0:NASym,0:NASym,0:2*NASym)
      REAL*8 R(0:Max4L,0:Max4L,0:Max4L)
      REAL*8 AuxR(0:Max4L)

      REAL*8 rInt
      REAL*8 CC(10)
      REAL*8 rf(0:10),fM(0:4,0:4)
      INTEGER LAmin,LAmax,LBmin,LBmax
      INTEGER IType,JType
      INTEGER LMNdx(0:Max2L,0:Max2L,0:Max2L)
      INTEGER IndL(Len2L)
      INTEGER IndM(Len2L)
      INTEGER IndN(Len2L)
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      LenSP(L)=L*(L+3)/2
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      DATA fM / 1.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,
     >          1.0D0, 2.0D0, 0.0D0, 0.0D0, 0.0D0,
     >          1.0D0, 3.0D0, 6.0D0, 0.0D0, 0.0D0,
     >          1.0D0, 4.0D0, 1.2D1, 2.4D1, 0.0D0,
     >          1.0D0, 5.0D0, 2.0D1, 6.0D1, 1.20D2 /
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      CALL MultipoleSetUp(MaxEll,Factorial,
     >                    FactOlm0,FactOlm1,FactOlm2,
     >                    FactMlm0,FactMlm1,FactMlm2)
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      Ltot=LAmax+LBmax
      Lmx =Ltot
      Ltot2=2*Ltot
      rInt=0.0D0
      iCnt=0
      Zmin=1.0D20
      iSA=LBegin(LAmin)
      iSB=LBegin(LBmin)
      iEA=LEnd(LAmax)  
      iEB=LEnd(LBmax)
      rf(0)=1.0D0
      DO I=1,Ltot
         rf(I)=rf(I-1)*DFLOAT(I)
      ENDDO
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
      I=len0*len0*(LenS+1)
      CALL RLoader(Cr,I,0.0D0)
      CALL RLoader(Sr,I,0.0D0)
      DO I=0,8
         Etmp1(I) = 0.0D0
         Etmp2(I) = 0.0D0
      ENDDO
c
c Find a good place for a multipole expansion.
c
      SumX=0.0D0
      SumY=0.0D0
      SumZ=0.0D0
      SumM=0.0D0
      DO IK=1,LenC
         Zeta = rBuf(1,IK)
         Px   = rBuf(2,IK)
         Py   = rBuf(3,IK)
         Pz   = rBuf(4,IK)
         rmom = rBuf(5,IK)*(rpi/Zeta)**(3.0D0/2.0D0)
         SumX = SumX+ABS(rmom)*Px
         SumY = SumY+ABS(rmom)*PY
         SumZ = SumZ+ABS(rmom)*PZ
         SumM = SumM+ABS(rmom)
         Zmin = MIN(Zmin,Zeta)
      ENDDO
      Rx = SumX/SumM
      Ry = SumY/SumM
      Rz = SumZ/SumM

      DO 5000 IK=1,LenC
         Zeta = rBuf(1,IK)
         Px   = rBuf(2,IK)
         Py   = rBuf(3,IK)
         Pz   = rBuf(4,IK)
         VSAB = rBuf(10,IK)
         IF(DABS(VSAB).LE.1.0D-15) GOTO 4990
         DO icc=1,NICase*NKCase
            CC(icc)=rBuf(5+icc,IK)
         ENDDO
         r1xZ=re1/Zeta
         rPxZ=(rpi*r1xZ)**(3.0D0/2.0D0)
         PRx = Px-Rx
         PRy = Py-Ry
         PRz = Pz-Rz
         PAx = Px-Ax
         PAy = Py-Ay
         PAz = Pz-Az
         PBx = Px-Bx
         PBy = Py-By
         PBz = Pz-Bz
         W   = Zeta*0.5D0
         Pre = Prev/(Zeta*Zeta*DSQRT(Zeta+Zeta))        
         PR  = DSQRT(PRx*PRx+PRy*PRy+PRz*PRz)
         CALL MD2TRRx(NASym,0,LAmax,LBmax,Zeta,
     >                PAx,PBx,PAy,PBy,PAz,PBz,Ex,Ey,Ez)
         jp = 0
         DO indA=iSA,iEA
            La=IndL(indA)
            Ma=IndM(indA)
            Na=IndN(indA)
            LtotA=La+Ma+Na
            ic1=LtotA-LAmin+1
            DO indB=iSB,iEB
               Lb=IndL(indB)
               Mb=IndM(indB)
               Nb=IndN(indB)
               LtotB=Lb+Mb+Nb
               ic2=LtotB-LBmin+1
               rccAB=CC(ic1+(ic2-1)*NICase)

               jp=jp+1
               IF(IType.GE.JType) THEN
                  jsp=jp
               ELSE
                  jsp=(indA-iSA+1)+(indB-iSB)*(iEA-iSA+1)
               ENDIF
 
               Ltot=La+Ma+Na+Lb+Mb+Nb
               DO LMN=1,Len1
                  EBra(LMN,jp)=0.0D0
                  EBraT(LMN)=0.0D0
               ENDDO
               DO Lab=0,La+Lb
                  EEx=VSAB*Ex(La,Lb,Lab)
               DO Mab=0,Ma+Mb
                  EExy=EEx*Ey(Ma,Mb,Mab)
               DO Nab=0,Na+Nb
                  LMN=LMNdx(Lab,Mab,Nab)
                  EBra(LMN,jp)=EExy*Ez(Na,Nb,Nab)
                  EBraT(LMN)=EBra(LMN,jp)*rccAB
                  IP=Lmx-Ltot+1
                  rtm=ABS(EBraT(LMN)*rPxZ*PR**IP)
                  Etmp1(Ltot)=MAX(Etmp1(Ltot),rtm)
                  Etmp2(Ltot)=MAX(Etmp2(Ltot),ABS(EBraT(LMN)*rPxZ))
               ENDDO
               ENDDO
               ENDDO
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
               LP=Ltot
               LR=Lmx
               LPR=Ltot+Lmx
               LenP=LenSP(LP)
               LenR=LenS
               LenPR=LenSP(LP+LR)             
               CALL HGTFToSP(Ltot,rpxz,Cp,Sp,EBraT)      
               CALL XLate(MaxEll,LP,LR,LPR,LenP,LenR,LenPR,
     >                    PRx,PRy,PRz,Cr(0,jsp),Sr(0,jsp),
     >                    Cpr,Spr,Cp,Sp,LegendreP,Cosine,Sine,RToTh,
     >                    FactOlm0,FactOlm1,FactOlm2)
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
            ENDDO
         ENDDO

         romg1=1.0D0
         romg2=-2.0D0*W
         DO J=0,Ltot2
            AuxR(J)=Pre*romg1/DFLOAT(2*J+1)
            romg1=romg1*romg2
         ENDDO

         DO J=LTot2,0,-1
            DO I0=LTot2-J,1,-1
               IF(I0-1.LE.0)THEN
                  R(I0,0,0)=0.0D0
               ELSE
                  R(I0,0,0)=DFLOAT(I0-1)*R(I0-2,0,0)
               ENDIF
               DO I1=LTot2-J-I0,1,-1
                  IF(I1-1.LE.0)THEN
                     R(I1,I0,0)=0.0D0
                  ELSE
                     R(I1,I0,0)=DFLOAT(I1-1)*R(I1-2,I0,0)
                  ENDIF
                  DO I2=LTot2-J-I0-I1,1,-1
                     IF(I2-1.LE.0)THEN
                        R(I2,I1,I0)=0.0D0
                     ELSE
                        R(I2,I1,I0)=DFLOAT(I2-1)*R(I2-2,I1,I0)
                     ENDIF
                  ENDDO
                  IF(I1-1.LE.0)THEN
                     R(0,I1,I0)=0.0D0
                     R(I1,0,I0)=0.0D0
                  ELSE
                     R(0,I1,I0)=DFLOAT(I1-1)*R(0,I1-2,I0)
                     R(I1,0,I0)=DFLOAT(I1-1)*R(I1-2,0,I0)
                  ENDIF
               ENDDO
               IF(I0-1.LE.0)THEN
                  R(0,I0,0)=0.0D0
                  R(0,0,I0)=0.0D0
               ELSE
                  R(0,I0,0)=DFLOAT(I0-1)*R(0,I0-2,0)
                  R(0,0,I0)=DFLOAT(I0-1)*R(0,0,I0-2)
               ENDIF
            ENDDO
            R(0,0,0)=AuxR(J)
         ENDDO

         jp=0
         sval = 0.0D0
         DO indA=iSA,iEA
            La=IndL(indA)
            Ma=IndM(indA)
            Na=IndN(indA)
            LtotA=La+Ma+Na
            ic1=LtotA-LAmin+1
            DO indB=iSB,iEB
               Lb=IndL(indB)
               Mb=IndM(indB)
               Nb=IndN(indB)
               LtotB=Lb+Mb+Nb
               ic2=LtotB-LBmin+1
               rccAB=CC(ic1+(ic2-1)*NICase)
               jp=jp+1

               r0=re0
               DO K=0,LTot
                  DO J=0,LTot-k
                     DO I=0,LTot-K-J
                        lmn1=LMNdx(I,J,K)
                        Phase=(-1.0D0)**(I+J+K)
                        TMP=Phase*EBra(lmn1,jp)
               DO N=0,LTot
                  DO M=0,LTot-N
                     DO L=0,LTot-N-M
                        lmn2=LMNdx(L,M,N)
                        r0=r0+TMP*EBra(lmn2,jp)*R(I+L,J+M,K+N)
                     ENDDO
                  ENDDO
               ENDDO
                     ENDDO
                  ENDDO
               ENDDO
               sval=MAX(sval,DSQRT(DABS(r0))*DABS(rccAB))
            ENDDO
         ENDDO

         IF(sval.GE.Thresh*1.0D-2) THEN
c         if(.true.) then
           iCnt=iCnt+1
           DO I=1,MXB
              rBufT(I,iCnt)=rBuf(I,IK)
           ENDDO
           rBufT(9,iCnt)=sval
         ENDIF

         rInt=rInt+sval

 4990    CONTINUE
 5000 CONTINUE  

      LenC=iCnt
      DO IK=1,iCnt
         BigInt=0.0D0
         DO I0=1,iCnt
            IF(rBufT(9,I0).GT.BigInt) THEN
               Iposs=I0
               BigInt=rBufT(9,I0)
            ENDIF
         ENDDO
         DO I=1,MXB
            rBuf(I,IK)=rBufT(I,Iposs)
         ENDDO
         rBuf(9,IK)=rBuf(9,IK)*DSQRT(DFLOAT(LenC))
         rBufT(9,IPoss)=-999.9D0
      ENDDO

      rBuf(1,LenC+1)=0.0D0
      DO I=0,Lmx
         rBuf(1,LenC+1)=rBuf(1,LenC+1)+fM(I,Lmx)*Etmp1(I)
         rBuf(2+I,LenC+1)=Etmp2(I)
      ENDDO
      rBuf(1,LenC+3)=Rx
      rBuf(2,LenC+3)=Ry
      rBuf(3,LenC+3)=Rz
      rBuf(8,LenC+3)=Zmin

      DO I=Lmx+3,MXB
         rBuf(I,LenC+1)=0.0D0
         rBuf(I,LenC+2)=0.0D0
      ENDDO

      RETURN
      END
