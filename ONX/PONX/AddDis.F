      SUBROUTINE AddDis(Num,NT,KLng,I1,I2,iPrm,iDis,iMom,SchT,
     &                  BfnT,BufT,VecT,
     &                  TmpBuf,PrmBuf,DisBuf,BfnBuf,VecBuf,
     >                  Ctmp,Stmp,CBuf,SBuf)
      IMPLICIT REAL*8 (a-h,o-z)
      INCLUDE "pmonx.inc"
      INTEGER Num,NT,iPrm,iDis
      INTEGER BfnT(500),BufT(Num)
      INTEGER BfnBuf(Num)
      REAL*8 SchT(Num)
      REAL*8 VecT(3,NT)
      REAL*8 TmpBuf(MXB,MXCNT*MXCNT+MInfo,NT)
      REAL*8 Ctmp(0:LenS,Len0*Len0,NT)
      REAL*8 Stmp(0:LenS,Len0*Len0,NT)
      REAL*8 PrmBuf(MXB,KLng+MInfo,Num)
      REAL*8 DisBuf(Num),VecBuf(3,Num)
      REAL*8 CBuf(0:I1,I2,Num)
      REAL*8 SBuf(0:I1,I2,Num)
      REAL*8 Etmp(MXB)

      DO I=1,MXB
         Etmp(I)=0.0D0
      ENDDO
      ZetaMin=1.0D20
      SumX=0.0D0
      SumY=0.0D0
      SumZ=0.0D0
      I0=0
c
c Loop over distributions in a sphere
c
      DO I=1,Num
         DisBuf(I)=SchT(I)
         ILoc     =BufT(I)
c
c Primitive information, exponents, centers, ect.
c
         DO IK=1,KLng
            DO J=1,MXB
               PrmBuf(J,IK,I)=TmpBuf(J,IK,ILoc)
            ENDDO
         ENDDO
         BfnBuf(I)=BfnT(ILoc)
         VecBuf(1,I)=VecT(1,ILoc)
         VecBuf(2,I)=VecT(2,ILoc)
         VecBuf(3,I)=VecT(3,ILoc)
c
c Moments of distributions (translated)
c
         DO J=1,I2
            DO K=0,I1
               I0=I0+1
               CBuf(K,J,I)=Ctmp(K,J,ILoc)
               SBuf(K,J,I)=Stmp(K,J,ILoc)
            ENDDO
         ENDDO
c
c MAC+PAC distribution (KLng+1) and sphere (KLng+2) info
c
         DO J=1,MXB
            PrmBuf(J,KLng+1,I)=TmpBuf(J,KLng+1,ILoc)
            PrmBuf(J,KLng+2,I)=MAX(Etmp(J),TmpBuf(J,KLng+1,ILoc))
            PrmBuf(J,KLng+3,I)=TmpBuf(J,KLng+3,ILoc)
            Etmp(J)=PrmBuf(J,KLng+2,I)
         ENDDO
         rDisM=0.0D0
         rSphM=0.0D0
         DO J=2,MXB-1
            rDisM=MAX(rDisM,PrmBuf(J,KLng+1,I))
            rSphM=MAX(rSphM,PrmBuf(J,KLng+2,I))
         ENDDO
         PrmBuf(MXB,KLng+1,I)=rDisM
         PrmBuf(MXB,KLng+2,I)=rSphM
c
c The smallest exponent in the sphere
c
         PrmBuf(9,KLng+3,I)=MIN(ZetaMin,PrmBuf(8,KLng+3,I))
         ZetaMin=PrmBuf(9,KLng+3,I)
c
c Find the sphere centers
c
         SumX=SumX+PrmBuf(1,KLng+3,I)
         SumY=SumY+PrmBuf(2,KLng+3,I)
         SumZ=SumZ+PrmBuf(3,KLng+3,I)
         PrmBuf(4,KLng+3,I)=SumX/FLOAT(I)
         PrmBuf(5,KLng+3,I)=SumY/FLOAT(I)
         PrmBuf(6,KLng+3,I)=SumZ/FLOAT(I)
      ENDDO
c
c Loop back over them to find the sphere radius
c
      DO I=1,Num
         Sx=PrmBuf(4,KLng+3,I)
         Sy=PrmBuf(5,KLng+3,I)
         Sz=PrmBuf(6,KLng+3,I)
         rad2=0.0D0
         DO J=1,I
            PSx=PrmBuf(1,KLng+3,J)-Sx
            PSy=PrmBuf(2,KLng+3,J)-Sy
            PSz=PrmBuf(3,KLng+3,J)-Sz
            r2=PSx*PSx+PSy*PSy+PSz*PSz
            rad2=MAX(r2,rad2)
         ENDDO
         PrmBuf(7,KLng+3,I)=DSQRT(rad2)
      ENDDO

      iPrm = iPrm + Num*MXB*(KLng+MInfo)
      iDis = iDis + Num
      iMom = iMom + I0

      RETURN
      END
