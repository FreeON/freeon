      SUBROUTINE RGen1C(Ltot,Kon,WK,CD,WR,WZ,R,
     >                  PrmBuf,DisBuf)
      IMPLICIT REAL*8 (a-h,o-z)
      IMPLICIT INTEGER (i-n)
      INCLUDE "pmsonx.inc"
      REAL*8 WK(4,Kon)
      REAL*8 CD(Kon,3)
      REAL*8 WR(6,Kon,Kon)
      REAL*8 WZ(3,Kon,Kon)
      REAL*8 R(Kon,Kon,*)
      REAL*8 PrmBuf(MXP,Kon+MInfo)
      REAL*8 DisBuf(MXC)

      DO IK=1,Kon
         CD(IK,1)=PrmBuf(6,IK)
         CD(IK,2)=PrmBuf(7,IK)
         CD(IK,3)=PrmBuf(8,IK)
         WK(1,IK)=PrmBuf(2,IK)-DisBuf(8)
         WK(2,IK)=PrmBuf(3,IK)-DisBuf(9)
         WK(3,IK)=PrmBuf(4,IK)-DisBuf(10)
         WK(4,IK)=0.5D0/PrmBuf(1,IK)
         DO IB=1,Kon
            Zeta = PrmBuf(1,IB)
            Px   = PrmBuf(2,IB)
            Py   = PrmBuf(3,IB)
            Pz   = PrmBuf(4,IB)
            Up   = PrmBuf(5,IB)
            Eta  = PrmBuf(1,IK)
            Qx   = PrmBuf(2,IK)
            Qy   = PrmBuf(3,IK)
            Qz   = PrmBuf(4,IK)
            Uq   = PrmBuf(5,IK)
            r1xZpE=1.0D0/(Zeta+Eta)
            Rkk=Up*Uq*DSQRT(r1xZpE)
            t1=(Zeta*Px+Eta*Qx)*r1xZpE
            t2=(Zeta*Py+Eta*Qy)*r1xZpE
            t3=(Zeta*Pz+Eta*Qz)*r1xZpE
            WR(1,IB,IK)=t1-Px
            WR(2,IB,IK)=t2-Py
            WR(3,IB,IK)=t3-Pz
            WR(4,IB,IK)=t1-Qx
            WR(5,IB,IK)=t2-Qy
            WR(6,IB,IK)=t3-Qz
            WZ(1,IB,IK)=Zeta*r1xZpE
            WZ(2,IB,IK)=Eta*r1xZpE
            WZ(3,IB,IK)=0.5D0*r1xZpE

            DO I=0,Ltot
               G=1.0D0/(1.0D0+2.0D0*DFLOAT(I))
               R(IB,IK,I+1)=G*Rkk
            ENDDO
         ENDDO
      ENDDO

      RETURN 
      END

