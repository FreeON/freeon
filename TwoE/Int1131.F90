! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (S S|P S) 
! ---------------------------------------------------------- 
   SUBROUTINE Int1131(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo,I, & 
                              OA,LDA,OB,LDB,OC,LDC,OD,LDD) 
      USE DerivedTypes
      USE GlobalScalars
      USE ONX2DataType
      USE GammaF0
      USE GammaF1
      IMPLICIT REAL(DOUBLE) (A,I,V,W)
      INTEGER        :: LBra,LKet
      REAL(DOUBLE)   :: PrmBufB(5,LBra),PrmBufK(5,LKet)
      TYPE(AtomInfo) :: ACInfo,BDInfo
      REAL(DOUBLE) :: I(*)
      REAL(DOUBLE)  :: Zeta,Eta,r1xZpE,HfxZpE,r1x2E,r1x2Z,ExZpE,ZxZpE,Omega,Up,Uq,Upq
      REAL(DOUBLE)  :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz,Qx,Qy,Qz,Px,Py,Pz,Wx,Wy,Wz
      REAL(DOUBLE)  :: QCx,QCy,QCz,PAx,PAy,PAz,PQx,PQy,PQz,WPx,WPy,WPz,WQx,WQy,WQz   
      REAL(DOUBLE)  :: T,ET,TwoT,InvT,SqInvT,ABx,ABy,ABz,CDx,CDy,CDz
      INTEGER       :: OA,LDA,OB,LDB,OC,LDC,OD,LDD,J,K,L
      I1Bar1=Zero
      I1Bar2=Zero
      I1Bar3=Zero
      I1Bar4=Zero
      Ax=ACInfo%Atm1X
      Ay=ACInfo%Atm1Y
      Az=ACInfo%Atm1Z
      Bx=ACInfo%Atm2X
      By=ACInfo%Atm2Y
      Bz=ACInfo%Atm2Z
      Cx=BDInfo%Atm1X
      Cy=BDInfo%Atm1Y
      Cz=BDInfo%Atm1Z
      Dx=BDInfo%Atm2X
      Dy=BDInfo%Atm2Y
      Dz=BDInfo%Atm2Z
      ABx=Ax-Bx
      ABy=Ay-By
      ABz=Az-Bz
      CDx=Cx-Dx
      CDy=Cy-Dy
      CDz=Cz-Dz
      DO J=1,LKet ! K^2 VRR |N0) loop 
         Eta=PrmBufK(1,J)
         Qx =PrmBufK(2,J)
         Qy =PrmBufK(3,J)
         Qz =PrmBufK(4,J)
         Uq =PrmBufK(5,J)
         QCx=Qx-Cx
         QCy=Qy-Cy
         QCz=Qz-Cz
         DO K=1,LBra ! K^2 VRR (M0| loop 
            Zeta=PrmBufB(1,K)
            Px  =PrmBufB(2,K)
            Py  =PrmBufB(3,K)
            Pz  =PrmBufB(4,K)
            Up  =PrmBufB(5,K)
            r1xZpE=One/(Zeta+Eta)
            Upq=SQRT(r1xZpE)*Up*Uq
            HfxZpE=Half/(Zeta+Eta)
            r1x2E=Half/Eta
            r1x2Z=Half/Zeta
            ExZpE=Eta*r1xZpE
            ZxZpE=Zeta*r1xZpE
            Omega=Eta*Zeta*r1xZpE
            Wx=(Zeta*Px+Eta*Qx)*r1xZpE
            Wy=(Zeta*Py+Eta*Qy)*r1xZpE
            Wz=(Zeta*Pz+Eta*Qz)*r1xZpE
            PAx=Px-Ax
            PAy=Py-Ay
            PAz=Pz-Az
            PQx=Px-Qx
            PQy=Py-Qy
            PQz=Pz-Qz
            WPx=Wx-Px
            WPy=Wy-Py
            WPz=Wz-Pz
            WQx=Wx-Qx
            WQy=Wy-Qy
            WQz=Wz-Qz
            T=Omega*(PQx*PQx+PQy*PQy+PQz*PQz)
            IF(T<Gamma_Switch)THEN
              L=AINT(T*Gamma_Grid)
              AuxR0=Upq*(F0_0(L)+T*(F0_1(L)+T*(F0_2(L)+T*(F0_3(L)+T*F0_4(L)))))
              AuxR1=Upq*(F1_0(L)+T*(F1_1(L)+T*(F1_2(L)+T*(F1_3(L)+T*F1_4(L)))))
            ELSE
              InvT=One/T
              SqInvT=DSQRT(InvT)
              AuxR0=+8.862269254527580D-01*Upq*SqInvT
              SqInvT=SqInvT*InvT
              AuxR1=+4.431134627263790D-01*Upq*SqInvT
            ENDIF
            I1Bar1=AuxR0+I1Bar1
            I1Bar2=AuxR0*QCx+AuxR1*WQx+I1Bar2
            I1Bar3=AuxR0*QCy+AuxR1*WQy+I1Bar3
            I1Bar4=AuxR0*QCz+AuxR1*WQz+I1Bar4
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! HRR 
      I((OA+1)*LDA+(OB+1)*LDB+(OC+1)*LDC+(OD+1)*LDD)=I1Bar2
      I((OA+1)*LDA+(OB+1)*LDB+(OC+2)*LDC+(OD+1)*LDD)=I1Bar3
      I((OA+1)*LDA+(OB+1)*LDB+(OC+3)*LDC+(OD+1)*LDD)=I1Bar4
   END SUBROUTINE Int1131
