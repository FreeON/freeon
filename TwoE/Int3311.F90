! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (P P|S S) 
! ---------------------------------------------------------- 
   SUBROUTINE Int3311(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, & 
                              OA,LDA,OB,LDB,OC,LDC,OD,LDD) 
      USE DerivedTypes
      USE GlobalScalars
      USE ShellPairStruct
      USE GammaF2
      IMPLICIT REAL(DOUBLE) (A,I,V,W)
      INTEGER        :: LBra,LKet
      REAL(DOUBLE)   :: PrmBufB(5,LBra),PrmBufK(5,LKet)
      TYPE(AtomInfo) :: ACInfo,BDInfo
      REAL(DOUBLE),DIMENSION(:) :: I
      REAL(DOUBLE)  :: Zeta,Eta,r1xZpE,HfxZpE,r1x2E,r1x2Z,ExZpE,ZxZpE,Omega,Up,Uq,Upq
      REAL(DOUBLE)  :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz,Qx,Qy,Qz,Px,Py,Pz,Wx,Wy,Wz
      REAL(DOUBLE)  :: QCx,QCy,QCz,PAx,PAy,PAz,PQx,PQy,PQz,WPx,WPy,WPz,WQx,WQy,WQz   
      REAL(DOUBLE)  :: T,ET,TwoT,InvT,SqInvT
      INTEGER       :: OA,LDA,OB,LDB,LC,LDC,OD,LDD,J,K,L
      I1Bar1=Zero
      I2Bar1=Zero
      I3Bar1=Zero
      I4Bar1=Zero
      I5Bar1=Zero
      I6Bar1=Zero
      I7Bar1=Zero
      I8Bar1=Zero
      I9Bar1=Zero
      I10Bar1=Zero
      Ax=ACInfo%Atm1X
      Ay=ACInfo%Atm1Y
      Az=ACInfo%Atm1Z
      Bx=BDInfo%Atm1X
      By=BDInfo%Atm1Y
      Bz=BDInfo%Atm1Z
      Cx=ACInfo%Atm2X
      Cy=ACInfo%Atm2Y
      Cz=ACInfo%Atm2Z
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
              ET=EXP(-T)
              TwoT=Two*T
              W2=(F2_0(L)+T*(F2_1(L)+T*(F2_2(L)+T*(F2_3(L)+T*F2_4(L)))))
              W1=+3.333333333333333D-01*(TwoT*W2+ET)
              W0=TwoT*W1+ET
              AuxR0=Upq*W0
              AuxR1=Upq*W1
              AuxR2=Upq*W2
            ELSE
              InvT=One/T
              SqInvT=DSQRT(InvT)
              AuxR0=+8.862269254527580D-01*Upq*SqInvT
              SqInvT=SqInvT*InvT
              AuxR1=+4.431134627263790D-01*Upq*SqInvT
              SqInvT=SqInvT*InvT
              AuxR2=+6.646701940895685D-01*Upq*SqInvT
            ENDIF
            V1=AuxR0*PAx
            V2=AuxR1*WPx
            V3=AuxR0*PAy
            V4=AuxR1*WPy
            V5=AuxR0*PAz
            V6=AuxR1*WPz
            V7=-(AuxR1*ExZpE)
            V8=AuxR0+V7
            V9=r1x2Z*V8
            V10=V3+V4
            V11=AuxR1*PAy
            V12=AuxR2*WPy
            V13=V11+V12
            V14=V5+V6
            V15=AuxR1*PAz
            V16=AuxR2*WPz
            V17=V15+V16
            I1Bar1=AuxR0+I1Bar1
            I2Bar1=V1+V2+I2Bar1
            I3Bar1=V3+V4+I3Bar1
            I4Bar1=V5+V6+I4Bar1
            W1=PAx*(V1+V2)+V9
            W2=WPx*(AuxR1*PAx+AuxR2*WPx)+I5Bar1
            I5Bar1=W1+W2
            I6Bar1=PAx*V10+V13*WPx+I6Bar1
            I7Bar1=PAy*V10+V9+V13*WPy+I7Bar1
            I8Bar1=PAx*V14+V17*WPx+I8Bar1
            I9Bar1=PAy*V14+V17*WPy+I9Bar1
            I10Bar1=PAz*V14+V9+V17*WPz+I10Bar1
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! HRR 
      I((OA+1)*LDA+(OB+1)*LDB+(OC+1)*LDC+(OD+1)*LDD)=ABx*I2Bar1+I5Bar1
      I((OA+2)*LDA+(OB+1)*LDB+(OC+1)*LDC+(OD+1)*LDD)=ABx*I3Bar1+I6Bar1
      I((OA+3)*LDA+(OB+1)*LDB+(OC+1)*LDC+(OD+1)*LDD)=ABx*I4Bar1+I8Bar1
      I((OA+1)*LDA+(OB+2)*LDB+(OC+1)*LDC+(OD+1)*LDD)=ABy*I2Bar1+I6Bar1
      I((OA+2)*LDA+(OB+2)*LDB+(OC+1)*LDC+(OD+1)*LDD)=ABy*I3Bar1+I7Bar1
      I((OA+3)*LDA+(OB+2)*LDB+(OC+1)*LDC+(OD+1)*LDD)=ABy*I4Bar1+I9Bar1
      I((OA+1)*LDA+(OB+3)*LDB+(OC+1)*LDC+(OD+1)*LDD)=ABz*I2Bar1+I8Bar1
      I((OA+2)*LDA+(OB+3)*LDB+(OC+1)*LDC+(OD+1)*LDD)=ABz*I3Bar1+I9Bar1
      I((OA+3)*LDA+(OB+3)*LDB+(OC+1)*LDC+(OD+1)*LDD)=I10Bar1+ABz*I4Bar1
   END SUBROUTINE Int3311
