! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (P P|P S) 
! ---------------------------------------------------------- 
   SUBROUTINE Int3331(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, & 
                              OA,LDA,OB,LDB,OC,LDC,OD,LDD,PBC,I) 
      USE DerivedTypes
      USE GlobalScalars
      USE ShellPairStruct
      USE GammaF3
      IMPLICIT REAL(DOUBLE) (A,I,V,W)
      INTEGER        :: LBra,LKet
      REAL(DOUBLE)   :: PrmBufB(7,LBra),PrmBufK(7,LKet)
      TYPE(SmallAtomInfo) :: ACInfo,BDInfo
      TYPE(PBCInfo) :: PBC
      REAL(DOUBLE) :: I(*)
      REAL(DOUBLE)  :: Zeta,Eta,r1xZpE,HfxZpE,r1x2E,r1x2Z,ExZpE,ZxZpE,Omega,Up,Uq,Upq
      REAL(DOUBLE)  :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz,Qx,Qy,Qz,Px,Py,Pz,Wx,Wy,Wz
      REAL(DOUBLE)  :: QCx,QCy,QCz,PAx,PAy,PAz,PQx,PQy,PQz,WPx,WPy,WPz,WQx,WQy,WQz   
      REAL(DOUBLE)  :: T,ET,TwoT,InvT,SqInvT,ABx,ABy,ABz,CDx,CDy,CDz
      INTEGER       :: OA,LDA,OB,LDB,OC,LDC,OD,LDD,J,K,L
      REAL(DOUBLE)  :: FPQx,FPQy,FPQz
      I1Bar1=0.0d0
      I2Bar1=0.0d0
      I3Bar1=0.0d0
      I4Bar1=0.0d0
      I5Bar1=0.0d0
      I6Bar1=0.0d0
      I7Bar1=0.0d0
      I8Bar1=0.0d0
      I9Bar1=0.0d0
      I10Bar1=0.0d0
      I1Bar2=0.0d0
      I2Bar2=0.0d0
      I3Bar2=0.0d0
      I4Bar2=0.0d0
      I5Bar2=0.0d0
      I6Bar2=0.0d0
      I7Bar2=0.0d0
      I8Bar2=0.0d0
      I9Bar2=0.0d0
      I10Bar2=0.0d0
      I1Bar3=0.0d0
      I2Bar3=0.0d0
      I3Bar3=0.0d0
      I4Bar3=0.0d0
      I5Bar3=0.0d0
      I6Bar3=0.0d0
      I7Bar3=0.0d0
      I8Bar3=0.0d0
      I9Bar3=0.0d0
      I10Bar3=0.0d0
      I1Bar4=0.0d0
      I2Bar4=0.0d0
      I3Bar4=0.0d0
      I4Bar4=0.0d0
      I5Bar4=0.0d0
      I6Bar4=0.0d0
      I7Bar4=0.0d0
      I8Bar4=0.0d0
      I9Bar4=0.0d0
      I10Bar4=0.0d0
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
            PAx=Px-Ax
            PAy=Py-Ay
            PAz=Pz-Az
            PQx=Px-Qx
            PQy=Py-Qy
            PQz=Pz-Qz
      ! Need to be improve...
            FPQx = PQx*PBC%InvBoxSh%D(1,1)+PQy*PBC%InvBoxSh%D(1,2)+PQz*PBC%InvBoxSh%D(1,3)
            FPQy = PQy*PBC%InvBoxSh%D(2,2)+PQz*PBC%InvBoxSh%D(2,3)
            FPQz = PQz*PBC%InvBoxSh%D(3,3)
            IF(PBC%AutoW%I(1)==1) FPQx = FPQx-ANINT(FPQx+sign(1d-9,FPQx))
            !IF(PBC%AutoW%I(1)==1) FPQx = FPQx-ANINT(ANINT(FPQx*1d9)*1d-9)
            IF(PBC%AutoW%I(2)==1) FPQy = FPQy-ANINT(ANINT(FPQy*1d9)*1d-9)
            IF(PBC%AutoW%I(3)==1) FPQz = FPQz-ANINT(ANINT(FPQz*1d9)*1d-9)
            PQx  = FPQx*PBC%BoxShape%D(1,1)+FPQy*PBC%BoxShape%D(1,2)+FPQz*PBC%BoxShape%D(1,3)
            PQy  = FPQy*PBC%BoxShape%D(2,2)+FPQz*PBC%BoxShape%D(2,3)
            PQz  = FPQz*PBC%BoxShape%D(3,3)
      !
            WPx = -Eta*PQx*r1xZpE
            WPy = -Eta*PQy*r1xZpE
            WPz = -Eta*PQz*r1xZpE
            WQx = Zeta*PQx*r1xZpE
            WQy = Zeta*PQy*r1xZpE
            WQz = Zeta*PQz*r1xZpE
            T=Omega*(PQx*PQx+PQy*PQy+PQz*PQz)
            IF(T<Gamma_Switch)THEN
              L=AINT(T*Gamma_Grid)
              ET=EXP(-T)
              TwoT=Two*T
              W3=(F3_0(L)+T*(F3_1(L)+T*(F3_2(L)+T*(F3_3(L)+T*F3_4(L)))))
              W2=+2.000000000000000D-01*(TwoT*W3+ET)
              W1=+3.333333333333333D-01*(TwoT*W2+ET)
              W0=TwoT*W1+ET
              AuxR0=Upq*W0
              AuxR1=Upq*W1
              AuxR2=Upq*W2
              AuxR3=Upq*W3
            ELSE
              InvT=One/T
              SqInvT=DSQRT(InvT)
              AuxR0=+8.862269254527580D-01*Upq*SqInvT
              SqInvT=SqInvT*InvT
              AuxR1=+4.431134627263790D-01*Upq*SqInvT
              SqInvT=SqInvT*InvT
              AuxR2=+6.646701940895685D-01*Upq*SqInvT
              SqInvT=SqInvT*InvT
              AuxR3=+1.661675485223921D+00*Upq*SqInvT
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
            V10=V1+V2
            V11=AuxR1*PAx
            V12=AuxR2*WPx
            V13=V11+V12
            V14=V3+V4
            V15=PAx*V14
            V16=AuxR1*PAy
            V17=AuxR2*WPy
            V18=V16+V17
            V19=V18*WPx
            V20=V5+V6
            V21=PAx*V20
            V22=AuxR1*PAz
            V23=AuxR2*WPz
            V24=V22+V23
            V25=V24*WPx
            V26=PAy*V20
            V27=V24*WPy
            V28=AuxR0*QCx
            V29=AuxR1*WQx
            V30=AuxR1*HfxZpE
            V31=QCx*V10
            V32=V13*WQx
            V33=QCx*V14
            V34=V18*WQx
            V35=QCx*V20
            V36=V24*WQx
            V37=HfxZpE*V13
            V38=AuxR2*HfxZpE
            V39=AuxR2*PAx
            V40=AuxR3*WPx
            V41=V39+V40
            V42=AuxR1*QCx
            V43=AuxR2*WQx
            V44=V42+V43
            V45=-(ExZpE*V44)
            V46=V28+V29+V45
            V47=r1x2Z*V46
            V48=HfxZpE*V18
            V49=V15+V19
            V50=PAx*V18
            V51=AuxR2*PAy
            V52=AuxR3*WPy
            V53=V51+V52
            V54=V53*WPx
            V55=V50+V54
            V56=HfxZpE*V24
            V57=V21+V25
            V58=PAx*V24
            V59=AuxR2*PAz
            V60=AuxR3*WPz
            V61=V59+V60
            V62=V61*WPx
            V63=V58+V62
            V64=V26+V27
            V65=PAy*V24
            V66=V61*WPy
            V67=V65+V66
            V68=AuxR0*QCy
            V69=AuxR1*WQy
            V70=QCy*V10
            V71=V13*WQy
            V72=QCy*V14
            V73=V18*WQy
            V74=QCy*V20
            V75=V24*WQy
            V76=AuxR1*QCy
            V77=AuxR2*WQy
            V78=V76+V77
            V79=-(ExZpE*V78)
            V80=V68+V69+V79
            V81=r1x2Z*V80
            V82=AuxR0*QCz
            V83=AuxR1*WQz
            V84=QCz*V10
            V85=V13*WQz
            V86=QCz*V14
            V87=V18*WQz
            V88=QCz*V20
            V89=V24*WQz
            V90=AuxR1*QCz
            V91=AuxR2*WQz
            V92=V90+V91
            V93=-(ExZpE*V92)
            V94=V82+V83+V93
            V95=r1x2Z*V94
            I1Bar1=AuxR0+I1Bar1
            I2Bar1=V1+V2+I2Bar1
            I3Bar1=V3+V4+I3Bar1
            I4Bar1=V5+V6+I4Bar1
            I5Bar1=PAx*V10+V9+V13*WPx+I5Bar1
            I6Bar1=V15+V19+I6Bar1
            I7Bar1=PAy*V14+V9+V18*WPy+I7Bar1
            I8Bar1=V21+V25+I8Bar1
            I9Bar1=V26+V27+I9Bar1
            I10Bar1=PAz*V20+V9+V24*WPz+I10Bar1
            I1Bar2=V28+V29+I1Bar2
            I2Bar2=V30+V31+V32+I2Bar2
            I3Bar2=V33+V34+I3Bar2
            I4Bar2=V35+V36+I4Bar2
            W1=PAx*(V30+V31+V32)+V37
            W2=V47+WPx*(QCx*V13+V38+V41*WQx)+I5Bar2
            I5Bar2=W1+W2
            I6Bar2=V48+QCx*V49+V55*WQx+I6Bar2
            W1=PAy*(V33+V34)+V47
            W2=WPy*(QCx*V18+V53*WQx)+I7Bar2
            I7Bar2=W1+W2
            I8Bar2=V56+QCx*V57+V63*WQx+I8Bar2
            I9Bar2=QCx*V64+V67*WQx+I9Bar2
            W1=PAz*(V35+V36)+V47
            W2=WPz*(QCx*V24+V61*WQx)+I10Bar2
            I10Bar2=W1+W2
            I1Bar3=V68+V69+I1Bar3
            I2Bar3=V70+V71+I2Bar3
            I3Bar3=V30+V72+V73+I3Bar3
            I4Bar3=V74+V75+I4Bar3
            W1=PAx*(V70+V71)+V81
            W2=WPx*(QCy*V13+V41*WQy)+I5Bar3
            I5Bar3=W1+W2
            I6Bar3=V37+QCy*V49+V55*WQy+I6Bar3
            W1=V48+PAy*(V30+V72+V73)
            W2=V81+WPy*(QCy*V18+V38+V53*WQy)+I7Bar3
            I7Bar3=W1+W2
            I8Bar3=QCy*V57+V63*WQy+I8Bar3
            I9Bar3=V56+QCy*V64+V67*WQy+I9Bar3
            W1=PAz*(V74+V75)+V81
            W2=WPz*(QCy*V24+V61*WQy)+I10Bar3
            I10Bar3=W1+W2
            I1Bar4=V82+V83+I1Bar4
            I2Bar4=V84+V85+I2Bar4
            I3Bar4=V86+V87+I3Bar4
            I4Bar4=V30+V88+V89+I4Bar4
            W1=PAx*(V84+V85)+V95
            W2=WPx*(QCz*V13+V41*WQz)+I5Bar4
            I5Bar4=W1+W2
            I6Bar4=QCz*V49+V55*WQz+I6Bar4
            W1=PAy*(V86+V87)+V95
            W2=WPy*(QCz*V18+V53*WQz)+I7Bar4
            I7Bar4=W1+W2
            I8Bar4=V37+QCz*V57+V63*WQz+I8Bar4
            I9Bar4=V48+QCz*V64+V67*WQz+I9Bar4
            W1=V56+PAz*(V30+V88+V89)
            W2=V95+WPz*(QCz*V24+V38+V61*WQz)+I10Bar4
            I10Bar4=W1+W2
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! HRR 
      I((OA+0)*LDA+(OB+0)*LDB+(OC+0)*LDC+OD*LDD)=I((OA+0)*LDA+(OB+0)*LDB+(OC+0)*LDC+OD*LDD)+ABx*I2Bar2+I5Bar2
      I((OA+1)*LDA+(OB+0)*LDB+(OC+0)*LDC+OD*LDD)=I((OA+1)*LDA+(OB+0)*LDB+(OC+0)*LDC+OD*LDD)+ABx*I3Bar2+I6Bar2
      I((OA+2)*LDA+(OB+0)*LDB+(OC+0)*LDC+OD*LDD)=I((OA+2)*LDA+(OB+0)*LDB+(OC+0)*LDC+OD*LDD)+ABx*I4Bar2+I8Bar2
      I((OA+0)*LDA+(OB+1)*LDB+(OC+0)*LDC+OD*LDD)=I((OA+0)*LDA+(OB+1)*LDB+(OC+0)*LDC+OD*LDD)+ABy*I2Bar2+I6Bar2
      I((OA+1)*LDA+(OB+1)*LDB+(OC+0)*LDC+OD*LDD)=I((OA+1)*LDA+(OB+1)*LDB+(OC+0)*LDC+OD*LDD)+ABy*I3Bar2+I7Bar2
      I((OA+2)*LDA+(OB+1)*LDB+(OC+0)*LDC+OD*LDD)=I((OA+2)*LDA+(OB+1)*LDB+(OC+0)*LDC+OD*LDD)+ABy*I4Bar2+I9Bar2
      I((OA+0)*LDA+(OB+2)*LDB+(OC+0)*LDC+OD*LDD)=I((OA+0)*LDA+(OB+2)*LDB+(OC+0)*LDC+OD*LDD)+ABz*I2Bar2+I8Bar2
      I((OA+1)*LDA+(OB+2)*LDB+(OC+0)*LDC+OD*LDD)=I((OA+1)*LDA+(OB+2)*LDB+(OC+0)*LDC+OD*LDD)+ABz*I3Bar2+I9Bar2
      I((OA+2)*LDA+(OB+2)*LDB+(OC+0)*LDC+OD*LDD)=I((OA+2)*LDA+(OB+2)*LDB+(OC+0)*LDC+OD*LDD)+I10Bar2+ABz*I4Bar2
      I((OA+0)*LDA+(OB+0)*LDB+(OC+1)*LDC+OD*LDD)=I((OA+0)*LDA+(OB+0)*LDB+(OC+1)*LDC+OD*LDD)+ABx*I2Bar3+I5Bar3
      I((OA+1)*LDA+(OB+0)*LDB+(OC+1)*LDC+OD*LDD)=I((OA+1)*LDA+(OB+0)*LDB+(OC+1)*LDC+OD*LDD)+ABx*I3Bar3+I6Bar3
      I((OA+2)*LDA+(OB+0)*LDB+(OC+1)*LDC+OD*LDD)=I((OA+2)*LDA+(OB+0)*LDB+(OC+1)*LDC+OD*LDD)+ABx*I4Bar3+I8Bar3
      I((OA+0)*LDA+(OB+1)*LDB+(OC+1)*LDC+OD*LDD)=I((OA+0)*LDA+(OB+1)*LDB+(OC+1)*LDC+OD*LDD)+ABy*I2Bar3+I6Bar3
      I((OA+1)*LDA+(OB+1)*LDB+(OC+1)*LDC+OD*LDD)=I((OA+1)*LDA+(OB+1)*LDB+(OC+1)*LDC+OD*LDD)+ABy*I3Bar3+I7Bar3
      I((OA+2)*LDA+(OB+1)*LDB+(OC+1)*LDC+OD*LDD)=I((OA+2)*LDA+(OB+1)*LDB+(OC+1)*LDC+OD*LDD)+ABy*I4Bar3+I9Bar3
      I((OA+0)*LDA+(OB+2)*LDB+(OC+1)*LDC+OD*LDD)=I((OA+0)*LDA+(OB+2)*LDB+(OC+1)*LDC+OD*LDD)+ABz*I2Bar3+I8Bar3
      I((OA+1)*LDA+(OB+2)*LDB+(OC+1)*LDC+OD*LDD)=I((OA+1)*LDA+(OB+2)*LDB+(OC+1)*LDC+OD*LDD)+ABz*I3Bar3+I9Bar3
      I((OA+2)*LDA+(OB+2)*LDB+(OC+1)*LDC+OD*LDD)=I((OA+2)*LDA+(OB+2)*LDB+(OC+1)*LDC+OD*LDD)+I10Bar3+ABz*I4Bar3
      I((OA+0)*LDA+(OB+0)*LDB+(OC+2)*LDC+OD*LDD)=I((OA+0)*LDA+(OB+0)*LDB+(OC+2)*LDC+OD*LDD)+ABx*I2Bar4+I5Bar4
      I((OA+1)*LDA+(OB+0)*LDB+(OC+2)*LDC+OD*LDD)=I((OA+1)*LDA+(OB+0)*LDB+(OC+2)*LDC+OD*LDD)+ABx*I3Bar4+I6Bar4
      I((OA+2)*LDA+(OB+0)*LDB+(OC+2)*LDC+OD*LDD)=I((OA+2)*LDA+(OB+0)*LDB+(OC+2)*LDC+OD*LDD)+ABx*I4Bar4+I8Bar4
      I((OA+0)*LDA+(OB+1)*LDB+(OC+2)*LDC+OD*LDD)=I((OA+0)*LDA+(OB+1)*LDB+(OC+2)*LDC+OD*LDD)+ABy*I2Bar4+I6Bar4
      I((OA+1)*LDA+(OB+1)*LDB+(OC+2)*LDC+OD*LDD)=I((OA+1)*LDA+(OB+1)*LDB+(OC+2)*LDC+OD*LDD)+ABy*I3Bar4+I7Bar4
      I((OA+2)*LDA+(OB+1)*LDB+(OC+2)*LDC+OD*LDD)=I((OA+2)*LDA+(OB+1)*LDB+(OC+2)*LDC+OD*LDD)+ABy*I4Bar4+I9Bar4
      I((OA+0)*LDA+(OB+2)*LDB+(OC+2)*LDC+OD*LDD)=I((OA+0)*LDA+(OB+2)*LDB+(OC+2)*LDC+OD*LDD)+ABz*I2Bar4+I8Bar4
      I((OA+1)*LDA+(OB+2)*LDB+(OC+2)*LDC+OD*LDD)=I((OA+1)*LDA+(OB+2)*LDB+(OC+2)*LDC+OD*LDD)+ABz*I3Bar4+I9Bar4
      I((OA+2)*LDA+(OB+2)*LDB+(OC+2)*LDC+OD*LDD)=I((OA+2)*LDA+(OB+2)*LDB+(OC+2)*LDC+OD*LDD)+I10Bar4+ABz*I4Bar4
   END SUBROUTINE Int3331
