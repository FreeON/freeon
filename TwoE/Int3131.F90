! ---------------------------------------------------------- 
! COMPUTES THE INTEGRAL CLASS (P S|P S) 
! ---------------------------------------------------------- 
   SUBROUTINE Int3131(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, & 
                              OA,LDA,OB,LDB,OC,LDC,OD,LDD,PBC,I) 
      USE DerivedTypes
      USE GlobalScalars
      USE ShellPairStruct
      USE GammaF2
      IMPLICIT REAL(DOUBLE) (A,I,V,W)
      INTEGER        :: LBra,LKet
      REAL(DOUBLE)   :: PrmBufB(5,LBra),PrmBufK(5,LKet)
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
      I1Bar2=0.0d0
      I2Bar2=0.0d0
      I3Bar2=0.0d0
      I4Bar2=0.0d0
      I1Bar3=0.0d0
      I2Bar3=0.0d0
      I3Bar3=0.0d0
      I4Bar3=0.0d0
      I1Bar4=0.0d0
      I2Bar4=0.0d0
      I3Bar4=0.0d0
      I4Bar4=0.0d0
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

!!$         write(*,*) 'Eta',Eta
!!$         write(*,*) 'Uq',Uq

         QCx=Qx-Cx
         QCy=Qy-Cy
         QCz=Qz-Cz
         DO K=1,LBra ! K^2 VRR (M0| loop 
            Zeta=PrmBufB(1,K)
            Px  =PrmBufB(2,K)
            Py  =PrmBufB(3,K)
            Pz  =PrmBufB(4,K)
            Up  =PrmBufB(5,K)

!!$            write(*,*) 'Zeta',Zeta
!!$            write(*,*) 'Up',Up

            r1xZpE=One/(Zeta+Eta)
            Upq=SQRT(r1xZpE)*Up*Uq
            HfxZpE=Half/(Zeta+Eta)
            r1x2E=Half/Eta
            r1x2Z=Half/Zeta
            ExZpE=Eta*r1xZpE
            ZxZpE=Zeta*r1xZpE
            Omega=Eta*Zeta*r1xZpE
            !Wx=(Zeta*Px+Eta*Qx)*r1xZpE
            !Wy=(Zeta*Py+Eta*Qy)*r1xZpE
            !Wz=(Zeta*Pz+Eta*Qz)*r1xZpE
            PAx=Px-Ax
            PAy=Py-Ay
            PAz=Pz-Az
            PQx=Px-Qx
            PQy=Py-Qy
            PQz=Pz-Qz

            !write(*,*) 'PQx',PQx

      ! Need to be improve...
            FPQx = PQx*PBC%InvBoxSh%D(1,1)+PQy*PBC%InvBoxSh%D(1,2)+PQz*PBC%InvBoxSh%D(1,3)
            FPQy = PQy*PBC%InvBoxSh%D(2,2)+PQz*PBC%InvBoxSh%D(2,3)
            FPQz = PQz*PBC%InvBoxSh%D(3,3)
            IF(PBC%AutoW%I(1)==1) FPQx = FPQx-ANINT(ANINT(FPQx*1d9)*1d-9)
            IF(PBC%AutoW%I(2)==1) FPQy = FPQy-ANINT(ANINT(FPQy*1d9)*1d-9)
            IF(PBC%AutoW%I(3)==1) FPQz = FPQz-ANINT(ANINT(FPQz*1d9)*1d-9)
            PQx  = FPQx*PBC%BoxShape%D(1,1)+FPQy*PBC%BoxShape%D(1,2)+FPQz*PBC%BoxShape%D(1,3)
            PQy  = FPQy*PBC%BoxShape%D(2,2)+FPQz*PBC%BoxShape%D(2,3)
            PQz  = FPQz*PBC%BoxShape%D(3,3)

            !write(*,*) 'PQx',PQx

WPx = -Eta*PQx*r1xZpE
WPy = -Eta*PQy*r1xZpE
WPz = -Eta*PQz*r1xZpE
WQx = Zeta*PQx*r1xZpE
WQy = Zeta*PQy*r1xZpE
WQz = Zeta*PQz*r1xZpE
      !
            !WPx=Wx-Px
            !WPy=Wy-Py
            !WPz=Wz-Pz

            !WPx=-0.944862994289462d0
            !write(*,*) 'WPx',WPx


            !WQx=Wx-Qx
            !WQy=Wy-Qy
            !WQz=Wz-Qz

            !WQx=0.944862994289462d0
            !write(*,*) 'WQx',WQx

            T=Omega*(PQx*PQx+PQy*PQy+PQz*PQz)

            !write(*,*) 'T',T

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
            V7=AuxR1*HfxZpE
            V8=V1+V2
            V9=AuxR1*PAx
            V10=AuxR2*WPx
            V11=V10+V9
            V12=V3+V4
            V13=AuxR1*PAy
            V14=AuxR2*WPy
            V15=V13+V14
            V16=V5+V6
            V17=AuxR1*PAz
            V18=AuxR2*WPz
            V19=V17+V18
            I1Bar1=AuxR0+I1Bar1
            I2Bar1=V1+V2+I2Bar1
            I3Bar1=V3+V4+I3Bar1
            I4Bar1=V5+V6+I4Bar1
            I1Bar2=AuxR0*QCx+AuxR1*WQx+I1Bar2
            I2Bar2=V7+QCx*V8+V11*WQx+I2Bar2
            I3Bar2=QCx*V12+V15*WQx+I3Bar2
            I4Bar2=QCx*V16+V19*WQx+I4Bar2
            I1Bar3=AuxR0*QCy+AuxR1*WQy+I1Bar3
            I2Bar3=QCy*V8+V11*WQy+I2Bar3
            I3Bar3=QCy*V12+V7+V15*WQy+I3Bar3
            I4Bar3=QCy*V16+V19*WQy+I4Bar3
            I1Bar4=AuxR0*QCz+AuxR1*WQz+I1Bar4
            I2Bar4=QCz*V8+V11*WQz+I2Bar4
            I3Bar4=QCz*V12+V15*WQz+I3Bar4
            I4Bar4=QCz*V16+V7+V19*WQz+I4Bar4
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! HRR 
      I((OA+0)*LDA+OB*LDB+(OC+0)*LDC+OD*LDD)=I((OA+0)*LDA+OB*LDB+(OC+0)*LDC+OD*LDD)+I2Bar2
      I((OA+1)*LDA+OB*LDB+(OC+0)*LDC+OD*LDD)=I((OA+1)*LDA+OB*LDB+(OC+0)*LDC+OD*LDD)+I3Bar2
      I((OA+2)*LDA+OB*LDB+(OC+0)*LDC+OD*LDD)=I((OA+2)*LDA+OB*LDB+(OC+0)*LDC+OD*LDD)+I4Bar2
      I((OA+0)*LDA+OB*LDB+(OC+1)*LDC+OD*LDD)=I((OA+0)*LDA+OB*LDB+(OC+1)*LDC+OD*LDD)+I2Bar3
      I((OA+1)*LDA+OB*LDB+(OC+1)*LDC+OD*LDD)=I((OA+1)*LDA+OB*LDB+(OC+1)*LDC+OD*LDD)+I3Bar3
      I((OA+2)*LDA+OB*LDB+(OC+1)*LDC+OD*LDD)=I((OA+2)*LDA+OB*LDB+(OC+1)*LDC+OD*LDD)+I4Bar3
      I((OA+0)*LDA+OB*LDB+(OC+2)*LDC+OD*LDD)=I((OA+0)*LDA+OB*LDB+(OC+2)*LDC+OD*LDD)+I2Bar4
      I((OA+1)*LDA+OB*LDB+(OC+2)*LDC+OD*LDD)=I((OA+1)*LDA+OB*LDB+(OC+2)*LDC+OD*LDD)+I3Bar4
      I((OA+2)*LDA+OB*LDB+(OC+2)*LDC+OD*LDD)=I((OA+2)*LDA+OB*LDB+(OC+2)*LDC+OD*LDD)+I4Bar4



!!$if(abs(I2Bar2).GT.1.0d-15)write(*,'(A,E22.15)') 'I2Bar2',I2Bar2
!!$if(abs(I3Bar2).GT.1.0d-15)write(*,'(A,E22.15)') 'I3Bar2',I3Bar2
!!$if(abs(I4Bar2).GT.1.0d-15)write(*,'(A,E22.15)') 'I4Bar2',I4Bar2
!!$if(abs(I2Bar3).GT.1.0d-15)write(*,'(A,E22.15)') 'I2Bar3',I2Bar3
!!$if(abs(I3Bar3).GT.1.0d-15)write(*,'(A,E22.15)') 'I3Bar3',I3Bar3
!!$if(abs(I4Bar3).GT.1.0d-15)write(*,'(A,E22.15)') 'I4Bar3',I4Bar3
!!$if(abs(I2Bar4).GT.1.0d-15)write(*,'(A,E22.15)') 'I2Bar4',I2Bar4
!!$if(abs(I3Bar4).GT.1.0d-15)write(*,'(A,E22.15)') 'I3Bar4',I3Bar4
!!$if(abs(I4Bar4).GT.1.0d-15)write(*,'(A,E22.15)') 'I4Bar4',I4Bar4

!!$write(*,'(A,9I4)') 'Off', &
!!$     & (OA+0)*LDA+OB*LDB+(OC+0)*LDC+OD*LDD, &
!!$     & (OA+1)*LDA+OB*LDB+(OC+0)*LDC+OD*LDD, &
!!$     & (OA+2)*LDA+OB*LDB+(OC+0)*LDC+OD*LDD, &
!!$     & (OA+0)*LDA+OB*LDB+(OC+1)*LDC+OD*LDD, &
!!$     & (OA+1)*LDA+OB*LDB+(OC+1)*LDC+OD*LDD, &
!!$     & (OA+2)*LDA+OB*LDB+(OC+1)*LDC+OD*LDD, &
!!$     & (OA+0)*LDA+OB*LDB+(OC+2)*LDC+OD*LDD, &
!!$     & (OA+1)*LDA+OB*LDB+(OC+2)*LDC+OD*LDD, &
!!$     & (OA+2)*LDA+OB*LDB+(OC+2)*LDC+OD*LDD

   END SUBROUTINE Int3131
