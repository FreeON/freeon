MODULE DEIntegrals
  USE DerivedTypes
  USE Indexing
  USE GlobalScalars
  USE ShellPairStruct
      USE GammaF0

  INTEGER, PARAMETER :: MxEll=1 ! Just through p functions for now, thanks very much!!
  INTEGER, PARAMETER :: LnEll=(2*MxEll+1)*(2*MxEll+2)*(2*MxEll+3)/6 
  INTEGER, PARAMETER :: NGrid=61
  INTEGER, DIMENSION(1:LnEll) :: XDex,YDex,ZDex
  REAL(DOUBLE), PARAMETER :: DECut=3D0
  REAL(DOUBLE),DIMENSION(1:NGrid) :: DEPointsSq,DEJacobian,Temp
  INTEGER, DIMENSION(1:NGrid) :: IPoint
  REAL(DOUBLE),DIMENSION(1:NGrid,0:2*MxEll,0:2*MxEll,0:2*MxEll) :: VRRx,VRRy,VRRz
  REAL(DOUBLE),DIMENSION(1:NGrid,0:2*MxEll,0:MxEll,0:2*MxEll,0:MxEll) :: HRRx,HRRy,HRRz
CONTAINS
  !
  SUBROUTINE DESetUp
    INTEGER      :: IGrid,L,M,N,LMN,I,NGrid2
    REAL(DOUBLE) :: h,X
    !---------------------------------------------------------------------
    ! Set up the Double Exponential grid mapping [0,1] -> (-Inf,Inf)
    NGrid2=FLOOR(Half*DBLE(NGrid))
    h=DECut/DBLE(NGrid2)
    IGrid=0
    DO I=NGrid2,1,-1
       IGrid=IGrid+1
       X=h*DBLE(I)
       DEPointsSq(IGrid)=(Half*(TANH(Half*Pi*SINH(X))+One))**2
       DEJacobian(IGrid)=h*(Pi*Cosh(x)/Cosh((Pi*Sinh(x))/2D0)**2)/4D0
       IGrid=IGrid+1
       X=-h*DBLE(I)
       DEPointsSq(IGrid)=(Half*(TANH(Half*Pi*SINH(X))+One))**2
       DEJacobian(IGrid)=h*(Pi*Cosh(x)/Cosh((Pi*Sinh(x))/2D0)**2)/4D0
    ENDDO
    X=Zero
    DEPointsSq(NGrid)=(Half*(TANH(Half*Pi*SINH(X))+One))**2
    DEJacobian(NGrid)=h*(Pi*Cosh(x)/Cosh((Pi*Sinh(x))/2D0)**2)/4D0
!
!    DO I=1,NGrid
!       WRITE(*,*)I,DEPointsSq(I),DEJacobian(I)
!    ENDDO

    ! Set up xyz indexing
    DO L=0,2*MxEll
       DO M=0,2*MxEll-L
          DO N=0,2*MxEll-L-M
             LMN=LMNDex(L,M,N)
             XDex(LMN)=L
             YDex(LMN)=M
             ZDex(LMN)=N
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE DESetUp
  !
  SUBROUTINE DEInts(PrmBufB,LBra,PrmBufK,LKet,ACInfo,BDInfo, &
                    OA,EllA,LDA,OB,EllB,LDB,OC,EllC,LDC,OD,EllD,LDD,PBC,INTGRL) 

      INTEGER             :: OffSet,AOFF,BOFF,COFF,LBra,LKet,EllBra,EllKet
      INTEGER             :: OA,EllA,LDA,OB,EllB,LDB,OC,EllC,LDC,OD,EllD,LDD
      INTEGER             :: IGrid,I,J,K,L,M,Ix,Jx,Kx,Lx,Iy,Jy,Ky,Ly,Iz,Jz,Kz,Lz
      INTEGER             :: Mp1,Ip1,Im1,Jm1,Km1,Kp1,Lm1,II,KK
      INTEGER             :: BegA,EndA,BegB,EndB,BegC,EndC,BegD,EndD
      REAL(DOUBLE)        :: PrmBufB(7,LBra),PrmBufK(7,LKet)
      REAL(DOUBLE)        :: INTGRL(*)  
      REAL(DOUBLE)        :: Zeta,Eta,r1xZpE,HfxZpE,r1x2E,r1x2Z,ExZpE,ZxZpE,Omega
      REAL(DOUBLE)        :: U2,UThrd,Up,Uq,Upq
      REAL(DOUBLE)        :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz,Qx,Qy,Qz,Px,Py,Pz
      REAL(DOUBLE)        :: QCx,QCy,QCz,PAx,PAy,PAz,PQx,PQy,PQz,WPx,WPy,WPz,WQx,WQy,WQz   
      REAL(DOUBLE)        :: T,ET,TwoT,InvT,SqInvT,ABx,ABy,ABz,CDx,CDy,CDz
      REAL(DOUBLE)        :: FPQx,FPQy,FPQz,PQ2X,PQ2Y,PQ2Z,Sum,AuxR0
      TYPE(SmallAtomInfo) :: ACInfo,BDInfo
      TYPE(PBCInfo)       :: PBC

      !------------------------------------------------------------------------------------------

      BegA=LBegin(EllA)
      EndA=LEnd(EllA)
      BegB=LBegin(EllB)
      EndB=LEnd(EllB)
      BegC=LBegin(EllC)
      EndC=LEnd(EllC)
      BegD=LBegin(EllD)
      EndD=LEnd(EllD)

      EllBra=EllA+EllB
      EllKet=EllC+EllD

      Ax=ACInfo%Atm1X
      Ay=ACInfo%Atm1Y
      Az=ACInfo%Atm1Z

      Cx=BDInfo%Atm1X
      Cy=BDInfo%Atm1Y
      Cz=BDInfo%Atm1Z

      Bx=ACInfo%Atm2X
      By=ACInfo%Atm2Y
      Bz=ACInfo%Atm2Z

      Dx=BDInfo%Atm2X
      Dy=BDInfo%Atm2Y
      Dz=BDInfo%Atm2Z

      ABx=Ax-Bx
      ABy=Ay-By
      ABz=Az-Bz

      CDx=Cx-Dx
      CDy=Cy-Dy
      CDz=Cz-Dz

      DO I=0,EllBra
         DO K=0,EllKet
            DO IGrid=1,NGrid
               HRRx(IGrid,I,0,K,0)=Zero
               HRRy(IGrid,I,0,K,0)=Zero
               HRRz(IGrid,I,0,K,0)=Zero
            ENDDO
         ENDDO
      ENDDO

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

            PAx=Px-Ax
            PAy=Py-Ay
            PAz=Pz-Az

            r1xZpE=One/(Zeta+Eta)
            Upq=SQRT(r1xZpE)*Up*Uq
            UThrd=Upq**(One/Three)

            HfxZpE=Half/(Zeta+Eta)
            R1x2E=Half/Eta
            R1x2Z=Half/Zeta

            ExZpE=Eta*r1xZpE
            ZxZpE=Zeta*r1xZpE

            Omega=Eta*Zeta*r1xZpE

            PQx=Px-Qx
            PQy=Py-Qy
            PQz=Pz-Qz

            FPQx = PQx*PBC%InvBoxSh%D(1,1)+PQy*PBC%InvBoxSh%D(1,2)+PQz*PBC%InvBoxSh%D(1,3)
            FPQy = PQy*PBC%InvBoxSh%D(2,2)+PQz*PBC%InvBoxSh%D(2,3)
            FPQz = PQz*PBC%InvBoxSh%D(3,3)
            IF(PBC%AutoW%I(1)==1) FPQx = FPQx-ANINT(ANINT(FPQx*1d9)*1d-9)
            IF(PBC%AutoW%I(2)==1) FPQy = FPQy-ANINT(ANINT(FPQy*1d9)*1d-9)
            IF(PBC%AutoW%I(3)==1) FPQz = FPQz-ANINT(ANINT(FPQz*1d9)*1d-9)

            PQx  = FPQx*PBC%BoxShape%D(1,1)+FPQy*PBC%BoxShape%D(1,2)+FPQz*PBC%BoxShape%D(1,3)
            PQy  = FPQy*PBC%BoxShape%D(2,2)+FPQz*PBC%BoxShape%D(2,3)
            PQz  = FPQz*PBC%BoxShape%D(3,3)

            WPx = -Eta*PQx*r1xZpE
            WPy = -Eta*PQy*r1xZpE
            WPz = -Eta*PQz*r1xZpE

            WQx = Zeta*PQx*r1xZpE
            WQy = Zeta*PQy*r1xZpE
            WQz = Zeta*PQz*r1xZpE

            PQ2x = PQx*PQx
            PQ2y = PQy*PQy
            PQ2z = PQz*PQz

            Sum=0
            DO IGrid=1,NGrid
               U2=DEPointsSq(IGrid)
               VRRx(IGrid,0,0,0)=UThrd*EXP(-Omega*U2*PQ2x)
               VRRy(IGrid,0,0,0)=UThrd*EXP(-Omega*U2*PQ2y)
               VRRz(IGrid,0,0,0)=UThrd*EXP(-Omega*U2*PQ2z)
               Sum=Sum+VRRx(IGrid,0,0,0)*VRRy(IGrid,0,0,0)*VRRz(IGrid,0,0,0)*DEJacobian(IGrid)
               ! +Upq*EXP(-Omega*U2*(PQ2x+PQ2y+PQ2z))*DEJacobian(IGrid)
!               WRITE(*,*)Upq*EXP(-Omega*U2*(PQ2x+PQ2y+PQ2z))*DEJacobian(IGrid)
            ENDDO

            T=Omega*(PQx*PQx+PQy*PQy+PQz*PQz)
            IF(T<Gamma_Switch)THEN
              L=AINT(T*Gamma_Grid)
              AuxR0=Upq*(F0_0(L)+T*(F0_1(L)+T*(F0_2(L)+T*(F0_3(L)+T*F0_4(L)))))
            ELSE
              InvT=One/T
              SqInvT=DSQRT(InvT)
              AuxR0=+8.862269254527580D-01*Upq*SqInvT
            ENDIF

            IF(ABS(Sum-AuxR0)>1D-12)THEN
               WRITE(*,*)' DIFF = ',ABS(Sum-AuxR0)
               WRITE(*,*)' Sum = ',Sum
               WRITE(*,*)' AuxR0 = ',AuxR0
               WRITE(*,*)' Upq = ',Upq
               WRITE(*,*)' T = ',T
               STOP
            ENDIF

            DO M=1,MAX(EllBra,EllKet) ! ??
               WRITE(*,*)'1 M = ',M
               DO IGrid=1,NGrid
                  U2=DEPointsSq(IGrid)
                  VRRx(IGrid,M,0,0)=U2*VRRx(IGrid,M-1,0,0)
                  VRRy(IGrid,M,0,0)=U2*VRRy(IGrid,M-1,0,0)
                  VRRz(IGrid,M,0,0)=U2*VRRz(IGrid,M-1,0,0)
               ENDDO
            ENDDO

            ! The p functions first
            DO M=0,EllBra-1
               WRITE(*,*)'2 M = ',M
               Mp1=M+1
               DO IGrid=1,NGrid               
                  VRRx(IGrid,M,1,0)=PAx*VRRx(IGrid,M,0,0)+WPx*VRRx(IGrid,Mp1,0,0)
                  VRRy(IGrid,M,1,0)=PAy*VRRy(IGrid,M,0,0)+WPy*VRRy(IGrid,Mp1,0,0)
                  VRRz(IGrid,M,1,0)=PAz*VRRz(IGrid,M,0,0)+WPz*VRRz(IGrid,Mp1,0,0)
               ENDDO
            ENDDO
            DO M=0,EllKet-1
               Mp1=M+1
               DO IGrid=1,NGrid               
                  VRRx(IGrid,M,0,1)=QCx*VRRx(IGrid,M,0,0)+WQx*VRRx(IGrid,Mp1,0,0)
                  VRRy(IGrid,M,0,1)=QCy*VRRy(IGrid,M,0,0)+WQy*VRRy(IGrid,Mp1,0,0)
                  VRRz(IGrid,M,0,1)=QCz*VRRz(IGrid,M,0,0)+WQz*VRRz(IGrid,Mp1,0,0)
               ENDDO
            ENDDO
            ! Then the d functions
            DO M=0,EllBra-2
               Mp1=M+1
               DO IGrid=1,NGrid               
                  VRRx(IGrid,M,2,0)=PAx*VRRx(IGrid,M,1,0)+WPx*VRRx(IGrid,Mp1,1,0)+R1x2Z*(VRRx(IGrid,M,0,0)-ExZpE*VRRx(IGrid,Mp1,0,0))
                  VRRy(IGrid,M,2,0)=PAy*VRRy(IGrid,M,1,0)+WPy*VRRy(IGrid,Mp1,1,0)+R1x2Z*(VRRy(IGrid,M,0,0)-ExZpE*VRRy(IGrid,Mp1,0,0))
                  VRRz(IGrid,M,2,0)=PAz*VRRz(IGrid,M,1,0)+WPz*VRRz(IGrid,Mp1,1,0)+R1x2Z*(VRRz(IGrid,M,0,0)-ExZpE*VRRz(IGrid,Mp1,0,0))
               ENDDO
            ENDDO
            DO M=0,EllKet-2
               Mp1=M+1
               DO IGrid=1,NGrid               
                  VRRx(IGrid,M,0,2)=QCx*VRRx(IGrid,M,0,1)+WQx*VRRx(IGrid,Mp1,0,1)+R1x2E*(VRRx(IGrid,M,0,0)-ZxZpE*VRRx(IGrid,Mp1,0,0))
                  VRRy(IGrid,M,0,2)=QCy*VRRy(IGrid,M,0,1)+WQy*VRRy(IGrid,Mp1,0,1)+R1x2E*(VRRy(IGrid,M,0,0)-ZxZpE*VRRy(IGrid,Mp1,0,0))
                  VRRz(IGrid,M,0,2)=QCz*VRRz(IGrid,M,0,1)+WQz*VRRz(IGrid,Mp1,0,1)+R1x2E*(VRRz(IGrid,M,0,0)-ZxZpE*VRRz(IGrid,Mp1,0,0))
               ENDDO
            ENDDO
            !
            DO II=0,EllBra
               DO KK=0,EllKet
                  DO IGrid=1,NGrid
                     HRRx(IGrid,II,0,KK,0)=HRRx(IGrid,I,0,KK,0)+VRRx(IGrid,0,II,KK)
                     HRRy(IGrid,II,0,KK,0)=HRRy(IGrid,I,0,KK,0)+VRRy(IGrid,0,II,KK)
                     HRRz(IGrid,II,0,KK,0)=HRRz(IGrid,I,0,KK,0)+VRRz(IGrid,0,II,KK)
                  ENDDO
               ENDDO
            ENDDO
            !
         ENDDO ! (M0| loop
      ENDDO ! |N0) loop
      ! Move Ell from A to B       
      DO I=1,EllA
         Ip1=I+1
         Im1=I-1
         DO J=1,EllC
            Jm1=J-1
            DO K=0,EllKet
               DO IGrid=1,NGrid
                  HRRx(IGrid,I,J,K,0)=HRRx(IGrid,Ip1,Jm1,K,0)+ABx*HRRx(IGrid,Im1,Jm1,K,0)
                  HRRy(IGrid,I,J,K,0)=HRRy(IGrid,Ip1,Jm1,K,0)+ABy*HRRy(IGrid,Im1,Jm1,K,0)
                  HRRz(IGrid,I,J,K,0)=HRRz(IGrid,Ip1,Jm1,K,0)+ABz*HRRz(IGrid,Im1,Jm1,K,0)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      ! Then, move some Ell from C to D 
      DO I=0,EllA
         DO J=0,EllC
            DO K=1,EllB
               Km1=K-1
               Kp1=K+1
               DO L=1,EllD
                  Lm1=L-1
                  DO IGrid=1,NGrid
                     HRRx(IGrid,I,J,K,L)=HRRx(IGrid,I,J,Kp1,Lm1)+CDx*HRRx(IGrid,I,J,Km1,Lm1)
                     HRRy(IGrid,I,J,K,L)=HRRy(IGrid,I,J,Kp1,Lm1)+CDy*HRRy(IGrid,I,J,Km1,Lm1)
                     HRRz(IGrid,I,J,K,L)=HRRz(IGrid,I,J,Kp1,Lm1)+CDz*HRRz(IGrid,I,J,Km1,Lm1)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      ! Fill in the integrals using numerics of the Laplace transform to achieve xyz coupling
      DO I=BegA,EndA
         Ix=XDex(I)
         Iy=YDex(I)
         Iz=ZDex(I)
         AOff=(OA+I-BegA)*LDA
         DO J=BegB,EndB
            Jx=XDex(J)
            Jy=YDex(J)
            Jz=ZDex(J)
            BOff=(OB+J-BegB)*LDB+AOff
            DO K=BegC,EndC
               Kx=XDex(K)
               Ky=YDex(K)
               Kz=ZDex(K)
               COff=(OC+K-BegC)*LDC+BOff
               DO L=BegD,EndD
                  Lx=XDex(L)
                  Ly=YDex(L)
                  Lz=ZDex(L)
                  OffSet=COff+(OD+L-BegD)*LDD 
                  Offset=1
                  INTGRL(OffSet)=Zero
                  DO IGrid=1,NGrid
                     INTGRL(OffSet)=INTGRL(OffSet)+HRRx(IGrid,Ix,Jx,Kx,Lx)*HRRy(IGrid,Iy,Jy,Ky,Ly)  &
                                                  *HRRz(IGrid,Iz,Jz,Kz,Lz)*DEJacobian(IGrid)
                  ENDDO
                  WRITE(*,*)' Int = ',INTGRL(OFFSET)
                  ! Probably need to do some renormalization of contraction coefficients here...
               ENDDO
            ENDDO
         ENDDO
      ENDDO
    END SUBROUTINE DEInts
!
  END MODULE DEIntegrals
