SUBROUTINE Int1166(N,IntCode,CBra,CKet,DisBufB,PrmBufB,DB,IB,SB,C,U)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE GammaF4
  USE InvExp
  IMPLICIT NONE
!--------------------------------------------------------------------------------
! Distribution buffer stuff
!--------------------------------------------------------------------------------
  TYPE(DBuf),INTENT(IN)    :: DB            ! ONX distribution buffers
  TYPE(IBuf),INTENT(INOUT) :: IB            ! ONX 2-e eval buffers
  TYPE(DSL),INTENT(IN)     :: SB            ! ONX distribution pointers
  INTEGER       :: N,IntCode,CBra,CKet
  REAL(DOUBLE)  :: DisBufB(DB%MAXC)
  REAL(DOUBLE)  :: PrmBufB(DB%MAXP,CBra+DB%MInfo)
!--------------------------------------------------------------------------------
! Temporary space for computing 2-e integrals
!--------------------------------------------------------------------------------
  REAL(DOUBLE)  :: U(48),C(N,31)
  REAL(DOUBLE)  :: Cx,Cy,Cz
  REAL(DOUBLE)  :: Zeta,Up,Px,Py,Pz,cp1,cp2,cp3
  REAL(DOUBLE)  :: Eta, Uq,Qx,Qy,Qz,cq1,cq2,cq3
  REAL(DOUBLE)  :: EtaQx,EtaQy,EtaQz
  REAL(DOUBLE)  :: r1xZpE,ZpE,Rkk,r1xZ,r1x2Z,rpxZ
  REAL(DOUBLE)  :: PAx,PAy,PAz,PQx,PQy,PQz
  REAL(DOUBLE)  :: QBx,QBy,QBz,WQx,WQy,WQz
  REAL(DOUBLE)  :: Wx,Wy,Wz,WPx,WPy,WPz
  REAL(DOUBLE)  :: T1,T2,T3,T4,TwoT
  REAL(DOUBLE)  :: R1,R2,R3,R4,R5,G0,G1,G2,G3,G4,ET
  REAL(DOUBLE)  :: Ts0,Ts1,Ts2,r1xE,r1x2E,rpxE,r1x2ZE
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  INTEGER       :: I,J,K,ME,MG,l
  INTEGER       :: I0,I1,I2

  IF(N.EQ.0) RETURN

  C=0.0D0
  Cx=DisBufB( 8)
  Cy=DisBufB( 9)
  Cz=DisBufB(10)

  DO I=1,N
    I0 = SB%SLPrm%I(I)
    I2 = SB%SLDis%I(I)-4
    DO J=1,CKet
      I1=I0+(J-1)*DB%MAXP
      Eta = DB%PrmBuf%D(I1)
      Qx  = DB%PrmBuf%D(I1+1)
      Qy  = DB%PrmBuf%D(I1+2)
      Qz  = DB%PrmBuf%D(I1+3)
      Uq  = DB%PrmBuf%D(I1+4)
      cq1 = DB%PrmBuf%D(I1+5)
      cq2 = DB%PrmBuf%D(I1+6)
      cq3 = DB%PrmBuf%D(I1+7)
      r1xE  = 1.0D0/Eta
      r1x2E = 0.5D0*r1xE
      EtaQx = Eta*Qx
      EtaQy = Eta*Qy
      EtaQz = Eta*Qz
      QBx = Qx-DB%DisBuf%D(I2+7)
      QBy = Qy-DB%DisBuf%D(I2+8)
      QBz = Qz-DB%DisBuf%D(I2+9)
      DO K=1,CBra
        Zeta = PrmBufB(1,K)
        Px   = PrmBufB(2,K)
        Py   = PrmBufB(3,K)
        Pz   = PrmBufB(4,K)
        Up   = PrmBufB(5,K)
        cp1  = PrmBufB(6,K)
        cp2  = PrmBufB(7,K)
        cp3  = PrmBufB(8,K)
        ZpE  = Zeta+Eta
        r1xZpE = 1.0D0/ZpE
        r1x2ZE = 0.5D0*r1xZpE
        r1xZ   = 1.0D0/Zeta
        r1x2Z  = 0.5D0*r1xZ
        rpxZ   = Eta*r1xZpE
        rpxE   = Zeta*r1xZpE
        Rkk = Up*Uq*DSQRT(r1xZpE)
        PAx = Px-Cx
        PAy = Py-Cy
        PAz = Pz-Cz
        Wx  = (Zeta*Px+EtaQx)*r1xZpE
        Wy  = (Zeta*Py+EtaQy)*r1xZpE
        Wz  = (Zeta*Pz+EtaQz)*r1xZpE
        WPx = Wx-Px
        WPy = Wy-Py
        WPz = Wz-Pz
        WQx = Wx-Qx
        WQy = Wy-Qy
        WQz = Wz-Qz
        PQx = Px-Qx
        PQy = Py-Qy
        PQz = Pz-Qz
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
        IF (T1<Gamma_Switch)THEN
          MG=AINT(T1*Gamma_Grid)
          ME=AINT(T1*Exp_Grid)
          T2=T1*T1
          T3=T2*T1
          T4=T2*T2
          TwoT=2.0D0*T1
          ET=Exp_0(ME)+T1*Exp_1(ME)+T2*Exp_2(ME)+T3*Exp_3(ME)+T4*Exp_4(ME)
          G4=F4_0(MG) +T1*F4_1(MG) +T2*F4_2(MG) +T3*F4_3(MG) +T4*F4_4(MG)
          G3=.142857142857143D+00*(TwoT*G4+ET)
          G2=.200000000000000D+00*(TwoT*G3+ET)
          G1=.333333333333333D+00*(TwoT*G2+ET)
          G0=TwoT*G1+ET
          R1=Rkk*G0
          R2=Rkk*G1
          R3=Rkk*G2
          R4=Rkk*G3
          R5=Rkk*G4
        ELSE
          T2=1.0D0/T1
          T3=DSQRT(T2)
          R1=Rkk*GammAss(0)*T3
          T3=T3*T2
          R2=Rkk*GammAss(1)*T3
          T3=T3*T2
          R3=Rkk*GammAss(2)*T3
          T3=T3*T2
          R4=Rkk*GammAss(3)*T3
          T3=T3*T2
          R5=Rkk*GammAss(4)*T3
        ENDIF

      U(1)=R1
      U(11)=R2
      U(12)=R3
      U(13)=R4
      U(5)=R5
      U(2)=QBx*U(1)+WQx*U(11)
      U(3)=QBy*U(1)+WQy*U(11)
      U(4)=QBz*U(1)+WQz*U(11)
      U(21)=QBx*U(11)+WQx*U(12)
      U(22)=QBy*U(11)+WQy*U(12)
      U(23)=QBz*U(11)+WQz*U(12)
      U(24)=QBx*U(12)+WQx*U(13)
      U(25)=QBy*U(12)+WQy*U(13)
      U(26)=QBz*U(12)+WQz*U(13)
      U(14)=QBx*U(13)+WQx*U(5)
      U(15)=QBy*U(13)+WQy*U(5)
      U(16)=QBz*U(13)+WQz*U(5)
      Ts0=r1x2E*(U(1)-rpxE*U(11))
      U(5)=QBx*U(2)+WQx*U(21)+ Ts0
      U(6)=QBx*U(3)+WQx*U(22)
      U(7)=QBy*U(3)+WQy*U(22)+ Ts0
      U(8)=QBx*U(4)+WQx*U(23)
      U(9)=QBy*U(4)+WQy*U(23)
      U(10)=QBz*U(4)+WQz*U(23)+ Ts0
      Ts1=r1x2E*(U(11)-rpxE*U(12))
      U(27)=QBx*U(21)+WQx*U(24)+ Ts1
      U(28)=QBx*U(22)+WQx*U(25)
      U(29)=QBy*U(22)+WQy*U(25)+ Ts1
      U(30)=QBx*U(23)+WQx*U(26)
      U(31)=QBy*U(23)+WQy*U(26)
      U(33)=QBz*U(23)+WQz*U(26)+ Ts1
      Ts2=r1x2E*(U(12)-rpxE*U(13))
      U(32)=QBx*U(24)+WQx*U(14)+ Ts2
      U(34)=QBx*U(25)+WQx*U(15)
      U(35)=QBy*U(25)+WQy*U(15)+ Ts2
      U(36)=QBx*U(26)+WQx*U(16)
      U(37)=QBy*U(26)+WQy*U(16)
      U(38)=QBz*U(26)+WQz*U(16)+ Ts2
      U(11)=QBx*U(5)+WQx*U(27)+2.0D0*r1x2E*(U(2)-rpxE*U(21))
      U(12)=QBx*U(6)+WQx*U(28)+r1x2E*(U(3)-rpxE*U(22))
      U(13)=QBx*U(7)+WQx*U(29)
      U(14)=QBy*U(7)+WQy*U(29)+2.0D0*r1x2E*(U(3)-rpxE*U(22))
      U(15)=QBx*U(8)+WQx*U(30)+r1x2E*(U(4)-rpxE*U(23))
      U(16)=QBx*U(9)+WQx*U(31)
      U(17)=QBy*U(9)+WQy*U(31)+r1x2E*(U(4)-rpxE*U(23))
      U(18)=QBx*U(10)+WQx*U(33)
      U(19)=QBy*U(10)+WQy*U(33)
      U(20)=QBz*U(10)+WQz*U(33)+2.0D0*r1x2E*(U(4)-rpxE*U(23))
      U(39)=QBx*U(27)+WQx*U(32)+2.0D0*r1x2E*(U(21)-rpxE*U(24))
      U(40)=QBx*U(28)+WQx*U(34)+r1x2E*(U(22)-rpxE*U(25))
      U(41)=QBx*U(29)+WQx*U(35)
      U(42)=QBy*U(29)+WQy*U(35)+2.0D0*r1x2E*(U(22)-rpxE*U(25))
      U(43)=QBx*U(30)+WQx*U(36)+r1x2E*(U(23)-rpxE*U(26))
      U(44)=QBx*U(31)+WQx*U(37)
      U(45)=QBy*U(31)+WQy*U(37)+r1x2E*(U(23)-rpxE*U(26))
      U(46)=QBx*U(33)+WQx*U(38)
      U(47)=QBy*U(33)+WQy*U(38)
      U(48)=QBz*U(33)+WQz*U(38)+2.0D0*r1x2E*(U(23)-rpxE*U(26))
      U(21)=QBx*U(11)+WQx*U(39)+3.0D0*r1x2E*(U(5)-rpxE*U(27))
      U(22)=QBx*U(12)+WQx*U(40)+2.0D0*r1x2E*(U(6)-rpxE*U(28))
      U(23)=QBx*U(13)+WQx*U(41)+r1x2E*(U(7)-rpxE*U(29))
      U(25)=QBy*U(14)+WQy*U(42)+3.0D0*r1x2E*(U(7)-rpxE*U(29))
      U(26)=QBx*U(15)+WQx*U(43)+2.0D0*r1x2E*(U(8)-rpxE*U(30))
      U(27)=QBx*U(16)+WQx*U(44)+r1x2E*(U(9)-rpxE*U(31))
      U(29)=QBy*U(17)+WQy*U(45)+2.0D0*r1x2E*(U(9)-rpxE*U(31))
      U(30)=QBx*U(18)+WQx*U(46)+r1x2E*(U(10)-rpxE*U(33))
      U(32)=QBy*U(19)+WQy*U(47)+r1x2E*(U(10)-rpxE*U(33))
      U(35)=QBz*U(20)+WQz*U(48)+3.0D0*r1x2E*(U(10)-rpxE*U(33))
      U(24)=QBx*U(14)+WQx*U(42)
      U(28)=QBx*U(17)+WQx*U(45)
      U(31)=QBx*U(19)+WQx*U(47)
      U(33)=QBx*U(20)+WQx*U(48)
      U(34)=QBy*U(20)+WQy*U(48)

      C(I,1)=C(I,1)+U(5)                                             
      C(I,2)=C(I,2)+U(6)                                             
      C(I,3)=C(I,3)+U(7)                                             
      C(I,4)=C(I,4)+U(8)                                             
      C(I,5)=C(I,5)+U(9)                                             
      C(I,6)=C(I,6)+U(10)                                            
      C(I,7)=C(I,7)+U(11)                                            
      C(I,8)=C(I,8)+U(12)                                            
      C(I,9)=C(I,9)+U(13)                                            
      C(I,10)=C(I,10)+U(14)                                          
      C(I,11)=C(I,11)+U(15)                                          
      C(I,12)=C(I,12)+U(16)                                          
      C(I,13)=C(I,13)+U(17)                                          
      C(I,14)=C(I,14)+U(18)                                          
      C(I,15)=C(I,15)+U(19)                                          
      C(I,16)=C(I,16)+U(20)                                          
      C(I,17)=C(I,17)+U(21)                                          
      C(I,18)=C(I,18)+U(22)                                          
      C(I,19)=C(I,19)+U(23)                                          
      C(I,20)=C(I,20)+U(24)                                          
      C(I,21)=C(I,21)+U(25)                                          
      C(I,22)=C(I,22)+U(26)                                          
      C(I,23)=C(I,23)+U(27)                                          
      C(I,24)=C(I,24)+U(28)                                          
      C(I,25)=C(I,25)+U(29)                                          
      C(I,26)=C(I,26)+U(30)                                          
      C(I,27)=C(I,27)+U(31)                                          
      C(I,28)=C(I,28)+U(32)                                          
      C(I,29)=C(I,29)+U(33)                                          
      C(I,30)=C(I,30)+U(34)                                          
      C(I,31)=C(I,31)+U(35)                           


      END DO ! K
    END DO ! J
  END DO ! I

END SUBROUTINE Int1166
