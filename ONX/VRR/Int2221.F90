SUBROUTINE Int2221(N,IntCode,CBra,CKet,DisBufB,PrmBufB,DB,IB,SB,C,U)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
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
  REAL(DOUBLE)  :: U(62),C(N,68)
  REAL(DOUBLE)  :: Cx,Cy,Cz
  REAL(DOUBLE)  :: Zeta,Up,Px,Py,Pz,cp1,cp2,cp3
  REAL(DOUBLE)  :: Eta, Uq,Qx,Qy,Qz,cq1,cq2,cq3
  REAL(DOUBLE)  :: EtaQx,EtaQy,EtaQz
  REAL(DOUBLE)  :: r1xZpE,ZpE,Rkk,r1xZ,r1x2Z,rpxZ
  REAL(DOUBLE)  :: PAx,PAy,PAz,PQx,PQy,PQz
  REAL(DOUBLE)  :: QBx,QBy,QBz,WQx,WQy,WQz
  REAL(DOUBLE)  :: Wx,Wy,Wz,WPx,WPy,WPz
  REAL(DOUBLE)  :: T1,T2,T3,TwoT
  REAL(DOUBLE)  :: R1,R2,R3,R4,G0,G1,G2,G3,ET
  REAL(DOUBLE)  :: Ts0,Ts1,r1xE,r1x2E,rpxE,r1x2ZE
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  INTEGER       :: I,J,K,M,l
  INTEGER       :: I0,I1,I2

  IF(N.EQ.0) RETURN

  C=0.0D0

  Cx=DisBufB( 8)
  Cy=DisBufB( 9)
  Cz=DisBufB(10)

  IF(IntCode.EQ.02020201) THEN

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

        IF (T1<IB%Switch) THEN
          M=INT(T1*IB%Grid)
          T2=T1*T1
          T3=T2*T1
          TwoT=2.0D0*T1
          ET=IB%ET%D(0,M)+T1*IB%ET%D(1,M)+T2*IB%ET%D(2,M)+T3*IB%ET%D(3,M)
          G3=IB%GT%D(0,M)+T1*IB%GT%D(1,M)+T2*IB%GT%D(2,M)+T3*IB%GT%D(3,M)
          G2=.200000000000000D+00*(TwoT*G3+ET)
          G1=.333333333333333D+00*(TwoT*G2+ET)
          G0=TwoT*G1+ET
          R1=Rkk*G0
          R2=Rkk*G1
          R3=Rkk*G2
          R4=Rkk*G3
        ELSE
          T2=1.0D0/T1
          T3=DSQRT(T2)
          R1=Rkk*IB%GammaA%D(1)*T3
          T3=T3*T2
          R2=Rkk*IB%GammaA%D(2)*T3
          T3=T3*T2
          R3=Rkk*IB%GammaA%D(3)*T3
          T3=T3*T2
          R4=Rkk*IB%GammaA%D(4)*T3
        ENDIF

      U(1)=R1
      U(6)=R2
      U(5)=R3
      U(2)=R4
      U(18)=QBx*U(1)+WQx*U(6)
      U(35)=QBy*U(1)+WQy*U(6)
      U(52)=QBz*U(1)+WQz*U(6)
      U(8)=QBx*U(6)+WQx*U(5)
      U(9)=QBy*U(6)+WQy*U(5)
      U(11)=QBz*U(6)+WQz*U(5)
      U(7)=QBx*U(5)+WQx*U(2)
      U(10)=QBy*U(5)+WQy*U(2)
      U(12)=QBz*U(5)+WQz*U(2)
      U(2)=PAx*U(1)+WPx*U(6)
      U(3)=PAy*U(1)+WPy*U(6)
      U(4)=PAz*U(1)+WPz*U(6)
      U(13)=PAx*U(6)+WPx*U(5)
      U(14)=PAy*U(6)+WPy*U(5)
      U(15)=PAz*U(6)+WPz*U(5)
      U(19)=PAx*U(18)+WPx*U(8)+r1x2ZE*U(6)
      U(20)=PAy*U(18)+WPy*U(8)
      U(21)=PAz*U(18)+WPz*U(8)
      U(36)=PAx*U(35)+WPx*U(9)
      U(37)=PAy*U(35)+WPy*U(9)+r1x2ZE*U(6)
      U(38)=PAz*U(35)+WPz*U(9)
      U(53)=PAx*U(52)+WPx*U(11)
      U(54)=PAy*U(52)+WPy*U(11)
      U(55)=PAz*U(52)+WPz*U(11)+r1x2ZE*U(6)
      U(16)=PAx*U(8)+WPx*U(7)+r1x2ZE*U(5)
      U(17)=PAy*U(8)+WPy*U(7)
      U(26)=PAz*U(8)+WPz*U(7)
      U(28)=PAx*U(9)+WPx*U(10)
      U(29)=PAy*U(9)+WPy*U(10)+r1x2ZE*U(5)
      U(30)=PAz*U(9)+WPz*U(10)
      U(31)=PAx*U(11)+WPx*U(12)
      U(32)=PAy*U(11)+WPy*U(12)
      U(33)=PAz*U(11)+WPz*U(12)+r1x2ZE*U(5)
      U(22)=PAx*U(19)+WPx*U(16)+r1x2Z*(U(18)-rpxZ*U(8))+r1x2ZE*U(13)
      U(23)=PAx*U(20)+WPx*U(17)+r1x2ZE*U(14)
      U(25)=PAx*U(21)+WPx*U(26)+r1x2ZE*U(15)
      U(41)=PAy*U(37)+WPy*U(29)+r1x2Z*(U(35)-rpxZ*U(9))+r1x2ZE*U(14)
      U(43)=PAy*U(38)+WPy*U(30)+r1x2ZE*U(15)
      U(61)=PAz*U(55)+WPz*U(33)+r1x2Z*(U(52)-rpxZ*U(11))+r1x2ZE*U(15)
      U(5)=PAx*U(2)+WPx*U(13)+r1x2Z*(U(1)-rpxZ*U(6))
      U(7)=PAy*U(3)+WPy*U(14)+r1x2Z*(U(1)-rpxZ*U(6))
      U(10)=PAz*U(4)+WPz*U(15)+r1x2Z*(U(1)-rpxZ*U(6))
      U(24)=PAy*U(20)+WPy*U(17)+r1x2Z*(U(18)-rpxZ*U(8))
      U(27)=PAz*U(21)+WPz*U(26)+r1x2Z*(U(18)-rpxZ*U(8))
      U(39)=PAx*U(36)+WPx*U(28)+r1x2Z*(U(35)-rpxZ*U(9))
      U(44)=PAz*U(38)+WPz*U(30)+r1x2Z*(U(35)-rpxZ*U(9))
      U(56)=PAx*U(53)+WPx*U(31)+r1x2Z*(U(52)-rpxZ*U(11))
      U(58)=PAy*U(54)+WPy*U(32)+r1x2Z*(U(52)-rpxZ*U(11))
      U(6)=PAx*U(3)+WPx*U(14)
      U(8)=PAx*U(4)+WPx*U(15)
      U(9)=PAy*U(4)+WPy*U(15)
      U(26)=PAy*U(21)+WPy*U(26)
      U(40)=PAx*U(37)+WPx*U(29)
      U(42)=PAx*U(38)+WPx*U(30)
      U(57)=PAx*U(54)+WPx*U(32)
      U(59)=PAx*U(55)+WPx*U(33)
      U(60)=PAy*U(55)+WPy*U(33)

      C(I,1)=C(I,1)+U(1)*cp1*cq1
      C(I,2)=C(I,2)+U(2)*cp3*cq1
      C(I,3)=C(I,3)+U(3)*cp3*cq1
      C(I,4)=C(I,4)+U(4)*cp3*cq1
      C(I,5)=C(I,5)+U(2)*cp2*cq1
      C(I,6)=C(I,6)+U(5)*cq1
      C(I,7)=C(I,7)+U(1)*cp3*cq1
      C(I,8)=C(I,8)+U(2)*cq1
      C(I,9)=C(I,9)+U(3)*cp2*cq1
      C(I,10)=C(I,10)+U(6)*cq1
      C(I,11)=C(I,11)+U(7)*cq1
      C(I,12)=C(I,12)+U(3)*cq1
      C(I,13)=C(I,13)+U(4)*cp2*cq1
      C(I,14)=C(I,14)+U(8)*cq1
      C(I,15)=C(I,15)+U(9)*cq1
      C(I,16)=C(I,16)+U(10)*cq1
      C(I,17)=C(I,17)+U(4)*cq1
      C(I,18)=C(I,18)+U(18)*cp1
      C(I,19)=C(I,19)+U(19)*cp3
      C(I,20)=C(I,20)+U(20)*cp3
      C(I,21)=C(I,21)+U(21)*cp3
      C(I,22)=C(I,22)+U(19)*cp2
      C(I,23)=C(I,23)+U(22)
      C(I,24)=C(I,24)+U(18)*cp3
      C(I,25)=C(I,25)+U(19)
      C(I,26)=C(I,26)+U(20)*cp2
      C(I,27)=C(I,27)+U(23)
      C(I,28)=C(I,28)+U(24)
      C(I,29)=C(I,29)+U(20)
      C(I,30)=C(I,30)+U(21)*cp2
      C(I,31)=C(I,31)+U(25)
      C(I,32)=C(I,32)+U(26)
      C(I,33)=C(I,33)+U(27)
      C(I,34)=C(I,34)+U(21)
      C(I,35)=C(I,35)+U(35)*cp1
      C(I,36)=C(I,36)+U(36)*cp3
      C(I,37)=C(I,37)+U(37)*cp3
      C(I,38)=C(I,38)+U(38)*cp3
      C(I,39)=C(I,39)+U(36)*cp2
      C(I,40)=C(I,40)+U(39)
      C(I,41)=C(I,41)+U(35)*cp3
      C(I,42)=C(I,42)+U(36)
      C(I,43)=C(I,43)+U(37)*cp2
      C(I,44)=C(I,44)+U(40)
      C(I,45)=C(I,45)+U(41)
      C(I,46)=C(I,46)+U(37)
      C(I,47)=C(I,47)+U(38)*cp2
      C(I,48)=C(I,48)+U(42)
      C(I,49)=C(I,49)+U(43)
      C(I,50)=C(I,50)+U(44)
      C(I,51)=C(I,51)+U(38)
      C(I,52)=C(I,52)+U(52)*cp1
      C(I,53)=C(I,53)+U(53)*cp3
      C(I,54)=C(I,54)+U(54)*cp3
      C(I,55)=C(I,55)+U(55)*cp3
      C(I,56)=C(I,56)+U(53)*cp2
      C(I,57)=C(I,57)+U(56)
      C(I,58)=C(I,58)+U(52)*cp3
      C(I,59)=C(I,59)+U(53)
      C(I,60)=C(I,60)+U(54)*cp2
      C(I,61)=C(I,61)+U(57)
      C(I,62)=C(I,62)+U(58)
      C(I,63)=C(I,63)+U(54)
      C(I,64)=C(I,64)+U(55)*cp2
      C(I,65)=C(I,65)+U(59)
      C(I,66)=C(I,66)+U(60)
      C(I,67)=C(I,67)+U(61)
      C(I,68)=C(I,68)+U(55)




      END DO ! K
    END DO ! J
  END DO ! I

  ELSEIF(IntCode.EQ.03020301) THEN
  ELSEIF(IntCode.EQ.03030301) THEN
  ELSEIF(IntCode.EQ.03020201) THEN
  ELSEIF(IntCode.EQ.03030201) THEN
  ELSEIF(IntCode.EQ.02020301) THEN
  ELSE
    CALL Halt('Illegal IntCode in Int2221')
  ENDIF
END SUBROUTINE Int2221
