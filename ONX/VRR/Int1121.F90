SUBROUTINE Int1121(N,IntCode,CBra,CKet,DisBufB,PrmBufB,DB,IB,SB,C,U)
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
  REAL(DOUBLE)  :: U(4),C(N,4)
  REAL(DOUBLE)  :: Cx,Cy,Cz
  REAL(DOUBLE)  :: Zeta,Up,Px,Py,Pz,cp1
  REAL(DOUBLE)  :: Eta, Uq,Qx,Qy,Qz,cq1
  REAL(DOUBLE)  :: EtaQx,EtaQy,EtaQz
  REAL(DOUBLE)  :: r1xZpE,ZpE,Rkk
  REAL(DOUBLE)  :: PAx,PAy,PAz,PQx,PQy,PQz
  REAL(DOUBLE)  :: QBx,QBy,QBz,WQx,WQy,WQz
  REAL(DOUBLE)  :: Wx,Wy,Wz,WPx,WPy,WPz
  REAL(DOUBLE)  :: T1,T2,T3
  REAL(DOUBLE)  :: R1,R2,G0,G1,ET
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  INTEGER       :: I,J,K,M
  INTEGER       :: I0,I1,I2

  IF(N.EQ.0) RETURN

  IF(IntCode.EQ.1010201) THEN

  C=0.0D0
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
        ZpE  = Zeta+Eta
        r1xZpE = 1.0D0/ZpE
        Rkk  = Up*Uq*DSQRT(r1xZpE)
        Wx  = (Zeta*Px+EtaQx)*r1xZpE
        Wy  = (Zeta*Py+EtaQy)*r1xZpE
        Wz  = (Zeta*Pz+EtaQz)*r1xZpE
        WQx = Wx - Qx
        WQy = Wy - Qy
        WQz = Wz - Qz
        PQx = Px - Qx
        PQy = Py - Qy
        PQz = Pz - Qz

        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
        IF (T1<IB%Switch) THEN
          M=INT(T1*IB%Grid)
          T2=T1*T1
          T3=T2*T1
          ET=IB%ET%D(0,M)+T1*IB%ET%D(1,M)+T2*IB%ET%D(2,M)+T3*IB%ET%D(3,M)
          G1=IB%GT%D(0,M)+T1*IB%GT%D(1,M)+T2*IB%GT%D(2,M)+T3*IB%GT%D(3,M)
          G0=2.0D0*T1*G1+ET
          R1=Rkk*G0
          R2=Rkk*G1
        ELSE
          T2=1.0D0/T1
          T3=DSQRT(T2)
          R1=Rkk*IB%GammaA%D(1)*T3
          T3=T3*T2
          R2=Rkk*IB%GammaA%D(2)*T3
        ENDIF
  
        U(1)=R1
        U(4)=R2
        U(2)=QBx*U(1)+WQx*U(4)
        U(3)=QBy*U(1)+WQy*U(4)
        U(4)=QBz*U(1)+WQz*U(4)
        C(I,1)=C(I,1)+U(1)*cq1
        C(I,2)=C(I,2)+U(2)
        C(I,3)=C(I,3)+U(3)
        C(I,4)=C(I,4)+U(4)

      END DO ! K
    END DO ! J
  END DO ! I

  ELSEIF(IntCode.EQ.1010301) THEN

  C=0.0D0
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
        ZpE  = Zeta+Eta
        r1xZpE = 1.0D0/ZpE
        Rkk  = Up*Uq*DSQRT(r1xZpE)
        Wx  = (Zeta*Px+EtaQx)*r1xZpE
        Wy  = (Zeta*Py+EtaQy)*r1xZpE
        Wz  = (Zeta*Pz+EtaQz)*r1xZpE
        WQx = Wx - Qx
        WQy = Wy - Qy
        WQz = Wz - Qz
        PQx = Px - Qx
        PQy = Py - Qy
        PQz = Pz - Qz

        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
        IF (T1<IB%Switch) THEN
          M=INT(T1*IB%Grid)
          T2=T1*T1
          T3=T2*T1
          ET=IB%ET%D(0,M)+T1*IB%ET%D(1,M)+T2*IB%ET%D(2,M)+T3*IB%ET%D(3,M)
          G1=IB%GT%D(0,M)+T1*IB%GT%D(1,M)+T2*IB%GT%D(2,M)+T3*IB%GT%D(3,M)
          G0=2.0D0*T1*G1+ET
          R1=Rkk*G0
          R2=Rkk*G1
        ELSE
          T2=1.0D0/T1
          T3=DSQRT(T2)
          R1=Rkk*IB%GammaA%D(1)*T3
          T3=T3*T2
          R2=Rkk*IB%GammaA%D(2)*T3
        ENDIF

        U(1)=R1
        U(4)=R2
        U(2)=QBx*U(1)+WQx*U(4)
        U(3)=QBy*U(1)+WQy*U(4)
        U(4)=QBz*U(1)+WQz*U(4)
        C(I,1)=C(I,1)+U(2)
        C(I,2)=C(I,2)+U(3)
        C(I,3)=C(I,3)+U(4)

      END DO ! K
    END DO ! J
  END DO ! I

  ELSE
    CALL Halt('Illegal IntCode in Int1121')
  ENDIF

END SUBROUTINE Int1121
