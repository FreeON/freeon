SUBROUTINE Int2111(N,CBra,CKet,DisBufB,PrmBufB,DB,IB,SB,C,U)
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
  INTEGER       :: N,CBra,CKet
  REAL(DOUBLE)  :: DisBufB(DB%MAXC)
  REAL(DOUBLE)  :: PrmBufB(DB%MAXP,CBra+DB%MInfo)
!--------------------------------------------------------------------------------
! Temporary space for computing 2-e integrals
!--------------------------------------------------------------------------------
  REAL(DOUBLE)  :: U(4),C(N,4)
  REAL(DOUBLE)  :: Cx,Cy,Cz
  REAL(DOUBLE)  :: Zeta,Up,Px,Py,Pz,cp1
  REAL(DOUBLE)  :: Eta, Uq,Qx,Qy,Qz
  REAL(DOUBLE)  :: EtaQx,EtaQy,EtaQz
  REAL(DOUBLE)  :: r1xZpE,ZpE,Rkk
  REAL(DOUBLE)  :: PAx,PAy,PAz,PQx,PQy,PQz
  REAL(DOUBLE)  :: Wx,Wy,Wz,WPx,WPy,WPz
  REAL(DOUBLE)  :: T1,T2,T3
  REAL(DOUBLE)  :: R1,R2,G0,G1,ET
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  INTEGER       :: I,J,K,M
  INTEGER       :: I0,I1

  IF(N.EQ.0) RETURN

  C=0.0D0
  Cx=DisBufB( 8)
  Cy=DisBufB( 9)
  Cz=DisBufB(10)
  DO I=1,N
    I0 = SB%SLPrm%I(I)
    DO J=1,CKet
      I1=I0+(J-1)*DB%MAXP
      Eta = DB%PrmBuf%D(I1)
      Qx  = DB%PrmBuf%D(I1+1)
      Qy  = DB%PrmBuf%D(I1+2)
      Qz  = DB%PrmBuf%D(I1+3)
      Uq  = DB%PrmBuf%D(I1+4)
      EtaQx = Eta*Qx
      EtaQy = Eta*Qy
      EtaQz = Eta*Qz
     
      DO K=1,CBra
        Zeta = PrmBufB(1,K)
        Px   = PrmBufB(2,K)
        Py   = PrmBufB(3,K)
        Pz   = PrmBufB(4,K)
        Up   = PrmBufB(5,K)
        cp1  = PrmBufB(6,K)
        ZpE  = Zeta+Eta
        r1xZpE = 1.0D0/ZpE
        Rkk  = Up*Uq*DSQRT(r1xZpE)
        PAx = Px-Cx
        PAy = Py-Cy
        PAz = Pz-Cz
        Wx  = (Zeta*Px+EtaQx)*r1xZpE
        Wy  = (Zeta*Py+EtaQy)*r1xZpE
        Wz  = (Zeta*Pz+EtaQz)*r1xZpE
        WPx = Wx - Px
        WPy = Wy - Py
        WPz = Wz - Pz
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
        U(2)=PAx*U(1)+WPx*U(4)
        U(3)=PAy*U(1)+WPy*U(4)
        U(4)=PAz*U(1)+WPz*U(4)
        C(I,1)=C(I,1)+U(1)*cp1
        C(I,2)=C(I,2)+U(2)
        C(I,3)=C(I,3)+U(3)
        C(I,4)=C(I,4)+U(4)

      END DO ! K
    END DO ! J
  END DO ! I

END SUBROUTINE Int2111
