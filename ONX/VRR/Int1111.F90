SUBROUTINE Int1111(N,IntCode,CBra,CKet,DisBufB,PrmBufB,DB,IB,SB,C,U,PBC)
  USE DerivedTypes
  USE GlobalScalars
  USE GammaF0
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
  REAL(DOUBLE)  :: C(N),U(1)
!--------------------------------------------------------------------------------
! Temporary space for computing 2-e integrals
!--------------------------------------------------------------------------------
  REAL(DOUBLE)  :: Zeta,Up,PQx,PQy,PQz
  REAL(DOUBLE)  :: Eta, Uq,Qx,Qy,Qz
  REAL(DOUBLE)  :: r1xZpE,TwoE
  REAL(DOUBLE)  :: T1
  REAL(DOUBLE)  :: R1
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  INTEGER       :: I,J,K,MG,M
  INTEGER       :: I0,I1 
!--------------------------------------------------------------------------------
! Periodic Variables
!--------------------------------------------------------------------------------
  TYPE(PBCInfo) :: PBC
  REAL(DOUBLE)  :: FPQx,FPQy,FPQz
!
  IF(N.EQ.0) RETURN
!
  DO I=1,N
    I0 = SB%SLPrm%I(I)
    C(I)=0.0D0
    DO J=1,CKet
      I1=I0+(J-1)*DB%MAXP
      Eta = DB%PrmBuf%D(I1)
      Qx  = DB%PrmBuf%D(I1+1)
      Qy  = DB%PrmBuf%D(I1+2)
      Qz  = DB%PrmBuf%D(I1+3)
      Uq  = DB%PrmBuf%D(I1+4)
      TwoE= 0.0D0
      DO K=1,CBra
        Zeta = PrmBufB(1,K)
        PQx  = PrmBufB(2,K)-Qx
        PQy  = PrmBufB(3,K)-Qy
        PQz  = PrmBufB(4,K)-Qz
#ifdef PERIODIC
        FPQx = PQx*PBC%InvBoxSh(1,1)+PQy*PBC%InvBoxSh(1,2)+PQz*PBC%InvBoxSh(1,3)
        FPQy = PQy*PBC%InvBoxSh(2,2)+PQz*PBC%InvBoxSh(2,3)
        FPQz = PQz*PBC%InvBoxSh(3,3)
        IF(PBC%AutoW(1)) FPQx = FPQx-ANINT(FPQx)
        IF(PBC%AutoW(2)) FPQy = FPQy-ANINT(FPQy)
        IF(PBC%AutoW(3)) FPQz = FPQz-ANINT(FPQz)
        PQx  = FPQx*PBC%BoxShape(1,1)+FPQy*PBC%BoxShape(1,2)+FPQz*PBC%BoxShape(1,3)
        PQy  = FPQy*PBC%BoxShape(2,2)+FPQz*PBC%BoxShape(2,3)
        PQz  = FPQz*PBC%BoxShape(3,3)
#endif
        Up   = PrmBufB(5,K)
        r1xZpE = 1.0D0/(Zeta+Eta)
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
        IF (T1<Gamma_Switch)THEN
          MG=AINT(T1*Gamma_Grid)
          R1=F0_0(MG)+T1*(F0_1(MG)+T1*(F0_2(MG)+T1*(F0_3(MG)+T1*F0_4(MG))))
        ELSE
          R1=GammAss(0)/DSQRT(T1)
        ENDIF
        TwoE=TwoE+R1*Up*DSQRT(r1xZpE)
      END DO ! K
      C(I)=C(I)+TwoE*Uq
    END DO ! J
  END DO ! I
END SUBROUTINE Int1111
