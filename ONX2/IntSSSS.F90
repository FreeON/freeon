SUBROUTINE IntSSSS5(PrmBufB,LBra,PrmBufK,LKet,C,PBC)
  USE DerivedTypes
  USE GlobalScalars
  USE GammaF0
  USE ONX2DataType
  IMPLICIT NONE
!--------------------------------------------------------------------------------
! Distribution buffer stuff
!--------------------------------------------------------------------------------
  INTEGER :: LBra,LKet
  REAL(DOUBLE) :: PrmBufB(7,LBra)
  REAL(DOUBLE) :: PrmBufK(7,LKet)
  REAL(DOUBLE) :: C(1)
!--------------------------------------------------------------------------------
! Temporary space for computing 2-e integrals
!--------------------------------------------------------------------------------
  REAL(DOUBLE)  :: Zeta,Up,PQx,PQy,PQz
  REAL(DOUBLE)  :: Eta,Uq,Qx,Qy,Qz
  REAL(DOUBLE)  :: r1xZpE,TwoE
  REAL(DOUBLE)  :: T1
  REAL(DOUBLE)  :: R1
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  INTEGER       :: J,K,Mg
!--------------------------------------------------------------------------------
! Periodic Variables
!--------------------------------------------------------------------------------
  TYPE(PBCInfo) :: PBC
  REAL(DOUBLE)  :: FPQx,FPQy,FPQz
!
  DO J=1,LKet
     Eta = PrmBufK(1,J)
     Qx  = PrmBufK(2,J)
     Qy  = PrmBufK(3,J)
     Qz  = PrmBufK(4,J)
     Uq  = PrmBufK(5,J)
     TwoE= 0.0D0
     DO K=1,LBra
        Zeta = PrmBufB(1,K)
        PQx  = PrmBufB(2,K)-Qx
        PQy  = PrmBufB(3,K)-Qy
        PQz  = PrmBufB(4,K)-Qz
        Up   = PrmBufB(5,K)
!
        FPQx = PQx*PBC%InvBoxSh%D(1,1)+PQy*PBC%InvBoxSh%D(1,2)+PQz*PBC%InvBoxSh%D(1,3)
        FPQy = PQy*PBC%InvBoxSh%D(2,2)+PQz*PBC%InvBoxSh%D(2,3)
        FPQz = PQz*PBC%InvBoxSh%D(3,3)
        IF(PBC%AutoW%I(1)==1) FPQx=FPQx-DNINT(FPQx-SIGN(1.D0,FPQx)*1.D-14)
        IF(PBC%AutoW%I(2)==1) FPQy=FPQy-DNINT(FPQy-SIGN(1.D0,FPQy)*1.D-14)
        IF(PBC%AutoW%I(3)==1) FPQz=FPQz-DNINT(FPQz-SIGN(1.D0,FPQz)*1.D-14)
        PQx  = FPQx*PBC%BoxShape%D(1,1)+FPQy*PBC%BoxShape%D(1,2)+FPQz*PBC%BoxShape%D(1,3)
        PQy  = FPQy*PBC%BoxShape%D(2,2)+FPQz*PBC%BoxShape%D(2,3)
        PQz  = FPQz*PBC%BoxShape%D(3,3)
!
        r1xZpE = 1.0D0/(Zeta+Eta)
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
!
        IF(T1<Gamma_Switch) THEN
           MG=AINT(T1*Gamma_Grid)
           R1=F0_0(MG)+T1*(F0_1(MG)+T1*(F0_2(MG)+T1*(F0_3(MG)+T1*F0_4(MG))))
        ELSE
           R1=GammAss(0)/DSQRT(T1)
        ENDIF
        TwoE=TwoE+R1*Up*DSQRT(r1xZpE)
     END DO ! K
     C(1)=C(1)+TwoE*Uq
  END DO ! J

END SUBROUTINE IntSSSS5


