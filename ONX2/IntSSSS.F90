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
        INCLUDE 'ERIMIC.Inc'
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


