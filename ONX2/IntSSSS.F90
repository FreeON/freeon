
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

  !write(*,*) 'intssss'

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
        IF(PBC%AutoW%I(1)==1) FPQx = FPQx-ANINT(FPQx)
        IF(PBC%AutoW%I(2)==1) FPQy = FPQy-ANINT(FPQy)
        IF(PBC%AutoW%I(3)==1) FPQz = FPQz-ANINT(FPQz)
        PQx  = FPQx*PBC%BoxShape%D(1,1)+FPQy*PBC%BoxShape%D(1,2)+FPQz*PBC%BoxShape%D(1,3)
        PQy  = FPQy*PBC%BoxShape%D(2,2)+FPQz*PBC%BoxShape%D(2,3)
        PQz  = FPQz*PBC%BoxShape%D(3,3)
!
        r1xZpE = 1.0D0/(Zeta+Eta)
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE

        !write(*,*) 'r1xZpE',r1xZpE
!write(*,*) 'Zeta',Zeta
!write(*,*) 'Eta',Eta
        !write(*,*) 'T1',T1

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


!!$SUBROUTINE IntSSSS5N(PrmBufB,LBra,PrmBufK,LKet,OffSet,C)
!!$  USE DerivedTypes
!!$  USE GlobalScalars
!!$  USE GammaF0
!!$  USE ONX2DataType
!!$  IMPLICIT NONE
!!$!--------------------------------------------------------------------------------
!!$! Distribution buffer stuff
!!$!--------------------------------------------------------------------------------
!!$  INTEGER :: LBra,LKet
!!$  REAL(DOUBLE) :: PrmBufB(5,LBra)
!!$  REAL(DOUBLE) :: PrmBufK(5,LKet)
!!$  REAL(DOUBLE), DIMENSION(3,3,3,3) :: C
!!$  TYPE(ONX2OffSt), INTENT(IN) :: OffSet
!!$!--------------------------------------------------------------------------------
!!$! Temporary space for computing 2-e integrals
!!$!--------------------------------------------------------------------------------
!!$  REAL(DOUBLE)  :: Zeta,Up,PQx,PQy,PQz
!!$  REAL(DOUBLE)  :: Eta,Uq,Qx,Qy,Qz
!!$  REAL(DOUBLE)  :: r1xZpE,TwoE
!!$  REAL(DOUBLE)  :: T1
!!$  REAL(DOUBLE)  :: R1
!!$!--------------------------------------------------------------------------------
!!$! Misc. internal variables
!!$!--------------------------------------------------------------------------------
!!$  INTEGER       :: J,K,Mg
!!$!
!!$!  C(1,1,1,1)=0.0D0
!!$  DO J=1,LKet
!!$     Eta = PrmBufK(1,J)
!!$     Qx  = PrmBufK(2,J)
!!$     Qy  = PrmBufK(3,J)
!!$     Qz  = PrmBufK(4,J)
!!$     Uq  = PrmBufK(5,J)
!!$     TwoE= 0.0D0
!!$     DO K=1,LBra
!!$        Zeta = PrmBufB(1,K)
!!$        PQx  = PrmBufB(2,K)-Qx
!!$        PQy  = PrmBufB(3,K)-Qy
!!$        PQz  = PrmBufB(4,K)-Qz
!!$        Up   = PrmBufB(5,K)
!!$        r1xZpE = 1.0D0/(Zeta+Eta)
!!$        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
!!$        IF(T1<Gamma_Switch) THEN
!!$           MG=AINT(T1*Gamma_Grid)
!!$           R1=F0_0(MG)+T1*(F0_1(MG)+T1*(F0_2(MG)+T1*(F0_3(MG)+T1*F0_4(MG))))
!!$        ELSE
!!$           R1=GammAss(0)/DSQRT(T1)
!!$        ENDIF
!!$        TwoE=TwoE+R1*Up*DSQRT(r1xZpE)*Uq
!!$     END DO ! K
!!$     !C(1)=C(1)+TwoE*Uq
!!$  END DO ! J
!!$  C(OffSet%A,OffSet%B,OffSet%C,OffSet%D)=TwoE
!!$END SUBROUTINE IntSSSS5N
!!$
!!$
!!$
!!$SUBROUTINE IntSSSS5N_2(PrmBufB,LBra,PrmBufK,LKet,OffSet,C)
!!$  USE DerivedTypes
!!$  USE GlobalScalars
!!$  USE GammaF0
!!$  USE ONX2DataType
!!$  IMPLICIT NONE
!!$!--------------------------------------------------------------------------------
!!$! Distribution buffer stuff
!!$!--------------------------------------------------------------------------------
!!$  INTEGER :: LBra,LKet
!!$  REAL(DOUBLE) :: PrmBufB(5,LBra)
!!$  REAL(DOUBLE) :: PrmBufK(5,LKet)
!!$  REAL(DOUBLE), DIMENSION(81) :: C
!!$  TYPE(ONX2OffSt), INTENT(IN) :: OffSet
!!$!--------------------------------------------------------------------------------
!!$! Temporary space for computing 2-e integrals
!!$!--------------------------------------------------------------------------------
!!$  REAL(DOUBLE)  :: Zeta,Up,PQx,PQy,PQz
!!$  REAL(DOUBLE)  :: Eta,Uq,Qx,Qy,Qz
!!$  REAL(DOUBLE)  :: r1xZpE,TwoE
!!$  REAL(DOUBLE)  :: T1
!!$  REAL(DOUBLE)  :: R1
!!$!--------------------------------------------------------------------------------
!!$! Misc. internal variables
!!$!--------------------------------------------------------------------------------
!!$  INTEGER       :: J,K,Mg
!!$!
!!$!  C(1,1,1,1)=0.0D0
!!$  DO J=1,LKet
!!$     Eta = PrmBufK(1,J)
!!$     Qx  = PrmBufK(2,J)
!!$     Qy  = PrmBufK(3,J)
!!$     Qz  = PrmBufK(4,J)
!!$     Uq  = PrmBufK(5,J)
!!$     TwoE= 0.0D0
!!$     DO K=1,LBra
!!$        Zeta = PrmBufB(1,K)
!!$        PQx  = PrmBufB(2,K)-Qx
!!$        PQy  = PrmBufB(3,K)-Qy
!!$        PQz  = PrmBufB(4,K)-Qz
!!$        Up   = PrmBufB(5,K)
!!$        r1xZpE = 1.0D0/(Zeta+Eta)
!!$        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
!!$        IF(T1<Gamma_Switch) THEN
!!$           MG=AINT(T1*Gamma_Grid)
!!$           R1=F0_0(MG)+T1*(F0_1(MG)+T1*(F0_2(MG)+T1*(F0_3(MG)+T1*F0_4(MG))))
!!$        ELSE
!!$           R1=GammAss(0)/DSQRT(T1)
!!$        ENDIF
!!$        TwoE=TwoE+R1*Up*DSQRT(r1xZpE)*Uq
!!$     END DO ! K
!!$     !C(1)=C(1)+TwoE*Uq
!!$  END DO ! J
!!$
!!$  C(     OffSet%A               + &
!!$       & OffSet%B*(OffSet%B-1)*3+ &
!!$       & OffSet%C*(OffSet%C-1)*9+ &
!!$       & OffSet%D*(OffSet%D-1)*27)=TwoE
!!$END SUBROUTINE IntSSSS5N_2
!!$
!!$
!!$
!!$SUBROUTINE IntSSSS6(PrmBufB,LB,PrmBufK,LK,C)
!!$  USE DerivedTypes
!!$  USE GlobalScalars
!!$  USE GammaF0
!!$  USE ONX2DataType
!!$  IMPLICIT NONE
!!$!--------------------------------------------------------------------------------
!!$! Distribution buffer stuff
!!$!--------------------------------------------------------------------------------
!!$  INTEGER :: LB,LK
!!$  REAL(DOUBLE) :: PrmBufB(5*LB)
!!$  REAL(DOUBLE) :: PrmBufK(5*LK)
!!$  REAL(DOUBLE), DIMENSION( 1) :: C
!!$!--------------------------------------------------------------------------------
!!$! Temporary space for computing 2-e integrals
!!$!--------------------------------------------------------------------------------
!!$  REAL(DOUBLE)  :: Zeta,InvZeta,Up,PQx,PQy,PQz
!!$  REAL(DOUBLE)  :: Eta,InvEta,Uq,Qx,Qy,Qz
!!$  REAL(DOUBLE)  :: r1xZpE,TwoE
!!$  REAL(DOUBLE)  :: T1
!!$  REAL(DOUBLE)  :: R1
!!$  REAL(DOUBLE)  :: Za,Zb,Zc,Zd
!!$!--------------------------------------------------------------------------------
!!$! Misc. internal variables
!!$!--------------------------------------------------------------------------------
!!$  INTEGER       :: J,K,Mg,JJ,KK
!!$!
!!$  C(1)=0.0D0
!!$  DO J=1,LK
!!$     JJ=(J-1)*5
!!$     Eta = PrmBufK(1+JJ)
!!$     Qx  = PrmBufK(2+JJ)
!!$     Qy  = PrmBufK(3+JJ)
!!$     Qz  = PrmBufK(4+JJ)
!!$     Uq  = PrmBufK(5+JJ)
!!$     TwoE= 0.0D0
!!$     DO K=1,LB
!!$        KK=(K-1)*5
!!$        Zeta = PrmBufB(1+KK)
!!$        PQx  = PrmBufB(2+KK)-Qx
!!$        PQy  = PrmBufB(3+KK)-Qy
!!$        PQz  = PrmBufB(4+KK)-Qz
!!$        Up   = PrmBufB(5+KK)
!!$        r1xZpE = 1.0D0/(Zeta+Eta)
!!$        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
!!$        IF (T1<Gamma_Switch)THEN
!!$           MG=AINT(T1*Gamma_Grid)
!!$           R1=F0_0(MG)+T1*(F0_1(MG)+T1*(F0_2(MG)+T1*(F0_3(MG)+T1*F0_4(MG))))
!!$        ELSE
!!$           R1=GammAss(0)/SQRT(T1)
!!$        ENDIF
!!$        TwoE=TwoE+R1*Up*SQRT(r1xZpE)
!!$     END DO ! K
!!$     C(1)=C(1)+TwoE*Uq
!!$  END DO ! J
!!$END SUBROUTINE IntSSSS6
!!$
!!$
!!$
!!$
!!$
!!$SUBROUTINE Int1111(N,IntCode,CBra,CKet,DisBufB,PrmBufB,DB,IB,SB,C,U,PBC)
!!$  USE DerivedTypes
!!$  USE GlobalScalars
!!$  USE GammaF0
!!$use ONX2DataType
!!$  IMPLICIT NONE
!!$!--------------------------------------------------------------------------------
!!$! Distribution buffer stuff
!!$!--------------------------------------------------------------------------------
!!$  TYPE(DBuf),INTENT(IN)    :: DB            ! ONX distribution buffers
!!$  TYPE(IBuf),INTENT(INOUT) :: IB            ! ONX 2-e eval buffers
!!$  TYPE(DSL),INTENT(IN)     :: SB            ! ONX distribution pointers
!!$  INTEGER       :: N,IntCode,CBra,CKet
!!$  REAL(DOUBLE)  :: DisBufB(DB%MAXC)
!!$  REAL(DOUBLE)  :: PrmBufB(DB%MAXP,CBra+DB%MInfo)
!!$  REAL(DOUBLE)  :: C(N),U(1)
!!$!--------------------------------------------------------------------------------
!!$! Temporary space for computing 2-e integrals
!!$!--------------------------------------------------------------------------------
!!$  REAL(DOUBLE)  :: Zeta,Up,PQx,PQy,PQz
!!$  REAL(DOUBLE)  :: Eta, Uq,Qx,Qy,Qz
!!$  REAL(DOUBLE)  :: r1xZpE,TwoE
!!$  REAL(DOUBLE)  :: T1
!!$  REAL(DOUBLE)  :: R1
!!$!--------------------------------------------------------------------------------
!!$! Misc. internal variables
!!$!--------------------------------------------------------------------------------
!!$  INTEGER       :: I,J,K,MG,M
!!$  INTEGER       :: I0,I1 
!!$!--------------------------------------------------------------------------------
!!$! Periodic Variables
!!$!--------------------------------------------------------------------------------
!!$  TYPE(PBCInfo) :: PBC
!!$  REAL(DOUBLE)  :: FPQx,FPQy,FPQz
!!$!
!!$  IF(N.EQ.0) RETURN
!!$!
!!$  DO I=1,N
!!$    I0 = SB%SLPrm%I(I)
!!$    C(I)=0.0D0
!!$    DO J=1,CKet
!!$      I1=I0+(J-1)*DB%MAXP
!!$      Eta = DB%PrmBuf%D(I1)
!!$      Qx  = DB%PrmBuf%D(I1+1)
!!$      Qy  = DB%PrmBuf%D(I1+2)
!!$      Qz  = DB%PrmBuf%D(I1+3)
!!$      Uq  = DB%PrmBuf%D(I1+4)
!!$      TwoE= 0.0D0
!!$      DO K=1,CBra
!!$        Zeta = PrmBufB(1,K)
!!$        PQx  = PrmBufB(2,K)-Qx
!!$        PQy  = PrmBufB(3,K)-Qy
!!$        PQz  = PrmBufB(4,K)-Qz
!!$!
!!$        FPQx = PQx*PBC%InvBoxSh%D(1,1)+PQy*PBC%InvBoxSh%D(1,2)+PQz*PBC%InvBoxSh%D(1,3)
!!$        FPQy = PQy*PBC%InvBoxSh%D(2,2)+PQz*PBC%InvBoxSh%D(2,3)
!!$        FPQz = PQz*PBC%InvBoxSh%D(3,3)
!!$        IF(PBC%AutoW%I(1)==1) FPQx = FPQx-ANINT(FPQx)
!!$        IF(PBC%AutoW%I(2)==1) FPQy = FPQy-ANINT(FPQy)
!!$        IF(PBC%AutoW%I(3)==1) FPQz = FPQz-ANINT(FPQz)
!!$        PQx  = FPQx*PBC%BoxShape%D(1,1)+FPQy*PBC%BoxShape%D(1,2)+FPQz*PBC%BoxShape%D(1,3)
!!$        PQy  = FPQy*PBC%BoxShape%D(2,2)+FPQz*PBC%BoxShape%D(2,3)
!!$        PQz  = FPQz*PBC%BoxShape%D(3,3)
!!$!
!!$        Up   = PrmBufB(5,K)
!!$        r1xZpE = 1.0D0/(Zeta+Eta)
!!$        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
!!$        IF (T1<Gamma_Switch)THEN
!!$          MG=AINT(T1*Gamma_Grid)
!!$          R1=F0_0(MG)+T1*(F0_1(MG)+T1*(F0_2(MG)+T1*(F0_3(MG)+T1*F0_4(MG))))
!!$        ELSE
!!$          R1=GammAss(0)/SQRT(T1)
!!$        ENDIF
!!$        TwoE=TwoE+R1*Up*SQRT(r1xZpE)
!!$      END DO ! K
!!$      C(I)=C(I)+TwoE*Uq
!!$    END DO ! J
!!$  END DO ! I
!!$END SUBROUTINE Int1111
