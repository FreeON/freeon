
#if 0
SUBROUTINE IntSSSS4(ShlPAC2,ShlPBD2,C)
  USE DerivedTypes
  USE GlobalScalars
  USE GammaF0
  USE ONX2DataType
  IMPLICIT NONE
!--------------------------------------------------------------------------------
! Distribution buffer stuff
!--------------------------------------------------------------------------------
  TYPE(ShellPair2), POINTER :: ShlPAC2,ShlPBD2
  REAL(DOUBLE), DIMENSION( 1) :: C
!--------------------------------------------------------------------------------
! Temporary space for computing 2-e integrals
!--------------------------------------------------------------------------------
  REAL(DOUBLE)  :: Zeta,InvZeta,Up,PQx,PQy,PQz
  REAL(DOUBLE)  :: Eta,InvEta,Uq,Qx,Qy,Qz
  REAL(DOUBLE)  :: r1xZpE,TwoE
  REAL(DOUBLE)  :: T1
  REAL(DOUBLE)  :: R1
  REAL(DOUBLE)  :: Za,Zb,Zc,Zd
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  INTEGER       :: J,K,Mg
!
  C(1)=0.0D0
  DO J=1,ShlPBD2%L
     Eta = ShlPBD2%Expt(J)
     Qx  = ShlPBD2%DCx(J)
     Qy  = ShlPBD2%Dcy(J)
     Qz  = ShlPBD2%DCz(J)
     Uq  = ShlPBD2%U(J)
!!$     Eta = ShlPBD2%Cst(1,J)
!!$     Qx  = ShlPBD2%Cst(2,J)
!!$     Qy  = ShlPBD2%Cst(3,J)
!!$     Qz  = ShlPBD2%Cst(4,J)
!!$     Uq  = ShlPBD2%Cst(5,J)
     TwoE= 0.0D0
     DO K=1,ShlPAC2%L
        Zeta = ShlPAC2%Expt(K)
        PQx  = ShlPAC2%DCx(K)-Qx
        PQy  = ShlPAC2%DCy(K)-Qy
        PQz  = ShlPAC2%DCz(K)-Qz
        Up   = ShlPAC2%U(K)
!!$        Zeta = ShlPAC2%Cst(1,K)
!!$        PQx  = ShlPAC2%Cst(2,K)-Qx
!!$        PQy  = ShlPAC2%Cst(3,K)-Qy
!!$        PQz  = ShlPAC2%Cst(4,K)-Qz
!!$        Up   = ShlPAC2%Cst(5,K)
        r1xZpE = 1.0D0/(Zeta+Eta)
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
        IF (T1<Gamma_Switch)THEN
           MG=AINT(T1*Gamma_Grid)
           R1=F0_0(MG)+T1*(F0_1(MG)+T1*(F0_2(MG)+T1*(F0_3(MG)+T1*F0_4(MG))))
        ELSE
           R1=GammAss(0)/SQRT(T1)
        ENDIF
        TwoE=TwoE+R1*Up*SQRT(r1xZpE)
     END DO ! K
     C(1)=C(1)+TwoE*Uq
  END DO ! J
END SUBROUTINE IntSSSS4



SUBROUTINE IntTest(PrmBufB,LBra,PrmBufK,LKet,C)
  USE DerivedTypes
  USE GlobalScalars
  USE GammaF0
  USE ONX2DataType
  IMPLICIT NONE
!--------------------------------------------------------------------------------
! Distribution buffer stuff
!--------------------------------------------------------------------------------
  INTEGER :: LBra,LKet
  REAL(DOUBLE) :: PrmBufB(6,LBra)
  REAL(DOUBLE) :: PrmBufK(6,LKet)
  REAL(DOUBLE) :: C(NA,NB,NC,ND)
!--------------------------------------------------------------------------------
! Temporary space for computing 2-e integrals
!--------------------------------------------------------------------------------
  REAL(DOUBLE)  :: Zeta,InvZeta,Up,PQx,PQy,PQz
  REAL(DOUBLE)  :: Eta,InvEta,Uq,Qx,Qy,Qz
  REAL(DOUBLE)  :: r1xZpE,TwoE
  REAL(DOUBLE)  :: T1
  REAL(DOUBLE)  :: R1
  REAL(DOUBLE)  :: Za,Zb,Zc,Zd
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  INTEGER       :: J,K,Mg
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
        r1xZpE = 1.0D0/(Zeta+Eta)
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
        !VRR
     END DO ! K
  END DO ! J
  !HRR                           !(ac|bd)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  C(OffA,OffB  ,OffC,OffD  )=....!(ss|pxpx)
  C(OffA,OffB+1,OffC,OffD  )=....!(ss|pypx)
  C(OffA,OffB+2,OffC,OffD  )=....!(ss|pzpx)
  !                        
  C(OffA,OffB  ,OffC,OffD+1)=....!(ss|pxpy)
  C(OffA,OffB+1,OffC,OffD+1)=....!(ss|pypy)
  C(OffA,OffB+2,OffC,OffD+1)=....!(ss|pzpy)
  !                        
  C(OffA,OffB  ,OffC,OffD+2)=....!(ss|pxpz)
  C(OffA,OffB+1,OffC,OffD+2)=....!(ss|pypz)
  C(OffA,OffB+2,OffC,OffD+2)=....!(ss|pzpz)
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  C(OffA,OffB,OffC  ,OffD)=....!(sdxx|ss)
  C(OffA,OffB,OffC+1,OffD)=....!(sdxy|ss)
  C(OffA,OffB,OffC+2,OffD)=....!(sdxz|ss)
  C(OffA,OffB,OffC+3,OffD)=....!(sdyy|ss)
  C(OffA,OffB,OffC+4,OffD)=....!(sdyz|ss)
  C(OffA,OffB,OffC+5,OffD)=....!(sdzz|ss)
  !
  !
END SUBROUTINE IntTest


!!$SUBROUTINE IntSSSS(ZAAry,ACoef,NA,ZBAry,BCoef,NB,ZCAry,CCoef,NC,ZDAry,DCoef,ND,&
!!$     &             ACoor,BCoor,CCoor,DCoor,AC2,BD2,C)
!!$  USE DerivedTypes
!!$  USE GlobalScalars
!!$  USE GammaF0
!!$  IMPLICIT NONE
!!$!--------------------------------------------------------------------------------
!!$! Distribution buffer stuff
!!$!--------------------------------------------------------------------------------
!!$  INTEGER                     :: NA,NB,NC,ND
!!$  REAL(DOUBLE), DIMENSION(NA) :: ZAAry,ACoef
!!$  REAL(DOUBLE), DIMENSION(NB) :: ZBAry,BCoef
!!$  REAL(DOUBLE), DIMENSION(NC) :: ZCAry,CCoef
!!$  REAL(DOUBLE), DIMENSION(ND) :: ZDAry,DCoef
!!$  REAL(DOUBLE), DIMENSION( 1) :: C
!!$  REAL(DOUBLE), DIMENSION( 3) :: ACoor,BCoor,CCoor,DCoor
!!$  REAL(DOUBLE) :: AC2,BD2
!!$!--------------------------------------------------------------------------------
!!$! Temporary space for computing 2-e integrals
!!$!--------------------------------------------------------------------------------
!!$  REAL(DOUBLE)  :: Zeta,Up,PQx,PQy,PQz
!!$  REAL(DOUBLE)  :: Eta, Uq,Qx,Qy,Qz
!!$  REAL(DOUBLE)  :: r1xZpE,TwoE
!!$  REAL(DOUBLE)  :: T1
!!$  REAL(DOUBLE)  :: R1
!!$  REAL(DOUBLE)  :: Za,Zb,Zc,Zd
!!$!--------------------------------------------------------------------------------
!!$! Misc. internal variables
!!$!--------------------------------------------------------------------------------
!!$  INTEGER       :: IA,IB,IC,ID,Mg
!!$!
!!$  TwoE= 0.0D0
!!$  DO ID=1,ND
!!$     Zd=ZDAry(ID)
!!$     DO IB=1,NB
!!$        Zb=ZBAry(IB)
!!$        Eta=Zd+Zb
!!$        Qx=(Zd*DCoor(1)+Zb*BCoor(1))/Eta
!!$        Qy=(Zd*DCoor(2)+Zb*BCoor(2))/Eta
!!$        Qz=(Zd*DCoor(3)+Zb*BCoor(3))/Eta
!!$        Uq=5.914967172796D0*EXP(-Zd*Zb/Eta*BD2)*DCoef(ID)*BCoef(IB)/Eta
!!$        DO IC=1,NC
!!$           Zc=ZCAry(IC)
!!$           DO IA=1,NA
!!$              Za=ZAAry(IA)
!!$              Zeta=Za+Zc
!!$              PQx=(Zc*CCoor(1)+Za*ACoor(1))/Zeta-Qx
!!$              PQy=(Zc*CCoor(2)+Za*ACoor(2))/Zeta-Qy
!!$              PQz=(Zc*CCoor(3)+Za*ACoor(3))/Zeta-Qz
!!$              Up=5.914967172796D0*EXP(-Za*Zc/Zeta*AC2)*ACoef(IA)*CCoef(IC)/Zeta
!!$              r1xZpE = 1.0D0/(Zeta+Eta)
!!$              T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
!!$              IF(T1<Gamma_Switch)THEN
!!$                 MG=AINT(T1*Gamma_Grid)
!!$                 R1=F0_0(MG)+T1*(F0_1(MG)+T1*(F0_2(MG)+T1*(F0_3(MG)+T1*F0_4(MG))))
!!$              ELSE
!!$                 R1=GammAss(0)/SQRT(T1)
!!$              ENDIF
!!$              TwoE=TwoE+R1*Up*Uq*SQRT(r1xZpE)
!!$           ENDDO
!!$        ENDDO
!!$        
!!$     ENDDO
!!$  ENDDO
!!$  C(1)=TwoE
!!$END SUBROUTINE IntSSSS



!!$SUBROUTINE IntSSSS2(C)
!!$!SUBROUTINE IntSSSS2(ShlPAC,ShlPBD,C)
!!$  USE DerivedTypes
!!$  USE GlobalScalars
!!$  USE GammaF0
!!$  USE ONX2DataType
!!$  IMPLICIT NONE
!!$!--------------------------------------------------------------------------------
!!$! Distribution buffer stuff
!!$!--------------------------------------------------------------------------------
!!$!  TYPE(ShellPair) :: ShlPAC,ShlPBD
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
!!$  INTEGER       :: IA,IB,IC,ID,Mg
!!$!
!!$  TwoE= 0.0D0
!!$  DO ID=1,ShlPBD%N2
!!$     Zd=ShlPBD%Expt2(ID)
!!$     DO IB=1,ShlPBD%N1
!!$        Zb=ShlPBD%Expt1(IB)
!!$        Eta=Zd+Zb
!!$        !InvEta=1.0d0/Eta
!!$        Qx=1.0d0!(Zd*ShlPBD%At2x+Zb*ShlPBD%At1x)*InvEta
!!$        Qy=1.0d0!(Zd*ShlPBD%At2y+Zb*ShlPBD%At1y)*InvEta
!!$        Qz=1.0d0!(Zd+ShlPBD%At2z*Zb*ShlPBD%At1z)*InvEta
!!$        Uq=5.914967172796D0!*EXP(-Zd*Zb*InvEta*ShlPBD%R12)*ShlPBD%Coef1(ID)*ShlPBD%Coef2(IB)*InvEta
!!$        DO IC=1,ShlPAC%N2
!!$           Zc=ShlPAC%Expt2(IC)
!!$           DO IA=1,ShlPAC%N1
!!$              Za=ShlPAC%Expt1(IA)
!!$              Zeta=Za+Zc
!!$              !InvZeta=1.0d0/Zeta
!!$              PQx=1.0d0!(Zc*ShlPAC%At2x+Za*ShlPAC%At1x)*InvZeta-Qx
!!$              PQy=1.0d0!(Zc*ShlPAC%At2y+Za*ShlPAC%At1y)*InvZeta-Qy
!!$              PQz=1.0d0!(Zc*ShlPAC%At2z+Za*ShlPAC%At1z)*InvZeta-Qz
!!$              Up=5.914967172796D0!*EXP(-Za*Zc*InvZeta*ShlPAC%R12)*ShlPAC%Coef1(IA)*ShlPAC%Coef2(IC)*InvZeta
!!$              r1xZpE = 1.0D0/(Zeta+Eta)
!!$              T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
!!$              IF(T1<Gamma_Switch)THEN
!!$                 MG=AINT(T1*Gamma_Grid)
!!$                 R1=F0_0(MG)+T1*(F0_1(MG)+T1*(F0_2(MG)+T1*(F0_3(MG)+T1*F0_4(MG))))
!!$              ELSE
!!$                 R1=GammAss(0)/SQRT(T1)
!!$              ENDIF
!!$              TwoE=TwoE+R1*Up*Uq*SQRT(r1xZpE)
!!$              !TwoE=TwoE*za*zc!+R1*Up*Uq*SQRT(r1xZpE)
!!$           ENDDO
!!$        ENDDO
!!$        
!!$     ENDDO
!!$  ENDDO
!!$  C(1)=TwoE ! *Uq*Qx*Qy*Qz
!!$END SUBROUTINE IntSSSS2


SUBROUTINE IntSSSS3(ShlPAC2,ShlPBD2,C)
!SUBROUTINE IntSSSS2(ShlPAC2,ShlPBD2,C)
  USE DerivedTypes
  USE GlobalScalars
  USE GammaF0
  USE ONX2DataType
  IMPLICIT NONE
!--------------------------------------------------------------------------------
! Distribution buffer stuff
!--------------------------------------------------------------------------------
  TYPE(ShellPair2) :: ShlPAC2,ShlPBD2
  REAL(DOUBLE), DIMENSION(1) :: C
!--------------------------------------------------------------------------------
! Temporary space for computing 2-e integrals
!--------------------------------------------------------------------------------
  REAL(DOUBLE)  :: Zeta,InvZeta,Up,PQx,PQy,PQz
  REAL(DOUBLE)  :: Eta,InvEta,Uq,Qx,Qy,Qz
  REAL(DOUBLE)  :: r1xZpE,TwoE
  REAL(DOUBLE)  :: T1
  REAL(DOUBLE)  :: R1
  REAL(DOUBLE)  :: Za,Zb,Zc,Zd
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  INTEGER       :: J,K,Mg
!
  C(1)=0.0D0
  DO J=1,ShlPBD2%L
!!$     Eta = ShlPBD2%Expt(J)
!!$     Qx  = ShlPBD2%DCx(J)
!!$     Qy  = ShlPBD2%Dcy(J)
!!$     Qz  = ShlPBD2%DCz(J)
!!$     Uq  = ShlPBD2%U(J)
     Eta = ShlPBD2%Cst(1,J)
     Qx  = ShlPBD2%Cst(2,J)
     Qy  = ShlPBD2%Cst(3,J)
     Qz  = ShlPBD2%Cst(4,J)
     Uq  = ShlPBD2%Cst(5,J)
     TwoE= 0.0D0
     DO K=1,ShlPAC2%L
!!$        Zeta = ShlPAC2%Expt(K)
!!$        PQx  = ShlPAC2%DCx(K)-Qx
!!$        PQy  = ShlPAC2%DCy(K)-Qy
!!$        PQz  = ShlPAC2%DCz(K)-Qz
!!$        Up   = ShlPAC2%U(K)
        Zeta = ShlPAC2%Cst(1,K)
        PQx  = ShlPAC2%Cst(2,K)-Qx
        PQy  = ShlPAC2%Cst(3,K)-Qy
        PQz  = ShlPAC2%Cst(4,K)-Qz
        Up   = ShlPAC2%Cst(5,K)
        r1xZpE = 1.0D0/(Zeta+Eta)
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
        IF (T1<Gamma_Switch)THEN
           MG=AINT(T1*Gamma_Grid)
           R1=F0_0(MG)+T1*(F0_1(MG)+T1*(F0_2(MG)+T1*(F0_3(MG)+T1*F0_4(MG))))
        ELSE
           R1=GammAss(0)/SQRT(T1)
        ENDIF
        TwoE=TwoE+R1*Up*SQRT(r1xZpE)
     END DO ! K
     C(1)=C(1)+TwoE*Uq
  END DO ! J
END SUBROUTINE IntSSSS3
#endif


SUBROUTINE IntSSSS5(PrmBufB,LBra,PrmBufK,LKet,C)
  USE DerivedTypes
  USE GlobalScalars
  USE GammaF0
  USE ONX2DataType
  IMPLICIT NONE
!--------------------------------------------------------------------------------
! Distribution buffer stuff
!--------------------------------------------------------------------------------
  INTEGER :: LBra,LKet
  REAL(DOUBLE) :: PrmBufB(5,LBra)
  REAL(DOUBLE) :: PrmBufK(5,LKet)
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
!
  C(1)=0.0D0
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
        r1xZpE = 1.0D0/(Zeta+Eta)
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
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


SUBROUTINE IntSSSS5N(PrmBufB,LBra,PrmBufK,LKet,OffSet,C)
  USE DerivedTypes
  USE GlobalScalars
  USE GammaF0
  USE ONX2DataType
  IMPLICIT NONE
!--------------------------------------------------------------------------------
! Distribution buffer stuff
!--------------------------------------------------------------------------------
  INTEGER :: LBra,LKet
  REAL(DOUBLE) :: PrmBufB(5,LBra)
  REAL(DOUBLE) :: PrmBufK(5,LKet)
  REAL(DOUBLE), DIMENSION(3,3,3,3) :: C
  TYPE(OffSt), INTENT(IN) :: OffSet
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
!
!  C(1,1,1,1)=0.0D0
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
        r1xZpE = 1.0D0/(Zeta+Eta)
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
        IF(T1<Gamma_Switch) THEN
           MG=AINT(T1*Gamma_Grid)
           R1=F0_0(MG)+T1*(F0_1(MG)+T1*(F0_2(MG)+T1*(F0_3(MG)+T1*F0_4(MG))))
        ELSE
           R1=GammAss(0)/DSQRT(T1)
        ENDIF
        TwoE=TwoE+R1*Up*DSQRT(r1xZpE)*Uq
     END DO ! K
     !C(1)=C(1)+TwoE*Uq
  END DO ! J
  C(OffSet%A,OffSet%B,OffSet%C,OffSet%D)=TwoE
END SUBROUTINE IntSSSS5N



SUBROUTINE IntSSSS5N_2(PrmBufB,LBra,PrmBufK,LKet,OffSet,C)
  USE DerivedTypes
  USE GlobalScalars
  USE GammaF0
  USE ONX2DataType
  IMPLICIT NONE
!--------------------------------------------------------------------------------
! Distribution buffer stuff
!--------------------------------------------------------------------------------
  INTEGER :: LBra,LKet
  REAL(DOUBLE) :: PrmBufB(5,LBra)
  REAL(DOUBLE) :: PrmBufK(5,LKet)
  REAL(DOUBLE), DIMENSION(81) :: C
  TYPE(OffSt), INTENT(IN) :: OffSet
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
!
!  C(1,1,1,1)=0.0D0
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
        r1xZpE = 1.0D0/(Zeta+Eta)
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
        IF(T1<Gamma_Switch) THEN
           MG=AINT(T1*Gamma_Grid)
           R1=F0_0(MG)+T1*(F0_1(MG)+T1*(F0_2(MG)+T1*(F0_3(MG)+T1*F0_4(MG))))
        ELSE
           R1=GammAss(0)/DSQRT(T1)
        ENDIF
        TwoE=TwoE+R1*Up*DSQRT(r1xZpE)*Uq
     END DO ! K
     !C(1)=C(1)+TwoE*Uq
  END DO ! J

  C(     OffSet%A               + &
       & OffSet%B*(OffSet%B-1)*3+ &
       & OffSet%C*(OffSet%C-1)*9+ &
       & OffSet%D*(OffSet%D-1)*27)=TwoE
END SUBROUTINE IntSSSS5N_2



SUBROUTINE IntSSSS6(PrmBufB,LB,PrmBufK,LK,C)
  USE DerivedTypes
  USE GlobalScalars
  USE GammaF0
  USE ONX2DataType
  IMPLICIT NONE
!--------------------------------------------------------------------------------
! Distribution buffer stuff
!--------------------------------------------------------------------------------
  INTEGER :: LB,LK
  REAL(DOUBLE) :: PrmBufB(5*LB)
  REAL(DOUBLE) :: PrmBufK(5*LK)
  REAL(DOUBLE), DIMENSION( 1) :: C
!--------------------------------------------------------------------------------
! Temporary space for computing 2-e integrals
!--------------------------------------------------------------------------------
  REAL(DOUBLE)  :: Zeta,InvZeta,Up,PQx,PQy,PQz
  REAL(DOUBLE)  :: Eta,InvEta,Uq,Qx,Qy,Qz
  REAL(DOUBLE)  :: r1xZpE,TwoE
  REAL(DOUBLE)  :: T1
  REAL(DOUBLE)  :: R1
  REAL(DOUBLE)  :: Za,Zb,Zc,Zd
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  INTEGER       :: J,K,Mg,JJ,KK
!
  C(1)=0.0D0
  DO J=1,LK
     JJ=(J-1)*5
     Eta = PrmBufK(1+JJ)
     Qx  = PrmBufK(2+JJ)
     Qy  = PrmBufK(3+JJ)
     Qz  = PrmBufK(4+JJ)
     Uq  = PrmBufK(5+JJ)
     TwoE= 0.0D0
     DO K=1,LB
        KK=(K-1)*5
        Zeta = PrmBufB(1+KK)
        PQx  = PrmBufB(2+KK)-Qx
        PQy  = PrmBufB(3+KK)-Qy
        PQz  = PrmBufB(4+KK)-Qz
        Up   = PrmBufB(5+KK)
        r1xZpE = 1.0D0/(Zeta+Eta)
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
        IF (T1<Gamma_Switch)THEN
           MG=AINT(T1*Gamma_Grid)
           R1=F0_0(MG)+T1*(F0_1(MG)+T1*(F0_2(MG)+T1*(F0_3(MG)+T1*F0_4(MG))))
        ELSE
           R1=GammAss(0)/SQRT(T1)
        ENDIF
        TwoE=TwoE+R1*Up*SQRT(r1xZpE)
     END DO ! K
     C(1)=C(1)+TwoE*Uq
  END DO ! J
END SUBROUTINE IntSSSS6





SUBROUTINE Int1111(N,IntCode,CBra,CKet,DisBufB,PrmBufB,DB,IB,SB,C,U,PBC)
  USE DerivedTypes
  USE GlobalScalars
  USE GammaF0
use ONX2DataType
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
        Up   = PrmBufB(5,K)
        r1xZpE = 1.0D0/(Zeta+Eta)
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
        IF (T1<Gamma_Switch)THEN
          MG=AINT(T1*Gamma_Grid)
          R1=F0_0(MG)+T1*(F0_1(MG)+T1*(F0_2(MG)+T1*(F0_3(MG)+T1*F0_4(MG))))
        ELSE
          R1=GammAss(0)/SQRT(T1)
        ENDIF
        TwoE=TwoE+R1*Up*SQRT(r1xZpE)
      END DO ! K
      C(I)=C(I)+TwoE*Uq
    END DO ! J
  END DO ! I
END SUBROUTINE Int1111
