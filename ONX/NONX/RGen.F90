SUBROUTINE RGen(N,Ltot,CBra,CKet,CB,CK,DisBufB,PrmBufB,R,DB,IB,SB)
  USE DerivedTypes
  USE GlobalScalars
  USE ONXParameters
  USE PrettyPrint
  USE InvExp
  USE GammaFunctions
  IMPLICIT NONE
!--------------------------------------------------------------------------------
! Distribution buffer stuff
!--------------------------------------------------------------------------------
  TYPE(DBuf),INTENT(IN)    :: DB            ! ONX distribution buffers
  TYPE(IBuf),INTENT(INOUT) :: IB            ! ONX 2-e eval buffers
  TYPE(DSL),INTENT(IN)     :: SB            ! ONX distribution pointers
  INTEGER       :: N,Ltot,CBra,CKet
  REAL(DOUBLE)  :: DisBufB(DB%MAXC)
  REAL(DOUBLE)  :: PrmBufB(DB%MAXP,CBra+DB%MInfo)
  REAL(DOUBLE)  :: CB(CBra,3)
  REAL(DOUBLE)  :: CK(N,CKet,3)
!--------------------------------------------------------------------------------
! Temporary space for computing 2-e integrals
!--------------------------------------------------------------------------------
  REAL(DOUBLE)  :: R(N,CBra,CKet,Ltot+1)
  REAL(DOUBLE)  :: Zeta,Up,Px,Py,Pz
  REAL(DOUBLE)  :: Eta, Uq,Qx,Qy,Qz
  REAL(DOUBLE)  :: Cx,Cy,Cz,PQx,PQy,PQz
  REAL(DOUBLE)  :: r1xZpE,Rkk,Tx,Ty,Tz
  REAL(DOUBLE)  :: T,T1,T2,T3,TS,TwoT,EX,GR(0:20)
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  INTEGER       :: I,J,K,M,Ind,IG
  INTEGER       :: I0,I1,I2

  IF(N.EQ.0) RETURN

  Cx=DisBufB( 8)
  Cy=DisBufB( 9)
  Cz=DisBufB(10)
  DO I=1,CBra
    CB(I,1)=PrmBufB(6,I)
    CB(I,2)=PrmBufB(7,I)
    CB(I,3)=PrmBufB(8,I)
  ENDDO ! I,CBra

  Ind=0
  DO I=1,CKet
    I0=(I-1)*DB%MAXP
    DO J=1,N
      I1        = SB%SLPrm%I(J)+I0
      CK(J,I,1) = DB%PrmBuf%D(I1+5)
      CK(J,I,2) = DB%PrmBuf%D(I1+6)
      CK(J,I,3) = DB%PrmBuf%D(I1+7)
    END DO
    DO K=1,CBra
      Zeta   = PrmBufB(1,K)
      Px     = PrmBufB(2,K)
      Py     = PrmBufB(3,K)
      Pz     = PrmBufB(4,K)
      Up     = PrmBufB(5,K)
      DO J=1,N
        Ind    = Ind+1
        I1     = SB%SLPrm%I(J)+I0
        I2     = SB%SLDis%I(J)-4
        Eta    = DB%PrmBuf%D(I1)
        Qx     = DB%PrmBuf%D(I1+1)
        Qy     = DB%PrmBuf%D(I1+2)
        Qz     = DB%PrmBuf%D(I1+3)
        Uq     = DB%PrmBuf%D(I1+4)
        r1xZpE = 1.0D0/(Zeta+Eta)
        Rkk    = Up*Uq*DSQRT(r1xZpE)
        Tx     = (Zeta*Px+Eta*Qx)*r1xZpE
        Ty     = (Zeta*Py+Eta*Qy)*r1xZpE
        Tz     = (Zeta*Pz+Eta*Qz)*r1xZpE

!        IB%WR%D( 1,Ind) = Px-Cx
!        IB%WR%D( 2,Ind) = Qx-DB%DisBuf%D(I2+7)
!        IB%WR%D( 3,Ind) = Py-Cy
!        IB%WR%D( 4,Ind) = Qy-DB%DisBuf%D(I2+8)
!        IB%WR%D( 5,Ind) = Pz-Cz
!        IB%WR%D( 6,Ind) = Qz-DB%DisBuf%D(I2+9)
!        IB%WR%D( 7,Ind) = Tx-Px
!        IB%WR%D( 8,Ind) = Tx-Qx
!        IB%WR%D( 9,Ind) = Ty-Py
!        IB%WR%D(10,Ind) = Ty-Qy
!        IB%WR%D(11,Ind) = Tz-Pz
!        IB%WR%D(12,Ind) = Tz-Qz
!        IB%WZ%D( 1,Ind) = Half/Eta
!        IB%WZ%D( 2,Ind) = Half/Zeta
!        IB%WZ%D( 3,Ind) = Zeta*r1xZpE
!        IB%WZ%D( 4,Ind) = Eta*r1xZpE
!        IB%WZ%D( 5,Ind) = Half*r1xZpE

        IB%WR%D(Ind, 1) = Px-Cx
        IB%WR%D(Ind, 2) = Qx-DB%DisBuf%D(I2+7)
        IB%WR%D(Ind, 3) = Py-Cy
        IB%WR%D(Ind, 4) = Qy-DB%DisBuf%D(I2+8)
        IB%WR%D(Ind, 5) = Pz-Cz
        IB%WR%D(Ind, 6) = Qz-DB%DisBuf%D(I2+9)
        IB%WR%D(Ind, 7) = Tx-Px
        IB%WR%D(Ind, 8) = Tx-Qx
        IB%WR%D(Ind, 9) = Ty-Py
        IB%WR%D(Ind,10) = Ty-Qy
        IB%WR%D(Ind,11) = Tz-Pz
        IB%WR%D(Ind,12) = Tz-Qz
        IB%WZ%D(Ind, 1) = Half/Eta
        IB%WZ%D(Ind, 2) = Half/Zeta
        IB%WZ%D(Ind, 3) = Zeta*r1xZpE
        IB%WZ%D(Ind, 4) = Eta*r1xZpE
        IB%WZ%D(Ind, 5) = Half*r1xZpE

        PQx=Px-Qx
        PQy=Py-Qy
        PQz=Pz-Qz
        T=(PQx*PQx+PQy*PQy+PQz*PQz)*Eta*Zeta*r1xZpE
        IF (T<Gamma_Switch)THEN 
          TwoT=2.0D0*T
          GR(LTot)=GammaF(LTot,T) 
          EX      =EXP(-T)
!          EX      =EXPInv(T)
          DO IG=Ltot-1,0,-1
            GR(IG)=(TwoT*GR(IG+1)+EX)/FLOAT(IG*2+1)
          END DO
          DO IG=1,LTot+1
            R(J,K,I,IG)=Rkk*GR(IG-1)
          END DO
        ELSE
          T1=1.0D0/T
          TS=DSQRT(T1)
          DO IG=1,LTot+1
            R(J,K,I,IG)=Rkk*IB%GammaA%D(IG)*TS
            TS=TS*T1
          END DO
        END IF ! Switch
      END DO ! J,N
    END DO ! K,CBra
  END DO ! I,CKet
END SUBROUTINE RGen
