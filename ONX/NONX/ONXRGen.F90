MODULE ONXRGen
!H=================================================================================
!H MODULE ONXRGen
!H This MODULE contains:
!H  o SUB RGen1C
!H  o SUB RGen
!H---------------------------------------------------------------------------------
  USE DerivedTypes
  USE GlobalScalars
  USE ONXParameters
  USE GammaFunctions
  !
  IMPLICIT NONE
  PRIVATE
  !
!---------------------------------------------------------------------------------
! PUBLIC DECLARATIONS
!---------------------------------------------------------------------------------
  PUBLIC :: RGen1C
  PUBLIC :: RGen
  !
CONTAINS
  !
  SUBROUTINE RGen1C(LDis,iB,Kon,CD,WR,WZ,R,TBufP,TBufC)
!--------------------------------------------------------------------------------
! Temporary space for computing 2-e integrals
!--------------------------------------------------------------------------------
    INTEGER                :: LDis,iB,Kon
    TYPE(DBL_RNK2)         :: TBufC,WR,WZ
    TYPE(DBL_RNK3)         :: TBufP
    REAL(DOUBLE)           :: CD(Kon,3)
    REAL(DOUBLE)           :: R(Kon,Kon,2*LDis+1)
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
    REAL(DOUBLE)           :: Cx,Cy,Cz
    REAL(DOUBLE)           :: Px,Py,Pz,Up,Zeta
    REAL(DOUBLE)           :: Qx,Qy,Qz,Uq,Eta
    REAL(DOUBLE)           :: Tx,Ty,Tz,r1xZpE,Rkk
    INTEGER                :: I,J,K,Ind
    !
    Cx=TBufC%D( 8,iB)
    Cy=TBufC%D( 9,iB)
    Cz=TBufC%D(10,iB)
    Ind=0
    DO I=1,Kon
       CD(I,1)=TBufP%D(6,I,iB)
       CD(I,2)=TBufP%D(7,I,iB)
       CD(I,3)=TBufP%D(8,I,iB)
       Eta = TBufP%D(1,I,iB)
       Qx  = TBufP%D(2,I,iB)
       Qy  = TBufP%D(3,I,iB)
       Qz  = TBufP%D(4,I,iB)
       Uq  = TBufP%D(5,I,iB)
       DO J=1,Kon
          Ind    = Ind+1
          Zeta   = TBufP%D(1,J,iB)
          Px     = TBufP%D(2,J,iB)
          Py     = TBufP%D(3,J,iB)
          Pz     = TBufP%D(4,J,iB)
          Up     = TBufP%D(5,J,iB)
          r1xZpE = 1.0D0/(Zeta+Eta)
          Rkk    = Up*Uq*DSQRT(r1xZpE)
          Tx     = (Zeta*Px+Eta*Qx)*r1xZpE
          Ty     = (Zeta*Py+Eta*Qy)*r1xZpE
          Tz     = (Zeta*Pz+Eta*Qz)*r1xZpE
          WR%D(Ind, 1) = Px-Cx
          WR%D(Ind, 2) = Qx-Cx
          WR%D(Ind, 3) = Py-Cy
          WR%D(Ind, 4) = Qy-Cy
          WR%D(Ind, 5) = Pz-Cz
          WR%D(Ind, 6) = Qz-Cz
          WR%D(Ind, 7) = Tx-Px
          WR%D(Ind, 8) = Tx-Qx
          WR%D(Ind, 9) = Ty-Py
          WR%D(Ind,10) = Ty-Qy
          WR%D(Ind,11) = Tz-Pz
          WR%D(Ind,12) = Tz-Qz
          WZ%D(Ind, 1) = Half/Eta
          WZ%D(Ind, 2) = Half/Zeta
          WZ%D(Ind, 3) = Zeta*r1xZpE
          WZ%D(Ind, 4) = Eta*r1xZpE
          WZ%D(Ind, 5) = Half*r1xZpE
          DO K=0,2*LDis
             R(J,I,K+1)=Rkk/(1.0D0+2.0D0*DBLE(K))
          ENDDO
       END DO ! J
    END DO ! I
    ! 
  END SUBROUTINE RGen1C
  !
  !
  SUBROUTINE RGen(N,Ltot,CBra,CKet,CB,CK,DisBufB,PrmBufB,R,DB,IB,SB,PBC)      !per
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
!--------------------------------------------------------------------------------
! Periodic Variables
!--------------------------------------------------------------------------------
    TYPE(PBCInfo) :: PBC                 !per
    REAL(DOUBLE)  :: FPQx,FPQy,FPQz      !per

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

             PQx=Px-Qx                                                                           
             PQy=Py-Qy                                                                           
             PQz=Pz-Qz                                                                           
!                                                                                                !per
             FPQx = PQx*PBC%InvBoxSh%D(1,1)+PQy*PBC%InvBoxSh%D(1,2)+PQz*PBC%InvBoxSh%D(1,3)      !per
             FPQy = PQy*PBC%InvBoxSh%D(2,2)+PQz*PBC%InvBoxSh%D(2,3)                              !per
             FPQz = PQz*PBC%InvBoxSh%D(3,3)                                                      !per
             IF(PBC%AutoW%I(1)==1) FPQx = FPQx-ANINT(FPQx)                                       !per
             IF(PBC%AutoW%I(2)==1) FPQy = FPQy-ANINT(FPQy)                                       !per
             IF(PBC%AutoW%I(3)==1) FPQz = FPQz-ANINT(FPQz)                                       !per
             PQx  = FPQx*PBC%BoxShape%D(1,1)+FPQy*PBC%BoxShape%D(1,2)+FPQz*PBC%BoxShape%D(1,3)   !per
             PQy  = FPQy*PBC%BoxShape%D(2,2)+FPQz*PBC%BoxShape%D(2,3)                            !per
             PQz  = FPQz*PBC%BoxShape%D(3,3)                                                     !per
!                                                                                                !per
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
!old             PQx=Px-Qx
!old             PQy=Py-Qy
!old             PQz=Pz-Qz
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
!            R(J,K,I,IG)=Rkk*IB%GammaA%D(IG)*TS
                   R(J,K,I,IG)=Rkk*GammAss(IG-1)*TS
                   TS=TS*T1
                END DO
             END IF ! Switch
          END DO ! J,N
       END DO ! K,CBra
    END DO ! I,CKet
  END SUBROUTINE RGen
  !
END MODULE ONXRGen



