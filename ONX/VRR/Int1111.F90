SUBROUTINE Int1111(N,CBra,CKet,DisBufB,PrmBufB,DB,IB,SB,W1)
  USE DerivedTypes
  USE GlobalScalars
  USE ONXParameters
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
  REAL(DOUBLE)  :: W1(N)
!--------------------------------------------------------------------------------
! Temporary space for computing 2-e integrals
!--------------------------------------------------------------------------------
  REAL(DOUBLE)  :: Zeta,Up,PQx,PQy,PQz
  REAL(DOUBLE)  :: Eta, Uq,Qx,Qy,Qz
  REAL(DOUBLE)  :: r1xZpE,TwoE
  REAL(DOUBLE)  :: T1,T2,T3
  REAL(DOUBLE)  :: R1
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  INTEGER       :: I,J,K,M
  INTEGER       :: I0,I1

  IF(N.EQ.0) RETURN

  DO I=1,N
    I0 = SB%SLPrm%I(I)
    W1(I)=0.0D0
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
        Up   = PrmBufB(5,K)
        r1xZpE = 1.0D0/(Zeta+Eta)
        T1=(PQx*PQx+PQy*PQy+PQz*PQz)*Zeta*Eta*r1xZpE
        IF (T1<IB%Switch) THEN
          M=INT(T1*IB%Grid)
          T2=T1*T1
          T3=T2*T1
          R1=IB%GT%D(0,M)+T1*IB%GT%D(1,M)+T2*IB%GT%D(2,M)+T3*IB%GT%D(3,M)
        ELSE
          R1=IB%GammaA%D(1)/DSQRT(T1)
        ENDIF
        TwoE=TwoE+R1*Up*DSQRT(r1xZpE)
      END DO ! K
      W1(I)=W1(I)+TwoE*Uq
    END DO ! J
  END DO ! I
END SUBROUTINE Int1111
