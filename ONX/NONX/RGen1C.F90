SUBROUTINE RGen1C(LDis,iB,Kon,WK,CD,WR,WZ,R,TBufP,TBufC)
  USE DerivedTypes
  USE GlobalScalars
  USE ONXParameters
  IMPLICIT NONE
!--------------------------------------------------------------------------------
! Temporary space for computing 2-e integrals
!--------------------------------------------------------------------------------
  INTEGER                :: LDis,iB,Kon
  TYPE(DBL_RNK2)         :: WK,CD,TBufC
  TYPE(DBL_RNK3)         :: WR,WZ,TBufP
  REAL(DOUBLE)           :: R(Kon,Kon,2*LDis+1)
!--------------------------------------------------------------------------------
! Misc. internal variables
!--------------------------------------------------------------------------------
  REAL(DOUBLE)           :: Px,Py,Pz,Up,Zeta
  REAL(DOUBLE)           :: Qx,Qy,Qz,Uq,Eta
  REAL(DOUBLE)           :: Tx,Ty,Tz,r1xZpE,Rkk
  INTEGER                :: I,J,K

  DO I=1,Kon
    CD%D(I,1)=TBufP%D(6,I,iB)
    CD%D(I,2)=TBufP%D(7,I,iB)
    CD%D(I,3)=TBufP%D(8,I,iB)
    WK%D(1,I)=TBufP%D(2,I,iB)-TBufC%D( 8,iB)
    WK%D(2,I)=TBufP%D(3,I,iB)-TBufC%D( 9,iB)
    WK%D(3,I)=TBufP%D(4,I,iB)-TBufC%D(10,iB)
    WK%D(4,I)=0.5D0/TBufP%D(1,I,iB)
    Eta = TBufP%D(1,I,iB)
    Qx  = TBufP%D(2,I,iB)
    Qy  = TBufP%D(3,I,iB)
    Qz  = TBufP%D(4,I,iB)
    Uq  = TBufP%D(5,I,iB)
    DO J=1,Kon
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
      WR%D(1,J,I)=Tx-Px
      WR%D(2,J,I)=Ty-Py
      WR%D(3,J,I)=Tz-Pz
      WR%D(4,J,I)=Tx-Qx
      WR%D(5,J,I)=Ty-Qy
      WR%D(6,J,I)=Tz-Qz
      WZ%D(1,J,I)= Zeta*r1xZpE
      WZ%D(2,J,I)=  Eta*r1xZpE
      WZ%D(3,J,I)=0.5D0*r1xZpE
      DO K=0,2*LDis
        R(J,I,K+1)=Rkk/(1.0D0+2.0D0*DFLOAT(I))
      ENDDO
    END DO ! J
  END DO ! I

END SUBROUTINE RGen1C
 
