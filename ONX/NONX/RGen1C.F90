SUBROUTINE RGen1C(LDis,iB,Kon,CD,WR,WZ,R,TBufP,TBufC)
  USE DerivedTypes
  USE GlobalScalars
  USE ONXParameters
  IMPLICIT NONE
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
      Ind=Ind+1
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
      WR%D( 1,Ind) = Px-Cx
      WR%D( 2,Ind) = Qx-Cx
      WR%D( 3,Ind) = Py-Cy
      WR%D( 4,Ind) = Qy-Cy
      WR%D( 5,Ind) = Pz-Cz
      WR%D( 6,Ind) = Qz-Cz
      WR%D( 7,Ind) = Tx-Px
      WR%D( 8,Ind) = Tx-Qx
      WR%D( 9,Ind) = Ty-Py
      WR%D(10,Ind) = Ty-Qy
      WR%D(11,Ind) = Tz-Pz
      WR%D(12,Ind) = Tz-Qz
      WZ%D( 1,Ind) = Half/Eta
      WZ%D( 2,Ind) = Half/Zeta
      WZ%D( 3,Ind) = Zeta*r1xZpE
      WZ%D( 4,Ind) = Eta*r1xZpE
      WZ%D( 5,Ind) = Half*r1xZpE
      DO K=0,2*LDis
        R(J,I,K+1)=Rkk/(1.0D0+2.0D0*DFLOAT(K))
      ENDDO
    END DO ! J
  END DO ! I

END SUBROUTINE RGen1C
 
