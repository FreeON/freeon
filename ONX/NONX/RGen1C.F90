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

END SUBROUTINE RGen1C
 
