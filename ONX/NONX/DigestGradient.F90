  SUBROUTINE DigestGradient(N,NA,NB,NC,ND,L1,L2,L3,L4,IntSwitch,AtA,AtC,AtD, &
                            NBFA,RS,SB,DB,D,SubInd,NTmp,Dcd,Dab,XFrc,W)
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  IMPLICIT NONE
  TYPE(DSL),INTENT(IN)         :: SB
  TYPE(DBuf),INTENT(IN)        :: DB
  TYPE(BCSR),INTENT(IN)        :: D
  TYPE(INT_RNK2),INTENT(IN)    :: SubInd
  INTEGER,INTENT(IN)           :: IntSwitch
  INTEGER,INTENT(IN)           :: L1,L2,L3,L4,NA,NB,NC,ND,N
  INTEGER,INTENT(IN)           :: AtA,AtC,AtD,RS
  INTEGER                      :: IA,IB,IC,ID
  INTEGER                      :: I,I1,Ind,Ioff
  INTEGER                      :: IndexA,NBFA
  INTEGER                      :: IndexB,AtB,NBFB,CS
  TYPE(INT_VECT),INTENT(INOUT) :: NTmp
  TYPE(DBL_RNK2),INTENT(INOUT) :: XFrc
  REAL(DOUBLE),INTENT(IN)      :: W(N,L2,L1,L4,L3,9)
  REAL(DOUBLE),INTENT(IN)      :: Dcd(NC,ND)
  REAL(DOUBLE),INTENT(IN)      :: Dab(N,NA,NB)
  REAL(DOUBLE)                 :: Dabcd
  REAL(DOUBLE)                 :: FAx,FAy,FAz
  REAL(DOUBLE)                 :: FBx,FBy,FBz
  REAL(DOUBLE)                 :: FCx,FCy,FCz
  REAL(DOUBLE)                 :: FDx,FDy,FDz

  DO I=1,N
    I1        = SB%SLDis%I(I)-3
    IndexB    = ABS(INT(DB%DisBuf%D(I1)))
    AtB       = SubInd%I(1,IndexB)
    NBFB      = SubInd%I(2,IndexB)
    CS        = SubInd%I(3,IndexB)
    NTmp%I(I) = AtB
    CALL GetAdrB(AtA,AtB,Ind,D,0)
    Ioff = D%BlkPt%I(Ind)
    CALL GetSubBlk(NBFA,NBFB,NA,NB,RS,CS,D%MTrix%D(Ioff),Dab(I,1,1))
  END DO

  write(*,*) "IntSwitch=",IntSwitch

  SELECT CASE (IntSwitch)
  CASE(11)
    DO IB=1,L3
    DO ID=1,L4
    DO IA=1,L1
    DO IC=1,L2
    DO I=1,N
      AtB = NTmp%I(I)
      Dabcd = -Dcd(IC,ID)*Dab(I,IA,IB)
      FAx=Dabcd*W(I,IC,IA,ID,IB,1)

    write(*,*) "Dabcd= ",Dabcd
    write(*,*) "W=",W(I,IC,IA,ID,IB,1)
!    call halt('enough')

      FAy=Dabcd*W(I,IC,IA,ID,IB,2)
      FAz=Dabcd*W(I,IC,IA,ID,IB,3)
      FCx=Dabcd*W(I,IC,IA,ID,IB,4)
      FCy=Dabcd*W(I,IC,IA,ID,IB,5)
      FCz=Dabcd*W(I,IC,IA,ID,IB,6)
      FBx=Dabcd*W(I,IC,IA,ID,IB,7)
      FBy=Dabcd*W(I,IC,IA,ID,IB,8)
      FBz=Dabcd*W(I,IC,IA,ID,IB,9)
      FDx=FAx+FCx+FBx                ! translational invariance
      FDy=FAy+FCy+FBy
      FDz=FAz+FCz+FBz
      XFrc%D(1,AtA)=XFrc%D(1,AtA)+FAx
      XFrc%D(2,AtA)=XFrc%D(2,AtA)+FAy
      XFrc%D(3,AtA)=XFrc%D(3,AtA)+FAz
      XFrc%D(1,AtC)=XFrc%D(1,AtC)+FCx
      XFrc%D(2,AtC)=XFrc%D(2,AtC)+FCy
      XFrc%D(3,AtC)=XFrc%D(3,AtC)+FCz
      XFrc%D(1,AtB)=XFrc%D(1,AtB)+FBx
      XFrc%D(2,AtB)=XFrc%D(2,AtB)+FBy
      XFrc%D(3,AtB)=XFrc%D(3,AtB)+FBz
      XFrc%D(1,AtD)=XFrc%D(1,AtD)-FDx
      XFrc%D(2,AtD)=XFrc%D(2,AtD)-FDy
      XFrc%D(3,AtD)=XFrc%D(3,AtD)-FDz
    END DO ! I
    END DO ! IC
    END DO ! IA
    END DO ! ID
    END DO ! IB
  CASE(01)
    DO IB=1,L3
    DO ID=1,L4
    DO IC=1,L1
    DO IA=1,L2
    DO I=1,N
      AtB = NTmp%I(I)
      Dabcd = -Dcd(IC,ID)*Dab(I,IA,IB)
      FAx=Dabcd*W(I,IA,IC,ID,IB,1)
      FAy=Dabcd*W(I,IA,IC,ID,IB,2)
      FAz=Dabcd*W(I,IA,IC,ID,IB,3)
      FCx=Dabcd*W(I,IA,IC,ID,IB,4)
      FCy=Dabcd*W(I,IA,IC,ID,IB,5)
      FCz=Dabcd*W(I,IA,IC,ID,IB,6)
      FBx=Dabcd*W(I,IA,IC,ID,IB,7)
      FBy=Dabcd*W(I,IA,IC,ID,IB,8)
      FBz=Dabcd*W(I,IA,IC,ID,IB,9)
      FDx=FAx+FCx+FBx                ! translational invariance
      FDy=FAy+FCy+FBy
      FDz=FAz+FCz+FBz
      XFrc%D(1,AtA)=XFrc%D(1,AtA)+FAx
      XFrc%D(2,AtA)=XFrc%D(2,AtA)+FAy
      XFrc%D(3,AtA)=XFrc%D(3,AtA)+FAz
      XFrc%D(1,AtC)=XFrc%D(1,AtC)+FCx
      XFrc%D(2,AtC)=XFrc%D(2,AtC)+FCy
      XFrc%D(3,AtC)=XFrc%D(3,AtC)+FCz
      XFrc%D(1,AtB)=XFrc%D(1,AtB)+FBx
      XFrc%D(2,AtB)=XFrc%D(2,AtB)+FBy
      XFrc%D(3,AtB)=XFrc%D(3,AtB)+FBz
      XFrc%D(1,AtD)=XFrc%D(1,AtD)-FDx
      XFrc%D(2,AtD)=XFrc%D(2,AtD)-FDy
      XFrc%D(3,AtD)=XFrc%D(3,AtD)-FDz
    END DO ! I
    END DO ! IA
    END DO ! IC
    END DO ! ID
    END DO ! IB
  CASE(10)
    DO ID=1,L3
    DO IB=1,L4
    DO IA=1,L1
    DO IC=1,L2
    DO I=1,N
      AtB = NTmp%I(I)
      Dabcd = -Dcd(IC,ID)*Dab(I,IA,IB)
      FAx=Dabcd*W(I,IC,IA,IB,ID,1)
      FAy=Dabcd*W(I,IC,IA,IB,ID,2)
      FAz=Dabcd*W(I,IC,IA,IB,ID,3)
      FCx=Dabcd*W(I,IC,IA,IB,ID,4)
      FCy=Dabcd*W(I,IC,IA,IB,ID,5)
      FCz=Dabcd*W(I,IC,IA,IB,ID,6)
      FBx=Dabcd*W(I,IC,IA,IB,ID,7)
      FBy=Dabcd*W(I,IC,IA,IB,ID,8)
      FBz=Dabcd*W(I,IC,IA,IB,ID,9)
      FDx=FAx+FCx+FBx                ! translational invariance
      FDy=FAy+FCy+FBy
      FDz=FAz+FCz+FBz
      XFrc%D(1,AtA)=XFrc%D(1,AtA)+FAx
      XFrc%D(2,AtA)=XFrc%D(2,AtA)+FAy
      XFrc%D(3,AtA)=XFrc%D(3,AtA)+FAz
      XFrc%D(1,AtC)=XFrc%D(1,AtC)+FCx
      XFrc%D(2,AtC)=XFrc%D(2,AtC)+FCy
      XFrc%D(3,AtC)=XFrc%D(3,AtC)+FCz
      XFrc%D(1,AtB)=XFrc%D(1,AtB)+FBx
      XFrc%D(2,AtB)=XFrc%D(2,AtB)+FBy
      XFrc%D(3,AtB)=XFrc%D(3,AtB)+FBz
      XFrc%D(1,AtD)=XFrc%D(1,AtD)-FDx
      XFrc%D(2,AtD)=XFrc%D(2,AtD)-FDy
      XFrc%D(3,AtD)=XFrc%D(3,AtD)-FDz
    END DO ! I
    END DO ! IC
    END DO ! IA
    END DO ! IB
    END DO ! ID
  CASE(00)
    DO ID=1,L3
    DO IB=1,L4
    DO IC=1,L1
    DO IA=1,L2
    DO I=1,N
      AtB = NTmp%I(I)
      Dabcd = -Dcd(IC,ID)*Dab(I,IA,IB)
      FAx=Dabcd*W(I,IA,IC,IB,ID,1)
      FAy=Dabcd*W(I,IA,IC,IB,ID,2)
      FAz=Dabcd*W(I,IA,IC,IB,ID,3)
      FCx=Dabcd*W(I,IA,IC,IB,ID,4)
      FCy=Dabcd*W(I,IA,IC,IB,ID,5)
      FCz=Dabcd*W(I,IA,IC,IB,ID,6)
      FBx=Dabcd*W(I,IA,IC,IB,ID,7)
      FBy=Dabcd*W(I,IA,IC,IB,ID,8)
      FBz=Dabcd*W(I,IA,IC,IB,ID,9)
      FDx=FAx+FCx+FBx                ! translational invariance
      FDy=FAy+FCy+FBy
      FDz=FAz+FCz+FBz
      XFrc%D(1,AtA)=XFrc%D(1,AtA)+FAx
      XFrc%D(2,AtA)=XFrc%D(2,AtA)+FAy
      XFrc%D(3,AtA)=XFrc%D(3,AtA)+FAz
      XFrc%D(1,AtC)=XFrc%D(1,AtC)+FCx
      XFrc%D(2,AtC)=XFrc%D(2,AtC)+FCy
      XFrc%D(3,AtC)=XFrc%D(3,AtC)+FCz
      XFrc%D(1,AtB)=XFrc%D(1,AtB)+FBx
      XFrc%D(2,AtB)=XFrc%D(2,AtB)+FBy
      XFrc%D(3,AtB)=XFrc%D(3,AtB)+FBz
      XFrc%D(1,AtD)=XFrc%D(1,AtD)-FDx
      XFrc%D(2,AtD)=XFrc%D(2,AtD)-FDy
      XFrc%D(3,AtD)=XFrc%D(3,AtD)-FDz
    END DO ! I
    END DO ! IA
    END DO ! IC
    END DO ! IB
    END DO ! ID
  CASE DEFAULT
    CALL Halt(' Illegal IntSwitch in XForce:DigestGradient')
  END SELECT
  END SUBROUTINE DigestGradient
