MODULE TreeWalk
  USE DerivedTypes
  USE GlobalScalars   
  USE GlobalObjects
  USE ProcessControl
  USE Indexing
  USE McMurchie
  USE InOut
  USE Thresholding
  USE AtomPairs
  USE BraBloks
  USE CubeTree 
  IMPLICIT NONE
  LOGICAL PrintFlag
!-----------------------------------------------------
! Globals
  REAL(DOUBLE),DIMENSION(HGLen)  :: Ket
  REAL(DOUBLE)                   :: PExtent
  TYPE(BBox)                     :: PBox
  TYPE(PrimPair)                 :: Prim
!----------!
  CONTAINS !
     RECURSIVE SUBROUTINE KxcWalk(Cube)
       TYPE(CubeNode), POINTER   :: Cube
!
       REAL(DOUBLE)              :: T
!-------------------------------------------------------------------------------------
!      Does PBox intersect this Cubes BBox? 
       T=ABS(PBox%Center(1)-Cube%Box%Center(1))
       IF(T>PBox%Half(1)+Cube%Box%Half(1))RETURN
       T=ABS(PBox%Center(2)-Cube%Box%Center(2))
       IF(T>PBox%Half(2)+Cube%Box%Half(2))RETURN
       T=ABS(PBox%Center(3)-Cube%Box%Center(3))
       IF(T>PBox%Half(3)+Cube%Box%Half(3))RETURN
       IF(Cube%Leaf)THEN
!         Evaluate primitives on this cube
          CALL KxcCube(Cube)
       ELSE
!         Keep on truckin      
          CALL KxcWalk(Cube%Descend)
          CALL KxcWalk(Cube%Descend%Travrse)
       ENDIF
     END SUBROUTINE KxcWalk
!====================================================================================================
!
!====================================================================================================
     SUBROUTINE KxcCube(Cube)
       TYPE(CubeNode), POINTER               :: Cube
       REAL(DOUBLE), DIMENSION(0:HGEll+1)    :: LambdaX,LambdaY,LambdaZ
       REAL(DOUBLE)                          :: Z,X,RPx,RPy,RPz,RP2,Xpt,TwoZ,  &
                                                Co,RL1,Dist, Wght,dEdRho,    &
                                                dEdAbsGradRho2,GradBraRhoDot,&
                                                GradRhoX,GradRhoY,GradRhoZ,  &
                                                PrimDist,GradPrimDistX,      &
                                                GradPrimDistY,GradPrimDistZ 
       INTEGER                               :: IC,IQ,JC,JQ,KC,KQ,NQ,Ell,    &
                                                I,J,L,M,N,L1,L2,LMN,LMNLen

       REAL(DOUBLE), DIMENSION(50) :: Gx,Gy,Gz
       REAL(DOUBLE)                :: GradX,GradY,GradZ
       INTEGER                     :: IA,IB,EllA,EllB,LMNA,LMNB

!
#ifdef EXPLICIT_SOURCE
       INCLUDE 'ExplicitBraElements.Inc'
#else
       DO I=1,NGrid
          RPx=Cube%Grid(I,1)-Prim%P(1)
          RPy=Cube%Grid(I,2)-Prim%P(2)
          RPz=Cube%Grid(I,3)-Prim%P(3)
          RP2=RPx**2+RPy**2+RPz**2
          Xpt=EXP(-Prim%Zeta*RP2)
          IF(Xpt>1.D-14)THEN
          TwoZ=Two*Prim%Zeta
          LambdaX(0)=One
          LambdaY(0)=One
          LambdaZ(0)=One
          LambdaX(1)=TwoZ*RPx 
          Lambday(1)=TwoZ*RPy
          LambdaZ(1)=TwoZ*RPz
          DO L=2,Prim%Ell+1
             L1=L-1
             L2=L-2
             RL1=DBLE(L1)
             LambdaX(L)=TwoZ*(RPx*LambdaX(L1)-RL1*LambdaX(L2))
             LambdaY(L)=TwoZ*(RPy*LambdaY(L1)-RL1*LambdaY(L2))
             LambdaZ(L)=TwoZ*(RPz*LambdaZ(L1)-RL1*LambdaZ(L2))
          ENDDO
          DO L=0,Prim%Ell
             DO M=0,Prim%Ell-L
                DO N=0,Prim%Ell-M-L
                   LMN=LMNDex(L,M,N)
                   PrimDist     = LambdaX(L  )*LambdaY(M  )*LambdaZ(N  )*Xpt
                   GradPrimDistX=-LambdaX(L+1)*LambdaY(M  )*LambdaZ(N  )*Xpt
                   GradPrimDistY=-LambdaX(L  )*LambdaY(M+1)*LambdaZ(N  )*Xpt
                   GradPrimDistZ=-LambdaX(L  )*LambdaY(M  )*LambdaZ(N+1)*Xpt
!                  Kxc(LMN)=Kxc(LMN)+w_i*(dE/dRho-2dE/d(GradRho)^2 Phi_a*Phi_B GradRho.Grad(Phi_a*Phi_B)
                   Ket(LMN)=Ket(LMN)+Cube%Wght(I)                  &
                                    *(Cube%Vals(I,1)*PrimDist       & 
                                    + Cube%Vals(I,2)*(              &
                                      Cube%Vals(I,3)*GradPrimDistX  &
                                    + Cube%Vals(I,4)*GradPrimDistY  &
                                    + Cube%Vals(I,5)*GradPrimDistZ))
!                   Ket(LMN)=Ket(LMN)+Cube%Wght(I)*Cube%Vals(I,1)*PrimDist
                ENDDO
             ENDDO
          ENDDO
          ENDIF
       ENDDO
#endif
     END SUBROUTINE KxcCube
END MODULE
