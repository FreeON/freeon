!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
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
  REAL(DOUBLE),DIMENSION(HGLen*2):: Ket !<<< SPIN We may need at most 2 rhos
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
       IMPLICIT NONE
       TYPE(CubeNode), POINTER               :: Cube
       REAL(DOUBLE), DIMENSION(0:HGEll+1)    :: LambdaX,LambdaY,LambdaZ
       REAL(DOUBLE)                          :: Z,X,RPx,RPy,RPz,RP2,Xpt,TwoZ,  &
                                                Co,RL1,Dist, Wght,dEdRho,    &
                                                dEdAbsGradRho2,GradBraRhoDot,&
                                                GradRhoX,GradRhoY,GradRhoZ,  &
                                                PrimDist,GradPrimDistX,      &
                                                GradPrimDistY,GradPrimDistZ, &
                                                dedrhoa,dedrhob,dedabsgradrho2aa,&
                                                dedabsgradrho2bb,dedabsgradrho2ab,&
                                                gradrhoxa,gradrhoya,gradrhoza,&
                                                gradrhoxb,gradrhoyb,gradrhozb
       INTEGER                               :: IC,IQ,JC,JQ,KC,KQ,NQ,Ell,    &
                                                I,J,L,M,N,L1,L2,LMN,LMNLen

       REAL(DOUBLE), DIMENSION(50) :: Gx,Gy,Gz
       REAL(DOUBLE)                :: GradX,GradY,GradZ,GRhoAG,GRhoBG
       INTEGER                     :: IA,IB,EllA,EllB,LMNA,LMNB,OffSDen,I1,I2,LMN1,LMN2,SDBeg,I0,iSDen

!
       SDBeg=0
       IF(NSDen.GT.1)SDBeg=1!<<< SPIN We don't need to integrate the rho_tot part.
       OffSDen=LHGTF(Prim%Ell)

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
          !
          IF(NSDen.EQ.1)THEN
             DO L=0,Prim%Ell
             DO M=0,Prim%Ell-L
             DO N=0,Prim%Ell-M-L
                PrimDist     = LambdaX(L  )*LambdaY(M  )*LambdaZ(N  )*Xpt
                GradPrimDistX=-LambdaX(L+1)*LambdaY(M  )*LambdaZ(N  )*Xpt
                GradPrimDistY=-LambdaX(L  )*LambdaY(M+1)*LambdaZ(N  )*Xpt
                GradPrimDistZ=-LambdaX(L  )*LambdaY(M  )*LambdaZ(N+1)*Xpt
                ! Kxc(LMN)=Kxc(LMN)+w_i*( dE/dRho*Phi_A*Phi_B +2*dE/d(GradRho)^2* GradRho.Grad(Phi_A*Phi_B) )
                LMN=LMNDex(L,M,N)
                Ket(LMN)=Ket(LMN)+Cube%Wght(I)                   &
                                 *(Cube%Vals(I,1)*PrimDist       &
                                 + Cube%Vals(I,2)*(              &
                                   Cube%Vals(I,3)*GradPrimDistX  &
                                 + Cube%Vals(I,4)*GradPrimDistY  &
                                 + Cube%Vals(I,5)*GradPrimDistZ))
             ENDDO
             ENDDO
             ENDDO
          ELSEIF(NSDen.EQ.3)THEN
             DO L=0,Prim%Ell
             DO M=0,Prim%Ell-L
             DO N=0,Prim%Ell-M-L
                PrimDist     = LambdaX(L  )*LambdaY(M  )*LambdaZ(N  )*Xpt
                GradPrimDistX=-LambdaX(L+1)*LambdaY(M  )*LambdaZ(N  )*Xpt
                GradPrimDistY=-LambdaX(L  )*LambdaY(M+1)*LambdaZ(N  )*Xpt
                GradPrimDistZ=-LambdaX(L  )*LambdaY(M  )*LambdaZ(N+1)*Xpt
                ! dEdRho(2*NGrid)
                ! dEdAbsGradRho2(3*NGrid)
                LMN1=LMNDex(L,M,N)
                LMN2=LMN1+OffSDen
                I1=I+  NGrid
                I2=I+2*NGrid
                !GradRho_a.Grad(Phi_A*Phi_B)
                GRhoAG=Cube%Vals(I ,3)*GradPrimDistX+Cube%Vals(I ,4)*GradPrimDistY+Cube%Vals(I ,5)*GradPrimDistZ
                !GradRho_b.Grad(Phi_A*Phi_B)
                GRhoBG=Cube%Vals(I1,3)*GradPrimDistX+Cube%Vals(I1,4)*GradPrimDistY+Cube%Vals(I1,5)*GradPrimDistZ
                !Kxc_a(LMN)=Kxc_a(LMN)+w_i*( dE/dRho_a*Phi_A*Phi_B
                !           +2*dE/d(GradRho_aa)^2*GradRho_a.Grad(Phi_A*Phi_B)
                !           +  dE/d(GradRho_ab)^2*GradRho_b.Grad(Phi_A*Phi_B) )
                Ket(LMN1)=Ket(LMN1)+Cube%Wght(I)*( &
                     & Cube%Vals(I ,1)*PrimDist + 2D0*Cube%Vals(I ,2)*GRhoAG + Cube%Vals(I2,2)*GRhoBG )
                !Kxc_b(LMN)=Kxc_b(LMN)+w_i*( dE/dRho_b*Phi_A*Phi_B
                !           +2*dE/d(GradRho_bb)^2*GradRho_b.Grad(Phi_A*Phi_B)
                !           +  dE/d(GradRho_ab)^2*GradRho_a.Grad(Phi_A*Phi_B) )
                Ket(LMN2)=Ket(LMN2)+Cube%Wght(I)*( &
                     & Cube%Vals(I1,1)*PrimDist + 2D0*Cube%Vals(I1,2)*GRhoBG + Cube%Vals(I2,2)*GRhoAG )
             ENDDO
             ENDDO
             ENDDO
          ELSE
             CALL Halt('Wrong NSDen in KxcCube!')
          ENDIF
          ENDIF
       ENDDO
#endif
     END SUBROUTINE KxcCube
END MODULE
