MODULE KxcGen
  USE DerivedTypes
  USE GlobalScalars   
  USE GlobalObjects
  USE ProcessControl
  USE Indexing
  USE McMurchie
  USE InOut
  USE Thresholding
  USE AtomPairs
  USE BraKetBloks
  USE CubeTree 
  IMPLICIT NONE
  LOGICAL PrintFlag
!---------------------------------------
!
! Global primitive 
  TYPE(PrimPair)                           :: Prim

   TYPE(BSET)                    :: BS              !  Global basis set
   TYPE(CRDS)                    :: GM              !  Global molecular geometry
   TYPE(DBL_RNK4)                :: MD             


!----------!
  CONTAINS !
!=============================================================================================
!
!=============================================================================================
    SUBROUTINE MakeKxc(Kxc,CubeRoot)
      TYPE(BCSR)             :: Kxc
      TYPE(CubeNode),POINTER :: CubeRoot
!
      TYPE(AtomPair)         :: Pair
      INTEGER                :: AtA,AtB
      INTEGER                :: JP,K,NA,NB,NAB,P,Q,R,I,J
#ifdef PERIODIC        
      INTEGER                :: NCA,NCB
      REAL(DOUBLE)           :: Ax,Ay,Az,Bx,By,Bz
      LOGICAL                :: NotMakeBlock
#endif     
!---------------------------------------------- 
!     Allocations 
      CALL New(MD,(/3,BS%NASym,BS%NASym,2*BS%NASym/),(/1,0,0,0/))
      CALL New(Kxc)
!     Initialize the matrix and associated indecies
      P=1; 
      R=1; 
      Kxc%RowPt%I(1)=1
      CALL SetEq(Kxc%MTrix,Zero)
!     Loop over atom pairs
      Kxc%NAtms=NAtoms
      DO AtA=1,NAtoms
         DO AtB=1,NAtoms
            IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
               NAB = Pair%NA*Pair%NB
               IF(AtB<=AtA)THEN
!              Compute only the lower triangle of symmetric Kxc
#ifdef PERIODIC
                  Ax = Pair%A(1)
                  Ay = Pair%A(2)
                  Az = Pair%A(3)
                  Bx = Pair%B(1)
                  By = Pair%B(2)           
                  Bz = Pair%B(3)
                  NotMakeBlock = .TRUE.
                  DO NCA = 1,CS%NCells
                     Pair%A(1) = Ax+CS%CellCarts%D(1,NCA)
                     Pair%A(2) = Ay+CS%CellCarts%D(2,NCA)
                     Pair%A(3) = Az+CS%CellCarts%D(3,NCA)
                     DO NCB = 1,CS%NCells
                        Pair%B(1) = Bx+CS%CellCarts%D(1,NCB)
                        Pair%B(2) = By+CS%CellCarts%D(2,NCB)
                        Pair%B(3) = Bz+CS%CellCarts%D(3,NCB)
                        Pair%AB2  = (Pair%A(1)-Pair%B(1))**2 &
                                  + (Pair%A(2)-Pair%B(2))**2 &
                                  + (Pair%A(3)-Pair%B(3))**2
                        IF(TestAtomPair(Pair)) THEN
                           Kxc%MTrix%D(R:R+NAB-1)=Kxc%MTrix%D(R:R+NAB-1)+KxcBlock(Pair,CubeRoot)
                           NotMakeBlock = .FALSE.
                        ENDIF
                     ENDDO
                  ENDDO
#else
                  Kxc%MTrix%D(R:R+NAB-1)=KxcBlock(Pair,CubeRoot)
#endif
               ENDIF
               Kxc%ColPt%I(P)=AtB
               Kxc%BlkPt%I(P)=R
               R=R+NAB
               P=P+1 
               Kxc%RowPt%I(AtA+1)=P        
               IF(R>MaxNon0.OR.P>MaxBlks) &
                  CALL Halt(' BCSR dimensions blown in Kxc ')
#ifdef PERIODIC
               IF(NotMakeBlock) &
                  CALL Halt(' Making a Zero Block in KxcGen ')
#endif
            ENDIF
         ENDDO
      ENDDO
      Kxc%NBlks=P-1
      Kxc%NNon0=R-1
!     Fill the upper triangle of Kxc
      DO I=1,NAtoms
         DO JP=Kxc%RowPt%I(I),Kxc%RowPt%I(I+1)-1
            J=Kxc%ColPt%I(JP)
            IF(I>J)THEN
               DO K=Kxc%RowPt%I(J),Kxc%RowPt%I(J+1)-1
                  IF(Kxc%ColPt%I(K)==I)THEN
                     Q=Kxc%BlkPt%I(K)                     
                     EXIT
                  ENDIF
               ENDDO               
               P=Kxc%BlkPt%I(JP)
               NA=BS%BFKnd%I(GM%AtTyp%I(I))
               NB=BS%BFKnd%I(GM%AtTyp%I(J))
               NAB=NA*NB
               CALL XPose(NA,NB,Kxc%MTrix%D(P:P+NAB-1),Kxc%MTrix%D(Q:Q+NAB-1))
            ENDIF
         ENDDO
      ENDDO
!
    END SUBROUTINE MakeKxc
!=============================================================================================
!
!=============================================================================================
     SUBROUTINE XPose(M,N,A,AT)
       INTEGER                      :: I,J,M,N,IDex,JDex
       REAL(DOUBLE), DIMENSION(M*N) :: A,AT
       DO I=1,N
          DO J=1,M
             IDex=(I-1)*M+J
             JDex=(J-1)*N+I
             AT(JDex)=A(IDex)
          ENDDO
       ENDDO
     END SUBROUTINE XPose
!=======================================================================================
!
!=======================================================================================
     FUNCTION KxcBlock(Pair,CubeRoot) RESULT(Kvct)
       TYPE(AtomPair)                           :: Pair
       TYPE(CubeNode), POINTER                  :: CubeRoot
!
       REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)  :: Kblk
       REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB)  :: Kvct
       REAL(DOUBLE)                             :: ZetaA,ZetaB,EtaAB,EtaIn,    &
                                                   XiAB,ExpAB,CA,CB,CC,Ov,     &
                                                   PAx,PAy,PAz,PBx,PBy,PBz,    &
                                                   MDx,MDxy,MDxyz,Amp2,MaxAmp
       INTEGER                                  :: KA,KB,CFA,CFB,PFA,PFB,      &
                                                   IndexA,IndexB,              &
                                                   StartLA,StartLB,            &
                                                   StopLA,StopLB
       INTEGER                                  :: I,J,MaxLA,MaxLB,IA,IB,    &
                                                   LMNA,LMNB,LA,LB,MA,MB,    &
                                                   NA,NB,LAB,MAB,NAB,LMN,EllA,EllB
!-------------------------------------------------------------------------------------- 
       KBlk=Zero
       KA=Pair%KA
       KB=Pair%KB
!
       IndexA=0                  
       DO CFA=1,BS%NCFnc%I(KA)    
       DO CFB=1,BS%NCFnc%I(KB) 
          IndexA  = CFBlokDex(BS,CFA,KA)                
          IndexB  = CFBlokDex(BS,CFB,KB)                  
          StartLA = BS%LStrt%I(CFA,KA)        
          StartLB = BS%LStrt%I(CFB,KB)
          StopLA  = BS%LStop%I(CFA,KA)
          StopLB  = BS%LStop%I(CFB,KB)
          MaxLA   = BS%ASymm%I(2,CFA,KA)
          MaxLB   = BS%ASymm%I(2,CFB,KB)
          DO PFA=1,BS%NPFnc%I(CFA,KA)          
          DO PFB=1,BS%NPFnc%I(CFB,KB)
             ZetaA=BS%Expnt%D(PFA,CFA,KA)
             ZetaB=BS%Expnt%D(PFB,CFB,KB)
             EtaAB=ZetaA+ZetaB 
             EtaIn=One/EtaAB
             XiAB =ZetaA*ZetaB*EtaIn
!            Primitive thresholding                      
             IF(TestPrimPair(XiAB,Pair%AB2))THEN
!               Set primitive values
                ExpAB=EXP(-XiAB*Pair%AB2)
                Prim%P(1)=(ZetaA*Pair%A(1)+ZetaB*Pair%B(1))*EtaIn
                Prim%P(2)=(ZetaA*Pair%A(2)+ZetaB*Pair%B(2))*EtaIn
                Prim%P(3)=(ZetaA*Pair%A(3)+ZetaB*Pair%B(3))*EtaIn
                Prim%Ell=MaxLA+MaxLB
                Prim%Z=EtaAB
                PAx=Prim%P(1)-Pair%A(1)
                PAy=Prim%P(2)-Pair%A(2)
                PAz=Prim%P(3)-Pair%A(3)
                PBx=Prim%P(1)-Pair%B(1)
                PBy=Prim%P(2)-Pair%B(2)
                PBz=Prim%P(3)-Pair%B(3)
!               McMurchie Davidson prefactors for HG Primitives
                CALL MD2TRR(BS%NASym,0,MaxLA,MaxLB,EtaAB,MD%D,PAx,PBx,PAy,PBy,PAz,PBz) 
!               Primitive coefficients in a HG representationj
                CALL SetBraBlok(CFA,PFA,KA,CFB,PFB,KB,ExpAB,BS,MD)
!               Find the max amplitude for this primitive distribution
                MaxAmp=Zero
                IA = IndexA
                DO LMNA=StartLA,StopLA
                   IA=IA+1
                   IB=IndexB
                   EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)                         
                   DO LMNB=StartLB,StopLB
                      IB=IB+1
                      Amp2=Zero
                      EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)                         
                      DO LMN=1,LHGTF(EllA+EllB)
                         Amp2=Amp2+HGBra%D(LMN,IA,IB)**2
                      ENDDO
                      MaxAmp=MAX(MaxAmp,SQRT(Amp2))
                   ENDDO
                ENDDO
                MaxAmp=LOG(MaxAmp)
!               Evaluate this primitives Ket contribution to Kxc_ab
                Prim%Ket=Zero
                Prim%Box%BndBox(1,:)=Prim%P(1)
                Prim%Box%BndBox(2,:)=Prim%P(2)
                Prim%Box%BndBox(3,:)=Prim%P(3)
                Prim%Extent=GaussianExtent(Prim%Z,MaxAmp)
                Prim%Box=ExpandBox(Prim%Box,Prim%Extent)
!               Walk the CubeTree
                CALL PrimKxcWalk(CubeRoot)!,Prim)
!               Contract <Bra|Ket> bloks to compute matrix elements of Kxc
                IA = IndexA
                DO LMNA=StartLA,StopLA
                   IA=IA+1
                   IB=IndexB
                   EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)                         
                   DO LMNB=StartLB,StopLB
                      IB=IB+1
                      EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)                         
                      DO LMN=1,LHGTF(EllA+EllB)
                         Kblk(IA,IB)=Kblk(IA,IB)+HGBra%D(LMN,IA,IB)*Prim%Ket(LMN)
                      ENDDO
                  ENDDO
                  ENDDO
             ENDIF !End primitive thresholding
          ENDDO 
          ENDDO
       ENDDO
       ENDDO
       Kvct=BlockToVect(Pair%NA,Pair%NB,Kblk)
     END FUNCTION KxcBlock
!====================================================================================================
!
!====================================================================================================
     RECURSIVE SUBROUTINE PrimKxcWalk(Cube)!,Prim)
       TYPE(CubeNode), POINTER   :: Cube
!       TYPE(PrimPair)            :: Prim
!
       REAL(DOUBLE)              :: T
!-------------------------------------------------------------------------------------
!        Is the BBox containing this primitive distribution outside this nodes BBox? 
!
       T=ABS(Prim%Box%Center(1)-Cube%Box%Center(1))
       IF(T>Prim%Box%Half(1)+Cube%Box%Half(1))RETURN
       T=ABS(Prim%Box%Center(2)-Cube%Box%Center(2))
       IF(T>Prim%Box%Half(2)+Cube%Box%Half(2))RETURN
       T=ABS(Prim%Box%Center(3)-Cube%Box%Center(3))
       IF(T>Prim%Box%Half(3)+Cube%Box%Half(3))RETURN

!       IF(BoxInBox(Prim%Box,Cube%Box))THEN
          IF(Cube%Leaf)THEN
             CALL KxcCube(Cube)!,Prim)
          ELSE
             CALL PrimKxcWalk(Cube%Descend)!,Prim)
             CALL PrimKxcWalk(Cube%Descend%Travrse)!,Prim)
          ENDIF
!       ENDIF
     END SUBROUTINE PrimKxcWalk
!====================================================================================================
!
!====================================================================================================
     SUBROUTINE KxcCube(Cube)!,Prim)
       TYPE(CubeNode), POINTER               :: Cube
!       TYPE(PrimPair)                        :: Prim
! 
       REAL(DOUBLE), DIMENSION(0:MaxEll+1)   :: LambdaX,LambdaY,LambdaZ
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
          Dist=Z*RP2
          Xpt=EXP(-Prim%Z*RP2)
          IF(Xpt>1.D-14)THEN
          TwoZ=Two*Prim%Z
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
!                  Prim_Kxc(LMN)=Prim_Kxc(LMN)+w_i*(dE/dRho-2dE/d(GradRho)^2 Phi_a*Phi_B GradRho.Grad(Phi_a*Phi_B)
                   Prim%Ket(LMN)=Prim%Ket(LMN)+Cube%Wght(I)        &
                                   *(Cube%Vals(I,1)*PrimDist       & 
                                    +Cube%Vals(I,2)*(              &
                                     Cube%Vals(I,3)*GradPrimDistX  &
                                    +Cube%Vals(I,4)*GradPrimDistY  &
                                    +Cube%Vals(I,5)*GradPrimDistZ))
                ENDDO
             ENDDO
          ENDDO
          ENDIF
       ENDDO
#endif
     END SUBROUTINE KxcCube
END MODULE KxcGen
