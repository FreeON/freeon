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
  USE CubeTree 
  IMPLICIT NONE
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
      INTEGER                :: AtA,AtB,JP,K,NA,NB,NAB,P,Q,R
#ifdef PERIODIC        
      INTEGER                :: NCA,NCB
      REAL(DOUBLE)           :: Ax,Ay,Az,Bx,By,Bz
      LOGICAL                :: NotMakeBlock
#endif     
!---------------------------------------------- 
!        Allocations 
!
      CALL New(MD,(/3,BS%NASym,BS%NASym,2*BS%NASym/),(/1,0,0,0/))
      CALL New(Kxc)
!
!        Initialize the matrix and associated indecies
!
      P=1; 
      R=1; 
      Kxc%RowPt%I(1)=1
      CALL SetEq(Kxc%MTrix,Zero)
!
!        Loop over atom pairs
!
      Kxc%NAtms=NAtoms
!
      DO AtA=1,NAtoms
         DO AtB=1,NAtoms
            IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
               NAB = Pair%NA*Pair%NB
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
!
!     Fill the upper triangle of Kxc
!
!!$      DO AtA=1,NAtoms
!!$         DO JP=Kxc%RowPt%I(AtA),Kxc%RowPt%I(AtA+1)-1
!!$            AtB=Kxc%ColPt%I(JP)
!!$            IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
!!$               NA  = Pair%NA
!!$               NB  = Pair%NB
!!$               NAB = Pair%NA*Pair%NB
!!$               WRITE(*,*) 'AtA = ',AtA,' AtB = ',AtB
!!$               IF(AtA>AtB)THEN
!!$                  DO K=Kxc%RowPt%I(AtB),Kxc%RowPt%I(AtB+1)-1
!!$                     IF(Kxc%ColPt%I(K)==AtA)THEN
!!$                        Q=Kxc%BlkPt%I(K)                     
!!$                        EXIT
!!$                     ENDIF
!!$                  ENDDO
!!$                  P=Kxc%BlkPt%I(JP)
!!$                  CALL XPose(NA,NB,Kxc%MTrix%D(P:P+NAB-1),Kxc%MTrix%D(Q:Q+NAB-1))
!!$               ENDIF
!!$            ENDIF
!!$         ENDDO
!!$      ENDDO
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
       TYPE(PrimPair)                           :: Prim
       REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)  :: Kblk
       REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB)  :: Kvct
       REAL(DOUBLE)                             :: ZetaA,ZetaB,EtaAB,EtaIn,    &
                                                   XiAB,ExpAB,CA,CB,CC,Ov,     &
                                                   PAx,PAy,PAz,PBx,PBy,PBz,    &
                                                   MDx,MDxy,MDxyz,Amp2,MaxAmp
       INTEGER                                  :: KA,KB
       INTEGER                                  :: CFA,CFB,PFA,PFB,  I,J,    &
                                                   IndexA,IndexB,        &
                                                   StartLA,StartLB,      &
                                                   StopLA,StopLB,        &
                                                   MaxLA,MaxLB,IA,IB,    &
                                                   LMNA,LMNB,LA,LB,MA,MB,&
                                                   NA,NB,LAB,MAB,NAB,LMN
!-------------------------------------------------------------------------------------- 
       KBlk=Zero
       KA=Pair%KA
       KB=Pair%KB
!
       IndexA=0                  
       DO CFA=1,BS%NCFnc%I(KA)    
          IndexA  = IBloDex(BS,CFA,KA)                
          StartLA = BS%LStrt%I(CFA,KA)        
          StopLA  = BS%LStop%I(CFA,KA)
          MaxLA   = BS%ASymm%I(2,CFA,KA)
          DO CFB=1,BS%NCFnc%I(KB) 
             IndexB  = IBloDex(BS,CFB,KB)                  
             StartLB = BS%LStrt%I(CFB,KB)
             StopLB  = BS%LStop%I(CFB,KB)
             MaxLB   = BS%ASymm%I(2,CFB,KB)
             DO PFA=1,BS%NPFnc%I(CFA,KA)          
                DO PFB=1,BS%NPFnc%I(CFB,KB)
!
!                  Set primitive values
!
                   ZetaA=BS%Expnt%D(PFA,CFA,KA)
                   ZetaB=BS%Expnt%D(PFB,CFB,KB)
                   EtaAB=ZetaA+ZetaB 
                   EtaIn=One/EtaAB
                   XiAB =ZetaA*ZetaB*EtaIn
                      
                   IF(TestPrimPair(XiAB,Pair%AB2))THEN
!
                      ExpAB=EXP(-XiAB*Pair%AB2)
                      Prim%P(1)=(ZetaA*Pair%A(1)+ZetaB*Pair%B(1))*EtaIn
                      Prim%P(2)=(ZetaA*Pair%A(2)+ZetaB*Pair%B(2))*EtaIn
                      Prim%P(3)=(ZetaA*Pair%A(3)+ZetaB*Pair%B(3))*EtaIn
!
                      Prim%Ell=MaxLA+MaxLB
                      Prim%Z=EtaAB
                      PAx=(Prim%P(1)-Pair%A(1))
                      PAy=(Prim%P(2)-Pair%A(2))
                      PAz=(Prim%P(3)-Pair%A(3))
                      PBx=(Prim%P(1)-Pair%B(1))
                      PBy=(Prim%P(2)-Pair%B(2))
                      PBz=(Prim%P(3)-Pair%B(3))
                      CALL MD2TRR(BS%NASym,0,MaxLA,MaxLB,EtaAB,MD%D, &
                              PAx,PBx,PAy,PBy,PAz,PBz) 
!
!                     Find the max amplitude for this primitive distribution
!
                      MaxAmp=Zero
                      DO LMNB=StartLB,StopLB
                         LB=BS%LxDex%I(LMNB)
                         MB=BS%LyDex%I(LMNB)
                         NB=BS%LzDex%I(LMNB)
                         CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                         DO LMNA=StartLA,StopLA
                            LA=BS%LxDex%I(LMNA)
                            MA=BS%LyDex%I(LMNA)
                            NA=BS%LzDex%I(LMNA)
                            CC=ExpAB*CB*BS%CCoef%D(LMNA,PFA,CFA,KA)
                            Amp2=Zero
                            DO LAB=0,LA+LB
                               MDx=CC*MD%D(1,LA,LB,LAB)
                               DO MAB=0,MA+MB
                                  MDxy=MDx*MD%D(2,MA,MB,MAB)
                                  DO NAB=0,NA+NB
                                     LMN=LMNdex(LAB,MAB,NAB)
                                     MDxyz=MDxy*MD%D(3,NA,NB,NAB)
                                     Amp2=Amp2+MDxyz**2
                                  ENDDO
                               ENDDO
                            ENDDO
                            MaxAmp=MAX(MaxAmp,SQRT(Amp2))
                         ENDDO
                      ENDDO
                      MaxAmp=LOG(MaxAmp)
!
!                     Evaluate this primitives contribution to Kxc_ab
!
                      Prim%Bra=Zero
                      Prim%Box%BndBox(1,:)=Prim%P(1)
                      Prim%Box%BndBox(2,:)=Prim%P(2)
                      Prim%Box%BndBox(3,:)=Prim%P(3)
                      Prim%Extent=GaussianExtent(Prim%Z,MaxAmp)
                      Prim%Box=ExpandBox(Prim%Box,Prim%Extent)
!
                      CALL PrimKxcWalk(CubeRoot,Prim)
!
!                     Compute matrix elements of Kxc
!
                      IB=IndexB
                      DO LMNB=StartLB,StopLB
                         LB=BS%LxDex%I(LMNB)
                         MB=BS%LyDex%I(LMNB)
                         NB=BS%LzDex%I(LMNB)
                         CB=BS%CCoef%D(LMNB,PFB,CFB,KB)
                         IB=IB+1
                         IA=IndexA
                         DO LMNA=StartLA,StopLA
                            IA=IA+1
                            LA=BS%LxDex%I(LMNA)
                            MA=BS%LyDex%I(LMNA)
                            NA=BS%LzDex%I(LMNA)
                            CC=ExpAB*CB*BS%CCoef%D(LMNA,PFA,CFA,KA)
                            DO LAB=0,LA+LB
                               MDx=CC*MD%D(1,LA,LB,LAB)
                               DO MAB=0,MA+MB
                                  MDxy=MDx*MD%D(2,MA,MB,MAB)
                                  DO NAB=0,NA+NB
                                     LMN=LMNdex(LAB,MAB,NAB)
                                     MDxyz=MDxy*MD%D(3,NA,NB,NAB)
                                     Kblk(IA,IB)=Kblk(IA,IB)+MDxyz*Prim%Bra(LMN) &
                                          *(-1.0D0)**(LAB+MAB+NAB)
                                  ENDDO
                               ENDDO
                            ENDDO
                         ENDDO
                      ENDDO
!
                   ENDIF
!
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       Kvct=BlockToVect(Pair%NA,Pair%NB,Kblk)
     END FUNCTION KxcBlock
!====================================================================================================
!
!====================================================================================================
     RECURSIVE SUBROUTINE PrimKxcWalk(Cube,Prim)
       TYPE(CubeNode), POINTER   :: Cube
       TYPE(PrimPair)            :: Prim
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
       IF(Cube%Leaf)THEN
          CALL KxcCube(Cube,Prim)
       ELSE
          CALL PrimKxcWalk(Cube%Descend,Prim)
          CALL PrimKxcWalk(Cube%Descend%Travrse,Prim)
       ENDIF
     END SUBROUTINE PrimKxcWalk
!====================================================================================================
!
!====================================================================================================
     SUBROUTINE KxcCube(Cube,Prim)
       TYPE(CubeNode), POINTER               :: Cube
       TYPE(PrimPair)                        :: Prim
! 
       REAL(DOUBLE), DIMENSION(0:MaxEll+1)   :: LambdaX,LambdaY,LambdaZ
       REAL(DOUBLE)                          :: Z,PRx,PRy,PRz,PR2,Xpt,TwoZ,  &
                                                Co,RL1,Dist, Wght,dEdRho,    &
                                                dEdAbsGradRho2,GradBraRhoDot,&
                                                GradRhoX,GradRhoY,GradRhoZ,  &
                                                PrimDist,GradPrimDistX,      &
                                                GradPrimDistY,GradPrimDistZ 
       INTEGER                               :: IC,IQ,JC,JQ,KC,KQ,NQ,Ell,    &
                                                I,L,M,N,L1,L2,LMN,LMNLen
!
       INCLUDE 'ExplicitBraElements.Inc'
       RETURN
       DO I=1,NGrid
          PRx=(Prim%P(1)-Cube%Grid(1,I))
          PRy=(Prim%P(2)-Cube%Grid(2,I))
          PRz=(Prim%P(3)-Cube%Grid(3,I))
          PR2=PRx**2+PRy**2+PRz**2
          Xpt=EXP(-Prim%Z*PR2)
          TwoZ=Two*Prim%Z
          LambdaX(0)=Xpt
          LambdaY(0)=One
          LambdaZ(0)=One
          LambdaX(1)=TwoZ*PRx*Xpt
          Lambday(1)=TwoZ*PRy
          LambdaZ(1)=TwoZ*PRz
          DO L=2,Prim%Ell+1
             L1=L-1
             L2=L-2
             RL1=DBLE(L1)
             LambdaX(L)=TwoZ*(PRx*LambdaX(L1)-RL1*LambdaX(L2))
             LambdaY(L)=TwoZ*(PRy*LambdaY(L1)-RL1*LambdaY(L2))
             LambdaZ(L)=TwoZ*(PRz*LambdaZ(L1)-RL1*LambdaZ(L2))
          ENDDO
          DO L=0,Prim%Ell
             DO M=0,Prim%Ell-L
                DO N=0,Prim%Ell-M-L
                   LMN=LMNDex(L,M,N)
                   PrimDist      =LambdaX(L  )*LambdaY(M  )*LambdaZ(N  )
                   GradPrimDistX=-LambdaX(L+1)*LambdaY(M  )*LambdaZ(N  )
                   GradPrimDistY=-LambdaX(L  )*LambdaY(M+1)*LambdaZ(N  )
                   GradPrimDistZ=-LambdaX(L  )*LambdaY(M  )*LambdaZ(N+1)
!
!   Prim_Kxc(LMN)=Prim_Kxc(LMN)+w_i*(  (dE/dRho) - 2(dE/d(GradRho)^2) Phi_a*Phi_B
!   GradRho.Grad(Phi_a*Phi_B)
!
                   Prim%Bra(LMN)=Prim%Bra(LMN)+Cube%Wght(I)   & 
                                *(Cube%Vals(I,1)*PrimDist       & 
                                -Two*Cube%Vals(I,2)             & 
                                *(Cube%Vals(I,3)*GradPrimDistX  & 
                                +Cube%Vals(I,4)*GradPrimDistY   & 
                                +Cube%Vals(I,5)*GradPrimDistZ))
                ENDDO
             ENDDO
          ENDDO
!
       ENDDO
     END SUBROUTINE KxcCube
!     
END MODULE KxcGen
