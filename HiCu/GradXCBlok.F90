MODULE GradXCBlock
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
!--------------------------------------------------------------------------
! Primative Atom Pair Type
!
  TYPE PrimPair
     INTEGER                        :: Ell  
     REAL(DOUBLE)                   :: Z
     REAL(DOUBLE),DIMENSION(3)      :: P
     REAL(DOUBLE),DIMENSION(dHGLen) :: Ket
     REAL(DOUBLE)                   :: Extent
     TYPE(BBox)                     :: Box
  ENDTYPE PrimPair
!---------------------------------------
! Global primitive 
!
  TYPE(PrimPair)                 :: Prim
!----------!
  CONTAINS !
!=======================================================================================
!
!=======================================================================================
     FUNCTION GradXCBlok(BS,MD,AtA,AtB,DiffAtm,P,Pair,CubeRoot) RESULT(dKvct)
       TYPE(AtomPair)                           :: Pair
       TYPE(BSET),                 INTENT(IN)    :: BS
       TYPE(DBL_RNK4),             INTENT(INOUT) :: MD
       TYPE(CubeNode), POINTER                  :: CubeRoot
!
!       TYPE(PrimPair)                           :: Prim
       REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)  :: P
       REAL(DOUBLE),DIMENSION(3)                :: dKvct
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
                                                   NA,NB,LAB,MAB,NAB,LMN,EllA,EllB, &
                                                   DiffAtm,AtA,AtB
!-------------------------------------------------------------------------------------- 
       dKvct=Zero
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
!               McMurchie Davidson prefactors for dervatives of HG Primitives
                CALL MD2TRR(BS%NASym+2,-1,MaxLA+1,MaxLB+1,EtaAB,MD%D, &
                     PAx,PBx,PAy,PBy,PAz,PBz) 
!               Compute the gradient of the Bra in <Bra|Ket>
                CALL dTwoE(DiffAtm,AtA,AtB,CFA,PFA,KA,CFB,PFB,KB,ExpAB,BS,MD,Pair)
!               Find the max amplitude for contribution involving this primitive distribution
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
                      DO LMN=1,LHGTF(EllA+EllB+1)
                         Amp2=Amp2+dHGBradX%D(LMN,IA,IB)**2 &
                                  +dHGBradY%D(LMN,IA,IB)**2 &
                                  +dHGBradZ%D(LMN,IA,IB)**2
                      ENDDO
                      MaxAmp=MAX(MaxAmp,SQRT(Amp2))
                   ENDDO
                ENDDO
                MaxAmp=LOG(MaxAmp)
!               Evaluate this primitives Ket contribution to Kxc_ab
!               Defined in /MondoMods/DerivedTypes.F90
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
                      DO LMN=1,LHGTF(EllA+EllB+1)
                         dKvct(1)=dKvct(1) + Two*P(IA,IB)*dHGBradX%D(LMN,IA,IB)*Prim%Ket(LMN)
                         dKvct(2)=dKvct(2) + Two*P(IA,IB)*dHGBradY%D(LMN,IA,IB)*Prim%Ket(LMN)
                         dKvct(3)=dKvct(3) + Two*P(IA,IB)*dHGBradZ%D(LMN,IA,IB)*Prim%Ket(LMN)
                      ENDDO
                  ENDDO
                  ENDDO
             ENDIF !End primitive thresholding
          ENDDO 
          ENDDO
       ENDDO
       ENDDO
!       Kvct=BlockToVect(Pair%NA,Pair%NB,Kblk)
     END FUNCTION GradXCBlok
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
       IF(Cube%Leaf)THEN
          CALL GradKxcCube(Cube)!,Prim)
       ELSE
          CALL PrimKxcWalk(Cube%Descend)!,Prim)
          CALL PrimKxcWalk(Cube%Descend%Travrse)!,Prim)
       ENDIF
     END SUBROUTINE PrimKxcWalk
!====================================================================================================
!
!====================================================================================================
     SUBROUTINE GradKxcCube(Cube)!,Prim)
       TYPE(CubeNode), POINTER               :: Cube
!       TYPE(PrimPair)                        :: Prim
! 
       REAL(DOUBLE), DIMENSION(0:HGEll+2)   :: LambdaX,LambdaY,LambdaZ
       REAL(DOUBLE)                          :: Z,RPx,RPy,RPz,RP2,Xpt,TwoZ,  &
                                                Co,RL1,Dist, Wght,dEdRho,    &
                                                dEdAbsGradRho2,GradBraRhoDot,&
                                                GradRhoX,GradRhoY,GradRhoZ,  &
                                                PrimDist,GradPrimDistX,      &
                                                GradPrimDistY,GradPrimDistZ 
       INTEGER                               :: IC,IQ,JC,JQ,KC,KQ,NQ,Ell,    &
                                                I,L,M,N,L1,L2,LMN,LMNLen

       REAL(DOUBLE)                :: GradX,GradY,GradZ
       INTEGER                     :: IA,IB,EllA,EllB,LMNA,LMNB

!
!#ifdef EXPLICIT_SOURCE
!       INCLUDE 'ExplicitBraElements.Inc'
!#else
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
          DO L=2,Prim%Ell+2
             L1=L-1
             L2=L-2
             RL1=DBLE(L1)
             LambdaX(L)=TwoZ*(RPx*LambdaX(L1)-RL1*LambdaX(L2))
             LambdaY(L)=TwoZ*(RPy*LambdaY(L1)-RL1*LambdaY(L2))
             LambdaZ(L)=TwoZ*(RPz*LambdaZ(L1)-RL1*LambdaZ(L2))
          ENDDO
          DO L=0,Prim%Ell+1
             DO M=0,Prim%Ell-L+1
                DO N=0,Prim%Ell-M-L+1
                   IF (L+M+N.LE.Prim%Ell+1) THEN
                      LMN=LMNDex(L,M,N)
                      PrimDist     = LambdaX(L  )*LambdaY(M  )*LambdaZ(N  )*Xpt
                      GradPrimDistX=-LambdaX(L+1)*LambdaY(M  )*LambdaZ(N  )*Xpt
                      GradPrimDistY=-LambdaX(L  )*LambdaY(M+1)*LambdaZ(N  )*Xpt
                      GradPrimDistZ=-LambdaX(L  )*LambdaY(M  )*LambdaZ(N+1)*Xpt
                      Prim%Ket(LMN)=Prim%Ket(LMN)+Cube%Wght(I)*(        &
                           Cube%Vals(I,1)*PrimDist       & 
                           +Cube%Vals(I,2)*(              &
                           Cube%Vals(I,3)*GradPrimDistX  &
                           +Cube%Vals(I,4)*GradPrimDistY  &
                           +Cube%Vals(I,5)*GradPrimDistZ))
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
          ENDIF
       ENDDO
!#endif
     END SUBROUTINE GradKxcCube
END MODULE GradXCBlock
