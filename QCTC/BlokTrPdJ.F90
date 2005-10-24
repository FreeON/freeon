!    FAST O(Lg N), BLOKWISE ACCUMULATION OF Tr{P_(A,B).dJ^T_(A,B)/dA}
!    Author: Matt Challacombe 
!------------------------------------------------------------------------------
MODULE BlokTrPdJ
  USE DerivedTypes
  USE GlobalScalars   
  USE GlobalObjects
  USE ProcessControl
  USE Indexing
  USE InOut
  USE Thresholding
  USE Globals
  USE AtomPairs
  USE BraBloks
  USE PoleTree
  USE TreeWalk
  USE PBCFarField
  IMPLICIT NONE
  LOGICAL PrintFlag
!----------!
  CONTAINS !
!=======================================================================================
!
!=======================================================================================
     FUNCTION TrPdJ(Pair,P,GMLoc) RESULT(Vck)
       TYPE(AtomPair)                             :: Pair
       TYPE(CRDS)                                 :: GMLoc
!
       REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)    :: P
       REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB,15) :: dJ
       REAL(DOUBLE),DIMENSION(3)                  :: PTmp,PShft
       REAL(DOUBLE),DIMENSION(0:SPLen)            :: SPBraC,SPBraS 
       REAL(DOUBLE)                               :: ZetaA,ZetaB,EtaAB,EtaIn,    &
                                                     XiAB,ExpAB,CA,CB,CC,Ov,     &
                                                     PAx,PAy,PAz,PBx,PBy,PBz,    &
                                                     MDx,MDxy,MDxyz,Amp2,MaxAmp, &
                                                     Pab,JNorm,Tau,OmegaMin,     &
                                                     PExtent,PStrength,Ext
       INTEGER                                    :: KA,KB,CFA,CFB,PFA,PFB,AtA,AtB,    &
                                                     IndexA,IndexB,              &
                                                     StartLA,StartLB,            &
                                                     StopLA,StopLB
       INTEGER                                    :: I,J,K,L,MaxLA,MaxLB,IA,IB,  &
                                                     LMNA,LMNB,LA,LB,MA,MB,      &
                                                     NA,NB,LAB,MAB,NAB,LM,LMN,   &
                                                     EllA,EllB,EllAB,LenHGTF,    &
                                                     LenSP,KI 
       REAL(DOUBLE), EXTERNAL                     :: BlkTrace_2 
       INTEGER                                    :: NC
       REAL(DOUBLE)                               :: Px,Py,Pz,HGFac,PiZ
!
       REAL(DOUBLE),DIMENSION(15)                 :: Vck
       REAL(DOUBLE),DIMENSION(3)                  :: nlm
       REAL(DOUBLE),DIMENSION(0:SPLen,3)          :: SPKetC_L,SPKetS_L
       REAL(DOUBLE),DIMENSION(1:HGLen,3)          :: HGKet_L
       REAL(DOUBLE),DIMENSION(1:HGLen)            :: HGKet_old
       REAL(DOUBLE),DIMENSION(0:SPLen)            :: SPKetC_old,SPKetS_old
!----------------------------------------------------------------------------------------- 
       Prim%A=Pair%A
       Prim%B=Pair%B
       Prim%KA=Pair%KA
       Prim%KB=Pair%KB
       Prim%AB2=Pair%AB2
!----------------------------------
       KA=Prim%KA
       KB=Prim%KB
       dJ=Zero
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
!----------------------------------
          Prim%CFA=CFA             
          Prim%CFB=CFB
          Prim%Ell=MaxLA+MaxLB+1 
!---------------------------------- 
          DO PFA=1,BS%NPFnc%I(CFA,KA)          
          DO PFB=1,BS%NPFnc%I(CFB,KB)
!----------------------------------
             Prim%ZA=BS%Expnt%D(PFA,CFA,KA)
             Prim%ZB=BS%Expnt%D(PFB,CFB,KB)
             Prim%Zeta=Prim%ZA+Prim%ZB
             Prim%Xi=Prim%ZA*Prim%ZB/Prim%Zeta
             IF(TestPrimPair(Prim%Xi,Prim%AB2))THEN
                Prim%PFA=PFA 
                Prim%PFB=PFB
                MaxAmp=SetBraBlok(Prim,BS,Gradients_O=.FALSE.)
!               Compute MAC and PAC
                DP2       = Zero
                PrimWCoef = Zero
                IA = IndexA
                DO LMNA=StartLA,StopLA
                   IA=IA+1
                   IB=IndexB
                   EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)                         
                   DO LMNB=StartLB,StopLB
                      IB=IB+1
                      EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)  
                      Pab     = P(IA,IB)
                      EllAB   = EllA+EllB+1
                      LenHGTF = LHGTF(EllAB)
                      LenSP   = LSP(EllAB)
                      PiZ     = (Pi/Prim%Zeta)**(ThreeHalves)
                      DO K=1,3
!                        Strength for MAC
                         CALL HGToSP_Direct(EllAB,LenHGTF,LenSP,PiZ,Pab*dHGBra%D(1:LenHGTF,IA,IB,K), &
                                            SPBraC(0:LenSP),SPBraS(0:LenSP))
                         PStrength = Zero
                         DO L=0,EllAB
                            PStrength = PStrength+FudgeFactorial(L,SPEll+1)*Unsold0(L,SPBraC,SPBraS)
                         ENDDO
                         PStrength=(PStrength/TauMAC)**(Two/DBLE(SPEll+2))
                         IF(DP2<PStrength) THEN
                            DP2=PStrength
                         ENDIF
!                        Strength for PAC
                         PrimWCoef = MAX(PrimWCoef,NodeWeight(EllAB,Prim%Zeta,Pab*dHGBra%D(1:LenHGTF,IA,IB,K)))
                      ENDDO
                   ENDDO
                ENDDO
                DP2=MIN(1.D10,DP2)
                IF(PrimWCoef>Zero .AND. DP2>Zero)THEN
!                  Initialize <KET|
                   CALL SetKet(Prim,PExtent)
!                  Initialize
                   LenSP=LSP(Prim%Ell)
                   LenHGTF=LHGTF(Prim%Ell)
                   HGKet_L(1:LenHGTF,1:3) = Zero
                   SPKetC_L(0:LenSP,1:3) = Zero
                   SPKetS_L(0:LenSP,1:3) = Zero
                   PTmp=Prim%P
!                  Sum over cells
                   DO NC=1,CS_IN%NCells 
                      Prim%P=PTmp+CS_IN%CellCarts%D(:,NC) 
                      PBox%Center=Prim%P
!                     Calculate the fractional coordinates
                      nlm =  AtomToFrac(GMLoc,CS_IN%CellCarts%D(:,NC))
!                     Store the Old Stuff
                      HGKet_old(1:LenHGTF)  = HGKet(1:LenHGTF)
                      SPKetC_old(0:LenSP) = SPKetC(0:LenSP)
                      SPKetS_old(0:LenSP) = SPKetS(0:LenSP)
!                     Walk the walk
                      CALL JWalk(PoleRoot)
!                     Acumulate the Lattice Forces
                      DO I=1,3
                         HGKet_L(1:LenHGTF,I)= HGKet_L(1:LenHGTF,I)+ nlm(I)*(HGKet(1:LenHGTF)- HGKet_old(1:LenHGTF))
                         SPKetC_L(0:LenSP,I) = SPKetC_L(0:LenSP,I) + nlm(I)*(SPKetC(0:LenSP) - SPKetC_old(0:LenSP))
                         SPKetS_L(0:LenSP,I) = SPKetS_L(0:LenSP,I) + nlm(I)*(SPKetS(0:LenSP) - SPKetS_old(0:LenSP))
                      ENDDO
                   ENDDO
                   Prim%P=PTmp
!                  Contract <Bra|Ket> bloks to compute matrix elements of J : SameAtom=.FALSE.
                   IA = IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1
                      IB=IndexB
                      EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)                         
                      DO LMNB=StartLB,StopLB
                         IB=IB+1
                         EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)                         
                         EllAB  = EllA+EllB+1
                         LenHGTF= LHGTF(EllAB)
                         LenSP  = LSP(EllAB)
                         PiZ    = (Pi/Prim%Zeta)**(ThreeHalves)
                         DO K=1,3
                            KI = K
                            DO LMN=1,LenHGTF
                               dJ(IA,IB,KI) = dJ(IA,IB,KI)+Phase%D(LMN)*dHGBra%D(LMN,IA,IB,K)*HGKet(LMN)
                            ENDDO
                            CALL HGToSP_Direct(EllAB,LenHGTF,LenSP,PiZ,dHGBra%D(1:LenHGTF,IA,IB,K), &
                                               SPBraC(0:LenSP),SPBraS(0:LenSP))
                            DO LM=0,LenSP
                               dJ(IA,IB,KI) = dJ(IA,IB,KI)+SPBraC(LM)*SPKetC(LM)+SPBraS(LM)*SPKetS(LM)
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO                   
!                  Calculate the FarField Multipole Contribution to the Matrix Element 
!                  Contract the Primative MM  with the density MM
                   IF(GMLoc%PBC%Dimen > 0) THEN
                      IA=IndexA
                      DO LMNA=StartLA,StopLA
                         IA=IA+1                    
                         IB=IndexB
                         DO LMNB=StartLB,StopLB  
                            IB=IB+1                      
                            DO K=1,3
                               KI = K
                               dJ(IA,IB,KI) = dJ(IA,IB,KI)+CTraxFF(Prim,dHGBra%D(:,IA,IB,K),GMLoc) 
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDIF
!                  Contract <Bra|Ket> bloks to compute matrix elements of J : SameAtom=.TRUE.
                   MaxAmp=SetBraBlok(Prim,BS,Gradients_O=.TRUE.)
                   IA = IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1
                      IB=IndexB
                      EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)                         
                      DO LMNB=StartLB,StopLB
                         IB=IB+1
                         EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)   
                         EllAB  = EllA+EllB+1
                         LenHGTF= LHGTF(EllAB)
                         LenSP  = LSP(EllAB)
                         PiZ    = (Pi/Prim%Zeta)**(ThreeHalves)
                         DO K=1,3
                            KI = K+3
                            DO LMN=1,LenHGTF
                               dJ(IA,IB,KI) = dJ(IA,IB,KI)+Phase%D(LMN)*dHGBra%D(LMN,IA,IB,K)*HGKet(LMN)
                            ENDDO
                            CALL HGToSP_Direct(EllAB,LenHGTF,LenSP,PiZ,dHGBra%D(1:LenHGTF,IA,IB,K), &
                                               SPBraC(0:LenSP),SPBraS(0:LenSP))
                            DO LM=0,LenSP
                               dJ(IA,IB,KI) = dJ(IA,IB,KI)+SPBraC(LM)*SPKetC(LM)+SPBraS(LM)*SPKetS(LM)
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO                   
!                  Calculate the FarField Multipole Contribution to the Matrix Element 
!                  Contract the Primative MM  with the density MM
                   IF(GMLoc%PBC%Dimen > 0) THEN
                      IA=IndexA
                      DO LMNA=StartLA,StopLA
                         IA=IA+1                    
                         IB=IndexB
                         DO LMNB=StartLB,StopLB  
                            IB=IB+1                      
                            DO K=1,3
                               KI = K+3
                               dJ(IA,IB,KI) = dJ(IA,IB,KI)+CTraxFF(Prim,dHGBra%D(:,IA,IB,K),GMLoc)  
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDIF
!                  Inner Sum J Box Forces
                   IA = IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1
                      IB=IndexB
                      EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)
                      DO LMNB=StartLB,StopLB
                         IB=IB+1
                         EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)
                         EllAB  = EllA+EllB+1
                         LenHGTF= LHGTF(EllAB)
                         LenSP  = LSP(EllAB)
                         PiZ    = (Pi/Prim%Zeta)**(ThreeHalves)
                         DO K=1,3
                            CALL HGToSP_Direct(EllAB,LenHGTF,LenSP,PiZ,dHGBra%D(1:LenHGTF,IA,IB,K), &
                                               SPBraC(0:LenSP),SPBraS(0:LenSP))
                            DO I=1,3
                               KI = K + 3*(I+1)
                               IF(GMLoc%PBC%AutoW%I(I)==1 .AND. GMLoc%PBC%AutoW%I(K)==1) THEN 
                                  DO LMN=1,LenHGTF
                                     dJ(IA,IB,KI)=dJ(IA,IB,KI)+Phase%D(LMN)*dHGBra%D(LMN,IA,IB,K)*HGKet_L(LMN,I)
                                  ENDDO
                                  DO LM=0,LenSP
                                     dJ(IA,IB,KI)=dJ(IA,IB,KI)+SPBraC(LM)*SPKetC_L(LM,I)+SPBraS(LM)*SPKetS_L(LM,I)
                                  ENDDO
                               ENDIF
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
!
                ENDIF
             ENDIF
          ENDDO 
          ENDDO
       ENDDO
       ENDDO
!
       DO K=1,15
          Vck(K)=BlkTrace_2(Pair%NA,Pair%NB,P,TRANSPOSE(dJ(:,:,K)))
       ENDDO
!
     END FUNCTION TrPdJ
!====================================================================================================
!
!====================================================================================================
    FUNCTION dNukE(GMLoc,At) RESULT(Vck)
       TYPE(CRDS)                        :: GMLoc
       REAL(DOUBLE)                      :: Tau,NukeCo,NukePole,PExtent
       REAL(DOUBLE),DIMENSION(4)         :: dBra
       REAL(DOUBLE),DIMENSION(0:SPLen)   :: SPBraC,SPBraS 
       INTEGER                           :: At,LenSP,LenHGTF,I,J,K,LM,LMN
       INTEGER                           :: NC,KI
       REAL(DOUBLE),DIMENSION(3)         :: PTmp,nlm
!
       REAL(DOUBLE),DIMENSION(15)        :: Vck
       REAL(DOUBLE),DIMENSION(0:SPLen,3) :: SPKetC_L,SPKetS_L
       REAL(DOUBLE),DIMENSION(1:HGLen,3) :: HGKet_L
       REAL(DOUBLE),DIMENSION(1:HGLen)   :: HGKet_old
       REAL(DOUBLE),DIMENSION(0:SPLen)   :: SPKetC_old,SPKetS_old
!---------------------------------------------------------------------------------------------
!      Initialize |dBRA>
       NukeCo=-GMLoc%AtNum%D(At)*(NuclearExpnt/Pi)**(ThreeHalves)
       DO K=1,3
          dHGBra%D(1:4,1,1,K)=Zero
       ENDDO
       dHGBra%D(2,1,1,1)=NukeCo
       dHGBra%D(3,1,1,2)=NukeCo
       dHGBra%D(4,1,1,3)=NukeCo
!      Initialize the primitive          
       Prim%Ell=1
       Prim%P=GMLoc%Carts%D(:,At)
       Prim%Zeta=NuclearExpnt
!      Set the MAC
       DP2=((FudgeFactorial(1,SPEll+1)*ABS(GMLoc%AtNum%D(At)))/TauMAC)**(Two/DBLE(SPEll+3))
       DP2=MIN(1.D10,DP2)
!      Set the PAC
       PrimWCoef = ABS(Three*GMLoc%AtNum%D(At)) 
!      Initialize <KET|
       CALL SetKet(Prim,PExtent)
       HGKet_L(:,:)  = Zero
       SPKetC_L(:,:) = Zero
       SPKetS_L(:,:) = Zero
!      KlumsyBl
       LenSP  =LSP(Prim%Ell)
       LenHGTF=LHGTF(Prim%Ell)
       PTmp=GMLoc%Carts%D(:,At)
       DO NC=1,CS_IN%NCells
!         Set atomic "primitive"
          Prim%P=PTmp+CS_IN%CellCarts%D(:,NC)
!         Calculate the fractional coordinates
          nlm = AtomToFrac(GMLoc,CS_IN%CellCarts%D(:,NC))
!         Store the Old Stuff
          HGKet_old(1:LenHGTF)= HGKet(1:LenHGTF)
          SPKetC_old(0:LenSP) = SPKetC(0:LenSP)
          SPKetS_old(0:LenSP) = SPKetS(0:LenSP)
!         Walk the walk
          CALL VWalk(PoleRoot)
!         Acumulate the Lattice Forces
          DO I=1,3
             HGKet_L(1:LenHGTF,I)= HGKet_L(1:LenHGTF,I)+ nlm(I)*(HGKet(1:LenHGTF)- HGKet_old(1:LenHGTF))
             SPKetC_L(0:LenSP,I) = SPKetC_L(0:LenSP,I) + nlm(I)*(SPKetC(0:LenSP) - SPKetC_old(0:LenSP))
             SPKetS_L(0:LenSP,I) = SPKetS_L(0:LenSP,I) + nlm(I)*(SPKetS(0:LenSP) - SPKetS_old(0:LenSP))
          ENDDO
       ENDDO
!      Reset the Atomic Coordinates
       Prim%P=PTmp
!      Init bra xforms
       Vck=Zero
       DO K=1,3
          DO LMN=1,LenHGTF
             Vck(K)=Vck(K)+Phase%D(LMN)*dHGBra%D(LMN,1,1,K)*HGKet(LMN)
          ENDDO
          CALL HGToSP(Prim%Zeta,1,dHGBra%D(:,1,1,K),SPBraC,SPBraS)
          DO LM=0,LenSP
             Vck(K)=Vck(K)+SPBraC(LM)*SPKetC(LM)+SPBraS(LM)*SPKetS(LM)
          ENDDO
       ENDDO
!      Add in the Far Field, Dipole and Quadripole Correction
       IF(GMLoc%PBC%Dimen>0) THEN
          DO K=1,3
             Vck(K)=Vck(K)+CTraxFF(Prim,dHGBra%D(:,1,1,K),GMLoc)
          ENDDO
       ENDIF
!      Inner Sum Nuc Box Force
       DO K=1,3
          CALL HGToSP(Prim%Zeta,1,dHGBra%D(:,1,1,K),SPBraC,SPBraS)
          DO I=1,3
             KI = K + 3*(I+1)
             IF(GMLoc%PBC%AutoW%I(I)==1 .AND. GMLoc%PBC%AutoW%I(K)==1) THEN
                DO LMN=1,LenHGTF
                   Vck(KI)=Vck(KI)+Phase%D(LMN)*dHGBra%D(LMN,1,1,K)*HGKet_L(LMN,I)
                ENDDO
                DO LM=0,LenSP
                   Vck(KI)=Vck(KI)+SPBraC(LM)*SPKetC_L(LM,I)+SPBraS(LM)*SPKetS_L(LM,I)
                ENDDO
             ENDIF
          ENDDO
       ENDDO  
!
     END FUNCTION dNukE
END MODULE BlokTrPdJ


