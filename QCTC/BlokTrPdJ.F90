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
#ifdef PERIODIC
  USE PBCFarField
#endif
  IMPLICIT NONE
  LOGICAL PrintFlag
!----------!
  CONTAINS !
!=======================================================================================
!
!=======================================================================================
     FUNCTION TrPdJ(Pair,P,GMLoc) RESULT(Vck)
       TYPE(AtomPair)                           :: Pair
       TYPE(CRDS)                               :: GMLoc
!
       REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)  :: P
       REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB,3):: dJ
       REAL(DOUBLE),DIMENSION(3)                :: PTmp,Vck
       REAL(DOUBLE),DIMENSION(0:SPLen)          :: SPBraC,SPBraS 
       REAL(DOUBLE)                             :: ZetaA,ZetaB,EtaAB,EtaIn,    &
                                                   XiAB,ExpAB,CA,CB,CC,Ov,     &
                                                   PAx,PAy,PAz,PBx,PBy,PBz,    &
                                                   MDx,MDxy,MDxyz,Amp2,MaxAmp, &
                                                   Pab,JNorm,Tau,OmegaMin,     &
                                                   PExtent,PStrength,Ext
       INTEGER                                  :: KA,KB,CFA,CFB,PFA,PFB,AtA,ATB,    &
                                                   IndexA,IndexB,              &
                                                   StartLA,StartLB,            &
                                                   StopLA,StopLB
       INTEGER                                  :: I,J,K,L,MaxLA,MaxLB,IA,IB,  &
                                                   LMNA,LMNB,LA,LB,MA,MB,    &
                                                   NA,NB,LAB,MAB,NAB,LM,LMN, &
                                                   Ell,EllA,EllB,EllAB,LenAB,&
                                                   HGLenEll,SPLenEll
       REAL(DOUBLE), EXTERNAL                   :: BlkTrace_2 
#ifdef PERIODIC
       INTEGER                                  :: NC
       REAL(DOUBLE)                             :: Px,Py,Pz
#endif
!----------------------------------------------------------------------------------------- 
!
       Prim%A=Pair%A
       Prim%B=Pair%B
       Prim%KA=Pair%KA
       Prim%KB=Pair%KB
       Prim%AB2=Pair%AB2
!----------------------------------
       KA=Prim%KA
       KB=Prim%KB
       dJ=Zero
!       PNorm=Zero
!       DO IA=1,Pair%NA; DO IB=1,Pair%NB
!          PNorm=PNorm+P(IA,IB)**2
!       ENDDO; ENDDO
!       PNorm=SQRT(PNorm)
!----------------------------------
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
                MaxAmp=SetBraBlok(Prim,BS,Gradients_O=Pair%SameAtom)
!---------------------------------------------------------------------------------------------
!               Compute maximal HG extent (for PAC) and Unsold esitmiate (for MAC)
!               Looping over all angular symmetries and Directions
                DP2=Zero
                PExtent=Zero
                IA = IndexA
                DO LMNA=StartLA,StopLA
                   IA=IA+1
                   IB=IndexB
                   EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)                         
                   DO LMNB=StartLB,StopLB
                      IB=IB+1
                      EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)    
                      Pab=P(IA,IB)
                      EllAB=EllA+EllB+1
                      LenAB=LHGTF(EllAB)
!                     Extent (for PAC)
                      DO K=1,3
                         Ext=Extent(EllAB,Prim%Zeta,Pab*dHGBra%D(1:LenAB,IA,IB,K),TauPAC,Potential_O=.TRUE.)
                         PExtent=MAX(PExtent,Ext)
!                        Strength (for MAC)
                         CALL HGToSP(Prim,Pab*dHGBra%D(1:LenAB,IA,IB,K),SPBraC,SPBraS)
                         DO L=0,EllA+EllB+1
                            PStrength = FudgeFactorial(L,SPEll+1)*Unsold0(L,SPBraC,SPBraS)
                            DP2       = MAX(DP2,(PStrength/TauMAC)**(Two/DBLE(SPELL+L+2)))
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
                DP2=MIN(1.D10,DP2)
!-----------------------------------------------------------------------------------------
                IF(PExtent>Zero)THEN ! Evaluate this primitives Ket contribution
!                  Initialize <KET|
                   CALL SetKet(Prim,PExtent)
#ifdef PERIODIC
#ifdef WRAPDIST
!                  WRAP the center of d Phi_A(r) Phi_B(r+R) back into the box
                   CALL AtomCyclic(GM,Prim%P)
#endif
                   PTmp=Prim%P
!                  Sum over cells
                   DO NC=1,CS_IN%NCells
                      Prim%P=PTmp+CS_IN%CellCarts%D(:,NC)
                      PBox%Center=Prim%P
!                     Walk the walk
                      CALL JWalk(PoleRoot)
                   ENDDO
                   Prim%P=PTmp
#else
!                  Walk the walk
                   CALL JWalk(PoleRoot)
#endif
!                  Contract <Bra|Ket> bloks to compute matrix elements of J
                   IA = IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1
                      IB=IndexB
                      EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)                         
                      DO LMNB=StartLB,StopLB
                         IB=IB+1
                         EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)                         
                         Ell=EllA+EllB+1
                         SPLenEll=LSP(Ell)
                         HGLenEll=LHGTF(Ell)
                         DO K=1,3
                            DO LMN=1,HGLenEll
                               dJ(IA,IB,K)=dJ(IA,IB,K)+Phase%D(LMN)*dHGBra%D(LMN,IA,IB,K)*HGKet(LMN)
                            ENDDO
                            CALL HGToSP(Prim,dHGBra%D(:,IA,IB,K),SPBraC,SPBraS)
                            DO LM=0,SPLenEll
                               dJ(IA,IB,K)=dJ(IA,IB,K)+SPBraC(LM)*SPKetC(LM)+SPBraS(LM)*SPKetS(LM)
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
#ifdef PERIODIC
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
                               dJ(IA,IB,K)=dJ(IA,IB,K) + CTraxFF_Grad(Prim,dHGBra%D(:,IA,IB,K),GMLoc) 
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDIF
#endif
             ENDIF
             ENDIF
          ENDDO 
          ENDDO
       ENDDO
       ENDDO
!
       DO K=1,3
          Vck(K)=BlkTrace_2(Pair%NA,Pair%NB,P,TRANSPOSE(dJ(:,:,K)))
       ENDDO
       IF(.NOT.Pair%SameAtom)THEN
          Vck=Four*Vck
       ELSE
          Vck=Two*Vck
       ENDIF
     END FUNCTION TrPdJ
!====================================================================================================
!
!====================================================================================================
    FUNCTION dNukE(GMLoc,At) RESULT(Vct)
       TYPE(CRDS)                      :: GMLoc
       REAL(DOUBLE)                    :: Tau,NukeCo,NukePole,PExtent
       REAL(DOUBLE),DIMENSION(4)       :: dBra
       REAL(DOUBLE),DIMENSION(3)       :: Vct
       REAL(DOUBLE),DIMENSION(0:SPLen) :: SPBraC,SPBraS 
       INTEGER                         :: At,SPLenEll,HGLenEll,K,LM,LMN
#ifdef PERIODIC
       INTEGER                         :: NC
       REAL(DOUBLE),DIMENSION(3)       :: PTmp
#endif
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
       PExtent=Extent(1,NuclearExpnt,dHGBra%D(:,1,1,1),TauPAC)
!      Initialize <KET|
       CALL SetKet(Prim,PExtent)
!      Klumsy
       SPLenEll=LSP(1)
       HGLenEll=LHGTF(1)
#ifdef PERIODIC
       PTmp=GMLoc%Carts%D(:,At)
       DO NC=1,CS_IN%NCells
!         Set atomic "primitive"
          Prim%P=PTmp+CS_IN%CellCarts%D(:,NC)
!         Walk the walk
          CALL VWalk(PoleRoot)
       ENDDO
!      Reset the Atomic Coordinates
       Prim%P=PTmp
!      Init bra xforms
       Vct=Zero
       DO K=1,3
          DO LMN=1,HGLenEll
             Vct(K)=Vct(K)+Phase%D(LMN)*dHGBra%D(LMN,1,1,K)*HGKet(LMN)
          ENDDO
          CALL HGToSP(Prim,dHGBra%D(:,1,1,K),SPBraC,SPBraS)
          DO LM=0,SPLenEll
             Vct(K)=Vct(K)+SPBraC(LM)*SPKetC(LM)+SPBraS(LM)*SPKetS(LM)
          ENDDO
       ENDDO
!      Add in the Far Field, Dipole and Quadripole Correction
       IF(GMLoc%PBC%Dimen>0) THEN
          DO K=1,3
             Vct(K)=Vct(K)+CTraxFF_Grad(Prim,dHGBra%D(:,1,1,K),GMLoc)
          ENDDO
       ENDIF
#else
!      Walk the walk
       CALL VWalk(PoleRoot)
!      <BRA|KET> 
       Vct=Zero
       DO K=1,3
          DO LMN=1,HGLenEll
             Vct(K)=Vct(K)+Phase%D(LMN)*dHGBra%D(LMN,1,1,K)*HGKet(LMN)
          ENDDO
          CALL HGToSP(Prim,dHGBra%D(:,1,1,K),SPBraC,SPBraS)
          DO LM=0,SPLenEll
             Vct(K)=Vct(K)+SPBraC(LM)*SPKetC(LM)+SPBraS(LM)*SPKetS(LM)
          ENDDO
       ENDDO
#endif
     END FUNCTION dNukE
END MODULE BlokTrPdJ
