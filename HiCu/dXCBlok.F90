!    COMPUTE THE EXCHANGE CORRELATION CONTRIBUTION TO THE TOTAL FORCE IN O(N)
!    USING HIERARCHICAL CUBATURE AND K-D BINARY TREES
!    Author: Matt Challacombe
!------------------------------------------------------------------------------
MODULE dXCBlok
  USE DerivedTypes
  USE GlobalScalars   
  USE GlobalObjects
  USE ProcessControl
  USE Indexing
  USE McMurchie
  USE InOut
  USE Thresholding
  USE HiCuThresholds
  USE AtomPairs
  USE BraBloks
  USE CubeTree 
  USE TreeWalk
  IMPLICIT NONE
!----------!
  CONTAINS !
!=======================================================================================
!
!=======================================================================================
     FUNCTION dXC(Pair,P) RESULT(Vck)
       TYPE(AtomPair)                           :: Pair
       REAL(DOUBLE),DIMENSION(3)                :: Vck
       REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)  :: P
       REAL(DOUBLE)                             :: ZetaA,ZetaB,EtaAB,EtaIn,    &
                                                   XiAB,ExpAB,CA,CB,CC,Ov,     &
                                                   PAx,PAy,PAz,PBx,PBy,PBz,    &
                                                   MDx,MDxy,MDxyz,Amp2,MaxAmp, &
                                                   Pab,NewEx
       INTEGER                                  :: KA,KB,CFA,CFB,PFA,PFB,      &
                                                   IndexA,IndexB,              &
                                                   StartLA,StartLB,            &
                                                   StopLA,StopLB,AtA,AtB
       INTEGER                                  :: I,J,K,MaxLA,MaxLB,IA,IB,  &
                                                   LMNA,LMNB,LA,LB,MA,MB,    &
                                                   NA,NB,LAB,MAB,NAB,LMN,    &
                                                   EllA,EllB,EllAB,LenAB,NC
!-------------------------------------------------------------------------------------- 
       Prim%A=Pair%A
       Prim%B=Pair%B
       Prim%KA=Pair%KA
       Prim%KB=Pair%KB
       KA=Prim%KA
       KB=Prim%KB
       Prim%AB2=Pair%AB2
       Vck=Zero
!----------------------------------
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
!
#ifdef PERIODIC
             CALL SetKxcCell(GM,Prim%Ell,Prim%ZA,Prim%ZB,ExtraEll_O=2) 
             DO NC=1,CS_Kxc%NCells
                Prim%A   = Pair%A + CS_Kxc%CellCarts%D(1:3,NC)
                Prim%B   = Pair%B + CS_Kxc%CellCarts%D(4:6,NC)
                Prim%AB2 = (Prim%A(1)-Prim%B(1))**2 &
                         + (Prim%A(2)-Prim%B(2))**2 &
                         + (Prim%A(3)-Prim%B(3))**2
#endif
             IF(TestPrimPair(Prim%Ell,Prim%Xi,Prim%AB2)) THEN
                Prim%PFA=PFA 
                Prim%PFB=PFB
!               Set primitive values, find distributions wheight
                MaxAmp=SetBraBlok(Prim,BS,Gradients_O=Pair%SameAtom)
!-----------------------------------------------------------------------------------------
!               Find the maximum extent of this primitive
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
                      DO K=1,3
                         PExtent=MAX(PExtent,Extent(EllAB,Prim%Zeta,  &
                                     Pab*dHGBra%D(1:LenAB,IA,IB,K),   &
                                     TauRho,ExtraEll_O=2))
                      ENDDO
                   ENDDO
                ENDDO
                IF(PExtent>Zero)THEN
                   Ket=Zero
!                  Set BBox for this primitive
                   PBox%BndBox(:,1)=Prim%P
                   PBox%BndBox(:,2)=Prim%P
                   PBox=ExpandBox(PBox,PExtent)
!                  Walk the walk
                   CALL KxcWalk(CubeRoot)
!                  Contract <Bra|Ket> bloks to compute matrix elements of Kxc
                   IA = IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1
                      IB=IndexB
                      EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)                         
                      DO LMNB=StartLB,StopLB
                         IB=IB+1
                          EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)                         
                          Pab=P(IA,IB)
                          LenAB=LHGTF(EllA+EllB+1)
                          DO K=1,3
                             DO LMN=1,LenAB
                                Vck(K)=Vck(K)+Pab*dHGBra%D(LMN,IA,IB,K)*Ket(LMN)
                             ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
!---------------------------------------------------------------------------
             ENDIF 
#ifdef PERIODIC
             ENDDO
#endif
          ENDDO 
          ENDDO
       ENDDO
       ENDDO
       IF(.NOT.Pair%SameAtom)THEN
          Vck=Four*Vck
       ELSE
          Vck=Two*Vck
       ENDIF
     END FUNCTION dXC

!====================================================================================================
!
!====================================================================================================
END MODULE dXCBlok
