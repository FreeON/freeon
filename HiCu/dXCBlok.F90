MODULE dXCBlok
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
                                                   Pab,PNorm
       INTEGER                                  :: KA,KB,CFA,CFB,PFA,PFB,      &
                                                   IndexA,IndexB,              &
                                                   StartLA,StartLB,            &
                                                   StopLA,StopLB,AtA,AtB,HGLenAB
       INTEGER                                  :: I,J,MaxLA,MaxLB,IA,IB,    &
                                                   LMNA,LMNB,LA,LB,MA,MB,    &
                                                   NA,NB,LAB,MAB,NAB,LMN,EllA,EllB,K
!-------------------------------------------------------------------------------------- 
       Prim%A=Pair%A
       Prim%B=Pair%B
       Prim%KA=Pair%KA
       Prim%KB=Pair%KB
       KA=Prim%KA
       KB=Prim%KB
       Prim%AB2=Pair%AB2
       Vck=Zero
       PNorm=Zero
       DO IA=1,Pair%NA; DO IB=1,Pair%NB
          PNorm=PNorm+P(IA,IB)**2
       ENDDO; ENDDO
       PNorm=SQRT(PNorm)
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
             IF(TestPrimPair(Prim%Xi,Prim%AB2))THEN
                Prim%PFA=PFA 
                Prim%PFB=PFB
!               Set primitive values, find distributions wheight
                MaxAmp=PNorm*SetBraBlok(Prim,BS,Gradients_O=Pair%SameAtom)
!---------------------------------------------------------------------------
                MaxAmp=LOG(MaxAmp)
!               Evaluate this primitives Ket contribution to Kxc_ab
                Ket=Zero
                PBox%BndBox(1,:)=Prim%P(1)
                PBox%BndBox(2,:)=Prim%P(2)
                PBox%BndBox(3,:)=Prim%P(3)
                PExtent=GaussianExtent(Prim%Zeta,MaxAmp)
                PBox=ExpandBox(PBox,PExtent)
!               Walk the walk
                CALL KxcWalk(CubeRoot)
!               Contract <Bra|Ket> bloks to compute matrix elements of Kxc
                IA = IndexA
                DO LMNA=StartLA,StopLA
                   IA=IA+1
                   IB=IndexB
                   EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)                         
                   DO LMNB=StartLB,StopLB
                      IB=IB+1
                      EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)                         
                      Pab=P(IA,IB)
                      HGLenAB=LHGTF(EllA+EllB+1)
                      DO K=1,3
                         DO LMN=1,HGLenAB
                            Vck(K)=Vck(K)+Pab*dHGBra%D(LMN,IA,IB,K)*Ket(LMN)
                         ENDDO
                      ENDDO
                   ENDDO
                ENDDO
!---------------------------------------------------------------------------
             ENDIF 
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
