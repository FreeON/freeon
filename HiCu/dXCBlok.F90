!------------------------------------------------------------------------------
!--  This code is part of the MondoSCF suite of programs for linear scaling 
!    electronic structure theory and ab initio molecular dynamics.
!
!--  Copyright (c) 2001, the Regents of the University of California.  
!    This SOFTWARE has been authored by an employee or employees of the 
!    University of California, operator of the Los Alamos National Laboratory 
!    under Contract No. W-7405-ENG-36 with the U.S. Department of Energy.  
!    The U.S. Government has rights to use, reproduce, and distribute this 
!    SOFTWARE.  The public may copy, distribute, prepare derivative works 
!    and publicly display this SOFTWARE without charge, provided that this 
!    Notice and any statement of authorship are reproduced on all copies.  
!    Neither the Government nor the University makes any warranty, express 
!    or implied, or assumes any liability or responsibility for the use of 
!    this SOFTWARE.  If SOFTWARE is modified to produce derivative works, 
!    such modified SOFTWARE should be clearly marked, so as not to confuse 
!    it with the version available from LANL.  The return of derivative works
!    to the primary author for integration and general release is encouraged. 
!    The first publication realized with the use of MondoSCF shall be
!    considered a joint work.  Publication of the results will appear
!    under the joint authorship of the researchers nominated by their
!    respective institutions. In future publications of work performed
!    with MondoSCF, the use of the software shall be properly acknowledged,
!    e.g. in the form "These calculations have been performed using MondoSCF, 
!    a suite of programs for linear scaling electronic structure theory and
!    ab initio molecular dynamics", and given appropriate citation.  
!------------------------------------------------------------------------------
!    Author: Matt Challacombe
!    COMPUTE THE EXCHANGE CORRELATION CONTRIBUTION TO THE TOTAL FORCE IN O(N)
!    USING HIERARCHICAL CUBATURE INVOLVING K-D BINARY TREES
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
                                                   Pab
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
!       PNorm=Zero
!       DO IA=1,Pair%NA; DO IB=1,Pair%NB
!          PNorm=PNorm+P(IA,IB)**2
!       ENDDO; ENDDO
!       PNorm=SQRT(PNorm)
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
!               Evaluate this primitives Ket contribution to Kxc_ab
                Ket=Zero
                PBox%BndBox(1,:)=Prim%P(1)
                PBox%BndBox(2,:)=Prim%P(2)
                PBox%BndBox(3,:)=Prim%P(3)
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
                      DO K=1,3
                         PExtent=MAX(PExtent, & 
                                     Extent(EllA+EllB,Prim%Zeta,Pab*dHGBra%D(:,IA,IB,K),TauRho,ExtraEll_O=2))
                      ENDDO
                   ENDDO
                ENDDO
                IF(PExtent>Zero)THEN
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
                          HGLenAB=LHGTF(EllA+EllB+1)
                          DO K=1,3
                             DO LMN=1,HGLenAB
                                Vck(K)=Vck(K)+Pab*dHGBra%D(LMN,IA,IB,K)*Ket(LMN)
                             ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
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
