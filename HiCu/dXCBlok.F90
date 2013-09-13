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
     FUNCTION dXC(Pair,P,NSMat) RESULT(Vck)
       TYPE(AtomPair)                           :: Pair
       INTEGER::NSMat
       REAL(DOUBLE),DIMENSION(6)                :: Vck
       REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB*NSMat)  :: P
       REAL(DOUBLE)                             :: ZetaA,ZetaB,EtaAB,EtaIn,    &
                                                   XiAB,ExpAB,CA,CB,CC,Ov,     &
                                                   PAx,PAy,PAz,PBx,PBy,PBz,    &
                                                   MDx,MDxy,MDxyz,Amp2,MaxAmp, &
                                                   Pab,NewEx,Pab1,Pab2
       INTEGER                                  :: KA,KB,CFA,CFB,PFA,PFB,      &
                                                   IndexA,IndexB,              &
                                                   StartLA,StartLB,            &
                                                   StopLA,StopLB,AtA,AtB
       INTEGER                                  :: I,J,K,MaxLA,MaxLB,IA,IB,  &
                                                   LMNA,LMNB,LA,LB,MA,MB,    &
                                                   NA,NB,LAB,MAB,NAB,LMN,    &
                                                   EllA,EllB,EllAB,LenAB,KI
       LOGICAL                                  :: SameAtom
!--------------------------------------------------------------------------------------
       SameAtom=.FALSE.
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
             IF(TestPrimPair(Prim%Xi,Prim%AB2))THEN
                Prim%PFA=PFA
                Prim%PFB=PFB
!               Set primitive values, find distributions wheight
                PExtent=SetBraBlok(Prim,BS,Gradients_O=.TRUE.,Tau_O=TauRho,ExtraEll_O=1)
                PExtent=MAX(PExtent,SetBraBlok(Prim,BS,Gradients_O=.FALSE.,Tau_O=TauRho,ExtraEll_O=1))
                PBox%BndBox(:,1)=Prim%P
                PBox%BndBox(:,2)=Prim%P
                PBox=ExpandBox(PBox,PExtent)
!               Quick check to see if primitive touches the grid
                IF(PExtent>Zero.AND.(.NOT.BoxOutSideBox(PBox,CubeRoot%Box)))THEN
!-----------------------------------------------------------------------------------------
!                  Find the extent of this primitive again, now factoring in the density matrix
                   PExtent=Zero
                   IA = IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1
                      IB=IndexB
                      EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)
                      DO LMNB=StartLB,StopLB
                         IB=IB+1
                         EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)
                         Pab=P(IA,IB) ! <<< SPIN !We must run over all the NSMat here!
                         EllAB=EllA+EllB+1
                         LenAB=LHGTF(EllAB)
                         DO K=1,3
                            PExtent=MAX(PExtent,Extent(EllAB,Prim%Zeta,  &
                                        Pab*dHGBra%D(1:LenAB,IA,IB,K),   &
                                        TauRho,ExtraEll_O=1))
                         ENDDO
                      ENDDO
                   ENDDO
                   Ket=Zero
!                  Set BBox for this primitive
                   PBox%BndBox(:,1)=Prim%P
                   PBox%BndBox(:,2)=Prim%P
                   PBox=ExpandBox(PBox,PExtent)
!                  Walk the walk
                   CALL KxcWalk(CubeRoot)
!                  Contract <Bra|Ket> bloks to compute matrix elements of Kxc:  SameAtom=.FALSE.
                   IA = IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1
                      IB=IndexB
                      EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)
                      DO LMNB=StartLB,StopLB
                         IB=IB+1
                          EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)
                          LenAB=LHGTF(EllA+EllB+1)
                          !Q&D... I know!
                          IF(NSMat.EQ.1)THEN
                             Pab=P(IA,IB)
                             DO K=1,3
                                KI = K
                                DO LMN=1,LenAB
                                   Vck(KI)=Vck(KI)+Pab*dHGBra%D(LMN,IA,IB,K)*Ket(LMN)
                                ENDDO
                             ENDDO
                          ELSEIF(NSMat.EQ.2)THEN
                             Pab1=P(IA,IB)
                             Pab2=P(IA,IB+Pair%NB)
                             DO K=1,3
                                KI = K
                                DO LMN=1,LenAB
                                   Vck(KI)=Vck(KI)+dHGBra%D(LMN,IA,IB,K)*( Pab1*Ket(LMN) + Pab2*Ket(LMN+LenAB) )
                                ENDDO
                             ENDDO
                          ELSE
                             CALL Halt('dXC: you are to impatient!')
                          ENDIF
                      ENDDO
                   ENDDO
!                  Contract <Bra|Ket> bloks to compute matrix elements of Kxc:  SameAtom=.TRUE.
                   PExtent=SetBraBlok(Prim,BS,Gradients_O=.TRUE.,Tau_O=TauRho,ExtraEll_O=1)
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
                          !Q&D... I know!
                          IF(NSMat.EQ.1)THEN
                             DO K=1,3
                                KI = K+3
                                DO LMN=1,LenAB
                                   Vck(KI)=Vck(KI)+Pab*dHGBra%D(LMN,IA,IB,K)*Ket(LMN)
                                ENDDO
                             ENDDO
                          ELSEIF(NSMat.EQ.2)THEN
                             Pab1=P(IA,IB)
                             Pab2=P(IA,IB+Pair%NB)
                             DO K=1,3
                                KI = K+3
                                DO LMN=1,LenAB
                                   Vck(KI)=Vck(KI)+dHGBra%D(LMN,IA,IB,K)*( Pab1*Ket(LMN) + Pab2*Ket(LMN+LenAB) )
                                ENDDO
                             ENDDO
                          ELSE
                             CALL Halt('dXC: you are to impatient!')
                          ENDIF

                      ENDDO
                   ENDDO
!
                ENDIF
!---------------------------------------------------------------------------
             ENDIF
          ENDDO
          ENDDO
       ENDDO
       ENDDO
!
     END FUNCTION dXC

!====================================================================================================
!
!====================================================================================================
END MODULE dXCBlok
