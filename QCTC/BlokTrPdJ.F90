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
!    FAST O(Lg N), BLOKWISE ACCUMULATION OF Tr{P_(A,B).dJ^T_(A,B)/dA}
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
     FUNCTION TrPdJ(Pair,P) RESULT(Vck)
       TYPE(AtomPair)                           :: Pair
!
       REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)  :: P
       REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB,3):: dJ
       REAL(DOUBLE),DIMENSION(3)                :: Vck
       REAL(DOUBLE),DIMENSION(0:SPLen)          :: SPBraC,SPBraS 
       REAL(DOUBLE)                             :: ZetaA,ZetaB,EtaAB,EtaIn,    &
                                                   XiAB,ExpAB,CA,CB,CC,Ov,     &
                                                   PAx,PAy,PAz,PBx,PBy,PBz,    &
                                                   MDx,MDxy,MDxyz,Amp2,MaxAmp, &
                                                   PNorm,JNorm,Tau,OmegaMin
       INTEGER                                  :: KA,KB,CFA,CFB,PFA,PFB,AtA,ATB,    &
                                                   IndexA,IndexB,              &
                                                   StartLA,StartLB,            &
                                                   StopLA,StopLB
       INTEGER                                  :: I,J,K,MaxLA,MaxLB,IA,IB,  &
                                                   LMNA,LMNB,LA,LB,MA,MB,    &
                                                   NA,NB,LAB,MAB,NAB,LM,LMN, &
                                                   Ell,EllA,EllB,HGLenEll,SPLenEll
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
       PNorm=Zero
       DO IA=1,Pair%NA; DO IB=1,Pair%NB
          PNorm=PNorm+P(IA,IB)**2
       ENDDO; ENDDO
       PNorm=SQRT(PNorm)
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
                MaxAmp=PNorm*SetBraBlok(Prim,BS,Gradients_O=Pair%SameAtom)
!                MaxAmp=SetBraBlok(Prim,BS,Gradients_O=Pair%SameAtom)
!-----------------------------------------------------------------------------------------
!               Evaluate this primitives Ket contribution to J_ab
                HGKet=Zero
                SPKetC=Zero
                SPKetS=Zero
!               Set MAC and PAC parameters
                Tau=Thresholds%TwoE
                DP2=(MaxAmp/Tau)**(Two/DBLE(SPEll+1))
                PoleSwitch=PFunk(Prim%Ell+MaxL,Tau/MaxAmp)
!               Zero the Acumulators
                HGKet=Zero
                SPKetC=Zero
                SPKetS=Zero 
#ifdef PERIODIC
!               Fold The Primative Back into the BOX
                CALL AtomCyclic(GM,Prim%P)
                Px = Prim%P(1)
                Py = Prim%P(2)
                Pz = Prim%P(3)
!               Sum over cells
                DO NC=1,CSMM1%NCells
                   Prim%P(1)=Px+CSMM1%CellCarts%D(1,NC)
                   Prim%P(2)=Py+CSMM1%CellCarts%D(2,NC)
                   Prim%P(3)=Pz+CSMM1%CellCarts%D(3,NC)
!                  Walk the walk
                   CALL JWalk(PoleRoot)
                ENDDO
                Prim%P(1)=Px
                Prim%P(2)=Py
                Prim%P(3)=Pz
#else
!               Walk the walk
                CALL JWalk(PoleRoot)
#endif
!               Contract <Bra|Ket> bloks to compute matrix elements of J
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
!               Calculate the FarField Multipole Contribution to the Matrix Element 
!               Contract the Primative MM  with the density MM
                IF(Dimen > 0) THEN
                   IA=IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1                    
                      IB=IndexB
                      DO LMNB=StartLB,StopLB  
                         IB=IB+1                      
                         DO K=1,3
                            dJ(IA,IB,K)=dJ(IA,IB,K) + CTraxFF(Prim,dHGBra%D(:,IA,IB,K)) &
                                                    + CTraxQ(Prim,dHGBra%D(:,IA,IB,K))                                 
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
#endif
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
    FUNCTION dNukE(At) RESULT(Vct)
       REAL(DOUBLE)                    :: Tau,NukeCo,NukePole
       REAL(DOUBLE),DIMENSION(4)       :: dBra
       REAL(DOUBLE),DIMENSION(3)       :: Vct
       REAL(DOUBLE),DIMENSION(0:SPLen) :: SPBraC,SPBraS 
       INTEGER                         :: At,SPLenEll,HGLenEll,K,LM,LMN
#ifdef PERIODIC
       INTEGER                         :: NC
       REAL(DOUBLE),DIMENSION(3)       :: QC
#endif
!---------------------------------------------------------------------------------------------
       NukeCo   =-GM%AtNum%I(At)*(NuclearExpnt/Pi)**(ThreeHalves)
       NukePole =-GM%AtNum%I(At)
!      Set MAC and PAC parameters
       Tau      =Thresholds%TwoE
       DP2      =(ABS(NukePole)/Tau)**(Two/DBLE(SPEll+1))
       PoleSwitch=Gamma_Switch
!      Set atomic "primitive"
       Prim%Ell=1
       Prim%Zeta=NuclearExpnt
       SPLenEll=LSP(Prim%Ell)
       HGLenEll=LHGTF(Prim%Ell)
!      Zero accumulators
       HGKet=Zero
       SPKetC=Zero
       SPKetS=Zero
#ifdef PERIODIC
       QC(:)  = GM%Carts%D(:,At)
       DO NC=1,CSMM1%NCells
!         Set atomic "primitive"
          Prim%P(:)=QC(:)+CSMM1%CellCarts%D(:,NC)
!         Walk the walk
          CALL VWalk(PoleRoot)
       ENDDO
!      Reset the Atomic Coordinates
       Prim%P = QC
!      Init bra xforms
       Vct=Zero
       dHGBra%D=Zero       
       dHGBra%D(2,1,1,1)=NukeCo
       dHGBra%D(3,1,1,2)=NukeCo
       dHGBra%D(4,1,1,3)=NukeCo
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
       IF(Dimen > 0) THEN
          DO K=1,3
             Vct(K)=Vct(K) + CTraxFF(Prim,dHGBra%D(:,1,1,K)) &
                           + CTraxQ(Prim,dHGBra%D(:,1,1,K))
          ENDDO
       ENDIF
#else
       Prim%P=GM%Carts%D(:,At) 
!      Walk the walk
       CALL VWalk(PoleRoot)
!      Init bra xforms
       Vct=Zero
       dHGBra%D=Zero       
       dHGBra%D(2,1,1,1)=NukeCo
       dHGBra%D(3,1,1,2)=NukeCo
       dHGBra%D(4,1,1,3)=NukeCo
       DO K=1,3
          DO LMN=1,HGLenEll
             Vct(K)=Vct(K) + Phase%D(LMN)*dHGBra%D(LMN,1,1,K)*HGKet(LMN)
          ENDDO
          CALL HGToSP(Prim,dHGBra%D(:,1,1,K),SPBraC,SPBraS)
          DO LM=0,SPLenEll
             Vct(K)=Vct(K)+SPBraC(LM)*SPKetC(LM)+SPBraS(LM)*SPKetS(LM)
          ENDDO
       ENDDO
#endif
     END FUNCTION dNukE
END MODULE BlokTrPdJ
