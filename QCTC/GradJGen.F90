!==============================================================================
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe and C.J. Tymczak
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!--  COMPUTE THE COULOMB MATRIX IN O(N Lg N) CPU TIME
!==============================================================================
MODULE GradJGen
  USE DerivedTypes
  USE GlobalScalars   
  USE GlobalObjects
  USE ProcessControl
  USE Indexing
  USE InOut
  USE Thresholding
  USE Globals
  USE AtomPairs
  USE BraKetBloks
  USE PoleTree
  USE MondoPoles 
  IMPLICIT NONE
  LOGICAL PrintFlag
!---------------------------------------
  TYPE(PrimPair)  :: Prim
  REAL(DOUBLE)    :: DP2
  REAL(DOUBLE)    :: PoleSwitch
!----------!
  CONTAINS !
!=============================================================================================
!
!=============================================================================================
     FUNCTION GradJBlock(BS,MD,AtA,AtB,DiffAtm,P,Pair,PoleRoot) RESULT(JPvck)
       TYPE(AtomPair)                           :: Pair
       TYPE(BSET),                 INTENT(IN)    :: BS
       TYPE(DBL_RNK4),             INTENT(INOUT) :: MD
       TYPE(PoleNode), POINTER                  :: PoleRoot
!
       REAL(DOUBLE), EXTERNAL                  :: BlkTrace_2
       REAL(DOUBLE),DIMENSION(3)                :: JPVck
       REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)  :: JBlk,JXP,JYP,JZP,P
       REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB)  :: Jvct
       REAL(DOUBLE),DIMENSION(0:SPLen+1)        :: SPBraC,SPBraS 
       REAL(DOUBLE)                             :: ZetaA,ZetaB,EtaAB,EtaIn,    &
                                                   XiAB,ExpAB,CA,CB,CC,Ov,     &
                                                   PAx,PAy,PAz,PBx,PBy,PBz,    &
                                                   MDx,MDxy,MDxyz,Amp2,MaxAmp, &
                                                   Tau,OmegaMin,Sum
       INTEGER                                  :: KA,KB,CFA,CFB,PFA,PFB,      &
                                                   IndexA,IndexB,              &
                                                   StartLA,StartLB,            &
                                                   StopLA,StopLB
       INTEGER                                  :: I,J,MaxLA,MaxLB,IA,IB,    &
                                                   LMNA,LMNB,LA,LB,MA,MB,    &
                                                   NA,NB,LAB,MAB,NAB,LM,LMN, &
                                                   Ell,EllA,EllB,AtA,AtB,DiffAtm
!-------------------------------------------------------------------------------------- 
       JBlk=Zero
       JXP=Zero
       JYP=Zero
       JZP=Zero
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
                Prim%Ell=MaxLA+MaxLB+1
                Prim%Zeta=EtaAB
                PAx=Prim%P(1)-Pair%A(1)
                PAy=Prim%P(2)-Pair%A(2)
                PAz=Prim%P(3)-Pair%A(3)
                PBx=Prim%P(1)-Pair%B(1)
                PBy=Prim%P(2)-Pair%B(2)
                PBz=Prim%P(3)-Pair%B(3)
!               McMurchie Davidson prefactors for HG Primitives
                CALL MD2TRR(BS%NASym+2,-1,MaxLA+1,MaxLB+1,EtaAB,MD%D,PAx,PBx,PAy,PBy,PAz,PBz) 
!               Derivative of Primitive coefficients in a HG representation
                CALL dTwoE(DiffAtm,AtA,AtB,CFA,PFA,KA,CFB,PFB,KB,ExpAB,BS,MD,Pair)
!               Find the max amplitude for trace involving this primitive distribution
                MaxAmp=Zero
                Sum=Zero
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
                      Sum=Sum+ABS(P(IA,IB))*SQRT(Amp2)
                      MaxAmp=MAX(MaxAmp,SQRT(Amp2))
                   ENDDO
                ENDDO
!               MaxAmp=Sum ??????????????
!               Evaluate this primitives Ket contribution to J_ab
                Prim%HGKet=Zero
                Prim%SPKetC=Zero!FarFC
                Prim%SPKetS=Zero!FarFS
                Prim%Box%BndBox(1,:)=Prim%P(1)
                Prim%Box%BndBox(2,:)=Prim%P(2)
                Prim%Box%BndBox(3,:)=Prim%P(3)
                Prim%Box%Center=Prim%P
!               Set MAC and PAC parameters
                Tau=Thresholds%TwoE
                DP2=(MaxAmp/Tau)**(Two/DBLE(SPEll+1))
                PoleSwitch=PFunk(Prim%Ell+PoleRoot%Ell,Tau/MaxAmp)
!
!               WRITE(*,22)Prim%Ell+PoleRoot%Ell,Tau/MaxAmp,PoleSwitch
!           22 FORMAT(' Ell = ',I2,' Acc = ',D12.3,' PoleSwitch = ',F12.6)
!
!               Walk the walk
                CALL JWalk(PoleRoot)
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
                      DO LMN=1,LHGTF(Ell)
!                         JBlk(IA,IB)=JBlk(IA,IB)+Phase%D(LMN)*HGBra%D(LMN,IA,IB)*Prim%HGKet(LMN)
                         JXP(IA,IB) = JXP(IA,IB) + Two*Phase%D(LMN)*dHGBradX%D(LMN,IA,IB)*Prim%HGKet(LMN)
                         JYP(IA,IB) = JYP(IA,IB) + Two*Phase%D(LMN)*dHGBradY%D(LMN,IA,IB)*Prim%HGKet(LMN)
                         JZP(IA,IB) = JZP(IA,IB) + Two*Phase%D(LMN)*dHGBradZ%D(LMN,IA,IB)*Prim%HGKet(LMN)
                      ENDDO
!                      CALL HGToSP(Prim,HGBra%D(:,IA,IB),SPBraC,SPBraS)
!                      DO LM=0,LSP(Ell)
!                         JBlk(IA,IB)=JBlk(IA,IB)+SPBraC(LM)*Prim%SPKetC(LM)+SPBraS(LM)*Prim%SPKetS(LM)
!                      ENDDO

                      CALL HGToSP(Prim,dHGBradX%D(:,IA,IB),SPBraC,SPBraS)
                      DO LM=0,LSP(Ell)
                         JXP(IA,IB)=JXP(IA,IB)+Two*(SPBraC(LM)*Prim%SPKetC(LM)+SPBraS(LM)*Prim%SPKetS(LM))
                      ENDDO

                      CALL HGToSP(Prim,dHGBradY%D(:,IA,IB),SPBraC,SPBraS)
                      DO LM=0,LSP(Ell)
                         JYP(IA,IB)=JYP(IA,IB)+Two*(SPBraC(LM)*Prim%SPKetC(LM)+SPBraS(LM)*Prim%SPKetS(LM))
                      ENDDO

                      CALL HGToSP(Prim,dHGBradZ%D(:,IA,IB),SPBraC,SPBraS)
                      DO LM=0,LSP(Ell)
                         JZP(IA,IB)=JZP(IA,IB)+Two*(SPBraC(LM)*Prim%SPKetC(LM)+SPBraS(LM)*Prim%SPKetS(LM))
                      ENDDO

                  ENDDO
                  ENDDO
             ENDIF !End primitive thresholding
          ENDDO 
          ENDDO
       ENDDO
       ENDDO
!       WRITE(*,*) AtA,AtB,DiffAtm,JZP
!       Jvct=BlockToVect(Pair%NA,Pair%NB,Jblk)
       JPVck(1)=  BlkTrace_2(Pair%NB,Pair%NA,P,JXP(1,1))
       JPVck(2)=  BlkTrace_2(Pair%NB,Pair%NA,P,JYP(1,1))
       JPVck(3)=  BlkTrace_2(Pair%NB,Pair%NA,P,JZP(1,1))

!      WRITE(*,*) AtA,AtB,DiffAtm,Sqrt(JPVck(1)**2+JPVck(2)*2+JPVck(3)**2)
     END FUNCTION GradJBlock
!====================================================================================================
!
!====================================================================================================
     FUNCTION PFunk(Ell,Ack)
        INTEGER                    :: Ell
        REAL(DOUBLE)               :: PFunk,X,Ack,MinAcc,MaxAcc
        REAL(DOUBLE),DIMENSION(30) :: W
        X=-LOG(Ack)
        INCLUDE "PFunk.Inc"
     END FUNCTION PFunk    
!====================================================================================================
!
!====================================================================================================
     RECURSIVE SUBROUTINE JWalk(Q)
       TYPE(PoleNode)                   :: Q
       REAL(DOUBLE)                     :: PQ2,PQ,PQxy,PQx,PQy,PQx2,PQy2,PQz,Omega, &
                                           CoTan,OneOvPQ,OneOvPQxy,RS,SQ,PQToThMnsL, &
                                           TT,TwoC,COne,SOne,CTwo,STwo  
       REAL(DOUBLE),DIMENSION(0:2*SPEll):: CoFact
       REAL(DOUBLE),DIMENSION(0:SPLen2) :: ALegendreP
       REAL(DOUBLE),DIMENSION(50)       :: W
       INTEGER                          :: Ell,LCode
!---------------------------------------------------------------------------------------------------
!      PAC: Exp[-W_pq PQ^2]/2 < Tau/Amp
       PQx=Prim%P(1)-Q%Box%Center(1)
       PQy=Prim%P(2)-Q%Box%Center(2)
       PQz=Prim%P(3)-Q%Box%Center(3)
       PQ2=PQx*PQx+PQy*PQy+PQz*PQz
       Omega=Prim%Zeta*Q%Zeta/(Prim%Zeta+Q%Zeta)
!       IF(Omega*PQ2>10.0)THEN
       IF(Omega*PQ2>PoleSwitch)THEN
!         MAC: (d_q/|PQ|)^(p+1) < Tau/Amp
          IF(PQ2>Q%D2*DP2.OR.Q%Leaf)THEN 
!          IF(Q%Leaf)THEN 
!            Evaluate multipoles
             Ell=Prim%Ell+Q%Ell
             LCode=100*Prim%Ell+Q%Ell 
#ifdef EXPLICIT_SOURCE
             INCLUDE "IrRegulars.Inc"
             INCLUDE "CTraX.Inc"
#else
             CALL CTrax(Prim,Q)
#endif
          ELSE
!            Keep walking
             CALL JWalk(Q%Descend)
             CALL JWalk(Q%Descend%Travrse)
          ENDIF
       ELSEIF(Q%Leaf)THEN
!         Compute an ERI
          CALL MDERI(Prim,Q)
       ELSE
!         Keep walking
          CALL JWalk(Q%Descend)
          CALL JWalk(Q%Descend%Travrse)
       ENDIF
     END SUBROUTINE JWalk
END MODULE GradJGen
