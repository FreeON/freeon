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
MODULE JGen
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
    SUBROUTINE MakeJ(J)
      TYPE(BCSR)             :: J
      TYPE(AtomPair)         :: Pair
      INTEGER                :: AtA,AtB
      INTEGER                :: JP,K,NA,NB,NAB,P,Q,R,I,L
#ifdef PERIODIC        
      INTEGER                :: NCA,NCB
      REAL(DOUBLE)           :: Ax,Ay,Az,Bx,By,Bz
      LOGICAL                :: NotMakeBlock
#endif     
!---------------------------------------------- 
!     Initialize the matrix and associated indecies
      P=1; 
      R=1; 
      J%RowPt%I(1)=1
      CALL SetEq(J%MTrix,Zero)
!     Loop over atom pairs
      J%NAtms=NAtoms
      DO AtA=1,NAtoms
         DO AtB=1,NAtoms
            IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
               NAB = Pair%NA*Pair%NB
               IF(AtB<=AtA)THEN
!              Compute only the lower triangle of symmetric J
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
                           J%MTrix%D(R:R+NAB-1)=J%MTrix%D(R:R+NAB-1)+Two*JBlock(Pair,PoleRoot)
                           NotMakeBlock = .FALSE.
                        ENDIF
                     ENDDO
                  ENDDO
#else
                  J%MTrix%D(R:R+NAB-1)=Two*JBlock(Pair,PoleRoot)
#endif
               ENDIF
               J%ColPt%I(P)=AtB
               J%BlkPt%I(P)=R
               R=R+NAB
               P=P+1 
               J%RowPt%I(AtA+1)=P        
               IF(R>MaxNon0.OR.P>MaxBlks) &
                  CALL Halt(' BCSR dimensions blown in J ')
#ifdef PERIODIC
               IF(NotMakeBlock) &
                  CALL Halt(' Making a Zero Block in JGen ')
#endif
            ENDIF
         ENDDO
      ENDDO
      J%NBlks=P-1
      J%NNon0=R-1
!     Fill the upper triangle of J
      DO I=1,NAtoms
         DO JP=J%RowPt%I(I),J%RowPt%I(I+1)-1
            L=J%ColPt%I(JP)
            IF(I>L)THEN
               DO K=J%RowPt%I(L),J%RowPt%I(L+1)-1
                  IF(J%ColPt%I(K)==I)THEN
                     Q=J%BlkPt%I(K)                     
                     EXIT
                  ENDIF
               ENDDO               
               P=J%BlkPt%I(JP)
               NA=BS%BFKnd%I(GM%AtTyp%I(I))
               NB=BS%BFKnd%I(GM%AtTyp%I(L))
               NAB=NA*NB
               CALL XPose(NA,NB,J%MTrix%D(P:P+NAB-1),J%MTrix%D(Q:Q+NAB-1))
            ENDIF
         ENDDO
      ENDDO
!
    END SUBROUTINE MakeJ
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
     FUNCTION JBlock(Pair,PoleRoot) RESULT(Jvct)
       TYPE(AtomPair)                           :: Pair
       TYPE(PoleNode), POINTER                  :: PoleRoot
!
       REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)  :: JBlk
       REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB)  :: Jvct
       REAL(DOUBLE),DIMENSION(0:SPLen)          :: SPBraC,SPBraS 
       REAL(DOUBLE)                             :: ZetaA,ZetaB,EtaAB,EtaIn,    &
                                                   XiAB,ExpAB,CA,CB,CC,Ov,     &
                                                   PAx,PAy,PAz,PBx,PBy,PBz,    &
                                                   MDx,MDxy,MDxyz,Amp2,MaxAmp,Tau,OmegaMin
       INTEGER                                  :: KA,KB,CFA,CFB,PFA,PFB,      &
                                                   IndexA,IndexB,              &
                                                   StartLA,StartLB,            &
                                                   StopLA,StopLB
       INTEGER                                  :: I,J,MaxLA,MaxLB,IA,IB,    &
                                                   LMNA,LMNB,LA,LB,MA,MB,    &
                                                   NA,NB,LAB,MAB,NAB,LM,LMN,Ell,EllA,EllB
!-------------------------------------------------------------------------------------- 
       JBlk=Zero
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
                Prim%Zeta=EtaAB
                PAx=Prim%P(1)-Pair%A(1)
                PAy=Prim%P(2)-Pair%A(2)
                PAz=Prim%P(3)-Pair%A(3)
                PBx=Prim%P(1)-Pair%B(1)
                PBy=Prim%P(2)-Pair%B(2)
                PBz=Prim%P(3)-Pair%B(3)
!               McMurchie Davidson prefactors for HG Primitives
                CALL MD2TRR(BS%NASym,0,MaxLA,MaxLB,EtaAB,MD%D,PAx,PBx,PAy,PBy,PAz,PBz) 
!               Primitive coefficients in a HG representationj
                CALL SetBraBlok(CFA,PFA,KA,CFB,PFB,KB,ExpAB,BS,MD)!,Phase_O=.TRUE.)
!               Find the max amplitude for this primitive distribution
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
                      DO LMN=1,LHGTF(EllA+EllB)
                         Amp2=Amp2+HGBra%D(LMN,IA,IB)**2
                      ENDDO
                      MaxAmp=MAX(MaxAmp,SQRT(Amp2))
                   ENDDO
                ENDDO
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
                PoleSwitch=MAX(1.D-2,-LOG(Two*Tau/MaxAmp))
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
                      Ell=EllA+EllB
                      DO LMN=1,LHGTF(Ell)
                         JBlk(IA,IB)=JBlk(IA,IB)+Phase%D(LMN)*HGBra%D(LMN,IA,IB)*Prim%HGKet(LMN)

                      ENDDO
                      CALL HGToSP(Prim,HGBra%D(:,IA,IB),SPBraC,SPBraS)
                      DO LM=0,LSP(Ell)
                         JBlk(IA,IB)=JBlk(IA,IB)+SPBraC(LM)*Prim%SPKetC(LM) &
                                                +SPBraS(LM)*Prim%SPKetS(LM)
                      ENDDO
                  ENDDO
                  ENDDO
             ENDIF !End primitive thresholding
          ENDDO 
          ENDDO
       ENDDO
       ENDDO
       Jvct=BlockToVect(Pair%NA,Pair%NB,Jblk)
     END FUNCTION JBlock
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
       IF(Omega*PQ2>PoleSwitch)THEN
!         MAC: (d_q/|PQ|)^(p+1) < Tau/Amp
          IF(PQ2>Q%D2*DP2.OR.Q%Leaf)THEN 
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
END MODULE JGen
