!    Authors: Matt Challacombe and C.J. Tymczak
!    COMPUTE THE COULOMB MATRIX IN O(N Lg N) CPU TIME
!-------------------------------------------------------------------------------
MODULE JGen
  USE DerivedTypes
  USE GlobalScalars   
  USE GlobalObjects
  USE ProcessControl
  USE Indexing
  USE InOut
  USE QCTCThresholds
  USE Globals
  USE AtomPairs
  USE BraBloks
  USE PoleTree
  USE TreeWalk
#ifdef PERIODIC
  USE PBCFarField
#endif
#ifdef PARALLEL
  USE FastMatrices
#endif
  IMPLICIT NONE
  LOGICAL PrintFlag
!-------------------------------------------------------------------------------
  CONTAINS
!===============================================================================

!===============================================================================
    SUBROUTINE MakeJ(J)
#ifdef PARALLEL
      TYPE(FastMat),POINTER     :: J
#else
      TYPE(BCSR)                :: J
#endif
      TYPE(DBL_RNK2)            :: Temp
      TYPE(AtomPair)            :: Pair
      INTEGER                   :: AtA,AtB
      INTEGER                   :: JP,K,NA,NB,NAB,P,Q,R,I1,I2,I3,L,I
#ifdef PERIODIC 
      INTEGER                   :: NC
      REAL(DOUBLE),DIMENSION(3) :: B
#endif    
!------------------------------------------------------------------------------- 
!     Initialize the matrix and associated indecies
#ifdef PARALLEL
#else
      P=1
      R=1
      J%RowPt%I(1)=1
      CALL SetEq(J%MTrix,Zero)
      J%NAtms= NAtoms
#endif
      ! Loop over atom pairs
#ifdef PARALLEL
      DO AtA=Beg%I(MyId),End%I(MyId)
#else
      DO AtA=1,NAtoms            
#endif
         DO AtB=1,NAtoms
            IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
               NAB = Pair%NA*Pair%NB
#ifndef PARALLEL
               ! Compute only the lower triangle of symmetric J
               IF(AtB<=AtA)THEN  
#endif
#ifdef PERIODIC
                  B = Pair%B
                  DO NC=1,CS_OUT%NCells
                     Pair%B   = B+CS_OUT%CellCarts%D(:,NC)
                     Pair%AB2 = (Pair%A(1)-Pair%B(1))**2+(Pair%A(2)-Pair%B(2))**2+(Pair%A(3)-Pair%B(3))**2
                     IF(TestAtomPair(Pair)) THEN
#ifdef PARALLEL
                        CALL AddFASTMATBlok(J,AtA,AtB,Two*JBlock(Pair,PoleRoot))
#else
                        J%MTrix%D(R:R+NAB-1)=J%MTrix%D(R:R+NAB-1)+Two*JBlock(Pair,PoleRoot)
#endif
                     ENDIF
                  ENDDO
#else

#ifdef PARALLEL
                  CALL AddFASTMATBlok(J,AtA,AtB,Two*JBlock(Pair,PoleRoot))
#else
                  J%MTrix%D(R:R+NAB-1)=Two*JBlock(Pair,PoleRoot)
#endif

#endif

#ifndef PARALLEL
               ENDIF
#endif

#ifdef PARALLEL
#else
               J%ColPt%I(P)=AtB
               J%BlkPt%I(P)=R
               R=R+NAB
               P=P+1  
               J%RowPt%I(AtA+1)=P        
               IF(R>MaxNon0.OR.P>MaxBlks) &
                  CALL Halt(' BCSR dimensions blown in J ')
#endif
            ENDIF
         ENDDO
      ENDDO
#ifdef PARALLEL
#else
      J%NBlks=P-1
      J%NNon0=R-1

      ! Fill the upper triangle of J
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
#endif

    END SUBROUTINE MakeJ
!===============================================================================

!===============================================================================
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
!===============================================================================

!===============================================================================
#ifdef PARALLEL
     FUNCTION JBlock(Pair,PoleRoot) RESULT(JBlk)
#else
     FUNCTION JBlock(Pair,PoleRoot) RESULT(Jvct)
#endif
       TYPE(AtomPair)                           :: Pair
       TYPE(PoleNode), POINTER                  :: PoleRoot

       REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)  :: JBlk
#ifndef PARALLEL
       REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB)  :: Jvct
#endif
       REAL(DOUBLE),DIMENSION(0:SPLen)          :: SPBraC,SPBraS 
       REAL(DOUBLE),DIMENSION(3)                :: PTmp
       REAL(DOUBLE)                             :: ZetaA,ZetaB,EtaAB,EtaIn,    &
                                                   XiAB,ExpAB,CA,CB,CC,Ov,     &
                                                   PAx,PAy,PAz,PBx,PBy,PBz,    &
                                                   MDx,MDxy,MDxyz,Amp2,MaxAmp, &
                                                   Tau,OmegaMin,Px,Py,Pz
       REAL(DOUBLE)                             :: PExtent,EX
       REAL(DOUBLE)                             :: PStrength,Error
       INTEGER                                  :: KA,KB,CFA,CFB,PFA,PFB,      &
                                                   IndexA,IndexB,StartLA,      &
                                                   StartLB,StopLA,StopLB
       INTEGER                                  :: I,J,MaxLA,MaxLB,IA,IB,LMNA, &
                                                   LMNB,LA,LB,MA,MB,NA,NB,LAB, &
                                                   MAB,NAB,LM,LMN,Ell,EllA,    &
                                                   EllB,NC,L,M
!------------------------------------------------------------------------------- 
       JBlk=Zero
       KA=Pair%KA
       KB=Pair%KB
       Prim%A=Pair%A
       Prim%B=Pair%B
       Prim%AB2=Pair%AB2
       Prim%KA=Pair%KA
       Prim%KB=Pair%KB
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
             Prim%CFA=CFA             
             Prim%CFB=CFB
             Prim%Ell=MaxLA+MaxLB
             DO PFA=1,BS%NPFnc%I(CFA,KA)          
                DO PFB=1,BS%NPFnc%I(CFB,KB)
                   Prim%ZA=BS%Expnt%D(PFA,CFA,KA)
                   Prim%ZB=BS%Expnt%D(PFB,CFB,KB)
                   Prim%Zeta=Prim%ZA+Prim%ZB
                   Prim%Xi=Prim%ZA*Prim%ZB/Prim%Zeta
                   IF(TestPrimPair(Prim%Xi,Prim%AB2))THEN
                      Prim%PFA=PFA 
                      Prim%PFB=PFB
                      MaxAmp=SetBraBlok(Prim,BS)
!-------------------------------------------------------------------------------
!                     Compute maximal HG extent (for PAC) and Unsold esitmiate (for MAC)
!                     looping over all angular symmetries

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
!                           Extent (for PAC)
                            EX=Extent(EllA+EllB,Prim%Zeta,HGBra%D(:,IA,IB),TauPAC,ExtraEll_O=0,Potential_O=.TRUE.)
                            PExtent=MAX(PExtent,EX)
!                           Strength (for MAC)
                            CALL HGToSP(Prim%Zeta,EllA+EllB,HGBra%D(:,IA,IB),SPBraC,SPBraS)
!
                            PStrength = Zero
                            DO L=0,EllA+EllB
                               PStrength = PStrength+FudgeFactorial(L,SPEll+1)*Unsold0(L,SPBraC,SPBraS)
                            ENDDO
                            PStrength = (PStrength/TauMAC)**(Two/DBLE(SPEll+2))
                            IF(DP2 < PStrength) THEN
                               DP2   = PStrength
                            ENDIF
!
                         ENDDO
                      ENDDO 
                      DP2=MIN(1.D10,DP2)
!-------------------------------------------------------------------------------
!                     If finite compute ...
                      IF(PExtent>Zero.AND.DP2>Zero)THEN
!                        Initialize <KET|
                         CALL SetKet(Prim,PExtent)
#ifdef PERIODIC  
!                        WRAP the center of Phi_A(r) Phi_B(r+R) back into the box
                         CALL AtomCyclic(GM,Prim%P)
                         PTmp=Prim%P
!                        Sum over cells
                         DO NC=1,CS_IN%NCells
                            Prim%P=PTmp+CS_IN%CellCarts%D(:,NC)
                            PBox%Center=Prim%P
!                           Walk the walk
                            CALL JWalk(PoleRoot) 
                         ENDDO
                         Prim%P=PTmp
#else
!                        Walk the walk
                         CALL JWalk(PoleRoot)
#endif
!-------------------------------------------------------------------------------
!                        <BRA|KET>
                         IA = IndexA
                         DO LMNA=StartLA,StopLA
                            IA=IA+1
                            IB=IndexB
                            EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)
                            DO LMNB=StartLB,StopLB
                               IB=IB+1
                               EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)
                               Ell=EllA+EllB
!                              Near field
                               DO LMN=1,LHGTF(Ell)
                                  JBlk(IA,IB)=JBlk(IA,IB)+Phase%D(LMN)*HGBra%D(LMN,IA,IB)*HGKet(LMN)
                               ENDDO
!                              Far field
                               CALL HGToSP(Prim%Zeta,Ell,HGBra%D(:,IA,IB),SPBraC,SPBraS)
                               DO LM=0,LSP(Ell)
                                  JBlk(IA,IB)=JBlk(IA,IB)+SPBraC(LM)*SPKetC(LM)+SPBraS(LM)*SPKetS(LM)
                               ENDDO
                            ENDDO
                         ENDDO
#ifdef PERIODIC
!                       Calculate the FarField Multipole Contribution to the Matrix Element
!                       Contract the Primative MM  with the density MM
                        IF(GM%PBC%Dimen > 0) THEN
                           IA = IndexA
                           DO LMNA=StartLA,StopLA
                              IA=IA+1
                              IB=IndexB
                              DO LMNB=StartLB,StopLB  
                                 IB=IB+1
                                 ! ANOTHER PROBLEM WITH COPY IN COPY OUT HERE 
                                 JBlk(IA,IB)=JBlk(IA,IB)+CTraxFF(Prim,HGBra%D(:,IA,IB),GM)
                              ENDDO
                           ENDDO
                        ENDIF
#endif
                      ENDIF                     
                   ENDIF !End primitive thresholding
                ENDDO
             ENDDO
          ENDDO
       ENDDO
#ifdef PARALLEL
#else
       Jvct=BlockToVect(Pair%NA,Pair%NB,Jblk)
#endif

     END FUNCTION JBlock 
END MODULE JGen

