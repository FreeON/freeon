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
  USE PBCFarField
#ifdef PARALLEL
  USE FastMatrices
  USE ParallelQCTC
#endif
  IMPLICIT NONE
  LOGICAL PrintFlag
 
#ifdef PARALLEL
  INTEGER :: PrIndex,AbsIndex
#endif

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
      INTEGER                   :: NC
      REAL(DOUBLE),DIMENSION(3) :: B  
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
      DO AtA=1,NAtoms             
         DO AtB=1,NAtoms
            IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
               NAB = Pair%NA*Pair%NB
#ifndef PARALLEL
               ! Compute only the lower triangle of symmetric J
               IF(AtB<=AtA)THEN  
#endif
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

      IF(PrIndex /= TotPrCount) THEN
        WRITE(*,*) 'PrIndex = ',PrIndex
        WRITE(*,*) 'TotPrCount = ',TotPrCount
        STOP 'ERR: Serious counting problem I!'
      ENDIF
      IF(AbsIndex /= GlobalCount) THEN
        WRITE(*,*) 'AbsIndex = ',AbsIndex
        WRITE(*,*) 'GlobalCount = ',GlobalCount
        STOP 'ERR: Serious counting problem II!'
      ENDIF

      DO PrIndex = 1, TotPrCount
        ETimer(2) = ETimer(2) + PosTimePair%D(4,PrIndex)
      ENDDO
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
       REAL(DOUBLE)                             :: PStrength,Error,PiZ
       INTEGER                                  :: KA,KB,CFA,CFB,PFA,PFB,      &
                                                   IndexA,IndexB,StartLA,      &
                                                   StartLB,StopLA,StopLB
       INTEGER                                  :: I,J,MaxLA,MaxLB,IA,IB,LMNA, &
                                                   LMNB,LA,LB,MA,MB,NA,NB,LAB, &
                                                   MAB,NAB,LM,LMN,EllAB,EllA,    &
                                                   EllB,NC,L,M,LenHGTF,LenSP
#ifdef PARALLEL
       REAL(DOUBLE)                             :: ZA,ZB,T1,T2
       LOGICAL                                  :: FirstIterQ,OtherIterQ
#endif
       REAL(DOUBLE),EXTERNAL :: MondoTimer
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
#ifdef PARALLEL
                      ZA = Prim%ZA
                      ZB = Prim%ZB
                      Prim%P=(ZA*Prim%A+ZB*Prim%B)/Prim%Zeta
                      AbsIndex = AbsIndex + 1
                      IF(GQLineLoc == 0) then
                        FirstIterQ = AbsIndex >= BegPrInd%I(MyID) &
                                  .AND.AbsIndex <= EndPrInd%I(MyID)
                      ELSE
                        OtherIterQ = Prim%P(1) >= LC%D(1,MyID+1) .AND. &
                                   Prim%P(1) < RC%D(1,MyID+1) .AND. &
                                   Prim%P(2) >= LC%D(2,MyID+1) .AND. &
                                   Prim%P(2) < RC%D(2,MyID+1) .AND. &
                                   Prim%P(3) >= LC%D(3,MyID+1) .AND. &
                                   Prim%P(3) < RC%D(3,MyID+1)
                      ENDIF
                      IF(FirstIterQ .OR. OtherIterQ) THEN

                      PrIndex = PrIndex + 1
                      T1 = MondoTimer()
                      PosTimePair%D(1:3,PrIndex) = Prim%P(1:3)

#endif
                      Prim%PFA=PFA 
                      Prim%PFB=PFB
                      MaxAmp=SetBraBlok(Prim,BS)
#ifdef NewPAC
!                     Compute MAC and PAC
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
                            EllAB   = EllA+EllB
                            LenHGTF = LHGTF(EllAB)
                            LenSP   = LSP(EllAB)
                            PiZ     = (Pi/Prim%Zeta)**(ThreeHalves)
!                           Strength for MAC
                            CALL HGToSP_Direct(EllAB,LenHGTF,LenSP,PiZ,HGBra%D(1:LenHGTF,IA,IB), &
                                               SPBraC(0:LenSP),SPBraS(0:LenSP))
                            PStrength = Zero
                            DO L=0,EllAB
                               PStrength = PStrength+FudgeFactorial(L,SPEll+1)*Unsold0(L,SPBraC,SPBraS)
                            ENDDO
                            PStrength=(PStrength/TauMAC)**(Two/DBLE(SPEll+2))
                            IF(DP2<PStrength) THEN
                               DP2=PStrength
                            ENDIF
!                           Strength for PAC
                            IF(EllAB==0) THEN
                               PrimBeta = Prim%Zeta
                            ELSE
                               PrimBeta = GFactor*Prim%Zeta
                            ENDIF
                            PrimWCoef = PrimWCoef+MaxCoef(EllAB,Prim%Zeta,HGBra%D(1:LenHGTF,IA,IB)) &
                                                  *((Pi/PrimBeta)**(ThreeHalves))
                         ENDDO
                      ENDDO
                      DP2=MIN(1.D10,DP2)
                      IF(PrimWCoef>Zero .AND. DP2>Zero)THEN
#else
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
                            EllAB=EllA+EllB
                            LenSP=LSP(EllAB)
                            LenHGTF=LHGTF(EllAB)
                            PiZ=(Pi/Prim%Zeta)**(ThreeHalves)
!                           Extent (for PAC)
                            EX=Extent(EllAB,Prim%Zeta,HGBra%D(1:LenHGTF,IA,IB),TauPAC,ExtraEll_O=0,Potential_O=.TRUE.)
                            PExtent=MAX(PExtent,EX)
!                           Strength (for MAC)
                            CALL HGToSP_Direct(EllAB,LenHGTF,LenSP,PiZ,HGBra%D(1:LenHGTF,IA,IB), &
                                               SPBraC(0:LenSP),SPBraS(0:LenSP))
                            PStrength = Zero
                            DO L=0,EllAB
                               PStrength = PStrength+FudgeFactorial(L,SPEll+1)*Unsold0(L,SPBraC,SPBraS)
                            ENDDO
                            PStrength=(PStrength/TauMAC)**(Two/DBLE(SPEll+2))
                            IF(DP2<PStrength) THEN
                               DP2=PStrength
                            ENDIF
                         ENDDO
                      ENDDO 
                      DP2=MIN(1.D10,DP2)
!-------------------------------------------------------------------------------
!                     If finite compute ...
                      IF(PExtent>Zero.AND.DP2>Zero)THEN
#endif
!                        Initialize <KET|
                         CALL SetKet(Prim,PExtent)
!
!                        WRAP the center of Phi_A(r) Phi_B(r+R) back into the box
!                        wrapping has to be turned off for now inorder for the
!                        lattice forces to be correct
!                        CALL AtomCyclic(GM,Prim%P)
!
                         PTmp(1)=Prim%P(1)
                         PTmp(2)=Prim%P(2)
                         PTmp(3)=Prim%P(3)
!                        Sum over cells
                         DO NC=1,CS_IN%NCells
                            Prim%P(1)=PTmp(1)+CS_IN%CellCarts%D(1,NC)
                            Prim%P(2)=PTmp(2)+CS_IN%CellCarts%D(2,NC)
                            Prim%P(3)=PTmp(3)+CS_IN%CellCarts%D(3,NC)
                            PBox%Center(1)=Prim%P(1)
                            PBox%Center(2)=Prim%P(2)
                            PBox%Center(3)=Prim%P(3)
!                           Walk the walk
                            CALL JWalk(PoleRoot) 
                         ENDDO
                         Prim%P(1)=PTmp(1)
                         Prim%P(2)=PTmp(2)
                         Prim%P(3)=PTmp(3)
!-------------------------------------------------------------------------------
!                        <BRA|KET>
                         IA = IndexA
                         DO LMNA=StartLA,StopLA
                            IA=IA+1
                            IB=IndexB
                            EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)
                            DO LMNB=StartLB,StopLB
                               IB=IB+1
                               EllB   = BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)
                               EllAB  = EllA+EllB
                               LenHGTF= LHGTF(EllAB)
                               LenSP  = LSP(EllAB)
                               PiZ    = (Pi/Prim%Zeta)**(ThreeHalves)
!                              Near field
                               DO LMN=1,LenHGTF
                                  JBlk(IA,IB)=JBlk(IA,IB)+Phase%D(LMN)*HGBra%D(LMN,IA,IB)*HGKet(LMN)
                               ENDDO
!                              Far field
                               CALL HGToSP_Direct(EllAB,LenHGTF,LenSP,PiZ,HGBra%D(1:LenHGTF,IA,IB), &
                                                  SPBraC(0:LenSP),SPBraS(0:LenSP))
                               DO LM=0,LenSP
                                  JBlk(IA,IB)=JBlk(IA,IB)+SPBraC(LM)*SPKetC(LM)+SPBraS(LM)*SPKetS(LM)
                               ENDDO
                            ENDDO
                         ENDDO
!                        Calculate the FarField Multipole Contribution to the Matrix Element
!                        Contract the Primative MM  with the density MM
                         IF(GM%PBC%Dimen > 0) THEN
                            IA = IndexA
                            DO LMNA=StartLA,StopLA
                               IA=IA+1
                               IB=IndexB
                               DO LMNB=StartLB,StopLB  
                                  IB=IB+1
                                  JBlk(IA,IB)=JBlk(IA,IB)+CTraxFF(Prim,HGBra%D(:,IA,IB),GM)
                               ENDDO
                            ENDDO
                         ENDIF
                      ENDIF
#ifdef PARALLEL
                      T2 = MondoTimer()
                      PosTimePair%D(4,PrIndex) = T2-T1
                      ENDIF
#endif

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

