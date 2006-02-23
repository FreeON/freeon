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
  IMPLICIT NONE
  LOGICAL PrintFlag
  ! Some extra temporary ds
  TYPE(ScalarHerm) :: TempHerm
  REAL(DOUBLE),ALLOCATABLE,DIMENSION(:) :: NNearCount
  !-------------------------------------------------------------------------------
CONTAINS
  !===============================================================================

  !===============================================================================
  SUBROUTINE MakeJ(J)
    TYPE(BCSR)                :: J
    TYPE(DBL_RNK2)            :: Temp
    TYPE(AtomPair)            :: Pair
    INTEGER                   :: AtA,AtB
    INTEGER                   :: JP,K,NA,NB,NAB,P,Q,R,I1,I2,I3,L,I
    INTEGER                   :: NC
    REAL(DOUBLE),DIMENSION(3) :: B  
    !------------------------------------------------------------------------------- 
    !    
    ALLOCATE(TempHerm%Coef(1:HGLen))
    !
    !     Initialize the matrix and associated indecies
    P=1
    R=1
    J%RowPt%I(1)=1
    CALL SetEq(J%MTrix,Zero)
    J%NAtms= NAtoms
    ! Loop over atom pairs
    DO AtA=1,NAtoms             
       DO AtB=1,NAtoms
          IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN            
             NAB = Pair%NA*Pair%NB
             IF(AtB<=AtA)THEN  
                B = Pair%B
                DO NC=1,CS_OUT%NCells
                   Pair%B   = B+CS_OUT%CellCarts%D(:,NC)
                   Pair%AB2 = (Pair%A(1)-Pair%B(1))**2+(Pair%A(2)-Pair%B(2))**2+(Pair%A(3)-Pair%B(3))**2
                   IF(TestAtomPair(Pair)) THEN
                      J%MTrix%D(R:R+NAB-1)=J%MTrix%D(R:R+NAB-1)+Two*JBlock(Pair,PoleRoot)
                   ENDIF
                ENDDO
             ENDIF
             J%ColPt%I(P)=AtB
             J%BlkPt%I(P)=R
             R=R+NAB
             P=P+1  
             J%RowPt%I(AtA+1)=P        
             IF(R>MaxNon0.OR.P>MaxBlks)THEN
                CALL Halt(' BCSR dimensions blown in J ')
             ENDIF
          ENDIF
       ENDDO
    ENDDO
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
    DEALLOCATE(TempHerm%Coef)
    !
  END SUBROUTINE MakeJ
  !===============================================================================
  !
  !===============================================================================
  FUNCTION JBlock(Pair,PoleRoot) RESULT(Jvct)
    TYPE(AtomPair)                            :: Pair
    TYPE(PoleNode), POINTER                   :: PoleRoot    
    TYPE(QCPrim)                              :: QP
    TYPE(PAC)                                 :: TmpPAC
    REAL(DOUBLE),DIMENSION(Pair%NA*Pair%NB)   :: Jvct
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)   :: JBlk
    REAL(DOUBLE),DIMENSION(0:SPLen)           :: SPBraC,SPBraS 
    REAL(DOUBLE),DIMENSION(3)                 :: PTmp
    REAL(DOUBLE)                              :: ZetaA,ZetaB,EtaAB,EtaIn,  X
    REAL(DOUBLE)                              :: XiAB,ExpAB,CA,CB,CC,Ov
    REAL(DOUBLE)                              :: PAx,PAy,PAz,PBx,PBy,PBz
    REAL(DOUBLE)                              :: MDx,MDxy,MDxyz,Amp2,MaxAmp
    REAL(DOUBLE)                              :: Tau,OmegaMin,Px,Py,Pz
    REAL(DOUBLE)                              :: BraEstimate
    REAL(DOUBLE)                              :: PStrength,Error,PiZ
    INTEGER                                   :: KA,KB,CFA,CFB,PFA,PFB
    INTEGER                                   :: IndexA,IndexB,StartLA
    INTEGER                                   :: StartLB,StopLA,StopLB
    INTEGER                                   :: I,J,MaxLA,MaxLB,IA,IB,LMNA
    INTEGER                                   :: LMNB,LA,LB,MA,MB,NA,NB,LAB
    INTEGER                                   :: MAB,NAB,LM,LMN,EllAB,EllA,LenAB
    INTEGER                                   :: EllB,NC,L,M,LenHGTF,LenSP ,NNearTmp
    REAL(DOUBLE),EXTERNAL                     :: MondoTimer
    !------------------------------------------------------------------------------- 
    JBlk=Zero
    KA=Pair%KA
    KB=Pair%KB
    QP%Prim%A=Pair%A
    QP%Prim%B=Pair%B
    QP%Prim%AB2=Pair%AB2
    QP%Prim%KA=Pair%KA
    QP%Prim%KB=Pair%KB
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
          QP%Prim%CFA=CFA             
          QP%Prim%CFB=CFB
          QP%Prim%Ell=MaxLA+MaxLB
          TempHerm%Ell=QP%Prim%Ell
          LenAB=LHGTF(QP%Prim%Ell)

          DO PFA=1,BS%NPFnc%I(CFA,KA)          
             DO PFB=1,BS%NPFnc%I(CFB,KB)
                QP%Prim%ZA=BS%Expnt%D(PFA,CFA,KA)
                QP%Prim%ZB=BS%Expnt%D(PFB,CFB,KB)
                QP%Prim%Zeta=QP%Prim%ZA+QP%Prim%ZB
                QP%Prim%Xi=QP%Prim%ZA*QP%Prim%ZB/QP%Prim%Zeta
                IF(TestPrimPair(QP%Prim%Xi,QP%Prim%AB2))THEN
                   QP%Prim%PFA=PFA 
                   QP%Prim%PFB=PFB
                   MaxAmp=SetBraBlok(QP%Prim,BS)
                   TempHerm%Zeta=QP%Prim%Zeta
                   TempHerm%Coef(1:LenAB)=Zero
                   IA = IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1
                      IB=IndexB
                      EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)                         
                      DO LMNB=StartLB,StopLB
                         IB=IB+1
                         EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)  
                         LenHGTF=LHGTF(EllA+EllB)
                         DO I=1,LenHGTF
                            TempHerm%Coef(I)=MAX(TempHerm%Coef(I),ABS(HGBra%D(I,IA,IB)))
                         ENDDO
                     ENDDO
                   ENDDO
                   !
                   CALL SetSerialPAC(QP%PAC,TempHerm)                   
                   ! Third level of primitive screening, based on Holders and 
                   ! Cramers inequalities.  Maybe this argues finally for a pair list...

!                   IF(QP%PAC%Wght/QP%PAC%Zeta<DBLE(NEl)*TauTwo)CYCLE
                   !
                   CALL SetSerialMAC(QP%MAC,TempHerm)                   
                   ! The integral estimate (ab|ab)^(1/2)
                   QP%IHalf=Estimate(QP%Prim%Ell,QP%Prim%Zeta,TempHerm%Coef(1:LenAB))
                   IF(QP%IHalf<TauTwo*1D-5)CYCLE
#ifdef PAC_DEBUG
                   DO L=1,LenAB
                      ERRBRA(L)=TempHerm%Coef(L)
                   ENDDO
#endif
#ifdef MAC_DEBUG
                   CALL HGToSP_Bra(QP%Prim%Zeta,QP%Prim%Ell,TempHerm%Coef,SPErrorBraC,SPErrorBraS)
#endif
                   ! Initialize <KET|
                   CALL DBL_VECT_EQ_DBL_SCLR(HGLen,HGKet(1),Zero)
                   CALL DBL_VECT_EQ_DBL_SCLR(SPLen+1,SPKetC(0),Zero)
                   CALL DBL_VECT_EQ_DBL_SCLR(SPLen+1,SPKetS(0),Zero)   
                  ! Sum over cells
                   PTmp=QP%Prim%P
                   DO NC=1,CS_IN%NCells
                      QP%Prim%P=PTmp+CS_IN%CellCarts%D(:,NC)
                      NNearTmp=NNearAv
                      NNearAv=0
                      CALL JWalk2(QP,PoleRoot) 
                      NNearCount(NC)=NNearCount(NC)+NNearAv
                      NNearAv=NNearTmp+NNearAv
                      NPrim=NPrim+1
                   ENDDO
                   QP%Prim%P=PTmp 
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
                         PiZ    = (Pi/QP%Prim%Zeta)**(ThreeHalves)
                         !  Near field
                         DO LMN=1,LenHGTF
                            JBlk(IA,IB)=JBlk(IA,IB)+Phase%D(LMN)*HGBra%D(LMN,IA,IB)*HGKet(LMN)
                         ENDDO
                         !  Far field
                         CALL HGToSP_Direct(EllAB,LenHGTF,LenSP,PiZ,HGBra%D(1:LenHGTF,IA,IB), &
                                            SPBraC(0:LenSP),SPBraS(0:LenSP))
                         DO LM=0,LenSP
                            JBlk(IA,IB)=JBlk(IA,IB)+SPBraC(LM)*SPKetC(LM)+SPBraS(LM)*SPKetS(LM)
                         ENDDO
                      ENDDO
                   ENDDO
                   ! Calculate the FarField Multipole Contribution to the Matrix Element
                   ! Contract the Primative MM  with the density MM
                   IF(GM%PBC%Dimen > 0) THEN
                      IA = IndexA
                      DO LMNA=StartLA,StopLA
                         IA=IA+1
                         IB=IndexB
                         DO LMNB=StartLB,StopLB  
                            IB=IB+1
                            JBlk(IA,IB)=JBlk(IA,IB)+CTraxFF(QP%Prim,HGBra%D(:,IA,IB),GM)
                         ENDDO
                      ENDDO
                   ENDIF
                ENDIF !End primitive thresholding
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    Jvct=BlockToVect(Pair%NA,Pair%NB,Jblk)
  END FUNCTION JBlock
END MODULE JGen

