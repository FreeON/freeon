!    FAST O(Lg N), BLOKWISE ACCUMULATION OF Tr{P_(A,B).dJ^T_(A,B)/dA}
!    Author: Matt Challacombe 
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
  USE PBCFarField
  IMPLICIT NONE
  LOGICAL PrintFlag
  ! Some extra temporary ds
  TYPE(ScalarHerm) :: TempHerm
  REAL(DOUBLE),ALLOCATABLE,DIMENSION(:) :: NNearCount
  !
CONTAINS 
  !=======================================================================================
  !
  !=======================================================================================
  FUNCTION TrPdJ(Pair,P,GMLoc) RESULT(Vck)
    TYPE(AtomPair)                             :: Pair
    TYPE(CRDS)                                 :: GMLoc
    TYPE(QCPrim)                               :: QP
    !
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)    :: P
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB,15) :: dJ
    REAL(DOUBLE),DIMENSION(3)                  :: PTmp,PShft
    REAL(DOUBLE),DIMENSION(0:SPLen)            :: SPBraC,SPBraS 
    REAL(DOUBLE)                               :: ZetaA,ZetaB,EtaAB,EtaIn,    &
         XiAB,ExpAB,CA,CB,CC,Ov,     &
         PAx,PAy,PAz,PBx,PBy,PBz,    &
         MDx,MDxy,MDxyz,Amp2,MaxAmp, &
         Pab,JNorm,Tau,OmegaMin,     &
         PExtent,PStrength,Ext
    INTEGER                                    :: KA,KB,CFA,CFB,PFA,PFB,AtA,AtB,    &
         IndexA,IndexB,              &
         StartLA,StartLB,            &
         StopLA,StopLB
    INTEGER                                    :: I,J,K,L,MaxLA,MaxLB,IA,IB,  &
         LMNA,LMNB,LA,LB,MA,MB,      &
         NA,NB,LAB,MAB,NAB,LM,LMN,   &
         EllA,EllB,EllAB,LenHGTF,    &
         LenSP,KI,LenAB
    REAL(DOUBLE), EXTERNAL                     :: BlkTrace_2 
    INTEGER                                    :: NC,NNearTmp
    REAL(DOUBLE)                               :: Px,Py,Pz,HGFac,PiZ
    !
    REAL(DOUBLE),DIMENSION(15)                 :: Vck
    REAL(DOUBLE),DIMENSION(3)                  :: nlm
    REAL(DOUBLE),DIMENSION(0:SPLen,3)          :: SPKetC_L,SPKetS_L
    REAL(DOUBLE),DIMENSION(1:HGLen,3)          :: HGKet_L
    REAL(DOUBLE),DIMENSION(1:HGLen)            :: HGKet_old
    REAL(DOUBLE),DIMENSION(0:SPLen)            :: SPKetC_old,SPKetS_old
    !----------------------------------------------------------------------------------------- 
    QP%Prim%A=Pair%A
    QP%Prim%B=Pair%B
    QP%Prim%KA=Pair%KA
    QP%Prim%KB=Pair%KB
    QP%Prim%AB2=Pair%AB2
    !----------------------------------
    KA=QP%Prim%KA
    KB=QP%Prim%KB
    dJ=Zero
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
          QP%Prim%CFA=CFA             
          QP%Prim%CFB=CFB
          QP%Prim%Ell=MaxLA+MaxLB+1 
          TempHerm%Ell=QP%Prim%Ell
          LenAB=LHGTF(QP%Prim%Ell)
          !---------------------------------- 
          DO PFA=1,BS%NPFnc%I(CFA,KA)          
             DO PFB=1,BS%NPFnc%I(CFB,KB)
                !----------------------------------
                QP%Prim%ZA=BS%Expnt%D(PFA,CFA,KA)
                QP%Prim%ZB=BS%Expnt%D(PFB,CFB,KB)
                QP%Prim%Zeta=QP%Prim%ZA+QP%Prim%ZB
                QP%Prim%Xi=QP%Prim%ZA*QP%Prim%ZB/QP%Prim%Zeta
                IF(TestPrimPair(QP%Prim%Xi,QP%Prim%AB2))THEN
                   QP%Prim%PFA=PFA 
                   QP%Prim%PFB=PFB
                   !
                   MaxAmp=SetBraBlok(QP%Prim,BS,Gradients_O=.FALSE.)

!!$                   QP%Prim%P=AtomToFrac(GM,QP%Prim%P)
!!$                   CALL FracCyclic(GM,QP%Prim%P)
!!$                   QP%Prim%P=FracToAtom(GM,QP%Prim%P)

                   !
                   TempHerm%Coef(:)=Zero
                   TempHerm%Zeta=QP%Prim%Zeta

                   IA = IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1
                      IB=IndexB
                      EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)                         
                      DO LMNB=StartLB,StopLB
                         IB=IB+1
                         EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)  
                         EllAB   = EllA+EllB+1
                         LenHGTF = LHGTF(EllAB)
                         DO K=1,3
                            DO I=1,LenHGTF
                               TempHerm%Coef(I)=MAX(TempHerm%Coef(I),ABS(P(IA,IB)*dHGBra%D(I,IA,IB,K)))
                            ENDDO
                         ENDDO
                      ENDDO
                   ENDDO
                   !
                   CALL SetSerialPAC(QP%PAC,TempHerm)                   
                   CALL SetSerialMAC(QP%MAC,TempHerm)                   
                   ! The integral estimate (ab|ab)^(1/2)
                   QP%IHalf=Estimate(QP%Prim%Ell,QP%Prim%Zeta,TempHerm%Coef(1:LenAB))
                   IF(QP%IHalf<TauTwo*1D-5)CYCLE
                   ! Initialize <KET|
                   LenSP=LSP(QP%Prim%Ell)
                   LenHGTF=LHGTF(QP%Prim%Ell)
                   CALL DBL_VECT_EQ_DBL_SCLR(LenHGTF,HGKet(1),Zero)
                   CALL DBL_VECT_EQ_DBL_SCLR(LenSP+1,SPKetC(0),Zero)
                   CALL DBL_VECT_EQ_DBL_SCLR(LenSP+1,SPKetS(0),Zero)   
                   DO K=1,3
                      CALL DBL_VECT_EQ_DBL_SCLR(LenHGTF,HGKet_L(1,K),Zero)   
                      CALL DBL_VECT_EQ_DBL_SCLR(LenSP+1,SPKetC_L(0,K),Zero)   
                      CALL DBL_VECT_EQ_DBL_SCLR(LenSP+1,SPKetS_L(0,K),Zero)   
                   ENDDO
                   !
                   PTmp=QP%Prim%P
                   !
                   DO NC=1,CS_IN%NCells 
                      QP%Prim%P=PTmp+CS_IN%CellCarts%D(:,NC) 
                      !                     Calculate the fractional coordinates
                      nlm =  AtomToFrac(GMLoc,CS_IN%CellCarts%D(:,NC))
                      !                     Store the Old Stuff
                      HGKet_old(1:LenHGTF)= HGKet(1:LenHGTF)
                      SPKetC_old(0:LenSP) = SPKetC(0:LenSP)
                      SPKetS_old(0:LenSP) = SPKetS(0:LenSP)
                      !                     Walk the walk
                      NNearTmp=NNearAv
                      NNearAv=0
                      CALL JWalk2(QP,PoleRoot) 
                      NNearCount(NC)=NNearCount(NC)+NNearAv
                      NNearAv=NNearTmp+NNearAv
                      NPrim=NPrim+1
                      !                     Acumulate the Lattice Forces
                      DO I=1,3
                         HGKet_L(1:LenHGTF,I)= HGKet_L(1:LenHGTF,I)+ nlm(I)*(HGKet(1:LenHGTF)- HGKet_old(1:LenHGTF))
                         SPKetC_L(0:LenSP,I) = SPKetC_L(0:LenSP,I) + nlm(I)*(SPKetC(0:LenSP) - SPKetC_old(0:LenSP))
                         SPKetS_L(0:LenSP,I) = SPKetS_L(0:LenSP,I) + nlm(I)*(SPKetS(0:LenSP) - SPKetS_old(0:LenSP))
                      ENDDO
                   ENDDO
                   !
                   QP%Prim%P=PTmp
                   ! Contract <Bra|Ket> bloks to compute matrix elements of J : SameAtom=.FALSE.
                   IA = IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1
                      IB=IndexB
                      EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)                         
                      DO LMNB=StartLB,StopLB
                         IB=IB+1
                         EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)                         
                         EllAB  = EllA+EllB+1
                         LenHGTF= LHGTF(EllAB)
                         LenSP  = LSP(EllAB)
                         PiZ    = (Pi/QP%Prim%Zeta)**(ThreeHalves)
                         DO K=1,3
                            KI = K
                            DO LMN=1,LenHGTF
                               dJ(IA,IB,KI) = dJ(IA,IB,KI)+Phase%D(LMN)*dHGBra%D(LMN,IA,IB,K)*HGKet(LMN)
                            ENDDO
                            CALL HGToSP_Direct(EllAB,LenHGTF,LenSP,PiZ,dHGBra%D(1:LenHGTF,IA,IB,K), &
                                 SPBraC(0:LenSP),SPBraS(0:LenSP))
                            DO LM=0,LenSP
                               dJ(IA,IB,KI) = dJ(IA,IB,KI)+SPBraC(LM)*SPKetC(LM)+SPBraS(LM)*SPKetS(LM)
                            ENDDO

                            dJ(IA,IB,KI) = dJ(IA,IB,KI)+CTraxFF(QP%Prim,dHGBra%D(1:1,IA,IB,K),GMLoc) 

                         ENDDO
                      ENDDO
                   ENDDO
                   ! Contract <Bra|Ket> bloks to compute matrix elements of J : SameAtom=.TRUE.
                   MaxAmp=SetBraBlok(QP%Prim,BS,Gradients_O=.TRUE.)
                   !
                   IA = IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1
                      IB=IndexB
                      EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)                         
                      DO LMNB=StartLB,StopLB
                         IB=IB+1
                         EllB=BS%LxDex%I(LMNB)+BS%LyDex%I(LMNB)+BS%LzDex%I(LMNB)   
                         EllAB  = EllA+EllB+1
                         LenHGTF= LHGTF(EllAB)
                         LenSP  = LSP(EllAB)
                         PiZ    = (Pi/QP%Prim%Zeta)**(ThreeHalves)
                         DO K=1,3
!                            CALL HGToSP77(EllAB,PiZ,dHGBra%D(1,IA,IB,K),SPBraC(0),SPBraS(0))
                            CALL HGToSP_Direct(EllAB,LenHGTF,LenSP,PiZ,dHGBra%D(1:LenHGTF,IA,IB,K), &
                                               SPBraC(0:LenSP),SPBraS(0:LenSP))                            

                            KI = K+3
                            DO LMN=1,LenHGTF
                               dJ(IA,IB,KI) = dJ(IA,IB,KI)+Phase%D(LMN)*dHGBra%D(LMN,IA,IB,K)*HGKet(LMN)
                            ENDDO

                            DO LM=0,LenSP
                               dJ(IA,IB,KI) = dJ(IA,IB,KI)+SPBraC(LM)*SPKetC(LM)+SPBraS(LM)*SPKetS(LM)
                            ENDDO

                            dJ(IA,IB,KI) = dJ(IA,IB,KI)+CTraxFF(QP%Prim,dHGBra%D(:,IA,IB,K),GMLoc)  

                            DO I=1,3
                               KI=K+3*(I+1)
                               IF(GMLoc%PBC%AutoW%I(I)==1 .AND. GMLoc%PBC%AutoW%I(K)==1) THEN 
                                  DO LMN=1,LenHGTF
                                     dJ(IA,IB,KI)=dJ(IA,IB,KI)+Phase%D(LMN)*dHGBra%D(LMN,IA,IB,K)*HGKet_L(LMN,I)
                                  ENDDO
                                  DO LM=0,LenSP
                                     dJ(IA,IB,KI)=dJ(IA,IB,KI)+SPBraC(LM)*SPKetC_L(LM,I)+SPBraS(LM)*SPKetS_L(LM,I)
                                  ENDDO
                               ENDIF
                            ENDDO

                         ENDDO
                      ENDDO
                   ENDDO

                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    DO K=1,15
       Vck(K)=BlkTrace_2(Pair%NA,Pair%NB,P,TRANSPOSE(dJ(:,:,K)))
    ENDDO


!    STOP
!!$
    !
  END FUNCTION TrPdJ
  !====================================================================================================
  !
  !====================================================================================================
  FUNCTION dNukE(GMLoc,At) RESULT(Vck)
    TYPE(CRDS)                        :: GMLoc
    TYPE(QCPrim)                      :: QP
    REAL(DOUBLE)                      :: Tau,NukeCo,NukePole,PExtent
    REAL(DOUBLE),DIMENSION(4)         :: dBra
    REAL(DOUBLE),DIMENSION(0:SPLen)   :: SPBraC,SPBraS 
    INTEGER                           :: At,LenSP,LenHGTF,I,J,K,LM,LMN
    INTEGER                           :: NC,KI
    REAL(DOUBLE),DIMENSION(3)         :: PTmp,nlm
    !
    REAL(DOUBLE),DIMENSION(15)        :: Vck
    REAL(DOUBLE),DIMENSION(0:SPLen,3) :: SPKetC_L,SPKetS_L
    REAL(DOUBLE),DIMENSION(1:HGLen,3) :: HGKet_L
    REAL(DOUBLE),DIMENSION(1:HGLen)   :: HGKet_old
    REAL(DOUBLE),DIMENSION(0:SPLen)   :: SPKetC_old,SPKetS_old
    !---------------------------------------------------------------------------------------------
    ! Initialize (dBRA|
    NukeCo=-GMLoc%AtNum%D(At)*(NuclearExpnt/Pi)**(ThreeHalves)
    DO K=1,3
       dHGBra%D(1:4,1,1,K)=Zero
    ENDDO
    dHGBra%D(2,1,1,1)=NukeCo
    dHGBra%D(3,1,1,2)=NukeCo
    dHGBra%D(4,1,1,3)=NukeCo
    ! Initialize the primitive          
    QP%Prim%Ell=1
    QP%Prim%P=GMLoc%Carts%D(:,At)
    QP%Prim%Zeta=NuclearExpnt
    QP%PAC%Zeta=QP%Prim%Zeta
    QP%PAC%Wght=GMLoc%AtNum%D(At)
    QP%MAC%O(0)=GMLoc%AtNum%D(At)
    QP%MAC%Delta=Zero
    QP%IHalf=ABS(NukeCo)  !! ??????????
    ! Initialize the |KET)
    HGKet(1:4)=Zero
    SPKetC(0:3)=Zero
    SPKetS(0:3)=Zero
    !
    HGKet_L(:,:)  = Zero
    SPKetC_L(:,:) = Zero
    SPKetS_L(:,:) = Zero
    !

    PTmp=GMLoc%Carts%D(:,At)

    DO NC=1,CS_IN%NCells
       !         Set atomic "primitive"
       QP%Prim%P=PTmp+CS_IN%CellCarts%D(:,NC)
       !         Calculate the fractional coordinates
       nlm = AtomToFrac(GMLoc,CS_IN%CellCarts%D(:,NC))
       !         Store the Old Stuff
       HGKet_old(1:4)= HGKet(1:4)
       SPKetC_old(0:3) = SPKetC(0:3)
       SPKetS_old(0:3) = SPKetS(0:3)
       !         Walk the walk
       CALL JWalk2(QP,PoleRoot,Nucular_O=.TRUE.)
       !         Acumulate the Lattice Forces
       DO I=1,3
          HGKet_L(1:4,I)= HGKet_L(1:4,I)+ nlm(I)*(HGKet(1:4)- HGKet_old(1:4))
          SPKetC_L(0:3,I) = SPKetC_L(0:3,I) + nlm(I)*(SPKetC(0:3) - SPKetC_old(0:3))
          SPKetS_L(0:3,I) = SPKetS_L(0:3,I) + nlm(I)*(SPKetS(0:3) - SPKetS_old(0:3))
       ENDDO
    ENDDO
    ! Reset the Atomic Coordinates
    QP%Prim%P=PTmp
    ! Init bra xforms
    Vck=Zero
    DO K=1,3
       DO LMN=1,4
          Vck(K)=Vck(K)+Phase%D(LMN)*dHGBra%D(LMN,1,1,K)*HGKet(LMN)
       ENDDO
       CALL HGToSP(QP%Prim%Zeta,1,dHGBra%D(:,1,1,K),SPBraC,SPBraS)
       DO LM=0,3
          Vck(K)=Vck(K)+SPBraC(LM)*SPKetC(LM)+SPBraS(LM)*SPKetS(LM)
       ENDDO
    ENDDO
    ! Add in the Far Field, Dipole and Quadripole Correction
    IF(GMLoc%PBC%Dimen>0) THEN
       DO K=1,3
          Vck(K)=Vck(K)+CTraxFF(QP%Prim,dHGBra%D(:,1,1,K),GMLoc)
       ENDDO
    ENDIF
    !      Inner Sum Nuc Box Force
    DO K=1,3
       CALL HGToSP(QP%Prim%Zeta,1,dHGBra%D(:,1,1,K),SPBraC,SPBraS)
       DO I=1,3
          KI = K + 3*(I+1)
          IF(GMLoc%PBC%AutoW%I(I)==1 .AND. GMLoc%PBC%AutoW%I(K)==1) THEN
             DO LMN=1,4
                Vck(KI)=Vck(KI)+Phase%D(LMN)*dHGBra%D(LMN,1,1,K)*HGKet_L(LMN,I)
             ENDDO
             DO LM=0,3
                Vck(KI)=Vck(KI)+SPBraC(LM)*SPKetC_L(LM,I)+SPBraS(LM)*SPKetS_L(LM,I)
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    !
  END FUNCTION dNukE
END MODULE BlokTrPdJ


