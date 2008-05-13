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
  TYPE(DBL_RNK2)::LatFrc_DWRAP
  !
CONTAINS 
  !=======================================================================================
  !
  !=======================================================================================
  !=======================================================================================
  !
  !=======================================================================================
  SUBROUTINE TrPdJ2(Pair,P,GMLoc,AtomForce,LattForce,NoWrap_O)
    TYPE(AtomPair)                             :: Pair
    TYPE(CRDS)                                 :: GMLoc
    TYPE(QCPrim)                               :: QP
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB)    :: P
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB,3)  :: HGGradA,HGGradP
    REAL(DOUBLE),DIMENSION(Pair%NA,Pair%NB,3,3):: WrapFrc

    REAL(DOUBLE)                               :: dWdL,Rx,Ry,Rz
    REAL(DOUBLE),DIMENSION(1:3,1:3)            :: LattForce,Delta,tmp
    LOGICAL, OPTIONAL                          :: NoWrap_O
    LOGICAL                                    :: NoWrap
    REAL(DOUBLE),DIMENSION(3)                  :: PTmp,PTmp2
    REAL(DOUBLE),DIMENSION(0:SPLen)            :: SPBraC,SPBraS 
    REAL(DOUBLE)                               :: ZetaA,ZetaB,EtaP,EtaIn,    &
         XiP,ExpP,CA,CB,CC,Ov,     &
         PAx,PAy,PAz,PBx,PBy,PBz,    &
         MDx,MDxy,MDxyz,Amp2,MaxAmp, &
         Pab,JNorm,Tau,OmegaMin,     &
         PExtent,PStrength,Ext,LTmp
    INTEGER                                    :: KA,KB,CFA,CFB,PFA,PFB,AtA,AtB,    &
         IndexA,IndexB,              &
         StartLA,StartLB,            &
         StopLA,StopLB
    INTEGER                                    :: I,J,K,L,MaxLA,MaxLB,IA,IB,  &
         LMNA,LMNB,LA,LB,MA,MB,      &
         NA,NB,LP,MP,NP,LM,LMN,      &
         EllA,EllB,EllP,LenHGTF,     &
         LenSP,KI,LenD,EllG,LenG,LenGSP
    REAL(DOUBLE), EXTERNAL                     :: BlkTrace_2 
    INTEGER                                    :: M,N,NC,NNearTmp
    REAL(DOUBLE)                               :: Px,Py,Pz,HGFac,PiZ,TFA,TFP,FFA,FFP, &
         dPdA,dPdB,E000,EIJ0
    !
    REAL(DOUBLE),DIMENSION(3)                  :: nlm,dAdL,dBdL,dCdL,dPwdL,dTFdPC
    REAL(DOUBLE),DIMENSION(1:3)                :: AtomForce,AtomA,AtomP,dEDdP,PBCNearFieldA,PBCNearFieldP
    REAL(DOUBLE),DIMENSION(1:3)                :: TinFoilA
    REAL(DOUBLE),DIMENSION(0:SPLen,3)          :: SPKetC_L,SPKetS_L
    REAL(DOUBLE),DIMENSION(1:HGLen,3)          :: HGKet_L,BraGradA,BraGradP,dMcMrchA
    REAL(DOUBLE),DIMENSION(1:HGLen)            :: EIJ,HGKet_old
    REAL(DOUBLE),DIMENSION(0:SPLen)            :: SPKetC_old,SPKetS_old
    REAL(DOUBLE),DIMENSION(6,3)                :: dDdP
    LOGICAL :: PP
    !----------------------------------------------------------------------------------------- 
    IF(PRESENT(NoWrap_O))THEN
       NoWrap=NoWrap_O
    ELSE
       NoWrap=.FALSE.
    ENDIF
    !
    QP%Prim%A=Pair%A
    QP%Prim%B=Pair%B
    QP%Prim%KA=Pair%KA
    QP%Prim%KB=Pair%KB
    QP%Prim%AB2=Pair%AB2
    ! Basic assumption about P(1,1) may be wrong here!! NEED MORE SOPHISTICATED ADDRESSING LATER
    !----------------------------------
    KA=QP%Prim%KA
    KB=QP%Prim%KB
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
          EllP=MaxLA+MaxLB
          EllG=EllP+1
          !
          QP%Prim%CFA=CFA             
          QP%Prim%CFB=CFB
          QP%Prim%Ell=EllG
          !
          LenD=LHGTF(EllP) ! Length of the HGTF distribution
          LenG=LHGTF(EllG) ! Length of the HGTF distributions gradient
          LenGSP=LSP(ELLG) ! Length of the SP distributions gradient
          !---------------------------------- 
          DO PFA=1,BS%NPFnc%I(CFA,KA)          
             DO PFB=1,BS%NPFnc%I(CFB,KB)
                !----------------------------------
                QP%Prim%ZA=BS%Expnt%D(PFA,CFA,KA)
                QP%Prim%ZB=BS%Expnt%D(PFB,CFB,KB)
                QP%Prim%Zeta=QP%Prim%ZA+QP%Prim%ZB
                QP%Prim%Xi=QP%Prim%ZA*QP%Prim%ZB/QP%Prim%Zeta
                IF(TestPrimPair(QP%Prim%Xi,QP%Prim%AB2))THEN

!!$                   WRITE(*,*)' ================================================='
!!$                   WRITE(*,*)' A = ',Pair%A(1),' B = ',Pair%B(1)
!!$                   WRITE(*,*)' CFA=', CFA,' CFB = ',CFB
!!$                   WRITE(*,*)' PFA=',PFA,' PFB = ',PFB

                   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                   ! Construct the local primitive here, then wrap and center it
                   QP%Prim%PFA=PFA 
                   QP%Prim%PFB=PFB
                   PiZ=(Pi/QP%Prim%Zeta)**(ThreeHalves)
                   QP%Prim%P=(QP%Prim%ZA*QP%Prim%A+QP%Prim%ZB*QP%Prim%B)/QP%Prim%Zeta
                   CALL PWrap(GMLoc,QP%Prim,.NOT.NoWrap)                   !
                   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                   ! Compute the HGTF McMurchie Davidson expansion coefficients and their 
                   ! derivatives with respect to nuclear center A (HGBra and ECoD)
                   CALL SetBraBlocks(QP%Prim,BS,CompPrim_O=.FALSE.,DerivativeE_O=.TRUE.)
                   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                   ! Also compute derivatives of the HGTF with respect to A (Gradients=FALSE)
                   ! and do not recompute the primitive 
                   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                   CALL SetBraBlocks(QP%Prim,BS,Gradients_O=.FALSE.,CompPrim_O=.FALSE.)
                   DO K=1,3
                      dMcMrchA(1:LenD,K)=0D0
                      BraGradA(1:LenG,K)=0D0
                      IA = IndexA
                      DO LMNA=StartLA,StopLA
                         IA=IA+1
                         IB=IndexB
                         EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)                         
                         DO LMNB=StartLB,StopLB
                            IB=IB+1
                            dMcMrchA(1:LenD,K)=dMcMrchA(1:LenD,K)+ECoD%D(1:LenD,IA,IB,K)*P(IA,IB)
                            BraGradA(1:LenG,K)=BraGradA(1:LenG,K)+dHGBra%D(1:LenG,IA,IB,K)*P(IA,IB)
                         ENDDO
                      ENDDO
                   ENDDO
                   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                   ! Now compute derivatives of the HGTF with respect to distribution center P (Gradients=TRUE)
                   ! while refraining from recomputing the primitive info
                   CALL SetBraBlocks(QP%Prim,BS,Gradients_O=.TRUE.,CompPrim_O=.FALSE.)
                   EIJ(1:LenG)=0D0
                   BraGradP=0D0
                   IA = IndexA
                   DO LMNA=StartLA,StopLA
                      IA=IA+1
                      IB=IndexB
                      EllA=BS%LxDex%I(LMNA)+BS%LyDex%I(LMNA)+BS%LzDex%I(LMNA)                         
                      DO LMNB=StartLB,StopLB
                         IB=IB+1
                         EIJ(1:LenG)=EIJ(1:LenG)+HGBra%D(1:LenG,IA,IB)*P(IA,IB)*PiZ
                         DO K=1,3
                            BraGradP(1:LenG,K)=BraGradP(1:LenG,K)+P(IA,IB)*dHGBra%D(1:LenG,IA,IB,K)
                         ENDDO
                      ENDDO
                   ENDDO
                   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                   TempHerm%Coef(:)=Zero
                   TempHerm%Ell=QP%Prim%Ell
                   TempHerm%Zeta=QP%Prim%Zeta
                   DO K=1,3
                      DO I=1,LenG
                         TempHerm%Coef(I)=MAX(TempHerm%Coef(I),MAX(ABS(BraGradP(I,K)),ABS(BraGradA(I,K))))
                      ENDDO
                   ENDDO

!!!!!!!!!!!                   CALL SetSerialPAC(QP%PAC,TempHerm)                   


                   CALL SetSerialMAC(QP%MAC,TempHerm)                   
                   ! The integral estimate (ab|ab)^(1/2)
!!!!!!!                   QP%IHalf=Estimate(QP%Prim%Ell,QP%Prim%Zeta,TempHerm%Coef(1:LenG))
!!!!!!!!                   IF(QP%IHalf<TauTwo*1D-5)CYCLE
#ifdef PAC_DEBUG
                   DO L=1,LenG
                      ERRBRA(L)=TempHerm%Coef(L)
                   ENDDO
#endif
#ifdef MAC_DEBUG
                   CALL HGToSP_Bra(QP%Prim%Zeta,QP%Prim%Ell,TempHerm%Coef,SPErrorBraC,SPErrorBraS)
#endif
                   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
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
                   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                   PTmp=QP%Prim%Pw
                   DO NC=1,CS_IN%NCells 
#ifdef MAC_DEBUG
                      SPKetC=0D0
                      SPKetS=0D0
                      SPErrorKetS=Zero
                      SPErrorKetC=Zero
#endif
                      QP%Prim%Pw=PTmp+CS_IN%CellCarts%D(:,NC) 
                      ! Calculate the fractional coordinates
                      nlm=AtomToFrac(GMLoc,CS_IN%CellCarts%D(:,NC))
                      ! Store the Old Stuff
                      HGKet_old(1:LenHGTF)= HGKet(1:LenHGTF)
                      SPKetC_old(0:LenSP) = SPKetC(0:LenSP)
                      SPKetS_old(0:LenSP) = SPKetS(0:LenSP)
                      ! Walk the walk
                      NNearTmp=NNearAv
                      NNearAv=0
                      CALL JWalk2(QP,PoleRoot) 
                      NNearCount(NC)=NNearCount(NC)+NNearAv
                      NNearAv=NNearTmp+NNearAv
                      NPrim=NPrim+1
                      ! Acumulate the Lattice Forces
                      DO I=1,3
                         SPKetC_L(0:LenSP,I)=SPKetC_L(0:LenSP,I)+nlm(I)*(SPKetC(0:LenSP)-SPKetC_old(0:LenSP))
                         SPKetS_L(0:LenSP,I)=SPKetS_L(0:LenSP,I)+nlm(I)*(SPKetS(0:LenSP)-SPKetS_old(0:LenSP))
                         HGKet_L(1:LenHGTF,I)= HGKet_L(1:LenHGTF,I)+ nlm(I)*(HGKet(1:LenHGTF)-HGKet_old(1:LenHGTF))
                      ENDDO
                   ENDDO
                   QP%Prim%Pw=PTmp
                   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                   PBCNearFieldA=0D0
                   PBCNearFieldP=0D0
                   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                   !
                   DO K=1,3
                      DO LMN=1,LenG
                         PBCNearFieldP(K)=PBCNearFieldP(K)+Phase%D(LMN)*BraGradP(LMN,K)*HGKet(LMN)
                      ENDDO
                      CALL HGToSP_Direct(EllG,LenG,LenGSP,PiZ,BraGradP(1:LenG,K),SPBraC(0:LenGSP),SPBraS(0:LenGSP))                            
                      DO LM=0,LenGSP
                         PBCNearFieldP(K)=PBCNearFieldP(K)+SPBraC(LM)*SPKetC(LM)+SPBraS(LM)*SPKetS(LM)
                      ENDDO
                      DO I=1,3
                         IF(GMLoc%PBC%AutoW%I(I)==1.AND.GMLoc%PBC%AutoW%I(K)==1) THEN                             
                            DO LMN=1,LenG                               
                               LattForce(K,I)=LattForce(K,I)+2D0*Phase%D(LMN)*BraGradP(LMN,K)*HGKet_L(LMN,I)
                            ENDDO
                            DO LM=0,LenGSP
                               LattForce(K,I)=LattForce(K,I)+2D0*(SPBraC(LM)*SPKetC_L(LM,I)+SPBraS(LM)*SPKetS_L(LM,I))
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                   !
                   DO K=1,3
                      DO LMN=1,LenG
                         PBCNearFieldA(K)=PBCNearFieldA(K)+Phase%D(LMN)*BraGradA(LMN,K)*HGKet(LMN)
                      ENDDO
                      CALL HGToSP_Direct(EllG,LenG,LenGSP,PiZ,BraGradA(1:LenG,K),SPBraC(0:LenGSP),SPBraS(0:LenGSP))
                      DO LM=0,LenGSP
                         PBCNearFieldA(K)=PBCNearFieldA(K)+SPBraC(LM)*SPKetC(LM)+SPBraS(LM)*SPKetS(LM)
                      ENDDO
                   ENDDO
                   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                   dPdA=QP%Prim%ZA/QP%Prim%Zeta
                   dPdB=QP%Prim%ZB/QP%Prim%Zeta
                   !
                   dTFdPC=DerivTF(QP%Prim,EIJ,GMLoc)
                   !
                   dAdL=AtomToFrac(GMLoc,Pair%A) 
                   dBdL=AtomToFrac(GMLoc,Pair%B) 
                   dPwdL=AtomToFrac(GMLoc,QP%Prim%Pw)
                   dCdL=AtomToFrac(GMLoc,GM%PBC%CellCenter%D)
                   !
                   DO K=1,3
                      ! Tin foil terms for the lattice forces
                      QP%Prim%Ell=QP%Prim%Ell-1
                      TinFoilA(K)=dTinFoil(QP%Prim,GMLoc,dMcMrchA(:,K),PiZ)
                      QP%Prim%Ell=QP%Prim%Ell+1
                      ! Near and far field (COULD SAVE A FACTOR OF 2 HERE IN BRA TRANSLATION....)
                      AtomA(K)=PBCNearFieldA(K)+FarField(QP%Prim,BraGradA(:,K),GMLoc)  
                      AtomP(K)=PBCNearFieldP(K)+FarField(QP%Prim,BraGradP(:,K),GMLoc)  
                   ENDDO
                   !
                   IF(Pair%SameAtom) THEN
                      DO K=1,3
                         AtomForce(K)=AtomForce(K)+4D0*AtomP(K) &
                              +Four*( TinFoilA(K)+dPdA*dTFdPC(K) )*2D0 ! <<<<<<<<<<<<
                      ENDDO
                   ELSE
                      DO K=1,3
                         AtomForce(K)=AtomForce(K)+8D0*AtomA(K) &
                              +Four*( TinFoilA(K)+dPdA*dTFdPC(K) )*2D0 ! <<<<<<<<<<<<<
                      ENDDO
                   ENDIF
                   !=========================================================================================
                   LattForce=LattForce+Four*LaticeForce(GMLoc,dAdL,AtomA)
                   LattForce=LattForce+Four*LaticeForce(GMLoc,dBdL,AtomP-AtomA)
                   ! Add in the tin foil part
                   LattForce=LattForce+Four*(LaticeForce(GMLoc,dAdL,TinFoilA)-LaticeForce(GM,dBdL,TinFoilA))
                   !=========================================================================================
                   ! Notice that dPw/dL corresponds to the wrapped and implicitly 
                   ! centered distribution.
                   LattForce=LattForce+Four*LaticeForce(GMLoc,dPwdL,dTFdPC) 
                   !
                   PTmp=(QP%Prim%A*QP%Prim%ZA+QP%Prim%B*QP%Prim%ZB)/QP%Prim%Zeta
                   PTmp2=AtomToFrac(GM,PTmp)
                   DO I=1,3
                      IF(GM%PBC%AutoW%I(I)==1)THEN
                         dWdL=MODULO(PTmp2(I),1D0)-PTmp2(I)
                         DO K=1,3
                            ! (dE/dPwrap)*(dPwrap/dL) 
                            LattForce(K,I)=LattForce(K,I)+4D0*AtomP(K)*dWdL 
                         ENDDO
                      ENDIF
                   ENDDO


                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO


!!$    CALL Print_LatForce(GMLoc,Lattforce,'END END',Unit_O=6) 
!!$    WRITE(*,*)'========================================='

    !
  END SUBROUTINE TrPdJ2
  !====================================================================================================
!!$  ! Wrapped Lattice Force term involving t
!!$  !====================================================================================================
!!$  FUNCTION WLFCellCenter(C,dEdP,P,L)
!!$    REAL(DOUBLE), DIMENSION(3) :: C,dEdP,P
!!$    REAL(DOUBLE), DIMENSION(3,3) :: L, WLFCellCenter
!!$    REAL(DOUBLE) :: p1,p2,p3,c1,c2,c3,a11,a21,a22,a31,a32,a33
!!$    REAL(DOUBLE) :: dEdP1,dEdP2,dEdP3
!!$    REAL(DOUBLE) :: Modu,Mod1,X
!!$    !
!!$    Modu(X)=MODULO(X,1D0)
!!$    Mod1(X)=MODULO(X,1D0)-X
!!$    !
!!$    c1=C(1);c2=C(2);c3=C(3)
!!$    p1=P(1);p2=P(2);p3=p(3)
!!$    dEdP1=dEdP(1);dEdP2=dEdP(2);dEdP3=dEdP(3)
!!$    !
!!$    a11=L(1,1)
!!$    a21=L(2,1)
!!$    a31=L(3,1)
!!$    a22=L(2,2)
!!$    a32=L(3,2)
!!$    a33=L(3,3)
!!$    !
!!$    WLFCellCenter=Zero
!!$    !
!!$    IF(A11==0D0)THEN
!!$       RETURN
!!$    ELSE
!!$       WLFCellCenter(1,1)=dEdP1*(-(c1/a11) + modu(p1/a11))
!!$       WLFCellCenter(2,1)=-((dEdP2*(c1 + p1 - a11*mod1(p1/a11)))/a11)
!!$       WLFCellCenter(3,1)=-((dEdP3*(c1 + p1 - a11*mod1(p1/a11)))/a11)
!!$    ENDIF
!!$    !
!!$    IF(A22==0D0)THEN
!!$       RETURN
!!$    ELSE
!!$       WLFCellCenter(2,2)=(dEdP2*(a21*(c1 + p1) - a11*(c2 + p2) + a11*a22*mod1((-(a21*p1) + a11*p2)/(a11*a22))))/(a11*a22)
!!$       WLFCellCenter(3,2)=(dEdP3*(a21*(c1 + p1) - a11*(c2 + p2) + a11*a22*mod1((-(a21*p1) + a11*p2)/(a11*a22))))/(a11*a22)
!!$    ENDIF
!!$    !
!!$    IF(A33==0D0)THEN
!!$       RETURN
!!$    ELSE
!!$       WLFCellCenter(3,3)=(dEdP3*(a32*(-(a21*(c1 + p1)) + a11*(c2 + p2)) + a22*(a31*(c1 + p1) - a11*(c3 + p3)) +  &
!!$     -      a11*a22*a33*mod1((-(a22*a31*p1) + a21*a32*p1 - a11*a32*p2 + a11*a22*p3)/(a11*a22*a33))))/(a11*a22*a33)
!!$    ENDIF    
!!$    !
!!$  END FUNCTION WLFCellCenter
!!$
!!$
!!$
!!$
!!$



  FUNCTION WrapLattForce(dEdP,P,L)
    REAL(DOUBLE), DIMENSION(3) :: dEdP,P
    REAL(DOUBLE), DIMENSION(3,3) :: L, WrapLattForce
    REAL(DOUBLE) :: p1,p2,p3,a11,a21,a22,a31,a32,a33,u1,u2,u3
    REAL(DOUBLE) :: Mod1,X
    Mod1(X)=MODULO(X,1D0)-X
    !
!!$    WRITE(*,*)'1 L = ',L(:,1)
!!$    WRITE(*,*)'2 L = ',L(:,2)
!!$    WRITE(*,*)'3 L = ',L(:,3)
    !
    p1=P(1)
    p2=P(2)
    p3=P(3)
    !
    a11=L(1,1)
    a21=L(2,1)
    a31=L(3,1)
    a22=L(2,2)
    a32=L(3,2)
    a33=L(3,3)
    !
    WrapLattForce=Zero
    !
    IF(A11==0D0)THEN
       RETURN
    ELSE
       u1=Mod1(p1/a11)
       WrapLattForce(1,1)=dEdP(1)*u1
       WrapLattForce(2,1)=dEdP(2)*u1
       WrapLattForce(3,1)=dEdP(3)*u1
    ENDIF
    !
    IF(A22==0D0)THEN
       RETURN
    ELSE
       u2=Mod1((a11*p2-a21*p1)/(a11*a22))
       WrapLattForce(2,2)=dEdP(2)*u2
       WrapLattForce(3,2)=dEdP(3)*u2
    ENDIF
    !
    IF(A33==0D0)THEN
       RETURN
    ELSE
       u3=Mod1((a11*a22*p3-a22*a31*p1+a21*a32*p1-a11*a32*p2)/(a11*a22*a33))
       WrapLattForce(3,3)=dEdP(3)*u3
    ENDIF    
    !
  END FUNCTION WrapLattForce
  !====================================================================================================
  !
  !====================================================================================================
  SUBROUTINE dNukE(GMLoc,At,AtomForce,LattForce,NoWrap_O)
    TYPE(CRDS)                        :: GMLoc
    TYPE(QCPrim)                      :: QP

    REAL(DOUBLE),DIMENSION(1:3)       :: AtomForce,AtomFrc,TinF,dTFdPC, &
                                         dAdL,dCdL
    REAL(DOUBLE),DIMENSION(1:3,1:3)   :: LattForce,CellFrc,tmp
 
    REAL(DOUBLE)                      :: Tau,NukeCo,NukePole,PExtent,E000,PiZ
    REAL(DOUBLE),DIMENSION(4)         :: dBra
    REAL(DOUBLE),DIMENSION(0:SPLen)   :: SPBraC,SPBraS 
    INTEGER                           :: At,LenSP,LenHGTF,I,J,K,LM,LMN
    INTEGER                           :: NC,KI
    REAL(DOUBLE),DIMENSION(1)         :: EIJ
    REAL(DOUBLE),DIMENSION(3)         :: PTmp,nlm
    !
    REAL(DOUBLE),DIMENSION(15)        :: Vck
    REAL(DOUBLE),DIMENSION(0:SPLen,3) :: SPKetC_L,SPKetS_L
    REAL(DOUBLE),DIMENSION(1:HGLen,3) :: HGKet_L
    REAL(DOUBLE),DIMENSION(1:HGLen)   :: HGKet_old
    REAL(DOUBLE),DIMENSION(0:SPLen)   :: SPKetC_old,SPKetS_old
    LOGICAL, OPTIONAL         :: NoWrap_O
    LOGICAL                   :: NoWrap
    !---------------------------------------------------------------------------------------------
    IF(PRESENT(NoWrap_O))THEN
       NoWrap=NoWrap_O
    ELSE
       NoWrap=.FALSE.
    ENDIF
    ! Initialize (dBRA|
    !
    PiZ=(Pi/NuclearExpnt)**(ThreeHalves)    
    EIJ(1)=-GMLoc%AtNum%D(At)
    NukeCo=-GMLoc%AtNum%D(At)/PiZ

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
    CALL PWrap(GM,QP%Prim,.NOT.NoWrap)

!!!!!!!!!!!!    QP%PAC%Zeta=QP%Prim%Zeta
!!!!!!!!    QP%PAC%Wght=GMLoc%AtNum%D(At)
    QP%MAC%O(0)=GMLoc%AtNum%D(At)
    QP%MAC%Delta=Zero
!!!!    QP%IHalf=ABS(NukeCo)  !! ??????????
    ! Initialize the |KET)
    CellFrc=Zero
    AtomFrc=Zero
    HGKet(1:4)=Zero
    SPKetC(0:3)=Zero
    SPKetS(0:3)=Zero
    HGKet_L(:,:)  = Zero
    SPKetC_L(:,:) = Zero
    SPKetS_L(:,:) = Zero
    !
    PTmp=QP%Prim%Pw
    DO NC=1,CS_IN%NCells
#ifdef MAC_DEBUG
       SPKetC=0D0
       SPKetS=0D0
       SPErrorKetS=Zero
       SPErrorKetC=Zero
#endif
       ! Set atomic "primitive"
       QP%Prim%Pw=PTmp+CS_IN%CellCarts%D(:,NC)
       ! Calculate the fractional coordinates
       nlm = AtomToFrac(GMLoc,CS_IN%CellCarts%D(:,NC))
       ! Store the Old Stuff
       HGKet_old(1:4)= HGKet(1:4)
       SPKetC_old(0:3) = SPKetC(0:3)
       SPKetS_old(0:3) = SPKetS(0:3)
       ! Walk the walk
       CALL JWalk2(QP,PoleRoot,Nucular_O=.TRUE.)
       ! Acumulate the Lattice Forces
       DO I=1,3
          HGKet_L(1:4,I)= HGKet_L(1:4,I)+ nlm(I)*(HGKet(1:4)- HGKet_old(1:4))
          SPKetC_L(0:3,I) = SPKetC_L(0:3,I) + nlm(I)*(SPKetC(0:3) - SPKetC_old(0:3))
          SPKetS_L(0:3,I) = SPKetS_L(0:3,I) + nlm(I)*(SPKetS(0:3) - SPKetS_old(0:3))
       ENDDO
    ENDDO

    ! Reset the Atomic Coordinates
    QP%Prim%Pw=PTmp
    ! Init bra xforms
    DO K=1,3
       DO LMN=1,4
          AtomFrc(K)=AtomFrc(K)+Phase%D(LMN)*dHGBra%D(LMN,1,1,K)*HGKet(LMN)
       ENDDO
       CALL HGToSP(QP%Prim%Zeta,1,dHGBra%D(1:4,1,1,K),SPBraC(0:2),SPBraS(0:2))
       DO LM=0,2
          AtomFrc(K)=AtomFrc(K)+SPBraC(LM)*SPKetC(LM)+SPBraS(LM)*SPKetS(LM)
       ENDDO
    ENDDO

    !      Inner Sum Nuc Box Force
    DO K=1,3
       CALL HGToSP(QP%Prim%Zeta,1,dHGBra%D(:,1,1,K),SPBraC,SPBraS)
       DO I=1,3
          KI = K + 3*(I+1)
          IF(GMLoc%PBC%AutoW%I(I)==1 .AND. GMLoc%PBC%AutoW%I(K)==1) THEN
             DO LMN=1,4
                CellFrc(K,I)=CellFrc(K,I)+Phase%D(LMN)*dHGBra%D(LMN,1,1,K)*HGKet_L(LMN,I)
             ENDDO
             DO LM=0,2
                CellFrc(K,I)=CellFrc(K,I)+SPBraC(LM)*SPKetC_L(LM,I)+SPBraS(LM)*SPKetS_L(LM,I)
             ENDDO
          ENDIF
       ENDDO
    ENDDO
    DO K=1,3
       TinF(K)=TinFoil(QP%Prim,GMLoc,dHGBra%D(1:4,1,1,K),PiZ)  
       AtomFrc(K)=AtomFrc(K)+FarField(QP%Prim,dHGBra%D(1:4,1,1,K),GMLoc)  
    ENDDO
    !
    dAdL=AtomToFrac(GMLoc,GMLoc%Carts%D(:,At))
    dCdL=AtomToFrac(GMLoc,GMLoc%PBC%CellCenter%D)
    dTFdPC=DerivTF(QP%Prim,EIJ,GMLoc)
    ! 
    AtomForce=AtomForce+Two*dTFdPC !<<<<<<<<<<<
    AtomForce=AtomForce+Two*AtomFrc
    !
    tmp=Two*LaticeForce(GM,dAdL,AtomFrc)    
!!$    write(*,*)' dAdL = ',dAdL
!!$    write(*,*)' AtomFrc = ',AtomFrc
!!$    STOP

    LattForce=LattForce+CellFrc+tmp
    tmp=Two*LaticeForce(GM,dAdL-dCdL,dTFdPC)
    LattForce=LattForce+tmp
    !
  END SUBROUTINE dNukE

END MODULE BlokTrPdJ


