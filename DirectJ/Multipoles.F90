!==============================================================================
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe and  C. J. Tymczak
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
!    MULTIPOLES
!==============================================================================
MODULE Multipoles
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  USE Moments
  USE MMML
  IMPLICIT NONE
#ifdef PERIODIC
!---------------------------------------------------------------------------
! Global variables
!---------------------------------------------------------------------------
  INTEGER               :: Dimen
  INTEGER               :: BigL
  INTEGER               :: BigLM
  TYPE(SpMoments)       :: MM_Rho
  TYPE(SpMoments)       :: MM_Tensor
  TYPE(SpMoments)       :: MM_TenRho
!
  REAL(DOUBLE)          :: CellCenterX
  REAL(DOUBLE)          :: CellCenterY
  REAL(DOUBLE)          :: CellCenterZ
!
  REAL(DOUBLE)          :: Volume
  REAL(DOUBLE)          :: QQFac
  REAL(DOUBLE)          :: DDFac
  TYPE(QMoments)        :: Q_Rho
!
CONTAINS
!========================================================================================
! Set up the nessesary information to calculate the multipoles
!========================================================================================
  SUBROUTINE MMSetup(MaxL,GM,Rho)
    TYPE(CRDS)                                :: GM
    TYPE(HGRho)                               :: Rho
    INTEGER                                   :: MaxL,MaxLM,I,J,K
    INTEGER                                   :: IMin,JMin,KMin
!
    MaxLM = LSP(MaxL)
    BigL  = MaxL
    BigLM = MaxLM
!
!   Calculate the Factorials and Allocate Memory
!
    CALL FactorialSetUpF90(BigL)
!
!   Calculate the Dimension
!
    Dimen = 0
    DO I = 1,3
       IF(GM%AutoW(I)) Dimen = Dimen+1
    ENDDO
!
!   Calculate the Box Volume
!
    Volume = One
    DO I=1,3
       IF(GM%AutoW(I)) THEN
          Volume = Volume*GM%BoxShape%D(I,I)
       ENDIF
    ENDDO
!
!   Calculate the Center of the Cell
!
    CellCenterX = Zero
    CellCenterY = Zero
    CellCenterZ = Zero
    DO I = 1,3
       IF(GM%AutoW(I)) THEN
          CellCenterX = CellCenterX+Half*GM%BoxShape%D(I,1)
          CellCenterY = CellCenterY+Half*GM%BoxShape%D(I,2)
          CellCenterZ = CellCenterZ+Half*GM%BoxShape%D(I,3)
       ENDIF
    ENDDO
!
!   Calculate the Size of the Box Needed  for the Direct J
!
    CALL BoxBounds(GM,Rho,IMin,JMin,KMin)
!
!   Create the Inner Box Cell Set for DirectJ
!
    CALL New_CellSet_Cube(CSMM1,GM%BoxShape%D,(/IMin,JMin,KMin/))
! 
!   Intitialize The MM Tensor
!  
    CALL New_SpMoments(MM_Tensor,BigL)
    MM_Tensor%CMMat%D = Zero
    MM_Tensor%SMMat%D = Zero
    MM_Tensor%CenterX = Zero
    MM_Tensor%CenterY = Zero
    MM_Tensor%CenterZ = Zero
!
    CALL New_MM_Tensor(GM,IMin,JMin,KMin)
!
!   Set Up the Factors For the the Cartesian Multipole Correction to the Matrix Elements
!
!      Phi_EW(R) = Phi_MM(R) - (4*Pi/3*V)*(R-R0).d + (4*Pi/3*V)*Tr(Q)  (3-Dimension
!
    IF(Dimen==0) THEN
       QQFac = Zero
       DDFac = Zero      
    ELSEIF(Dimen==1) THEN
       QQFac = Zero
       DDFac = Zero
    ELSEIF(Dimen==2) THEN
    ELSEIF(Dimen==3) THEN
       QQFac = Four*Pi/(Three*Volume)
       DDFac = Four*Pi/(Three*Volume)
    ENDIF
  END SUBROUTINE MMSetup
!========================================================================================
! Calculate the MM_Tensor
!========================================================================================
  SUBROUTINE New_MM_Tensor(GM,IMin,JMin,KMin)
    TYPE(CRDS)                         :: GM 
    INTEGER                            :: IMin,JMin,KMin
    INTEGER                            :: I,J,L,M,LM,NC
    INTEGER                            :: IMax,JMax,KMax,IJKMax,LSwitch
    REAL(DOUBLE)                       :: PQx,PQy,PQz
    REAL(DOUBLE)                       :: CFac,SFac,Beta,Beta0,Rad,R2B,ExpFac
    REAL(DOUBLE),DIMENSION(0:BigLM)    :: Cpq,Spq
!
    LSwitch = 8
    Beta0 = Pi/(Volume**(1.D0/3.D0))
!
    IF(Dimen==1) THEN
    ELSEIF(Dimen==2) THEN
    ELSEIF(Dimen==3) THEN
!
!      Sum the Real Space
!
       IJKMax = 32
       IMax = IJKMax+IMin
       JMax = IJKMax+KMin
       KMax = IJKMax+JMin
       CALL New_CellSet_Cube(CSMM2,GM%BoxShape%D,(/IMax,JMax,KMax/),(/IMin+1,IMin+1,IMin+1/))
       DO NC = 1,CSMM2%NCells
          PQx = CSMM2%CellCarts%D(1,NC)
          PQy = CSMM2%CellCarts%D(2,NC)
          PQz = CSMM2%CellCarts%D(3,NC)
          Rad  = SQRT(PQx*PQx+PQy*PQy+PQz*PQz)
          CALL IrRegularF90(BigL,BigLM,PQx,PQy,PQz,Cpq,Spq)
          DO L = 3,BigL
             IF(L .LE. LSwitch) THEN
                Beta = Beta0
                R2B  = Beta*Beta*Rad*Rad
                CFac = GScript(L,R2B)
             ELSE
                CFac = One
             ENDIF
             DO M = 0,L
                LM = LTD(L)+M
                MM_Tensor%CMMat%D(LM)=MM_Tensor%CMMat%D(LM)+Cpq(LM)*CFac
                MM_Tensor%SMMat%D(LM)=MM_Tensor%SMMat%D(LM)+Spq(LM)*CFac
             ENDDO
          ENDDO
       ENDDO
       CALL Delete_CellSet(CSMM2)
!
!      Sum the Reciprical Space 
!
       ExpFac = (Pi*Pi)/(Beta0*Beta0)
       CALL New_CellSet_Cube(CSMM2,GM%InvBoxSh%D,(/IMax,JMax,KMax/),(/1,1,1/))
       DO NC = 1,CSMM2%NCells
          PQx = CSMM2%CellCarts%D(1,NC)
          PQy = CSMM2%CellCarts%D(2,NC)
          PQz = CSMM2%CellCarts%D(3,NC)
          Rad = SQRT(PQx*PQx+PQy*PQy+PQz*PQz)
          CALL IrRegularF90(BigL,BigLM,PQx,PQy,PQz,Cpq,Spq)
          DO L = 3,BigL
             IF(L .LE. LSwitch) THEN             
                CFac = FT_FScriptC(L,ExpFac,Rad)/Volume
                SFac = FT_FScriptS(L,ExpFac,Rad)/Volume
             ELSE
                CFac = Zero
                SFac = Zero
             ENDIF
             DO M = 0,L
                LM = LTD(L)+M
                MM_Tensor%CMMat%D(LM)=MM_Tensor%CMMat%D(LM)+Cpq(LM)*CFac-Spq(LM)*SFac
                MM_Tensor%SMMat%D(LM)=MM_Tensor%SMMat%D(LM)+Spq(LM)*CFac+Cpq(LM)*SFac
             ENDDO
          ENDDO
       ENDDO
       CALL Delete_CellSet(CSMM2)
!
!      Substract the inner boxes
!
       CALL New_CellSet_Cube(CSMM2,GM%BoxShape%D,(/IMin,JMin,KMin/),(/1,1,1/))  
       DO NC = 1,CSMM2%NCells
          PQx=CSMM2%CellCarts%D(1,NC)
          PQy=CSMM2%CellCarts%D(2,NC)
          PQz=CSMM2%CellCarts%D(3,NC)
          Rad  = SQRT(PQx*PQx+PQy*PQy+PQz*PQz)
          CALL IrRegularF90(BigL,BigLM,PQx,PQy,PQz,Cpq,Spq)
          DO L = 3,BigL
             IF(L .LE. LSwitch) THEN
                Beta = Beta0
                R2B  = Beta*Beta*Rad*Rad
                CFac = FScript(L,R2B)
             ELSE
                CFac = Zero
             ENDIF
             DO M = 0,L
                LM = LTD(L)+M
                MM_Tensor%CMMat%D(LM)=MM_Tensor%CMMat%D(LM)-Cpq(LM)*CFac
                MM_Tensor%SMMat%D(LM)=MM_Tensor%SMMat%D(LM)-Spq(LM)*CFac
             ENDDO
          ENDDO
       ENDDO
       CALL Delete_CellSet(CSMM2)
!
!      Filter out the Zero Elements
!
       DO L = 0,BigL
          DO M = 0,L
             LM = LTD(L)+M
             IF(ABS(MM_Tensor%CMMat%D(LM)) .LT. 1.D-15) MM_Tensor%CMMat%D(LM)=Zero
             IF(ABS(MM_Tensor%SMMat%D(LM)) .LT. 1.D-15) MM_Tensor%SMMat%D(LM)=Zero
          ENDDO
       ENDDO
!
    ENDIF
!
!    CALL PPrint_SpMoments(MM_Tensor,'MM_Tensor','MM_Tensor',6)
!
  END SUBROUTINE New_MM_Tensor
!========================================================================================
! Calculate the Multipole Tensor of the local density
!========================================================================================
  SUBROUTINE CalMMRho(Rho)
    TYPE(HGRho)                     :: Rho
    INTEGER                         :: LP,LQ,LPQ,LenP,LenQ,LenPQ,LKet
    INTEGER                         :: zq,iq,iadd,jadd,NQ,OffQ,OffR,LM,I
    REAL(DOUBLE)                    :: Zeta,Px,Py,Pz,PQx,PQy,PQz
    REAL(DOUBLE),DIMENSION(0:BigLM) :: Cq,Sq,Cpq,Spq
!
!   Initialize
!
    CALL New_SpMoments(MM_Rho,BigL)
    MM_Rho%CMMat%D = Zero
    MM_Rho%SMMat%D = Zero
    MM_Rho%CenterX = CellCenterX
    MM_Rho%CenterY = CellCenterY
    MM_Rho%CenterZ = CellCenterZ
!
!   Calculate the Moments of the Density, Centered in the Box
!
    DO zq=1,Rho%NExpt
       NQ    = Rho%NQ%I(zq)
       Zeta  = Rho%Expt%D(zq)
       OffQ  = Rho%OffQ%I(zq)
       OffR  = Rho%OffR%I(zq)
!
       LQ    = Rho%Lndx%I(zq)
       LP    = BigL-LQ
       LPQ   = LP+LQ
       LenP  = LSP(LP)  
       LenQ  = LSP(LQ)  
       LenPQ = LSP(LPQ) 
       LKet  = LHGTF(LQ)
!
       IF(NQ /= 0) THEN
          DO iq = 1,NQ
             iadd = Rho%OffQ%I(zq)+iq
             jadd = Rho%OffR%I(zq)+(iq-1)*LKet+1
!
             Px = Rho%Qx%D(iadd) - CellCenterX
             Py = Rho%Qy%D(iadd) - CellCenterY
             Pz = Rho%Qz%D(iadd) - CellCenterZ
!
             CALL HGTFToSP(LQ,Zeta,Cq,Sq,Rho%Co%D(jadd:))
             IF(NoTranslate(Px,Py,Pz))THEN
                DO LM = 0,LenQ
                   MM_Rho%CMMat%D(LM) =  MM_Rho%CMMat%D(LM) + Cq(LM)
                   MM_Rho%SMMat%D(LM) =  MM_Rho%SMMat%D(LM) + Sq(LM)
                ENDDO
             ELSE
                CALL XLate(BigL,LP,LQ,LPQ,LenP,LenQ,LenPQ,Px,Py,Pz, &
                           MM_Rho%CMMat%D,MM_Rho%SMMat%D,           &
                           Cpq,Spq,                                 &
                           Cq,Sq,                                   &
                           LegendreP%D,Cosine%D,Sine%D,RToTh%D,     &
                           FactOlm0%D,FactOlm1%D,FactOlm2%D)
             ENDIF
          ENDDO
       ENDIF
    ENDDO
!
!    CALL PPrint_SpMoments(MM_Rho,'MM_Rho','MM_Rho',6)
!
!   Zero the Monopole piece to maintain charge neutrality
!
    MM_Rho%CMMat%D(0) = Zero
    MM_Rho%SMMat%D(0) = Zero
!
!   Contract the Density Multipoles with the Multipole Tensor
!
    CALL New_SpMoments(MM_TenRho,BigL)
    MM_TenRho%CMMat%D = Zero
    MM_TenRho%SMMat%D = Zero
    MM_TenRho%CenterX = CellCenterX
    MM_TenRho%CenterY = CellCenterY
    MM_TenRho%CenterZ = CellCenterZ
!
!   Calculate MM_TenRho(l,m) = Sum_l'm' M(l+l',m+m') * rho(l'.m')
!
    LQ     = BigL/2
    LP     = BigL/2
    LPQ    = LP+LQ
    LenP   = LSP(LP)  
    LenQ   = LSP(LQ)  
    LenPQ  = LSP(LPQ)
!
    CALL CTraX(BigL,LP,LQ,LPQ,LenP,LenQ,LenPQ,PQx,PQy,PQz, &
               MM_TenRho%CMMat%D,MM_TenRho%SMMat%D,        &
               MM_Tensor%CMMat%D,MM_Tensor%SMMat%D,        &
               MM_Rho%CMMat%D,MM_Rho%SMMat%D,              &
               LegendreP%D,Cosine%D,Sine%D,RToTh%D,        &
               FactMlm0%D,FactMlm1%D,FactMlm2%D,.FALSE.)
!
!    CALL PPrint_SpMoments(MM_TenRho,'MM_TenRho','MM_TenRho',6)
!
  END SUBROUTINE CalMMRho
!========================================================================================
! Calculate the Cartesian Multipole Tensor of the local density
!========================================================================================
  SUBROUTINE CalQRho(Rho)
    TYPE(HGRho)                     :: Rho
    INTEGER                         :: zq,iq,iadd,jadd,NQ,OffQ,OffR,LQ,LP,LenP,LenQ
    REAL(DOUBLE)                    :: RX,RY,RZ,Expt
!
    LQ   = 2
    LenQ = LHGTF(LP)
    CALL New_QMoments(Q_Rho,LQ)
    Q_Rho%QMat%D  = Zero
    Q_Rho%CenterX = CellCenterX
    Q_Rho%CenterY = CellCenterY
    Q_Rho%CenterZ = CellCenterZ
!
    DO zq=1,Rho%NExpt
       NQ    = Rho%NQ%I(zq)
       Expt  = Rho%Expt%D(zq)
       OffQ  = Rho%OffQ%I(zq)
       OffR  = Rho%OffR%I(zq)
       LP    = Rho%Lndx%I(zq) 
       LenP  = LHGTF(LP)
       IF(NQ > 0) THEN
          DO iq = 1,NQ
             iadd = Rho%OffQ%I(zq)+iq
             jadd = Rho%OffR%I(zq)+(iq-1)*LenP+1
             RX = Rho%Qx%D(iadd)-CellCenterX
             RY = Rho%Qy%D(iadd)-CellCenterY
             RZ = Rho%Qz%D(iadd)-CellCenterZ 
             CALL HGTFToCarts(LQ,LP,Expt,RX,RY,RZ,Q_Rho%QMat%D,Rho%Co%D(jadd:)) 
          ENDDO
       ENDIF
    ENDDO
!
    CALL PPrint_QMoments(Q_Rho,'Q_Rho','Q_Rho',6)
!
  END SUBROUTINE CalQRho
!========================================================================================
! Contract the Primative with the Density
!========================================================================================
  FUNCTION CTraxBraKet(LBra,Expt,Qx,Qy,Qz,Coef)
!
    INTEGER                             :: L,M,LM,LBra,LQ,LP,LPQ,LenQ,LenP,LenPQ,LenBra
    REAL(DOUBLE)                        :: Expt,Qx,Qy,Qz,DiPole,QuadPole
    REAL(DOUBLE)                        :: PQx,PQy,PQz,CTraxBraKet
    REAL(DOUBLE),DIMENSION(:)           :: Coef
    REAL(DOUBLE),DIMENSION(0:BigLM)     :: Cpq,Spq,Cq,Sq,Cp,Sp
!
    CTraxBraKet = Zero
!
    PQx   = Qx-CellCenterX
    PQy   = Qy-CellCenterY
    PQz   = Qz-CellCenterZ
!
    LQ    = LBra
    LP    = BigL-LQ
    LPQ   = LQ+LP
    LenQ  = LSP(LQ)
    LenP  = LSP(LP)
    LenPQ = LSP(LPQ)
    LenBra= LHGTF(LQ)
!
    Cp(:) = Zero
    Sp(:) = Zero
    Cq(:) = Zero
    Sq(:) = Zero
!
    CALL HGTFToSP(LQ,Expt,Cq,Sq,Coef)
    IF(.NOT. NoTranslate(PQx,PQy,PQz)) THEN
       CALL XLate(BigL,LP,LQ,LPQ,LenP,LenQ,LenPQ,PQx,PQy,PQz, &
                  Cp,Sp,                                   &
                  Cpq,Spq,                                 &
                  Cq,Sq,                                   &
                  LegendreP%D,Cosine%D,Sine%D,RToTh%D,     &
                  FactOlm0%D,FactOlm1%D,FactOlm2%D)
    ENDIF
!
    DO L = 0,BigL
       DO M = 0,L
          LM = LTD(L)+M
          CTraxBraKet = CTraxBraKet + (Cp(LM)*MM_TenRho%CMMat%D(LM) &
               + Sp(LM)*MM_TenRho%SMMat%D(LM))
       ENDDO
    ENDDO
!
  END FUNCTION CTraxBraKet
!========================================================================================
! Contract the Primative with the Density
!========================================================================================
  FUNCTION QTraxBraKet(LBra,Expt,Qx,Qy,Qz,Coef)
!
    INTEGER                             :: LQ,LBra
    REAL(DOUBLE)                        :: Expt,Qx,Qy,Qz,DiPole,QuadPole
    REAL(DOUBLE)                        :: RX,RY,RZ,QTraxBraKet
    REAL(DOUBLE),DIMENSION(:)           :: Coef
    REAL(DOUBLE),DIMENSION(1:10)        :: QTemp
!
    RX   = Qx-CellCenterX
    RY   = Qy-CellCenterY
    RZ   = Qz-CellCenterZ
!
    LQ    = 2
    QTemp = Zero
    CALL HGTFToCarts(LQ,LBra,Expt,RX,RY,RZ,QTemp,Coef)
!
    DiPole   = DDFac*(QTemp(LMNDex(1,0,0))*Q_Rho%QMat%D(LMNDex(1,0,0)) + &
                      QTemp(LMNDex(0,1,0))*Q_Rho%QMat%D(LMNDex(0,1,0)) + &
                      QTemp(LMNDex(0,0,1))*Q_Rho%QMat%D(LMNDex(0,0,1)))
!
    QuadPole = QQFac*QTemp(LMNDex(0,0,0))*(Q_Rho%QMat%D(LMNDex(2,0,0)) + &
                                           Q_Rho%QMat%D(LMNDex(0,2,0)) + &
                                           Q_Rho%QMat%D(LMNDex(0,0,2)))
!
    QTraxBraKet = -DiPole+QuadPole
!
  END FUNCTION QTraxBraKet
!========================================================================================
! Calculate the Box Bounds Needed for the Direct Sum
!========================================================================================
  SUBROUTINE BoxBounds(GM,Rho,I,J,K) 
    TYPE(CRDS)                      :: GM
    TYPE(HGRho)                     :: Rho
    INTEGER                         :: I,J,K,iadd,iq,zq,NQ
    REAL(DOUBLE)                    :: XMin,YMin,ZMin,X,Y,Z
!
    XMin = Zero
    YMin = Zero
    ZMin = Zero
    DO zq=1,Rho%NExpt
       NQ  = Rho%NQ%I(zq)
       DO iq = 1,NQ
          iadd = Rho%OffQ%I(zq)+iq
          X = Rho%Qx%D(iadd)-CellCenterX
          Y = Rho%Qy%D(iadd)-CellCenterY
          Z = Rho%Qz%D(iadd)-CellCenterZ
          IF(X > XMin) XMin = X
          IF(Y > YMin) YMin = Y
          IF(Z > ZMin) ZMin = Z
       ENDDO
    ENDDO
!
    I  = XMin/GM%BoxShape%D(1,1)+1
    J  = YMin/GM%BoxShape%D(2,2)+1
    K  = ZMin/GM%BoxShape%D(3,3)+1
!
    IF(.NOT. GM%AutoW(1)) I = 0
!    IF(.NOT. GM%AutoW(2)) J = 0
!    IF(.NOT. GM%AutoW(3)) K = 0
!
  END SUBROUTINE BoxBounds
!========================================================================================
! Calculate the Box Bounds Needed for the Direct Sum
!========================================================================================
  FUNCTION NoTranslate(X,Y,Z) 
    REAL(DOUBLE)           :: X,Y,Z
    REAL(DOUBLE),PARAMETER :: TOL=1.0D-15
    LOGICAL                :: NoTranslate
!
    NoTranslate = (ABS(X).LT.TOL) .AND. (ABS(Y).LT.TOL) .AND. (ABS(Z).LT.TOL)
!
  END FUNCTION NoTranslate
!
#endif
!
END MODULE Multipoles










