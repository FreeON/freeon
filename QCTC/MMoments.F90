MODULE MMoments
  USE DerivedTypes
  USE GlobalScalars
  USE PrettyPrint
  USE McMurchie
  USE AtomPairs
  IMPLICIT NONE
!-------------------------------------------
! TYPE MMom, Holds Info on Multipole Moments
!-------------------------------------------
  TYPE Moments
     INTEGER            :: Alloc
     INTEGER            :: MaxL
     INTEGER            :: MaxLM
     REAL(DOUBLE)       :: CenterX
     REAL(DOUBLE)       :: CenterY
     REAL(DOUBLE)       :: CenterZ    
     TYPE(DBL_VECT)     :: CMMat
     TYPE(DBL_VECT)     :: SMMat
  ENDTYPE Moments
!---------------------------------------------------------------------------
! Global variables
!---------------------------------------------------------------------------
  INTEGER               :: BigL
  INTEGER               :: BigLM
  TYPE(Moments)         :: MM_Rho
  TYPE(Moments)         :: MM_Tensor
  TYPE(Moments)         :: MM_TenRho
  TYPE(DBL_VECT)        :: Factorial
  TYPE(DBL_VECT)        :: FactOlm0
  TYPE(DBL_VECT)        :: FactOlm1
  TYPE(DBL_RNK2)        :: FactOlm2
  TYPE(DBL_VECT)        :: FactMlm0
  TYPE(DBL_VECT)        :: FactMlm1
  TYPE(DBL_RNK2)        :: FactMlm2
!
  TYPE(DBL_VECT)        :: Cosine
  TYPE(DBL_VECT)        :: Sine
  TYPE(DBL_VECT)        :: RToTh
  TYPE(DBL_RNK2)        :: LegendreP
!
  REAL(DOUBLE)          :: CellCenterX
  REAL(DOUBLE)          :: CellCenterY
  REAL(DOUBLE)          :: CellCenterZ
!
CONTAINS
!========================================================================================
! ALLOCATE  new moments
!========================================================================================
  SUBROUTINE New_Moments(A,N_O)
    TYPE(Moments)                   :: A
    INTEGER,OPTIONAL                :: N_O
!
   IF(AllocQ(A%Alloc)) THEN
       CALL Delete_Moments(A)       
       IF(PRESENT(N_O)) THEN
          A%MaxL  = N_O
          A%MaxLM = LSP(N_O)
       ELSE
          A%MaxL  = 0        
          A%MaxLM = 0
       ENDIF
       CALL New(A%CMMat,A%MaxLM,0)  
       CALL New(A%SMMat,A%MaxLM,0)     
    ELSE
       A%Alloc=ALLOCATED_TRUE
       IF(PRESENT(N_O)) THEN
          A%MaxL  = N_O
          A%MaxLM = LSP(N_O)
       ELSE
          A%MaxL  = 0
          A%MaxLM = 0          
       ENDIF
       CALL New(A%CMMat,A%MaxLM,0)  
       CALL New(A%SMMat,A%MaxLM,0) 
    ENDIF
!
  END SUBROUTINE New_Moments
!========================================================================================
! Delete  moments
!========================================================================================
  SUBROUTINE Delete_Moments(A)
    TYPE(Moments)                   :: A
!
    IF(AllocQ(A%Alloc)) THEN
       A%Alloc=ALLOCATED_FALSE
       A%MaxL  = 0
       A%MaxLM = 0
       CALL Delete(A%CMMat)
       CALL Delete(A%SMMat)
    ENDIF
!
  END SUBROUTINE Delete_Moments
!========================================================================================
! Print  moments
!========================================================================================
  SUBROUTINE PPrint_Moments(A,Name,FileName_O,Unit_O)
    TYPE(Moments)                      :: A
    CHARACTER(LEN=*)                   :: Name   
    CHARACTER(LEN=*),OPTIONAL          :: FileName_O
    INTEGER,OPTIONAL                   :: Unit_O
    INTEGER                            :: OutU
    INTEGER                            :: L,M,LM
!
    IF(.NOT. AllocQ(A%Alloc)) THEN
       CALL Halt(' Moments not allocated in PPrint_Moments')
    ENDIF
    IF(PRESENT(Unit_O)) THEN
       OutU=Unit_O
    ELSE
       OutU=Out
    ENDIF
    IF(PRESENT(FileName_O) .AND. OutU /= 6) THEN
       CALL OpenASCII(FileName_O,OutU)
    ELSEIF(OutU /= 6) THEN
       CALL OpenASCII(OutFile,OutU)
    ENDIF
!
    WRITE(OutU,5)
    WRITE(OutU,10) Name
    WRITE(OutU,11) A%MaxL,A%MaxLM
    WRITE(OutU,12) A%CenterX,A%CenterY,A%CenterZ
    WRITE(OutU,6)
!
    DO L = 0,A%MaxL
       WRITE(OutU,20) L
       DO M = 0,L
          LM = LTD(L)+M
          WRITE(OutU,21) L,M,A%CMMat%D(LM),L,M,A%SMMat%D(LM)
       ENDDO
    ENDDO
    WRITE(OutU,5)
    RETURN
!
5   FORMAT(1x,'=========================================================================')
6   FORMAT(1x,'=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=')
10  FORMAT(1x,A)
11  FORMAT(1x,'MaxL = ',I2,' MaxLM = ',I3)
12  FORMAT(1x,'CellCenter = (',D16.8,',',D16.8,',',D16.8,')')
20  FORMAT(1x,'L = ',I2)
21  FORMAT(2x,'CMat(',I2,',',I2,') = ',D16.8,2x,'SMat(',I2,',',I2,') = ',D16.8)
!
  END SUBROUTINE PPrint_Moments
!========================================================================================
! Set up the nessesary information to calculate the multipoles
!========================================================================================
  SUBROUTINE MMSetup(MaxL,GM,Rho)
    TYPE(CRDS)                          :: GM
    TYPE(HGRho)                         :: Rho
    INTEGER                             :: MaxL,MaxLM,I,J,K,NC,L,M,LM
    INTEGER                             :: IMin,JMin,KMin,IJKMax
    REAL(DOUBLE)                        :: PQx,PQy,PQz,RMAG,VAL
    REAL(DOUBLE),DIMENSION(0:MaxL*(MaxL+3)/2) :: Cpq,Spq
!
    MaxLM = LSP(MaxL)
    BigL  = MaxL
    BigLM = MaxLM
!
!   Calculate the Factorials and Allocate Memory
!
    CALL New(Factorial,2*MaxL,0)
!
    CALL New(FactOlm0 ,MaxL  ,0)
    CALL New(FactOlm1 ,MaxL  ,0)
    CALL New(FactMlm0 ,MaxL  ,0)
    CALL New(FactMlm1 ,MaxL  ,0)
!
    CALL New(FactOlm2 ,(/MaxL,MaxL/),(/0,0/))
    CALL New(FactMlm2 ,(/MaxL,MaxL/),(/0,0/))
    CALL MultipoleSetUp(MaxL,Factorial%D,FactOlm0%D,FactOlm1%D,FactOlm2%D, &
                                         FactMlm0%D,FactMlm1%D,FactMlm2%D)
!
    CALL New(Cosine   ,MaxL,0)
    CALL New(Sine     ,MaxL,0)
    CALL New(RToTh    ,MaxL,0)
    CALL New(LegendreP,(/MaxL,MaxL/),(/0,0/))
! 
!   Intitialize The MM Tensor
!  
    CALL New_Moments(MM_Tensor,MaxL)
    MM_Tensor%CMMat%D = Zero
    MM_Tensor%SMMat%D = Zero
    MM_Tensor%CenterX = Zero
    MM_Tensor%CenterY = Zero
    MM_Tensor%CenterZ = Zero
!
#ifdef PERIODIC
!
!   Calculate the Center of the Cell
!
    DO I = 1,3
       CellCenterX = CellCenterX+Half*GM%BoxShape%D(I,1)
       CellCenterY = CellCenterY+Half*GM%BoxShape%D(I,2)
       CellCenterZ = CellCenterZ+Half*GM%BoxShape%D(I,3)
    ENDDO
!
!   Calculate the Size of the Box Needed  for the Direct J
!
    CALL BoxBounds(GM,Rho,IMin,JMin,KMin)
!
!   Create the Inner Box Cell Set for DirectJ
!
    CALL New_CellSet_Cube(CSMM1,GM,(/IMin,JMin,KMin/))

!
    IF(BoxIsCube(GM)) THEN
       OPEN(UNIT=99,FILE='/earth/save/tymczak/MONDO/TwoE/Multipoles_Cube.dat',STATUS='OLD')
       DO I = 1,100
          READ(99,*) L,M,VAL
          IF(L .LE. MaxL) THEN
             LM = LTD(L)+M
             MM_Tensor%CMMat%D(LM) = VAL
          ENDIF
       ENDDO
       CLOSE(99)
!
!      Substract the inner boxes and normailze to box size
!
       CALL New_CellSet_Cube(CSMM2,GM,(/IMin,JMin,KMin/),(/1,1,1/)) 
       VAL = One
       DO L = 0,MaxL
          VAL  = VAL*GM%BoxShape%D(1,1)
          DO M = 0,L
             LM = LTD(L)+M
             MM_Tensor%CMMat%D(LM)=MM_Tensor%CMMat%D(LM)/VAL
          ENDDO
       ENDDO
       DO NC = 1,CSMM2%NCells
          PQx=CSMM2%CellCarts%D(1,NC)
          PQy=CSMM2%CellCarts%D(2,NC)
          PQz=CSMM2%CellCarts%D(3,NC)
          CALL IrRegular(MaxL,MaxL,MaxLM,PQx,PQy,PQz,Cpq,Spq, &
                         LegendreP%D,Cosine%D,Sine%D,RToTh%D, &
                         FactMlm0%D,FactMlm1%D,FactMlm2%D)
          DO L = 0,MaxL
             DO M = 0,L
                LM = LTD(L)+M
                IF(MM_Tensor%CMMat%D(LM) .NE. Zero) THEN
                   MM_Tensor%CMMat%D(LM)=MM_Tensor%CMMat%D(LM)-Cpq(LM)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       CALL Delete_CellSet(CSMM2)
    ELSE
!
!      For Now do a Direct Summation
!
       IJKMax = 64
       CALL New_CellSet_Cube(CSMM2,GM,(/IJKMax,IJKMax,IJKMax/),(/IMin,JMin,KMin/))
!
       DO NC = 1,CSMM2%NCells
          PQx=CSMM2%CellCarts%D(1,NC)
          PQy=CSMM2%CellCarts%D(2,NC)
          PQz=CSMM2%CellCarts%D(3,NC)
          CALL IrRegular(MaxL,MaxL,MaxLM,PQx,PQy,PQz,Cpq,Spq, &
                         LegendreP%D,Cosine%D,Sine%D,RToTh%D, &
                         FactMlm0%D,FactMlm1%D,FactMlm2%D)
          DO J = 1,MaxLM
             MM_Tensor%CMMat%D(J)=MM_Tensor%CMMat%D(J)+Cpq(J)
             MM_Tensor%SMMat%D(J)=MM_Tensor%SMMat%D(J)+Spq(J)
          ENDDO
       ENDDO
       CALL Delete_CellSet(CSMM2)
    ENDIF
#endif
  END SUBROUTINE MMSetup
!========================================================================================
!
!========================================================================================
  SUBROUTINE CalMMRho(Rho)
    TYPE(HGRho)                     :: Rho
!
    INTEGER                         :: LP,LQ,LPQ,LenP,LenQ,LenPQ,LKet
    INTEGER                         :: zq,iq,iadd,jadd,NQ,OffQ,OffR,LM,I
    INTEGER                         :: L1,L2,L12,M1,M2,M12,LM1,LM2,LM12
    REAL(DOUBLE)                    :: Zeta,Px,Py,Pz
    REAL(DOUBLE)                    :: CM1,CM2,CM12,SM1,SM2,SM12,Csum,Ssum
    REAL(DOUBLE),DIMENSION(0:BigLM) :: Cq,Sq,Cpq,Spq
!
!   Initialize
!
    CALL New_Moments(MM_Rho,BigL)
    MM_Rho%CMMat%D = Zero
    MM_Rho%SMMat%D = Zero
    MM_Rho%CenterX = CellCenterX
    MM_Rho%CenterY = CellCenterY
    MM_Rho%CenterZ = CellCenterZ
!
!   Calculate the Moments of the Density, Centered
!   in the Box and Contracted with the Multipole Tensor
!
    DO zq=1,Rho%NExpt
       NQ    = Rho%NQ%I(zq)
       Zeta  = (Rho%Expt%D(zq))**(-threehalf)
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
             jadd = Rho%OffR%I(zq)+(iq-1)*LKet
!
             Px = Rho%Qx%D(iadd)-MM_Rho%CenterX
             Py = Rho%Qy%D(iadd)-MM_Rho%CenterY
             Pz = Rho%Qz%D(iadd)-MM_Rho%CenterZ
!
             IF(Px == Zero .AND. Py == Zero .AND. Pz == Zero) THEN
                CALL HGTFToSP(LQ,Zeta,Cq,Sq,Rho%Co%D(jadd+1))
                DO LM = 0,LenQ
                   MM_Rho%CMMat%D(LM) =  MM_Rho%CMMat%D(LM) + Cq(LM)
                   MM_Rho%SMMat%D(LM) =  MM_Rho%SMMat%D(LM) + Sq(LM)
                ENDDO
             ELSE
!
                CALL HGTFToSP(LQ,Zeta,Cq,Sq,Rho%Co%D(jadd+1))
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
!   Contract the Density Multipoles with the Multipole Tensor
!
    CALL New_Moments(MM_TenRho,BigL)
    MM_TenRho%CMMat%D = Zero
    MM_TenRho%SMMat%D = Zero
    MM_TenRho%CenterX = CellCenterX
    MM_TenRho%CenterY = CellCenterY
    MM_TenRho%CenterZ = CellCenterZ
!
!   Calculate Sum_l'm' M(l+l',m+m') * rho(l'.m')
!
    DO L1 = 0,BigL
       DO M1 = 0,L1
          LM1 = LTD(L1)+M1
          Csum = Zero
          Ssum = Zero
          DO L2 = 0,BigL
             DO M2 = -L2,L2
                LM2 = LTD(L2)+ABS(M2)
                L12 = L1+L2
                M12 = M1+M2
                LM12 = LTD(L12)+ABS(M12)
                IF(L12 .LE. BigL) THEN 
                   IF(M12 .GE. -L12 .OR. M12 .LE. L12) THEN
                      IF(M2 < 0) THEN
                         CM2 = (-One)**ABS(M2  )
                         SM2 = (-One)**ABS(M2+1)
                      ELSE
                         CM2 = One
                         SM2 = One
                      ENDIF
                      IF(M12 < 0) THEN
                         CM12 = (-One)**ABS(M12  )
                         SM12 = (-One)**ABS(M12+1)
                      ELSE
                         CM12 = One
                         SM12 = One
                      ENDIF
!
                      Csum = Csum                                              &
                        +CM12*CM2*MM_Tensor%CMMat%D(LM12)*MM_Rho%CMMat%D(LM2)  &
                        +SM12*SM2*MM_Tensor%SMMat%D(LM12)*MM_Rho%SMMat%D(LM2) 
                      Ssum = Ssum                                              &
                        +CM12*SM2*MM_Tensor%CMMat%D(LM12)*MM_Rho%SMMat%D(LM2)  &
                        -SM12*CM2*MM_Tensor%SMMat%D(LM12)*MM_Rho%CMMat%D(LM2) 
!                     
                   ENDIF
                ENDIF
             ENDDO
          ENDDO
          MM_TenRho%CMMat%D(LM1) = Csum
          MM_TenRho%SMMat%D(LM1) = Ssum
       ENDDO
    ENDDO

!
  END SUBROUTINE CalMMRho
!========================================================================================
!
!========================================================================================
  FUNCTION CTraxBraKet(LBra,Expt,Qx,Qy,Qz,Coef)
!
    INTEGER                             :: L,M,LM,LBra,LQ,LP,LPQ,LenQ,LenP,LenPQ
    REAL(DOUBLE)                        :: Expt,Qx,Qy,Qz
    REAL(DOUBLE)                        :: Zeta,Px,Py,Pz,CTraxBraKet
    REAL(DOUBLE),DIMENSION(:)           :: Coef
    REAL(DOUBLE),DIMENSION(0:BigLM)     :: Cpq,Spq,Cq,Sq,Cp,Sp
!
    INTEGER                             :: L1,L2,L12,M1,M2,M12,LM1,LM2,LM12
    REAL(DOUBLE)                        :: CM1,CM2,CM12,SM1,SM2,SM12
!
    Px   = Qx-MM_Rho%CenterX
    Py   = Qy-MM_Rho%CenterY
    Pz   = Qz-MM_Rho%CenterZ
    Zeta = Pi3*Expt**(-ThreeHalf)
!
    LQ    = LBra
    LP    = BigL-LQ
    LPQ   = LQ+LP
    LenQ  = LHGTF(LQ)
    LenP  = LHGTF(LP)
    LenPQ = LHGTF(LPQ)
    Cp(:) = Zero
    Sp(:) = Zero
    Cq(:) = Zero
    Sq(:) = Zero
!
    IF(Px == Zero .AND. Py == Zero .AND. Pz == Zero) THEN
       CALL HGTFToSP(LQ,Zeta,Cp,Sp,Coef)
    ELSE
       CALL HGTFToSP(LQ,Zeta,Cq,Sq,Coef)
       CALL XLate(BigL,LP,LQ,LPQ,LenP,LenQ,LenPQ,Px,Py,Pz, &
                  Cp,Sp,                                   &
                  Cpq,Spq,                                 &
                  Cq,Sq,                                   &
                  LegendreP%D,Cosine%D,Sine%D,RToTh%D,     &
                  FactOlm0%D,FactOlm1%D,FactOlm2%D)
    ENDIF
!
!   Calculate Sum_lm,l'm'  (-1)^l * rho_ab(l,m) * M(l+l',m+m') * rho(l'.m')
!
    CTraxBraKet = Zero
    DO L1 = 0,BigL
       DO M1 = -L1,L1
          LM1 = LTD(L1)+ABS(M1)
          CTraxBraKet = CTraxBraKet +((-One)**L1)*(          &
                      + Cp(LM1)*MM_TenRho%CMMat%D(LM1)       &
                      - Sp(LM1)*MM_TenRho%SMMat%D(LM1))
       ENDDO
    ENDDO
!
  END FUNCTION CTraxBraKet
#ifdef PERIODIC
!========================================================================================
!
!========================================================================================
  FUNCTION BoxIsCube(GM) 
    TYPE(CRDS)                          :: GM
    LOGICAL                             :: BoxIsCube
    REAL(DOUBLE)                        :: AA,BB,CC,AB,AC,BC,BoxIsCubeTol
!
    BoxIsCube    = .FALSE.
    BoxIsCubeTol = 1.0D-12
!
!   Test to see if latice vectors form a cube
!
    AA = GM%BoxShape%D(1,1)**2
    BB = GM%BoxShape%D(1,2)**2+GM%BoxShape%D(2,2)**2
    CC = GM%BoxShape%D(1,3)**2+GM%BoxShape%D(2,3)**2+GM%BoxShape%D(3,3)**2
    AB = GM%BoxShape%D(1,1)*GM%BoxShape%D(1,2)
    AC = GM%BoxShape%D(1,1)*GM%BoxShape%D(1,3)
    BC = GM%BoxShape%D(1,2)*GM%BoxShape%D(1,3)+GM%BoxShape%D(2,2)*GM%BoxShape%D(2,3)
!
    IF(AB < BoxIsCubeTol .AND. AC < BoxIsCubeTol .AND. BC < BoxIsCubeTol) THEN
       IF( ABS(AA-BB) < BoxIsCubeTol .AND. ABS(BB-CC) < BoxIsCubeTol) THEN
          BoxIsCube = .TRUE.
       ENDIF
    ENDIF 
!
!   Test the PBC
!
    IF(.NOT. (GM%AutoW(1) .AND. GM%AutoW(2) .AND. GM%AutoW(3))) THEN
       BoxIsCube = .FALSE.
    ENDIF
!
  END FUNCTION BoxIsCube
!========================================================================================
!
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
          X = Rho%QRx%D(iadd)-CellCenterX
          Y = Rho%QRy%D(iadd)-CellCenterY
          Z = Rho%QRz%D(iadd)-CellCenterZ
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
  END SUBROUTINE BoxBounds
#endif
!
!
!
END MODULE MMoments



















