MODULE PBCFarField
  USE Derivedtypes
  USE GlobalScalars   
  USE GlobalObjects
  USE ProcessControl
  USE Indexing
  USE Parse
  USE InOut
  USE Macros
  USE Thresholding
  USE BoundingBox
  USE McMurchie
  USE MondoPoles
  USE PoleTree
  USE SpecFun
  USE Globals
  USE PFFTensor 
  IMPLICIT NONE
!====================================================================================
! Global Array intermediates for multipole opperations
!====================================================================================
  REAL(DOUBLE)                      :: DFac,QFac,Quadripole
  REAL(DOUBLE), DIMENSION(3)        :: CellCenter,Dipole
  REAL(DOUBLE), DIMENSION(0:FFLen)  :: RhoC,RhoS
  CONTAINS
#ifdef PERIODIC
!====================================================================================
!   Setup the PBCFarField Matrix. 
!====================================================================================
    SUBROUTINE PBCFarFieldSetUp(FFL)
      INTEGER                         :: FFL,I,J,K,IMin,KMin,JMin,NC,LM
      REAL(DOUBLE),DIMENSION(3)       :: PQ
      REAL(DOUBLE)                    :: E_PFF,E_DP,E_QP
      REAL(DOUBLE), DIMENSION(0:FFLen):: FarFC,FarFS
!
!     Calculate the Dimension
!
      Dimen = 0
      DO I = 1,3; IF(GM%AutoW(I)) Dimen = Dimen+1;ENDDO
!
!     Calculate the Box Volume
!
      Volume = One
      DO I=1,3;IF(GM%AutoW(I)) Volume = Volume*GM%BoxShape%D(I,I);ENDDO
!
!     Calculate the Dipole and Quadripole Corections factors
!
      IF(Dimen < 3) THEN
         DFac = Zero
         QFac = Zero
      ELSEIF(Dimen==3) THEN
         DFac = -(Four*Pi)/(Three*Volume)
         QFac = (Two*Pi)/(Three*Volume)
      ENDIF
!
!     Calculate the Center of the Cell
!
      CellCenter=Zero
      DO I = 1,3
         IF(GM%AutoW(I)) THEN
            CellCenter(1) = CellCenter(1)+Half*GM%BoxShape%D(I,1)
            CellCenter(2) = CellCenter(2)+Half*GM%BoxShape%D(I,2)
            CellCenter(3) = CellCenter(3)+Half*GM%BoxShape%D(I,3)
         ENDIF
      ENDDO
!
!     Calculate the Size of the Box Needed  for the Direct J
!
      CALL BoxBounds(IMin,JMin,KMin)
!
!     Generate the Cells for the Inner Box
!  
      CALL New_CellSet_Cube(CSMM1,GM%BoxShape%D,(/IMin,JMin,KMin/))
!
!     Calculate the Box Moments
!
      CALL RhoToSP(RhoC,RhoS)
!
!     Calculate the PBC FarField Tensor
!
      CALL PFFTensor(IMin,JMin,KMin)
!
!     CALCULATE THE DIPOLE AND QUADRIPOLE MOMENTS
!
      CALL RhoToQ2(Quadripole,Dipole(1),Dipole(2),Dipole(3))
!
!     Calculate PFF    energy
!
      E_PFF  = Zero
      FarFC  = Zero
      FarFS  = Zero
      CALL CTraX77(FFELL,FFEll,FarFC,FarFS,TensorC,TensorS,RhoC,RhoS)
      DO LM = 0,LSP(FFEll)
         E_PFF = E_PFF + RhoC(LM)*FarFC(LM)+RhoS(LM)*FarFS(LM)
      ENDDO
      E_PFF = Two*E_PFF
!
!     Calculate Dipole energy
!
      E_DP = Two*DFAC*(Dipole(1)**2+Dipole(2)**2+Dipole(3)**2)
!
!     Calculate Quad   energy
!
      E_QP = ABS(GM%NElec)*QFAC*Quadripole
!
      WRITE(*,*) 'CellCenterX = ',CellCenter(1)
      WRITE(*,*) 'CellCenterY = ',CellCenter(2)
      WRITE(*,*) 'CellCenterZ = ',CellCenter(3)
      WRITE(*,*) 'IMin        = ',IMin
      WRITE(*,*) 'JMin        = ',JMin
      WRITE(*,*) 'KMin        = ',KMin
      WRITE(*,*) 'D(1)        = ',Dipole(1)
      WRITE(*,*) 'D(2)        = ',Dipole(2)
      WRITE(*,*) 'D(3)        = ',Dipole(3)
      WRITE(*,*) 'Quadripole  = ',Quadripole
      WRITE(*,*) ' '
      WRITE(*,*) 'Correction to the Energy: PFF        = ',E_PFF
      WRITE(*,*) 'Correction to the Energy: Dipole     = ',E_DP
      WRITE(*,*) 'Correction to the Energy: Quadripole = ',E_QP
!
    END SUBROUTINE PBCFarFieldSetUp
!====================================================================================
!   Calculate the FarField Component of the J matrix
!====================================================================================
   FUNCTION CTraxFF(Prim,HGBra)
      TYPE(PrimPair)                   :: Prim
      INTEGER                          :: LM
      REAL(DOUBLE)                     :: PiZ,CTraxFF
      REAL(DOUBLE), DIMENSION(3)       :: PQ
      REAL(DOUBLE), DIMENSION(1:)      :: HGBra
      REAL(DOUBLE), DIMENSION(0:FFLen) :: SPBraC,SPBraS,SPKetC,SPKetS,FarFC,FarFS
!
      CTraxFF = Zero
!      IF(.TRUE.) RETURN
!
!     Transform <Bra| coefficients from HG to SP
!
      PiZ=(Pi/Prim%Zeta)**(ThreeHalves)
      CALL HGToSP_Gen(Prim%Ell,PiZ,HGBra,SPBraC,SPBraS)   
      PQ = CellCenter-Prim%P
!
      IF(.NOT. NoTranslate(PQ)) THEN
         SPKetC = Zero
         SPKetS = Zero
         FarFC  = Zero
         FarFS  = Zero
         CALL Regular(FFELL,PQ)
         CALL XLate77(FFELL,FFEll,SPKetC,SPKetS,Cpq,Spq,RhoC,RhoS)
         CALL CTraX77(FFELL,FFEll,FarFC,FarFS,TensorC,TensorS,SPKetC,SPKetS)
      ELSE
         FarFC  = Zero
         FarFS  = Zero
         CALL CTraX77(FFELL,FFEll,FarFC,FarFS,TensorC,TensorS,RhoC,RhoS)
      ENDIF
!
!     Contrax the <Bra|Ket> FF corection
!
      DO LM = 0,LSP(Prim%Ell)
         CTraxFF = CTraxFF + SPBraC(LM)*FarFC(LM)+SPBraS(LM)*FarFS(LM)
      ENDDO
!
!     INCLUDE THE DIPOLE correction to FarFC and FarFS
!
      SELECT CASE(Prim%Ell)
      CASE (0)
         CTraxFF = CTraxFF + (-PiZ*HGBra(1)*PQ(1))*(DFAC*Dipole(1)) &
                           + (-PiZ*HGBra(1)*PQ(2))*(DFAC*Dipole(2)) &
                           + (-PiZ*HGBra(1)*PQ(3))*(DFAC*Dipole(3))
      CASE(1:)
         CTraxFF = CTraxFF + (-PiZ*HGBra(1)*PQ(1)+PiZ*HGBra(2))*(DFAC*Dipole(1)) &
                           + (-PiZ*HGBra(1)*PQ(2)+PiZ*HGBra(3))*(DFAC*Dipole(2)) &
                           + (-PiZ*HGBra(1)*PQ(3)+PiZ*HGBra(4))*(DFAC*Dipole(3))
      END SELECT
!
    END FUNCTION CTraxFF
!====================================================================================
!   Calculate the Q of the Coulomb Matrix
!====================================================================================
    FUNCTION CTraxQ(Prim,HGBra)
      TYPE(PrimPair)                   :: Prim
      REAL(DOUBLE)                     :: PiZ,CTraxQ
      REAL(DOUBLE), DIMENSION(1:)      :: HGBra
!
      CTraxQ = Zero
!      IF(.TRUE.) RETURN
!
      PiZ=(Pi/Prim%Zeta)**(ThreeHalves)
!
!     INCLUDE THE QUADRIPOLE correction 
!
      CTraxQ = CTraxQ + PiZ*HGBra(1)*QFAC*Quadripole
!
    END FUNCTION CTraxQ
!========================================================================================
! Calculate Quadripole = Q(2,0,0)+Q(0,2,0)+Q(0,0,2) of the local density
!========================================================================================
    SUBROUTINE RhoToQ2(Q2,DX,DY,DZ)
      INTEGER                  :: zq,iq,iadd,jadd,NQ,OffQ,OffR,LQ,LenQ
      REAL(DOUBLE)             :: Q2,RX,RY,RZ,R2,Expt,PiExpt,Expt75,DX,DY,DZ,SUM
!
      Q2 = Zero
      DX = Zero
      DY = Zero
      DZ = Zero
      DO zq=1,Rho%NExpt
         NQ     = Rho%NQ%I(zq)
         Expt   = Rho%Expt%D(zq)
         PiExpt = (Pi/Expt)**(threehalves)
         Expt75 = 0.75D0/Expt
         OffQ   = Rho%OffQ%I(zq)
         OffR   = Rho%OffR%I(zq)
         LQ     = Rho%Lndx%I(zq) 
         LenQ   = LHGTF(LQ) 
         IF(NQ > 0) THEN
            DO iq = 1,NQ
               iadd = Rho%OffQ%I(zq)+iq
               jadd = Rho%OffR%I(zq)+(iq-1)*LenQ+1
               RX   = Rho%Qx%D(iadd)-CellCenter(1)
               RY   = Rho%Qy%D(iadd)-CellCenter(2)
               RZ   = Rho%Qz%D(iadd)-CellCenter(3)
               R2   = Half*(RX*RX+RY*RY+RZ*RZ)+Expt75
!
               SELECT CASE(LQ)
               CASE (0)
                  Q2 = Q2 + PiExpt*(Rho%Co%D(jadd)*R2)
!
                  DX = DX + PiExpt*(Rho%Co%D(jadd)*RX)
                  DY = DY + PiExpt*(Rho%Co%D(jadd)*RY)
                  DZ = DZ + PiExpt*(Rho%Co%D(jadd)*RZ)
               CASE(1)
                  Q2 = Q2+ PiExpt*(Rho%Co%D(jadd)*R2 + &
                       RX*Rho%Co%D(jadd+1) + RY*Rho%Co%D(jadd+2) + RZ*Rho%Co%D(jadd+3))
!
                  DX = DX + PiExpt*(Rho%Co%D(jadd)*RX+Rho%Co%D(jadd+1))
                  DY = DY + PiExpt*(Rho%Co%D(jadd)*RY+Rho%Co%D(jadd+2))
                  DZ = DZ + PiExpt*(Rho%Co%D(jadd)*RZ+Rho%Co%D(jadd+3))
               CASE(2:)
                  Q2 = Q2+ PiExpt*(Rho%Co%D(jadd)*R2 + &
                       RX*Rho%Co%D(jadd+1) + RY*Rho%Co%D(jadd+2) + RZ*Rho%Co%D(jadd+3) + &
                       Rho%Co%D(jadd+4) + Rho%Co%D(jadd+6) + Rho%Co%D(jadd+9))
!
                  DX = DX + PiExpt*(Rho%Co%D(jadd)*RX+Rho%Co%D(jadd+1))
                  DY = DY + PiExpt*(Rho%Co%D(jadd)*RY+Rho%Co%D(jadd+2))
                  DZ = DZ + PiExpt*(Rho%Co%D(jadd)*RZ+Rho%Co%D(jadd+3))
               END SELECT
!
            ENDDO
         ENDIF
      ENDDO
!
    END SUBROUTINE RhoToQ2  
!========================================================================================
! Calculate the SP Moments of Rho_Loc
!========================================================================================
  SUBROUTINE RhoToSP(Cp,Sp)
      INTEGER                         :: LP,LQ,LPQ,LenP,LenQ,LenPQ,LKet
      INTEGER                         :: zq,iq,iadd,jadd,NQ,OffQ,OffR,LM
      REAL(DOUBLE)                    :: Zeta,PiZ
      REAL(DOUBLE),DIMENSION(3)       :: PQ
      REAL(DOUBLE),DIMENSION(0:FFLen) :: Cp,Sp,Cq,Sq
!
      Cp  = Zero
      Sp  = Zero
      DO zq=1,Rho%NExpt
         NQ    = Rho%NQ%I(zq)
         Zeta  = Rho%Expt%D(zq)
         OffQ  = Rho%OffQ%I(zq)
         OffR  = Rho%OffR%I(zq)
!
         LQ    = Rho%Lndx%I(zq)
         LP    = FFEll
         LPQ   = LP+LQ
         LenQ  = LSP(LQ)  
         LenP  = LSP(LP)  
         LenPQ = LSP(LPQ) 
         LKet  = LHGTF(LQ)
!
         IF(NQ > 0) THEN
            DO iq = 1,NQ
               iadd = Rho%OffQ%I(zq)+iq
               jadd = Rho%OffR%I(zq)+(iq-1)*LKet+1
               PQ(1)  = Rho%Qx%D(iadd)-CellCenter(1)
               PQ(2)  = Rho%Qy%D(iadd)-CellCenter(2)
               PQ(3)  = Rho%Qz%D(iadd)-CellCenter(3)
!
               PiZ=(Pi/Zeta)**(ThreeHalves)
               CALL HGToSP_Gen(LQ,PiZ,Rho%Co%D(jadd:),Cq,Sq) 
!
               IF(NoTranslate(PQ))THEN
                  DO LM = 0,LenQ
                     Cp(LM) =  Cp(LM) + Cq(LM)   
                     Sp(LM) =  Sp(LM) + Sq(LM)
                  ENDDO
               ELSE
                  CALL Regular(LPQ,PQ)
                  CALL XLate90(LP,LQ,Cp,Sp,Cpq,Spq,Cq,Sq)
               ENDIF
            ENDDO
         ENDIF
      ENDDO
!
    END SUBROUTINE RhoToSP
!========================================================================================
! Calculate the Box Bounds Needed for the Direct Sum
!========================================================================================
  SUBROUTINE BoxBounds(I,J,K) 
    INTEGER                         :: I,J,K,iadd,iq,zq,NQ
    REAL(DOUBLE)                    :: PacMin,MacMin,BoxR,RMin
!
!   Box Size
!
    BoxR = Two*SQRT(CellCenter(1)**2+CellCenter(2)**2+CellCenter(3)**2)
!
!   PAC Distance (Approx)
!
    PacMin = PFunk(2*BS%NAsym,Thresholds%TwoE)
    PacMin = SQRT(Two*PacMin/Rho%Expt%D(1))
!
!   MAC DISTANCE
!
    MacMin = ( (BoxR**(DBLE(FFEll+1)))/Thresholds%TwoE )**(One/DBLE(FFEll+2))
!
!   Inner Box Distance: MAC AND PAC
!
    RMin  = MAX(PacMin,MacMin)
    I  = RMin/GM%BoxShape%D(1,1)
    J  = RMin/GM%BoxShape%D(2,2)
    K  = RMin/GM%BoxShape%D(3,3)
!
!   Threshold
!
    I = MIN(I,5)
    J = MIN(J,5)
    K = MIN(K,5)
    I = MAX(I,1)
    J = MAX(J,1)
    K = MAX(K,1)
!
!   Impose BCs
!
    IF(.NOT. GM%AutoW(1)) I = 0
    IF(.NOT. GM%AutoW(2)) J = 0
    IF(.NOT. GM%AutoW(3)) K = 0
!
    WRITE(*,*) 'BoxR   = ',BoxR
    WRITE(*,*) 'PacMin = ',PacMin
    WRITE(*,*) 'MacMin = ',MacMin
!
  END SUBROUTINE BoxBounds
!========================================================================================
! Calculate the Box Bounds Needed for the Direct Sum
!========================================================================================
  FUNCTION NoTranslate(X) 
    REAL(DOUBLE),DIMENSION(3) :: X
    REAL(DOUBLE),PARAMETER    :: TOL=1.0D-14
    LOGICAL                   :: NoTranslate
!
    NoTranslate = (ABS(X(1)).LT.TOL) .AND. (ABS(X(2)).LT.TOL) .AND. (ABS(X(3)).LT.TOL)
!
  END FUNCTION NoTranslate
!========================================================================================
! Print Fields
!========================================================================================
  SUBROUTINE Print_SP(Ell,C,S,Tag_O,Pre_O)
    CHARACTER(LEN=*),OPTIONAL        :: Tag_O,Pre_O 
    CHARACTER(LEN=DEFAULT_CHR_LEN)   :: Line
    INTEGER                          :: Ell,L,M,LMDx,IPre
    REAL(DOUBLE)                     :: Cpp,Spp
    REAL(DOUBLE),DIMENSION(0:)       :: C,S
    WRITE(*,*)'======================================================' 
    IF(PRESENT(Tag_O)) WRITE(*,*) Tag_O
    IF(PRESENT(Pre_O)) THEN
       IF(Pre_O == 'Short') IPre = 0
       IF(Pre_O == 'Med'  ) IPre = 1
       IF(Pre_O == 'Long' ) IPre = 2       
    ELSE
       IPre = 1
    ENDIF
!
    IF(IPre == 0) THEN
       DO l=0,Ell
          DO m=0,l
             lmdx=LTD(l)+m
             Cpp = C(lmdx)
             Spp = S(lmdx)
             IF(ABS(Cpp) .LT. 1.D-12) Cpp = Zero 
             IF(ABS(Spp) .LT. 1.D-12) Spp = Zero 
             Line='L = '//TRIM(IntToChar(L))//' M = '//TRIM(IntToChar(M)) &
                  //' Cq = '//TRIM(DblToShrtChar(Cpp)) &
                  //' Sq = '//TRIM(DblToShrtChar(Spp))
             WRITE(*,*)TRIM(Line)
          ENDDO
       ENDDO
    ELSEIF(IPre == 1) THEN
       DO l=0,Ell
          DO m=0,l
             lmdx=LTD(l)+m
             Cpp = C(lmdx)
             Spp = S(lmdx)
             IF(ABS(Cpp) .LT. 1.D-12) Cpp = Zero 
             IF(ABS(Spp) .LT. 1.D-12) Spp = Zero 
             Line='L = '//TRIM(IntToChar(L))//' M = '//TRIM(IntToChar(M)) &
                  //' Cq = '//TRIM(DblToMedmChar(Cpp)) &
                  //' Sq = '//TRIM(DblToMedmChar(Spp))
             WRITE(*,*)TRIM(Line)
          ENDDO
       ENDDO
    ELSEIF(IPre == 2) THEN
       DO l=0,Ell
          DO m=0,l
             lmdx=LTD(l)+m
             Cpp = C(lmdx)
             Spp = S(lmdx)
             IF(ABS(Cpp) .LT. 1.D-12) Cpp = Zero 
             IF(ABS(Spp) .LT. 1.D-12) Spp = Zero 
             Line='L = '//TRIM(IntToChar(L))//' M = '//TRIM(IntToChar(M)) &
                  //' Cq = '//TRIM(DblToChar(Cpp)) &
                  //' Sq = '//TRIM(DblToChar(Spp))
             WRITE(*,*)TRIM(Line)
          ENDDO
       ENDDO
    ENDIF
!
  END SUBROUTINE Print_SP
!========================================================================================
! Print Cells in Two D
!========================================================================================
  SUBROUTINE Print_Occ(CS)
    TYPE(CellSet)         :: CS
    INTEGER               :: NC,IMAX,JMAX,KMAX,I,J
    REAL(DOUBLE)          :: XMAX,YMAX,ZMAX,X,Y,Z
    LOGICAL               :: InCell
    INTEGER,DIMENSION(20) :: OCC
!    
    XMAX = Zero
    YMAX = Zero
    ZMAx = Zero
    DO NC=1,CS%NCells
       XMAX = MAX(XMAX,CS%CellCarts%D(1,NC))
       YMAX = MAX(YMAX,CS%CellCarts%D(2,NC))
       ZMAX = MAX(ZMAX,CS%CellCarts%D(3,NC))
    ENDDO
    IMAX  = XMAX/GM%BoxShape%D(1,1)
    JMAX  = YMAX/GM%BoxShape%D(2,2)
    KMAX  = ZMAX/GM%BoxShape%D(3,3)
    DO I=-IMAX,IMAX
       WRITE(*,*) ' '
       DO J=-JMAX,JMAX
          X  = I*GM%BoxShape%D(1,1)
          Y  = I*GM%BoxShape%D(1,2)+J*GM%BoxShape%D(2,2)
          Z  = I*GM%BoxShape%D(1,3)+J*GM%BoxShape%D(2,3)
          IF(InCell_CellSet(CS,X,Y,Z)) THEN
             OCC(J+JMAX+1) = 1
          ELSE
             OCC(J+JMAX+1) = 0
          ENDIF
       ENDDO                
       WRITE(*,'20(I2)') (OCC(J),J=1,2*JMAX+1)
    ENDDO
!
  END SUBROUTINE Print_Occ
#endif
!
!
!
END MODULE PBCFarField

