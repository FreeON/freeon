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
  IMPLICIT NONE
!====================================================================================
! Global Array intermediates for multipole opperations
!====================================================================================
  INTEGER                           :: Dimen
  REAL(DOUBLE)                      :: Volume,DFac,QFac,Quadripole
  REAL(DOUBLE), DIMENSION(3)        :: CellCenter,Dipole
  REAL(DOUBLE), DIMENSION(0:FFLen)  :: RhoC,RhoS
  REAL(DOUBLE), DIMENSION(0:FFLen2) :: TensorC,TensorS
  CONTAINS
#ifdef PERIODIC
!====================================================================================
!   Setup the PBCFarField Matrix. 
!====================================================================================
    SUBROUTINE PBCFarFieldSetUp(FFL)
      INTEGER                         :: FFL,I,J,K,IMin,KMin,JMin,NC
      REAL(DOUBLE),DIMENSION(3)       :: PQ
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
!!
!   Calculate the Center of the Cell
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
!      CALL PBCTensor(IMin,JMin,KMin)
!
!     Calculate the Quadripole Moment
!
      CALL RhoToQ2(Quadripole,Dipole(1),Dipole(2),Dipole(3))
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
      WRITE(*,*) 'Correction to the Energy: Dipole     = ',DFAC*(Dipole(1)**2+Dipole(2)**2+Dipole(3)**2)
      WRITE(*,*) 'Correction to the Energy: Quadripole = ',Half*QFAC*Quadripole*NAtoms*GM%AtNum%I(1)
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
      REAL(DOUBLE), DIMENSION(:)       :: HGBra
      REAL(DOUBLE), DIMENSION(0:FFLen) :: SPBraC,SPBraS,SPKetC,SPKetS,FarFC,FarFS
!
      CTraxFF = Zero
!
!     Transform <Bra| coefficients from HG to SP
!
      PiZ=(Pi/Prim%Zeta)**(ThreeHalves)
      CALL HGToSP_Gen(Prim%Ell,PiZ,HGBra,SPBraC,SPBraS)   
      PQ = CellCenter-Prim%P
!
!!$      IF(.NOT. NoTranslate(PQ)) THEN
!!$         SPKetC = Zero
!!$         SPKetS = Zero
!!$         FarFC  = Zero
!!$         FarFS  = Zero
!!$         CALL Regular(FFELL,PQ)
!!$         CALL XLate77(FFELL,FFEll,SPKetC,SPKetS,Cpq,Spq,RhoC,RhoS)
!!$         CALL CTraX77(FFELL,FFEll,FarFC,FarFS,TensorC,TensorS,SPKetC,SPKetS)
!!$      ELSE
!!$         FarFC  = Zero
!!$         FarFS  = Zero
!!$         CALL CTraX77(FFELL,FFEll,FarFC,FarFS,TensorC,TensorS,RhoC,RhoS)
!!$      ENDIF
!!$!
!!$!     Contrax the <Bra|Ket>
!!$!
!!$      DO LM = 0,LSP(Prim%Ell)
!!$         CTraxFF = CTraxFF + SPBraC(LM)*FarFC(LM)+SPBraS(LM)*FarFS(LM)
!!$      ENDDO
!!
!     Include the Dipole correction to FarFC and FarFS
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
!========================================================================================
! Calculate Quadripole = Q(2,0,0)+Q(0,2,0)+Q(0,0,2) of the local density
!========================================================================================
    SUBROUTINE RhoToQ2(Q2,DX,DY,DZ)
      INTEGER                  :: zq,iq,iadd,jadd,NQ,OffQ,OffR,LQ,LenQ
      REAL(DOUBLE)             :: Q2,RX,RY,RZ,R2,Expt,PiExpt,Expt75,DX,DY,DZ
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
    I = MIN(I,1)
    J = MIN(J,1)
    K = MIN(K,1)
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
!========================================================================================
! Calculate the PBCTensor
!========================================================================================
  SUBROUTINE PBCTensor(IMin,JMin,KMin,MaxL_O)
    INTEGER,OPTIONAL                  :: MaxL_O
    INTEGER                           :: IMin,JMin,KMin
    INTEGER                           :: I,J,K,L,M,LM,NC,LMax
    INTEGER                           :: IMax,JMax,KMax,IJKMax,LSwitch
    REAL(DOUBLE),DIMENSION(3)         :: PQ
    REAL(DOUBLE)                      :: CFac,SFac,BetaSq,Rad,RadSq,ExpFac
!
    IF(PRESENT(MaxL_O)) THEN
       LMax = MaxL_O
    ELSE
       LMax = FFEll2
    ENDIF
!
    TensorC = Zero
    TensorS = Zero
    IF(Dimen==0) THEN
       TensorC = Zero
       TensorS = Zero
       RETURN
    ENDIF
!
!   One Dimension
!
    IF(Dimen==1) THEN
       IF(GM%AutoW(1)) THEN
          NC = IMin
          CALL IrRegular(LMax,(/GM%BoxShape%D(1,1),Zero,Zero/))
       ELSEIF(GM%AutoW(2)) THEN
          NC = JMin
          CALL IrRegular(LMax,(/Zero,GM%BoxShape%D(2,2),Zero/))
       ELSEIF(GM%AutoW(3)) THEN  
          NC = KMin       
          CALL IrRegular(LMax,(/Zero,Zero,GM%BoxShape%D(3,3)/))
       ENDIF
       DO L=1,LMax
          DO M = 0,L
             LM = LTD(L)
             TensorC(LM) = Cpq(LM)*RZeta(L+1,NC)
          ENDDO
       ENDDO
!
!   Two Dimension
!
    ELSEIF(Dimen==2) THEN
       LSwitch = 8
       BetaSq  = (Pi*Pi)/Volume
       IJKMax  = 32
       IMax    = IJKMax+IMin
       JMax    = IJKMax+KMin
       KMax    = IJKMax+JMin
!
!      Sum the Real Space
!
       IF(.NOT. GM%AutoW(1)) THEN
          CALL New_CellSet_Cube(CSMM2,GM%BoxShape%D,(/1,JMax,KMax/),(/1,JMin+1,KMin+1/))
       ELSEIF(.NOT. GM%AutoW(2)) THEN
          CALL New_CellSet_Cube(CSMM2,GM%BoxShape%D,(/IMax,1,KMax/),(/IMin+1,1,KMin+1/))
       ELSEIF(.NOT. GM%AutoW(3)) THEN
          CALL New_CellSet_Cube(CSMM2,GM%BoxShape%D,(/IMax,JMax,1/),(/IMin+1,JMin+1,1/))
       ENDIF
       DO NC = 1,CSMM2%NCells
          PQ(:) = CSMM2%CellCarts%D(:,NC)
          RadSq = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
          CALL IrRegular(LMax,PQ)
          DO L = 1,LMax
             IF(L .LE. LSwitch) THEN
                CFac = GScript(L,RadSq)
             ELSE
                CFac = One
             ENDIF
             DO M = 0,L
                LM = LTD(L)+M
                TensorC(LM)=TensorC(LM)+Cpq(LM)*CFac
                TensorS(LM)=TensorS(LM)+Spq(LM)*CFac
             ENDDO
          ENDDO
       ENDDO
       CALL Delete_CellSet(CSMM2)
!
!      Sum the Reciprical Space 
!
       ExpFac = (Pi*Pi)/(BetaSq)
       IF(.NOT. GM%AutoW(1)) THEN
          CALL New_CellSet_Cube(CSMM2,GM%InvBoxSh%D,(/1,JMax,KMax/),(/1,1,1/))
       ELSEIF(.NOT. GM%AutoW(2)) THEN
          CALL New_CellSet_Cube(CSMM2,GM%InvBoxSh%D,(/IMax,1,KMax/),(/1,1,1/))
       ELSEIF(.NOT. GM%AutoW(3)) THEN
          CALL New_CellSet_Cube(CSMM2,GM%InvBoxSh%D,(/IMax,JMax,1/),(/1,1,1/))
       ENDIF
       DO NC = 1,CSMM2%NCells
          PQ(:) = CSMM2%CellCarts%D(:,NC)
          Rad   = SQRT(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
          CALL IrRegular(LMax,PQ)
          DO L = 1,LMax
             IF(L .LE. LSwitch) THEN             
                CFac = FT_FScriptC(L,ExpFac,Rad)/Volume
                SFac = FT_FScriptS(L,ExpFac,Rad)/Volume
             ELSE
                CFac = Zero
                SFac = Zero
             ENDIF
             DO M = 0,L
                LM = LTD(L)+M
                TensorC(LM)=TensorC(LM)+Cpq(LM)*CFac-Spq(LM)*SFac
                TensorS(LM)=TensorS(LM)+Spq(LM)*CFac+Cpq(LM)*SFac
             ENDDO
          ENDDO
       ENDDO
       CALL Delete_CellSet(CSMM2)
!
!      Substract the inner boxes
!
       IF(.NOT. GM%AutoW(1)) THEN
          CALL New_CellSet_Cube(CSMM2,GM%BoxShape%D,(/1,JMin,KMin/),(/1,1,1/))
       ELSEIF(.NOT. GM%AutoW(2)) THEN
          CALL New_CellSet_Cube(CSMM2,GM%BoxShape%D,(/IMin,1,KMin/),(/1,1,1/))
       ELSEIF(.NOT. GM%AutoW(3)) THEN
          CALL New_CellSet_Cube(CSMM2,GM%BoxShape%D,(/IMin,JMin,1/),(/1,1,1/))
       ENDIF
       DO NC = 1,CSMM2%NCells
          PQ(:) = CSMM2%CellCarts%D(:,NC)
          RadSq = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
          CALL IrRegular(LMax,PQ)
          DO L = 1,LMax
             IF(L .LE. LSwitch) THEN
                CFac = FScript(L,RadSq)
             ELSE
                CFac = Zero
             ENDIF
             DO M = 0,L
                LM = LTD(L)+M
                TensorC(LM)=TensorC(LM)-Cpq(LM)*CFac
                TensorS(LM)=TensorS(LM)-Spq(LM)*CFac
             ENDDO
          ENDDO
       ENDDO
       CALL Delete_CellSet(CSMM2)
!
!   Three Dimension
!
    ELSEIF(Dimen==3) THEN
       LSwitch = 8
       BetaSq  = Pi*Pi/(Volume**(Two/Three))
!
!      Sum the Real Space
!
       IJKMax = 32
       IMax = IJKMax+IMin
       JMax = IJKMax+KMin
       KMax = IJKMax+JMin
       CALL New_CellSet_Cube(CSMM2,GM%BoxShape%D,(/IMax,JMax,KMax/),(/IMin+1,JMin+1,KMin+1/))
       DO NC = 1,CSMM2%NCells
          PQ(:) = CSMM2%CellCarts%D(:,NC)
          RadSq = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
          CALL IrRegular(LMax,PQ)
          DO L = 1,LMax
             IF(L .LE. LSwitch) THEN
                CFac = GScript(L,RadSq)
             ELSE
                CFac = One
             ENDIF
             DO M = 0,L
                LM = LTD(L)+M
                TensorC(LM)=TensorC(LM)+Cpq(LM)*CFac
                TensorS(LM)=TensorS(LM)+Spq(LM)*CFac
             ENDDO
          ENDDO
       ENDDO
       CALL Delete_CellSet(CSMM2)
!
!      Sum the Reciprical Space 
!
       ExpFac = (Pi*Pi)/(BetaSq)
       CALL New_CellSet_Cube(CSMM2,GM%InvBoxSh%D,(/IMax,JMax,KMax/),(/1,1,1/))
       DO NC = 1,CSMM2%NCells
          PQ(:) = CSMM2%CellCarts%D(:,NC)
          Rad   = SQRT(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
          CALL IrRegular(LMax,PQ)
          DO L = 1,LMax
             IF(L .LE. LSwitch) THEN        
                CFac = FT_FScriptC(L,ExpFac,Rad)/Volume
                SFac = FT_FScriptS(L,ExpFac,Rad)/Volume
             ELSE
                CFac = Zero
                SFac = Zero
             ENDIF
             DO M = 0,L
                LM = LTD(L)+M
                TensorC(LM)=TensorC(LM)+Cpq(LM)*CFac-Spq(LM)*SFac
                TensorS(LM)=TensorS(LM)+Spq(LM)*CFac+Cpq(LM)*SFac
             ENDDO
          ENDDO
       ENDDO
       CALL Delete_CellSet(CSMM2)
!
!      Substract the inner boxes
!
       CALL New_CellSet_Cube(CSMM2,GM%BoxShape%D,(/IMin,JMin,KMin/),(/1,1,1/))  
       DO NC = 1,CSMM2%NCells
          PQ(:) = CSMM2%CellCarts%D(:,NC)
          RadSq = BetaSq*(PQ(1)*PQ(1)+PQ(2)*PQ(2)+PQ(3)*PQ(3))
          CALL IrRegular(LMax,PQ)
          DO L = 1,LMax
             IF(L .LE. LSwitch) THEN
                CFac = FScript(L,RadSq)
             ELSE
                CFac = Zero
             ENDIF 
             DO M = 0,L
                LM = LTD(L)+M
                TensorC(LM)=TensorC(LM)-Cpq(LM)*CFac
                TensorS(LM)=TensorS(LM)-Spq(LM)*CFac
             ENDDO
          ENDDO
       ENDDO
       CALL Delete_CellSet(CSMM2)
    ENDIF
!
!   Filter out the Zero Elements
!
    DO L = 0,LMax
       DO M = 0,L
          LM = LTD(L)+M
          IF(ABS(TensorC(LM)) .LT. 1.D-14) TensorC(LM) = Zero
          IF(ABS(TensorS(LM)) .LT. 1.D-14) TensorS(LM) = Zero
       ENDDO
    ENDDO
!
  END SUBROUTINE PBCTensor
!========================================================================================
! FT_FSCriptC
!========================================================================================
  FUNCTION FT_FScriptC(L,ExpFac,R)
    INTEGER                    :: L,IFac,Isgn
    REAL(DOUBLE)               :: ExpFac,R,FT_FScriptC,Fac
    REAL(DOUBLE),DIMENSION(41) :: Sfac = (/2.00000000000000000D0 ,1.61199195401646964D0 , &
                    1.39399011673806825D0 ,1.24836057075206815D0 ,1.14180213615392840D0 , &
                    1.05927697236749128D0 ,9.92824332104422014D-1,9.37767832227969827D-1, &
                    8.91149652589504141D-1,8.50992321055495993D-1,8.15915401864279526D-1, &
                    7.84921233738737833D-1,7.57267966537096089D-1,7.32390674279586654D-1, &
                    7.09850372573553877D-1,6.89299957250757043D-1,6.70460791874183787D-1, &
                    6.53106213773379367D-1,6.37049661171763433D-1,6.22135962734544435D-1, &
                    6.08234838286273008D-1,5.95235975457893168D-1,5.83045248969098145D-1, &
                    5.71581781316674666D-1,5.60775631817133548D-1,5.50565960943842365D-1, &
                    5.40899558418448758D-1,5.31729652703953292D-1,5.23014940361368490D-1, &
                    5.14718788772630577D-1,5.06808576734283345D-1,4.99255145565447538D-1, &
                    4.92032339458264325D-1,4.85116618392624313D-1,4.78486730436807416D-1, & 
                    4.72123432945047398D-1,4.66009254246335764D-1,4.60128289044821961D-1, &
                    4.54466022030378857D-1,4.49009175209461305D-1,4.43745575272018022D-1 /) 
!
    IFac = 1+(-1)**L
    IF(IFac == 0) THEN
       FT_FScriptC = Zero
    ELSE
       Isgn = (-1)**(L/2)
       Fac  = ExpFac/DBLE(2*L-1)
       FT_FScriptC = Isgn*(Sfac(L)*R*EXP(-Fac*R*R))**(2*L-1)   
    ENDIF
! 
  END FUNCTION FT_FScriptC
!========================================================================================
! FT_FSCriptS
!========================================================================================
  FUNCTION FT_FScriptS(L,ExpFac,R)
    INTEGER                    :: L,IFac,Isgn
    REAL(DOUBLE)               :: ExpFac,R,FT_FScriptS,Fac
    REAL(DOUBLE),DIMENSION(41) :: Sfac = (/2.00000000000000000D0 ,1.61199195401646964D0 , &
                    1.39399011673806825D0 ,1.24836057075206815D0 ,1.14180213615392840D0 , &
                    1.05927697236749128D0 ,9.92824332104422014D-1,9.37767832227969827D-1, &
                    8.91149652589504141D-1,8.50992321055495993D-1,8.15915401864279526D-1, &
                    7.84921233738737833D-1,7.57267966537096089D-1,7.32390674279586654D-1, &
                    7.09850372573553877D-1,6.89299957250757043D-1,6.70460791874183787D-1, &
                    6.53106213773379367D-1,6.37049661171763433D-1,6.22135962734544435D-1, &
                    6.08234838286273008D-1,5.95235975457893168D-1,5.83045248969098145D-1, &
                    5.71581781316674666D-1,5.60775631817133548D-1,5.50565960943842365D-1, &
                    5.40899558418448758D-1,5.31729652703953292D-1,5.23014940361368490D-1, &
                    5.14718788772630577D-1,5.06808576734283345D-1,4.99255145565447538D-1, &
                    4.92032339458264325D-1,4.85116618392624313D-1,4.78486730436807416D-1, & 
                    4.72123432945047398D-1,4.66009254246335764D-1,4.60128289044821961D-1, &
                    4.54466022030378857D-1,4.49009175209461305D-1,4.43745575272018022D-1 /) 
!
    IFac = 1-(-1)**L
    IF(IFac == 0) THEN
       FT_FScriptS = Zero
    ELSE
       Isgn = (-1)**((L-1)/2)
       Fac  = ExpFac/DBLE(2*L-1)
       FT_FScriptS = Isgn*(Sfac(L)*R*EXP(-Fac*R*R))**(2*L-1)
    ENDIF
!
  END FUNCTION FT_FScriptS
!========================================================================================
! FT_FSCriptC
!========================================================================================
  FUNCTION GScript(L,R)
    INTEGER                    :: L,LS
    REAL(DOUBLE)               :: R,SqrtR,GScript,XSUM
    REAL(DOUBLE),DIMENSION(41) :: Sfac = (/2.0000000000000000D0 ,1.333333333333333333D0 ,&
                       7.3029674334022148D-1,5.3412580601946498D-1,4.2897258944050779D-1,&
                       3.6130260767122454D-1,3.1338173978619282D-1,2.7736675692301273D-1,&
                       2.4916970717917522D-1,2.2642030439700098D-1,2.0763707534036204D-1,&
                       1.9184081786310919D-1,1.7835560611619894D-1,1.6669836456773905D-1,&
                       1.5651386447776487D-1,1.4753458765116469D-1,1.3955494670757578D-1,&
                       1.3241416731800685D-1,1.2598459485542503D-1,1.2016350807025758D-1,&
                       1.1486726370711862D-1,1.1002702835105538D-1,1.0558561445605698D-1,&
                       1.0149509929347097D-1,9.7715008595692035D-2,9.4210913822889116D-2,&
                       9.0953336661706548D-2,8.7916884656620846D-2,8.5079562763772593D-2,&
                       8.2422220248006543D-2,7.9928102738676689D-2,7.7582486742601370D-2,&
                       7.5372379364895518D-2,7.3286270006162048D-2,7.1313923796257465D-2,&
                       6.9446208774383175D-2,6.7674950532187841D-2,6.5992809342867331D-2,&
                       6.4393175806982839D-2,6.2870081828999622D-2,6.1418124351705448D-2 /)
!
    SqrtR      = SQRT(R)
    IF(L == 0) THEN
       GScript = (One-ERF(SqrtR))
    ELSE
       XSUM = SFAC(1)*EXP(-R)
       DO LS = 2,L
          XSUM  = XSUM+(Sfac(LS)*R*EXP(-R/DBLE(LS-1)))**(LS-1)
       ENDDO 
       GScript = (One-ERF(SqrtR))+(SqrtR/SqrtPi)*XSUM
    ENDIF
!
  END FUNCTION GScript
!========================================================================================
! FT_FSCriptS
!========================================================================================
  FUNCTION FScript(L,R)
    INTEGER                    :: L,LS
    REAL(DOUBLE)               :: R,SqrtR,FScript,XSUM
    REAL(DOUBLE),DIMENSION(41) :: Sfac = (/2.0000000000000000D0 ,1.333333333333333333D0 ,&
                       7.3029674334022148D-1,5.3412580601946498D-1,4.2897258944050779D-1,&
                       3.6130260767122454D-1,3.1338173978619282D-1,2.7736675692301273D-1,&
                       2.4916970717917522D-1,2.2642030439700098D-1,2.0763707534036204D-1,&
                       1.9184081786310919D-1,1.7835560611619894D-1,1.6669836456773905D-1,&
                       1.5651386447776487D-1,1.4753458765116469D-1,1.3955494670757578D-1,&
                       1.3241416731800685D-1,1.2598459485542503D-1,1.2016350807025758D-1,&
                       1.1486726370711862D-1,1.1002702835105538D-1,1.0558561445605698D-1,&
                       1.0149509929347097D-1,9.7715008595692035D-2,9.4210913822889116D-2,&
                       9.0953336661706548D-2,8.7916884656620846D-2,8.5079562763772593D-2,&
                       8.2422220248006543D-2,7.9928102738676689D-2,7.7582486742601370D-2,&
                       7.5372379364895518D-2,7.3286270006162048D-2,7.1313923796257465D-2,&
                       6.9446208774383175D-2,6.7674950532187841D-2,6.5992809342867331D-2,&
                       6.4393175806982839D-2,6.2870081828999622D-2,6.1418124351705448D-2 /)
!
    SqrtR      = SQRT(R)
    IF(L == 0) THEN
       FScript = ERF(SqrtR)
    ELSE
       XSUM = SFAC(1)*EXP(-R)
       DO LS = 2,L
          XSUM  = XSUM+(Sfac(LS)*R*EXP(-R/DBLE(LS-1)))**(LS-1)
       ENDDO 
       FScript = ERF(SqrtR)-(SqrtR/SqrtPi)*XSUM
    ENDIF
!
  END FUNCTION FScript
!========================================================================================
! Rieman Zeta Function
!========================================================================================
  FUNCTION RZeta(N,M)
    INTEGER                    :: N,M,I
    REAL(DOUBLE)               :: RZeta,RSum
    REAL(DOUBLE),DIMENSION(56) :: RZ = (/ 0.00000000000000000D0, 1.64493406684822644D0, & 
                   1.20205690315959429D0, 1.08232323371113819D0, 1.03692775514336993D0, & 
                   1.01734306198444914D0, 1.00834927738192283D0, 1.00407735619794434D0, &
                   1.00200839282608221D0, 1.00099457512781809D0, 1.00049418860411946D0, &
                   1.00024608655330805D0, 1.00012271334757849D0, 1.00006124813505870D0, &
                   1.00003058823630702D0, 1.00001528225940865D0, 1.00000763719763790D0, &
                   1.00000381729326500D0, 1.00000190821271655D0, 1.00000095396203387D0, &
                   1.00000047693298679D0, 1.00000023845050273D0, 1.00000011921992597D0, &
                   1.00000005960818905D0, 1.00000002980350351D0, 1.00000001490155483D0, &
                   1.00000000745071179D0, 1.00000000372533402D0, 1.00000000186265972D0, &
                   1.00000000093132743D0, 1.00000000046566291D0, 1.00000000023283118D0, &
                   1.00000000011641550D0, 1.00000000005820772D0, 1.00000000002910385D0, &
                   1.00000000001455192D0, 1.00000000000727596D0, 1.00000000000363798D0, &
                   1.00000000000181899D0, 1.00000000000090949D0, 1.00000000000045475D0, &
                   1.00000000000022737D0, 1.00000000000011369D0, 1.00000000000005684D0, & 
                   1.00000000000002842D0, 1.00000000000001421D0, 1.00000000000000711D0, &
                   1.00000000000000355D0, 1.00000000000000178D0, 1.00000000000000089D0, &
                   1.00000000000000044D0, 1.00000000000000022D0, 1.00000000000000011D0, &
                   1.00000000000000006D0, 1.00000000000000003D0, 1.00000000000000001D0  /)

!
    IF(N .LE. 56) THEN
       RZeta = RZ(N)
    ELSE
       RZeta = One
    ENDIF
!
    DO I=2,M
       RSum = RSum + One/(DBLE(I)**N)
    ENDDO
    RZeta = RZ(N) - RSum
!
  END FUNCTION RZeta
#endif
!
!
!
END MODULE PBCFarField

