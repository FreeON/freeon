MODULE PBCFarField
#ifdef PERIODIC
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
  USE PFFT
  IMPLICIT NONE
!====================================================================================
! Global Array intermediates for multipole opperations
!====================================================================================
  REAL(DOUBLE)                      :: DFac,QFac,Quadripole
  REAL(DOUBLE), DIMENSION(3)        :: CellCenter,Dipole
  REAL(DOUBLE), DIMENSION(0:FFLen)  :: RhoC,RhoS
  CONTAINS
!====================================================================================
!   Setup the PBCFarField Matrix. 
!====================================================================================
    SUBROUTINE PBCFarFieldSetUp(FFL,Q)
      INTEGER                         :: FFL,I,J,K,NC,LM
      REAL(DOUBLE),DIMENSION(3)       :: PQ
      REAL(DOUBLE)                    :: E_PFF,E_DP,E_QP,MACDist,PACDist,RMIN
      REAL(DOUBLE), DIMENSION(0:FFLen):: FarFC,FarFS
      TYPE(PoleNode)                  :: Q
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
      CellCenter = Zero
      DO I=1,3
         IF(GM%AutoW(I)) THEN
            DO J=1,3
               IF(GM%AutoW(J)) THEN
                  CellCenter(I) =  CellCenter(I) + Half*GM%BoxShape%D(I,J)
               ENDIF
            ENDDO
         ELSE
            CellCenter(I) = Q%Box%Center(I)
         ENDIF
      ENDDO

!     Calculate the Box Moments
!
      CALL RhoToSP(RhoC,RhoS)
!
!     Calculate the Dipole and Quadripole Moments
!    
      CALL RhoToPoles(Dipole,Quadripole)
!
!     Calculate the Size of the Box Needed  for the Direct J and Generate the Cells for the Inner Box
!
      CALL BoxBounds(Q,Dimen,MACDist,PACDist,RMIN)
!
!     Calculate the PBC FarField Tensor
!
      CALL PFFTensor(RMIN) 
!
!     Calculate PFF  energy
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
      CALL OpenASCII(OutFile,Out)
      WRITE(Out,*) '==================================Periodic================================' 
      WRITE(Out,*) 'Volume = ',Volume
      WRITE(Out,*) 'DFac   = ',DFac
      WRITE(Out,*) 'QFac   = ',QFac
      WRITE(Out,*) 'CellCenterX = ',CellCenter(1)
      WRITE(Out,*) 'CellCenterY = ',CellCenter(2)
      WRITE(Out,*) 'CellCenterZ = ',CellCenter(3)
      WRITE(Out,*) 'Pac Dist    = ',PACDist
      WRITE(Out,*) 'Mac Dist    = ',MACDist
      WRITE(Out,*) 'Num Cells   = ',CSMM1%NCells
      WRITE(Out,*) 'D(1)        = ',Dipole(1)
      WRITE(Out,*) 'D(2)        = ',Dipole(2)
      WRITE(Out,*) 'D(3)        = ',Dipole(3)
      WRITE(Out,*) 'Quadripole  = ',Quadripole
      WRITE(Out,*) ' '
      WRITE(Out,*) 'Correction to the Energy: PFF        = ',E_PFF
      WRITE(Out,*) 'Correction to the Energy: Dipole     = ',E_DP
      WRITE(Out,*) 'Correction to the Energy: Quadripole = ',E_QP  
      WRITE(Out,*) '=========================================================================='
      CLOSE(Out)
!
      WRITE(*,*) '==================================Periodic================================'  
      WRITE(*,*) 'Volume = ',Volume
      WRITE(*,*) 'DFac   = ',DFac
      WRITE(*,*) 'QFac   = ',QFac
      WRITE(*,*) 'CellCenterX = ',CellCenter(1)
      WRITE(*,*) 'CellCenterY = ',CellCenter(2)
      WRITE(*,*) 'CellCenterZ = ',CellCenter(3)
      WRITE(*,*) 'Pac Dist    = ',PACDist
      WRITE(*,*) 'Mac Dist    = ',MACDist
      WRITE(*,*) 'Num Cells   = ',CSMM1%NCells
      WRITE(*,*) 'D(1)        = ',Dipole(1)
      WRITE(*,*) 'D(2)        = ',Dipole(2)
      WRITE(*,*) 'D(3)        = ',Dipole(3)
      WRITE(*,*) 'Quadripole  = ',Quadripole
      WRITE(*,*) ' '
      WRITE(*,*) 'Correction to the Energy: PFF        = ',E_PFF
      WRITE(*,*) 'Correction to the Energy: Dipole     = ',E_DP
      WRITE(*,*) 'Correction to the Energy: Quadripole = ',E_QP  
      WRITE(*,*) '=========================================================================='
!
    END SUBROUTINE PBCFarFieldSetUp
!====================================================================================
!   Calculate the FarField Component of the J matrix
!====================================================================================
   FUNCTION CTraxFF(Prim,HGBra)
      TYPE(PrimPair)                   :: Prim
      INTEGER                          :: LM
      REAL(DOUBLE)                     :: PiZ,CTraxFF
      REAL(DOUBLE), DIMENSION(3)       :: PQ,HGDipole
      REAL(DOUBLE), DIMENSION(1:)      :: HGBra
      REAL(DOUBLE), DIMENSION(0:FFLen) :: SPBraC,SPBraS,SPKetC,SPKetS,FarFC,FarFS
!
      CTraxFF = Zero
!
!     Transform <Bra| coefficients from HG to SP
!
      PiZ=(Pi/Prim%Zeta)**(ThreeHalves)
      CALL HGToSP_Gen(Prim%Ell,PiZ,HGBra,SPBraC,SPBraS)   
      PQ = -Prim%P+CellCenter
!
      IF(.NOT. NoTranslate(PQ)) THEN
         SPKetC = Zero
         SPKetS = Zero
         FarFC  = Zero
         FarFS  = Zero
         CALL Regular(FFELL,PQ(1),PQ(2),PQ(3))
         CALL XLate77(FFELL,FFEll,SPKetC,SPKetS,Cpq,Spq,RhoC,RhoS)
         CALL CTraX77(Prim%ELL,FFEll,FarFC,FarFS,TensorC,TensorS,SPKetC,SPKetS)
      ELSE
         FarFC  = Zero
         FarFS  = Zero
         CALL CTraX77(Prim%ELL,FFEll,FarFC,FarFS,TensorC,TensorS,RhoC,RhoS)
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
      HGDipole = CalculateDiPole(Prim%Ell,PiZ,-PQ(1),-PQ(2),-PQ(3),HGBra(1:))
      CTraxFF  = CTraxFF + DFAC*(Dipole(1)*HGDipole(1)+Dipole(2)*HGDipole(2)+Dipole(3)*HGDipole(3))  
!
!     INCLUDE THE QUADRIPOLE correction 
!
      CTraxFF = CTraxFF + PiZ*HGBra(1)*QFAC*Quadripole
!
    END FUNCTION CTraxFF
!========================================================================================
! Calculate Quadripole = Q(2,0,0)+Q(0,2,0)+Q(0,0,2) of the local density
!========================================================================================
    SUBROUTINE RhoToPoles(Dip,Q2)
      INTEGER                   :: zq,iq,iadd,jadd,NQ,OffQ,OffR,LQ,LenQ
      REAL(DOUBLE)              :: RX,RY,RZ,R2,Expt,PiExpt
      REAL(DOUBLE)              :: Q2
      REAL(DOUBLE),DIMENSION(3) :: Dip
!
      Q2  = Zero
      Dip = Zero
      DO zq=1,Rho%NExpt
         NQ     = Rho%NQ%I(zq)
         Expt   = Rho%Expt%D(zq)
         PiExpt = (Pi/Expt)**(threehalves)
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
!
               Dip = Dip + CalculateDiPole(LQ,PiExpt,RX,RY,RZ,Rho%Co%D(jadd:jadd+LenQ+1))
               Q2  = Q2  + CalculateQuPole(LQ,Expt,PiExpt,RX,RY,RZ,Rho%Co%D(jadd:jadd+LenQ+1)) 
!
            ENDDO
         ENDIF
      ENDDO
!
    END SUBROUTINE RhoToPoles
!========================================================================================
!
!========================================================================================
    FUNCTION CalculateDiPole(LQ,PiExpt,RX,RY,RZ,Coef)
      INTEGER                   :: LQ
      REAL(DOUBLE)              :: RX,RY,RZ,PiExpt
      REAL(DOUBLE),DIMENSION(3) :: CalculateDiPole
      REAL(DOUBLE),DIMENSION(:) :: Coef
!
      SELECT CASE(LQ)
      CASE (0)
         CalculateDiPole(1) = -PiExpt*Coef(1)*RX
         CalculateDiPole(2) = -PiExpt*Coef(1)*RY
         CalculateDiPole(3) = -PiExpt*Coef(1)*RZ
      CASE(1:)
         CalculateDiPole(1) = -PiExpt*(Coef(2)+Coef(1)*RX)
         CalculateDiPole(2) = -PiExpt*(Coef(3)+Coef(1)*RY)
         CalculateDiPole(3) = -PiExpt*(Coef(4)+Coef(1)*RZ)
      END SELECT
!
    END FUNCTION CalculateDiPole
!========================================================================================
!
!========================================================================================
    FUNCTION CalculateQuPole(LQ,Expt,PiExpt,RX,RY,RZ,Coef)
      INTEGER                   :: LQ
      REAL(DOUBLE)              :: R2,RX,RY,RZ,PiExpt,Expt
      REAL(DOUBLE)              :: CalculateQuPole
      REAL(DOUBLE),DIMENSION(:) :: Coef
!
      R2 = Half/Expt + RX*RX + RY*RY + RZ*RZ 
      SELECT CASE(LQ)
      CASE (0)
         CalculateQuPole = PiExpt*Coef(1)*R2
      CASE(1)
         CalculateQuPole = PiExpt*(Coef(1)*R2 +                                         &
                                   Two*Coef(2)*RX + Two*Coef(3)*RY + Two*Coef(4)*RZ)
      CASE(2:)
         CalculateQuPole = PiExpt*(Coef(1)*R2 +                                         &
                                   Two*Coef(2)*RX + Two*Coef(3)*RY + Two*Coef(4)*RZ +   &
                                   Two*Coef(5)    + Two*Coef(7)    + Two*Coef(10))
      END SELECT
!  
    END FUNCTION CalculateQuPole
!========================================================================================
! Calculate the SP Moments of Rho_Loc
!========================================================================================
  SUBROUTINE RhoToSP(Cp,Sp)
      INTEGER                         :: LP,LQ,LPQ,LenP,LenQ,LenPQ,LKet,I
      INTEGER                         :: zq,iq,iadd,jadd,NQ,OffQ,OffR,LM
      REAL(DOUBLE)                    :: Zeta,PiZ,SUM1,SUM2
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
         LPQ   = MAX(LP,LQ)
         LenQ  = LSP(LQ)  
         LenP  = LSP(LP)  
         LenPQ = LSP(LPQ) 
         LKet  = LHGTF(LQ)
!
         IF(NQ > 0) THEN
            DO iq=1,NQ
               iadd = Rho%OffQ%I(zq)+iq
               jadd = Rho%OffR%I(zq)+(iq-1)*LKet+1
               PQ(1)  = Rho%Qx%D(iadd)-CellCenter(1)
               PQ(2)  = Rho%Qy%D(iadd)-CellCenter(2)
               PQ(3)  = Rho%Qz%D(iadd)-CellCenter(3)
               PiZ=(Pi/Zeta)**(ThreeHalves)
               CALL HGToSP_Gen(LQ,PiZ,Rho%Co%D(jadd:jadd+LKet-1),Cq(0:SPLen),Sq(0:SPLen)) 
               IF(NoTranslate(PQ))THEN
                  DO LM = 0,LenQ
                     Cp(LM) =  Cp(LM) + Cq(LM)   
                     Sp(LM) =  Sp(LM) + Sq(LM)
                  ENDDO
               ELSE
                  CALL Regular(LPQ,PQ(1),PQ(2),PQ(3))
                  CALL XLate77(LP,LQ,Cp,Sp,Cpq,Spq,Cq,Sq)
               ENDIF
            ENDDO
         ENDIF
      ENDDO
!
    END SUBROUTINE RhoToSP
!========================================================================================
! Calculate the Box Bounds Needed for the Direct Sum
!========================================================================================
  SUBROUTINE BoxBounds(Q,Dimen,MACDist,PACDist,Radius) 
    INTEGER                          :: I,J,K,LP,Dimen
    INTEGER                          :: IRmin,IRmax
    TYPE(PoleNode)                   :: Q
    REAL(DOUBLE)                     :: Px,Py,Pz,O_LP,O_FFEll,NFac,Dist,MACDist,PACDist
    REAL(DOUBLE)                     :: Radd,Radius,A0,B0,C0
!
    IF(Dimen==0) THEN
       IRmin = 0
       IRMax = 2
    ELSEIF(Dimen==1) THEN
       IRmin = 1
       IRMax = 50
    ELSEIF(Dimen==2) THEN
       IRmin = 7
       IRMax = 500
    ELSEIF(Dimen==3) THEN
       IRmin = 26
       IRMax = 5000
    ENDIF
!
!   PAC Distance (From PoleRoot)
!
    Px=Half*(Q%Box%BndBox(1,2)-Q%Box%BndBox(1,1))
    Py=Half*(Q%Box%BndBox(2,2)-Q%Box%BndBox(2,1))   
    Pz=Half*(Q%Box%BndBox(3,2)-Q%Box%BndBox(3,1))
    PACDist = SQRT(Px*Px+Py*Py+Pz*Pz)
!
!   MAC DISTANCE
!
    MACDist = Zero
    DO LP=0,BS%NASym 
       O_LP     = UnsoldO(LP,RhoC,RhoS)/(Thresholds%TwoE)
       O_FFELL  = UnsoldO(FFEll,RhoC,RhoS)
       NFac     = FudgeFactorial(LP,FFELL)
       Dist     = (O_LP*NFac*O_FFEll)**(One/DBLE(FFEll+LP+2))
       MACDist = MAX(MACDist,Dist)
    ENDDO
    MACDist = MACDist
!
!   Generate the Cells for the Inner Box
!  
    Radius = PACDist+MACDist
    CALL New_CellSet_Sphere(CSMM1,GM%AutoW,GM%BoxShape%D,Radius)
    DO 
       IF(CSMM1%NCells .LT. IRmax) THEN
          EXIT
       ELSE
          Radius = 0.99*Radius
          CALL Delete_CellSet(CSMM1)
          CALL New_CellSet_Sphere(CSMM1,GM%AutoW,GM%BoxShape%D,Radius)
       ENDIF
    ENDDO
!
    DO 
       IF(CSMM1%NCells .GT. IRmin) THEN
          EXIT
       ELSE
          Radius = 1.01*Radius
          CALL Delete_CellSet(CSMM1)
          CALL New_CellSet_Sphere(CSMM1,GM%AutoW,GM%BoxShape%D,Radius)
       ENDIF
    ENDDO
!
  END SUBROUTINE BoxBounds
!========================================================================================
! If QP is < TOL, do not translate
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
       WRITE(*,100) (OCC(J),J=1,2*JMAX+1)
 100   FORMAT(20I3)
    ENDDO
!
  END SUBROUTINE Print_Occ
#endif
END MODULE PBCFarField

