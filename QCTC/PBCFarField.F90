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
  USE AtomPairs
  USE MondoPoles
  USE PoleTree
  USE SpecFun
  USE Globals
  USE PFFT
  IMPLICIT NONE
!
  REAL(DOUBLE)                        :: E_PFF,E_DP,E_QP,QRPole
  REAL(DOUBLE)                        :: MACDist,PACDist,RDist
!
  CONTAINS
!====================================================================================
!   Setup the PBCFarField Matrix. 
!====================================================================================
    SUBROUTINE PBCFarFieldSetUp(Q)
      TYPE(PoleNode)                  :: Q
      INTEGER                         :: I,J,K,NC,LM
      REAL(DOUBLE),DIMENSION(3)       :: PQ
!
!     Calculate the Box Moments
!
      CALL RhoToSP()
!
!     Calculate the Size of the Box Needed  for the Direct J and Generate the Cells for the Inner Box
!
      CALL BoxBounds(Q,MACDist,PACDist,RDist)
!
!     Calculate the PBC FarField Tensor
!
      CALL PFFTensor(FFEll2,RDist) 
!
!     Calculate PFF  energy
!
      E_PFF  = Zero
      FarFC  = Zero
      FarFS  = Zero
      CALL CTraX77(FFELL,FFEll,PFFBraC,PFFBraS,TensorC,TensorS,RhoC,RhoS)
      DO LM = 0,LSP(FFEll)
         E_PFF = E_PFF + RhoC(LM)*PFFBraC(LM)+RhoS(LM)*PFFBraS(LM)
      ENDDO
      E_PFF = Two*E_PFF
!
!     Calculate Dipole energy
!
      E_DP = Two*GM%PBC%DipoleFAC*(RhoPoles%DPole%D(1)**2+RhoPoles%DPole%D(2)**2+RhoPoles%DPole%D(3)**2)
!
!     Calculate the Radial Quadrupole
!
      QRPole = RhoPoles%QPole%D(1)+RhoPoles%QPole%D(2)+RhoPoles%QPole%D(3)
!
!     Calculate Quadrupole energy
!
      E_QP = ABS(GM%NElec)*GM%PBC%QupoleFAC*QRPole
!
    END SUBROUTINE PBCFarFieldSetUp
!---------------------------------------------------------------------------------------------- 
!   Print out Periodic Info
!----------------------------------------------------------------------------------------------  
    SUBROUTINE Print_Periodic(Unit_O)
      INTEGER,OPTIONAL    :: Unit_O
      INTEGER             :: Unit
!
      IF(PRESENT(Unit_O)) THEN
         Unit = Unit_O
      ELSE
         Unit=Out
      ENDIF
!
      IF(PrintFlags%Key>DEBUG_MINIMUM) THEN
         IF(GM%PBC%Dimen > 0) THEN
            CALL OpenASCII(OutFile,Unit)  
            WRITE(Unit,100)
            WRITE(Unit,101) QRPole
            WRITE(Unit,102) RhoPoles%DPole%D(1),RhoPoles%DPole%D(2),RhoPoles%DPole%D(3)
            WRITE(Unit,103) PACDist,MACDist,RDist
            WRITE(Unit,104) CS_OUT%NCells
            WRITE(Unit,110) CS_IN%NCells
            WRITE(Unit,105)
            WRITE(Unit,106) E_PFF,E_DP,E_QP
            WRITE(Unit,107)       
            CLOSE(Unit)
         ENDIF
      ENDIF
100   FORMAT('================================Periodic Information==============================')
101   FORMAT(' Radial Quadrupole  = ',F12.6)
102   FORMAT(' Dipole Moment      = (',F12.6,',',F12.6,',',F12.6,')')
103   FORMAT(' PAC Distance       = ',F6.2,'  MAC Distance = ',F6.2,'  Total Distance = ',F6.2)
104   FORMAT(' Outer No. of Cells = ',I4)
110   FORMAT(' Inner No. of Cells = ',I4)
105   FORMAT(' Correction to the Energy:')
106   FORMAT('   PFF = ',E14.6,'  Dipole = ',E14.6,'  Quadrupole =',E14.6)
107   FORMAT('==================================================================================')
!
    END SUBROUTINE Print_Periodic
!====================================================================================
!   Calculate the FarField Component of the J matrix
!====================================================================================
   FUNCTION CTraxFF(Prim,HGBra)
      TYPE(PrimPair)                   :: Prim
      INTEGER                          :: LM
      REAL(DOUBLE)                     :: PiZ,CTraxFF
      REAL(DOUBLE), DIMENSION(3)       :: PQ,HGDipole
      REAL(DOUBLE), DIMENSION(1:)      :: HGBra
!
      CTraxFF = Zero
!
!     Transform <Bra| coefficients from HG to SP
!
      PiZ=(Pi/Prim%Zeta)**(ThreeHalves)
      CALL HGToSP_Gen(Prim%Ell,PiZ,HGBra,PFFBraC,PFFBraS)   
      PQ = -Prim%P+GM%PBC%CellCenter
!
      IF(.NOT. NoTranslate(PQ)) THEN
         PFFKetC = Zero
         PFFKetS = Zero
         FarFC  = Zero
         FarFS  = Zero
         CALL Regular(FFELL,PQ(1),PQ(2),PQ(3))
         CALL XLate77(FFELL,FFEll,PFFKetC,PFFKetS,Cpq,Spq,RhoC,RhoS)
         CALL CTraX77(Prim%ELL,FFEll,FarFC,FarFS,TensorC,TensorS,PFFKetC,PFFKetS)
      ELSE
         FarFC  = Zero
         FarFS  = Zero
         CALL CTraX77(Prim%ELL,FFEll,FarFC,FarFS,TensorC,TensorS,RhoC,RhoS)
      ENDIF
!
!     Contrax the <Bra|Ket> FF corection
!
      DO LM = 0,LSP(Prim%Ell)
         CTraxFF = CTraxFF + PFFBraC(LM)*FarFC(LM)+PFFBraS(LM)*FarFS(LM)
      ENDDO
!
!     INCLUDE THE DIPOLE correction to FarFC and FarFS
!
      HGDipole = CalculateDiPole(Prim%Ell,Prim%Zeta,-PQ(1),-PQ(2),-PQ(3),HGBra(1:))
      CTraxFF  = CTraxFF + GM%PBC%DipoleFAC*(HGDipole(1)*RhoPoles%DPole%D(1) &
                                           + HGDipole(2)*RhoPoles%DPole%D(2) &
                                           + HGDipole(3)*RhoPoles%DPole%D(3) )
!
!     INCLUDE THE QUADRIPOLE correction 
!
      CTraxFF = CTraxFF + PiZ*HGBra(1)*GM%PBC%QupoleFAC*QRPole
!
    END FUNCTION CTraxFF
!========================================================================================
! Calculate the SP Moments of Rho_Loc
!========================================================================================
  SUBROUTINE RhoToSP()
      INTEGER                         :: LP,LQ,LPQ,LenP,LenQ,LenPQ,LKet,I
      INTEGER                         :: zq,iq,iadd,jadd,NQ,OffQ,OffR,LM
      REAL(DOUBLE)                    :: Zeta,PiZ
      REAL(DOUBLE),DIMENSION(3)       :: PQ
!
      RhoC  = Zero
      RhoS  = Zero
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
               PQ(1)  = Rho%Qx%D(iadd)-GM%PBC%CellCenter(1)
               PQ(2)  = Rho%Qy%D(iadd)-GM%PBC%CellCenter(2)
               PQ(3)  = Rho%Qz%D(iadd)-GM%PBC%CellCenter(3)
               PiZ=(Pi/Zeta)**(ThreeHalves)
               CALL HGToSP_Gen(LQ,PiZ,Rho%Co%D(jadd:jadd+LKet-1),PFFKetC(0:SPLen),PFFKetS(0:SPLen)) 
               IF(NoTranslate(PQ))THEN
                  DO LM = 0,LenQ
                     RhoC(LM) =  RhoC(LM) + PFFKetC(LM)   
                     RhoS(LM) =  RhoS(LM) + PFFKetS(LM)
                  ENDDO
               ELSE
                  CALL Regular(LPQ,PQ(1),PQ(2),PQ(3))
                  CALL XLate77(LP,LQ,RhoC,RhoS,Cpq,Spq,PFFKetC,PFFKetS)
               ENDIF
            ENDDO
         ENDIF
      ENDDO
!
    END SUBROUTINE RhoToSP
!========================================================================================
! Calculate the Box Bounds Needed for the Direct Sum
!========================================================================================
  SUBROUTINE BoxBounds(Q,MACDist,PACDist,Radius) 
    TYPE(PoleNode)                   :: Q
    INTEGER                          :: I,J,K,LP
    INTEGER                          :: IRmin,IRmax
    REAL(DOUBLE)                     :: Px,Py,Pz,O_LP,O_FFEll,NFac,Dist
    REAL(DOUBLE)                     :: MACDist,PACDist,Radius
!
    IF(GM%PBC%Dimen==0) THEN
       IRmin = 0
       IRmax = 2
    ELSEIF(GM%PBC%Dimen==1) THEN
       IRmin = 1
       IRmax = 10
    ELSEIF(GM%PBC%Dimen==2) THEN
       IRmin = 7
       IRmax = 100
    ELSEIF(GM%PBC%Dimen==3) THEN
       IRmin = 26
       IRmax = 1000
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
    CALL New_CellSet_Sphere(CS_IN,GM%PBC%AutoW,GM%PBC%BoxShape,Radius)
    DO 
       IF(CS_IN%NCells .LT. IRmax) THEN
          EXIT
       ELSE
          Radius = 0.99D0*Radius
          CALL Delete_CellSet(CS_IN)
          CALL New_CellSet_Sphere(CS_IN,GM%PBC%AutoW,GM%PBC%BoxShape,Radius)
       ENDIF
    ENDDO
!
    DO 
       IF(CS_IN%NCells .GT. IRmin) THEN
          EXIT
       ELSE
          Radius = 1.01D0*Radius
          CALL Delete_CellSet(CS_IN)
          CALL New_CellSet_Sphere(CS_IN,GM%PBC%AutoW,GM%PBC%BoxShape,Radius)
       ENDIF
    ENDDO
!
  END SUBROUTINE BoxBounds
!========================================================================================
! If QP is < TOL, do not translate
!========================================================================================
  FUNCTION NoTranslate(X) 
    REAL(DOUBLE),DIMENSION(3) :: X
    REAL(DOUBLE),PARAMETER    :: TOL=1.0D-12
    LOGICAL                   :: NoTranslate
!
    NoTranslate = (ABS(X(1)).LT.TOL) .AND. (ABS(X(2)).LT.TOL) .AND. (ABS(X(3)).LT.TOL)
!
  END FUNCTION NoTranslate
!
#endif
END MODULE PBCFarField

