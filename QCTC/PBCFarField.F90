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
  USE PFFTen
  IMPLICIT NONE
!  
  INTEGER                             :: MaxELL
  REAL(DOUBLE)                        :: E_PFF,E_DP
  REAL(DOUBLE)                        :: MACDist,PACDist,BDist,RDist
  REAL(DOUBLE),DIMENSION(3)           :: BOXDist
! 
  CONTAINS
!====================================================================================
!   Setup the PBCFarField Matrix. 
!====================================================================================
    SUBROUTINE PBCFarFieldSetUp(Q,GMLoc)
      TYPE(PoleNode)                  :: Q
      TYPE(CRDS)                      :: GMLoc
      INTEGER                         :: I,J,K,NC,L,LM,LP
      REAL(DOUBLE)                    :: Radius,OL,NL,FAC,Px,Py,Pz
      REAL(DOUBLE),DIMENSION(3)       :: PQ
!
!     Get the new Tensor and MaxEll
!
      CALL Get(MaxEll,'MaxEll'//CurGeom)
      CALL New(TensorC,LSP(2*MaxEll),0)
      CALL New(TensorS,LSP(2*MaxEll),0) 
      CALL Get(TensorC,'PFFTensorC'//CurGeom)
      CALL Get(TensorS,'PFFTensorS'//CurGeom)
!
!     Get CS_IN and the Inner Cell Radius
!
      CALL Get(Radius,'CS_IN%Radius'//CurBase//CurGeom)
      CALL Get_CellSet(CS_IN,'CS_IN'//CurBase//CurGeom)
!
!     Calculate the Box Moments
!
      CALL RhoToSP(GMLoc)
!
!     Calculate BDist
!
      BDist = SQRT(BOXDist(1)**2+BOXDist(2)**2+BOXDist(3)**2)

      WRITE(*,*) 
!
!     PACDist (From PoleRoot)
!
      Px=Half*(Q%Box%BndBox(1,2)-Q%Box%BndBox(1,1))
      Py=Half*(Q%Box%BndBox(2,2)-Q%Box%BndBox(2,1))   
      Pz=Half*(Q%Box%BndBox(3,2)-Q%Box%BndBox(3,1))
      PACDist = SQRT(Px*Px+Py*Py+Pz*Pz)
!
!     If Not Over Riden Calculate MaxEll (Still Needs Work)
!
      IF(.NOT. GMLoc%PBC%PFFOvRide) THEN
         LP=BS%NAsym
         DO L = GMLoc%PBC%PFFMaxEll,FFELL
            OL = Unsold2(L/2,L,RhoC,RhoS)
            NL = FudgeFactorial(LP,L)
            FAC = OL*NL/((Radius-BDist)**(DBLE(L)+1))
            IF(FAC .LT. TauMAC) THEN
               MaxEll = L
               EXIT
            ENDIF
         ENDDO
      ENDIF
      IF(MaxEll == FFELL) THEN
         WRITE(*,*) ' *** WARNING *** MaxEll > FFEll *** WARNING *** '
      ENDIF
      RDist = SQRT(BDist**2+Radius**2)
!
!     Calculate PFF  energy
!
      E_PFF  = Zero
      FarFC  = Zero
      FarFS  = Zero
      CALL CTraX77(MaxELL,MaxELL,PFFBraC,PFFBraS,TensorC%D,TensorS%D,RhoC,RhoS)
      DO LM = 0,LSP(MaxELL)
         E_PFF = E_PFF + RhoC(LM)*PFFBraC(LM)+RhoS(LM)*PFFBraS(LM)
      ENDDO
      E_PFF = Two*E_PFF
!
!     Calculate Dipole energy
!
      IF(GMLoc%PBC%Dimen == 2) THEN
         IF(.NOT. GMLoc%PBC%AutoW(1) ) E_DP = Two*GMLoc%PBC%DipoleFAC*RhoPoles%DPole%D(1)**2
         IF(.NOT. GMLoc%PBC%AutoW(2) ) E_DP = Two*GMLoc%PBC%DipoleFAC*RhoPoles%DPole%D(2)**2         
         IF(.NOT. GMLoc%PBC%AutoW(3) ) E_DP = Two*GMLoc%PBC%DipoleFAC*RhoPoles%DPole%D(3)**2
      ELSEIF(GMLoc%PBC%Dimen == 3) THEN 
         E_DP = Two*GMLoc%PBC%DipoleFAC*(RhoPoles%DPole%D(1)**2+RhoPoles%DPole%D(2)**2+RhoPoles%DPole%D(3)**2)
      ENDIF
!
!     Output
!
      WRITE(*,*) 'GM%PBC%Dimen  = ',GMLoc%PBC%Dimen 
      WRITE(*,*) 'MaxEll        = ',MaxEll
      WRITE(*,*) 'MaxLay        = ',Radius/MaxBoxDim(GMLoc)
      WRITE(*,*) 'CS_IN%NCells  = ',CS_IN%NCells
      WRITE(*,*) 'CS_OUT%NCells = ',CS_OUT%NCells
      WRITE(*,*) 'PACDist       = ',PACDist
      WRITE(*,*) 'BOXDist       = ',BDist
      WRITE(*,*) 'RDist         = ',RDist      
      WRITE(*,*) '|Dipole|      = ',SQRT(RhoPoles%DPole%D(1)**2+RhoPoles%DPole%D(2)**2+RhoPoles%DPole%D(3)**2)
      WRITE(*,*)
      WRITE(*,*) 'Epsilon       = ',GMLoc%PBC%Epsilon
      WRITE(*,*) 'DipoleFAC     = ',GMLoc%PBC%DipoleFAC
      WRITE(*,*) 'E_PFF         = ',E_PFF
      WRITE(*,*) 'E_DP          = ',E_DP     
      WRITE(*,*)
!
!!$!
!!$!     Calculate the Size of the Box Needed  for the Direct J and Generate the Cells for the Inner Box
!!$!
!!$      CALL BoxBounds(GMLoc)
!!$!
!!$!     Calculate the MACDist, PACDist and BDist
!!$!
!!$      CALL CalMPB(Q,GMLoc)
!!$!
!!$!     Calculate the PBC FarField Tensor
!!$!
!!$      CALL PFFTensor(MaxELL+LP,GMLoc,Args) 
!
    END SUBROUTINE PBCFarFieldSetUp
!---------------------------------------------------------------------------------------------- 
!   Print out Periodic Info
!----------------------------------------------------------------------------------------------  
    SUBROUTINE Print_Periodic(Unit_O,GMLoc)
      INTEGER,OPTIONAL    :: Unit_O
      INTEGER             :: Unit
      TYPE(CRDS)          :: GMLoc
!
      IF(PRESENT(Unit_O)) THEN
         Unit=Unit_O
      ELSE
         Unit=Out
      ENDIF
!
      IF(PrintFlags%Key>DEBUG_MINIMUM) THEN
         IF(GMLoc%PBC%Dimen > 0) THEN
            CALL OpenASCII(OutFile,Unit)  
            WRITE(Unit,100)
            WRITE(Unit,101) MaxELL,GMLoc%PBC%Dimen
            WRITE(Unit,102) RhoPoles%DPole%D(1),RhoPoles%DPole%D(2),RhoPoles%DPole%D(3)
            WRITE(Unit,103) PACDIst,BDist,RDist
            WRITE(Unit,108) SQRT(RDist**2-BDist**2)/MaxBoxDim(GMLoc)
            WRITE(Unit,104) CS_OUT%NCells
            WRITE(Unit,110) CS_IN%NCells
            WRITE(Unit,105)
            WRITE(Unit,106) E_PFF,E_DP
            WRITE(Unit,107)       
            CLOSE(Unit)
         ENDIF
      ENDIF
100   FORMAT('========================================Periodic Information======================================')
101   FORMAT(' MaxEll             = ',I3,'     Dimension = ',I2)
102   FORMAT(' Dipole Moment      = (',F12.6,',',F12.6,',',F12.6,')')
103   FORMAT(' PAC Distance  = ',F6.2,' BOX Distance = ',F6.2,' Total Distance = ',F6.2)
108   FORMAT(' No. of Layers = ',F6.2)
104   FORMAT(' Outer No. of Cells = ',I4)
110   FORMAT(' Inner No. of Cells = ',I4)
105   FORMAT(' Correction to the Energy:')
106   FORMAT('   PFF = ',E14.6,'  Dipole = ',E14.6)
107   FORMAT('=========================================END======================================================')
!
    END SUBROUTINE Print_Periodic
!====================================================================================
!   Calculate the FarField Component of the J matrix
!====================================================================================
   FUNCTION CTraxFF(Prim,HGBra,GMLoc)
      TYPE(PrimPair)                   :: Prim
      INTEGER                          :: LM
      REAL(DOUBLE)                     :: PiZ,CTraxFF
      REAL(DOUBLE), DIMENSION(3)       :: PQ,HGDipole
      REAL(DOUBLE), DIMENSION(1:)      :: HGBra
      TYPE(CRDS)                       :: GMLoc
!
      CTraxFF = Zero
!
!     Transform <Bra| coefficients from HG to SP
!
      PiZ=(Pi/Prim%Zeta)**(ThreeHalves)
      CALL HGToSP_Gen(Prim%Ell,PiZ,HGBra,PFFBraC,PFFBraS)   
      PQ = -Prim%P+GMLoc%PBC%CellCenter
!
      IF(.NOT. NoTranslate(PQ)) THEN
         PFFKetC = Zero
         PFFKetS = Zero
         FarFC  = Zero
         FarFS  = Zero
         CALL Regular(MaxELL,PQ(1),PQ(2),PQ(3))
         CALL XLate77(MaxELL,MaxELL,PFFKetC,PFFKetS,Cpq,Spq,RhoC,RhoS)
         CALL CTraX77(Prim%ELL,MaxELL,FarFC,FarFS,TensorC%D,TensorS%D,PFFKetC,PFFKetS)
      ELSE
         FarFC  = Zero
         FarFS  = Zero
         CALL CTraX77(Prim%ELL,MaxELL,FarFC,FarFS,TensorC%D,TensorS%D,RhoC,RhoS)
      ENDIF
!
!     Contrax the <Bra|Ket> FF corection
!
      DO LM = 0,LSP(Prim%Ell)
         CTraxFF = CTraxFF + PFFBraC(LM)*FarFC(LM)+PFFBraS(LM)*FarFS(LM)
      ENDDO
!
!     Include the Dipole correction to FarFC and FarFS
!
      IF(GMLoc%PBC%Dimen == 1) THEN
         RETURN
      ELSEIF(GMLoc%PBC%Dimen == 2) THEN
         HGDipole = CalculateDiPole(Prim%Ell,Prim%Zeta,-PQ(1),-PQ(2),-PQ(3),HGBra(1:))
         IF(.NOT. GMLoc%PBC%AutoW(1)) CTraxFF  = CTraxFF + GMLoc%PBC%DipoleFAC*HGDipole(1)*RhoPoles%DPole%D(1)         
         IF(.NOT. GMLoc%PBC%AutoW(2)) CTraxFF  = CTraxFF + GMLoc%PBC%DipoleFAC*HGDipole(2)*RhoPoles%DPole%D(2)
         IF(.NOT. GMLoc%PBC%AutoW(3)) CTraxFF  = CTraxFF + GMLoc%PBC%DipoleFAC*HGDipole(3)*RhoPoles%DPole%D(3)
         RETURN
      ELSEIF(GMLoc%PBC%Dimen == 3) THEN
         HGDipole = CalculateDiPole(Prim%Ell,Prim%Zeta,-PQ(1),-PQ(2),-PQ(3),HGBra(1:))
         CTraxFF  = CTraxFF + GMLoc%PBC%DipoleFAC*(HGDipole(1)*RhoPoles%DPole%D(1) &
                                              + HGDipole(2)*RhoPoles%DPole%D(2) &
                                              + HGDipole(3)*RhoPoles%DPole%D(3) )
         RETURN
      ENDIF
!
    END FUNCTION CTraxFF
!====================================================================================
!   Calculate the FarField Component of the JForce matrix
!====================================================================================
   FUNCTION CTraxFF_Grad(Prim,HGBra,GMLoc)
      TYPE(PrimPair)                   :: Prim
      INTEGER                          :: LM
      REAL(DOUBLE)                     :: PiZ,CTraxFF_Grad
      REAL(DOUBLE), DIMENSION(3)       :: PQ,HGDipole
      REAL(DOUBLE), DIMENSION(1:)      :: HGBra
      TYPE(CRDS)                       :: GMLoc
!
      CTraxFF_Grad = Zero
!
!     Transform <Bra| coefficients from HG to SP
!
      PiZ=(Pi/Prim%Zeta)**(ThreeHalves)
      CALL HGToSP_Gen(Prim%Ell,PiZ,HGBra,PFFBraC,PFFBraS)   
      PQ = -Prim%P+GMLoc%PBC%CellCenter
!
      IF(.NOT. NoTranslate(PQ)) THEN 
         PFFKetC = Zero
         PFFKetS = Zero
         FarFC  = Zero
         FarFS  = Zero
         CALL Regular(MaxELL,PQ(1),PQ(2),PQ(3))
         CALL XLate77(MaxELL,MaxELL,PFFKetC,PFFKetS,Cpq,Spq,RhoC,RhoS)
         CALL CTraX77(Prim%ELL,MaxELL,FarFC,FarFS,TensorC%D,TensorS%D,PFFKetC,PFFKetS)
      ELSE
         FarFC  = Zero
         FarFS  = Zero
         CALL CTraX77(Prim%ELL,MaxELL,FarFC,FarFS,TensorC%D,TensorS%D,RhoC,RhoS)
      ENDIF
!
!     Contrax the <Bra|Ket> FF corection
!
      DO LM = 0,LSP(Prim%Ell)
         CTraxFF_Grad = CTraxFF_Grad + PFFBraC(LM)*FarFC(LM)+PFFBraS(LM)*FarFS(LM)
      ENDDO
!
!     Include the Dipole and Quadripole Correction to FarFC and FarFS
!
      IF(GMLoc%PBC%Dimen == 1) THEN
         RETURN
      ELSEIF(GMLoc%PBC%Dimen == 2) THEN
         HGDipole = CalculateDiPole(Prim%Ell,Prim%Zeta,-PQ(1),-PQ(2),-PQ(3),HGBra(1:))
         IF(.NOT. GMLoc%PBC%AutoW(1)) CTraxFF_Grad  = CTraxFF_Grad + GMLoc%PBC%DipoleFAC*HGDipole(1)*RhoPoles%DPole%D(1)
         IF(.NOT. GMLoc%PBC%AutoW(2)) CTraxFF_Grad  = CTraxFF_Grad + GMLoc%PBC%DipoleFAC*HGDipole(2)*RhoPoles%DPole%D(2)
         IF(.NOT. GMLoc%PBC%AutoW(3)) CTraxFF_Grad  = CTraxFF_Grad + GMLoc%PBC%DipoleFAC*HGDipole(3)*RhoPoles%DPole%D(3)
         RETURN
      ELSEIF(GMLoc%PBC%Dimen == 3) THEN
         HGDipole = CalculateDiPole(Prim%Ell,Prim%Zeta,-PQ(1),-PQ(2),-PQ(3),HGBra(1:))
         CTraxFF_Grad  = CTraxFF_Grad + GMLoc%PBC%DipoleFAC*(HGDipole(1)*RhoPoles%DPole%D(1) &
                                                        + HGDipole(2)*RhoPoles%DPole%D(2) &
                                                        + HGDipole(3)*RhoPoles%DPole%D(3) )
         RETURN
      ENDIF
!
    END FUNCTION CTraxFF_Grad
!========================================================================================
! Calculate the SP Moments of Rho_Loc
!========================================================================================
  SUBROUTINE RhoToSP(GMLoc)
      INTEGER                         :: LP,LQ,LPQ,LenP,LenQ,LenPQ,LKet,I
      INTEGER                         :: zq,iq,iadd,jadd,NQ,OffQ,OffR,LM
      REAL(DOUBLE)                    :: Zeta,PiZ
      REAL(DOUBLE),DIMENSION(3)       :: PQ
      TYPE(CRDS)                      :: GMLoc
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
         LP    = MaxELL
         LPQ   = MAX(LP,LQ)
         LenQ  = LSP(LQ)  
         LenP  = LSP(LP)  
         LenPQ = LSP(LPQ) 
         LKet  = LHGTF(LQ)
!
         IF(NQ > 0) THEN
            DO iq=1,NQ
               iadd   = Rho%OffQ%I(zq)+iq
               jadd   = Rho%OffR%I(zq)+(iq-1)*LKet+1
               PQ(1)  = Rho%Qx%D(iadd)-GMLoc%PBC%CellCenter(1)
               PQ(2)  = Rho%Qy%D(iadd)-GMLoc%PBC%CellCenter(2)
               PQ(3)  = Rho%Qz%D(iadd)-GMLoc%PBC%CellCenter(3)
               PiZ=(Pi/Zeta)**(ThreeHalves)               
               BOXDist(1) = MAX(BOXDist(1),ABS(PQ(1)))
               BOXDist(2) = MAX(BOXDist(2),ABS(PQ(2)))
               BOXDist(3) = MAX(BOXDist(3),ABS(PQ(3)))
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
      RhoC(0) = Zero
      RhoS(0) = Zero
!
    END SUBROUTINE RhoToSP
!========================================================================================
! Calculate the Box Bounds Needed for the Direct Sum
!========================================================================================
  SUBROUTINE BoxBounds(GMLoc)
    INTEGER                          :: I
    INTEGER                          :: IRmin,IRmax
    REAL(DOUBLE)                     :: Radius
    TYPE(CRDS)                       :: GMLoc
!        
    IF(GMLoc%PBC%Dimen==0) THEN
       IRmin = 1
       IRmax = 1
    ELSEIF(GMLoc%PBC%Dimen==1) THEN
       IRmin = 3
       IRmax = MAX(  50,CS_OUT%NCells+1)
    ELSEIF(GMLoc%PBC%Dimen==2) THEN
       IRmin = 9
       IRmax = MAX( 500,CS_OUT%NCells+1)
    ELSEIF(GMLoc%PBC%Dimen==3) THEN
       IRmin = 27
       IRmax = MAX(5000,CS_OUT%NCells+1)
    ENDIF
!
!   Generate the Cells for the Inner Box
!  
    Radius = RDist
    CALL New_CellSet_Sphere(CS_IN,GMLoc%PBC%AutoW,GMLoc%PBC%BoxShape,Radius)
!
!   Rescale the box size  so that the box is within the bounds of IRmin and IRmax
!
    IF(CS_IN%NCells < IRmin) THEN
       DO I=1,1000
          Radius = 1.001D0*Radius
          CALL Delete_CellSet(CS_IN)
          CALL New_CellSet_Sphere(CS_IN,GMLoc%PBC%AutoW,GMLoc%PBC%BoxShape,Radius)
          IF(CS_IN%NCells > IRmin) THEN
             EXIT
          ENDIF
       ENDDO
       RDist = Radius
       RETURN
    ELSEIF(CS_IN%NCells > IRMax) THEN
       DO I=1,1000
          Radius = 0.999D0*Radius
          CALL Delete_CellSet(CS_IN)
          CALL New_CellSet_Sphere(CS_IN,GMLoc%PBC%AutoW,GMLoc%PBC%BoxShape,Radius)
          IF(CS_IN%NCells < IRmax) THEN
             EXIT
          ENDIF
       ENDDO
       RDist = Radius
       RETURN
    ENDIF
    RDist = Radius
!
  END SUBROUTINE BoxBounds
!========================================================================================
! Calculate the Box Bounds Needed for the Direct Sum
!========================================================================================
!========================================================================================
! Calculate the Box Bounds Needed for the Direct Sum
!========================================================================================
  SUBROUTINE CalMPB(Q,GMLoc) 
    TYPE(PoleNode)                   :: Q
    INTEGER                          :: I,J,K,LP
    REAL(DOUBLE)                     :: Px,Py,Pz,O_FFEll,NFac,Dist,Mx,My,Mz
    REAL(DOUBLE)                     :: Radius
    TYPE(CRDS)                       :: GMLoc
!
!   PAC Distance (From PoleRoot)
!
    Px=Half*(Q%Box%BndBox(1,2)-Q%Box%BndBox(1,1))
    Py=Half*(Q%Box%BndBox(2,2)-Q%Box%BndBox(2,1))   
    Pz=Half*(Q%Box%BndBox(3,2)-Q%Box%BndBox(3,1))
    PACDist = SQRT(Px*Px+Py*Py+Pz*Pz)
!
!   MAC Distance
!
    LP      = 0
    MACDist = Zero
    DO I=-1,1
       DO J=-1,1
          DO K=-1,1
             Mx = DBLE(I)*BOXDist(1)
             My = DBLE(J)*BOXDist(2)
             Mz = DBLE(K)*BOXDist(3)
             IF(.NOT. GMLoc%PBC%AutoW(1)) Mx = Zero
             IF(.NOT. GMLoc%PBC%AutoW(2)) My = Zero
             IF(.NOT. GMLoc%PBC%AutoW(3)) Mz = Zero
             CALL Regular(MaxELL,Mx,My,Mz)
             CALL XLate77(MaxELL,MaxELL,FarFC,FarFS,Cpq,Spq,RhoC,RhoS)
!
             O_FFEll  = Unsold2(MaxEll-2,MaxELL,FarFC,FarFS)
             NFac     = FudgeFactorial(LP,MaxELL)
             Dist     = (NFac*O_FFEll/TauMAC)**(One/DBLE(MaxELL+LP+1))
             MACDist  = MAX(MACDist,Dist)
!
          ENDDO
       ENDDO
    ENDDO
!
!   BOX Distance
!
    BDist = SQRT(BOXDist(1)**2+BOXDist(2)**2+BOXDist(3)**2)
!
!   Total Distance
!
    RDist = MAX(PACDist,SQRT(MACDIst**2+BDist**2))
!
  END SUBROUTINE CalMPB
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

