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
  REAL(DOUBLE)                        :: PDist,BDist,RDist
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
      REAL(DOUBLE)                    :: Layers,OL,NL,FAC,Px,Py,Pz
      REAL(DOUBLE),DIMENSION(3)       :: PQ
!
!     Get CS_IN and the Closest Cell Distance
!
      CALL Get(RDist,'CS_IN%Radius'//CurBase//CurGeom)
      CALL Get_CellSet(CS_IN,'CS_IN'//CurBase//CurGeom)
      IF(GMLoc%PBC%Dimen==0) RETURN
!
!     Get the new Tensor and MaxEll
!
      CALL Get(MaxEll,'MaxEll'//CurBase//CurGeom)
      CALL New(TensorC,LSP(2*MaxEll),0)
      CALL New(TensorS,LSP(2*MaxEll),0) 
      CALL Get(TensorC,'PFFTensorC'//CurBase//CurGeom)
      CALL Get(TensorS,'PFFTensorS'//CurBase//CurGeom)
!
!     Initialize Arrays
!
      CALL New(TenRhoC,LSP(MaxEll),0)
      CALL New(TenRhoS,LSP(MaxEll),0)
!
      CALL New(PFFBraC,LSP(MaxEll),0)
      CALL New(PFFBraS,LSP(MaxEll),0)
      CALL New(PFFKetC,LSP(MaxEll),0)
      CALL New(PFFKetS,LSP(MaxELL),0)
!
!     Initailize RhoC and RhoS
!
      CALL New(RhoC   ,LSP(FFEll),0)  
      CALL New(RhoS   ,LSP(FFEll),0) 
!
!     Calculate the Box Moments and the BDist
!
      CALL RhoToSP(GMLoc)
!
!     Contract TenRho 
!
      TenRhoC%D=Zero
      TenRhoS%D=Zero
      CALL CTraX77(MaxELL,MaxEll,TenRhoC%D,TenRhoS%D,TensorC%D,TensorS%D,RhoC%D,RhoS%D)  
!
!     PACDist (From PoleRoot)
!
      Px=Half*(Q%Box%BndBox(1,2)-Q%Box%BndBox(1,1))
      Py=Half*(Q%Box%BndBox(2,2)-Q%Box%BndBox(2,1))   
      Pz=Half*(Q%Box%BndBox(3,2)-Q%Box%BndBox(3,1))
      PDist = SQRT(Px*Px+Py*Py+Pz*Pz)
!
!     If Not Over Riden Calculate MaxEll 
!
      CALL CalMaxEll(GMLoc)
!
!     Calculate PFF  energy
!
      DO LM = 0,LSP(MaxELL)
         E_PFF = E_PFF + RhoC%D(LM)*TenRhoC%D(LM)+RhoS%D(LM)*TenRhoS%D(LM)
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
!!$      Layers = SQRT(CS_IN%CellCarts%D(1,1)**2+CS_IN%CellCarts%D(2,1)**2+CS_IN%CellCarts%D(3,1)**2)/MaxBoxDim(GMLoc)
!!$      WRITE(*,*)
!!$      WRITE(*,*) 'GM%PBC%Dimen  = ',GMLoc%PBC%Dimen 
!!$      WRITE(*,*) 'MaxEll        = ',MaxEll
!!$      WRITE(*,*) 'MaxLay        = ',Layers
!!$      WRITE(*,*) 'CS_IN%NCells  = ',CS_IN%NCells
!!$      WRITE(*,*) 'CS_OUT%NCells = ',CS_OUT%NCells
!!$      WRITE(*,*) 'PACDist       = ',PDist
!!$      WRITE(*,*) 'BOXDist       = ',BDist
!!$      WRITE(*,*) 'RDist         = ',RDist      
!!$      WRITE(*,*) '|Dipole|      = ',SQRT(RhoPoles%DPole%D(1)**2+RhoPoles%DPole%D(2)**2+RhoPoles%DPole%D(3)**2)
!!$      WRITE(*,*)
!!$      WRITE(*,*) 'Epsilon       = ',GMLoc%PBC%Epsilon
!!$      WRITE(*,*) 'DipoleFAC     = ',GMLoc%PBC%DipoleFAC
!!$      WRITE(*,*) 'E_PFF         = ',E_PFF
!!$      WRITE(*,*) 'E_DP          = ',E_DP     
!!$      WRITE(*,*)
!
    END SUBROUTINE PBCFarFieldSetUp
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
      CtraxFF   = Zero   
      PFFBraC%D = Zero
      PFFBraS%D = Zero
!
!     Transform <Bra| coefficients from HG to SP
!
      PiZ=(Pi/Prim%Zeta)**(ThreeHalves)
      CALL HGToSP_Gen(Prim%Ell,PiZ,HGBra,PFFKetC%D,PFFKetS%D)   
      PQ = Prim%P-GMLoc%PBC%CellCenter   
!
!     Contract
!
      IF(NoTranslate(PQ)) THEN
         DO LM = 0,LSP(Prim%Ell)
            CTraxFF = CTraxFF + PFFKetC%D(LM)*TenRhoC%D(LM)+PFFKetS%D(LM)*TenRhoS%D(LM)
         ENDDO
      ELSE
         CALL Regular(MaxELL,PQ(1),PQ(2),PQ(3))
         CALL XLate77(MaxELL,Prim%Ell,PFFBraC%D,PFFBraS%D,Cpq,Spq,PFFKetC%D,PFFKetS%D)
         DO LM = 0,LSP(MaxEll)
            CTraxFF = CTraxFF + PFFBraC%D(LM)*TenRhoC%D(LM)+PFFBraS%D(LM)*TenRhoS%D(LM)
         ENDDO
      ENDIF
!
!     Include the Dipole correction to FarFC and FarFS
!      
      PQ = Prim%P-GMLoc%PBC%CellCenter 
      IF(GMLoc%PBC%Dimen == 1) THEN
         RETURN
      ELSEIF(GMLoc%PBC%Dimen == 2) THEN
         HGDipole = CalculateDiPole(Prim%Ell,Prim%Zeta,PQ(1),PQ(2),PQ(3),HGBra(1:))
         IF(.NOT. GMLoc%PBC%AutoW(1)) CTraxFF  = CTraxFF + GMLoc%PBC%DipoleFAC*HGDipole(1)*RhoPoles%DPole%D(1)         
         IF(.NOT. GMLoc%PBC%AutoW(2)) CTraxFF  = CTraxFF + GMLoc%PBC%DipoleFAC*HGDipole(2)*RhoPoles%DPole%D(2)
         IF(.NOT. GMLoc%PBC%AutoW(3)) CTraxFF  = CTraxFF + GMLoc%PBC%DipoleFAC*HGDipole(3)*RhoPoles%DPole%D(3)
         RETURN
      ELSEIF(GMLoc%PBC%Dimen == 3) THEN
         HGDipole = CalculateDiPole(Prim%Ell,Prim%Zeta,PQ(1),PQ(2),PQ(3),HGBra(1:))
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
      CtraxFF_Grad = Zero   
      PFFBraC%D    = Zero
      PFFBraS%D    = Zero
!
!     Transform <Bra| coefficients from HG to SP
!
      PiZ=(Pi/Prim%Zeta)**(ThreeHalves)
      CALL HGToSP_Gen(Prim%Ell,PiZ,HGBra,PFFKetC%D,PFFKetS%D)   
      PQ = Prim%P-GMLoc%PBC%CellCenter   
!
!     Contract
!      
      IF(NoTranslate(PQ)) THEN
         DO LM = 0,LSP(Prim%Ell)
            CTraxFF_Grad = CTraxFF_Grad + PFFKetC%D(LM)*TenRhoC%D(LM)+PFFKetS%D(LM)*TenRhoS%D(LM)
         ENDDO
      ELSE
         CALL Regular(MaxELL,PQ(1),PQ(2),PQ(3))
         CALL XLate77(MaxELL,Prim%Ell,PFFBraC%D,PFFBraS%D,Cpq,Spq,PFFKetC%D,PFFKetS%D)
         DO LM = 0,LSP(MaxEll)
            CTraxFF_Grad = CTraxFF_Grad + PFFBraC%D(LM)*TenRhoC%D(LM)+PFFBraS%D(LM)*TenRhoS%D(LM)
         ENDDO
      ENDIF
!
!     Include the Dipole and Quadripole Correction to FarFC and FarFS
!
      PQ = Prim%P-GMLoc%PBC%CellCenter 
      IF(GMLoc%PBC%Dimen == 1) THEN
         RETURN
      ELSEIF(GMLoc%PBC%Dimen == 2) THEN
         HGDipole = CalculateDiPole(Prim%Ell,Prim%Zeta,PQ(1),PQ(2),PQ(3),HGBra(1:))
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
!---------------------------------------------------------------------------------------------- 
!   Print out Periodic Info
!----------------------------------------------------------------------------------------------  
    SUBROUTINE Print_Periodic(GMLoc,Prog,Unit_O)
      CHARACTER(LEN=*)    :: Prog
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
      INTEGER,OPTIONAL    :: Unit_O
      INTEGER             :: Unit
      TYPE(CRDS)          :: GMLoc
      REAL(DOUBLE)        :: Layers
!      IF(GMLoc%PBC%Dimen==0.OR.PrintFlags%Key<=DEBUG_MINIMUM)RETURN

#ifdef PARALLEL
  IF(MyID == ROOT) THEN
#endif
      Unit=OpenPU(Unit_O=Unit_O)
      Mssg=ProcessName(Prog,TRIM(IntToChar(GMLoc%PBC%Dimen))//'-D periodics') &
         //'MxL = '//TRIM(IntToChar(MaxEll))                                  &
         //', PAC Cells = '//TRIM(IntToChar(CS_OUT%NCells))                   &
         //', MAC Cells = '//TRIM(IntToChar(CS_IN%NCells))                    
      WRITE(Unit,*)TRIM(Mssg)
      Mssg=ProcessName(Prog,'FF Energy')//'<FF> = '//TRIM(DblToChar(E_PFF))
      WRITE(Unit,*)TRIM(Mssg)

#ifdef LDJFLDJFLDJF
!
      IF(PRESENT(Unit_O)) THEN
         Unit=Unit_O
      ELSE
         Unit=Out
      ENDIF
!     
      Layers = SQRT(CS_IN%CellCarts%D(1,1)**2+CS_IN%CellCarts%D(2,1)**2+CS_IN%CellCarts%D(3,1)**2)/MaxBoxDim(GMLoc)



         IF(GMLoc%PBC%Dimen > 0) THEN
            CALL OpenASCII(OutFile,Unit)  
            WRITE(Unit,100)
            WRITE(Unit,101) MaxELL,GMLoc%PBC%Dimen
            WRITE(Unit,102) RhoPoles%DPole%D(1),RhoPoles%DPole%D(2),RhoPoles%DPole%D(3)
            WRITE(Unit,103) PDIst,BDist,RDist
            WRITE(Unit,108) Layers
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
103   FORMAT(' PAC Distance       = ',F6.2,' BOX Distance = ',F6.2,' Cell Distance = ',F6.2)
108   FORMAT(' No. of Layers      = ',F6.2)
104   FORMAT(' Outer No. of Cells = ',I4)
110   FORMAT(' Inner No. of Cells = ',I4)
105   FORMAT(' Correction to the Energy:')
106   FORMAT('   PFF = ',E14.6,'  Dipole = ',E14.6)
107   FORMAT('=========================================END======================================================')
#endif

      CALL ClosePU(Unit)
#ifdef PARALLEL
  ENDIF
#endif

!
    END SUBROUTINE Print_Periodic
!========================================================================================
! Calculate the SP Moments of Rho_Loc
!========================================================================================
  SUBROUTINE RhoToSP(GMLoc)
      INTEGER                         :: LP,LQ,LPQ,LenP,LenQ,LenPQ,LKet,I
      INTEGER                         :: zq,iq,iadd,jadd,NQ,OffQ,OffR,LM
      REAL(DOUBLE)                    :: Zeta,PiZ,Dist
      REAL(DOUBLE),DIMENSION(3)       :: PQ
      TYPE(CRDS)                      :: GMLoc
!
      BDist   = Zero
      RhoC%D  = Zero
      RhoS%D  = Zero
      DO zq=1,Rho%NExpt
         NQ    = Rho%NQ%I(zq)
         Zeta  = Rho%Expt%D(zq)
         OffQ  = Rho%OffQ%I(zq)
         OffR  = Rho%OffR%I(zq)
!
         LQ    = Rho%Lndx%I(zq)
         LP    = FFELL
!
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
               Dist = SQRT(PQ(1)**2+PQ(2)**2+PQ(3)**2)
               IF(Dist > BDist) THEN
                  BDist      = Dist
                  BOXDist(1) = PQ(1)
                  BOXDist(2) = PQ(2)
                  BOXDist(3) = PQ(3)
               ENDIF
               PFFKetC%D = Zero
               PFFKetS%D = Zero
               CALL HGToSP_Gen(LQ,PiZ,Rho%Co%D(jadd:jadd+LKet-1),PFFKetC%D,PFFKetS%D) 
               IF(NoTranslate(PQ)) THEN
                  DO LM = 0,LSP(LQ)
                     RhoC%D(LM) = RhoC%D(LM)+PFFKetC%D(LM)
                     RhoS%D(LM) = RhoS%D(LM)+PFFKetS%D(LM)
                  ENDDO
               ELSE
                  CALL Regular(LPQ,PQ(1),PQ(2),PQ(3))
                  CALL XLate77(LP,LQ,RhoC%D,RhoS%D,Cpq,Spq,PFFKetC%D,PFFKetS%D)
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      RhoC%D(0) = Zero
      RhoS%D(0) = Zero
!
    END SUBROUTINE RhoToSP
!=======================================================================================
! Calculate the Maxium Ell
!========================================================================================
    SUBROUTINE CalMaxEll(GMLoc)
      TYPE(CRDS)                            :: GMLoc
      INTEGER                               :: L,LP,ML1
      REAL(DOUBLE)                          :: Radius,Dist
      REAL(DOUBLE)                          :: PL,OL,FAC,NL,OLDist
      REAL(DOUBLE)                          :: PFFTau
!
      PFFTau = TauMAC
      ML1=FFELL+1
!
      Radius=SQRT(RDist**2-BDist**2)
!
      IF(BDist > Radius) THEN
         CALL Halt('BDist > Radius increase PFFMaxLay')
      ENDIF
!
      OLDist = (Unsold1(FFEll,RhoC%D,RhoS%D))**(One/DBLE(FFEll))
      Dist   = MAX(BDist,OLDist)
      DO L   = GMLoc%PBC%PFFMaxEll,FFELL
         PL  = Dist**DBLE(L)
         FAC = PL/(Radius**(DBLE(L+1)))
         IF(FAC .LT. PFFTau) THEN
            ML1 = L
            EXIT
         ENDIF
      ENDDO
!
      CALL OpenASCII(OutFile,Out)  
      IF(GM%PBC%PFFOvRide) THEN
         WRITE(Out,990) ML1
!         WRITE(*,990) ML1
      ELSE
         MaxEll = ML1
         IF(MaxEll > FFELL) THEN
            WRITE(Out,991)
!            WRITE(*,991)
            MaxEll = FFELL
         ENDIF
      ENDIF
      CLOSE(Out)
!
990   FORMAT(' OverRide is On: Optimal Ell ==> ',I3)
991   FORMAT(' *** WARNING *** MaxEll > FFEll *** WARNING *** ')
!
    END SUBROUTINE CalMaxEll
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
!!$!========================================================================================
!!$! Calculate the Box Bounds Needed for the Direct Sum
!!$!========================================================================================
!!$  SUBROUTINE BoxBounds(GMLoc)
!!$    INTEGER                          :: I
!!$    INTEGER                          :: IRmin,IRmax
!!$    REAL(DOUBLE)                     :: Radius
!!$    TYPE(CRDS)                       :: GMLoc
!!$!        
!!$    IF(GMLoc%PBC%Dimen==0) THEN
!!$       IRmin = 1
!!$       IRmax = 1
!!$    ELSEIF(GMLoc%PBC%Dimen==1) THEN
!!$       IRmin = 3
!!$       IRmax = MAX(  50,CS_OUT%NCells+1)
!!$    ELSEIF(GMLoc%PBC%Dimen==2) THEN
!!$       IRmin = 9
!!$       IRmax = MAX( 500,CS_OUT%NCells+1)
!!$    ELSEIF(GMLoc%PBC%Dimen==3) THEN
!!$       IRmin = 27
!!$       IRmax = MAX(5000,CS_OUT%NCells+1)
!!$    ENDIF
!!$!
!!$!   Generate the Cells for the Inner Box
!!$!  
!!$    Radius = RDist
!!$    CALL New_CellSet_Sphere(CS_IN,GMLoc%PBC%AutoW,GMLoc%PBC%BoxShape,Radius)
!!$!
!!$!   Rescale the box size  so that the box is within the bounds of IRmin and IRmax
!!$!
!!$    IF(CS_IN%NCells < IRmin) THEN
!!$       DO I=1,1000
!!$          Radius = 1.001D0*Radius
!!$          CALL Delete_CellSet(CS_IN)
!!$          CALL New_CellSet_Sphere(CS_IN,GMLoc%PBC%AutoW,GMLoc%PBC%BoxShape,Radius)
!!$          IF(CS_IN%NCells > IRmin) THEN
!!$             EXIT
!!$          ENDIF
!!$       ENDDO
!!$       RDist = Radius
!!$       RETURN
!!$    ELSEIF(CS_IN%NCells > IRMax) THEN
!!$       DO I=1,1000
!!$          Radius = 0.999D0*Radius
!!$          CALL Delete_CellSet(CS_IN)
!!$          CALL New_CellSet_Sphere(CS_IN,GMLoc%PBC%AutoW,GMLoc%PBC%BoxShape,Radius)
!!$          IF(CS_IN%NCells < IRmax) THEN
!!$             EXIT
!!$          ENDIF
!!$       ENDDO
!!$       RDist = Radius
!!$       RETURN
!!$    ENDIF
!!$    RDist = Radius
!!$!
!!$  END SUBROUTINE BoxBounds
!!$!========================================================================================
!!$! Calculate the Box Bounds Needed for the Direct Sum
!!$!========================================================================================
!!$  SUBROUTINE CalMPB(Q,GMLoc) 
!!$    TYPE(PoleNode)                   :: Q
!!$    INTEGER                          :: I,J,K,LP
!!$    REAL(DOUBLE)                     :: Px,Py,Pz,O_FFEll,NFac,Dist,Mx,My,Mz
!!$    REAL(DOUBLE)                     :: Radius
!!$    TYPE(CRDS)                       :: GMLoc
!!$!
!!$!   PAC Distance (From PoleRoot)
!!$!
!!$    Px=Half*(Q%Box%BndBox(1,2)-Q%Box%BndBox(1,1))
!!$    Py=Half*(Q%Box%BndBox(2,2)-Q%Box%BndBox(2,1))   
!!$    Pz=Half*(Q%Box%BndBox(3,2)-Q%Box%BndBox(3,1))
!!$    PACDist = SQRT(Px*Px+Py*Py+Pz*Pz)
!!$!
!!$!   MAC Distance
!!$!
!!$    LP      = 0
!!$    MACDist = Zero
!!$    DO I=-1,1
!!$       DO J=-1,1
!!$          DO K=-1,1
!!$             Mx = DBLE(I)*BOXDist(1)
!!$             My = DBLE(J)*BOXDist(2)
!!$             Mz = DBLE(K)*BOXDist(3)
!!$             IF(.NOT. GMLoc%PBC%AutoW(1)) Mx = Zero
!!$             IF(.NOT. GMLoc%PBC%AutoW(2)) My = Zero
!!$             IF(.NOT. GMLoc%PBC%AutoW(3)) Mz = Zero
!!$             CALL Regular(MaxELL,Mx,My,Mz)
!!$             CALL XLate77(MaxELL,MaxELL,FarFC,FarFS,Cpq,Spq,RhoC%D,RhoS%D)
!!$!
!!$             O_FFEll  = Unsold2(MaxEll-2,MaxELL,FarFC,FarFS)
!!$             NFac     = FudgeFactorial(LP,MaxELL)
!!$             Dist     = (NFac*O_FFEll/TauMAC)**(One/DBLE(MaxELL+LP+1))
!!$             MACDist  = MAX(MACDist,Dist)
!!$!
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
!!$!
!!$!   BOX Distance
!!$!
!!$    BDist = SQRT(BOXDist(1)**2+BOXDist(2)**2+BOXDist(3)**2)
!!$!
!!$!   Total Distance
!!$!
!!$    RDist = MAX(PACDist,SQRT(MACDIst**2+BDist**2))
!!$!
!!$  END SUBROUTINE CalMPB
