MODULE PBCFarField
  USE Derivedtypes
  USE GlobalScalars   
  USE GlobalObjects
  USE ProcessControl
  USE Indexing
  USE InOut
  USE Macros
  USE Thresholding
  USE BoundingBox
  USE MondoPoles
  USE PoleTree
  USE Globals
  IMPLICIT NONE
!  
  REAL(DOUBLE)                        :: E_PFF,E_DP
  REAL(DOUBLE)                        :: PDist,BDist,RDist
  REAL(DOUBLE),DIMENSION(3)           :: BOXDist
!
  ! Variables that contain the Periodic Far Field potential
  TYPE(DBL_VECT)                      :: TenRhoC,TenRhoS
  TYPE(DBL_RNK3)                      :: dTenRhoC,dTenRhoS
  ! Some other shit
  TYPE(DBL_VECT)                      :: PFFBraC,PFFBraS
  TYPE(DBL_VECT)                      :: PFFKetC,PFFKetS


  CONTAINS
!====================================================================================
!   Setup the PBCFarField Matrix. 
!====================================================================================
    SUBROUTINE PBCFarFuckingFieldSetUp(GMLoc,RhoLoc,MaxPFFFEll,EPFF,LatFrc)
      INTEGER                         :: MaxPFFFEll,MaxPFFFLen
      TYPE (HGRho)                    :: RhoLoc
      TYPE(PoleNode)                  :: Q
      TYPE(CRDS)                      :: GMLoc
      INTEGER                         :: I,J,K,NC,L,LM,LP
      REAL(DOUBLE)                    :: Layers,OL,NL,FAC,Px,Py,Pz
      REAL(DOUBLE),DIMENSION(3)       :: PQ
      TYPE(DBL_VECT)                  :: TensorC,TensorS
      TYPE(DBL_RNK3)                  :: dTensorC,dTensorS
      TYPE(DBL_VECT)                  :: RhoC,RhoS
      REAL(DOUBLE)                    :: EPFF
      REAL(DOUBLE),DIMENSION(3,3)     :: DivCV
      TYPE(DBL_RNK2)                  :: LatFrc_Dip,LatFrc_PFF
      TYPE(DBL_RNK2),OPTIONAL         :: LatFrc
      !
      IF(GMLoc%PBC%Dimen==0) RETURN
      !
      IF(MaxPFFFEll>FFEll) &
           STOP ' MaxPFFFEll to large in PBCFarFuckingFieldSetUp '
      LenPFFFEll=LSP(MaxPFFFEll)
      !Allocate and fill global PBC Far Field tensors
      CALL New(TensorC,LSP(2*MaxPFFFEll),0)
      CALL New(TensorS,LSP(2*MaxPFFFEll),0) 
      CALL New(dTensorC,(/LSP(2*MaxPFFFEll),3,3/),(/0,1,1/))
      CALL New(dTensorS,(/LSP(2*MaxPFFFEll),3,3/),(/0,1,1/))
      ! Get the script M tensors and their derivatives from HDF 
      CALL Get(TensorC,'PFFTensorC')
      CALL Get(TensorS,'PFFTensorS')
      CALL Get(dTensorC,'dPFFTensorC')
      CALL Get(dTensorS,'dPFFTensorS')
      ! 
      CALL New(TenRhoC,LenPFFFEll,0)
      CALL New(TenRhoS,LenPFFFEll,0)
      CALL New(dTenRhoC,(/LenPFFFEll,3,3/),(/0,1,1/))
      CALL New(dTenRhoS,(/LenPFFFEll,3,3/),(/0,1,1/))
      ! Allocate local tensors that hold the density translated to the cell center
      CALL New(RhoC,LenPFFFEll,0)  
      CALL New(RhoS,LenPFFFEll,0) 
      ! Compute the multipole tensors of the density translated to the cell center
      CALL HGToSP_HGRho(MaxPFFFEll,LenPFFFEll,RhoLoc,GMLoc%PBC%CellCenter%D,RhoC%D(0),RhoS%D(0))
      ! Contract expansion of the density with the PFF tensors and put in TenRho 
      TenRhoC%D=Zero
      TenRhoS%D=Zero
      CALL CTraX77(MaxPFFFEll,MaxPFFFEll,TenRhoC%D,TenRhoS%D,TensorC%D,TensorS%D,RhoC%D,RhoS%D)  
      ! Contract expansion of the density with the derivative PFF tensors and put in dTenRho 
      dTenRhoC%D=Zero
      dTenRhoS%D=Zero
      DO I=1,3
         DO J=1,3
            IF(GMLoc%PBC%AutoW%I(I)==1 .AND. GMLoc%PBC%AutoW%I(J)==1) THEN
               CALL CTraX77(MaxPFFFEll,MaxPFFFEll,dTenRhoC%D(:,I,J),dTenRhoS%D(:,I,J), &
                    dTensorC%D(:,I,J),dTensorS%D(:,I,J),RhoC%D,RhoS%D)  
            ENDIF
         ENDDO
      ENDDO
      ! Here is where we reset MaxPFFFEll based on cell-cell MAC
      ! but punt for now...
      !      MaxPFFFEll=CalMaxPFFFEll(GMLoc)
      !
      !
      !
      !  Calculate energy of the crystal field
      EPFF=Zero
      DO LM=0,LenPFFFEll
         EPFF=EPFF+Two*RhoC%D(LM)*TenRhoC%D(LM)+Two*RhoS%D(LM)*TenRhoS%D(LM)
      ENDDO
      WRITE(*,*)' E_PFF = ',E_PFF
      IF(PRESENT(LatFrc))THEN
         ! The tin foil (dipole) contribution to the energy and lattice forces
         IF(GMLoc%PBC%Dimen == 1) THEN
            E_DP = 0.0D0
         ELSEIF(GMLoc%PBC%Dimen == 2) THEN
            IF(GMLoc%PBC%AutoW%I(1)==0 ) E_DP = Two*GMLoc%PBC%DipoleFAC*RhoPoles%DPole%D(1)**2
            IF(GMLoc%PBC%AutoW%I(2)==0 ) E_DP = Two*GMLoc%PBC%DipoleFAC*RhoPoles%DPole%D(2)**2         
            IF(GMLoc%PBC%AutoW%I(3)==0 ) E_DP = Two*GMLoc%PBC%DipoleFAC*RhoPoles%DPole%D(3)**2
         ELSEIF(GMLoc%PBC%Dimen == 3) THEN 
            E_DP = Two*GMLoc%PBC%DipoleFAC*(RhoPoles%DPole%D(1)**2+RhoPoles%DPole%D(2)**2+RhoPoles%DPole%D(3)**2)
         ELSE
            CALL Halt('PBCFarField: Unknown dimension <'//IntToChar(GMLoc%PBC%Dimen)//'>')
         ENDIF
         CALL New(LatFrc_Dip,(/3,3/))
         LatFrc_Dip%D = Zero
         IF(GMLoc%PBC%Dimen > 0) THEN
            DivCV=DivCellVolume(GMLoc%PBC%BoxShape%D,GMLoc%PBC%AutoW%I)
            DO I=1,3
               DO J=1,3
                  IF(GMLoc%PBC%AutoW%I(I)==1.AND.GMLoc%PBC%AutoW%I(J)==1)THEN
                     LatFrc_Dip%D(I,J)=LatFrc_Dip%D(I,J)-E_DP*DivCV(I,J)/GMLoc%PBC%CellVolume
                  ENDIF
               ENDDO
            ENDDO
         ENDIF

         ! The crystal field contribution to the lattice forces
         CALL New(LatFrc_PFF,(/3,3/))
         LatFrc_PFF%D = Zero
         IF(GMLoc%PBC%Dimen > 0) THEN
            DO I=1,3
               DO J=1,3
                  IF(GMLoc%PBC%AutoW%I(I)==1 .AND. GMLoc%PBC%AutoW%I(J)==1) THEN
                     DO K=0,LenPFFFEll
                        LatFrc_PFF%D(I,J)=LatFrc_PFF%D(I,J)-Two*(RhoC%D(K)*dTenRhoC%D(K,I,J)+RhoS%D(K)*dTenRhoS%D(K,I,J))
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
         PrintFlags%Key=DEBUG_MAXIMUM	
         PrintFlags%MM=DEBUG_FRC
         CALL Print_LatForce(GMLoc,LatFrc_PFF%D,'J  PFF   Lattice Force')
         CALL Print_LatForce(GMLoc,LatFrc_PFF%D,'J  PFF   Lattice Force',Unit_O=6)

         ! Initialize the lattice forces with the tin foil and crystal field terms
         LatFrc%D=LatFrc_Dip%D+LatFrc_PFF%D
         CALL Delete(LatFrc_Dip)
         CALL Delete(LatFrc_PFF)
      ENDIF

      !
    END SUBROUTINE PBCFarFuckingFieldSetUp
!====================================================================================
!   Calculate the FarField Component of the J matrix
!====================================================================================
    FUNCTION CTraxFF(Prim,HGBra,GMLoc) RESULT(CTFF)
      TYPE(PrimPair)                   :: Prim
      INTEGER                          :: LM
      REAL(DOUBLE)                     :: PiZ,CTFF
      REAL(DOUBLE), DIMENSION(3)       :: PQ,HGDipole 
      REAL(DOUBLE), DIMENSION(1:)      :: HGBra
      TYPE(CRDS)                       :: GMLoc
      REAL(DOUBLE),DIMENSION(0:FFLen)  :: KetC,KetS
      REAL(DOUBLE),DIMENSION(0:FFLen)  :: BraC,BraS

      CTFF=Zero
      IF(GMLoc%PBC%Dimen==0) RETURN
!     Transform <Bra| coefficients from HG to SP
      PiZ=(Pi/Prim%Zeta)**(ThreeHalves)
      CALL HGToSP_Gen(Prim%Ell,PiZ,HGBra,KetC,KetS)   
      PQ=Prim%P-GMLoc%PBC%CellCenter%D   
!     Contract
      IF(NoTranslate(PQ)) THEN
         DO LM = 0,LSP(Prim%Ell)
            CTFF=CTFF+KetC(LM)*TenRhoC%D(LM)+KetS(LM)*TenRhoS%D(LM)
         ENDDO
      ELSE
         BraC(0:LenPFFFEll)=Zero 
         BraS(0:LenPFFFEll)=Zero 
         CALL Regular(MaxPFFFELL,PQ(1),PQ(2),PQ(3))
         CALL XLate77(MaxPFFFEll,Prim%Ell,BraC(0),BraS(0),Cpq(0),Spq(0),KetC(0),KetS(0))
         DO LM=0,LenPFFFEll
            CTFF=CTFF+BraC(LM)*TenRhoC%D(LM)+BraS(LM)*TenRhoS%D(LM)
         ENDDO
      ENDIF
!     Include the Dipole correction to FarFC and FarFS
      PQ=Prim%P-GMLoc%PBC%CellCenter%D 
      IF(GMLoc%PBC%Dimen==2) THEN
         HGDipole = CalculateDiPole(Prim%Ell,Prim%Zeta,PQ(1),PQ(2),PQ(3),HGBra(1:))
         IF(GMLoc%PBC%AutoW%I(1)==0) CTFF  = CTFF + GMLoc%PBC%DipoleFAC*HGDipole(1)*RhoPoles%DPole%D(1)         
         IF(GMLoc%PBC%AutoW%I(2)==0) CTFF  = CTFF + GMLoc%PBC%DipoleFAC*HGDipole(2)*RhoPoles%DPole%D(2)
         IF(GMLoc%PBC%AutoW%I(3)==0) CTFF  = CTFF + GMLoc%PBC%DipoleFAC*HGDipole(3)*RhoPoles%DPole%D(3)
      ELSEIF(GMLoc%PBC%Dimen==3) THEN
         HGDipole = CalculateDiPole(Prim%Ell,Prim%Zeta,PQ(1),PQ(2),PQ(3),HGBra(1:))
         CTFF  = CTFF + GMLoc%PBC%DipoleFAC*(HGDipole(1)*RhoPoles%DPole%D(1)         &
                                           + HGDipole(2)*RhoPoles%DPole%D(2)         &
                                           + HGDipole(3)*RhoPoles%DPole%D(3) )
      ENDIF
    END FUNCTION CTraxFF
!---------------------------------------------------------------------------------------------- 
!   Print out Periodic Info
!----------------------------------------------------------------------------------------------  
    SUBROUTINE Print_Periodic(GMLoc,Prog,Unit_O)
      INTEGER,OPTIONAL                :: Unit_O
      INTEGER                         :: Unit
      TYPE(CRDS)                      :: GMLoc
      CHARACTER(LEN=*)                :: Prog
      CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Mssg
      REAL(DOUBLE)                    :: Layers
!
#ifdef PARALLEL
  IF(MyID==ROOT) THEN
#endif
      Layers = SQRT(CS_IN%CellCarts%D(1,1)**2+CS_IN%CellCarts%D(2,1)**2+CS_IN%CellCarts%D(3,1)**2)/MaxBoxDim(GMLoc)
      IF(GMLoc%PBC%Dimen > 0) THEN
         Unit=OpenPU(Unit_O=Unit_O)
         IF(PrintFlags%Key > DEBUG_MINIMUM) THEN
            Mssg=ProcessName(Prog,TRIM(IntToChar(GMLoc%PBC%Dimen))//'-D periodics')   &
                 //'MxL = '//TRIM(IntToChar(MaxPFFFEll))                                  &
                 //', PAC Cells = '//TRIM(IntToChar(CS_OUT%NCells))                   &
                 //', MAC Cells = '//TRIM(IntToChar(CS_IN%NCells))                    
            WRITE(Unit,*)TRIM(Mssg)
            Mssg=ProcessName(Prog,'FF Energy')//'<FF> = '//TRIM(DblToChar(E_PFF))
            WRITE(Unit,*)TRIM(Mssg)
         ENDIF
         CALL ClosePU(Unit)
      ENDIF
#ifdef PARALLEL
   ENDIF
#endif
!
100   FORMAT('========================================Periodic Information======================================')
101   FORMAT(' MaxPFFFEll         = ',I3,'     Dimension = ',I2)
102   FORMAT(' Dipole Moment      = (',F12.6,',',F12.6,',',F12.6,')')
103   FORMAT(' PAC Distance       = ',F6.2,' BOX Distance = ',F6.2,' Cell Distance = ',F6.2)
108   FORMAT(' No. of Layers      = ',F6.2)
104   FORMAT(' Outer No. of Cells = ',I4)
110   FORMAT(' Inner No. of Cells = ',I4)
105   FORMAT(' Correction to the Energy:')
106   FORMAT('   PFF = ',E14.6,'  Dipole = ',E14.6)
107   FORMAT('=========================================END======================================================')
!
    END SUBROUTINE Print_Periodic
!========================================================================================
! If QP is < TOL, do not translate
!========================================================================================
    FUNCTION NoTranslate(X) 
      REAL(DOUBLE),DIMENSION(3) :: X
      REAL(DOUBLE),PARAMETER    :: TOL=1.0D-12
      LOGICAL                   :: NoTranslate
      NoTranslate = (ABS(X(1)).LT.TOL) .AND. (ABS(X(2)).LT.TOL) .AND. (ABS(X(3)).LT.TOL)
    END FUNCTION NoTranslate

!
END MODULE PBCFarField

