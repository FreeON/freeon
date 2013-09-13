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
!  REAL(DOUBLE)                        :: E_PFF,E_DP,E_QP
  REAL(DOUBLE)                        :: PDist,BDist,RDist
  REAL(DOUBLE),DIMENSION(3)           :: BOXDist,NewPole
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
    SUBROUTINE PBCFarFieldSetUp(GMLoc,RhoLoc,Prog,MaxPFFFEll,ETotal,LatFrc)
      INTEGER                         :: MaxPFFFEll,MaxPFFFLen
      TYPE (HGRho)                    :: RhoLoc
      TYPE(PoleNode)                  :: Q
      TYPE(CRDS)                      :: GMLoc
      INTEGER                         :: I,J,K,M,NC,L,LM,LP,Ell2Use,LL,LDX,ID
      REAL(DOUBLE)                    :: Layers,OL,NL,FAC,Px,Py,Pz
      REAL(DOUBLE),DIMENSION(3)       :: PQ
      TYPE(DBL_VECT)                  :: TensorC,TensorS
      TYPE(DBL_RNK3)                  :: dTensorC,dTensorS
      TYPE(DBL_VECT)                  :: RhoC,RhoS
      REAL(DOUBLE)                    :: E_DP,E_QP,EPFF,DeltaPFF,DeltaFrc,TargetG,ETotal
      REAL(DOUBLE),DIMENSION(3,3)     :: DivCV
      TYPE(DBL_RNK2)                  :: LatFrc_DiP,LatFrc_QuP,LatFrc_Dlt,LatFrc_PFF
      TYPE(DBL_RNK2),OPTIONAL         :: LatFrc
      CHARACTER(LEN=*)                :: Prog
      CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Mssg

      ID(L)=L*(L+1)/2

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
      ! Compute the multipole tensors of the density translated at the cell center
      ! NOTE: the refrence cell center is at the origen, a policy that must be strictly
      ! adhered to.  This policy is administered by PWrap.
      CALL HGToSP_HGRho(MaxPFFFEll,LenPFFFEll,RhoLoc,(/0D0,0D0,0D0/),RhoC%D(0),RhoS%D(0))
      ! Contract expansion of the density with the PFF tensors and put in TenRho
      TenRhoC%D=Zero
      TenRhoS%D=Zero

!!$      DO L=0,8
!!$         DO M=0,L
!!$            LM=LTD(L)+M
!!$            WRITE(*,55)L,M,RhoC%D(LM),TensorC%D(LM)
!!$55          FORMAT(2(I3,","),2(D16.10,", "))
!!$         ENDDO
!!$      ENDDO

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
      ! Calculate the dipole and quadrupole contributions to the total energy.
      ! Note that the quadrupole term is zero for charge neutral systems, but
      ! that non-scf variation of the geometry (ie finite differences) can lead to variation in this term.
      E_DP=0D0
      E_QP=0D0
      DO I=1,3
         IF(GMLoc%PBC%AutoW%I(I)==1)E_DP=E_DP+Half*GMLoc%PBC%DipoleFAC*RhoPoles%DPole%D(I)**2
         IF(GMLoc%PBC%AutoW%I(I)==1)E_QP=E_QP+GMLoc%PBC%QupoleFAC*RhoPoles%MPole*RhoPoles%QPole%D(I)
      ENDDO
      IF(.NOT.PRESENT(LatFrc))THEN
         EPFF=Zero
         Ell2Use=MaxPFFFEll
         DO L=0,MaxPFFFEll
            DeltaPFF=0D0
            DO M=0,L
               LM=LTD(L)+M
               DeltaPFF=DeltaPFF+Two*(RhoC%D(LM)*TenRhoC%D(LM)+RhoS%D(LM)*TenRhoS%D(LM))
            ENDDO
            EPFF=EPFF+DeltaPFF
!!$            WRITE(*,33)L,EPFF+E_DP,EPFF,DELTAPFF
!!$33          FORMAT('Ell = ',I3,', ETotal = ',D12.6,', E_PFF = ',D20.12,', Delta = ',D12.6)

            IF((ABS(DeltaPFF/ETotal)<Thresholds%ETol*1D-2).AND.(L>4))THEN
               Ell2Use=L
               EXIT
            ENDIF

         ENDDO
         IF(Ell2Use==MaxPFFFEll.AND.ETotal.NE.1D2)THEN
            CALL Warn('Maximum angular symmetry in periodic far-field expansion of the energy is too low,'//RTRN//     &
                 'resulting in a relative error '//TRIM(DblToShrtChar(DeltaPFF/ETotal))//" in the total energy, "//RTRN// &
                 'which is larger than the requested accuracy goal of '//TRIM(DblToShrtChar(Thresholds%ETol))//'. '//RTRN// &
                 'Please consider increasing the value of PFFMaxEll from '//TRIM(IntToChar(MaxPFFFEll))//'. '//RTRN// &
                 'Note that increasing PFFMaxEll can actually decrease the cost of QCTC ')
         ENDIF
      ELSE
        ! The crystal field contribution to the lattice forces
         CALL New(LatFrc_PFF,(/3,3/))
         CALL New(LatFrc_Dlt,(/3,3/))
         LatFrc_PFF%D=0D0
         Ell2Use=MaxPFFFEll
         DO L=0,MaxPFFFEll
            LatFrc_Dlt%D=0D0
            DO I=1,3
               DO J=1,3
                  IF(GMLoc%PBC%AutoW%I(I)==1 .AND. GMLoc%PBC%AutoW%I(J)==1) THEN
                     DO M=0,L
                        LM=LTD(L)+M
                        LatFrc_Dlt%D(I,J)=LatFrc_Dlt%D(I,J)-Two*(RhoC%D(LM)*dTenRhoC%D(LM,I,J)+RhoS%D(LM)*dTenRhoS%D(LM,I,J))
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
            LatFrc_PFF%D=LatFrc_PFF%D+LatFrc_Dlt%D
            DeltaFrc=0D0
            DO I=1,3
               DO J=1,3
                  DeltaFrc=MAX(DeltaFrc,ABS(LatFrc_DLT%D(I,J)))
               ENDDO
            ENDDO
            ! Now here we have to be carefull.  Consider that the ABSOLUTE accuracy in the forces
            ! is approx the square root of the RELATIVE accuracy we want in the energy.
!!$            TargetG=SQRT(Thresholds%ETol)*1D-1 ! 1/2 fudge factor to be sure...
!!$            IF(DeltaFrc<TargetG.AND.(L>4))THEN
!!$               Ell2Use=L
!!$               EXIT
!!$            ENDIF

         ENDDO
         !
         EPFF=0D0
         DO L=0,Ell2Use
            DeltaPFF=0D0
            DO M=0,L
               LM=LTD(L)+M
               DeltaPFF=DeltaPFF+Two*(RhoC%D(LM)*TenRhoC%D(LM)+RhoS%D(LM)*TenRhoS%D(LM))
            ENDDO
            EPFF=EPFF+DeltaPFF
         ENDDO
         !
         IF(Ell2Use==MaxPFFFEll.AND.ETotal<1D200)THEN
!!$            CALL Warn(' Maximum angular symmetry in periodic far-field expansion of the gradient is too low,'                   //RTRN//  &
!!$                 '    resulting in an absolute error = '//TRIM(DblToShrtChar(DeltaFrc))//' in the lattice force'                //RTRN//  &
!!$                 '    that is larger than the requested absolute accuracy goal of '//TRIM(DblToShrtChar(TargetG))//'. '         //RTRN//  &
!!$                 '    The cooresponding relative error in the total energy is '//TRIM(DblToShrtChar(DeltaPFF/SQRT(ABS(ETotal))))//RTRN//  &
!!$                 '    Please consider increasing the value of PFFMaxEll from '//TRIM(IntToChar(MaxPFFFEll))//'. '               //RTRN//  &
!!$                 '    Note that increasing PFFMaxEll can actually decrease the cost of QCTC ')
         ENDIF

      ENDIF
      MaxPFFFEll=Ell2Use
      ! For a modern discussion of these issues, see JCP 126 p.124106 (2007)
      ! Here is the INTRINSIC energy due to the crystal Coulomb field
      ! This is the shape dependent Lorentz field + the shape dependent surface term

      CALL MondoLog(DEBUG_MAXIMUM,Prog,'PFFEll = '//TRIM(IntToChar(Ell2Use))// &
                  ', Lorentz Field = <'//TRIM(DblToChar(EPFF+E_DP)) &
                     //'> to within <'//TRIM(DblToShrtChar(ABS(DeltaPFF)))//'>')

!!$      Mssg=ProcessName(Prog,'Ell  = '//TRIM(IntToChar(Ell2Use)))
!!$      Mssg=TRIM(Mssg)//' DIPOLE 1= <'//TRIM(DblToChar(RhoPoles%DPole%D(1)))//'>'
!!$      Mssg=TRIM(Mssg)//' DIPOLE 2= <'//TRIM(DblToChar(RhoPoles%DPole%D(2)))//'>'
!!$      Mssg=TRIM(Mssg)//' DIPOLE 3= <'//TRIM(DblToChar(RhoPoles%DPole%D(3)))//'>'
!!$      Mssg=TRIM(Mssg)//' QUPOLE = <'//TRIM(DblToChar(RhoPoles%QPole%D(2)))//'>'
!!$      WRITE(*,*)TRIM(Mssg)
!!$      Mssg=ProcessName(Prog,'Ell  = '//TRIM(IntToChar(Ell2Use)))
!!$      Mssg=TRIM(Mssg)//' DiP_FoilEnergy = <'//TRIM(DblToChar(E_DP))//'>'
!!$      Mssg=TRIM(Mssg)//' QuP_FoilEnergy = <'//TRIM(DblToChar(E_QP))//'>'
!!$      WRITE(*,*)TRIM(Mssg)
       !
      IF(PRESENT(LatFrc))THEN
         CALL New(LatFrc_DiP,(/3,3/))
         CALL New(LatFrc_QuP,(/3,3/))
         LatFrc_DiP%D=Zero
         LatFrc_QuP%D=Zero
         IF(GMLoc%PBC%Dimen>0) THEN
          DivCV=DivCellVolume(GMLoc%PBC%BoxShape%D,GMLoc%PBC%AutoW%I)
            DO I=1,3
               DO J=1,3
                  IF(GMLoc%PBC%AutoW%I(I)==1.AND.GMLoc%PBC%AutoW%I(J)==1)THEN
                     LatFrc_Dip%D(I,J)=LatFrc_Dip%D(I,J)-(E_DP)*DivCV(I,J)/GMLoc%PBC%CellVolume
                     LatFrc_QuP%D(I,J)=LatFrc_Qup%D(I,J)-(E_QP)*DivCV(I,J)/GMLoc%PBC%CellVolume
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
         !
         ! Initialize the lattice forces with the tin foil and crystal field terms
         ! Note that both of these terms are shape dependent and that only their sum is invarient
         LatFrc%D=LatFrc_PFF%D+LatFrc_DiP%D
         !
!!$         PrintFlags%Key=DEBUG_MAXIMUM
!!$         PrintFlags%MM=DEBUG_FRC
!!$         CALL Print_LatForce(GMLoc,LatFrc_DiP%D,'J  Dipole Volume Lattice Force',Unit_O=6)
!!$!!         CALL Print_LatForce(GMLoc,LatFrc_QuP%D,'J  Qupole Volume Lattice Force',Unit_O=6)
!!$         CALL Print_LatForce(GMLoc,LatFrc_PFF%D,'J  PFF    Lattice Force',Unit_O=6)
!!$         !
!!$         CALL Print_LatForce(GMLoc,LatFrc%D,'J PFF+Dipole+QuadPole  Lattice Force',Unit_O=6)
!!$         !


         CALL Delete(LatFrc_DiP)
         CALL Delete(LatFrc_QuP)
         CALL Delete(LatFrc_PFF)
         CALL Delete(LatFrc_Dlt)
      ENDIF
      !
    END SUBROUTINE PBCFarFieldSetUp
    !====================================================================================
    !   Calculate the Extrinisc (unique) Far-Field component of the J matrix,
    !   including the non-unique surface (tin foil) and Lorentz field (Nijboer-Dewette) terms
    !====================================================================================
    FUNCTION CTraxFF(P,G,E,PiZ) RESULT(CTFF)
      TYPE(PrimPair)              :: P
      TYPE(CRDS)                  :: G
      REAL(DOUBLE), DIMENSION(:)  :: E
      REAL(DOUBLE)                :: PiZ,CTFF
      IF(G%PBC%Dimen==0)THEN
         CTFF=Zero
         RETURN
      ELSE
         CTFF=TinFoil(P,G,E,PiZ)
         CTFF=CTFF+FarField(P,E,G)
      ENDIF
    END FUNCTION CTraxFF
    !====================================================================================
    !   Calculate the FarField (Lorentz field) Component of the J matrix
    !====================================================================================
    FUNCTION FarField(Prim,HGBra,G)     RESULT(FF)
      TYPE(PrimPair)                   :: Prim
      INTEGER                          :: LM
      REAL(DOUBLE)                     :: PiZ,FF
      REAL(DOUBLE), DIMENSION(3)       :: PQ,HGDipole
      REAL(DOUBLE), DIMENSION(1:)      :: HGBra
      TYPE(CRDS)                       :: G
      REAL(DOUBLE),DIMENSION(0:FFLen)  :: KetC,KetS
      REAL(DOUBLE),DIMENSION(0:FFLen)  :: BraC,BraS
      !
      FF=Zero
      IF(G%PBC%Dimen==0) RETURN
      ! Transform <Bra| coefficients from HG to SP
      PiZ=(Pi/Prim%Zeta)**(ThreeHalves)
      CALL HGToSP_Gen(Prim%Ell,PiZ,HGBra,KetC,KetS)
      PQ=Prim%Pw
      ! Contract
      IF(NoTranslate(PQ)) THEN
         DO LM = 0,LSP(Prim%Ell)
            FF=FF+(KetC(LM)*TenRhoC%D(LM)+KetS(LM)*TenRhoS%D(LM))
         ENDDO
      ELSE
         BraC(0:LenPFFFEll)=Zero
         BraS(0:LenPFFFEll)=Zero
         CALL Regular(MaxPFFFELL,PQ(1),PQ(2),PQ(3))
         CALL XLate77(MaxPFFFEll,Prim%Ell,BraC(0),BraS(0),Cpq(0),Spq(0),KetC(0),KetS(0))
         DO LM=0,LenPFFFEll
            FF=FF+(BraC(LM)*TenRhoC%D(LM)+BraS(LM)*TenRhoS%D(LM))
         ENDDO
      ENDIF
    END FUNCTION FarField
    !----------------------------------------------------------------------------------------------
    ! This function adds in the Tin Foil (surface) term to the Coulomb matrix
    !----------------------------------------------------------------------------------------------
    FUNCTION TinFoil(P,G,E,PiZ) RESULT(TF)
      TYPE(PrimPair)                            :: P
      TYPE(CRDS)                                :: G
       REAL(DOUBLE), DIMENSION(1:)               :: E
      REAL(DOUBLE)                              :: PiZ,TF,MP
      REAL(DOUBLE), DIMENSION(3)                :: PQ,DP
      REAL(DOUBLE), DIMENSION(6)                :: QP
      INTEGER                                   :: I
      !
      TF=Zero
      PQ=P%Pw
      MP=PiZ*E(1)
      DP=DPole(P%Ell,PiZ,PQ(1),PQ(2),PQ(3),E)
      QP=QPole(P%Ell,P%Zeta,PiZ,PQ(1),PQ(2),PQ(3),E)
      DO I=1,3
         IF(G%PBC%AutoW%I(I)==1)THEN
            TF=TF+G%PBC%DipoleFAC*DP(I)*RhoPoles%DPole%D(I)*Half
            TF=TF+MP*G%PBC%QupoleFAC*RhoPoles%QPole%D(I)*Half
         ENDIF
      ENDDO
    END FUNCTION TinFoil
    !----------------------------------------------------------------------------------------------
    ! This function returns the derivative of the Tin Foil term with respect to
    ! change in McMurchie-Davidson basis function elements E, as outlined by Helgaker and Taylor
    ! in 1992 TCA publication.
    !----------------------------------------------------------------------------------------------
    FUNCTION dTinFoil(P,G,PdE,PiZ) RESULT(TF)
      TYPE(PrimPair)              :: P
      TYPE(CRDS)                  :: G
      REAL(DOUBLE)                :: PiZ,TF,dMP
      REAL(DOUBLE), DIMENSION(3)  :: PQ,dDP
      REAL(DOUBLE), DIMENSION(6)  :: dQP
      REAL(DOUBLE), DIMENSION(1:) :: PdE
      INTEGER                     :: I
      TF=Zero
      PQ=P%Pw
      ! Monopole derivative
      dMP=PiZ*PdE(1)
      ! Dipole derivative
      dDP=DPole(P%Ell,PiZ,PQ(1),PQ(2),PQ(3),PdE)
      !
      DO I=1,3
         IF(G%PBC%AutoW%I(I)==1)THEN
            TF=TF+G%PBC%DipoleFAC*dDP(I)*RhoPoles%DPole%D(I)*Half
            TF=TF+dMP*G%PBC%QupoleFAC*RhoPoles%QPole%D(I)*Half
         ENDIF
      ENDDO
    END FUNCTION DTinFoil
    !----------------------------------------------------------------------------------------------
    ! This function returns the derivative of the Tin Foil term with respect to
    ! change in the position-center variable (P-C), which is just the wrapped and
    ! re-centered position, Pw.
    !----------------------------------------------------------------------------------------------
    FUNCTION DerivTF(Prim,E,GMLoc) RESULT(TF)
      TYPE(PrimPair)                   :: Prim
      TYPE(CRDS)                       :: GMLoc
      REAL(DOUBLE)                     :: dRi2
      REAL(DOUBLE), DIMENSION(:)       :: E
      REAL(DOUBLE), DIMENSION(3)       :: PQ,TF,HGDipole
      INTEGER                          :: I
      TF=Zero
      ! Note sign flip because of inverted nuclear/electron sign (electrons are positive)
      DO I=1,3
         IF(GMLoc%PBC%AutoW%I(I)==1)THEN
            ! Dipole term
            TF(I)=TF(I)-Half*GMLoc%PBC%DipoleFAC*E(1)*RhoPoles%DPole%D(I)
         ENDIF
      ENDDO
    END FUNCTION DerivTF
!----------------------------------------------------------------------------------------------
!   Print out Periodic Info
!----------------------------------------------------------------------------------------------
!!$    SUBROUTINE Print_Periodic(GMLoc,Prog,Unit_O)
!!$      INTEGER,OPTIONAL                :: Unit_O
!!$      INTEGER                         :: Unit
!!$      TYPE(CRDS)                      :: GMLoc
!!$      CHARACTER(LEN=*)                :: Prog
!!$      CHARACTER(LEN=DEFAULT_CHR_LEN)  :: Mssg
!!$      REAL(DOUBLE)                    :: Layers
!!$!
!!$#ifdef PARALLEL
!!$  IF(MyID==ROOT) THEN
!!$#endif
!!$      Layers = SQRT(CS_IN%CellCarts%D(1,1)**2+CS_IN%CellCarts%D(2,1)**2+CS_IN%CellCarts%D(3,1)**2)/MaxBoxDim(GMLoc)
!!$      IF(GMLoc%PBC%Dimen > 0) THEN
!!$         Unit=OpenPU(Unit_O=Unit_O)
!!$         IF(PrintFlags%Key > DEBUG_MINIMUM) THEN
!!$            Mssg=ProcessName(Prog,TRIM(IntToChar(GMLoc%PBC%Dimen))//'-D periodics')   &
!!$                 //'MxL = '//TRIM(IntToChar(MaxPFFFEll))                                  &
!!$                 //', PAC Cells = '//TRIM(IntToChar(CS_OUT%NCells))                   &
!!$                 //', MAC Cells = '//TRIM(IntToChar(CS_IN%NCells))
!!$            WRITE(Unit,*)TRIM(Mssg)
!!$            Mssg=ProcessName(Prog,'FF Energy')//'<FF> = '//TRIM(DblToChar(E_PFF))
!!$            WRITE(Unit,*)TRIM(Mssg)
!!$         ENDIF
!!$         CALL ClosePU(Unit)
!!$      ENDIF
!!$#ifdef PARALLEL
!!$   ENDIF
!!$#endif
!!$!
!!$100   FORMAT('========================================Periodic Information======================================')
!!$101   FORMAT(' MaxPFFFEll         = ',I3,'     Dimension = ',I2)
!!$102   FORMAT(' Dipole Moment      = (',F12.6,',',F12.6,',',F12.6,')')
!!$103   FORMAT(' PAC Distance       = ',F6.2,' BOX Distance = ',F6.2,' Cell Distance = ',F6.2)
!!$108   FORMAT(' No. of Layers      = ',F6.2)
!!$104   FORMAT(' Outer No. of Cells = ',I4)
!!$110   FORMAT(' Inner No. of Cells = ',I4)
!!$105   FORMAT(' Correction to the Energy:')
!!$106   FORMAT('   PFF = ',E14.6,'  Dipole = ',E14.6)
!!$107   FORMAT('=========================================END======================================================')
!!$!
!!$    END SUBROUTINE Print_Periodic
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
