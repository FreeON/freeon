!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
!  LAY FUNCTIONALS DEFINED IN MMA/FUNCTIONALS ONTO A GRID
!  Author: Matt Challacombe
!  For reference #s see: http://www.dl.ac.uk/DFTlib/
!-----------------------------------------------------------
MODULE Functionals
  USE DerivedTypes
  USE ProcessControl
  USE Parse
  !INTERFACE ExcOnTheGrid
  !   MODULE PROCEDURE ExcOnTheGrid_ClSh
  !END INTERFACE
  !  Parsing keys for <Options.Model=>
  CHARACTER(LEN=*),  PARAMETER :: MODEL_OPTION='ModelChem'
  !  Coulomb only
  CHARACTER(LEN=*),  PARAMETER :: MODEL_Hartree       ='Hartree'
  !  Exchange only
  CHARACTER(LEN=*),  PARAMETER :: MODEL_ExactX       ='HF'
  CHARACTER(LEN=*),  PARAMETER :: MODEL_SD           ='SlaterDirac'
  CHARACTER(LEN=*),  PARAMETER :: MODEL_XA           ='XAlpha'
  CHARACTER(LEN=*),  PARAMETER :: MODEL_B88x         ='B88x'
  CHARACTER(LEN=*),  PARAMETER :: MODEL_PBEx         ='PBEx'
  CHARACTER(LEN=*),  PARAMETER :: MODEL_PW91x        ='PW91x'
  !  Pure exchange-correlation functionals
  CHARACTER(LEN=*),  PARAMETER :: MODEL_VWN3         ='VWN3xc'
  CHARACTER(LEN=*),  PARAMETER :: MODEL_VWN5         ='VWN5xc'
  CHARACTER(LEN=*),  PARAMETER :: MODEL_PW91PW91     ='PW91xc'
  CHARACTER(LEN=*),  PARAMETER :: MODEL_PW91LYP      ='PW91LYP'
  CHARACTER(LEN=*),  PARAMETER :: MODEL_BLYP         ='BLYPxc'
  CHARACTER(LEN=*),  PARAMETER :: MODEL_BPW91        ='BPW91xc'
  CHARACTER(LEN=*),  PARAMETER :: MODEL_PBEPBE       ='PBExc'
  CHARACTER(LEN=*),  PARAMETER :: MODEL_HCTH93       ='HCTH93xc'
  CHARACTER(LEN=*),  PARAMETER :: MODEL_HCTH120      ='HCTH120xc'
  CHARACTER(LEN=*),  PARAMETER :: MODEL_HCTH147      ='HCTH147xc'
  CHARACTER(LEN=*),  PARAMETER :: MODEL_HCTH407      ='HCTH407xc'
  CHARACTER(LEN=*),  PARAMETER :: MODEL_XLYP         ='XLYPxc'
  !  Hybrid exchange-correlation functionals
  CHARACTER(LEN=*),  PARAMETER :: MODEL_B3LYP_VWN3   ='B3LYP'
  CHARACTER(LEN=*),  PARAMETER :: MODEL_B3LYP_VWN5   ='B3LYP/VWN5'
  CHARACTER(LEN=*),  PARAMETER :: MODEL_PBE0         ='PBE0'
  CHARACTER(LEN=*),  PARAMETER :: MODEL_X3LYP        ='X3LYP'
  !-----------------------------------------------------------------------------------------------------
  !  Numerical keys for model chemistries
  INTEGER, PARAMETER :: NO_EXCHANGE      =      -1 ! Use no exchange (or correlation); this is bare Coulomb
  INTEGER, PARAMETER :: EXACT_EXCHANGE   =03045805 ! Use exact Hartree-Fock exchange
  INTEGER, PARAMETER :: HAS_DFT          =10000000 ! Any key above this has DFT
  INTEGER, PARAMETER :: HYBRID_B3LYP_VWN3=10268305 ! Use Beckes B3LYP model with VWN3 LSD correlation
  INTEGER, PARAMETER :: HYBRID_B3LYP_VWN5=10429495 ! Use Beckes B3LYP model with VWN5 LSD correlation
  INTEGER, PARAMETER :: HYBRID_X3LYP     =10685608 ! Xtended LYP hybrid
  INTEGER, PARAMETER :: HYBRID_PBE0      =14506981 ! Use Adamo and Barones PBE0 model
  INTEGER, PARAMETER :: HAS_HF           =20000000 ! Any key below this has exact exchange
  INTEGER, PARAMETER :: SD_EXCHANGE      =32098243 ! Slater-Dirac Uniform Electron Gas exchange
  INTEGER, PARAMETER :: XA_EXCHANGE      =33890824 ! X-Alpha LDA exchange
  INTEGER, PARAMETER :: B88_EXCHANGE     =34065852 ! Becke 88 exchange
  INTEGER, PARAMETER :: PW91_EXCHANGE    =38458583 ! Perdew Wang 91 exchange
  INTEGER, PARAMETER :: PBE_EXCHANGE     =36792335 ! Perdew Burke Enzerhoff exchange
  INTEGER, PARAMETER :: PURE_VWN3_LSD    =40022304 ! Slater exchange and VWN3 LSDA correlation
  INTEGER, PARAMETER :: PURE_VWN5_LSD    =40038402 ! Slater exchange and VWN5 LSDA correlation
  INTEGER, PARAMETER :: PURE_PBE_PBE     =42509434 ! Perdew Burke Enzerhoff GGA exchange-correlation
  INTEGER, PARAMETER :: PURE_B88_LYP     =43083038 ! B88 exchange with LYP correlation
  INTEGER, PARAMETER :: PURE_B88_PW91    =43083039 ! B88 exchange with PW91 correlation
  INTEGER, PARAMETER :: PURE_PW91_LYP    =42034802 ! PW91 exchange with LYP correlation
  INTEGER, PARAMETER :: PURE_PW91_PW91   =41243253 ! PW91 exchange with PW91 correlation
  INTEGER, PARAMETER :: PURE_HCTH93      =41243260 ! HCTH Handy's family functional
  INTEGER, PARAMETER :: PURE_HCTH120     =41243261 !
  INTEGER, PARAMETER :: PURE_HCTH147     =41243262 !
  INTEGER, PARAMETER :: PURE_HCTH407     =41243263 !
  INTEGER, PARAMETER :: PURE_XLYP        =10685609 ! Xtended LYP hybrid

  ! Spin model keys.
  CHARACTER(LEN=*), PARAMETER :: SPIN_MODEL_OPTION       = "SpinModel"
  CHARACTER(LEN=*), PARAMETER :: SPIN_MODEL_RESTRICTED   = "R"
  CHARACTER(LEN=*), PARAMETER :: SPIN_MODEL_UNRESTRICTED = "U"
  CHARACTER(LEN=*), PARAMETER :: SPIN_MODEL_GENERALIZED  = "G"

  INTEGER, PARAMETER          :: SPIN_MODEL_RESTRICTED_VALUE   = 1
  INTEGER, PARAMETER          :: SPIN_MODEL_UNRESTRICTED_VALUE = 2
  INTEGER, PARAMETER          :: SPIN_MODEL_GENERALIZED_VALUE  = 4

  !  Avoid under and over flows
  REAL(DOUBLE), PARAMETER :: NoNAN=1.D-30
  !  Global variable set to the current model chemistry
  INTEGER :: ModelChem
  INTEGER :: NSMat
  !  Global intermediates for optimized functional forms
  REAL(DOUBLE),DIMENSION(500) :: UTmp
  REAL(DOUBLE),DIMENSION(500) :: VTmp
  !------------
CONTAINS !
  !====================================================================================
  !     RETURN A STRING WITH THE FUNCTIONALS NAME/DISCRIPTION
  !====================================================================================
  FUNCTION FunctionalName(Key) RESULT(Name)
    INTEGER            :: Key
    CHARACTER(LEN=30)  :: Name
    SELECT CASE(Key)
    CASE(EXACT_EXCHANGE)
       Name='HFx'
    CASE(SD_EXCHANGE)
       Name='Slater-Dirac exchange'
    CASE(XA_EXCHANGE)
       Name='X-Alpha exchange'
    CASE(B88_EXCHANGE)
       Name='B88x'
    CASE(PW91_EXCHANGE)
       Name='PW91x'
    CASE(PBE_EXCHANGE)
       Name='PBEx'
    CASE(PURE_VWN3_LSD)
       Name='X-Alpha/VWN3c'
    CASE(PURE_VWN5_LSD)
       Name='X-Alpha/VWN5c'
    CASE(PURE_PW91_PW91)
       Name='PW91x/PW91c'
    CASE(PURE_B88_LYP)
       Name='B88x/LYPc'
    CASE(PURE_B88_PW91)
       Name='B88x/PW91c'
    CASE(PURE_PBE_PBE)
       Name='PBEx/PBEc'
    CASE(PURE_XLYP)
       Name='XLYP'
    CASE(HYBRID_PBE0)
       Name='PBE0'
    CASE(HYBRID_B3LYP_VWN3)
       Name='B3LYP(VWN3)'
    CASE(HYBRID_B3LYP_VWN5)
       Name='B3LYP(VWN5)'
    CASE(HYBRID_X3LYP)
       Name='X3LYP'
    CASE DEFAULT
       CALL Halt('Unknown functional key = '//TRIM(IntToChar(Key)))
    END SELECT
  END FUNCTION FunctionalName

  FUNCTION ExactXScale(Key) RESULT(Scalar)
    INTEGER      :: Key
    REAL(DOUBLE) :: Scalar
    IF(Key==EXACT_EXCHANGE)THEN
       Scalar=1.0D0
    ELSEIF(Key==HYBRID_PBE0)THEN
       Scalar=0.25D0
    ELSEIF(Key==HYBRID_B3LYP_VWN3.OR. &
         Key==HYBRID_B3LYP_VWN5)THEN
       Scalar=0.20D0
    ELSEIF(Key==HYBRID_X3LYP)THEN
       Scalar=0.218D0
    ELSE
       Scalar=0.0D0
    ENDIF
  END FUNCTION ExactXScale
  !
  FUNCTION HasDFT(Key)
    INTEGER  :: Key
    LOGICAL  :: HasDFT
    IF(Key>Has_DFT)THEN
       HasDFT=.TRUE.
    ELSE
       HasDFT=.FALSE.
    ENDIF
  END FUNCTION HasDFT
  !
  FUNCTION HasHF(Key)
    INTEGER  :: Key
    LOGICAL  :: HasHF
    IF(Key>0.AND.Key<Has_HF)THEN
       HasHF=.TRUE.
    ELSE
       HasHF=.FALSE.
    ENDIF
  END FUNCTION HasHF
  !====================================================================================
  !
  !====================================================================================
  SUBROUTINE ExcOnTheGrid_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam)
    INTEGER                    :: NGrid,NSDen
    REAL(DOUBLE)               :: a0,ax,ac,Ax1,Ax2
    REAL(DOUBLE), DIMENSION(1:NGrid)   :: Rho,AbsGradRho2
    REAL(DOUBLE), DIMENSION(1:NGrid)   :: E,dEdRho,dEdGam
    REAL(DOUBLE), DIMENSION(1:3*NGrid) :: Buf
    !---------------------------------------------------------------------------------------
    N2=2*NGrid
    N3=3*NGrid
    E=Zero
    dEdRho=Zero
    dEdGam=Zero
    !
    SELECT CASE(ModelChem)
    CASE(SD_EXCHANGE)
       CALL SDx_ClSh  (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One)
    CASE(XA_EXCHANGE)
       CALL XAx_ClSh  (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One)
    CASE(B88_EXCHANGE)
       CALL B88x_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One)
    CASE(PW91_EXCHANGE)
       CALL PW91x_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One)
    CASE(PBE_EXCHANGE)
       CALL PBEx_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One)
    CASE(PURE_VWN3_LSD)
       CALL XAx_ClSh  (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One)
       CALL VWN3c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One)
    CASE(PURE_VWN5_LSD)
       CALL XAx_ClSh  (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One)
       CALL VWN5c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One)
    CASE(PURE_PW91_PW91)
       Rho(1:NGrid)=2D0*Rho(1:NGrid)
       AbsGradRho2(1:NGrid)=4D0*AbsGradRho2(1:NGrid)
       CALL rks_xc_pw91(1,NGrid,Rho(1),AbsGradRho2(1),E(1),dEdRho(1),dEdGam(1),0D0,0D0,0D0)
    CASE(PURE_PW91_LYP)
       CALL PW91x_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One)
       CALL LYPc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One)
    CASE(PURE_B88_LYP)
       CALL B88x_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One)
       CALL LYPc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One)
    CASE(PURE_B88_PW91)
       Rho(1:NGrid)=2D0*Rho(1:NGrid)
       AbsGradRho2(1:NGrid)=4D0*AbsGradRho2(1:NGrid)
       ! Exchange
       CALL rks_x_b88(1,NGrid,Rho(1),AbsGradRho2(1),Buf(1),Buf(N1),Buf(N2),0D0,0D0,0D0)
       ! Copy the exchange in the mondo arrays !
       E     (1:NGrid)=Buf( 1: 1+NGrid-1)
       dEdRho(1:NGrid)=Buf(N1:N1+NGrid-1)
       dEdGam(1:NGrid)=Buf(N2:N2+NGrid-1)
       ! Correlation
       CALL rks_c_pw91(1,NGrid,Rho(1),AbsGradRho2(1),Buf(1),Buf(N1),Buf(N2),0D0,0D0,0D0)
       ! Add the correlation part !
       E     (1:NGrid)=Buf( 1: 1+NGrid-1)+     E(1:NGrid)
       dEdRho(1:NGrid)=Buf(N1:N1+NGrid-1)+dEdRho(1:NGrid)
       dEdGam(1:NGrid)=Buf(N2:N2+NGrid-1)+dEdGam(1:NGrid)
    CASE(PURE_HCTH93)
       Rho(1:NGrid)=2D0*Rho(1:NGrid)
       AbsGradRho2(1:NGrid)=4D0*AbsGradRho2(1:NGrid)
       CALL rks_xc_hcth(1,NGrid,Rho(1),AbsGradRho2(1),E(1),dEdRho(1),dEdGam(1),0D0,0D0,0D0)
    CASE(PURE_HCTH120)
       Rho(1:NGrid)=2D0*Rho(1:NGrid)
       AbsGradRho2(1:NGrid)=4D0*AbsGradRho2(1:NGrid)
       CALL rks_xc_hcth120(1,NGrid,Rho(1),AbsGradRho2(1),E(1),dEdRho(1),dEdGam(1),0D0,0D0,0D0)
    CASE(PURE_HCTH147)
       Rho(1:NGrid)=2D0*Rho(1:NGrid)
       AbsGradRho2(1:NGrid)=4D0*AbsGradRho2(1:NGrid)
       CALL rks_xc_hcth147(1,NGrid,Rho(1),AbsGradRho2(1),E(1),dEdRho(1),dEdGam(1),0D0,0D0,0D0)
    CASE(PURE_HCTH407)
       Rho(1:NGrid)=2D0*Rho(1:NGrid)
       AbsGradRho2(1:NGrid)=4D0*AbsGradRho2(1:NGrid)
       CALL rks_xc_hcth407(1,NGrid,Rho(1),AbsGradRho2(1),E(1),dEdRho(1),dEdGam(1),0D0,0D0,0D0)
    CASE(PURE_PBE_PBE)
       CALL PBEx_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One)
       CALL PBEc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One)
    CASE(PURE_XLYP)
       Ax1=0.722D0
       Ax2=0.347D0
       CALL XG04x_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam, One, Ax1, Ax2 )
       CALL LYPc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam, One )
    CASE(HYBRID_PBE0)
       ! PBE0 Model, Adamo and Barone: JCP 110, p.6158 (1999)
       ! E^{PBE0}_{xc}=E^{PBE}_{xc}+(1/4)(E^{HF}_x-E^{GGA}_x)
       !              =E^{PBE}_{c}+(3/4)E^{PBE}_x+(1/4)E^{HF}_x
       CALL PBEx_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,0.75D0)
       CALL PBEc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One)
    CASE(HYBRID_B3LYP_VWN3)
       !           Gaussianized B3LYP Model, uses VWN3 LSD correlation functional
       !           A. Becke: JCP 98, p.5648 (1993)
       !           R. H. Hertwig and W. Koch: CPL 268, p.345 (1997)
       !           P. J. Stevens et al, J. Phys. Chem. 98, p.11623 (1994)
       !           a0=0.2, ax=0.72, ac=0.81
       !           E^{B3LYP}_{xc}=(1-a0-ax)E^{LSD}_{x}+(1-ac)E^{VWN5}_c
       !                         +ax*E^{B88}_x+ac*E^{LYP}_c+a0*E^{HF}_x
       a0=0.20D0
       ax=0.72D0
       ac=0.81D0
       CALL XAx_ClSh  (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam, One-a0-ax )
       CALL B88x_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam, ax )
       CALL VWN3c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam, One-ac )
       CALL LYPc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam, ac )
    CASE(HYBRID_B3LYP_VWN5)
       !           Turbomolized B3LYP Model,
       !           A. Becke: JCP 98, p.5648 (1993), uses VWN5 LSD correlation functional
       !           R. H. Hertwig and W. Koch: CPL 268, p.345 (1997)
       !           a0=0.2, ax=0.72, ac=0.81
       !           E^{B3LYP}_{xc}=(1-a0-ax)E^{LSD}_{x}+(1-ac)E^{VWN5}_c
       !                         +ax*E^{B88}_x+ac*E^{LYP}_c+a0*E^{HF}_x
       a0=0.20D0
       ax=0.72D0
       ac=0.81D0
       CALL XAx_ClSh  (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam, One-a0-ax )
       CALL B88x_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam, ax )
       CALL VWN5c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam, One-ac )
       CALL LYPc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam, ac )
    CASE(HYBRID_X3LYP)
       !           Extended hybrid LYP functional
       !           Xu and Goddard, PNAS 101, p.2673 (2004)
       !           E^{X3LYP}_{xc}=(1-a0-ax)E^{LSD}_{x}+ac*E^{VWN3}_c
       !                         +ax*E^{XG04}_x+(1-ac)*E^{LYP}_c+a0*E^{HF}_x
       a0=0.218D0
       ax=0.709D0
       ac=0.129D0
       Ax1=0.765D0
       Ax2=0.235D0
       CALL XAx_ClSh  (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam, One-a0-ax )
       CALL XG04x_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam, ax, Ax1, Ax2)

       CALL VWN3c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam, ac )
       CALL LYPc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam, One-ac )
    CASE DEFAULT
       CALL Halt('Unknown functional requested ')
    END SELECT
  END SUBROUTINE ExcOnTheGrid_ClSh

  SUBROUTINE ExcOnTheGrid_OpnSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam)
    INTEGER                    :: NGrid,NSDen
    REAL(DOUBLE)               :: a0,ax,ac,Ax1,Ax2 , mxa,mxb,mna,mnb
    REAL(DOUBLE), DIMENSION(1:3*NGrid) :: Rho,AbsGradRho2
    REAL(DOUBLE), DIMENSION(1:3*NGrid) :: E,dEdRho,dEdGam
    REAL(DOUBLE), DIMENSION(1:6*NGrid) :: Buf
    INTEGER                    :: N1,N2,N3,N4,N5

    E=Zero
    dEdRho=Zero
    dEdGam=Zero
    !
    N1=  NGrid+1
    N2=2*NGrid+1
    N3=3*NGrid+1
    N4=4*NGrid+1
    N5=5*NGrid+1
    !
    DO I=1,2*NGrid
       IF(Rho(I)<1D-16)Rho(I)=Zero
    END DO
    !
    SELECT CASE(ModelChem)
    CASE(SD_EXCHANGE)
       CALL uks_x_lda(1,NGrid,Rho(1),Rho(N1),AbsGradRho2(1),AbsGradRho2(N1),AbsGradRho2(N2), &
            E(1),dEdRho(1),dEdRho(N1),dEdGam(1),dEdGam(N1),dEdGam(N2), &
            0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
            0D0,0D0,0D0,0D0,0D0,0D0,0D0)
    CASE(XA_EXCHANGE)
       IF(NSDen.NE.1)CALL Halt('This functional is not available for unrestricted! '//IntToChar(ModelChem))
    CASE(B88_EXCHANGE)
       IF(NSDen.NE.1)CALL Halt('This functional is not available for unrestricted! '//IntToChar(ModelChem))
    CASE(PW91_EXCHANGE)
       IF(NSDen.NE.1)CALL Halt('This functional is not available for unrestricted! '//IntToChar(ModelChem))
    CASE(PBE_EXCHANGE)
       IF(NSDen.NE.1)CALL Halt('This functional is not available for unrestricted! '//IntToChar(ModelChem))
    CASE(PURE_VWN3_LSD)
       IF(NSDen.NE.1)CALL Halt('This functional is not available for unrestricted! '//IntToChar(ModelChem))
    CASE(PURE_VWN5_LSD)
       IF(NSDen.NE.1)CALL Halt('This functional is not available for unrestricted! '//IntToChar(ModelChem))
    CASE(PURE_PW91_PW91)
       CALL uks_xc_pw91(1,NGrid,Rho(1),Rho(N1),AbsGradRho2(1),AbsGradRho2(N1),AbsGradRho2(N2), &
            E(1),dEdRho(1),dEdRho(N1),dEdGam(1),dEdGam(N1),dEdGam(N2), &
            0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
            0D0,0D0,0D0,0D0,0D0,0D0)
!!$
!!$               Buf=Zero
!!$               CALL uks_x_pw91(1,NGrid,Rho(1),Rho(N1),AbsGradRho2(1),AbsGradRho2(N1),AbsGradRho2(N2), &
!!$                       Buf(1),Buf(N1),Buf(N2),Buf(N3),Buf(N4),Buf(N5), &
!!$                       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
!!$                       0D0,0D0,0D0,0D0,0D0,0D0)
!!$               ! Copy the exchange in the mondo arrays !
!!$               E     ( 1: 1+NGrid-1)=Buf( 1: 1+NGrid-1)
!!$               dEdRho( 1: 1+NGrid-1)=Buf(N1:N1+NGrid-1)
!!$               dEdRho(N1:N1+NGrid-1)=Buf(N2:N2+NGrid-1)
!!$               dEdGam( 1: 1+NGrid-1)=Buf(N3:N3+NGrid-1)
!!$               dEdGam(N1:N1+NGrid-1)=Buf(N4:N4+NGrid-1)
!!$               dEdGam(N2:N2+NGrid-1)=Buf(N5:N5+NGrid-1)
!!$               ! Correlation
!!$               Buf=Zero
!!$               CALL uks_c_pw91(1,NGrid,Rho(1),Rho(N1),AbsGradRho2(1),AbsGradRho2(N1),AbsGradRho2(N2), &
!!$                       Buf(1),Buf(N1),Buf(N2),Buf(N3),Buf(N4),Buf(N5), &
!!$                       0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
!!$                       0D0,0D0,0D0,0D0,0D0,0D0)
!!$               ! Add the correlation part !
!!$               E     ( 1: 1+NGrid-1)=Buf( 1: 1+NGrid-1)+     E( 1: 1+NGrid-1)
!!$               dEdRho( 1: 1+NGrid-1)=Buf(N1:N1+NGrid-1)+dEdRho( 1: 1+NGrid-1)
!!$               dEdRho(N1:N1+NGrid-1)=Buf(N2:N2+NGrid-1)+dEdRho(N1:N1+NGrid-1)
!!$               dEdGam( 1: 1+NGrid-1)=Buf(N3:N3+NGrid-1)+dEdGam( 1: 1+NGrid-1)
!!$               dEdGam(N1:N1+NGrid-1)=Buf(N4:N4+NGrid-1)+dEdGam(N1:N1+NGrid-1)
!!$               dEdGam(N2:N2+NGrid-1)=Buf(N5:N5+NGrid-1)+dEdGam(N2:N2+NGrid-1)
 CASE(PURE_PW91_LYP)
    Buf=Zero
    CALL uks_x_pw91(1,NGrid,Rho(1),Rho(N1),AbsGradRho2(1),AbsGradRho2(N1),AbsGradRho2(N2), &
         Buf(1),Buf(N1),Buf(N2),Buf(N3),Buf(N4),Buf(N5), &
         0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
         0D0,0D0,0D0,0D0,0D0,0D0)
    ! Copy the exchange in the mondo arrays !
    E     ( 1: 1+NGrid-1)=Buf( 1: 1+NGrid-1)
    dEdRho( 1: 1+NGrid-1)=Buf(N1:N1+NGrid-1)
    dEdRho(N1:N1+NGrid-1)=Buf(N2:N2+NGrid-1)
    dEdGam( 1: 1+NGrid-1)=Buf(N3:N3+NGrid-1)
    dEdGam(N1:N1+NGrid-1)=Buf(N4:N4+NGrid-1)
    dEdGam(N2:N2+NGrid-1)=Buf(N5:N5+NGrid-1)
    ! Correlation
    Buf=Zero
    CALL uks_c_lyp(1,NGrid,Rho(1),Rho(N1),AbsGradRho2(1),AbsGradRho2(N1),AbsGradRho2(N2), &
         Buf(1),Buf(N1),Buf(N2),Buf(N3),Buf(N4),Buf(N5), &
         0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
         0D0,0D0,0D0,0D0,0D0,0D0)
    ! Add the correlation part !
    E     ( 1: 1+NGrid-1)=Buf( 1: 1+NGrid-1)+     E( 1: 1+NGrid-1)
    dEdRho( 1: 1+NGrid-1)=Buf(N1:N1+NGrid-1)+dEdRho( 1: 1+NGrid-1)
    dEdRho(N1:N1+NGrid-1)=Buf(N2:N2+NGrid-1)+dEdRho(N1:N1+NGrid-1)
    dEdGam( 1: 1+NGrid-1)=Buf(N3:N3+NGrid-1)+dEdGam( 1: 1+NGrid-1)
    dEdGam(N1:N1+NGrid-1)=Buf(N4:N4+NGrid-1)+dEdGam(N1:N1+NGrid-1)
    dEdGam(N2:N2+NGrid-1)=Buf(N5:N5+NGrid-1)+dEdGam(N2:N2+NGrid-1)
 CASE(PURE_B88_LYP)
    Buf=Zero
    CALL uks_x_b88(1,NGrid,Rho(1),Rho(N1),AbsGradRho2(1),AbsGradRho2(N1),AbsGradRho2(N2), &
         Buf(1),Buf(N1),Buf(N2),Buf(N3),Buf(N4),Buf(N5), &
         0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
         0D0,0D0,0D0,0D0,0D0,0D0)
    ! Copy the exchange in the mondo arrays !
    E     ( 1: 1+NGrid-1)=Buf( 1: 1+NGrid-1)
    dEdRho( 1: 1+NGrid-1)=Buf(N1:N1+NGrid-1)
    dEdRho(N1:N1+NGrid-1)=Buf(N2:N2+NGrid-1)
    dEdGam( 1: 1+NGrid-1)=Buf(N3:N3+NGrid-1)
    dEdGam(N1:N1+NGrid-1)=Buf(N4:N4+NGrid-1)
    dEdGam(N2:N2+NGrid-1)=Buf(N5:N5+NGrid-1)
    ! Correlation
    Buf=Zero
    CALL uks_c_lyp(1,NGrid,Rho(1),Rho(N1),AbsGradRho2(1),AbsGradRho2(N1),AbsGradRho2(N2), &
         Buf(1),Buf(N1),Buf(N2),Buf(N3),Buf(N4),Buf(N5), &
         0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
         0D0,0D0,0D0,0D0,0D0,0D0)
    ! Add the correlation part !
    E     ( 1: 1+NGrid-1)=Buf( 1: 1+NGrid-1)+     E( 1: 1+NGrid-1)
    dEdRho( 1: 1+NGrid-1)=Buf(N1:N1+NGrid-1)+dEdRho( 1: 1+NGrid-1)
    dEdRho(N1:N1+NGrid-1)=Buf(N2:N2+NGrid-1)+dEdRho(N1:N1+NGrid-1)
    dEdGam( 1: 1+NGrid-1)=Buf(N3:N3+NGrid-1)+dEdGam( 1: 1+NGrid-1)
    dEdGam(N1:N1+NGrid-1)=Buf(N4:N4+NGrid-1)+dEdGam(N1:N1+NGrid-1)
    dEdGam(N2:N2+NGrid-1)=Buf(N5:N5+NGrid-1)+dEdGam(N2:N2+NGrid-1)
 CASE(PURE_B88_PW91)
    Buf=Zero
    CALL uks_x_b88(1,NGrid,Rho(1),Rho(N1),AbsGradRho2(1),AbsGradRho2(N1),AbsGradRho2(N2), &
         Buf(1),Buf(N1),Buf(N2),Buf(N3),Buf(N4),Buf(N5), &
         0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
         0D0,0D0,0D0,0D0,0D0,0D0)
    ! Copy the exchange in the mondo arrays !
    E     ( 1: 1+NGrid-1)=Buf( 1: 1+NGrid-1)
    dEdRho( 1: 1+NGrid-1)=Buf(N1:N1+NGrid-1)
    dEdRho(N1:N1+NGrid-1)=Buf(N2:N2+NGrid-1)
    dEdGam( 1: 1+NGrid-1)=Buf(N3:N3+NGrid-1)
    dEdGam(N1:N1+NGrid-1)=Buf(N4:N4+NGrid-1)
    dEdGam(N2:N2+NGrid-1)=Buf(N5:N5+NGrid-1)
    ! Correlation
    Buf=Zero
    CALL uks_c_pw91(1,NGrid,Rho(1),Rho(N1),AbsGradRho2(1),AbsGradRho2(N1),AbsGradRho2(N2), &
         Buf(1),Buf(N1),Buf(N2),Buf(N3),Buf(N4),Buf(N5), &
         0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
         0D0,0D0,0D0,0D0,0D0,0D0)
    ! Add the correlation part !
    E     ( 1: 1+NGrid-1)=Buf( 1: 1+NGrid-1)+     E( 1: 1+NGrid-1)
    dEdRho( 1: 1+NGrid-1)=Buf(N1:N1+NGrid-1)+dEdRho( 1: 1+NGrid-1)
    dEdRho(N1:N1+NGrid-1)=Buf(N2:N2+NGrid-1)+dEdRho(N1:N1+NGrid-1)
    dEdGam( 1: 1+NGrid-1)=Buf(N3:N3+NGrid-1)+dEdGam( 1: 1+NGrid-1)
    dEdGam(N1:N1+NGrid-1)=Buf(N4:N4+NGrid-1)+dEdGam(N1:N1+NGrid-1)
    dEdGam(N2:N2+NGrid-1)=Buf(N5:N5+NGrid-1)+dEdGam(N2:N2+NGrid-1)
 CASE(PURE_HCTH93)
    CALL uks_xc_hcth(1,NGrid,Rho(1),Rho(N1),AbsGradRho2(1),AbsGradRho2(N1),AbsGradRho2(N2), &
         E(1),dEdRho(1),dEdRho(N1),dEdGam(1),dEdGam(N1),dEdGam(N2), &
         0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
         0D0,0D0,0D0,0D0,0D0,0D0)
 CASE(PURE_HCTH120)
    CALL uks_xc_hcth120(1,NGrid,Rho(1),Rho(N1),AbsGradRho2(1),AbsGradRho2(N1),AbsGradRho2(N2), &
         E(1),dEdRho(1),dEdRho(N1),dEdGam(1),dEdGam(N1),dEdGam(N2), &
         0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
         0D0,0D0,0D0,0D0,0D0,0D0)
 CASE(PURE_HCTH147)
    CALL uks_xc_hcth147(1,NGrid,Rho(1),Rho(N1),AbsGradRho2(1),AbsGradRho2(N1),AbsGradRho2(N2), &
         E(1),dEdRho(1),dEdRho(N1),dEdGam(1),dEdGam(N1),dEdGam(N2), &
         0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
         0D0,0D0,0D0,0D0,0D0,0D0)
 CASE(PURE_HCTH407)
    CALL uks_xc_hcth407(1,NGrid,Rho(1),Rho(N1),AbsGradRho2(1),AbsGradRho2(N1),AbsGradRho2(N2), &
         E(1),dEdRho(1),dEdRho(N1),dEdGam(1),dEdGam(N1),dEdGam(N2), &
         0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
         0D0,0D0,0D0,0D0,0D0,0D0)
 CASE(PURE_PBE_PBE)
    Buf=Zero
    CALL uks_x_pbe(1,NGrid,Rho(1),Rho(N1),AbsGradRho2(1),AbsGradRho2(N1),AbsGradRho2(N2), &
         Buf(1),Buf(N1),Buf(N2),Buf(N3),Buf(N4),Buf(N5), &
         0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
         0D0,0D0,0D0,0D0,0D0,0D0)
    ! Copy the exchange in the mondo arrays !
    E     ( 1: 1+NGrid-1)=Buf( 1: 1+NGrid-1)
    dEdRho( 1: 1+NGrid-1)=Buf(N1:N1+NGrid-1)
    dEdRho(N1:N1+NGrid-1)=Buf(N2:N2+NGrid-1)
    dEdGam( 1: 1+NGrid-1)=Buf(N3:N3+NGrid-1)
    dEdGam(N1:N1+NGrid-1)=Buf(N4:N4+NGrid-1)
    dEdGam(N2:N2+NGrid-1)=Buf(N5:N5+NGrid-1)
    ! Correlation
    Buf=Zero
    CALL uks_c_pbe(1,NGrid,Rho(1),Rho(N1),AbsGradRho2(1),AbsGradRho2(N1),AbsGradRho2(N2), &
         Buf(1),Buf(N1),Buf(N2),Buf(N3),Buf(N4),Buf(N5), &
         0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
         0D0,0D0,0D0,0D0,0D0,0D0)
    ! Add the correlation part !
    E     ( 1: 1+NGrid-1)=Buf( 1: 1+NGrid-1)+     E( 1: 1+NGrid-1)
    dEdRho( 1: 1+NGrid-1)=Buf(N1:N1+NGrid-1)+dEdRho( 1: 1+NGrid-1)
    dEdRho(N1:N1+NGrid-1)=Buf(N2:N2+NGrid-1)+dEdRho(N1:N1+NGrid-1)
    dEdGam( 1: 1+NGrid-1)=Buf(N3:N3+NGrid-1)+dEdGam( 1: 1+NGrid-1)
    dEdGam(N1:N1+NGrid-1)=Buf(N4:N4+NGrid-1)+dEdGam(N1:N1+NGrid-1)
    dEdGam(N2:N2+NGrid-1)=Buf(N5:N5+NGrid-1)+dEdGam(N2:N2+NGrid-1)
 CASE(PURE_XLYP)
    IF(NSDen.NE.1)CALL Halt('This functional is not available for unrestricted!'//IntToChar(ModelChem))
 CASE(HYBRID_PBE0)
    !           PBE0 Model, Adamo and Barone: JCP 110, p.6158 (1999)
    !           E^{PBE0}_{xc}=E^{PBE}_{xc}+(1/4)(E^{HF}_x-E^{GGA}_x)
    !                        =E^{PBE}_{c}+(3/4)E^{PBE}_x+(1/4)E^{HF}_x
    Buf=Zero
    CALL uks_x_pbe(1,NGrid,Rho(1),Rho(N1),AbsGradRho2(1),AbsGradRho2(N1),AbsGradRho2(N2), &
         Buf(1),Buf(N1),Buf(N2),Buf(N3),Buf(N4),Buf(N5), &
         0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
         0D0,0D0,0D0,0D0,0D0,0D0)
    ! Copy the exchange in the mondo arrays !
    E     ( 1: 1+NGrid-1)=Buf( 1: 1+NGrid-1)*0.75D0
    dEdRho( 1: 1+NGrid-1)=Buf(N1:N1+NGrid-1)*0.75D0
    dEdRho(N1:N1+NGrid-1)=Buf(N2:N2+NGrid-1)*0.75D0
    dEdGam( 1: 1+NGrid-1)=Buf(N3:N3+NGrid-1)*0.75D0
    dEdGam(N1:N1+NGrid-1)=Buf(N4:N4+NGrid-1)*0.75D0
    dEdGam(N2:N2+NGrid-1)=Buf(N5:N5+NGrid-1)*0.75D0
    ! Correlation
    Buf=Zero
    CALL uks_c_pbe(1,NGrid,Rho(1),Rho(N1),AbsGradRho2(1),AbsGradRho2(N1),AbsGradRho2(N2), &
         Buf(1),Buf(N1),Buf(N2),Buf(N3),Buf(N4),Buf(N5), &
         0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
         0D0,0D0,0D0,0D0,0D0,0D0)
    ! Add the correlation part !
    E     ( 1: 1+NGrid-1)=Buf( 1: 1+NGrid-1)+     E( 1: 1+NGrid-1)
    dEdRho( 1: 1+NGrid-1)=Buf(N1:N1+NGrid-1)+dEdRho( 1: 1+NGrid-1)
    dEdRho(N1:N1+NGrid-1)=Buf(N2:N2+NGrid-1)+dEdRho(N1:N1+NGrid-1)
    dEdGam( 1: 1+NGrid-1)=Buf(N3:N3+NGrid-1)+dEdGam( 1: 1+NGrid-1)
    dEdGam(N1:N1+NGrid-1)=Buf(N4:N4+NGrid-1)+dEdGam(N1:N1+NGrid-1)
    dEdGam(N2:N2+NGrid-1)=Buf(N5:N5+NGrid-1)+dEdGam(N2:N2+NGrid-1)
CASE(HYBRID_B3LYP_VWN3)
   !           Gaussianized B3LYP Model, uses VWN3 LSD correlation functional
   !           A. Becke: JCP 98, p.5648 (1993)
   !           R. H. Hertwig and W. Koch: CPL 268, p.345 (1997)
   !           P. J. Stevens et al, J. Phys. Chem. 98, p.11623 (1994)
   !           a0=0.2, ax=0.72, ac=0.81
   !           E^{B3LYP}_{xc}=(1-a0-ax)E^{LSD}_{x}+(1-ac)E^{VWN5}_c
   !                         +ax*E^{B88}_x+ac*E^{LYP}_c+a0*E^{HF}_x
 CALL uks_xc_b3lyp(1,NGrid,Rho(1),Rho(N1),AbsGradRho2(1),AbsGradRho2(N1),AbsGradRho2(N2), &
      E(1),dEdRho(1),dEdRho(N1),dEdGam(1),dEdGam(N1),dEdGam(N2), &
      0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0,0D0, &
      0D0,0D0,0D0,0D0,0D0,0D0)
CASE(HYBRID_B3LYP_VWN5)
 IF(NSDen.NE.1)CALL Halt('This functional is not available for unrestricted!'//IntToChar(ModelChem))
CASE(HYBRID_X3LYP)
 IF(NSDen.NE.1)CALL Halt('This functional is not available for unrestricted!'//IntToChar(ModelChem))
 !           Extended hybrid LYP functional
 !           Xu and Goddard, PNAS 101, p.2673 (2004)
 !           E^{X3LYP}_{xc}=(1-a0-ax)E^{LSD}_{x}+ac*E^{VWN3}_c
 !                         +ax*E^{XG04}_x+(1-ac)*E^{LYP}_c+a0*E^{HF}_x
CASE DEFAULT
 CALL Halt('Unknown functional requested ')
END SELECT
END SUBROUTINE ExcOnTheGrid_OpnSh

!====================================================================================
!     Closed shell version of X-alpha exchange
!====================================================================================
SUBROUTINE XAx_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,Scale)
INTEGER                        :: I,NGrid
REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2,E,dEdRho,dEdGam
REAL(DOUBLE)                   :: R,A,Scale,Ei,dEdRhoi,dEdGami,X,ASinh
ASinh(X)=LOG(X+SQRT(1.0D0+X*X))
!-------------------------------------------------------------------------------------
DO I=1,NGrid
 R=Rho(I)
 IF(NoNAN<R)THEN
    INCLUDE "XAx.Inc"
    E(I)=E(I)+Scale*Ei
    dEdRho(I)=dEdRho(I)+Scale*dEdRhoi
 ENDIF
ENDDO
END SUBROUTINE XAx_ClSh
!====================================================================================
!     Closed shell version of Slater-Dirac exchange
!====================================================================================
SUBROUTINE SDx_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,Scale)
INTEGER                        :: I,NGrid
REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2,E,dEdRho,dEdGam
REAL(DOUBLE)                   :: R,A,Scale,Ei,dEdRhoi,dEdGami,X,ASinh
ASinh(X)=LOG(X+SQRT(1.0D0+X*X))
!-------------------------------------------------------------------------------------
DO I=1,NGrid
 R=Rho(I)
 IF(NoNAN<R)THEN
    A=AbsGradRho2(I)
    INCLUDE "SDx.Inc"
    E(I)=E(I)+Scale*Ei
    dEdRho(I)=dEdRho(I)+Scale*dEdRhoi
    dEdGam(I)=dEdGam(I)+Scale*dEdGami
 ENDIF
ENDDO
END SUBROUTINE SDx_ClSh
!====================================================================================
!     Closed shell Becke 88 exchange engery functional
!     Physical Review A 38, p.3098 (1988), Equation (8).
!====================================================================================
SUBROUTINE B88x_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,Scale)
INTEGER                        :: I,NGrid
REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2,E,dEdRho,dEdGam
REAL(DOUBLE)                   :: R,A,Scale,Ei,dEdRhoi,dEdGami,X,ASinh
ASinh(X)=LOG(X+SQRT(1.0D0+X*X))
!-------------------------------------------------------------------------------------
DO I=1,NGrid
 R=Rho(I)
 IF(NoNAN<R)THEN
    A=AbsGradRho2(I)
    INCLUDE "B88x.Inc"
    E(I)=E(I)+Scale*Ei
    dEdRho(I)=dEdRho(I)+Scale*dEdRhoi
    dEdGam(I)=dEdGam(I)+Scale*dEdGami
 ENDIF
ENDDO
END SUBROUTINE B88x_ClSh
!====================================================================================
!     Closed shell version of the Perdew Wang 91 GGA-II
!     Physical Review B 46, p.6671 (1992)
!====================================================================================
SUBROUTINE PW91x_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,Scale)
INTEGER                        :: I,NGrid
REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2,E,dEdRho,dEdGam
REAL(DOUBLE)                   :: R,A,Scale,Ei,dEdRhoi,dEdGami,X,ASinh
ASinh(X)=LOG(X+SQRT(1.0D0+X*X))
!-------------------------------------------------------------------------------------
DO I=1,NGrid
 R=Rho(I)
 IF(NoNAN<R)THEN
    A=AbsGradRho2(I)
    INCLUDE "PW91x.Inc"
    E(I)=E(I)+Scale*Ei
    dEdRho(I)=dEdRho(I)+Scale*dEdRhoi
    dEdGam(I)=dEdGam(I)+Scale*dEdGami
 ENDIF
ENDDO
END SUBROUTINE PW91x_ClSh
!====================================================================================
!     Closed shell version of the Xu Goddard 04 GGA
!     PNAS, p.2673 (2004)
!====================================================================================
SUBROUTINE XG04x_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,Scale,Ax1,Ax2)
INTEGER                        :: I,NGrid
REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2,E,dEdRho,dEdGam
REAL(DOUBLE)                   :: R,A,Scale,Ei,dEdRhoi,dEdGami,X,ASinh,Ax1,Ax2
ASinh(X)=LOG(X+SQRT(1.0D0+X*X))
!-------------------------------------------------------------------------------------
DO I=1,NGrid
 R=Rho(I)
 IF(NoNAN<R)THEN
    A=AbsGradRho2(I)
    INCLUDE "XG04x.Inc"
    E(I)=E(I)+Scale*Ei
    dEdRho(I)=dEdRho(I)+Scale*dEdRhoi
    dEdGam(I)=dEdGam(I)+Scale*dEdGami
 ENDIF
ENDDO
END SUBROUTINE XG04x_ClSh
!====================================================================================
!     Closed shell version of the Perdew Burke Ernnzerhof exchange GGA
!     Physical Review Letters 77, p.3865 (1996)
!     E_x=\rho \epsilon^{PBE}_x
!====================================================================================
SUBROUTINE PBEx_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,Scale)
INTEGER                        :: I,NGrid
REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2,E,dEdRho,dEdGam
REAL(DOUBLE)                   :: R,A,Scale,Ei,dEdRhoi,dEdGami,X,ASinh
ASinh(X)=LOG(X+SQRT(1.0D0+X*X))
!-------------------------------------------------------------------------------------
DO I=1,NGrid
 R=Rho(I)
 IF(NoNAN<R)THEN
    A=AbsGradRho2(I)
    INCLUDE "PBEx.Inc"
    E(I)=E(I)+Scale*Ei
    dEdRho(I)=dEdRho(I)+Scale*dEdRhoi
    dEdGam(I)=dEdGam(I)+Scale*dEdGami
 ENDIF
ENDDO
END SUBROUTINE PBEx_ClSh
!====================================================================================
!     Closed shell version of the VWN3 LSD correlation functional
!====================================================================================
SUBROUTINE VWN3c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,Scale)
INTEGER                        :: I,NGrid
REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2,E,dEdRho,dEdGam
REAL(DOUBLE)                   :: R,A,Scale,Ei,dEdRhoi,dEdGami,X,ASinh
ASinh(X)=LOG(X+SQRT(1.0D0+X*X))
!-------------------------------------------------------------------------------------
DO I=1,NGrid
 R=Rho(I)
 IF(NoNAN<R)THEN
    A=AbsGradRho2(I)
    INCLUDE "VWN3c.Inc"
    E(I)=E(I)+Scale*Ei
    dEdRho(I)=dEdRho(I)+Scale*dEdRhoi
    dEdGam(I)=dEdGam(I)+Scale*dEdGami
 ENDIF
ENDDO
END SUBROUTINE VWN3c_ClSh
!====================================================================================
!     Closed shell version of the VWN5 LSD correlation functional
!====================================================================================
SUBROUTINE VWN5c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,Scale)
INTEGER                        :: I,NGrid
REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2,E,dEdRho,dEdGam
REAL(DOUBLE)                   :: R,A,Scale,Ei,dEdRhoi,dEdGami,X,ASinh
ASinh(X)=LOG(X+SQRT(1.0D0+X*X))
!-------------------------------------------------------------------------------------
DO I=1,NGrid
 R=Rho(I)
 IF(NoNAN<R)THEN
    A=AbsGradRho2(I)
    INCLUDE "VWN5c.Inc"
    E(I)=E(I)+Scale*Ei
    dEdRho(I)=dEdRho(I)+Scale*dEdRhoi
    dEdGam(I)=dEdGam(I)+Scale*dEdGami
 ENDIF
ENDDO
END SUBROUTINE VWN5c_ClSh
!====================================================================================
!     Closed shell Lee Yang Parr correlation engery functional
!     Physical Review B 37, p.785 (1988)
!     Rederived by Miehlich, Savin, Stoll, and Preuss
!     Chemical Physics Letters 157, P.200 (1989), Equation (2).
!====================================================================================
SUBROUTINE LYPc_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,Scale)
INTEGER                        :: I,NGrid
REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2,E,dEdRho,dEdGam
REAL(DOUBLE)                   :: R,A,Scale,Ei,dEdRhoi,dEdGami,X,ASinh
ASinh(X)=LOG(X+SQRT(1.0D0+X*X))
!-------------------------------------------------------------------------------------
DO I=1,NGrid
 R=Rho(I)
 IF(NoNAN<R)THEN
    A=AbsGradRho2(I)
    INCLUDE "LYPc.Inc"
    E(I)=E(I)+Scale*Ei
    dEdRho(I)=dEdRho(I)+Scale*dEdRhoi
    dEdGam(I)=dEdGam(I)+Scale*dEdGami
 ENDIF
ENDDO
END SUBROUTINE LYPc_ClSh
!====================================================================================
!     Closed shell version Perdew Wang LSD correlation functional, beyond RPA (p=1)
!     Perdew and Wang, Phys. Rev. B., p.13244 (1992)
!====================================================================================
SUBROUTINE PW91c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,Scale)
INTEGER                        :: I,NGrid
REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2,E,dEdRho,dEdGam
REAL(DOUBLE)                   :: R,A,Scale,Ei,dEdRhoi,dEdGami,X,ASinh
ASinh(X)=LOG(X+SQRT(1.0D0+X*X))
!-------------------------------------------------------------------------------------
DO I=1,NGrid
 R=Rho(I)
 IF(NoNAN<R)THEN
    A=(I)
    INCLUDE "PW91c.Inc"
    E(I)=E(I)+Scale*Ei
    dEdRho(I)=dEdRho(I)+Scale*dEdRhoi
    dEdGam(I)=dEdGam(I)+Scale*dEdGami
 ENDIF
ENDDO
END SUBROUTINE PW91c_ClSh
!====================================================================================
!     Closed shell Perdew Burke Ernnzerhof correlation GGA
!     Physical Review Letters 77, p.3865 (1996)
!====================================================================================
SUBROUTINE PBEc_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,Scale)
INTEGER                        :: I,NGrid
REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2,E,dEdRho,dEdGam
REAL(DOUBLE)                   :: R,A,Scale,Ei,dEdRhoi,dEdGami,X,ASinh
ASinh(X)=LOG(X+SQRT(1.0D0+X*X))
!------------------------------------------------------------------------------------
DO I=1,NGrid
 R=Rho(I)
 IF(NoNAN<R)THEN
    A=AbsGradRho2(I)
    INCLUDE "PBEc.Inc"
    E(I)=E(I)+Scale*Ei
    dEdRho(I)=dEdRho(I)+Scale*dEdRhoi
    dEdGam(I)=dEdGam(I)+Scale*dEdGami
 ENDIF
ENDDO
END SUBROUTINE PBEc_ClSh
END MODULE
