!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 2000, The University of California
!
MODULE Functionals
   USE DerivedTypes
   USE ProcessControl
   USE Parse
   IMPLICIT REAL(DOUBLE) (A-Z)
   INTERFACE ExcOnTheGrid
      MODULE PROCEDURE ExcOnTheGrid_ClSh
   END INTERFACE
!  Parsing keys for <Options.Model=>
   CHARACTER(LEN=9),  PARAMETER :: MODEL_OPTION='ModelChem'
!  Exchange only
   CHARACTER(LEN=2),  PARAMETER :: MODEL_ExactX       ='HF'
   CHARACTER(LEN=11), PARAMETER :: MODEL_SD           ='SlaterDirac'
   CHARACTER(LEN=6),  PARAMETER :: MODEL_XA           ='XAlpha'
   CHARACTER(LEN=4),  PARAMETER :: MODEL_B88x         ='B88x'
   CHARACTER(LEN=4),  PARAMETER :: MODEL_PBEx         ='PBEx'
   CHARACTER(LEN=5),  PARAMETER :: MODEL_PW91x        ='PW91x'
!  Pure exchange-correlation functionals
   CHARACTER(LEN=6),  PARAMETER :: MODEL_VWN3xc       ='VWN3xc'
   CHARACTER(LEN=6),  PARAMETER :: MODEL_VWN5xc       ='VWN5xc'
   CHARACTER(LEN=6),  PARAMETER :: MODEL_PW91xc       ='PW91xc'
   CHARACTER(LEN=6),  PARAMETER :: MODEL_BLYPxc       ='BLYPxc'
   CHARACTER(LEN=5),  PARAMETER :: MODEL_PBExc        ='PBExc'
!  Hybrid exchange-correlation functionals
   CHARACTER(LEN=10),  PARAMETER :: MODEL_B3LYP_VWN3xc='B3LYP/VWN3'
   CHARACTER(LEN=10),  PARAMETER :: MODEL_B3LYP_VWN5xc='B3LYP/VWN5'
   CHARACTER(LEN=10),  PARAMETER :: MODEL_B3LYP_PW91xc='B3LYP/PW91'
   CHARACTER(LEN=4),   PARAMETER :: MODEL_PBE0xc      ='PBE0'
!-----------------------------------------------------------------------------------------------------
!  Numerical keys for model chemistries
   INTEGER, PARAMETER :: EXACT_EXCHANGE   =03045805 ! Use exact Hartree-Fock exchange
   INTEGER, PARAMETER :: HAS_DFT          =10000000 ! Any key above this has DFT
   INTEGER, PARAMETER :: HYBRID_B3LYP_VWN3=10268305 ! Use Beckes B3LYP model with VWN3 LSD correlation
   INTEGER, PARAMETER :: HYBRID_B3LYP_VWN5=10429495 ! Use Beckes B3LYP model with VWN5 LSD correlation
   INTEGER, PARAMETER :: HYBRID_B3LYP_PW91=10685608 ! Use Beckes B3LYP model with PW91 LSD correlation
   INTEGER, PARAMETER :: HYBRID_PBE0      =14506981 ! Use Adamo and Barones PBE0 model
   INTEGER, PARAMETER :: HAS_HF           =20000000 ! Any key below this has exact exchange   
   INTEGER, PARAMETER :: SD_EXCHANGE      =32098243 ! Slater-Dirac Uniform Electron Gas exchange
   INTEGER, PARAMETER :: XA_EXCHANGE      =33890824 ! X-Alpha LDA exchange
   INTEGER, PARAMETER :: B88_EXCHANGE     =34065852 ! Becke 88 exchange
   INTEGER, PARAMETER :: PW91_EXCHANGE    =38458583 ! Perdew Wang 91 exchange
   INTEGER, PARAMETER :: PBE_EXCHANGE     =36792335 ! Perdew Burke Enzerhoff exchange
   INTEGER, PARAMETER :: PURE_VWN3_LSD    =40022304 ! Slater exchange and VWN3 LSDA correlation
   INTEGER, PARAMETER :: PURE_VWN5_LSD    =40038402 ! Slater exchange and VWN5 LSDA correlation
   INTEGER, PARAMETER :: PURE_PW91_LSD    =40230174 ! Slater exchange and PW92 beyond RPA LSDA correlation
   INTEGER, PARAMETER :: PURE_PBE_GGA     =42509434 ! Perdew Burke Enzerhoff GGA exchange-correlation
   INTEGER, PARAMETER :: PURE_BLYP_GGA    =43083038 ! B88 exchange with LYP correlation
!-----------------------------------------------------------------------------------------------------
!  Avoid under and over flows
   REAL(DOUBLE), PARAMETER :: NoNAN=1.D-30
!  Global variable set to the current model chemistry
   INTEGER :: ModelChem
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
         CASE(PURE_PW91_LSD)
            Name='PW91x/PW91c'
         CASE(PURE_BLYP_GGA)
            Name='B88x/LYPc'
         CASE(PURE_PBE_GGA)
            Name='PBEx/PBEc'  
         CASE(HYBRID_PBE0)
            Name='PBE0'
         CASE(HYBRID_B3LYP_VWN3)
            Name='B3LYP(VWN3)' 
         CASE(HYBRID_B3LYP_VWN5)
            Name='B3LYP(VWN5)'
         CASE(HYBRID_B3LYP_PW91)
            Name='B3LYP(PW91)'
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
                Key==HYBRID_B3LYP_VWN5.OR. &
                Key==HYBRID_B3LYP_PW91)THEN
            Scalar=0.20D0
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
         IF(Key<Has_HF)THEN
            HasHF=.TRUE.
         ELSE
            HasHF=.FALSE.
         ENDIF
      END FUNCTION HasHF
!====================================================================================
!
!====================================================================================
      SUBROUTINE ExcOnTheGrid_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam) 
         INTEGER                        :: NGrid
         REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2
         REAL(DOUBLE), DIMENSION(NGrid) :: E,dEdRho,dEdGam
!---------------------------------------------------------------------------------------
         E=Zero 
         dEdRho=Zero 
         dEdGam=Zero
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
         CASE(PURE_PW91_LSD)
            CALL Halt(' Functional not yet validated ')
            CALL PW91x_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One)
            CALL PW91c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One)
         CASE(PURE_BLYP_GGA)
            CALL B88x_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One) 
            CALL LYPc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One) 
         CASE(PURE_PBE_GGA)
            CALL Halt(' Functional not yet validated ')
            CALL PBEx_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One) 
            CALL PBEc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One) 
         CASE(HYBRID_PBE0)
            CALL Halt(' Functional not yet validated ')
!           PBE0 Model, Adamo and Barone: JCP 110, p.6158 (1999)
!           E^{PBE0}_{xc}=E^{PBE}_{xc}+(1/4)(E^{HF}_x-E^{GGA}_x)
!                        =E^{PBE}_{c}+(3/4)E^{PBE}_x+(1/4)E^{HF}_x
            CALL PBEx_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,0.75D0) 
            CALL PBEc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,One)
         CASE(HYBRID_B3LYP_VWN3)
!           Gaussianized B3LYP Model, uses VWN3 LSD correlation functional
!           A. Becke: JCP 98, p.5648 (1993)
!           R. H. Hertwig and W. Koch: CPL 268, p.345 (1997)
!           P. J. Stevens et al, J. Phys. Chem. 98, p.11623 (1994)
!           a0=0.2, ax=0.72, ac=0.81
!           E^{B3LYP}_{xc}=(1-a0-ax)E^{LSD}_{x}+(1-ac)E^{VWN3}_c 
!                         +ax*E^{B88}_x+ac*E^{LYP}_c+a0*E^{HF}_x                 
            CALL XAx_ClSh  (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,0.08D0)
            CALL B88x_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,0.72D0)
            CALL VWN3c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,0.19D0)
            CALL LYPc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,0.81D0)
         CASE(HYBRID_B3LYP_VWN5)
!           Turbomolized B3LYP Model, 
!           A. Becke: JCP 98, p.5648 (1993), uses VWN5 LSD correlation functional
!           R. H. Hertwig and W. Koch: CPL 268, p.345 (1997)
!           a0=0.2, ax=0.72, ac=0.81
!           E^{B3LYP}_{xc}=(1-a0-ax)E^{LSD}_{x}+(1-ac)E^{VWN5}_c 
!                         +ax*E^{B88}_x+ac*E^{LYP}_c+a0*E^{HF}_x                 
            CALL XAx_ClSh  (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,0.08D0)
            CALL B88x_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,0.72D0)
            CALL VWN5c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,0.19D0)
            CALL LYPc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,0.81D0)
         CASE(HYBRID_B3LYP_PW91)
!           A. Becke: JCP 98, p.5648 (1993), uses PW91 LSD correlation functional
!           a0=0.2, ax=0.72, ac=0.81
!           E^{B3LYP}_{xc}=(1-a0-ax)E^{LSD}_{x}+(1-ac)E^{VWN5}_c 
!                         +ax*E^{B88}_x+ac*E^{LYP}_c+a0*E^{HF}_x                 
            CALL XAx_ClSh  (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,0.08D0)
            CALL B88x_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,0.72D0)
            CALL PW91c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,0.19D0)
            CALL LYPc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdGam,0.81D0)
         CASE DEFAULT
            CALL Halt('Unknown functional requested ') 
         END SELECT
      END SUBROUTINE ExcOnTheGrid_ClSh
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
