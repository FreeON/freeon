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
   CHARACTER(LEN=5),  PARAMETER :: MODEL_OPTION='ModelChem'
!  Exchange only
   CHARACTER(LEN=2),  PARAMETER :: MODEL_ExactX       ='HF'
   CHARACTER(LEN=6),  PARAMETER :: MODEL_SlaterX      ='Slater'
   CHARACTER(LEN=4),  PARAMETER :: MODEL_B88x         ='B88x'
   CHARACTER(LEN=4),  PARAMETER :: MODEL_PBEx         ='PBEx'
   CHARACTER(LEN=5),  PARAMETER :: MODEL_PW91x        ='PW91x'
!  Pure exchange-correlation functionals
   CHARACTER(LEN=5),  PARAMETER :: MODEL_SVWN3xc      ='SVWN3'
   CHARACTER(LEN=5),  PARAMETER :: MODEL_SVWN5xc      ='SVWN5'
   CHARACTER(LEN=5),  PARAMETER :: MODEL_SPW92xc      ='SPW92'
   CHARACTER(LEN=4),  PARAMETER :: MODEL_BLYPxc       ='BLYP'
   CHARACTER(LEN=3),  PARAMETER :: MODEL_PBExc        ='PBE'
!  Hybrid exchange-correlation functionals
   CHARACTER(LEN=10),  PARAMETER :: MODEL_B3LYP_VWN3xc='B3LYP/VWN3'
   CHARACTER(LEN=10),  PARAMETER :: MODEL_B3LYP_VWN5xc='B3LYP/VWN5'
   CHARACTER(LEN=10),  PARAMETER :: MODEL_B3LYP_PW92xc='B3LYP/PW92'
   CHARACTER(LEN=4),   PARAMETER :: MODEL_PBE0xc       ='PBE0'
!  Numerical keys for model chemistries
   INTEGER, PARAMETER :: EXACT_EXCHANGE   =03045805 ! Use exact Hartree-Fock exchange
   INTEGER, PARAMETER :: HAS_DFT          =10000000 ! Any key above this has DFT
   INTEGER, PARAMETER :: HYBRID_B3LYP_VWN3=10268305 ! Use Beckes B3LYP model with VWN3 LSD correlation
   INTEGER, PARAMETER :: HYBRID_B3LYP_VWN5=10429495 ! Use Beckes B3LYP model with VWN5 LSD correlation
   INTEGER, PARAMETER :: HYBRID_B3LYP_PW92=10685608 ! Use Beckes B3LYP model with PW92 LSD correlation
   INTEGER, PARAMETER :: HYBRID_PBE0      =14506981 ! Use Adamo and Barones PBE0 model
   INTEGER, PARAMETER :: HAS_HF           =20000000 ! Any key below this has exact exchange   
   INTEGER, PARAMETER :: LDA_EXCHANGE     =32098243 ! Slater LDA exchange
   INTEGER, PARAMETER :: B88_EXCHANGE     =34065852 ! Becke 88 exchange
   INTEGER, PARAMETER :: PW91_EXCHANGE    =38458583 ! Perdew Wang 91 exchange
   INTEGER, PARAMETER :: PBE_EXCHANGE     =36792335 ! Perdew Burke Enzerhoff exchange
   INTEGER, PARAMETER :: PURE_SVWN3_LSD   =40022304 ! Slater exchange and VWN3 LSDA correlation
   INTEGER, PARAMETER :: PURE_SVWN5_LSD   =40038402 ! Slater exchange and VWN5 LSDA correlation
   INTEGER, PARAMETER :: PURE_SPW92_LSD   =40230174 ! Slater exchange and PW92 beyond RPA LSDA correlation
   INTEGER, PARAMETER :: PURE_PBE_GGA     =42509434 ! Perdew Burke Enzerhoff GGA exchange-correlation
   INTEGER, PARAMETER :: PURE_BLYP_GGA    =43083038 ! B88 exchange with LYP correlation
!  Avoid under and over flows
   REAL(DOUBLE), PARAMETER :: NoNAN=1.D-10
!  Global variable set to the current model chemistry
   INTEGER :: ModelChem
   CONTAINS
!====================================================================================
!
!====================================================================================
      FUNCTION FunctionalName(Key) RESULT(Name)
         INTEGER            :: Key
         CHARACTER(LEN=30)  :: Name
         SELECT CASE(Key)
         CASE(EXACT_EXCHANGE)
            Name='HF'
         CASE(LDA_EXCHANGE)
            Name='S'
         CASE(B88_EXCHANGE)
            Name='B88'
         CASE(PW91_EXCHANGE)
            Name='PW91x'
         CASE(PBE_EXCHANGE)
            Name='PBEx'
         CASE(PURE_SVWN3_LSD)
            Name='SVWN3'
         CASE(PURE_SVWN5_LSD)
            Name='SVWN5'
         CASE(PURE_SPW92_LSD)
            Name='SPW92'
         CASE(PURE_BLYP_GGA)
            Name='BLYP'
         CASE(PURE_PBE_GGA)
            Name='PBE'  
         CASE(HYBRID_PBE0)
            Name='PBE0'
         CASE(HYBRID_B3LYP_VWN3)
            Name='B3LYP(VWN3)' 
         CASE(HYBRID_B3LYP_VWN5)
            Name='B3LYP(VWN5)'
         CASE(HYBRID_B3LYP_PW92)
            Name='B3LYP(PW92)'
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
                Key==HYBRID_B3LYP_PW92)THEN
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
      SUBROUTINE ExcOnTheGrid_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2) 
         INTEGER                        :: NGrid
         REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2
         REAL(DOUBLE), DIMENSION(NGrid) :: E,dEdRho,dEdAbsGradRho2
         E=Zero 
         dEdRho=Zero 
         dEdAbsGradRho2=Zero
         SELECT CASE(ModelChem)         
         CASE(LDA_EXCHANGE)
            CALL LSDx_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,One)
         CASE(B88_EXCHANGE)
            CALL B88x_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,One) 
         CASE(PW91_EXCHANGE)
            CALL PW91x_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,One) 
         CASE(PBE_EXCHANGE)
            CALL PBEx_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,One) 
         CASE(PURE_SVWN3_LSD)
            CALL LSDx_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,One)
            CALL VWN3c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,One)
         CASE(PURE_SVWN5_LSD)
            CALL LSDx_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,One)
            CALL VWN5c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,One)
         CASE(PURE_SPW92_LSD)
            CALL LSDx_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,One)
            CALL PW92c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,One)
         CASE(PURE_BLYP_GGA)
            CALL B88x_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,One) 
            CALL LYPc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,One) 
         CASE(PURE_PBE_GGA)
            CALL PBEx_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,One) 
            CALL PBEc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,One) 
         CASE(HYBRID_PBE0)
!           PBE0 Model, Adamo and Barone: JCP 110, p.6158 (1999)
!           E^{PBE0}_{xc}=E^{PBE}_{xc}+(1/4)(E^{HF}_x-E^{GGA}_x)
!                        =E^{PBE}_{c}+(3/4)E^{PBE}_x+(1/4)E^{HF}_x
            CALL PBEx_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,0.75D0) 
            CALL PBEc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,One)
         CASE(HYBRID_B3LYP_VWN3)
!           Gaussianized B3LYP Model, uses VWN5 LSD correlation functional
!           A. Becke: JCP 98, p.5648 (1993)
!           R. H. Hertwig and W. Koch: CPL 268, p.345 (1997)
!           P. J. Stevens et al, J. Phys. Chem. 98, p.11623 (1994)
!           a0=0.2, ax=0.72, ac=0.81
!           E^{B3LYP}_{xc}=(1-a0-ax)E^{LSD}_{x}+(1-ac)E^{VWN3}_c 
!                         +ax*E^{B88}_x+ac*E^{LYP}_c+a0*E^{HF}_x                 
            CALL LSDx_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,0.08D0)
            CALL B88x_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,0.72D0)
            CALL VWN3c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,0.19D0)
            CALL LYPc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,0.81D0)
         CASE(HYBRID_B3LYP_VWN5)


! STO-3G
!    9     -75.2665126291      5.34E-09  000000  Convergence criterion met
! ----------------------------------------------------
! One-Electron    Energy =  -122.5946137216
! Total Coulomb   Energy =    47.4596205222
! Alpha Exchange  Energy =    -0.9119178430
! Beta  Exchange  Energy =    -0.9119178430
! DFT   Exchange  Energy =    -7.2024742957
! DFT Correlation Energy =    -0.4032173523

!           Turbomolized B3LYP Model, 
!           A. Becke: JCP 98, p.5648 (1993), uses VWN5 LSD correlation functional
!           R. H. Hertwig and W. Koch: CPL 268, p.345 (1997)
!           a0=0.2, ax=0.72, ac=0.81
!           E^{B3LYP}_{xc}=(1-a0-ax)E^{LSD}_{x}+(1-ac)E^{VWN5}_c 
!                         +ax*E^{B88}_x+ac*E^{LYP}_c+a0*E^{HF}_x                 
            CALL LSDx_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,0.08D0)
            CALL B88x_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,0.72D0)
            CALL VWN5c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,0.19D0)
            CALL LYPc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,0.81D0)
         CASE(HYBRID_B3LYP_PW92)
!           Mondoized B3LYP Model, 
!           A. Becke: JCP 98, p.5648 (1993), uses PW92 LSD correlation functional
!           a0=0.2, ax=0.72, ac=0.81
!           E^{B3LYP}_{xc}=(1-a0-ax)E^{LSD}_{x}+(1-ac)E^{VWN5}_c 
!                         +ax*E^{B88}_x+ac*E^{LYP}_c+a0*E^{HF}_x                 
            CALL LSDx_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,0.08D0)
            CALL B88x_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,0.72D0)
            CALL PW92c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,0.19D0)
            CALL LYPc_ClSh (NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,0.81D0)
         CASE DEFAULT
            CALL Halt('Unknown functional requested ') 
         END SELECT
      END SUBROUTINE ExcOnTheGrid_ClSh
!====================================================================================
!     Closed shell version of Slater exchange
!====================================================================================
      SUBROUTINE LSDx_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,Scale) 
         INTEGER                        :: NGrid
         REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2
         REAL(DOUBLE), DIMENSION(NGrid) :: E,dEdRho,dEdAbsGradRho2
         REAL(DOUBLE)                   :: R,A,Scale
         INTEGER                        :: I
!-------------------------------------------------------------------------------------
         DO I=1,NGrid
            R=Rho(I)
            IF(NoNAN<R)THEN
               R=Two*R
               R1=R**(1.d0/3.d0)
               E(I)=E(I)+Scale*(-7.38558766382022d-1*R*R1) 
               dEdRho(I)=dEdRho(I)+Scale*( -9.84745021842697d-1*R1)
            ENDIF
         ENDDO
      END SUBROUTINE LSDx_ClSh
!====================================================================================
!     Closed shell version Perdew Wang LSD correlation functional, beyond RPA (p=1) 
!     Perdew and Wang, Phys. Rev. B., p.13244 (1992)
!====================================================================================
      SUBROUTINE PW92c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,Scale) 
         INTEGER                        :: NGrid
         REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2
         REAL(DOUBLE), DIMENSION(NGrid) :: E,dEdRho,dEdAbsGradRho2
         REAL(DOUBLE)                   :: R,A,Scale
         INTEGER                        :: I
!-------------------------------------------------------------------------------------
         DO I=1,NGrid
            R=Rho(I)
            IF(NoNAN<R)THEN
               R=Two*R
               R1=R**(1.d0/3.d0)
               R2=R1**2
               R3=R**(1.d0/6.d0)
               R4=1.33793879983596d-1*R3
               R5=3.72010139328451d-1*R1
               R6=sqrt(R)
               R7=R4+R5+R6
               R8=3.17089566751485d-2+R7
               R9=1/R8
               R10=R2*R9
               R11=2.68814775092151d0*R10
               R12=1.d0+R11
               R13=log(R12)
               R14=1/R12
               R15=R8**2
               R16=1/R15
               R17=R3**2
               R18=R17**2
               E(I)=E(I)+Scale*((-8.24331979256531d-3-6.21814d-2*R1)*R13*R2)
               RTmp1=-6.21814d-2*R13-(5.49554652837688d-3*R13)/R1
               RTmp2=2.07275109683497d-2*R*R14*R16+2.74782332684713d-3*R14*R16*R2
               RTmp3=8.35763952795754d-2*R*R14*R16*R3+1.48069675134213d-2*R14*R16*R18*R3
               RTmp4=4.94128930292024d-4*R14*R16*R6
               RTmp5=-1.47728410403408d-2*R1*R14*R9-1.11435193706101d-1*R14*R2*R9
               dEdRho(I)=dEdRho(I)+Scale*(RTmp1+RTmp2+RTmp3+RTmp4+RTmp5)
            ENDIF
         ENDDO
      END SUBROUTINE PW92c_ClSh
!====================================================================================
!     Closed shell version of the VWN3 LSD correlation functional
!====================================================================================
      SUBROUTINE VWN3c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,Scale) 
         INTEGER                        :: NGrid
         REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2
         REAL(DOUBLE), DIMENSION(NGrid) :: E,dEdRho,dEdAbsGradRho2
         REAL(DOUBLE)                   :: R,A,Scale
         INTEGER                        :: I
!-------------------------------------------------------------------------------------
         DO I=1,NGrid
            R=Rho(I)
            IF(NoNAN<R)THEN
               R=Two*R

        R1=R**(1.d0/6.d0)
        R2=1.3072d1*R1
        R3=1.57524663579949d0+R2
        R4=1.0D0/R3
        R5=R1*R4
        R6=4.48998886414974d-2*R5
        R7=atan(R6)
        R8=2.05219729377415d1*R7
        R9=1.0D0/R1
        R10=7.87623317899743d-1*R9
        R11=4.09286d-1+R10
        R12=log(R11)
        R13=8.86274753549907d-3*R12
        R14=R**(1.d0/3.d0)
        R15=1.0D0/R14
        R16=6.203504908994d-1*R15
        R17=1.02958120115854d1*R9
        R18=R16+R17
        R19=4.27198d1+R18
        R20=log(R19)
        R21=-3.55220737677495d-2*R20
        R22=log(R)
        R23=-1.03635666666667d-2*R22
        R24=1.0D0/R19
        R25=R3**2
        R26=1.0D0/R25
        R27=R14*R26
        R28=2.01600000001887d-3*R27
        R29=1.d0+R28
        R30=1.0D0/R29
        E(I)=E(I)+Scale*(R*(-1.48448968239848d-2+R13+R21+R23+R8))
        RTmp1=-2.52084634906515d-2+R13+R21+R23+7.34537863319604d-3*R15*R24
        RTmp2=-2.0074981940802d0*R14*R26*R30+1.53572383268069d-1*R1*R30*R4
        RTmp3=R8-(1.16341776993626d-3*R9)/R11+6.095476562907d-2*R24*R9
        dEdRho(I)=dEdRho(I)+Scale*(RTmp1+RTmp2+RTmp3)
        
            ENDIF
         ENDDO
      END SUBROUTINE VWN3c_ClSh
!====================================================================================
!     Closed shell version of the VWN5 LSD correlation functional
!====================================================================================
      SUBROUTINE VWN5c_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,Scale) 
         INTEGER                        :: NGrid
         REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2
         REAL(DOUBLE), DIMENSION(NGrid) :: E,dEdRho,dEdAbsGradRho2
         REAL(DOUBLE)                   :: R,A,Scale
         INTEGER                        :: I
!-------------------------------------------------------------------------------------
         DO I=1,NGrid
            R=Rho(I)
            IF(NoNAN<R)THEN
               R=Two*R


        R1=R**(1.d0/6.d0)
        R2=3.72744d0*R1
        R3=1.57524663579949d0+R2
        R4=1.0D0/R3
        R5=R1*R4
        R6=6.15199081975908d0*R5
        R7=atan(R6)
        R8=3.8783294878113d-2*R7
        R9=1.0D0/R1
        R10=7.87623317899743d-1*R9
        R11=1.0498d-1+R10
        R12=log(R11)
        R13=1.93804554230887d-3*R12
        R14=R**(1.d0/3.d0)
        R15=1.0D0/R14
        R16=6.203504908994d-1*R15
        R17=2.93581866007222d0*R9
        R18=R16+R17
        R19=1.29352d1+R18
        R20=log(R19)
        R21=-3.20597227711544d-2*R20
        R22=log(R)
        R23=-1.03635666666667d-2*R22
        R24=1.0D0/R19
        R25=R3**2
        R26=1.0D0/R25
        R27=R14*R26
        R28=3.78469910464d1*R27
        R29=1.d0+R28
        R30=1.0D0/R29
        E(I)=E(I)+Scale*(R*(-1.48448968239848d-2+R13+R21+R23+R8))
        RTmp1=-2.52084634906515d-2+R13+R21+R23+6.62942158639478d-3*R15*R24
        RTmp2=-1.48224431058922d-1*R14*R26*R30+3.97657456750268d-2*R1*R30*R4
        RTmp3=R8-(2.54408310045687d-4*R9)/R11+1.56869220580496d-2*R24*R9
        dEdRho(I)=dEdRho(I)+Scale*(RTmp1+RTmp2+RTmp3)
            ENDIF
         ENDDO
      END SUBROUTINE VWN5c_ClSh
!====================================================================================
!     Closed shell version of the Perdew Burke Ernnzerhof exchange GGA
!     Physical Review Letters 77, p.3865 (1996)
!     E_x=\rho \epsilon^{PBE}_x
!====================================================================================
      SUBROUTINE PBEx_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,Scale) 
         INTEGER                        :: NGrid
         REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2
         REAL(DOUBLE), DIMENSION(NGrid) :: E,dEdRho,dEdAbsGradRho2
         REAL(DOUBLE)                   :: R,A,Scale
         INTEGER                        :: I
!-------------------------------------------------------------------------------------
         DO I=1,NGrid
            R=Rho(I)
            IF(NoNAN<R)THEN
               R=Two*R
               A=Four*AbsGradRho2(I)
               A=AbsGradRho2(I)
               R1=R**(1.d0/3.d0)
               R2=R*R1
               R3=R**2
               R4=R1**2
               R5=R3*R4
               R6=7.13182658760049d-3*A
               R7=R5+R6
               R8=R7**2
               R9=1.D0/R8
               R10=A**2
               R11=R3**2
               E(I)=E(I)+Scale*( &
                   -1.33236001455317d0*R2+(5.93801248171146d-1*R2)/(1.d0+ &
                   (7.13182658760049d-3*A)/R5))
               dEdRho(I)=dEdRho(I)+Scale*( &
                   (-9.03570152478593d-5*R1*R10-8.39954475162682d-3*A*R*R3 &
                   -9.84745021842697d-1*R*R11*R4)*R9)
               dEdAbsGradRho2(I)=dEdAbsGradRho2(I)+Scale*(  &
                   -4.23488752945734d-3*R11*R9)
            ENDIF
         ENDDO
      END SUBROUTINE PBEx_ClSh
!====================================================================================
!     Closed shell version of the Perdew Wang 91 GGA-II
!     Physical Review B 46, p.6671 (1992)
!====================================================================================
      SUBROUTINE PW91x_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,Scale)
         INTEGER                        :: NGrid
         REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2
         REAL(DOUBLE), DIMENSION(NGrid) :: E,dEdRho,dEdAbsGradRho2
         REAL(DOUBLE)                   :: R,A,ASinh,DASinh,X,Scale
         INTEGER                        :: I
!-------------------------------------------------------------------------------------
         DO I=1,NGrid
            R=Rho(I)
            IF(NoNAN<R)THEN
               R=Two*R
               A=Four*AbsGradRho2(I)
               R1=R**(1.d0/3.d0)
               R2=R*R1
               R3=1.D0/R2
               R4=A**2
               R5=R**2
               R6=R5**2
               R7=R*R6
               R8=R1*R7
               R9=1.D0/R8
               R10=R4*R9
               R11=2.72926271249799d-6*R10
               R12=sqrt(A)
               R13=R1**2
               R14=R13*R5
               R15=1.D0/R14
               R16=A*R15
               R17=1.58741971281381d0*R16
               R18=1.d0+R17
               R19=sqrt(R18)
               R20=R12*R3
               R21=1.2599284554346d0*R20
               R22=R19+R21
               R23=log(R22)
               R24=R12*R23*R3
               R25=3.17503393029564d-2*R24
               R26=R11+R25
               R27=1.d0+R26
               R28=1.D0/R27
               R29=-2.61211729852336d0*R16
               R30=exp(R29)
               R31=R4**2
               R32=1.D0/R19
               R33=R*R12
               R34=1.2599284554346d0*R33
               R35=R1*R5
               R36=R19*R35
               R37=R34+R36
               R38=1.D0/R37
               R39=2.72926271249799d-6*R4
               R40=R12*R23*R6
               R41=3.17503393029564d-2*R40
               R42=R39+R41+R8
               R43=R42**2
               R44=1.D0/R43
               R45=A*R4
               R46=R12*R45
               R47=R12*R4
               R48=R5*R6
               R49=R13*R48
               R50=R6**2
               R51=A*R12
               R52=R*R50
               R53=R1*R52
               R54=R5*R50
               R55=R13*R54
               R56=R50*R6
               R57=R23**2
               R58=1.2599284554346d0*R12
               R59=R19*R2
               R60=R58+R59
               R61=1.D0/R60
               E(I)=E(I)+Scale*( &
                   -7.38558766382022d-1*R2*R28-2.34494914278021d-2*R12*R23*R28 &
                   -5.29180144160952d-3*A*R28*R3+2.90923681150097d-3*A*R28*R3*R30)
               dEdRho(I)=dEdRho(I)+Scale*( &
                    6.96837811561903d-8*R12*R2*R30*R31*R38*R44-9.17066106037242d-8 &
                   *R14*R31*R32*R38*R44+8.77965787609615d-8*A*R30*R31*R32*R38      &
                   *R44+1.05724633077259d-7*R14*R30*R31*R32*R38*R44                &
                   -1.69310941515903d-5*R38*R44*R47*R49+2.55320899805982d-2*R30    &
                   *R38*R44*R47*R49-3.55616920479275d-4*R32*R38*R44*R47*R49        &
                   -3.41332387098355d-7*R23*R32*R38*R44*R47*R49                    &
                   +1.95505036851165d-4*R30*R32*R38*R44*R47*R49                    &
                   +6.434115496804d-4*R23*R30*R32*R38*R44*R47*R49                  &
                   -2.8225167781976d-4*R38*R4*R44*R50+1.55171538517025d-4          &
                   *R30*R38*R4*R44*R50+1.11869750935962d-2*R32*R38*R4*R44          &
                   *R50+1.41071409837667d-2*R30*R32*R38*R4*R44*R50                 &
                   +8.88972162239157d-3*R38*R44*R51*R53-4.88724032321053d-3        &
                   *R30*R38*R44*R51*R53-9.92644931945378d-2*R23*R32*R38*R44        &
                   *R51*R53-7.8785817374287d-2*A*R23*R38*R44*R55                   &
                   -1.55614792451289d0*A*R32*R38*R44*R55-3.87898241533463d-3       & 
                   *A*R30*R32*R38*R44*R55-1.24070827436718d0*R12*R38*R44           &
                   *R56-6.25319771408057d-2*R12*R23*R32*R38*R44*R56                &
                   -1.57584066983129d-3*R32*R38*R4*R44*R50*R57                     &
                   -1.25073821694718d-3*R38*R44*R51*R53*R57-9.92705745752647d-4    &
                   *A*R32*R38*R44*R55*R57-7.27871572454413d-8*R38*R44*R46*R6       &
                   +4.00156883434654d-8*R30*R38*R44*R46*R6+1.35459439975431d-7     &
                   *R32*R38*R44*R46*R6-5.41837759901725d-7*R23*R32*R38*R44*R46     & 
                   *R6+1.02136417741475d-3*R23*R30*R32*R38*R44*R46*R6              &  
                   -9.84745021842697d-1*R*R1*R32*R38*R44*R50*R6                    &
                   +1.07513596816659d-7*R38*R44*R45*R8-4.30054387266636d-7*R23     &
                   *R38*R44*R45*R8+8.10652519997608d-4*R23*R30*R38*R44*R45*R8      &
                   -2.13897381686571d-5*R32*R38*R44*R45*R8+3.21686384535586d-2     &
                   *R30*R32*R38*R44*R45*R8)
               dEdAbsGradRho2(I)=dEdAbsGradRho2(I)+Scale*(  &
                   -3.29237170353606d-8*R30*R31*R32*R44*R61+2.2926652650931d-8*R14  &
                   *R32*R44*R45*R61-3.33446243530055d-8*R14*R30*R32*R44*R45         &
                   *R61-2.61314179335714d-8*R2*R30*R44*R46*R61                      &
                   +1.0584437918241d-4*A*R44*R50*R61-1.0584437918241d-4*A*R23       &
                   *R44*R50*R61-5.81893269438842d-5*A*R30*R44*R50*R61               &
                   +5.81893269438842d-5*A*R23                                       &        
                   *R30*R44*R50*R61-8.39627848290337d-3*A*R32*R44*R50*R61           &
                   -2.98108793700238d-3*A*R30*R32*R44*R50*R61+5.0793282454771d-6    &
                   *R44*R49*R51*R61-9.57453374272434d-3*R30*R44*R49*R51*R61         &
                   +1.33356345179728d-4*R32*R44*R49*R51*R61-1.33260345445857d-4     &
                   *R23*R32*R44*R49*R51*R61-7.3314388819187d-5*R30*R32*R44*R49      &    
                   *R51*R61-1.67964942310963d-4*R23*R30*R32*R44*R49*R51*R61         &
                   -6.66729121679367d-3*R12*R44*R53*R61+3.6654302424079d-3*R12      &
                   *R30*R44*R53*R61-8.40082456474881d-5*R12*R23*R32*R44*R53         &
                   *R61+4.61846279389034d-5*R12*R23*R30*R32*R44*R53*R61             &
                   -5.29180144160952d-3*R32*R44*R55*R61+2.90923681150097d-3         &
                   *R30*R32*R44*R55*R61+1.81967893113603d-8*R44*R47*R6*R61          &
                   -1.00039220858663d-8*R30*R44*R47*R6*R61-5.07972899907867d-8      &
                   *R32*R44*R47*R6*R61+1.5239186997236d-7*R23*R32*R44*R47*R6        &
                   *R61-3.83011566530532d-4*R23*R30*R32*R44*R47*R6*R61              &
                   -4.03175988062471d-8*R4*R44*R61*R8+1.20952796418741d-7           &
                   *R23*R4*R44*R61*R8-3.03994694999103d-4*R23*R30*R4*R44*R61        &
                   *R8+6.41403290732583d-6*R32*R4*R44*R61*R8-1.20632354500487d-2    &
                   *R30*R32*R4*R44*R61*R8)
            ENDIF
         ENDDO
      END SUBROUTINE PW91x_ClSh
!====================================================================================
!     Closed shell Becke 88 exchange engery functional
!     Physical Review A 38, p.3098 (1988), Equation (8).
!====================================================================================
      SUBROUTINE B88x_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,Scale) 
         INTEGER                        :: NGrid
         REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2
         REAL(DOUBLE), DIMENSION(NGrid) :: E,dEdRho,dEdAbsGradRho2
         REAL(DOUBLE)                   :: R,A,ASinh,DASinh,X,Scale
         INTEGER                        :: I
!-------------------------------------------------------------------------------------
         DO I=1,NGrid
            R=Rho(I)
            IF(NoNAN<R)THEN
               R=Two*R
               A=Four*AbsGradRho2(I)
               R1=R**(1.d0/3.d0)
               R2=R*R1
               R3=sqrt(A)
               R4=R**2
               R5=R1**2
               R6=R4*R5
               R7=1.D0/R6
               R8=A*R7
               R9=1.5874010519682d0*R8
               R10=1.d0+R9
               R11=sqrt(R10)
               R12=1.D0/R2
               R13=R12*R3
               R14=1.25992104989487d0*R13
               R15=R11+R14
        R16=log(R15)
        R17=A**2
        R18=1.D0/R11
        R19=1.D0/R15
        R20=R4**2
        R21=R12*R16*R3
        R22=3.17500104573508d-2*R21
        R23=1.d0+R22
        R24=R23**2
        R25=1.D0/R24
        R26=R*R20
        R27=A*R3
        R28=1.D0/R23
        RTmp1=-7.38558766382022d-1*R2
        RTmp2=(-5.29166840955847d-3*A)/(R2+3.17500104573508d-2*R16*R3)
        E(I)=E(I)+Scale*(RTmp1+RTmp2)
        RTmp1=-1.96949004368539d0*R1-(5.6448d-4*R17*R19*R25)/R26
        RTmp2=(1.41111157588226d-2*A*R28)/(R1*R4)
        RTmp3=(-7.11200234244658d-4*R17*R18*R19*R25*R3)/(R1*R20*R4)
        RTmp4=(-4.48028072907505d-4*R16*R25*R27)/(R*R4*R5)
        dEdRho(I)=dEdRho(I)+Scale*Half*(RTmp1+RTmp2+RTmp3+RTmp4)
        RTmp1=(4.2336d-4*A*R19*R25)/R20
        RTmp2=(5.33400175683494d-4*R18*R19*R25*R27)/(R1*R26)
        RTmp3=-2.11666736382339d-2*R12*R28+3.36021054680628d-4*R16*R25*R3*R7
        dEdAbsGradRho2(I)=dEdAbsGradRho2(I)+Scale*Half*(RTmp1+RTmp2+RTmp3)
            ENDIF
         ENDDO
      END SUBROUTINE B88x_ClSh
!====================================================================================
!     Closed shell Lee Yang Parr correlation engery functional 
!     Physical Review B 37, p.785 (1988)
!     Rederived by Miehlich, Savin, Stoll, and Preuss
!     Chemical Physics Letters 157, P.200 (1989), Equation (2).
!====================================================================================
      SUBROUTINE LYPc_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,Scale) 
         INTEGER                        :: NGrid
         REAL(DOUBLE), DIMENSION(NGrid) :: Rho,AbsGradRho2
         REAL(DOUBLE), DIMENSION(NGrid) :: E,dEdRho,dEdAbsGradRho2
         REAL(DOUBLE)                   :: R,A,ASinh,DASinh,X,Scale
         INTEGER                        :: I
!         ASinh(X)=LOG( X+SQRT(1.0D0+X))
!         DASinh(X)=1.0D0/SQRT(1.0D0+X)
!-------------------------------------------------------------------------------------
         DO I=1,NGrid
            R=Rho(I)
            IF(NoNAN<R)THEN
               R=Two*R
               A=Four*AbsGradRho2(I)


        R1=R**(1.d0/3.d0)
        R2=1.d0*R1
        R3=3.49d-1+R2
        R4=R3**2
        R5=1.0D0/R4
        R6=R1**2
        R7=R*R6
        R8=1.0D0/R7
        R9=1.0D0/R1
        R10=-2.533d-1*R9
        R11=exp(R10)
        R12=R**2
        R13=R*R12
        R14=R3*R4
        R15=1.0D0/R14
        R16=R1*R12
        R17=1.0D0/R16
        R18=1.0D0/R12
        RTmp1=R11*R5
        RTmp2=2.7049d-4*A*(1.26730158948386d-1+R1)*(1.62763650771828d0+R1)
        RTmp3=3.79002888024841d-1*R1+R3*exp(2.533d-1*R9)
        RTmp4=RTmp2-4.918d-2*R13*(1.3227200792067d-1+RTmp3)
        E(I)=E(I)+Scale*(R8*RTmp1*RTmp4)
        RTmp1=3.27866666666667d-2*R*R15+1.24262413553745d-2*R*R11*R15
        RTmp2=-3.71960957402222d-5*A*R11*R15*R17-3.16359093111111d-4*A*R11*R15*R18
        RTmp3=-2.28850933333333d-2*R1*R5-5.49250430212703d-4*R11*R5
        RTmp4=-1.02472999337096d-2*R1*R11*R5+(4.71088552549915d-6*A*R11*R5)/R13
        RTmp5=-6.09879813888889d-4*A*R11*R17*R5-2.7049d-4*A*R11*R18*R5
        RTmp6=(-5.29233602080333d-5*A*R11*R5)/(R12*R6)
        RTmp7=1.14425466666667d-2*R15*R6
        RTmp8=4.33675823302569d-3*R11*R15*R6-8.19666666666667d-2*R5*R6
        RTmp9=-3.10656033884362d-2*R11*R5*R6-1.80326666666667d-4*A*R11*R15*R8
        RTmp10=RTmp1+RTmp2+RTmp3+RTmp4
        RTmp11=RTmp5+RTmp6+RTmp7+RTmp8+RTmp9
        dEdRho(I)=dEdRho(I)+Scale*(RTmp10+RTmp11)
        RTmp1=(1.08196d-3*R11*R5)/R
        RTmp2=(1.89815455866667d-3*R11*R5)/(R*R1)
        RTmp3=2.23176574441333d-4*R11*R5*R8
        dEdAbsGradRho2(I)=dEdAbsGradRho2(I)+Scale*0.5D0*(RTmp1+RTmp2+RTmp3)
            ENDIF
         ENDDO
      END SUBROUTINE LYPc_ClSh
!====================================================================================
!     Closed shell Perdew Burke Ernnzerhof correlation GGA
!     Physical Review Letters 77, p.3865 (1996)
!====================================================================================
      SUBROUTINE PBEc_ClSh(NGrid,Rho,AbsGradRho2,E,dEdRho,dEdAbsGradRho2,Scale) 
         INTEGER                         :: NGrid
         REAL(DOUBLE), DIMENSION(NGrid)  :: Rho,AbsGradRho2
         REAL(DOUBLE), DIMENSION(NGrid)  :: E,dEdRho,dEdAbsGradRho2
         REAL(DOUBLE)                    :: R,A,Scale
         INTEGER                         :: I
!------------------------------------------------------------------------------------
         DO I=1,NGrid
            R=Rho(I)
            IF(NoNAN<R)THEN
               R=Two*R
               A=Four*AbsGradRho2(I)
               A=AbsGradRho2(I)
               R1=R**(1.d0/3.d0)
               R2=R1**2
               R3=1.D0/R2
               R4=1.89700432574756d-1*R3
               R5=sqrt(R)
               R6=1.D0/R5
               R7=8.00428634999363d-1*R6
               R8=1.D0/R1
               R9=2.22556942115069d0*R8
               R10=R**(1.d0/6.d0)
               R11=1.D0/R10
               R12=5.98255043577108d0*R11
               R13=R12+R4+R7+R9
               R14=1.D0/R13
               R15=1.60819794986925d1*R14
               R16=1.d0+R15
               R17=LOG(R16)
               R18=1.32568899905202d-1*R8
               R19=1.d0+R18
               R20=2.00000058733626d0*R19
               R21=R16**R20
               R22=-1.d0+R21
               R23=1.D0/R22
               R24=R**2
               R25=R1*R24
               R26=1.D0/R25
               R27=A*R23*R26
               R28=1.36210788856759d-1*R27
               R29=1.d0+R28
               R30=A**2
               R31=R22**2
               R32=1.D0/R31
               R33=R24**2
               R34=R2*R33
               R35=1.D0/R34
               R36=R30*R32*R35
               R37=1.85533790009806d-2*R36
               R38=R28+R37
               R39=1.d0+R38
               R40=1.D0/R39
               R41=A*R26*R29*R40
               R42=1.36210788856759d-1*R41
               R43=1.d0+R42
               R44=LOG(R43)
               R45=R30**2
               R46=2.65137877672926d-1*R8
               R47=2.00000058733626d0+R46
               R48=R16**R47
               R49=-1.d0+R48
               R50=R49**2
               R51=R50**2
               R52=1.D0/R51
               R53=1.00000058733626d0+R46
               R54=R16**R53
               R55=1.D0/R50
               R56=R30*R35*R55
               R57=1.85533790009806d-2*R56
               R58=1.D0/R49
               R59=A*R26*R58
               R60=1.36210788856759d-1*R59
               R61=R57+R60
               R62=1.d0+R61
               R63=R62**2
               R64=1.D0/R63
               R65=1.d0+R60
               R66=1.D0/R62
               R67=A*R26*R65*R66
               R68=1.36210788856759d-1*R67
               R69=1.d0+R68
               R70=1.D0/R69
               R71=R13**2
               R72=1.D0/R71
               R73=R33**2
               R74=R24*R73
               R75=R*R73
               R76=R10**2
               R77=R76**2
               R78=R10*R77
               R79=R2*R75
               R80=1.D0/R79
               R81=R49*R50
               R82=1.D0/R81
               R83=A*R30
               R84=R*R24*R33
               R85=R1*R84
               R86=1.D0/R85
               R87=R*R33
               R88=R2*R87
               R89=1.D0/R88
               R90=R5*R87
               R91=1.D0/R90
               R92=R1*R87
               R93=1.D0/R92
               R94=R10*R87
               R95=1.D0/R94
               R96=1.D0/R87
               R97=R33*R78
               R98=1.D0/R97
               R99=1.D0/R16
               R100=R*R24
               R101=R100*R2
               R102=1.D0/R101
               E(I)=E(I)+Scale*(  &
                   -6.21814d-2*R*R17-8.24331979256531d-3*R17*R2+3.10906908696549d-2*R*R44)
               RTmp1=-6.21814d-2*R17+3.10906908696549d-2*R44+1.34595386591642d-3            &
                    *R30*R35*R58*R64*R70-9.88140423540045d-3*A*R26*R66*R70                  &
                    -2.69190773183283d-3*R30*R35*R58*R66*R70+(1.54774525355071d-4           & 
                    *R45*R52*R54*R64*R70*R72)/R74+(1.15423633791477d-5*R45*R52              & 
                    *R54*R64*R70*R72)/(R1*R74)+(3.65266885748814d-5*R45*R52                 &
                    *R54*R64*R70*R72)/(R10*R74)+(6.86453535918717d-4*R45*R52*R54            & 
                    *R64*R70*R72)/(R5*R75)+(3.66532236089851d-4*R45*R52*R54*R64             & 
                    *R70*R72)/(R75*R78)-5.49554652837688d-3*R17*R8-1.89172028888574d-6      &
                    *R17*R45*R48*R52*R64*R70*R80+5.1073535107926d-4*R45*R52*R54*R64*R70     &
                    *R72*R80+(4.99439843843647d-5*R45*R64*R70*R82)/(R1*R75)                 &
                    +(1.27108470731556d-4*R54*R64*R70*R72*R82*R83)/R73+(5.5000031352384d-4  &
                    *R55*R64*R70*R83)/R84+(7.55946215802992d-3*R54*R64*R70*R72*R82*R83)     &
                    /(R10*R84)+(1.70443024360391d-3*R54*R64*R70*R72*R82*R83)/(R2*R84)       &
                    +(4.03637890030099d-3*R54*R64*R70*R72*R82*R83)/(R5*R84)                 &
                    +(4.02244442765395d-4*R54*R64*R70*R72*R82*R83)/(R78*R84)                &
                    -2.08322736924507d-5*R17*R48*R64*R70*R82*R83*R86
               RTmp2=5.62439314131374d-3*R54*R64*R70*R72*R82*R83*R86+3.11058254632152d-4    &
                     *R30*R54*R55*R64*R70*R72*R89-3.11058254632152d-4*R30*R54*R55*R66*R70   &
                     *R72*R89+9.84367553019607d-4*R30*R54*R55*R64*R70*R72*R91               &
                     -9.84367553019607d-4*R30*R54*R55*R66*R70*R72*R91+4.17106030515763d-3   &
                     *R30*R54*R55*R64*R70*R72*R93-4.17106030515763d-3*R30*R54*R55*R66*R70   &
                     *R72*R93+9.87777579680979d-3*R30*R54*R55*R64*R70*R72*R95               &
                     -9.87777579680979d-3*R30*R54*R55*R66*R70*R72*R95-5.09804787792498d-5   &
                     *R17*R30*R48*R55*R64*R70*R96+5.09804787792498d-5*R17*R30*R48*R55*R66   &
                     *R70*R96+1.37639443212006d-2*R30*R54*R55*R64*R70*R72*R96               &
                     -1.37639443212006d-2*R30*R54*R55*R66*R70*R72*R96+1.84994209378905d-2   &
                     *R30*R54*R55*R64*R70*R72*R98-1.84994209378905d-2*R30*R54*R55*R66*R70   &
                     *R72*R98-(1.67655851053175d-2*R72*R99)/R-9.9709173929518d-1*R11*R72    &
                     *R99-2.24814051658038d-1*R3*R72*R99-5.32397672482608d-1*R6*R72*R99     &
                     -(5.3055971797244d-2*R72*R99)/R78-7.41856473716896d-1*R72*R8*R99
               dEdRho(I)=dEdRho(I)+Scale*(RTmp1+RTmp2)
               dEdAbsGradRho2(I)=dEdAbsGradRho2(I) +Scale*(                                 &
                     (-2.35714420081646d-4*R30*R55*R64*R70)/(R24*R33)-5.76837371107036d-4   &
                    *A*R102*R58*R64*R70+(4.23488752945733d-3*R66*R70)/(R*R1)                &
                    +1.15367474221407d-3*A*R102*R58*R66*R70-(2.14045647361563d-5*R64*R70    &
                    *R82*R83)/(R1*R73))
            ENDIF
         ENDDO
      END SUBROUTINE PBEc_ClSh
END MODULE
