MODULE ParseOptions
  USE InOut
  USE Parse
  USE MemMan
  USE Numerics
  USE OptionKeys
  USE Functionals
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalObjects
  USE ProcessControl
  USE GlobalCharacters
  USE ParsingConstants
  USE ControlStructures
#ifdef NAG
  USE F90_UNIX
#endif
  IMPLICIT NONE
CONTAINS 
  !============================================================================
  !  LOAD THE OPTIONS OBJECT
  !============================================================================
  SUBROUTINE LoadOptions(N,O)
    TYPE(FileNames) :: N
    TYPE(Options)   :: O
    !-------------------------------------------------------------------------!  
    CALL OpenASCII(N%IFile,Inp)
    ! Check for restart, and extract the current state from the restart HDF file
    CALL ParseRestart(N%M_PWD,N%NewFileID,O%Guess,N%RFile,O%RestartState)
    ! Parse print output options, load global object PrintFlags
    CALL ParsePrintFlags(O%PFlags,O%GeomPrint)
    ! Parse for accuracy levels and set thresholds
    CALL ParseThresholds(O%NThrsh,O%AccuracyLevels,O%Thresholds)
    ! Parse for SCF methods to use in solution of SCF equations and put to HDF
    CALL ParseSCFMethods(O%NMthds,O%Methods)
    ! Parse for model chemistries 
    CALL ParseModelChems(O%NModls,O%Models)
    ! Parse for gradient options.  
    CALL ParseGradients(O%NSteps,O%Coordinates,O%Grad,O%DoGDIIS,O%SteepStep)
    ! Parse for NEB options.
    CALL ParseNEB(O%NEBSpring,O%NEBClimb)
    CLOSE(UNIT=Inp,STATUS='KEEP')
  END SUBROUTINE LoadOptions
  !============================================================================
  !  PARSE THE METHODS TO USE IN SOLUTION OF THE SCF EQUATIONS
  !============================================================================
  SUBROUTINE ParseModelChems(NModls,Models)
    INTEGER                       :: NModls,I
    INTEGER,   DIMENSION(MaxSets) :: Models
    !-------------------------------------------------------------------------!
    NModls=0
    IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_ExactX,MaxSets,NLoc,Location))THEN
       NModls=NModls+NLoc
       DO I=1,NLoc
          Models(Location(I))=EXACT_EXCHANGE
       ENDDO
    ENDIF
    ! Slater-Dirac exchage
    IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_SD,MaxSets,NLoc,Location))THEN
       NModls=NModls+NLoc
       DO I=1,NLoc
          Models(Location(I))=SD_EXCHANGE
       ENDDO
    ENDIF
    ! X-alpha exchage
    IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_XA,MaxSets,NLoc,Location))THEN
       NModls=NModls+NLoc
       DO I=1,NLoc 
          Models(Location(I))=XA_EXCHANGE
       ENDDO
    ENDIF
    ! PW91 exchange
    IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_PW91x,MaxSets,NLoc,Location))THEN
       NModls=NModls+NLoc
       DO I=1,NLoc
          Models(Location(I))=PW91_EXCHANGE
       ENDDO
    ENDIF
    ! PBE exchange
    IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_PBEx,MaxSets,NLoc,Location))THEN
       NModls=NModls+NLoc
       DO I=1,NLoc
          Models(Location(I))=PBE_EXCHANGE
       ENDDO
    ENDIF
    ! B88 exchange
    IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_B88x,MaxSets,NLoc,Location))THEN
       NModls=NModls+NLoc
       DO I=1,NLoc 
          Models(Location(I))=B88_EXCHANGE
       ENDDO
    ENDIF
    ! Pure Slater exchange with VWN3 LSDA correlation 
    IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_VWN3xc,MaxSets,NLoc,Location))THEN
       NModls=NModls+NLoc
       DO I=1,NLoc
          Models(Location(I))=PURE_VWN3_LSD
       ENDDO
    ENDIF
    ! Pure Slater exchange with VWN5 LSDA correlation 
    IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_VWN5xc,MaxSets,NLoc,Location))THEN
       NModls=NModls+NLoc
       DO I=1,NLoc
          Models(Location(I))=PURE_VWN5_LSD
       ENDDO
    ENDIF
    ! Pure Slater exchange with PW91 LSDA correlation 
    IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_PW91xc,MaxSets,NLoc,Location))THEN
       NModls=NModls+NLoc
       DO I=1,NLoc
          Models(Location(I))=PURE_PW91_LSD
       ENDDO
    ENDIF
    ! Pure PBE GGA exchange-correlation 
    IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_PBExc,MaxSets,NLoc,Location))THEN
       NModls=NModls+NLoc
       DO I=1,NLoc
          Models(Location(I))=PURE_PBE_GGA
       ENDDO
    ENDIF
    ! Pure BLYP GGA exchange-correlation 
    IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_BLYPxc,MaxSets,NLoc,Location))THEN
       NModls=NModls+NLoc
       DO I=1,NLoc
          Models(Location(I))=PURE_BLYP_GGA
       ENDDO
    ENDIF
    ! Hybrid B3LYP/VWN3 exchange-correlation 
    IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_B3LYP_VWN3xc,MaxSets,NLoc,Location))THEN
       NModls=NModls+NLoc
       DO I=1,NLoc
          Models(Location(I))=HYBRID_B3LYP_VWN3
       ENDDO
    ENDIF
    ! Hybrid B3LYP/VWN5 exchange-correlation 
    IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_B3LYP_VWN5xc,MaxSets,NLoc,Location))THEN
       NModls=NModls+NLoc
       DO I=1,NLoc
          Models(Location(I))=HYBRID_B3LYP_VWN5
       ENDDO
    ENDIF
    ! Hybrid B3LYP/PW91 exchange-correlation 
    IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_B3LYP_PW91xc,MaxSets,NLoc,Location))THEN
       NModls=NModls+NLoc
       DO I=1,NLoc
          Models(Location(I))=HYBRID_B3LYP_PW91
       ENDDO
    ENDIF
    ! Hybrid PBE0 exchange-correlation 
    IF(OptKeyLocQ(Inp,MODEL_OPTION,MODEL_PBE0xc,MaxSets,NLoc,Location))THEN
       NModls=NModls+NLoc
       DO I=1,NLoc 
          Models(Location(I))=HYBRID_PBE0
       ENDDO
    ENDIF
    IF(NModls==0) &
         CALL MondoHalt(PRSE_ERROR,'Option '//MODEL_OPTION//' not set in input.'//RTRN   &
         //'Options include '//MODEL_ExactX//', '//MODEL_SD//', '//       &
         MODEL_XA//', '//MODEL_PW91x//', '//MODEL_PBEx//', '//           &
         MODEL_B88x//', '//MODEL_VWN3xc//', '//MODEL_VWN5xc//', '//      &
         MODEL_PW91xc//', '//MODEL_PBExc//', '//MODEL_BLYPxc//', '//     &
         MODEL_B3LYP_VWN3xc//', '//MODEL_B3LYP_PW91xc//', and '//MODEL_PBE0xc)
  END SUBROUTINE ParseModelChems
  !============================================================================
  !  PARSE THE METHODS TO USE IN SOLUTION OF THE SCF EQUATIONS
  !============================================================================
  SUBROUTINE ParseSCFMethods(NMthds,Methods)
    INTEGER                       :: NMthds,NLoc,I
    INTEGER,   DIMENSION(MaxSets) :: Methods
    !----------------------------------------------------------------------------
    NMthds=0
    IF(OptKeyLocQ(Inp,SCF_OPTION,SCF_SDMM,MaxSets,NLoc,Location))THEN
       NMthds=NMthds+NLoc
       DO I=1,NLoc 
          Methods(Location(I))=SDMM_R_SCF
       ENDDO
    ENDIF
    IF(OptKeyLocQ(Inp,SCF_OPTION,SCF_PM,MaxSets,NLoc,Location))THEN
       NMthds=NMthds+NLoc
       DO I=1,NLoc 
          Methods(Location(I))=PM_R_SCF
       ENDDO
    ENDIF
    IF(OptKeyLocQ(Inp,SCF_OPTION,SCF_SP2,MaxSets,NLoc,Location))THEN
       NMthds=NMthds+NLoc
       DO I=1,NLoc 
          Methods(Location(I))=SP2_R_SCF 
       ENDDO
    ENDIF
    IF(OptKeyLocQ(Inp,SCF_OPTION,SCF_SP4,MaxSets,NLoc,Location))THEN
       NMthds=NMthds+NLoc
       DO I=1,NLoc 
          Methods(Location(I))=SP4_R_SCF
       ENDDO
    ENDIF
    IF(OptKeyLocQ(Inp,SCF_OPTION,SCF_TS4,MaxSets,NLoc,Location))THEN
       NMthds=NMthds+NLoc
       DO I=1,NLoc 
          Methods(Location(I))=TS4_R_SCF 
       ENDDO
    ENDIF
    IF(OptKeyLocQ(Inp,SCF_OPTION,SCF_RHHF,MaxSets,NLoc,Location))THEN
       NMthds=NMthds+NLoc
       DO I=1,NLoc 
          Methods(Location(I))=RH_R_SCF 
       ENDDO
    ENDIF
    IF(NMthds==0) &
         CALL MondoHalt(PRSE_ERROR,'Option '//SCF_OPTION//' not set in input.'//RTRN   &
         //'Options include '//SCF_SDMM//', '//SCF_PM//', '//SCF_SP2//', '//SCF_TS4//', and '//SCF_RHHF)
  END SUBROUTINE ParseSCFMethods
  !===============================================================================================
  ! PARSE FOR ACCURACY LEVELS AND SET THRESHOLDS
  !===============================================================================================
  SUBROUTINE ParseThresholds(NThrsh,AccuracyLevels,Thresholds)
    INTEGER                       :: NThrsh,I
    INTEGER,   DIMENSION(MaxSets) :: AccuracyLevels
    TYPE(TOLS),DIMENSION(MaxSets) :: Thresholds
    !-----------------------------------------------------------------------------------------------
    NThrsh=0
    IF(OptKeyLocQ(Inp,ACCURACY_OPTION,ACCURACY_CHEEZY,MaxSets,NLoc,Location))THEN
       NThrsh=NThrsh+NLoc
       DO I=1,NLoc 
          AccuracyLevels(Location(I))=1
       ENDDO
    ENDIF
    IF(OptKeyLocQ(Inp,ACCURACY_OPTION,ACCURACY_GOOD,MaxSets,NLoc,Location))THEN
       NThrsh=NThrsh+NLoc
       DO I=1,NLoc 
          AccuracyLevels(Location(I))=2
       ENDDO
    ENDIF
    IF(OptKeyLocQ(Inp,ACCURACY_OPTION,ACCURACY_TIGHT,MaxSets,NLoc,Location))THEN
       NThrsh=NThrsh+NLoc
       DO I=1,NLoc 
          AccuracyLevels(Location(I))=3
       ENDDO
    ENDIF
    IF(OptKeyLocQ(Inp,ACCURACY_OPTION,ACCURACY_RETENTIVE,MaxSets,NLoc,Location))THEN
       NThrsh=NThrsh+NLoc
       DO I=1,NLoc 
          AccuracyLevels(Location(I))=4 
       ENDDO
    ENDIF
    IF(NThrsh==0) &
         CALL MondoHalt(PRSE_ERROR,'Option '//ACCURACY_OPTION//' not set in input.'//RTRN   &
         //'Options include '//ACCURACY_CHEEZY//', '//ACCURACY_GOOD//', '   &
         //ACCURACY_TIGHT//', and '//ACCURACY_RETENTIVE)
    ! Set thresholds
    DO I=1,NThrsh
       Thresholds(I)%Cube=CubeNeglect(AccuracyLevels(I))
       Thresholds(I)%Trix=TrixNeglect(AccuracyLevels(I))
       Thresholds(I)%Dist=DistNeglect(AccuracyLevels(I))
       Thresholds(I)%TwoE=TwoENeglect(AccuracyLevels(I))
       Thresholds(I)%ETol=ETol(AccuracyLevels(I))
       Thresholds(I)%DTol=DTol(AccuracyLevels(I))
    ENDDO
  END SUBROUTINE ParseThresholds
  !===============================================================================================
  ! PARSE THE INPUT FILE FOR GUESS OPTIONS, GET PREVIOUS CURRENT SCF STATE IF RESTART
  !===============================================================================================
  SUBROUTINE ParseRestart(M_PWD,NewFileID,Guess,RestartHDF,RestartState)
    CHARACTER(LEN=*)                 :: M_PWD
    CHARACTER(LEN=DCL)               :: RestartHDF
    TYPE(INT_VECT)                   :: RestartState
    INTEGER                          :: Guess,NewFileID
    !-----------------------------------------------------------------------------------------------         
    IF(OptKeyQ(Inp,GUESS_OPTION,GUESS_RESTART))THEN
       Guess=GUESS_EQ_RESTART
       IF(.NOT.OptCharQ(Inp,RESTART_INFO,RestartHDF))  &
            CALL MondoHalt(PRSE_ERROR,'Restart requested, but no HDF file specified.')
       ! Check for absolute path 
       IF(RestartHDF(1:1)/='/')  &
            CALL MondoHalt(PRSE_ERROR,'Please use absolute path to the restart HDF file')
       CALL New(RestartState,3)
       ! Open the old restart HDF file
       HDF_CurrentID=OpenHDF(RestartHDF)       
       ! Get the current state to restart from
       CALL Get(RestartState,'current_state')
       ! Now close the old file...
       CALL CloseHDF(HDF_CurrentID)
    ELSEIF(OptKeyQ(Inp,GUESS_OPTION,GUESS_CORE))THEN
       Guess=GUESS_EQ_CORE
    ELSEIF(OptKeyQ(Inp,GUESS_OPTION,GUESS_SUPER))THEN
       Guess=GUESS_EQ_SUPR
    ELSE
       CALL MondoHalt(PRSE_ERROR,'No guess specified in input file')
    ENDIF
  END SUBROUTINE ParseRestart
  !===============================================================================================
  ! PARSE THE PRINT OUT OPTIONS AND LOAD THE GLOBAL PRINT FLAGS OBJECT.
  !===============================================================================================
  SUBROUTINE ParsePrintFlags(PFlags,GeomPrint)
    TYPE(DEBG)         :: PFlags
    CHARACTER(LEN=3)   :: GeomPrint
    !-----------------------------------------------------------------------------------------------
    IF(OptKeyQ(Inp,GLOBAL_DEBUG,DBG_NONE) )THEN
       PFlags%Key=DEBUG_NONE
    ELSEIF(OptKeyQ(Inp,GLOBAL_DEBUG,DBG_MEDIUM) )THEN
       PFlags%Key=DEBUG_MEDIUM
    ELSEIF(OptKeyQ(Inp,GLOBAL_DEBUG,DBG_MAXIMUM) )THEN
       PFlags%Key=DEBUG_MAXIMUM
    ELSE
       PFlags%Key=DEBUG_MINIMUM
    ENDIF
    IF(OptKeyQ(Inp,GLOBAL_DEBUG,DBG_MATRICES) )THEN
       PFlags%Mat=DEBUG_MATRICES
    ELSEIF(OptKeyQ(Inp,GLOBAL_DEBUG,PLT_MATRICES) )THEN
       PFlags%Mat=PLOT_MATRICES
    ELSE
       PFlags%Mat=DEBUG_NONE
    ENDIF
    IF(OptKeyQ(Inp,GLOBAL_DEBUG,DBG_CHKSUMS))THEN
       PFlags%Chk=DEBUG_CHKSUMS
    ELSE
       PFlags%Chk=DEBUG_NONE
    ENDIF
    IF(OptKeyQ(Inp,GLOBAL_DEBUG,DBG_MMA_STYLE) )THEN
       PFlags%Fmt=DEBUG_MMASTYLE
    ELSEIF(OptKeyQ(Inp,GLOBAL_DEBUG,DBG_FLT_STYLE) )THEN
       PFlags%Fmt=DEBUG_FLTSTYLE
    ELSE
       PFlags%Fmt=DEBUG_DBLSTYLE
    ENDIF
    IF(OptKeyQ(Inp,GLOBAL_DEBUG,DBG_PRT_INTS))THEN
       PFlags%Int=DEBUG_INTEGRAL
    ELSE
       PFlags%Int=DEBUG_NONE
    ENDIF
    IF(OptKeyQ(Inp,GLOBAL_DEBUG,DBG_PRT_SETS))THEN
       PFlags%Set=DEBUG_BASISSET
    ELSE
       PFlags%Set=DEBUG_NONE
    ENDIF
    IF(OptKeyQ(Inp,GLOBAL_DEBUG,DBG_PRT_GEOP))THEN
       PFlags%GeOp=DEBUG_GEOP    
    ELSE
       PFlags%GeOp=DEBUG_NONE
    ENDIF
    IF(OptKeyQ(Inp,GLOBAL_DEBUG,DBG_PRT_MM))THEN
       PFlags%MM=DEBUG_MM    
    ELSE
       PFlags%MM=DEBUG_NONE
    ENDIF
!    GeomPrint='PDB' ! hard set for now
    GeomPrint='XYZ' ! hard set for now
  END SUBROUTINE ParsePrintFlags
  !===============================================================================================
  !
  !===============================================================================================
  SUBROUTINE ParseGradients(NSteps,Coordinates,Grad,DoGDIIS,SteepStep)
    INTEGER NSteps,Coordinates,Grad
    LOGICAL OneBase,DoGDIIS,SteepStep
    !-----------------------------------------------------------------------------------------------
    ! Default max geometry steps is 100
    IF(.NOT.OptIntQ(Inp,GRADIENTS,NSteps))THEN
       NSteps=100               
    ENDIF
    ! Macro gradient options
    IF(OptKeyQ(Inp,GRADIENTS,GRAD_FORCE))THEN
       Grad=GRAD_ONE_FORCE
       NSteps=1
    ELSEIF(OptKeyQ(Inp,GRADIENTS,GRAD_DYNAMICS))THEN
       ! Do molecular dynamics
       Grad=GRAD_DO_DYNAMICS    
    ELSEIF(OptKeyQ(Inp,GRADIENTS,GRAD_OPTIMIZE))THEN
       ! Go downhill in energy
       Grad=GRAD_GO_DOWNHILL
    ELSEIF(OptKeyQ(Inp,GRADIENTS,GRAD_TS_SEARCH))THEN
       ! Transition state search 
       Grad=GRAD_TS_SEARCH_NEB 
    ELSE
       ! Single point energy only
       Grad=GRAD_NO_GRAD       
       NSTeps=1
    ENDIF
    ! Use internal or Cartesian coordinates ?
    IF(OptKeyQ(Inp,GRADIENTS,GRAD_INTERNALS)) THEN
       ! Yes, use internals where available
       Coordinates=GRAD_INTS_OPT 
    ELSE
       ! Default is Cartesians
       Coordinates=GRAD_CART_OPT  
    ENDIF 
    ! Use GDIIS in optimization ?
    IF(OptKeyQ(Inp,GRADIENTS,GRAD_GDIIS))THEN
       ! Yep
       DoGDIIS=.TRUE.             
    ELSE
       ! Nope
       DoGDIIS=.FALSE.    
    ENDIF
    ! Use approximate second order methods in optimization
    IF(OptKeyQ(Inp,GRADIENTS,GRAD_APPRX_HESS))THEN
       ! Yeah, but only if using internal coordinates
       SteepStep=.FALSE.
    ELSE
       ! No, gradients only
       SteepStep=.TRUE.
    ENDIF    
  END SUBROUTINE ParseGradients
  !===============================================================================================
  !
  !===============================================================================================
  SUBROUTINE ParseNEB(NEBSpring,NEBClimb)
    !-----------------------------------------------------------------------------------------------
    REAL(DOUBLE) :: NEBSpring
    Logical      :: NEBClimb
    ! Set the spring constant between NEB images
    IF(.NOT.OptDblQ(Inp,NEB_SPRING,NEBSpring))THEN
       NEBSpring=1.0
    ENDIF
    ! Use the climbing image?
    IF(OptKeyQ(Inp,NEB_CLIMB,NEB_DO_CLIMB))THEN
       NEBClimb=.TRUE.
    ELSE
       NEBClimb=.FALSE.
    ENDIF    
  END SUBROUTINE ParseNEB
END MODULE ParseOptions
