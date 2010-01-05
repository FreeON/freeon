!------------------------------------------------------------------------------
! This code is part of the MondoSCF suite of programs for linear scaling
! electronic structure theory and ab initio molecular dynamics.
!
! Copyright (2004). The Regents of the University of California. This
! material was produced under U.S. Government contract W-7405-ENG-36
! for Los Alamos National Laboratory, which is operated by the University
! of California for the U.S. Department of Energy. The U.S. Government has
! rights to use, reproduce, and distribute this software.  NEITHER THE
! GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
! OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by the
! Free Software Foundation; either version 2 of the License, or (at your
! option) any later version. Accordingly, this program is distributed in
! the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
! the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
! PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
! While you may do as you like with this software, the GNU license requires
! that you clearly mark derivative software.  In addition, you are encouraged
! to return derivative works to the MondoSCF group for review, and possible
! disemination in future releases.
!------------------------------------------------------------------------------
! Authors : Hugh Nymeyer, CJ Tymczak and Matt Challacombe

MODULE ParseDynamics
  USE Parse
  USE InOut
  USE PrettyPrint
  USE DynamicsKeys
  USE ControlStructures
  USE OptionKeys
  USE MondoLogger

  IMPLICIT NONE

CONTAINS
  !========================================================================================
  ! PARSE FOR MOLECULAR DYNAMICS DIRECTIVES AND VALUES
  !========================================================================================
  SUBROUTINE LoadDynamics(N,O,D)
    TYPE(FileNames)    :: N
    TYPE(Options)      :: O
    TYPE(Dynamics)     :: D
    CHARACTER(LEN=DCL) :: Line
    !---------------------------------------------------------------------------------------!
    CALL OpenASCII(N%IFile,Inp)

    ! Initialize
    D%DoingMD        = .TRUE.
    D%Initial_Temp   = .FALSE.
    D%Temp_Scaling   = .FALSE.
    D%Const_Temp     = .FALSE.
    D%Const_Press    = .FALSE.
    D%Parallel_Rep   = .FALSE.
    D%MDGuess        = MD_DGuess

    ! Parse MD
    IF(O%Grad==GRAD_DO_DYNAMICS) THEN
      ! Parse the Density Matrix Projection Algoithm
      IF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMLinear)) THEN
        D%MDGuess=MD_DMLinear
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMTRBO)) THEN
        D%MDGuess=MD_DMTRBO
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMTRBO_Damp_dt3)) THEN
        D%MDGuess=MD_DMTRBO_Damp_dt3
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMTRBO_Damp_dt5)) THEN
        D%MDGuess=MD_DMTRBO_Damp_dt5
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMTRBO_Damp_dt7)) THEN
        D%MDGuess=MD_DMTRBO_Damp_dt7
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMTRBO_Damp_dt9)) THEN
        D%MDGuess=MD_DMTRBO_Damp_dt9
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMTRBO_Damp_dt11)) THEN
        D%MDGuess=MD_DMTRBO_Damp_dt11
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMSymplectic)) THEN
        D%MDGuess=MD_DMSymplectic
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_FMVerlet0)) THEN
        D%MDGuess=MD_FMVerlet0
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_FMVerlet1)) THEN
        D%MDGuess=MD_FMVerlet1
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP0)) THEN
        D%MDGuess=MD_DMP0
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP1)) THEN
        D%MDGuess=MD_DMP1
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP2)) THEN
        D%MDGuess=MD_DMP2
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP3)) THEN
        D%MDGuess=MD_DMP3
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP4)) THEN
        D%MDGuess=MD_DMP4
      ELSE
        CALL MondoLog(DEBUG_NONE, "LoadDynamics", 'In input file, no MD DM Projection (MDProjection) algorithm defined')
        CALL MondoLog(DEBUG_NONE, "LoadDynamics", "using DMLinear")
        D%MDGuess = MD_DMLinear
      ENDIF

      ! Parse MD Options: First MD Algorithmn
      IF(OptKeyQ(Inp,MD_AL_OPTION,MD_AL_VERLET))THEN
        D%MDAlgorithm = MD_AL_VERLET
      ELSEIF(OptKeyQ(Inp,MD_AL_OPTION,MD_AL_GEAR)) THEN
        D%MDAlgorithm = MD_AL_GEAR
      ELSEIF(OptKeyQ(Inp, MD_AL_OPTION, MD_AL_SYMPLECTIC)) THEN
        D%MDAlgorithm = MD_AL_SYMPLECTIC
      ELSE
        CALL MondoLog(DEBUG_NONE, "LoadDynamics", "No MD algorithm specified, setting to Verlet")
        D%MDAlgorithm = MD_AL_VERLET
      ENDIF
      ! MD Time Step
      IF(.NOT. OptDblQ(Inp,MD_TIME_STEP,D%DTime)) THEN
        CALL MondoHalt(PRSE_ERROR, "option "//MD_TIME_STEP//' not found in input.')
      ELSE
        ! Convert that to internal time units.
        D%DTime = D%DTime*FemtosecondsToInternalTime
      ENDIF
      ! MD MaxSteps
      IF(.NOT. OptIntQ(Inp,MD_MAX_STEP,D%MDMaxSteps)) THEN
        CALL MondoHalt(PRSE_ERROR,MD_MAX_STEP//' not found in input.')
      ENDIF
      ! Parse for Number of SCF cycles
      IF(.NOT. OptIntQ(Inp,MD_NUM_SCF,D%MDNumSCF)) THEN
        CALL MondoLog(DEBUG_NONE, "LoadDynamics", 'Number of SCF cycles for MD is not set')
        D%MDNumSCF = -1
      ENDIF

      ! Parse for Initial Temp, if any
      IF(OptDblQ(Inp,MD_INITIAL_TEMP,D%TempInit))THEN
        D%Initial_Temp = .TRUE.
      ELSE
        D%Initial_Temp = .FALSE.
        D%TempInit     = Zero
      ENDIF

      ! Parse for Temperature Scaling
      IF(.NOT.OptDblQ(Inp,MD_TARGET_TEMP,D%TargetTemp)) THEN
        D%TargetTemp = Zero
        D%RescaleInt = 1
      ENDIF

      IF(.NOT.OptIntQ(Inp,MD_TSCALE_STEPS,D%RescaleInt)) THEN
        D%RescaleInt = D%MDMaxSteps
      ENDIF

      ! Parse for Energy Scaling
      IF(OptDblQ(Inp,MD_TARGET_ETOTAL,D%TargetEtotal)) THEN
        D%Energy_Scaling_Set = .TRUE.
        ! Change units of the energy to internal units, i.e. Hartree.
        D%TargetEtotal = eV2au*D%TargetEtotal
        CALL MondoLog(DEBUG_NONE, "LoadDynamics", "thermostating Etotal to "//TRIM(DblToChar(D%TargetEtotal*au2eV))//" eV")
      ELSE
        D%Energy_Scaling_Set = .FALSE.
        D%TargetEtotal       = Zero
        D%RescaleInt         = 1
      ENDIF

      ! Parse for thermostat.
      IF(OptKeyQ(Inp, MD_THERMOSTAT, MD_THERM_SCALING)) THEN
        CALL MondoLog(DEBUG_NONE, "LoadDynamics", "Using Scaling thermostat")
        D%Thermostat = MD_THERM_SCALING
      ELSEIF(OptKeyQ(Inp, MD_THERMOSTAT, MD_THERM_BERENDSEN)) THEN
        CALL MondoLog(DEBUG_NONE, "LoadDynamics", "Using Berendsen thermostat")
        D%Thermostat = MD_THERM_BERENDSEN
      ELSEIF(OptKeyQ(Inp, MD_THERMOSTAT, MD_THERM_BERENDSEN_ETOT)) THEN
        CALL MondoLog(DEBUG_NONE, "LoadDynamics", "Using Berendsen thermostat for Etotal")
        D%Thermostat = MD_THERM_BERENDSEN_ETOT
      ELSEIF(OptKeyQ(Inp, MD_THERMOSTAT, MD_ACT)) THEN
        CALL MondoLog(DEBUG_NONE, "LoadDynamics", "Using Action Control Theory thermostat")
        D%Thermostat = MD_ACT
      ELSE
        CALL MondoLog(DEBUG_NONE, "LoadDynamics", "Unknown thermostat or thermostat not set")
        D%Thermostat = MD_THERM_UNSET
      ENDIF

      IF(OptDblQ(Inp, MD_ACT_MAX_FORCE_ERROR, D%ACTMaxForceError)) THEN
        CALL MondoLog(DEBUG_NONE, "LoadDynamics", "ACT Max Force Error = "//TRIM(DblToChar(D%ACTMaxForceError)))
      ELSE
        IF(D%Thermostat == MD_ACT) THEN
          D%ACTMaxForceError = 1.D-6
          CALL MondoLog(DEBUG_NONE, "LoadDynamics", "ACT Max Force Error defaults to "//TRIM(DblToChar(D%ACTMaxForceError)))
        ENDIF
      ENDIF

      IF(OptDblQ(Inp, MD_ACT_ALPHA, D%ACTAlpha)) THEN
        CALL MondoLog(DEBUG_NONE, "LoadDynamics", "ACT alpha = "//TRIM(DblToChar(D%ACTAlpha)))
      ELSE
        IF(D%Thermostat == MD_ACT) THEN
          D%ACTAlpha = 1.D0
          CALL MondoLog(DEBUG_NONE, "LoadDynamics", "ACT alpha defaults to "//TRIM(DblToChar(D%ACTAlpha)))
        ENDIF
      ENDIF

      IF(OptDblQ(Inp, MD_ACT_BETA, D%ACTBeta)) THEN
        CALL MondoLog(DEBUG_NONE, "LoadDynamics", "ACT beta = "//TRIM(DblToChar(D%ACTBeta)))
      ELSE
        IF(D%Thermostat == MD_ACT) THEN
          D%ACTBeta = 1.D-2
          CALL MondoLog(DEBUG_NONE, "LoadDynamics", "ACT beta defaults to "//TRIM(DblToChar(D%ACTBeta)))
        ENDIF
      ENDIF

      IF(OptDblQ(Inp, MD_BERENDSEN_TAU, D%BerendsenTau)) THEN
        D%BerendsenTau = D%BerendsenTau*FemtosecondsToInternalTime
        CALL MondoLog(DEBUG_NONE, "LoadDynamics", "Using tau = "//TRIM(DblToChar(D%BerendsenTau*InternalTimeToFemtoseconds))//" fs")
      ELSE
        CALL MondoLog(DEBUG_NONE, "LoadDynamics", "BerendsenTau not specified in input: Setting tau to "//MD_TIME_STEP)
        D%BerendsenTau = D%DTime
      ENDIF

      ! Parse MD damping.
      IF(.NOT.OptDblQ(Inp, MD_DAMPING, D%MDalpha)) THEN
        D%MDalpha = 0.0d0
      ENDIF

      ! Parse how many time steps damping should be used.
      IF(.NOT.OptIntQ(Inp, MD_DAMP_STEP, D%MDDampStep)) THEN
        ! We will set this to a very large number here, which means that we want
        ! damping for all time steps. If we want to be smarter about this and
        ! use MDMaxSteps for example, we need to make sure that the parameters
        ! MDDampStep and MDMaxSteps are set in the correct order in the input
        ! file.
        D%MDDampStep = 10000000
      ENDIF

      ! Parse Hybrid MC
    ELSEIF(O%Grad==GRAD_DO_HYBRIDMC) THEN
      ! Parse the Density Matrix Projection Algoithm
      IF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMLinear)) THEN
        D%MDGuess=MD_DMLinear
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMTRBO)) THEN
        D%MDGuess=MD_DMTRBO
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMSymplectic)) THEN
        D%MDGuess=MD_DMSymplectic
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_FMVerlet0)) THEN
        D%MDGuess=MD_FMVerlet0
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_FMVerlet1)) THEN
        D%MDGuess=MD_FMVerlet1
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP0)) THEN
        D%MDGuess=MD_DMP0
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP1)) THEN
        D%MDGuess=MD_DMP1
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP2)) THEN
        D%MDGuess=MD_DMP2
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP3)) THEN
        D%MDGuess=MD_DMP3
      ELSEIF(OptKeyQ(Inp,MD_PM_OPTION,MD_DMP4)) THEN
        D%MDGuess=MD_DMP4
      ELSE
        CALL MondoLog(DEBUG_NONE, "LoadDynamics", 'In input file, no MD DM Projection algorithm defined')
      ENDIF
      ! MD Algorithmn
      IF(OptKeyQ(Inp,MD_AL_OPTION,MD_AL_VERLET))THEN
        D%MDAlgorithm = MD_AL_VERLET
      ELSEIF(OptKeyQ(Inp,MD_AL_OPTION,MD_AL_GEAR)) THEN
        D%MDAlgorithm = MD_AL_GEAR
      ELSE
        CALL MondoLog(DEBUG_NONE, "LoadDynamics", "No MD algorithm specified, setting to Verlet")
        D%MDAlgorithm = MD_AL_VERLET
      ENDIF
      ! MD Time Step
      IF(.NOT. OptDblQ(Inp,MD_TIME_STEP,D%DTime)) THEN
        CALL MondoHalt(PRSE_ERROR,MD_TIME_STEP//' not found in input.')
      ELSE
        ! Convert that to internal time units.
        D%DTime = D%DTime*FemtosecondsToInternalTime
      ENDIF
      ! MD MaxSteps
      IF(.NOT. OptIntQ(Inp,MD_MAX_STEP,D%MDMaxSteps)) THEN
        CALL MondoHalt(PRSE_ERROR,MD_MAX_STEP//' not found in input.')
      ENDIF
      ! Parse Hybrid MC Options: First  MC MaxSteps
      IF(.NOT. OptIntQ(Inp,MC_MAX_STEP,D%MCMaxSteps)) THEN
        CALL MondoHalt(PRSE_ERROR,MC_MAX_STEP//' not found in input.')
      ENDIF
      ! MC Temp
      IF(.NOT. OptDblQ(Inp,MC_TEMP,D%MCTemp))THEN
        CALL MondoHalt(PRSE_ERROR,MC_TEMP//' not found in input.')
      ENDIF
    ENDIF

  END SUBROUTINE LoadDynamics
  !---------------------------------------------------------------------------------------!
END MODULE ParseDynamics
