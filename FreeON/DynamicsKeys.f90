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
MODULE DynamicsKeys

  IMPLICIT NONE

  ! DM Projection Method
  CHARACTER(LEN=*), PARAMETER :: MD_PM_OPTION           = 'MDProjection'

  ! DM Verlet
  CHARACTER(LEN=*), PARAMETER :: MD_DMLinear            = 'DMLinear'
  CHARACTER(LEN=*), PARAMETER :: MD_DMTRBO              = 'DMTRBO'
  CHARACTER(LEN=*), PARAMETER :: MD_DMTRBO_Damp_dt3     = 'DMTRBO_Damp_dt3'
  CHARACTER(LEN=*), PARAMETER :: MD_DMTRBO_Damp_dt5     = 'DMTRBO_Damp_dt5'
  CHARACTER(LEN=*), PARAMETER :: MD_DMTRBO_Damp_dt7     = 'DMTRBO_Damp_dt7'
  CHARACTER(LEN=*), PARAMETER :: MD_DMTRBO_Damp_dt9     = 'DMTRBO_Damp_dt9'
  CHARACTER(LEN=*), PARAMETER :: MD_DMTRBO_Damp_dt11    = 'DMTRBO_Damp_dt11'
  CHARACTER(LEN=*), PARAMETER :: MD_DMSymplectic        = 'DMSymplectic'

  ! FM Verlet
  CHARACTER(LEN=*), PARAMETER :: MD_FMVerlet0           = 'FMVerlet0'
  CHARACTER(LEN=*), PARAMETER :: MD_FMVerlet1           = 'FMVerlet1'

  ! Stencil Denisty Matrix Projectors
  CHARACTER(LEN=*), PARAMETER :: MD_DMP0                = 'DMProj0'
  CHARACTER(LEN=*), PARAMETER :: MD_DMP1                = 'DMProj1'
  CHARACTER(LEN=*), PARAMETER :: MD_DMP2                = 'DMProj2'
  CHARACTER(LEN=*), PARAMETER :: MD_DMP3                = 'DMProj3'
  CHARACTER(LEN=*), PARAMETER :: MD_DMP4                = 'DMProj4'

  ! Diagonal Guess
  CHARACTER(LEN=*), PARAMETER :: MD_DGuess              = 'DMDGuess'

  ! MD Algorithm
  CHARACTER(LEN=*), PARAMETER :: MD_AL_OPTION           = 'MDMethod'

  ! Velocity Velet
  CHARACTER(LEN=*), PARAMETER :: MD_AL_VERLET           = 'Verlet'

  ! Symplectic nth order.
  CHARACTER(LEN=*), PARAMETER :: MD_AL_SYMPLECTIC       = 'Symplectic'

  ! Gear PC
  CHARACTER(LEN=*), PARAMETER :: MD_AL_GEAR             = 'Gear'
  INTEGER,          PARAMETER :: GEAR_MD_AL             = 34235424

  ! MD Inputs
  CHARACTER(LEN=*), PARAMETER :: MD_TIME_STEP           = 'DeltaTime'
  CHARACTER(LEN=*), PARAMETER :: MD_MAX_STEP            = 'MaxMDStep'
  CHARACTER(LEN=*), PARAMETER :: MD_DAMP_STEP           = 'MDDampStep'
  CHARACTER(LEN=*), PARAMETER :: MD_NUM_SCF             = 'MDNumSCF'
  CHARACTER(LEN=*), PARAMETER :: MD_DAMPING             = 'MDDamping'

  ! Initial Temp
  CHARACTER(LEN=*), PARAMETER :: MD_INITIAL_TEMP        = 'InitialTemp'

  ! Temperature Rescaling
  CHARACTER(LEN=*), PARAMETER :: MD_TARGET_TEMP         = 'TargetTemp'
  CHARACTER(LEN=*), PARAMETER :: MD_TARGET_ETOTAL       = 'TargetEtotal'
  CHARACTER(LEN=*), PARAMETER :: MD_TSCALE_STEPS        = 'TempScalingSteps'

  ! Thermostats.
  CHARACTER(LEN=*), PARAMETER :: MD_THERMOSTAT           = 'Thermostat'
  CHARACTER(LEN=*), PARAMETER :: MD_THERM_UNSET          = 'ThermUnset'
  CHARACTER(LEN=*), PARAMETER :: MD_THERM_SCALING        = 'Scaling'
  CHARACTER(LEN=*), PARAMETER :: MD_THERM_BERENDSEN      = 'Berendsen'
  CHARACTER(LEN=*), PARAMETER :: MD_THERM_BERENDSEN_ETOT = 'BerendsenEtotal'
  CHARACTER(LEN=*), PARAMETER :: MD_BERENDSEN_TAU        = 'BerendsenTau'

  ! Action Control Theory.
  CHARACTER(LEN=*), PARAMETER :: MD_ACT                 = 'ActionControlTheory'
  CHARACTER(LEN=*), PARAMETER :: MD_ACT_MAX_FORCE_ERROR = "ACTMaxForceError"
  CHARACTER(LEN=*), PARAMETER :: MD_ACT_ALPHA           = "ACTAlpha"
  CHARACTER(LEN=*), PARAMETER :: MD_ACT_BETA            = "ACTBeta"

  ! MC Inputs
  CHARACTER(LEN=*), PARAMETER :: MC_MAX_STEP            = 'MaxMCStep'
  CHARACTER(LEN=*), PARAMETER :: MC_TEMP                = 'MCTemperature'

END MODULE DynamicsKeys
