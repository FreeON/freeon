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

#include "MondoConfig.h"

MODULE MDynamics
  USE Order
  USE SCFs
  USE InOut
  USE MatFunk
  USE Numerics
  USE AtomPairs
  USE ControlStructures
  USE GeomOptKeys
  USE PunchHDF
  USE GlobalScalars
  USE SetXYZ
  USE DynamicsKeys
  USE Utilities
  USE MondoLogger

  IMPLICIT NONE

  TYPE(DBL_VECT) :: MDTime,MDEkin,MDEpot,MDEtot,MDTemp,MDTave
  TYPE(DBL_RNK2) :: MDLinP

  REAL(DOUBLE) :: ACTTargetEnergy, ACTIntegratedEnergyError, ACTGamma
  REAL(DOUBLE) :: ACTActionTransfer, ACTLambda

  INTEGER :: MinMDGeo

CONTAINS
  !--------------------------------------------------------------
  ! Main MD Subroutine
  !--------------------------------------------------------------
  SUBROUTINE MD(C)
    TYPE(Controls)  :: C
    INTEGER         :: I,iSCF,iBAS,iGEO,iCLONE,iATS,iREMOVE
    INTEGER         :: iGEOBegin,iMDStep
    REAL(DOUBLE)    :: Temp, MDDeltaTime
    LOGICAL         :: NewECMD,OrthogDM
    INTEGER         :: oldState

    !--------------------------------------------------------------
    ! Do Molecular Dynamics:Loop over Time Steps

    ! Intitialize
    CALL MondoLog(DEBUG_NONE, "MD", "initializing MD simulation")
    C%Stat%Previous%I = (/0,1,1/)
    iGEO    = 1
    iMDStep = 1

    ! Some output.
    CALL MondoLog(DEBUG_MAXIMUM, "FreeON", "Variables in HDF recycled every "//TRIM(IntToChar(RecycleHDF))//" geometry steps.")

    ! Allocate some memory.
    CALL New(MDTime,C%Geos%Clones)
    CALL New(MDEkin,C%Geos%Clones)
    CALL New(MDEpot,C%Geos%Clones)
    CALL New(MDEtot,C%Geos%Clones)
    CALL New(MDTemp,C%Geos%Clones)
    CALL New(MDTave,C%Geos%Clones)
    CALL New(MDLinP,(/3,C%Geos%Clones/))

    ! Set some stuff to zero.
    MDTime%D = Zero
    MDEkin%D = Zero
    MDEpot%D = Zero
    MDEtot%D = Zero
    MDTemp%D = Zero
    MDTave%D = Zero
    MDLinP%D = Zero

    ! Initial Guess
    IF(C%Opts%Guess==GUESS_EQ_RESTART) THEN

      iGEOBegin = C%Stat%Current%I(3)

      ! Init from old hdf file.
      CALL MondoLog(DEBUG_NONE, "MD", "loading from restart hdf file "//TRIM(C%Nams%RFile))
      HDFFileID=OpenHDF(C%Nams%RFile)
      DO iCLONE=1,C%Geos%Clones
        HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
        CALL Get(iMDStep,"iMDStep")
        CALL Get(MDTime%D(iCLONE),"MDTime")
        CALL Get(MDDeltaTime, "MDDeltaTime", Tag_O = TRIM(IntToChar(iGEOBegin)))
        CALL Get(C%Dyns%MDGuess,"MDGuess")
        CALL CloseHDFGroup(HDF_CurrentID)
      ENDDO
      CALL CloseHDF(HDFFileID)

      ! Save to Current HDF
      CALL MondoLog(DEBUG_NONE, "MD", "saving to current hdf file "//TRIM(C%Nams%HFile))
      HDFFileID=OpenHDF(C%Nams%HFile)
      DO iCLONE=1,C%Geos%Clones
        HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
        CALL Put(iMDStep,"iMDStep")
        CALL Put(MDTime%D(iCLONE),"MDTime")
        CALL Put(MDDeltaTime, "MDDeltaTime", Tag_O = TRIM(IntToChar(iGEOBegin)))
        CALL Put(.TRUE.,"DoingMD")
        CALL Put(C%Dyns%MDGuess,"MDGuess")
        CALL CloseHDFGroup(HDF_CurrentID)
      ENDDO
      CALL CloseHDF(HDFFileID)

      ! Determine the Begin Step
      IF(iMDStep > iGEOBegin) THEN
        CALL MondoLog(DEBUG_NONE, "MD", "iMDStep ("//TRIM(IntToChar(iMDStep)) &
          //") > iGEOBegin ("//TRIM(IntToChar(iGEOBegin))//")")
        iMDStep = iGEOBegin
        IF(C%Geos%Clones > 1) THEN
          CALL MondoLog(DEBUG_NONE, "MD", "[FIXME] this might fail if you are using clones")
        ENDIF
        MDTime%D(1) = MDTime%D(1)-C%Dyns%DTime

        ! Since we are recalculating the last time step, we will add 1 step to
        ! MDMaxSteps to account for that additional step.
        CALL MondoLog(DEBUG_NONE, "MD", "I will recalculate the last time step")
        C%Dyns%MDMaxSteps = C%Dyns%MDMaxSteps+1
      ENDIF

      ! Determine iREMOVE and MinMDGeo
      CALL CalculateMDGeo(C%Dyns,iREMOVE,MinMDGeo)

      ! Print some stuff out.
      CALL MondoLog(DEBUG_NONE, "MD", "restarting MD calculation")
      CALL MondoLog(DEBUG_NONE, "MD", "iMDStep    = "//TRIM(IntToChar(iMDStep)))
      CALL MondoLog(DEBUG_NONE, "MD", "iGEOBegin  = "//TRIM(IntToChar(iGEOBegin)))
      CALL MondoLog(DEBUG_NONE, "MD", "MDTime     = "//TRIM(FltToChar(MDTime%D(1)*InternalTimeToFemtoseconds))//" fs")
      CALL MondoLog(DEBUG_NONE, "MD", "MinMDGeo   = "//TRIM(IntToChar(MinMDGeo)))
      CALL MondoLog(DEBUG_NONE, "MD", "iREMOVE    = "//TRIM(IntToChar(iREMOVE)))
      CALL MondoLog(DEBUG_NONE, "MD", "MDMaxSteps = "//TRIM(IntToChar(C%Dyns%MDMaxSteps)))

      ! copy .DOsave (.FOsave) and .DOPsave (.FOPsave) matrice to their new names
      CALL CopyRestart(C,iGEOBegin,MinMDGeo)
    ELSE
      ! Init the Time
      MDTime%D(:) = Zero
      MDDeltaTime = C%Dyns%DTime
      HDFFileID=OpenHDF(C%Nams%HFile)
      DO iCLONE=1,C%Geos%Clones
        HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
        CALL Put(MDTime%D(iCLONE),"MDTime")

        CALL Put(.TRUE.,"DoingMD")
        CALL Put(C%Dyns%MDGuess,"MDGuess")
        CALL CloseHDFGroup(HDF_CurrentID)
      ENDDO
      CALL CloseHDF(HDFFileID)

      ! Determine the Begin Step
      iGEOBegin = iMDStep

      ! Determine iREMOVE and MinMDGeo
      CALL CalculateMDGeo(C%Dyns,iREMOVE,MinMDGeo)

      ! Print some stuff out.
      CALL MondoLog(DEBUG_NONE, "MD", "new MD calculation")
      CALL MondoLog(DEBUG_NONE, "MD", "iMDStep     = "//TRIM(IntToChar(iMDStep)))
      CALL MondoLog(DEBUG_NONE, "MD", "iGEOBegin   = "//TRIM(IntToChar(iGEOBegin)))
      CALL MondoLog(DEBUG_NONE, "MD", "MDTime      = "//TRIM(FltToChar(MDTime%D(1)*InternalTimeToFemtoseconds))//" fs")
      CALL MondoLog(DEBUG_NONE, "MD", "MDDeltaTime = "//TRIM(FltToChar(C%Dyns%DTime*InternalTimeToFemtoseconds))//" fs")
      CALL MondoLog(DEBUG_NONE, "MD", "MinMDGeo    = "//TRIM(IntToChar(MinMDGeo)))
      CALL MondoLog(DEBUG_NONE, "MD", "iREMOVE     = "//TRIM(IntToChar(iREMOVE)))
      CALL MondoLog(DEBUG_NONE, "MD", "MDMaxSteps  = "//TRIM(IntToChar(C%Dyns%MDMaxSteps)))
    ENDIF

    ! Do MD
    DO iGEO = iGEOBegin, C%Dyns%MDMaxSteps+iGEOBegin-1

      ! Write time step size to hdf.
      HDFFileID = OpenHDF(C%Nams%HFile)
      DO iCLONE = 1, C%Geos%Clones
        HDF_CurrentID = OpenHDFGroup(HDFFileID, "Clone #"//TRIM(IntToChar(iCLONE)))
        CALL Put(MDDeltaTime, "MDDeltaTime", Tag_O = TRIM(IntToChar(iGEO)))
        CALL CloseHDFGroup(HDF_CurrentID)
      ENDDO
      CALL CloseHDF(HDFFileID)

      ! Initialize with temperature distribution.
      IF((iGEO == 1) .AND. C%Dyns%Initial_Temp) THEN
        CALL SetTempMaxBoltDist(C,C%Dyns%TempInit)
      ENDIF

      ! For iMDStep <= MinMDGeo make Diagonal or Core Guess
      IF(iGEO <= MinMDGeo) THEN
        CALL MondoLog(DEBUG_NONE, "MD", "iGEO <= MinMDGeo, guessing...")
        DO iBAS=1,C%Sets%NBSets
          IF(iBAS==1) THEN
            C%Opts%GuessToP2Use=SCF_SUPERPOSITION
          ELSE
            C%Opts%GuessTOP2Use=SCF_BASISSETSWITCH
          ENDIF
          CALL GeomArchive(iBAS,iGEO,C%Nams,C%Opts,C%Sets,C%Geos)
          CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
          CALL SCF(iBAS,iGEO,C)
        ENDDO

        ! Copy Last SCF to P(iGEO)_save
        CALL CopyMatrices(C,.TRUE.)

        ! Copy Last SCF to P(iGEO)_save_tilde
        CALL CopyMatrices(C,.FALSE.)
      ELSE
        ! Normal MD-SCF
        C%Opts%GuessTOP2Use = C%Dyns%MDGuess

        IF(C%Dyns%MDNumSCF < 0) THEN
          CALL MondoLog(DEBUG_NONE, "MD", "MDNumSCF not set")
        ELSE
#if defined MD_DEBUG
          C%Opts%MinSCF = C%Dyns%MDNumSCF
          C%Opts%MaxSCF = HAVE_MAX_SCF
#else
          C%Opts%MinSCF = C%Dyns%MDNumSCF
          C%Opts%MaxSCF = C%Dyns%MDNumSCF
#endif
          CALL MondoLog(DEBUG_NONE, "MD", "MDNumSCF set to "//TRIM(IntToChar(C%Dyns%MDNumSCF)))
        ENDIF
        CALL MondoLog(DEBUG_NONE, "MD", "MinSCF = "//TRIM(IntToChar(C%Opts%MinSCF)))
        CALL MondoLog(DEBUG_NONE, "MD", "MaxSCF = "//TRIM(IntToChar(C%Opts%MaxSCF)))

        iBAS=C%Sets%NBSets
        CALL GeomArchive(iBAS,iGEO,C%Nams,C%Opts,C%Sets,C%Geos)
        CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
        CALL SCF(iBAS,iGEO,C)

        ! Copy Last SCF to P(iGEO)_save
#if defined MD_DEBUG
        IF(C%Dyns%MDNumSCF < 0) THEN
          CALL MondoLog(DEBUG_NONE, "MD", "MDNumSCF not set")
        ELSE
          CALL MondoLog(DEBUG_NONE, "MD", "hardwiring density matrix to use result from "//TRIM(IntToChar(C%Dyns%MDNumSCF+1))//" SCFs")
          oldState = C%Stat%Current%I(1)
          C%Stat%Current%I(1) = C%Dyns%MDNumSCF
        ENDIF
#endif
        CALL CopyMatrices(C,.FALSE.)
#if defined MD_DEBUG
        IF(C%Dyns%MDNumSCF >= 0) THEN
          CALL MondoLog(DEBUG_NONE, "MD", "resetting back to "//TRIM(IntToChar(oldState)))
          C%Stat%Current%I(1) = oldState
        ENDIF
#endif
      ENDIF

      ! Calculate the Forces.
      !
      ! Hardwire force to use density after 2 SCFs.
#if defined MD_DEBUG
      IF(C%Dyns%MDNumSCF < 0) THEN
        CALL MondoLog(DEBUG_NONE, "MD", "MDNumSCF not set")
      ELSE
        IF(iGEO > MinMDGeo) THEN
          CALL MondoLog(DEBUG_NONE, "MD", "hardwiring force to use result from "//TRIM(IntToChar(C%Dyns%MDNumSCF+1))//" SCFs")
          oldState = C%Stat%Current%I(1)
          C%Stat%Current%I(1) = C%Dyns%MDNumSCF
        ENDIF
      ENDIF
#endif
      CALL Force(C%Sets%NBSets,iGEO,C%Nams,C%Opts,C%Stat,C%Geos,C%Sets,C%MPIs)
#if defined MD_DEBUG
      IF(C%Dyns%MDNumSCF >= 0) THEN
        IF(iGEO > MinMDGeo) THEN
          CALL MondoLog(DEBUG_NONE, "MD", "resetting back to "//TRIM(IntToChar(oldState)))
          C%Stat%Current%I(1) = oldState
        ENDIF
      ENDIF
#endif

      ! Move the Atoms, apply the thermostats and print some output.
      IF(C%Dyns%MDAlgorithm == MD_AL_VERLET) THEN
        IF(C%Dyns%Const_Temp .AND. .NOT.C%Dyns%Const_Press) THEN
          CALL Halt('Constant Temperature MD Not Implemented')
          CALL MDVerlet_NVT(C,iGEO)
        ELSEIF(C%Dyns%Const_Press .AND. .NOT.C%Dyns%Const_Temp) THEN
          CALL Halt('Constant Presure MD Not Implemented')
          CALL MDVerlet_NPH(C,iGEO)
        ELSEIF(C%Dyns%Const_Press.AND.C%Dyns%Const_Temp) THEN
          CALL Halt('Constant Presure/Constant Temperature MD Not Implemented')
          CALL MDVerlet_NPT(C,iGEO)
        ELSEIF(.NOT.C%Dyns%Const_Press .AND. .NOT.C%Dyns%Const_Temp) THEN
          CALL MDVerlet_NVE(C,iGEO)
        ENDIF
      ELSEIF(C%Dyns%MDAlgorithm == MD_AL_SYMPLECTIC) THEN
        CALL MDSymplectic_4th_Order_NVE(C, iGeo)
      ELSEIF(C%Dyns%MDAlgorithm == MD_AL_GEAR) THEN
        CALL Halt('Gear MD Not Implemented')
      ENDIF

      ! If Parallel Rep, switch Geometries
      IF(C%Dyns%Parallel_Rep) THEN
        CALL Halt('Parallel Replicate MD Not Implemented')
      ENDIF

      ! Refresh the Time and other Stuff
      HDFFileID=OpenHDF(C%Nams%HFile)
      DO iCLONE=1,C%Geos%Clones
        iMDStep = iMDStep+1

        IF(C%Dyns%MDAlgorithm == MD_AL_SYMPLECTIC) THEN
          MDTime%D(iCLONE) = MDTime%D(iCLONE)+C%Dyns%DTime/4.0D0
        ELSE
          MDTime%D(iCLONE) = MDTime%D(iCLONE)+C%Dyns%DTime
        ENDIF

        HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
        CALL Put(iMDStep,'iMDStep')
        CALL Put(MDTime%D(iCLONE),"MDTime")
        CALL CloseHDFGroup(HDF_CurrentID)
        CALL MondoLog(DEBUG_NONE, "MD", "done with time step, incrementing counters and putting to hdf")
        CALL MondoLog(DEBUG_NONE, "MD", "iMDStep = "//TRIM(IntToChar(iMDStep)))
        CALL MondoLog(DEBUG_NONE, "MD", "MDTime = "//TRIM(FltToChar(MDTime%D(1)*InternalTimeToFemtoseconds))//" fs")
        CALL MondoLog(DEBUG_NONE, "MD", "MDEtot = "//TRIM(DblToChar(MDEtot%D(iCLONE)))//" Hartree")
        CALL MondoLog(DEBUG_NONE, "MD", "MDEkin = "//TRIM(DblToChar(MDEkin%D(iCLONE)))//" Hartree")
        CALL MondoLog(DEBUG_NONE, "MD", "MDEpot = "//TRIM(DblToChar(MDEpot%D(iCLONE)))//" Hartree")
        CALL MondoLog(DEBUG_NONE, "MD", "T = "//TRIM(DblToChar(MDTemp%D(iCLONE)))//" K")
      ENDDO
      CALL CloseHDF(HDFFileID)

      ! Remove old Stuff from Scratch
      IF(C%Stat%Current%I(3)-iREMOVE > 0) THEN
        CALL MondoLog(DEBUG_NONE, "MD", "removing outdated densities")
        IF(doCleanScratch) THEN
          CALL CleanScratch(C,iGEO-iREMOVE-1,.TRUE.)
        ENDIF
      ENDIF

    ENDDO

    ! Clean up.
    CALL Delete(MDTime)
    CALL Delete(MDEkin)
    CALL Delete(MDEpot)
    CALL Delete(MDEtot)
    CALL Delete(MDTemp)
    CALL Delete(MDTave)
    CALL Delete(MDLinP)

  END SUBROUTINE MD

  !--------------------------------------------------------------
  ! The Velocity Verlet Algorithmn : NVE
  !--------------------------------------------------------------
  SUBROUTINE MDVerlet_NVE(C,iGEO)
    TYPE(Controls)            :: C
    INTEGER                   :: iGEO
    INTEGER                   :: iCLONE,iATS,AOut
    INTEGER                   :: numberConstrainedAtoms
    REAL(DOUBLE)              :: Mass,totalMass,dT,dT2,dTSq2,Time,Dist,averageTorqueNorm
    REAL(DOUBLE),DIMENSION(3) :: Pos,Vel,Acc,torque,PosSave,VelSave,centerOfMass
    REAL(DOUBLE)              :: v_scale

    ! Initialize
    dT    = C%Dyns%DTime
    dT2   = Half*dT
    dTSq2 = Half*dT*dT

    ! Count the number of constrained atoms.
    numberConstrainedAtoms = 0
    centerOfMass = Zero
    DO iCLONE = 1,C%Geos%Clones
      DO iATS = 1,C%Geos%Clone(iCLONE)%NAtms
        IF(C%Geos%Clone(iCLONE)%CConstrain%I(iATS) /= 0) THEN
          numberConstrainedAtoms = numberConstrainedAtoms + 1
          centerOfMass = C%Geos%Clone(iCLONE)%Carts%D(1:3, iATS)
        ENDIF
      ENDDO
    ENDDO

    ! Set the total torque to zero. In the case of no constrained atom, we take
    ! the center of mass as the center of rotation. In the case of 1 constrained
    ! atom, we take that atom position as the center of rotation. In the case of
    ! more than 1 constrained atom, we don't worry about the spinning ice-cube
    ! correction.
    IF(numberConstrainedAtoms == 0) THEN
      totalMass = Zero
      DO iATS = 1, C%Geos%Clone(1)%NAtms
        centerOfMass(1:3) = centerOfMass(1:3) + C%Geos%Clone(1)%AtMss%D(iATS)*C%Geos%Clone(1)%Carts%D(1:3, iAts)
        totalMass = totalMass + C%Geos%Clone(1)%AtMss%D(iATS)
      ENDDO
      centerOfMass(1:3) = centerOfMass(1:3)/totalMass
    ENDIF

    CALL MondoLog(DEBUG_NONE, "MD:Verlet", "CM = [ "// &
      TRIM(DblToChar(centerOfMass(1)))//" "// &
      TRIM(DblToChar(centerOfMass(2)))//" "// &
      TRIM(DblToChar(centerOfMass(3)))//" ]")

    IF(numberConstrainedAtoms < 2) THEN
      torque = Zero
      DO iATS = 1, C%Geos%Clone(1)%NAtms
        !CALL MondoLog(DEBUG_NONE, "MD:Verlet", "F("//TRIM(IntToChar(iATS))//") = [ " &
        !  //TRIM(DblToChar(-C%Geos%Clone(1)%Gradients%D(1, iATS)))//" " &
        !  //TRIM(DblToChar(-C%Geos%Clone(1)%Gradients%D(2, iATS)))//" " &
        !  //TRIM(DblToChar(-C%Geos%Clone(1)%Gradients%D(3, iATS)))//" ]")
        torque(1:3) = torque(1:3) - CROSS_PRODUCT(C%Geos%Clone(1)%Gradients%D(1:3, iATS), C%Geos%Clone(1)%Carts%D(1:3, iATS)-centerOfMass(1:3))
      ENDDO

      CALL MondoLog(DEBUG_NONE, "MD:Verlet", "spinning ice-cube correction: [ " &
        //TRIM(DblToChar(torque(1)))//" " &
        //TRIM(DblToChar(torque(2)))//" " &
        //TRIM(DblToChar(torque(3)))//" ]")

      averageTorqueNorm = VABS(torque/DBLE(C%Geos%Clone(1)%NAtms-numberConstrainedAtoms))
      CALL MondoLog(DEBUG_NONE, "MD:Verlet", "norm(average torque) = "//TRIM(DblToChar(averageTorqueNorm)))

      DO iATS = 1, C%Geos%Clone(1)%NAtms
        Acc(1:3) = CROSS_PRODUCT(torque, C%Geos%Clone(1)%Carts%D(1:3, iATS)-centerOfMass(1:3))
        IF(VABS(Acc) > 1D-12) THEN
          Acc(1:3) = Acc(1:3)/VABS(Acc)
          Acc(1:3) = Acc(1:3)*averageTorqueNorm/VABS(C%Geos%Clone(1)%Carts%D(1:3, iATS)-centerOfMass(1:3))
          C%Geos%Clone(1)%Gradients%D(1:3, iATS) = C%Geos%Clone(1)%Gradients%D(1:3, iATS)-Acc(1:3)
        ENDIF
      ENDDO

      ! Check torque.
      torque = Zero
      DO iATS = 1, C%Geos%Clone(1)%NAtms
        torque(1:3) = torque(1:3) - CROSS_PRODUCT(C%Geos%Clone(1)%Gradients%D(1:3, iATS), C%Geos%Clone(1)%Carts%D(1:3, iATS)-centerOfMass(1:3))
      ENDDO

      CALL MondoLog(DEBUG_NONE, "MD:Verlet", "total torque after spinning ice-cube correction: [ " &
        //TRIM(DblToChar(torque(1)))//" " &
        //TRIM(DblToChar(torque(2)))//" " &
        //TRIM(DblToChar(torque(3)))//" ]")
    ENDIF

    ! Set the Sum of the Forces equal to zero in case of no constraint atoms.
    IF(numberConstrainedAtoms == 0) THEN
      DO iCLONE = 1,C%Geos%Clones
        Acc(1:3) = Zero
        DO iATS = 1,C%Geos%Clone(iCLONE)%NAtms
          IF(C%Geos%Clone(iCLONE)%CConstrain%I(iATS) == 0) THEN
            Acc(1:3) = Acc(1:3)+C%Geos%Clone(iCLONE)%Gradients%D(1:3,iATS)
          ENDIF
        ENDDO

        CALL MondoLog(DEBUG_NONE, "MD:Verlet", "flying ice-cube correction: [ " &
          //TRIM(DblToChar(Acc(1)))//" " &
          //TRIM(DblToChar(Acc(2)))//" " &
          //TRIM(DblToChar(Acc(3)))//" ]")

        ! Make sure the total force on the system is zero.
        CALL MondoLog(DEBUG_NONE, "MD:Verlet", "setting total force to zero")
        DO iATS = 1,C%Geos%Clone(iCLONE)%NAtms
          IF(C%Geos%Clone(iCLONE)%CConstrain%I(iATS) == 0) THEN
            C%Geos%Clone(iCLONE)%Gradients%D(1:3,iATS) = C%Geos%Clone(iCLONE)%Gradients%D(1:3,iATS)-Acc(1:3)*C%Geos%Clone(iCLONE)%AtMss%D(iATS)/totalMass
          ENDIF
        ENDDO

        ! Check torque.
        torque = Zero
        DO iATS = 1, C%Geos%Clone(1)%NAtms
          torque(1:3) = torque(1:3) - CROSS_PRODUCT(C%Geos%Clone(1)%Gradients%D(1:3, iATS), C%Geos%Clone(1)%Carts%D(1:3, iATS)-centerOfMass(1:3))
        ENDDO

        CALL MondoLog(DEBUG_NONE, "MD:Verlet", "total torque after flying ice-cube correction: [ " &
          //TRIM(DblToChar(torque(1)))//" " &
          //TRIM(DblToChar(torque(2)))//" " &
          //TRIM(DblToChar(torque(3)))//" ]")
      ENDDO
    ENDIF

    ! Update the Velocity if not the first step
    IF(iGEO .NE. 1) THEN
      CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", "updating velocity")
      DO iCLONE = 1,C%Geos%Clones
        DO iATS = 1,C%Geos%Clone(iCLONE)%NAtms
          IF(C%Geos%Clone(iCLONE)%CConstrain%I(iATS) == 0) THEN
            Mass     =  C%Geos%Clone(iCLONE)%AtMss%D(iATS)
            Vel(1:3) =  C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS)
            Acc(1:3) = -C%Geos%Clone(iCLONE)%Gradients%D(1:3,iATS)/Mass

            ! Velocity: v(t) = v(t-dT/2)+a(t)dT/2
            Vel(1:3) = Vel(1:3) + Acc(1:3)*dT2
            C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS) = Vel(1:3)
          ENDIF
        ENDDO
      ENDDO
    ELSE
      CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", "first step, not updating velocity")
    ENDIF

    ! Reset the linear momentum, compute the kinetic energy and average kinectic energy.
    DO iCLONE = 1,C%Geos%Clones
      ! Calculate kinetic energy and temperature, update average temperature.
      CALL CalculateMDKin(C%Geos%Clone(iCLONE),MDEkin%D(iCLONE),MDTemp%D(iCLONE))
      MDTave%D(iCLONE) = (DBLE(iGEO-1)/DBLE(iGEO))*MDTave%D(iCLONE) +(One/DBLE(iGEO))*MDTemp%D(iCLONE)

      ! Store Potential and Total Energy
      CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", "ETotal = "//TRIM(DblToChar(C%Geos%Clone(iCLONE)%ETotal))//" Hartree")
      CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", "ETotalPerSCF = "//TRIM(DblVectToChar(C%Geos%Clone(iCLONE)%ETotalPerSCF, &
        (/ 0, C%Stat%Current%I(1) /)))//" Hartree")
      CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", "EKin = "//TRIM(DblToChar(MDEkin%D(iCLONE)))//" Hartree")
      CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", "T = "//TRIM(DblToChar(MDTemp%D(iCLONE)))//" K")
#if defined MD_DEBUG
      IF(C%Dyns%MDNumSCF < 0) THEN
        CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", "MDNumSCF not set")
      ELSE
        IF(iGEO > MinMDGeo) THEN
          CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", "hardwiring potential energy to use result from "//TRIM(IntToChar(C%Dyns%MDNumSCF+1))//" SCFs")
          MDEpot%D(iCLONE) = C%Geos%Clone(iCLONE)%ETotalPerSCF%D(C%Dyns%MDNumSCF)
        ELSE
          MDEpot%D(iCLONE) = C%Geos%Clone(iCLONE)%ETotal
        ENDIF
      ENDIF
#else
      MDEpot%D(iCLONE) = C%Geos%Clone(iCLONE)%ETotal
#endif

      ! Calculate total energy.
      MDEtot%D(iCLONE) = MDEpot%D(iCLONE) + MDEkin%D(iCLONE)
    ENDDO

    ! Rescaling thermostat.
    IF(C%Dyns%Temp_Scaling) THEN
      IF(MOD(iGEO,C%Dyns%RescaleInt)==0) THEN
        CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", 'Rescaling Temperature')
        CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", 'MD temperature     = '//TRIM(DblToChar(MDTemp%D(1)))//" K")
        CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", 'Target temperature = '//TRIM(DblToChar(C%Dyns%TargetTemp)))
        DO iCLONE = 1,C%Geos%Clones
          CALL RescaleVelocity(C%Geos%Clone(iCLONE),MDTemp%D(iCLONE),C%Dyns%TargetTemp)
        ENDDO
      ENDIF
    ENDIF

    ! Berendsen thermostat.
    IF(C%Dyns%Thermostat == MD_THERM_BERENDSEN) THEN
      CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", "Applying Berendsen thermostat")
      CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", "MD temperature     = "//TRIM(DblToChar(MDTemp%D(1)))//" K")
      CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", "Target temperature = "//TRIM(DblToChar(C%Dyns%TargetTemp))//" K")
      DO iCLONE = 1, C%Geos%Clones
        CALL BerendsenThermostat(C%Geos%Clone(iCLONE), MDTemp%D(iCLONE), C%Dyns%TargetTemp, C%Dyns%DTime, C%Dyns%BerendsenTau, v_scale)
        C%Dyns%BerendsenVScale = v_scale

        ! Store v_scale in hdf.
        HDFFileID=OpenHDF(C%Nams%HFile)
        HDF_CurrentID = OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
        CALL Put(v_scale, "v_scale", Tag_O=TRIM(IntToChar(iGEO)))
        CALL CloseHDFGroup(HDF_CurrentID)
        CALL CloseHDF(HDFFileID)
        CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", "putting v_scale("//TRIM(IntToChar(iGEO))//") to hdf = "//TRIM(DblToChar(v_scale)))
      ENDDO

    ! Berendsen thermostat for Etotal.
    ELSEIF(C%Dyns%Thermostat == MD_THERM_BERENDSEN_ETOT) THEN

      IF(iGEO == 1 .AND. (.NOT.C%Dyns%Energy_Scaling_Set)) THEN
        C%Dyns%TargetEtotal = MDEtot%D(1)
        CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", "Using total energy of first MD step as target energy = "// &
          TRIM(DblToChar(C%Dyns%TargetEtotal*au2eV))//" eV")
      ENDIF

      CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", "Applying Berendsen thermostat for E_total")
      CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", "MD E_total     = "//TRIM(DblToChar(MDEtot%D(1)*au2eV))//" eV")
      CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", "Target E_total = "//TRIM(DblToChar(C%Dyns%TargetEtotal*au2eV))//" eV")
      DO iCLONE = 1, C%Geos%Clones
        CALL BerendsenThermostatTotalEnergy(C%Geos%Clone(iCLONE), MDEpot%D(iCLONE), MDEkin%D(iCLONE), C%Dyns%TargetEtotal, &
          C%Dyns%DTime, C%Dyns%BerendsenTau, v_scale)
        C%Dyns%BerendsenVScale = v_scale

        ! Store v_scale in hdf.
        HDFFileID=OpenHDF(C%Nams%HFile)
        HDF_CurrentID = OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
        CALL Put(v_scale, "v_scale", Tag_O=TRIM(IntToChar(iGEO)))
        CALL CloseHDFGroup(HDF_CurrentID)
        CALL CloseHDF(HDFFileID)
        CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", "putting v_scale("//TRIM(IntToChar(iGEO))//") to hdf = "//TRIM(DblToChar(v_scale)))
      ENDDO

    ELSE

      DO iCLONE = 1, C%Geos%Clones
        ! Let's punch something so we don't die in P2Use.
        HDFFileID=OpenHDF(C%Nams%HFile)
        HDF_CurrentID = OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
        CALL Put(0D0, "v_scale", Tag_O=TRIM(IntToChar(iGEO)))
        CALL CloseHDFGroup(HDF_CurrentID)
        CALL CloseHDF(HDFFileID)
        CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", "putting dummy v_scale("//TRIM(IntToChar(iGEO))//") to hdf = "//TRIM(DblToChar(0D0)))
      ENDDO

    ENDIF

    ! Generate Output
    CALL MondoLog(DEBUG_NONE, "MD:Verlet_NVE", "writing coordinates to file")
    CALL OutputMD(C,iGEO)

    ! Update the Positions
    DO iCLONE=1,C%Geos%Clones
      DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
        IF(C%Geos%Clone(iCLONE)%CConstrain%I(iATS) == 0)THEN
          Mass     =  C%Geos%Clone(iCLONE)%AtMss%D(iATS)
          Pos(1:3) =  C%Geos%Clone(iCLONE)%Carts%D(1:3,iATS)
          Vel(1:3) =  C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS)
          Acc(1:3) = -C%Geos%Clone(iCLONE)%Gradients%D(1:3,iATS)/Mass

          ! Position: r(t+dT)= r(t)+v(t)*dT+a(t)*dT*dT/2
          Pos(1:3) = Pos(1:3) + Vel(1:3)*dT + Acc(1:3)*dTSq2

          ! Update.
          C%Geos%Clone(iCLONE)%Carts%D(1:3,iATS) = Pos(1:3)
        ENDIF
      ENDDO
    ENDDO

    ! Update the Velocity.
    DO iCLONE=1,C%Geos%Clones
      DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
        IF(C%Geos%Clone(iCLONE)%CConstrain%I(iATS) == 0)THEN
          Mass     =  C%Geos%Clone(iCLONE)%AtMss%D(iATS)
          Vel(1:3) =  C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS)
          Acc(1:3) = -C%Geos%Clone(iCLONE)%Gradients%D(1:3,iATS)/Mass

          ! Velocity: v(t) = v(t-dT/2)+a(t)dT/2
          Vel(1:3) = Vel(1:3) + Acc(1:3)*dT2
          C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS) = Vel(1:3)
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE MDVerlet_NVE

  !--------------------------------------------------------------
  ! Symplectic 4th order integration.
  !--------------------------------------------------------------
  SUBROUTINE MDSymplectic_4th_Order_NVE(C, iGEO)
    TYPE(Controls)            :: C
    INTEGER                   :: iGEO
    INTEGER                   :: iCLONE,iATS,AOut
    INTEGER                   :: numberUnconstrainedAtoms
    INTEGER                   :: m_step
    REAL(DOUBLE)              :: Mass,dT,dT2,dTSq2,Time,Dist
    REAL(DOUBLE),DIMENSION(3) :: Pos,Vel,Acc,PosSave,VelSave
    REAL(DOUBLE)              :: v_scale

    ! Initialize
    dT    = C%Dyns%DTime
    dT2   = Half*dT
    dTSq2 = Half*dT*dT

    ! Set the Sum of the Forces equal to zero
    numberUnconstrainedAtoms = 0
    DO iCLONE = 1,C%Geos%Clones
      Acc(1:3) = Zero
      DO iATS = 1,C%Geos%Clone(iCLONE)%NAtms
        IF(C%Geos%Clone(iCLONE)%CConstrain%I(iATS) == 0) THEN
          Acc(1:3) = Acc(1:3)+C%Geos%Clone(iCLONE)%Gradients%D(1:3,iATS)
          numberUnconstrainedAtoms = numberUnconstrainedAtoms + 1
        ENDIF
      ENDDO

      ! Calculate total acceleration.
      Acc(1:3) = Acc(1:3)/DBLE(numberUnconstrainedAtoms)

      CALL MondoLog(DEBUG_NONE, "MD:Symplectic", "flying ice-cube correction: [ " &
        //TRIM(DblToChar(Acc(1)))//" " &
        //TRIM(DblToChar(Acc(2)))//" " &
        //TRIM(DblToChar(Acc(3)))//" ]")

      IF(numberUnconstrainedAtoms == C%Geos%Clone(iCLONE)%NAtms) THEN
        ! Make sure the total force on the system is zero.
        DO iATS = 1,C%Geos%Clone(iCLONE)%NAtms
          IF(C%Geos%Clone(iCLONE)%CConstrain%I(iATS) == 0) THEN
            C%Geos%Clone(iCLONE)%Gradients%D(1:3,iATS) = C%Geos%Clone(iCLONE)%Gradients%D(1:3,iATS)-Acc(1:3)
          ENDIF
        ENDDO
      ENDIF
    ENDDO

    ! Get the symplectic counter.
    m_step = MOD(iGEO-1,4)+1

    IF(m_step == 1) THEN

      ! Reset the linear momentum, compute the kinectic energy and ave kinectic energy
      DO iCLONE=1,C%Geos%Clones
        ! Calculate kinectic energy and temp, update ave temp
        CALL CalculateMDKin(C%Geos%Clone(iCLONE),MDEkin%D(iCLONE),MDTemp%D(iCLONE))
        MDTave%D(iCLONE) = (DBLE(iGEO-1)/DBLE(iGEO))*MDTave%D(iCLONE)+(One/DBLE(iGEO))*MDTemp%D(iCLONE)

        ! Store Potential and Total Energy
        CALL MondoLog(DEBUG_NONE, "MD:Symplectic", "T = "//TRIM(DblToChar(MDTemp%D(iCLONE)))//" K")
        CALL MondoLog(DEBUG_NONE, "MD:Symplectic", "ETotal = "//TRIM(DblToChar(C%Geos%Clone(iCLONE)%ETotal)))
        CALL MondoLog(DEBUG_NONE, "MD:Symplectic", "ETotalPerSCF = "//TRIM(DblVectToChar(C%Geos%Clone(iCLONE)%ETotalPerSCF, (/ 0, C%Stat%Current%I(1) /))))
#if defined MD_DEBUG
        IF(C%Dyns%MDNumSCF < 0) THEN
          CALL MondoLog(DEBUG_NONE, "MD:Symplectic", "MDNumSCF not set")
        ELSE
          IF(iGEO > MinMDGeo) THEN
            CALL MondoLog(DEBUG_NONE, "MD:Symplectic", "hardwiring potential energy to use result from "//TRIM(IntToChar(C%Dyns%MDNumSCF+1))//" SCFs")
            MDEpot%D(iCLONE) = C%Geos%Clone(iCLONE)%ETotalPerSCF%D(C%Dyns%MDNumSCF)
          ELSE
            MDEpot%D(iCLONE) = C%Geos%Clone(iCLONE)%ETotal
          ENDIF
        ENDIF
#else
        MDEpot%D(iCLONE) = C%Geos%Clone(iCLONE)%ETotal
#endif

        ! Store Potential and Total Energy
        MDEtot%D(iCLONE) = MDEpot%D(iCLONE) + MDEkin%D(iCLONE)
      ENDDO

      ! Thermostats
      IF(C%Dyns%Temp_Scaling) THEN
        IF(MOD(iGEO,C%Dyns%RescaleInt)==0) THEN
          CALL MondoLog(DEBUG_NONE, "MD:Symlectic", 'Rescaling Temperature')
          CALL MondoLog(DEBUG_NONE, "MD:Symlectic", 'MD temperature     = '//TRIM(DblToChar(MDTemp%D(1))))
          CALL MondoLog(DEBUG_NONE, "MD:Symlectic", 'Target temperature = '//TRIM(DblToChar(C%Dyns%TargetTemp)))
          DO iCLONE=1,C%Geos%Clones
            CALL RescaleVelocity(C%Geos%Clone(iCLONE),MDTemp%D(iCLONE),C%Dyns%TargetTemp)
          ENDDO
        ENDIF
      ENDIF

      ! Berendsen thermostat.
      IF(C%Dyns%Thermostat == MD_THERM_BERENDSEN) then
        CALL MondoLog(DEBUG_NONE, "MD:Symlectic", "Applying Berendsen thermostat")
        CALL MondoLog(DEBUG_NONE, "MD:Symlectic", "MD temperature     = "//TRIM(DblToChar(MDTemp%D(1))))
        CALL MondoLog(DEBUG_NONE, "MD:Symlectic", "Target temperature = "//TRIM(DblToChar(C%Dyns%TargetTemp)))
        DO iCLONE = 1, C%Geos%Clones
          CALL BerendsenThermostat(C%Geos%Clone(iCLONE), MDTemp%D(iCLONE), C%Dyns%TargetTemp, C%Dyns%DTime, C%Dyns%BerendsenTau, v_scale)
          C%Dyns%BerendsenVScale = v_scale

          ! Store v_scale in hdf.
          HDFFileID=OpenHDF(C%Nams%HFile)
          HDF_CurrentID = OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
          CALL Put(v_scale, "v_scale", Tag_O=TRIM(IntToChar(iGEO)))
          CALL CloseHDFGroup(HDF_CurrentID)
          CALL CloseHDF(HDFFileID)
          CALL MondoLog(DEBUG_NONE, "MD:Symplectic", "putting v_scale("//TRIM(IntToChar(iGEO))//") to hdf = "//TRIM(DblToChar(v_scale)))
          CALL MondoLog(DEBUG_NONE, "MD:Symplectic", "a("//TRIM(IntToChar(m_step))//") = "//TRIM(FltToChar(Symplectic_4th_Order_a(m_step))))
          CALL MondoLog(DEBUG_NONE, "MD:Symplectic", "b("//TRIM(IntToChar(m_step))//") = "//TRIM(FltToChar(Symplectic_4th_Order_b(m_step))))
        ENDDO

      ELSE

        DO iCLONE = 1, C%Geos%Clones
          ! Let's punch something so we don't die in P2Use.
          HDFFileID=OpenHDF(C%Nams%HFile)
          HDF_CurrentID = OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
          CALL Put(0D0, "v_scale", Tag_O=TRIM(IntToChar(iGEO)))
          CALL CloseHDFGroup(HDF_CurrentID)
          CALL CloseHDF(HDFFileID)
          CALL MondoLog(DEBUG_NONE, "MD:Symplectic", "putting dummy v_scale("//TRIM(IntToChar(iGEO))//") to hdf = "//TRIM(DblToChar(0D0)))
        ENDDO

      ENDIF

    ENDIF

    ! Update the Velocity if not the first step.
    !IF(iGEO > 1) THEN
    DO iCLONE=1,C%Geos%Clones
      DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
        IF(C%Geos%Clone(iCLONE)%CConstrain%I(iATS) == 0)THEN
          Mass      =  C%Geos%Clone(iCLONE)%AtMss%D(iATS)
          Vel(1:3)  =  C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS)
          Acc(1:3)  = -C%Geos%Clone(iCLONE)%Gradients%D(1:3,iATS)/Mass

          ! Velocity: v(t) = v(t-dT/2)+a(t)dT/2
          Vel(1:3) = Vel(1:3) + Acc(1:3)*dT*Symplectic_4th_Order_b(m_step)
          C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS) = Vel(1:3)
        ENDIF
      ENDDO
    ENDDO
    !ENDIF

    IF(m_step == 1) THEN
      ! Generate Output.
      CALL MondoLog(DEBUG_NONE, "MD:Symplectic", "writing coordinates to file")
      CALL OutputMD(C,iGEO)
    ENDIF

    ! Update the Positions.
    DO iCLONE=1,C%Geos%Clones
      DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
        IF(C%Geos%Clone(iCLONE)%CConstrain%I(iATS) == 0)THEN
          Mass     =  C%Geos%Clone(iCLONE)%AtMss%D(iATS)
          Pos(1:3) =  C%Geos%Clone(iCLONE)%Carts%D(1:3,iATS)
          Vel(1:3) =  C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS)
          Acc(1:3) = -C%Geos%Clone(iCLONE)%Gradients%D(1:3,iATS)/Mass

          ! Position: r(t+dT)= r(t)+v(t)*dT+a(t)*dT*dT/2
          ! Position: r(t+dT)= r(t)+v(t+Tt)*dT*a_coeff(t)
          Pos(1:3) = Pos(1:3) + Vel(1:3)*dT*Symplectic_4th_Order_a(m_step)

          ! Update
          C%Geos%Clone(iCLONE)%Carts%D(1:3,iATS) = Pos(1:3)
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE MDSymplectic_4th_Order_NVE

  !--------------------------------------------------------------
  ! Set an initial Temperature
  !--------------------------------------------------------------
  SUBROUTINE SetTempMaxBoltDist(C,Temp0)
    TYPE(Controls)        :: C
    REAL(DOUBLE)          :: Temp0
    REAL(DOUBLE)          :: Temp,Kin,Mass,TVel,VX,VY,VZ
    INTEGER               :: iCLONE,iATS,I,J,Jmax

    Jmax = 20
    DO iCLONE=1,C%Geos%Clones
      DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
        IF(C%Geos%Clone(iCLONE)%CConstrain%I(iATS) == 0)THEN
          Mass  = C%Geos%Clone(iCLONE)%AtMss%D(iATS)
          TVel  = SQRT(Three*Temp0*KelvinToHartrees/Mass)
          VX = Zero
          VY = Zero
          VZ = Zero
          DO J=1,Jmax
            VX = VX+TVel*Random((/-One,One/))
            VY = VY+TVel*Random((/-One,One/))
            VZ = VZ+TVel*Random((/-One,One/))
          ENDDO
          C%Geos%Clone(iCLONE)%Velocity%D(1,iATS) = VX
          C%Geos%Clone(iCLONE)%Velocity%D(2,iATS) = VY
          C%Geos%Clone(iCLONE)%Velocity%D(3,iATS) = VZ
        ENDIF
      ENDDO
      CALL ResetMomentum(C%Geos%Clone(iCLONE),Zero,Zero,Zero)
      CALL CalculateMDKin(C%Geos%Clone(iCLONE),Kin,Temp)
      CALL RescaleVelocity(C%Geos%Clone(iCLONE),Temp,Temp0)
      CALL CalculateMDKin(C%Geos%Clone(iCLONE),Kin,Temp)

      MDEkin%D(iCLONE)  = Kin
      MDTemp%D(iCLONE) = Temp
      MDTave%D(iCLONE) = Temp
      MDEpot%D(iCLONE) = C%Geos%Clone(iCLONE)%ETotal
      MDEtot%D(iCLONE) = MDEpot%D(iCLONE) + MDEkin%D(iCLONE)

    ENDDO
  END SUBROUTINE SetTempMaxBoltDist

  !--------------------------------------------------------------
  ! Reset the Velocities so that PX=PX0,PY=PY0,PZ=PZ0
  !--------------------------------------------------------------
  SUBROUTINE ResetMomentum(GM,PX0,PY0,PZ0)
    TYPE(CRDS)            :: GM
    REAL(DOUBLE)          :: Mass,PX,PY,PZ,PX0,PY0,PZ0
    INTEGER               :: iATS,nATOMS
    ! Calculate Linear Momentum
    PX = Zero
    PY = Zero
    PZ = Zero
    nATOMS = 0
    DO iATS=1,GM%NAtms
      IF(GM%CConstrain%I(iATS) == 0)THEN
        nATOMS = nATOMS + 1
        Mass = GM%AtMss%D(iATS)
        PX  = PX + Mass*GM%Velocity%D(1,iATS)
        PY  = PY + Mass*GM%Velocity%D(2,iATS)
        PZ  = PZ + Mass*GM%Velocity%D(3,iATS)
      ENDIF
    ENDDO
    ! Reset Linear Momentum
    DO iATS=1,GM%NAtms
      IF(GM%CConstrain%I(iATS) == 0)THEN
        Mass = GM%AtMss%D(iATS)
        GM%Velocity%D(1,iATS) = GM%Velocity%D(1,iATS)-(PX-PX0)/(nATOMS*Mass)
        GM%Velocity%D(2,iATS) = GM%Velocity%D(2,iATS)-(PY-PY0)/(nATOMS*Mass)
        GM%Velocity%D(3,iATS) = GM%Velocity%D(3,iATS)-(PZ-PZ0)/(nATOMS*Mass)
      ENDIF
    ENDDO

  END SUBROUTINE ResetMomentum

  !--------------------------------------------------------------
  ! Rescale the Velocities to Temp
  !--------------------------------------------------------------
  SUBROUTINE RescaleVelocity(GM,Temp,Temp0)
    TYPE(CRDS)            :: GM
    REAL(DOUBLE)          :: Temp,Temp0,Scale
    INTEGER               :: iATS

    Scale = SQRT(Temp0/Temp)
    DO iATS=1,GM%NAtms
      IF(GM%CConstrain%I(iATS) == 0)THEN
        GM%Velocity%D(1,iATS) = Scale*GM%Velocity%D(1,iATS)
        GM%Velocity%D(2,iATS) = Scale*GM%Velocity%D(2,iATS)
        GM%Velocity%D(3,iATS) = Scale*GM%Velocity%D(3,iATS)
      ENDIF
    ENDDO

  END SUBROUTINE RescaleVelocity

  SUBROUTINE BerendsenThermostat(GM, T, T0, delta_t, tau, v_scale_O)

    ! Fix v_scale. We use Berendsen's weak coupling to an external heat bath, J.
    ! Chem. Phys. 81, 3684 (1984). tau is the time constant of thermal coupling,
    !
    ! dT/dt = 1/tau * (T_0 - T)
    !
    ! where T_0 is the temperature of the heat bath. Perfect coupling means that
    !
    ! tau = Delta_t.

    TYPE(CRDS), INTENT(INOUT)           :: GM
    REAL(DOUBLE), INTENT(IN)            :: T, T0, delta_t, tau
    REAL(DOUBLE), OPTIONAL, INTENT(OUT) :: v_scale_O
    REAL(DOUBLE)                        :: v_scale
    INTEGER                             :: i

    IF(T > 1.0D-4) THEN
      v_scale = SQRT(1 + delta_t/tau * (T0/T - 1))
    ELSE
      v_scale = 1.0D0
    ENDIF

    IF(PRESENT(v_scale_O)) THEN
      v_scale_O = v_scale
    ENDIF

    CALL MondoLog(DEBUG_NONE, "Berendsen Thermostat", "T       = "//TRIM(DblToChar(T))//" K")
    CALL MondoLog(DEBUG_NONE, "Berendsen Thermostat", "T0      = "//TRIM(DblToChar(T0))//" K")
    CALL MondoLog(DEBUG_NONE, "Berendsen Thermostat", "delta_t = "//TRIM(DblToChar(delta_t*InternalTimeToFemtoseconds))//" fs")
    CALL MondoLog(DEBUG_NONE, "Berendsen Thermostat", "tau     = "//TRIM(DblToChar(tau*InternalTimeToFemtoseconds))//" fs")
    CALL MondoLog(DEBUG_NONE, "Berendsen Thermostat", "v_scale = "//TRIM(FltToChar(v_scale)))

    DO i = 1, GM%NAtms
      IF(GM%CConstrain%I(i) == 0) THEN
        GM%Velocity%D(1,i) = v_scale*GM%Velocity%D(1,i)
        GM%Velocity%D(2,i) = v_scale*GM%Velocity%D(2,i)
        GM%Velocity%D(3,i) = v_scale*GM%Velocity%D(3,i)
      ENDIF
    ENDDO

  END SUBROUTINE BerendsenThermostat

  SUBROUTINE BerendsenThermostatTotalEnergy(GM, Epot, Ekin, E0, delta_t, tau, v_scale_O)

    ! Fix v_scale. We use Berendsen's weak coupling to an external heat bath, J.
    ! Chem. Phys. 81, 3684 (1984). tau is the time constant of thermal coupling,
    !
    ! dT/dt = 1/tau * (T_0 - T)
    !
    ! where T_0 is the temperature of the heat bath. Perfect coupling means that
    !
    ! tau = Delta_t.

    TYPE(CRDS), INTENT(INOUT)           :: GM
    REAL(DOUBLE), INTENT(IN)            :: Epot, Ekin, E0, delta_t, tau
    REAL(DOUBLE), OPTIONAL, INTENT(OUT) :: v_scale_O
    REAL(DOUBLE)                        :: v_scale, E, Ekin_target
    REAL(DOUBLE), SAVE                  :: EIntegrated = 0.D0
    INTEGER                             :: i

    ! Calculate v_scale.
    E = Epot + Ekin
    Ekin_target = E0-Epot

    IF(Ekin_target > 0 .AND. Ekin > 1.0D-4) THEN
      EIntegrated = EIntegrated + (Ekin_target/Ekin - 1)
      v_scale = SQRT(1 + delta_t/tau * ((Ekin_target/Ekin - 1) + 1/10.D0*EIntegrated))
    ELSE
      v_scale = 1.0D0
    ENDIF

    IF(PRESENT(v_scale_O)) THEN
      v_scale_O = v_scale
    ENDIF

    CALL MondoLog(DEBUG_NONE, "Berendsen Thermostat", "E       = "//TRIM(DblToChar(E))//" eV")
    CALL MondoLog(DEBUG_NONE, "Berendsen Thermostat", "E0      = "//TRIM(DblToChar(E0))//" eV")
    CALL MondoLog(DEBUG_NONE, "Berendsen Thermostat", "delta_t = "//TRIM(DblToChar(delta_t*InternalTimeToFemtoseconds))//" fs")
    CALL MondoLog(DEBUG_NONE, "Berendsen Thermostat", "tau     = "//TRIM(DblToChar(tau*InternalTimeToFemtoseconds))//" fs")
    CALL MondoLog(DEBUG_NONE, "Berendsen Thermostat", "v_scale = "//TRIM(FltToChar(v_scale)))

    DO i = 1, GM%NAtms
      IF(GM%CConstrain%I(i) == 0) THEN
        GM%Velocity%D(1,i) = v_scale*GM%Velocity%D(1,i)
        GM%Velocity%D(2,i) = v_scale*GM%Velocity%D(2,i)
        GM%Velocity%D(3,i) = v_scale*GM%Velocity%D(3,i)
      ENDIF
    ENDDO

  END SUBROUTINE BerendsenThermostatTotalEnergy

  !--------------------------------------------------------------
  ! Calculate the Temperature
  !--------------------------------------------------------------
  SUBROUTINE CalculateMDKin(GM,Kin,Temp)
    TYPE(CRDS)            :: GM
    REAL(DOUBLE)          :: Kin,Temp,Mass
    INTEGER               :: iATS

    Kin = Zero
    nATOMS = 0
    DO iATS=1,GM%NAtms
      IF(GM%CConstrain%I(iATS) == 0)THEN
        nATOMS = nATOMS + 1
        Mass = GM%AtMss%D(iATS)
        Kin  = Kin + Half*GM%AtMss%D(iATS) &
             *(GM%Velocity%D(1,iATS)**2 &
              +GM%Velocity%D(2,iATS)**2 &
              +GM%Velocity%D(3,iATS)**2)
      ENDIF
    ENDDO
    Temp = (Two/Three)*(Kin/DBLE(nATOMS))*HartreesToKelvin

  END SUBROUTINE CalculateMDKin

  !--------------------------------------------------------------
  ! The Verlet Algorithmn : NVT
  !--------------------------------------------------------------
  SUBROUTINE MDVerlet_NVT(C,iGEO)
    TYPE(Controls)            :: C
    INTEGER                   :: iGEO,iCLONE,iATS
    REAL(DOUBLE)              :: Mass,dT
    REAL(DOUBLE),DIMENSION(3) :: Pos,Vel,Acc
    !--------------------------------------------------------------
  END SUBROUTINE MDVerlet_NVT

  !--------------------------------------------------------------
  ! The Verlet Algorithmn : NPH
  !--------------------------------------------------------------
  SUBROUTINE MDVerlet_NPH(C,iGEO)
    TYPE(Controls)            :: C
    INTEGER                   :: iGEO,iCLONE,iATS
    REAL(DOUBLE)              :: Mass,dT
    REAL(DOUBLE),DIMENSION(3) :: Pos,Vel,Acc
    !--------------------------------------------------------------
  END SUBROUTINE MDVerlet_NPH

  !--------------------------------------------------------------
  ! The Verlet Algorithmn : NPT
  !--------------------------------------------------------------
  SUBROUTINE MDVerlet_NPT(C,iGEO)
    TYPE(Controls)            :: C
    INTEGER                   :: iGEO,iCLONE,iATS
    REAL(DOUBLE)              :: Mass,dT
    REAL(DOUBLE),DIMENSION(3) :: Pos,Vel,Acc
    !--------------------------------------------------------------
  END SUBROUTINE MDVerlet_NPT

  !--------------------------------------------------------------
  ! Output
  !--------------------------------------------------------------
  SUBROUTINE OutputMD(C,iGEO)
    TYPE(Controls)                 :: C
    INTEGER                        :: iGEO,iCLONE,iATS,I,J,K
    REAL(DOUBLE)                   :: Pressure,DMax,Entropy,Ftot
    REAL(DOUBLE),DIMENSION(3,3)    :: Latt,InvLatt,LattFrc
    REAL(DOUBLE),DIMENSION(3)      :: Pos,Vel
    LOGICAL, SAVE                  :: FirstTime = .TRUE.
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: Remark

    DO iCLONE=1,C%Geos%Clones
      CALL MondoLog(DEBUG_MAXIMUM, "MD:OutputMD", "writing coordinates, t = " &
        //TRIM(DblToMedmChar(MDTime%D(iCLONE)*InternalTimeToFemtoseconds))//" fs")
      Remark = 't = '//TRIM(DblToMedmChar(MDTime%D(iCLONE)*InternalTimeToFemtoseconds))//" fs, "// &
               "Ekin = "//TRIM(DblToChar(MDEkin%D(iCLONE)*au2eV))//" eV, "// &
               "Epot = "//TRIM(DblToChar(MDEpot%D(iCLONE)*au2eV))//" eV, "// &
               "Etot = "//TRIM(DblToChar(MDEtot%D(iCLONE)*au2eV))//" eV, "// &
               "T = "//TRIM(DblToMedmChar(MDTemp%D(iCLONE)))//" K"
      HDFFileID=OpenHDF(C%Nams%HFile)
      HDF_CurrentID = OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
      CALL Get(Entropy, "Entropy")
      CALL CloseHDFGroup(HDF_CurrentID)
      CALL CloseHDF(HDFFileID)
      IF(Entropy > Zero) THEN
        Ftot = MDEtot%D(iCLONE)-Entropy
        Remark = TRIM(Remark)//", Entropy = "//TRIM(DblToChar(Entropy*au2eV))// " eV, Ftot = " &
                 //TRIM(DblToChar(Ftot*au2eV))//" eV"
      ENDIF
      CALL PPrint(C%Geos%Clone(iCLONE), FileName_O=C%Nams%GFile, Unit_O=Geo, &
                  PrintGeom_O=C%Opts%GeomPrint, Clone_O=iCLONE, Remark_O=Remark, Gradients_O='Velocities')
    ENDDO

  END SUBROUTINE OutputMD

  !--------------------------------------------------------------
  ! Move the last P or D marices to other names
  !--------------------------------------------------------------
  SUBROUTINE CopyMatrices(C,Tilde)
    TYPE(Controls)                 :: C
    INTEGER                        :: iCLONE,iSCF,iBAS,iGEO,oldMyClone
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: chSCF,chBAS,chGEO,chCLONE
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: PoldFile,PnewFile
    LOGICAL                        :: Tilde

    iSCF = C%Stat%Current%I(1)
    iBAS = C%Stat%Current%I(2)
    iGEO = C%Stat%Current%I(3)

    CALL MondoLog(DEBUG_MAXIMUM, "CopyMatrices", "Current state = "//TRIM(IntVectToChar(C%Stat%Current)))

    IF(Tilde) THEN
      SELECT CASE(C%Dyns%MDGuess)
      CASE('DMLinear', &
           'DMTRBO', &
           'DMTRBO_Damp_dt3', &
           'DMTRBO_Damp_dt5', &
           'DMTRBO_Damp_dt7', &
           'DMTRBO_Damp_dt9', &
           'DMTRBO_Damp_dt11', &
           'DMSymplectic', &
           'DMProj0', &
           'DMProj1', &
           'DMProj2', &
           'DMProj3', &
           'DMProj4')

        oldMyClone = myClone
        DO iCLONE=1,C%Geos%Clones
          myClone = iCLONE
          CALL FileCopy(TrixFile("OrthoD",  Stats_O = C%Stat%Current%I, Offset_O = 1), &
                        TrixFile("DOPsave", Stats_O = C%Stat%Current%I))
        ENDDO
        myClone = oldMyClone
      CASE('FMVerlet0','FMVerlet1')
        DO iCLONE=1,C%Geos%Clones
          chSCF   = IntToChar(iSCF)
          chBAS   = IntToChar(iBAS)
          chGEO   = IntToChar(iGEO)
          chCLONE = IntToChar(iCLONE)

          PoldFile = TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)  &
               //'_Geom#'//TRIM(chGEO)  &
               //'_Base#'//TRIM(chBAS)  &
               //'_Cycl#'//TRIM(chSCF)  &
               //'_Clone#'//TRIM(chCLONE)//'.OrthoF'
          PnewFile = TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)  &
               //'_G#'//TRIM(chGEO)     &
               //'_C#'//TRIM(chCLONE)//'.FOPsave'
          CALL FileCopy(PoldFile, PnewFile)
        ENDDO
      CASE DEFAULT
        CALL MondoLog(DEBUG_MAXIMUM, "CopyMatrices", "in default case with tilde")
      ENDSELECT
    ELSE
      SELECT CASE(C%Dyns%MDGuess)
      CASE('DMLinear', &
           'DMTRBO', &
           'DMTRBO_Damp_dt3', &
           'DMTRBO_Damp_dt5', &
           'DMTRBO_Damp_dt7', &
           'DMTRBO_Damp_dt9', &
           'DMTRBO_Damp_dt11', &
           'DMSymplectic', &
           'DMProj0', &
           'DMProj1', &
           'DMProj2', &
           'DMProj3', &
           'DMProj4')

        oldMyClone = myClone
        DO iCLONE=1,C%Geos%Clones
          myClone = iCLONE
          CALL FileCopy(TrixFile("OrthoD", Stats_O = C%Stat%Current%I, Offset_O = 1), &
                        TrixFile("DOsave", Stats_O = C%Stat%Current%I))
        ENDDO
        myClone = oldMyClone
      CASE('FMVerlet0','FMVerlet1')
        DO iCLONE=1,C%Geos%Clones
          chSCF   = IntToChar(iSCF)
          chBAS   = IntToChar(iBAS)
          chGEO   = IntToChar(iGEO)
          chCLONE = IntToChar(iCLONE)
          PoldFile = TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)  &
               //'_Geom#'//TRIM(chGEO)  &
               //'_Base#'//TRIM(chBAS)  &
               //'_Cycl#'//TRIM(chSCF)  &
               //'_Clone#'//TRIM(chCLONE)//'.OrthoF'
          PnewFile = TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)  &
               //'_G#'//TRIM(chGEO)     &
               //'_C#'//TRIM(chCLONE)//'.FOsave'
          CALL FileCopy(PoldFile, PnewFile)
        ENDDO
      CASE DEFAULT
        CALL MondoLog(DEBUG_MAXIMUM, "CopyMatrices", "in default case without tilde")
      END SELECT
    ENDIF

  END SUBROUTINE CopyMatrices

  !--------------------------------------------------------------
  ! Copy the Restart Marices
  !--------------------------------------------------------------
  SUBROUTINE CopyRestart(C,iGEOBegin,MinMDGeo)
    TYPE(Controls)                 :: C
    INTEGER                        :: iCLONE,iGEO,iGEOBegin,MinMDGeo
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: chGEO,chCLONE
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: PoldFile,PnewFile
    INTEGER                        :: copyStatus
    INTEGER                        :: oldMyClone

    CALL MondoLog(DEBUG_NONE, "CopyRestart", "copying density matrix files between "//TRIM(IntToChar(iGEOBegin-MinMDGeo))//" and "//TRIM(IntToChar(iGEOBegin)))
    oldMyClone = MyClone
    DO iCLONE=1,C%Geos%Clones
      MyClone = iCLONE
      DO iGEO=iGEOBegin-MinMDGeo,iGEOBegin
        chCLONE = IntToChar(iCLONE)
        chGEO   = IntToChar(iGEO)
        SELECT CASE(C%Dyns%MDGuess)
        CASE('DMLinear', &
             'DMTRBO', &
             "DMTRBO_Damp_dt3", &
             "DMTRBO_Damp_dt5", &
             "DMTRBO_Damp_dt7", &
             "DMTRBO_Damp_dt9", &
             "DMTRBO_Damp_dt11", &
             'DMSymplectic')

          CALL FileCopy(TrixFile("DOsave", Name_O = C%Nams%RFile(1:LEN(TRIM(C%Nams%RFile))-4), &
                                           PWD_O = "", &
                                           Stats_O = (/ 0, 1, iGeo /)), &
                        TrixFile("DOsave", Stats_O = (/ 0, 1, iGeo /)))

          CALL FileCopy(TrixFile("DOPsave", Name_O = C%Nams%RFile(1:LEN(TRIM(C%Nams%RFile))-4), &
                                            PWD_O = "", &
                                            Stats_O = (/ 0, 1, iGeo /)), &
                        TrixFile("DOPsave", Stats_O = (/ 0, 1, iGeo /)))

        CASE('FMVerlet0','FMVerlet1')

          PoldFile = TRIM(C%Nams%RFile(1:LEN(TRIM(C%Nams%RFile))-4))//'_G#'//TRIM(chGEO)//'_C#'//TRIM(chCLONE)//'.FOsave'
          PnewFile = TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)//'_G#'//TRIM(chGEO)//'_C#'//TRIM(chCLONE)//'.FOsave'
          CALL MondoLog(DEBUG_NONE, "CopyRestart", TRIM(PoldFile)//" --> "//TRIM(PnewFile))
          CALL FileCopy(PoldFile, PnewFile)

          PoldFile = TRIM(C%Nams%RFile(1:LEN(TRIM(C%Nams%RFile))-4))//'_G#'//TRIM(chGEO)//'_C#'//TRIM(chCLONE)//'.FOPsave'
          PnewFile = TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)//'_G#'//TRIM(chGEO)//'_C#'//TRIM(chCLONE)//'.FOPsave'
          CALL MondoLog(DEBUG_NONE, "CopyRestart", TRIM(PoldFile)//" --> "//TRIM(PnewFile))
          CALL FileCopy(PoldFile, PnewFile)

        END SELECT
      ENDDO
    ENDDO
    MyClone = oldMyClone

  END SUBROUTINE CopyRestart

  SUBROUTINE OptimizeLattice(C)
    TYPE(Controls)  :: C
    INTEGER         :: iBAS,iGEO,iCLONE,I,J

    iBAS   = C%Sets%NBSets
    iCLONE = 1
    DO iGEO = 1,50
      ! Determine a New Lattice
      CALL OpenASCII("LatticeOpt.dat",77)
      WRITE(77,*) 'Geom   = ',iGEO
      WRITE(77,*) "Energy = ",C%Geos%Clone(iCLONE)%ETotal
      WRITE(77,*) 'Lattice'
      DO I=1,3
        WRITE(77,*) (C%Geos%Clone(iCLONE)%PBC%BoxShape%D(I,J)*AUToAngstroms,J=1,3)
      ENDDO
      WRITE(77,*) 'LatFrc'
      DO I=1,3
        WRITE(77,*) (C%Geos%Clone(iCLONE)%PBC%LatFrc%D(I,J)*AUToAngstroms,J=1,3)
      ENDDO
      CLOSE(77)

      WRITE(*,*) 'Geom   = ',iGEO
      WRITE(*,*) "Energy = ",C%Geos%Clone(iCLONE)%ETotal
      WRITE(*,*) 'Lattice'
      DO I=1,3
        WRITE(*,*) (C%Geos%Clone(iCLONE)%PBC%BoxShape%D(I,J)*AUToAngstroms,J=1,3)
      ENDDO
      WRITE(*,*) 'LatFrc'
      DO I=1,3
        WRITE(*,*) (C%Geos%Clone(iCLONE)%PBC%LatFrc%D(I,J)*AUToAngstroms,J=1,3)
      ENDDO

      DO I=1,3
        DO J=I,3
          C%Geos%Clone(iCLONE)%PBC%BoxShape%D(I,J) = C%Geos%Clone(iCLONE)%PBC%BoxShape%D(I,J)    &
               - 5.0D0*C%Geos%Clone(iCLONE)%PBC%LatFrc%D(I,J)
        ENDDO
      ENDDO

      ! Archive Geometry for next step
      CALL MkGeomPeriodic(C%Geos%Clone(iCLONE))
      C%Geos%Clone(1)%Carts%D = C%Geos%Clone(iCLONE)%Carts%D
      CALL GeomArchive(iBAS,iGEO+1,C%Nams,C%Opts,C%Sets,C%Geos)

      ! Evaluate energies at the new geometry
      CALL SCF(iBAS,iGEO+1,C)

      ! Calculate Force
      CALL Force(iBAS,iGEO+1,C%Nams,C%Opts,C%Stat,C%Geos,C%Sets,C%MPIs)
    ENDDO

  END SUBROUTINE OptimizeLattice

END MODULE MDynamics
