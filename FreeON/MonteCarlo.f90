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
MODULE MonteCarlo
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
  USE MDynamics

  IMPLICIT NONE

  TYPE(INT_VECT)      :: MDIter
  TYPE(DBL_VECT)      :: MCEtot0,MCTemp0
  TYPE(DBL_RNK2)      :: Carts
  TYPE(DBL_RNK3)      :: MCCarts0

  CONTAINS
!--------------------------------------------------------------
! Main MD Subroutine
!--------------------------------------------------------------
  SUBROUTINE  HybridMC(C)
    TYPE(Controls)  :: C
    INTEGER         :: iMC,iMD,iSTART,iSTATUS
    INTEGER         :: iSCF,iBAS,iGEO,iCLONE,iREMOVE,I
    REAL(DOUBLE)    :: DelEtot,ExpFac,RanFac

!--------------------------------------------------------------
!   Intitialize
    C%Stat%Previous%I=(/0,1,1/)
    iGEO      = 1
!
    WRITE(*,*) "Hybrid Monte-Carlo"
!
    CALL New(MDIter ,C%Geos%Clones)
    CALL New(MDTime ,C%Geos%Clones)
    CALL New(MDEkin ,C%Geos%Clones)
    CALL New(MDEpot ,C%Geos%Clones)
    CALL New(MDEtot ,C%Geos%Clones)
    CALL New(MDTemp ,C%Geos%Clones)
    CALL New(MDTave ,C%Geos%Clones)
    CALL New(MDLinP ,(/3,C%Geos%Clones/))
!
    CALL New(MCEtot0,C%Geos%Clones)
    CALL New(MCTemp0,C%Geos%Clones)
    CALL New(MCCarts0,(/3,C%Geos%Clone(1)%NAtms,C%Geos%Clones/))
    CALL New(Carts   ,(/3,C%Geos%Clone(1)%NAtms/))
!
    MDIter%I  = 1
    MDTime%D  = Zero
    MDEkin%D  = Zero
    MDEpot%D  = Zero
    MDEtot%D  = Zero
    MDTemp%D  = Zero
    MDTave%D  = Zero
    MDLinP%D  = Zero
!
    MCEtot0%D = Zero
    MCTemp0%D = Zero
    MCCarts0%D= Zero
!
!   Initial Guess
    IF(C%Opts%Guess==GUESS_EQ_RESTART) THEN
       WRITE(*,*) 'Geuss==RESTART'
!      Init the Time
       HDFFileID=OpenHDF(C%Nams%RFile)
       DO iCLONE=1,C%Geos%Clones
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
          CALL Get(MDIter%I(iCLONE)  ,"MDIter")
          CALL Get(C%Dyns%MDGeuss    ,"MDGeuss")
          CALL Get(MCEtot0%D(iCLONE) ,"MCEtot0")
          CALL Get(MCTemp0%D(iCLONE) ,"MCTemp0")
          CALL Get(Carts,"MCCarts0")
          MCCarts0%D(:,:,iCLONE)=Carts%D(:,:)
          CALL CloseHDFGroup(HDF_CurrentID)
       ENDDO
       CALL CloseHDF(HDFFileID)
!      Save to Current HDF
       HDFFileID=OpenHDF(C%Nams%HFile)
       DO iCLONE=1,C%Geos%Clones
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
          CALL Put(MDIter%I(iCLONE),"MDIter")
          CALL Put(MDTime%D(iCLONE),"MDTime")
          CALL Put(.TRUE.,"DoingMD")
          CALL Put(.TRUE.,"DoingHybridMC")
          CALL Put(C%Dyns%MDGeuss   ,"MDGeuss")
          CALL Put(MCEtot0%D(iCLONE),"MCEtot0")
          CALL Put(MCTemp0%D(iCLONE),"MCTemp0")
          Carts%D(:,:) = MCCarts0%D(:,:,iCLONE)
          CALL Put(Carts,"MCCarts0")
          CALL CloseHDFGroup(HDF_CurrentID)
       ENDDO
       CALL CloseHDF(HDFFileID)
!      Do The initial SCF
       iBAS=C%Sets%NBSets
       CALL GeomArchive(iBAS,iGEO,C%Nams,C%Opts,C%Sets,C%Geos)
       CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
       CALL SCF(iBAS,iGEO,C)
    ELSE
       WRITE(*,*) 'Geuss==SUPERPOS'
!      Do The initial SCF: Build the SCF
       DO iBAS=1,C%Sets%NBSets
          CALL GeomArchive(iBAS,iGEO,C%Nams,C%Opts,C%Sets,C%Geos)
          CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
          CALL SCF(iBAS,iGEO,C)
       ENDDO
!      Set up initial Temperature for the MD run for the MC Hybrid
       CALL SetTempMaxBoltDist(C,C%Dyns%MCTemp)
!      Store stuff for the MC, this step (Temp, Energy, Carts)
       DO iCLONE=1,C%Geos%Clones
          MCTemp0%D(iCLONE)     = C%Dyns%MCTemp
          MCEtot0%D(iCLONE)     = MDEtot%D(iCLONE)
          MCCarts0%D(:,:,iCLONE)= C%Geos%Clone(iCLONE)%Carts%D(:,:)
       ENDDO
!      Save to Current HDF
       HDFFileID=OpenHDF(C%Nams%HFile)
       DO iCLONE=1,C%Geos%Clones
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
          CALL Put(MDIter%I(iCLONE),"MDIter")
          CALL Put(MDTime%D(iCLONE),"MDTime")
          CALL Put(.TRUE.,"DoingMD")
          CALL Put(.TRUE.,"DoingHybridMC")
          CALL Put(C%Dyns%MDGeuss   ,"MDGeuss")
          CALL Put(MCEtot0%D(iCLONE),"MCEtot0")
          CALL Put(MCTemp0%D(iCLONE),"MCTemp0")
          Carts%D(:,:) = MCCarts0%D(:,:,iCLONE)
          CALL Put(Carts,"MCCarts0")
          CALL CloseHDFGroup(HDF_CurrentID)
       ENDDO
       CALL CloseHDF(HDFFileID)
    ENDIF
!   Determine iREMOVE
    SELECT CASE(C%Dyns%MDGeuss)
    CASE ('DMVerlet')
       iREMOVE = 4
    CASE ('FMVerlet')
       iREMOVE = 4
    CASE ('DMProj0')
       iREMOVE = 1
    CASE ('DMProj1')
       iREMOVE = 2
    CASE ('DMProj2')
       iREMOVE = 3
    CASE ('DMProj3')
       iREMOVE = 4
    CASE ('DMProj4')
       iREMOVE = 5
    CASE ('DMDGeuss')
       iREMOVE = 0
    END SELECT
!   Initialize MD and MC
    iBAS=C%Sets%NBSets
!    CALL RenameDensityMatrix(C,C%Stat%Current%I(1),C%Stat%Current%I(2),C%Stat%Current%I(3))
!    CALL CopyDensityMatrix(C,C%Stat%Current%I(3),C%Stat%Current%I(3)+C%Dyns%MDMaxSteps)
    CALL OutputMD(C,0)
    CALL OutputMC(C,0,1)
!   Do MC
    iGEO = 0
    Do iMC = 1,C%Dyns%MCMaxSteps
!      Do MD
       IF(iMC==1) THEN
          iSTART = 2
       ELSE
          iSTART = 1
       ENDIF
       DO iMD = iSTART,C%Dyns%MDMaxSteps-MDIter%I(1)+1
          iGEO = iGEO+1
!         Calculate the Force
          CALL Force(iBAS,iGEO,C%Nams,C%Opts,C%Stat,C%Geos,C%Sets,C%MPIs)
          IF(C%Dyns%MDAlgorithm == MD_AL_VERLET) THEN
             CALL MDVerlet_NVE(C,iGEO)
          ELSE
             CALL Halt('Other MD algorithms are Not Implimented for Hybrid MC')
          ENDIF
!         Generate Output
          CALL OutputMD(C,iMD)
!         Refresh the Time
          HDFFileID=OpenHDF(C%Nams%HFile)
          DO iCLONE=1,C%Geos%Clones
             MDIter%I(iCLONE) = iMD
             MDTime%D(iCLONE) = MDTime%D(iCLONE)+C%Dyns%DTime
             HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
             CALL Put(MDIter%I(iCLONE),"MDIter")
             CALL Put(MDTime%D(iCLONE),"MDTime")
             CALL CloseHDFGroup(HDF_CurrentID)
          ENDDO
          CALL CloseHDF(HDFFileID)
!         Archive Geometry for next step
          CALL GeomArchive(iBAS,iGEO+1,C%Nams,C%Opts,C%Sets,C%Geos)
!         Evaluate energies at the new geometry
          CALL SCF(iBAS,iGEO+1,C)
!         Store the Last P matrix
!          CALL RenameDensityMatrix(C,C%Stat%Current%I(1),C%Stat%Current%I(2),C%Stat%Current%I(3))
          IF(C%Stat%Current%I(3)-iREMOVE > 0) THEN
!             CALL RemoveDensityMatrix(C,C%Stat%Current%I(3)-iREMOVE)
          ENDIF
          IF(iMD < iREMOVE+1) THEN
!             CALL CopyDensityMatrix(C,C%Stat%Current%I(3),C%Stat%Current%I(3)+C%Dyns%MDMaxSteps)
          ENDIF
!         Remove old Stuff from Scratch
          CALL CleanScratch(C,iGEO)
       ENDDO
       MDIter%I  = 1
       MDTime%D  = Zero
!      Determine what to do Next
       DO iCLONE=1,C%Geos%Clones
          DelEtot =  MCEtot0%D(iCLONE)-MDEtot%D(iCLONE)
          ExpFac = EXP(-DelEtot/(KelvinToHartrees*C%Dyns%MCTemp))
          RanFac = Random((/Zero,One/))
          IF(RanFac < ExpFac) THEN
!            Accept Move, redue Max-Botz Temp Distribution
             CALL SetTempMaxBoltDist(C,C%Dyns%MCTemp)
!
             MCEtot0%D(iCLONE) = MDEtot%D(iCLONE)
             MCCarts0%D(:,:,iCLONE)= C%Geos%Clone(iCLONE)%Carts%D(:,:)
             DO I=1,iREMOVE
!                CALL RemoveDensityMatrix(C,C%Stat%Current%I(3)+I)
             ENDDO
             iSTATUS=1
          ELSE
!            Reject Move, redue Max-Botz Temp Distribution
             CALL SetTempMaxBoltDist(C,C%Dyns%MCTemp)
             MCEtot0%D(iCLONE) = MDEtot%D(iCLONE)
             C%Geos%Clone(iCLONE)%Carts%D(:,:)=MCCarts0%D(:,:,iCLONE)
!            Place old density matrices where they are needed
             DO I=1,iREMOVE
!                CALL   CopyDensityMatrix(C,C%Stat%Current%I(3)+I,C%Stat%Current%I(3)+I-iREMOVE)
!                CALL RemoveDensityMatrix(C,C%Stat%Current%I(3)+I)
             ENDDO
             iSTATUS=0
          ENDIF
       ENDDO
!      Save to Current HDF
       HDFFileID=OpenHDF(C%Nams%HFile)
       DO iCLONE=1,C%Geos%Clones
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
          CALL Put(MCEtot0%D(iCLONE),"MCEtot0")
          CALL Put(MCTemp0%D(iCLONE),"MCTemp0")
          Carts%D(:,:) = MCCarts0%D(:,:,iCLONE)
          CALL Put(Carts,"MCCarts0")
          CALL CloseHDFGroup(HDF_CurrentID)
       ENDDO
       CALL CloseHDF(HDFFileID)
!      Generate Output
       CALL OutputMC(C,iMC,iSTATUS)
    ENDDO
!
 END SUBROUTINE HybridMC
!--------------------------------------------------------------
! Output
!--------------------------------------------------------------
 SUBROUTINE OutputMC(C,iMC,iSTATUS)
   TYPE(Controls)                 :: C
   INTEGER                        :: iMC,iSTATUS,iCLONE,iATS,I,J,K
   CHARACTER(LEN=DEFAULT_CHR_LEN) :: File,Line
!
    DO iCLONE=1,C%Geos%Clones
!      Open File
       File = TRIM(C%Nams%SCF_NAME)//"_Clone#"//TRIM(IntToChar(iCLONE))//".MCO"
       CALL OpenASCII(File,Out)
!      Add Header
       IF(iMC==0) THEN
          Line = "##########################################################################"
          WRITE(Out,97) Line
          Line = "# MC Clone No. = "//TRIM(IntToChar(iCLONE))
          WRITE(Out,97) Line
          Line = "# MC Atoms No. = "//TRIM(IntToChar(C%Geos%Clone(iCLONE)%NAtms))
          WRITE(Out,97) Line
          Line = "# MD Algorithm = Verlet"
          WRITE(Out,97) Line
          Line = "# MDGeuss      = "//TRIM(C%Dyns%MDGeuss)
          WRITE(Out,97) Line
          Line = "# Minium SCF   = "//TRIM(IntToChar(C%Opts%MinSCF))
          WRITE(Out,97) Line
          Line = "# Maxium SCF   = "//TRIM(IntToChar(C%Opts%MaxSCF))
          WRITE(Out,97) Line
          Line = "# MD Time Step = "//TRIM(DblToMedmChar(C%Dyns%DTime))//" au"
          WRITE(Out,97) Line
          Line = "# MD Time Step = "//TRIM(DblToMedmChar(C%Dyns%DTime*InternalTimeToSeconds))//" sec"
          WRITE(Out,97) Line
          Line = "# MD Max Step  = "//TRIM(IntToChar(C%Dyns%MDMaxSteps))
          WRITE(Out,97) Line
          Line = "# MC Max Step  = "//TRIM(IntToChar(C%Dyns%MCMaxSteps))
          WRITE(Out,97) Line
          Line = "# MC Temper    = "//TRIM(DblToMedmChar(C%Dyns%MCTemp))
          WRITE(Out,97) Line
          Line = "##########################################################################"
          WRITE(Out,97) Line
          Line = " "
          WRITE(Out,97) Line
       ENDIF
       Line = "-------------------------------------------MC-"//TRIM(IntToChar(iMC))// &
              "-MC-------------------------------------------"
       WRITE(Out,97) Line
       WRITE(Out,98) "MC Potential   = ",MDEpot%D(iCLONE)
       IF(iSTATUS==0) THEN
          WRITE(Out,*) "Move Rejected"
       ELSE
          WRITE(Out,*) "Move Accepted"
       ENDIF
       WRITE(Out,85)
       DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
          WRITE(Out,99) TRIM(C%Geos%Clone(iCLONE)%AtNam%C(iATS)),  &
                        C%Geos%Clone(iCLONE)%Carts%D(1:3,iATS)
       ENDDO
       CLOSE(Out)
    ENDDO
!
80  FORMAT(a18,F16.10)
81  FORMAT(a18,F16.10,1x,F16.10)
82  FORMAT(a18,F16.10,1x,F16.10,1x,F16.10)
85  FORMAT("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
!
96  FORMAT(a18,3(F16.10,1x))
97  FORMAT(a128)
98  FORMAT(a18,F16.10)
99  FORMAT(a3,6(1x,F18.12))
!
 END SUBROUTINE OutputMC

!
!
!
END MODULE MonteCarlo
