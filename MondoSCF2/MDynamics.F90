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
  IMPLICIT NONE
  TYPE(DBL_VECT)      :: MDTime,MDKin,MDEpot,MDEtot,MDTemp,MDTave
  TYPE(DBL_RNK2)      :: MDLinP
  CONTAINS
!--------------------------------------------------------------
! Main MD Subroutine
!--------------------------------------------------------------
 SUBROUTINE  MD(C)  
    TYPE(Controls)  :: C
    INTEGER         :: iSCF,iBAS,iGEO,iCLONE,iREMOVE,I,DMPOrder
    REAL(DOUBLE)    :: Temp 
!--------------------------------------------------------------
!   Do Molecular Dynamics:Loop over Time Steps
!   Intitialize
    C%Stat%Previous%I=(/0,1,1/)
    iGEO      = 1 
    DMPOrder  = C%Opts%DMPOrder
    iREMOVE   = DMPOrder+1
    CALL New(MDTime ,C%Geos%Clones)
    CALL New(MDKin  ,C%Geos%Clones)
    CALL New(MDEpot ,C%Geos%Clones)
    CALL New(MDETot ,C%Geos%Clones)
    CALL New(MDTemp ,C%Geos%Clones)
    CALL New(MDTave ,C%Geos%Clones)
    CALL New(MDLinP ,(/3,C%Geos%Clones/))
!
    MDTime%D = Zero
    MDKin%D  = Zero
    MDEpot%D = Zero
    MDETot%D = Zero
    MDTemp%D = Zero
    MDTave%D = Zero
    MDLinP%D = Zero
!
!   Initial Guess     
    IF(C%Opts%Guess==GUESS_EQ_RESTART) THEN
!      Init the Time
       HDFFileID=OpenHDF(C%Nams%RFile)
       DO iCLONE=1,C%Geos%Clones
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
          CALL Get(MDTime%D(iCLONE),"MDTime")
          CALL Get(C%Opts%DMPOrder,"DMPOrder")
          CALL CloseHDFGroup(HDF_CurrentID)
       ENDDO   
       CALL CloseHDF(HDFFileID)
!      Save to Current HDF
       HDFFileID=OpenHDF(C%Nams%HFile)
       DO iCLONE=1,C%Geos%Clones
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
          CALL Put(MDTime%D(iCLONE),"MDTime")
          CALL Put(.TRUE.,"DoingMD")
          CALL Put(C%Opts%DMPOrder,"DMPOrder")
          CALL CloseHDFGroup(HDF_CurrentID)
       ENDDO
       CALL CloseHDF(HDFFileID)
!      Do The initial SCF
       iBAS=C%Sets%NBSets
       CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos) 
       CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
       CALL SCF(iBAS,iGEO,C)
    ELSE
!      Give an intial Maxwell Boltzman Temp
       IF(C%Dyns%Initial_Temp) THEN
          CALL SetTempMaxBoltDist(C,C%Dyns%TempInit)
       ENDIF
!      Init the Time
       MDTime%D(:) = Zero
       HDFFileID=OpenHDF(C%Nams%HFile)
       DO iCLONE=1,C%Geos%Clones
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
          CALL Put(MDTime%D(iCLONE),"MDTime")
          CALL Put(.TRUE.,'DoingMD')
          CALL Put( C%Opts%DMPOrder,'DMPOrder')
          CALL CloseHDFGroup(HDF_CurrentID)
       ENDDO
       CALL CloseHDF(HDFFileID)
!      Build the SCF 
       DO iBAS=1,C%Sets%NBSets
          CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos) 
          CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
          CALL SCF(iBAS,iGEO,C)
       ENDDO
    ENDIF
!   Initialize MD
    iBAS=C%Sets%NBSets
    CALL RenameDensityMatrix(C,C%Stat%Current%I(1),C%Stat%Current%I(2),C%Stat%Current%I(3))
    CALL OutputMD(C,0)
!   Do MD
    DO iGEO = 1,C%Dyns%MDMaxSteps
!      Calculate the Force
       CALL Force(iBAS,iGEO,C%Nams,C%Opts,C%Stat,C%Geos,C%Sets,C%MPIs)
!      Move the Atoms, apply the theromstats and output
       IF(C%Dyns%MDAlgorithm==VERLET_MD_AL) THEN
          IF(         C%Dyns%Const_Temp  .AND. .NOT.C%Dyns%Const_Press) THEN
             CALL Halt('Constant Temperature MD Not Implimented')
             CALL MDVerlet_NVT(C,iGEO)
          ELSEIF(     C%Dyns%Const_Press .AND. .NOT.C%Dyns%Const_Temp ) THEN
             CALL Halt('Constant Presure MD Not Implimented')
             CALL MDVerlet_NPH(C,iGEO)
          ELSEIF(     C%Dyns%Const_Press .AND.      C%Dyns%Const_Temp ) THEN
             CALL Halt('Constant Presure/Constant Temperature MD Not Implimented')
             CALL MDVerlet_NPT(C,iGEO)
          ELSEIF(.NOT.C%Dyns%Const_Press .AND. .NOT.C%Dyns%Const_Temp ) THEN
             CALL MDVerlet_NVE(C,iGEO)
          ENDIF
       ELSEIF(C%Dyns%MDAlgorithm==GEAR_MD_AL) THEN
          CALL Halt('Gear MD Not Implimented')
       ENDIF
!      If Parallel Rep, switch Geometries
       IF(C%Dyns%Parallel_Rep) THEN
          CALL Halt('Parallel Replicate MD Not Implimented')
       ENDIF
!      Generate Output
       CALL OutputMD(C,iGEO)
!      Theomstates
       IF(C%Dyns%Temp_Scaling) THEN
          IF(MOD(iGEO,C%Dyns%RescaleInt)==0) THEN
             WRITE(*,*) 'Rescaling Temperature' 
             DO iCLONE=1,C%Geos%Clones
                IF(.TRUE.) THEN
                   CALL RescaleVelocity(C%Geos%Clone(iCLONE),MDTemp%D(iCLONE),C%Dyns%TargetTemp)
                ELSE
                   CALL RescaleVelocity(C%Geos%Clone(iCLONE),MDTave%D(iCLONE),C%Dyns%TargetTemp)
                ENDIF
             ENDDO
          ENDIF   
       ENDIF
!      Refresh the Time
       HDFFileID=OpenHDF(C%Nams%HFile)
       DO iCLONE=1,C%Geos%Clones
          MDTime%D(iCLONE) = MDTime%D(iCLONE)+C%Dyns%DTime
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
          CALL Put(MDTime%D(iCLONE),"MDTime")
       ENDDO
!      Archive Geometry for next step
       CALL GeomArchive(iBAS,iGEO+1,C%Nams,C%Sets,C%Geos)     
!      Evaluate energies at the new geometry
       CALL SCF(iBAS,iGEO+1,C)
!      Store the Last P matrix
       CALL RenameDensityMatrix(C,C%Stat%Current%I(1),C%Stat%Current%I(2),C%Stat%Current%I(3))
       IF(C%Stat%Current%I(3)-iREMOVE > 0) THEN
          CALL RemoveDensityMatrix(C,C%Stat%Current%I(3)-iREMOVE-1)
       ENDIF
!      Remove old Stuff from Scratch
       CALL CleanScratch(C,iGEO)
    ENDDO
  END SUBROUTINE MD
!--------------------------------------------------------------
! Set an initial Temperature
!--------------------------------------------------------------
  SUBROUTINE SetTempMaxBoltDist(C,Temp0)
    TYPE(Controls)        :: C
    REAL(DOUBLE)          :: Temp0
    REAL(DOUBLE)          :: Temp,Kin,Mass,TVel,VX,VY,VZ
    INTEGER               :: iCLONE,iATS,I,J,Jmax
!
    Jmax = 20
    DO iCLONE=1,C%Geos%Clones
       DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
          IF(C%Geos%Clone(iCLONE)%CConstrain%I(iATS)==0)THEN
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
!
       MDKin%D(iCLONE)  = Kin
       MDTemp%D(iCLONE) = Temp
       MDTave%D(iCLONE) = Temp
       MDEpot%D(iCLONE) = C%Geos%Clone(iCLONE)%ETotal
       MDEtot%D(iCLONE) = MDEpot%D(iCLONE) + MDKin%D(iCLONE)
!
    ENDDO
!
  END SUBROUTINE SetTempMaxBoltDist
!--------------------------------------------------------------
! The Velocity Verlet Algorithmn : NVE
!--------------------------------------------------------------
  SUBROUTINE MDVerlet_NVE(C,iGEO)
    TYPE(Controls)            :: C
    INTEGER                   :: iGEO
    INTEGER                   :: iCLONE,iATS
    REAL(DOUBLE)              :: Mass,dT,dT2,dTSq2,Time,Dist
    REAL(DOUBLE),DIMENSION(3) :: Pos,Vel,Acc
!--------------------------------------------------------------
!   initialize
    dT    = C%Dyns%DTime
    dT2   = Half*dT
    dTSq2 = Half*dT*dT
!   Clone Loop
    DO iCLONE=1,C%Geos%Clones
!      Move The Atoms
       DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
          IF(C%Geos%Clone(iCLONE)%CConstrain%I(iATS)==0)THEN
             Mass      =  C%Geos%Clone(iCLONE)%AtMss%D(iATS)
             Pos(1:3)  =  C%Geos%Clone(iCLONE)%Carts%D(1:3,iATS)
             Vel(1:3)  =  C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS)
             Acc(1:3)  = -C%Geos%Clone(iCLONE)%Gradients%D(1:3,iATS)/Mass
!            Velocity: v(t) = v(t-dT/2)+a(t)dT/2
             Vel(1:3) = Vel(1:3) + Acc(1:3)*dT2
!            Position: r(t+dT)= r(t)+v(t)*dT+a(t)*dT*dT/2
             Pos(1:3) = Pos(1:3) + Vel(1:3)*dT + Acc(1:3)*dTSq2
!            Update
             C%Geos%Clone(iCLONE)%Carts%D(1:3,iATS)    = Pos(1:3) 
             C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS) = Vel(1:3)
          ENDIF
       ENDDO
!      Calculate Kinectic and  Temp, update Ave Temp
       CALL CalculateMDKin(C%Geos%Clone(iCLONE),MDKin%D(iCLONE),MDTemp%D(iCLONE))
       MDTave%D(iCLONE) = (DBLE(iGEO)/DBLE(iGEO+1))*MDTave%D(iCLONE) &
                         +(       One/DBLE(iGEO+1))*MDTemp%D(iCLONE)
!      Velocity: v(t+dT/2) = v(t)+a(t)dT/2
       DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
          IF(C%Geos%Clone(iCLONE)%CConstrain%I(iATS)==0)THEN
             Mass      =  C%Geos%Clone(iCLONE)%AtMss%D(iATS)
             Vel(1:3)  =  C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS)
             Acc(1:3)  = -C%Geos%Clone(iCLONE)%Gradients%D(1:3,iATS)/Mass
!            Velocity: v(t+dT/2= v(t) + a(t)*dT/2
             Vel(1:3) = Vel(1:3) + Acc(1:3)*dT2
!            Update
             C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS) = Vel(1:3)
          ENDIF
       ENDDO
!      Reset the Momentum to Zero    
       CALL ResetMomentum(C%Geos%Clone(iCLONE),Zero,Zero,Zero)
!      Store Potential and Total Energy
       MDEpot%D(iCLONE) = C%Geos%Clone(iCLONE)%ETotal
       MDEtot%D(iCLONE) = MDEpot%D(iCLONE) + MDKin%D(iCLONE)
!
       C%Geos%Clone(iCLONE)%Carts%D=C%Geos%Clone(iCLONE)%Carts%D
!
    ENDDO
!
    IF(.TRUE.) THEN
       CALL OpenASCII("EnergiesMD.dat",99)
       WRITE(99,'(F10.4,1x,F18.12,1x,F18.12,1x,F18.12)') MDTime%D(1),MDKin%D(1),MDEpot%D(1),MDEtot%D(1)
       WRITE(*,*) "Time = ",MDTime%D(1)," Temperature = ",MDTemp%D(1),' Ave Temp = ',MDTave%D(1)
    ENDIF
!
  END SUBROUTINE MDVerlet_NVE
!--------------------------------------------------------------
! Reset the Velocities so that PX=PX0,PY=PY0,PZ=PZ0
!--------------------------------------------------------------
  SUBROUTINE ResetMomentum(GM,PX0,PY0,PZ0)
    TYPE(CRDS)            :: GM
    REAL(DOUBLE)          :: Mass,PX,PY,PZ,PX0,PY0,PZ0
    INTEGER               :: iATS,nATOMS
!   Calculate Linear Momentum
    PX = Zero
    PY = Zero
    PZ = Zero
    nATOMS = 0
    DO iATS=1,GM%NAtms
       IF(GM%CConstrain%I(iATS)==0)THEN
          nATOMS = nATOMS + 1
          Mass = GM%AtMss%D(iATS)
          PX  = PX + Mass*GM%Velocity%D(1,iATS)
          PY  = PY + Mass*GM%Velocity%D(2,iATS)
          PZ  = PZ + Mass*GM%Velocity%D(3,iATS)
       ENDIF
    ENDDO
!   Reset Linear Momentum
    DO iATS=1,GM%NAtms
       IF(GM%CConstrain%I(iATS)==0)THEN
          Mass = GM%AtMss%D(iATS)
          GM%Velocity%D(1,iATS) = GM%Velocity%D(1,iATS)-(PX-PX0)/(nATOMS*Mass)
          GM%Velocity%D(2,iATS) = GM%Velocity%D(2,iATS)-(PY-PY0)/(nATOMS*Mass)
          GM%Velocity%D(3,iATS) = GM%Velocity%D(3,iATS)-(PZ-PZ0)/(nATOMS*Mass)
       ENDIF
    ENDDO
!
  END SUBROUTINE ResetMomentum
!--------------------------------------------------------------
! Rescale the Velocities to Temp
!--------------------------------------------------------------
  SUBROUTINE RescaleVelocity(GM,Temp,Temp0)
    TYPE(CRDS)            :: GM
    REAL(DOUBLE)          :: Temp,Temp0,Scale
    INTEGER               :: iATS
!  
    Scale = SQRT(Temp0/Temp)
    DO iATS=1,GM%NAtms
       IF(GM%CConstrain%I(iATS)==0)THEN
          GM%Velocity%D(1,iATS) = Scale*GM%Velocity%D(1,iATS)
          GM%Velocity%D(2,iATS) = Scale*GM%Velocity%D(2,iATS)
          GM%Velocity%D(3,iATS) = Scale*GM%Velocity%D(3,iATS)
       ENDIF
    ENDDO
!
  END SUBROUTINE RescaleVelocity
!--------------------------------------------------------------
! Calcualte the Temperature
!--------------------------------------------------------------
  SUBROUTINE CalculateMDKin(GM,Kin,Temp)
    TYPE(CRDS)            :: GM
    REAL(DOUBLE)          :: Kin,Temp,Mass
    INTEGER               :: iATS
!
    Kin    = Zero
    nATOMS = 0
    DO iATS=1,GM%NAtms
       IF(GM%CConstrain%I(iATS)==0)THEN
          nATOMS = nATOMS + 1
          Mass = GM%AtMss%D(iATS)
          Kin  = Kin + Half*GM%AtMss%D(iATS)*(GM%Velocity%D(1,iATS)**2  &
                                             +GM%Velocity%D(2,iATS)**2  &
                                             +GM%Velocity%D(3,iATS)**2)
       ENDIF
    ENDDO
    Temp = (Two/Three)*(Kin/DBLE(nATOMS))*HartreesToKelvin
!   
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
! The Verlet Algorithmn : NPH
!--------------------------------------------------------------
  SUBROUTINE MDVerlet_NPT(C,iGEO)
    TYPE(Controls)            :: C
    INTEGER                   :: iGEO,iCLONE,iATS
    REAL(DOUBLE)              :: Mass,dT
    REAL(DOUBLE),DIMENSION(3) :: Pos,Vel,Acc
!--------------------------------------------------------------
  END SUBROUTINE MDVerlet_NPT
!--------------------------------------------------------------
!--------------------------------------------------------------
! Output
!--------------------------------------------------------------
  SUBROUTINE OutputMD(C,iGEO)
    TYPE(Controls)                 :: C
    INTEGER                        :: iGEO,iCLONE,iATS,I,J,K
    REAL(DOUBLE)                   :: Pressure
    REAL(DOUBLE),DIMENSION(3,3)    :: Latt,InvLatt,LattFrc
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: File,Line
!
    DO iCLONE=1,C%Geos%Clones
!      Open File
       File = TRIM(C%Nams%SCF_NAME)//"_Clone#"//TRIM(IntToChar(iCLONE))//".MDO"
       CALL OpenASCII(File,Out)
!      Add Header
       IF(iGEO==0) THEN
          Line = "##########################################################################"
          WRITE(Out,97) Line
          Line = "# MD Clone No. = "//TRIM(IntToChar(iCLONE))
          WRITE(Out,97) Line
          Line = "# MD Atoms No. = "//TRIM(IntToChar(C%Geos%Clone(iCLONE)%NAtms))
          WRITE(Out,97) Line
          IF(    C%Dyns%MDAlgorithm==VERLET_MD_AL) THEN
             Line = "# MD Algorithm = Verlet"
          ELSEIF(C%Dyns%MDAlgorithm==GEAR_MD_AL) THEN
             Line = "# MD Algorithm = Gear"
          ENDIF
          WRITE(Out,97) Line
          Line = "# DMPOrder     = "//TRIM(IntToChar(C%Opts%DMPOrder))
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
          IF(C%Dyns%Temp_Scaling) THEN
             Line = "# Temperature  : Scaling is On: Target Temp = "//TRIM(IntToChar(INT(C%Dyns%TargetTemp)))//&
                    " K  Scaling Interval = "//TRIM(IntToChar(C%Dyns%RescaleInt))
             WRITE(Out,97) Line
          ELSE
             Line = "# Temperature  : Scaling is Off"
             WRITE(Out,97) Line
          ENDIF
          IF(C%Dyns%Const_Temp) THEN
          ELSE
             Line = "# Anderson T.S.: Scaling is Off"
             WRITE(Out,97) Line
          ENDIF
          IF(C%Dyns%Const_Press) THEN
          ELSE
             Line = "# Volume       : Scaling is Off"
             WRITE(Out,97) Line
          ENDIF
          IF(C%Dyns%Parallel_Rep) THEN
          ELSE
             Line = "# *** Parallel Replicates are Off ***"
             WRITE(Out,97) Line
          ENDIF
          Line = "##########################################################################"
          WRITE(Out,97) Line
          WRITE(Out,*)
!         Set to Zero Stuff that has not been Calculated
          MDEpot%D(iCLONE) = C%Geos%Clone(iCLONE)%ETotal
          MDEtot%D(iCLONE) = MDEpot%D(iCLONE) + MDKin%D(iCLONE)
          C%Geos%Clone(iCLONE)%PBC%LatFrc%D = Zero
       ENDIF
!      Add MD Timestep Configuration
       Line = "-------------------------------------------MD-"//TRIM(IntToChar(iGEO))// &
              "-MD-------------------------------------------"
       WRITE(Out,97) Line
       WRITE(Out,98) "MD Time        = ",MDTime%D(iCLONE)
       WRITE(Out,98) "MD Kinetic     = ",MDKin%D(iCLONE)
       WRITE(Out,98) "MD Potential   = ",MDEPot%D(iCLONE)
       WRITE(Out,98) "MD Total       = ",MDEtot%D(iCLONE)
       WRITE(Out,98) "MD Temperature = ",MDTemp%D(iCLONE)
       WRITE(Out,96) "MD L. Momentum = ",MDLinP%D(1:3,iCLONE)
!      Output Lattice and Lattice Forces
       Latt    = C%Geos%Clone(iCLONE)%PBC%BoxShape%D
       InvLatt = C%Geos%Clone(iCLONE)%PBC%InvBoxSh%D
       LattFrc = C%Geos%Clone(iCLONE)%PBC%LatFrc%D
       IF(C%Geos%Clone(iCLONE)%PBC%Dimen==1) THEN
          IF(C%Geos%Clone(iCLONE)%PBC%AutoW%I(1)==1) THEN
             WRITE(Out,85) 
             WRITE(Out,82) "       a       = ",Latt(1,1)
             WRITE(Out,85) 
             WRITE(Out,82) "     div(a)    = ",LattFrc(1,1)
             WRITE(Out,85) 
          ENDIF
          IF(C%Geos%Clone(iCLONE)%PBC%AutoW%I(2)==1) THEN
             WRITE(Out,85) 
             WRITE(Out,82) "       b       = ",Latt(2,2)
             WRITE(Out,85) 
             WRITE(Out,82) "     div(b)    = ",LattFrc(2,2)
             WRITE(Out,85) 
          ENDIF
          IF(C%Geos%Clone(iCLONE)%PBC%AutoW%I(3)==1) THEN
             WRITE(Out,85) 
             WRITE(Out,82) "       c       = ",Latt(3,3)
             WRITE(Out,85) 
             WRITE(Out,82) "     div(c)    = ",LattFrc(3,3)
             WRITE(Out,85) 
          ENDIF
       ELSEIF(C%Geos%Clone(iCLONE)%PBC%Dimen==2) THEN
          IF(C%Geos%Clone(iCLONE)%PBC%AutoW%I(1)==0) THEN
             WRITE(Out,85) 
             WRITE(Out,82) "       b       = ",Latt(2,2), Latt(2,3)
             WRITE(Out,82) "       c       = ",Latt(3,2), Latt(3,3)
             WRITE(Out,85) 
             WRITE(Out,82) "     div(b)    = ",LattFrc(2,2), LattFrc(2,3)
             WRITE(Out,82) "     div(c)    = ",LattFrc(3,2), LattFrc(3,3)
             WRITE(Out,85) 
          ENDIF
          IF(C%Geos%Clone(iCLONE)%PBC%AutoW%I(2)==0) THEN
             WRITE(Out,85) 
             WRITE(Out,82) "       a       = ",Latt(1,1), Latt(1,3)
             WRITE(Out,82) "       c       = ",Latt(3,1), Latt(3,3)
             WRITE(Out,85) 
             WRITE(Out,82) "     div(a)    = ",LattFrc(1,1), LattFrc(1,3)
             WRITE(Out,82) "     div(c)    = ",LattFrc(3,1), LattFrc(3,3)
             WRITE(Out,85) 
          ENDIF
          IF(C%Geos%Clone(iCLONE)%PBC%AutoW%I(3)==0) THEN
             WRITE(Out,85) 
             WRITE(Out,82) "       a       = ",Latt(1,1), Latt(1,2)
             WRITE(Out,82) "       b       = ",Latt(2,1), Latt(2,2)
             WRITE(Out,85) 
             WRITE(Out,82) "     div(a)    = ",LattFrc(1,1), LattFrc(1,2)
             WRITE(Out,82) "     div(b)    = ",LattFrc(2,1), LattFrc(2,2)
             WRITE(Out,85) 
          ENDIF
       ELSEIF(C%Geos%Clone(iCLONE)%PBC%Dimen==3) THEN
!         Compute The Pressure
          Pressure = (2.D0/3.D0)*MDKin%D(iCLONE)
          DO I=1,3
             DO J=1,3
                Pressure = Pressure + Latt(I,J)*LattFrc(I,J)
             ENDDO
          ENDDO
          Pressure=Pressure/(C%Geos%Clone(iCLONE)%PBC%CellVolume*GPaToAU)
!
          WRITE(Out,85) 
          WRITE(Out,82) "       a       = ",Latt(1,1), Latt(1,2), Latt(1,3)
          WRITE(Out,82) "       b       = ",Latt(2,1), Latt(2,2), Latt(2,3)
          WRITE(Out,82) "       c       = ",Latt(3,1), Latt(3,2), Latt(3,3)
          WRITE(Out,85) 
          WRITE(Out,82) "     div(a)    = ",LattFrc(1,1), LattFrc(1,2), LattFrc(1,3)
          WRITE(Out,82) "     div(b)    = ",LattFrc(2,1), LattFrc(2,2), LattFrc(2,3)
          WRITE(Out,82) "     div(c)    = ",LattFrc(3,1), LattFrc(3,2), LattFrc(3,3)
          WRITE(Out,85) 
          WRITE(Out,80) "    Pressure   = ",Pressure," GPa"
          WRITE(Out,85) 
       ENDIF
       DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
          WRITE(Out,99) TRIM(C%Geos%Clone(iCLONE)%AtNam%C(iATS)),  &
                        C%Geos%Clone(iCLONE)%Carts%D(1:3,iATS),    &
                        C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS)
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
99  FORMAT(a3,6(1x,F14.8))
  END SUBROUTINE OutputMD
!--------------------------------------------------------------
! Rename the Last Density Matrix, Remove old ones
!--------------------------------------------------------------
  SUBROUTINE RenameDensityMatrix(C,iSCF,iBAS,iGEO)
    TYPE(Controls)                 :: C
    INTEGER                        :: iCLONE,iSCF,iBAS,iGEO
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: chSCF,chBAS,chGEO,chCLONE
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: PoldFile,PnewFile
!      
!!$    chSCF = IntToChar(C%Stat%Current%I(1))
!!$    chBAS = IntToChar(C%Stat%Current%I(2))
!!$    chGEO = IntToChar(C%Stat%Current%I(3))
!
    chSCF = IntToChar(iSCF)
    chBAS = IntToChar(iBAS)
    chGEO = IntToChar(iGEO)
!
    DO iCLONE=1,C%Geos%Clones
       chCLONE = IntToChar(iCLONE)
       PoldFile = TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)  &
                                        //'_Geom#'//TRIM(chGEO)  &
                                        //'_Base#'//TRIM(chBAS)  &
                                        //'_Cycl#'//TRIM(chSCF)  &
                                        //'_Clone#'//TRIM(chCLONE)//'.D'
       PnewFile = TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)  &
                                        //'_G#'//TRIM(chGEO)     &
                                        //'_C#'//TRIM(chCLONE)//'.Dsave'
       CALL SYSTEM('/bin/cp -f  '//PoldFile//' '//PnewFile)
    ENDDO
!
  END SUBROUTINE RenameDensityMatrix
!--------------------------------------------------------------
! Rename the Last Density Matrix, Remove old ones
!--------------------------------------------------------------
  SUBROUTINE RemoveDensityMatrix(C,iGEO)
    TYPE(Controls)                 :: C
    INTEGER                        :: iCLONE,iGEO
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: chSCF,chBAS,chGEO,chCLONE
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: PoldFile,PnewFile
!
    DO iCLONE=1,C%Geos%Clones
       chCLONE = IntToChar(iCLONE)
       chGEO   = IntToChar(iGEO)
       PoldFile = TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)  &
            //'_G#'//TRIM(chGEO)     &
            //'_C#'//TRIM(chCLONE)//'.Dsave'
       CALL SYSTEM('/bin/rm -f  '//PoldFile)
    ENDDO
!
  END SUBROUTINE RemoveDensityMatrix
!--------------------------------------------------------------
! Rename the Last Density Matrix, Remove old ones
!--------------------------------------------------------------
  SUBROUTINE CopyDensityMatrix(C,iGEO,nGEO)
    TYPE(Controls)                 :: C
    INTEGER                        :: iCLONE,iGEO,nGEO
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: chGEO,chCLONE
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: PoldFile,PnewFile
!
    DO iCLONE=1,C%Geos%Clones
       chCLONE = IntToChar(iCLONE)
       chGEO   = IntToChar(iGEO)
       PoldFile = TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)  &
            //'_G#'//TRIM(chGEO)     &
            //'_C#'//TRIM(chCLONE)//'.Dsave'
       chGEO   = IntToChar(nGEO)
       PnewFile = TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)  &
            //'_G#'//TRIM(chGEO)     &
            //'_C#'//TRIM(chCLONE)//'.Dsave'
       CALL SYSTEM('/bin/cp -f  '//PoldFile//' '//PnewFile)
    ENDDO
!
  END SUBROUTINE CopyDensityMatrix
!--------------------------------------------------------------
!
!
!--------------------------------------------------------------
  SUBROUTINE OptimizeLattice(C)
    TYPE(Controls)  :: C
    INTEGER         :: iBAS,iGEO,iCLONE,I,J
!
    iBAS   = C%Sets%NBSets
    iCLONE = 1
    DO iGEO = 1,50
!      Determine a New Lattice
       CALL OpenASCII("LatticeOpt.dat",77)
       WRITE(77,*) 'Geom   = ',iGEO
       WRITE(77,*) "Energy = ",C%Geos%Clone(iCLONE)%ETotal
       WRITE(77,*) 'Lattice'
       DO I=1,3
          WRITE(77,*) (C%Geos%Clone(iCLONE)%PBC%BoxShape%D(I,J)*BohrsToAngstroms,J=1,3)
       ENDDO
       WRITE(77,*) 'LatFrc'
       DO I=1,3
          WRITE(77,*) (C%Geos%Clone(iCLONE)%PBC%LatFrc%D(I,J)*BohrsToAngstroms,J=1,3)
       ENDDO
       CLOSE(77)
!
       WRITE(*,*) 'Geom   = ',iGEO
       WRITE(*,*) "Energy = ",C%Geos%Clone(iCLONE)%ETotal
       WRITE(*,*) 'Lattice'
       DO I=1,3
          WRITE(*,*) (C%Geos%Clone(iCLONE)%PBC%BoxShape%D(I,J)*BohrsToAngstroms,J=1,3)
       ENDDO
       WRITE(*,*) 'LatFrc'
       DO I=1,3
          WRITE(*,*) (C%Geos%Clone(iCLONE)%PBC%LatFrc%D(I,J)*BohrsToAngstroms,J=1,3)
       ENDDO
!
       DO I=1,3
          DO J=I,3
             C%Geos%Clone(iCLONE)%PBC%BoxShape%D(I,J) = C%Geos%Clone(iCLONE)%PBC%BoxShape%D(I,J)    &
                                                      - 5.0D0*C%Geos%Clone(iCLONE)%PBC%LatFrc%D(I,J)
          ENDDO
       ENDDO
!      Archive Geometry for next step
       CALL MakeGMPeriodic(C%Geos%Clone(iCLONE))
       C%Geos%Clone(1)%Carts%D = C%Geos%Clone(iCLONE)%Carts%D
       CALL GeomArchive(iBAS,iGEO+1,C%Nams,C%Sets,C%Geos)     
!      Evaluate energies at the new geometry
       CALL SCF(iBAS,iGEO+1,C)
!      Calculate Force
       CALL Force(iBAS,iGEO+1,C%Nams,C%Opts,C%Stat,C%Geos,C%Sets,C%MPIs)
    ENDDO
!
  END SUBROUTINE OptimizeLattice
!
!
!
END MODULE MDynamics
