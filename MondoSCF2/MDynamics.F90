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
  TYPE(DBL_VECT)      :: MDTime,MDKin,MDEpot,MDEtot,MDTemp
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
    CALL New(MDLinP ,(/3,C%Geos%Clones/))
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
          CALL Put( C%Opts%DMPOrder,"DMPOrder")
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
    iBAS=C%Sets%NBSets
!   Store the Last P matrix
    CALL RenameDensityMatrix(C,iREMOVE)
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
!      Refresh the Time
       HDFFileID=OpenHDF(C%Nams%HFile)
       DO iCLONE=1,C%Geos%Clones
          MDTime%D(iCLONE) = MDTime%D(iCLONE)+C%Dyns%DTime
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
          CALL Put(MDTime%D(iCLONE),"MDTime")
       ENDDO
       IF(.FALSE.) THEN
          C%Opts%Guess=GUESS_EQ_SUPR
          DO iBAS=1,C%Sets%NBSets
             CALL GeomArchive(iBAS,iGEO+1,C%Nams,C%Sets,C%Geos) 
             CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
             CALL SCF(iBAS,iGEO+1,C)
          ENDDO
          iBAS=C%Sets%NBSets
       ELSE
!         Archive Geometry for next step
          CALL GeomArchive(iBAS,iGEO+1,C%Nams,C%Sets,C%Geos)     
!         Evaluate energies at the new geometry
          CALL SCF(iBAS,iGEO+1,C)
       ENDIF
!      Store the Last P matrix
       CALL RenameDensityMatrix(C,iREMOVE)
!      Remove old Stuff from Scratch
       CALL CleanScratch(C,iGEO)
    ENDDO
  END SUBROUTINE MD
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
       MDKin%D(iCLONE)      = Zero
       MDLinP%D(1:3,iCLONE) = Zero
       DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
          Mass      =  C%Geos%Clone(iCLONE)%AtMss%D(iATS)
          Pos(1:3)  =  C%Geos%Clone(iCLONE)%AbCarts%D(1:3,iATS)
          Vel(1:3)  =  C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS)
          Acc(1:3)  = -C%Geos%Clone(iCLONE)%Gradients%D(1:3,iATS)/Mass
!         Velocity: v(t) = v(t-dT/2(t)*dT/2
          Vel(1:3) = Vel(1:3) + Acc(1:3)*dT2
!         Calculate EKin
          MDKin%D(iCLONE) = MDKin%D(iCLONE) + Half*Mass*(Vel(1)**2+Vel(2)**2+Vel(3)**2)
!         Calculate Linear Momentum
          MDLinP%D(1:3,iCLONE) = MDLinP%D(1:3,iCLONE)+Vel(1:3)*Mass
!         Position: r(t+dT)= r(t)+v(t)*dT+a(t)*dT*dT/2
          Pos(1:3) = Pos(1:3) + Vel(1:3)*dT + Acc(1:3)*dTSq2
!         Velocity: v(t+dT/2= v(t) + a(t)*dT/2
          Vel(1:3) = Vel(1:3) + Acc(1:3)*dT2
!         Update
          C%Geos%Clone(iCLONE)%AbCarts%D(1:3,iATS)  = Pos(1:3) 
          C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS) = Vel(1:3)
       ENDDO
!
       MDEpot%D(iCLONE) = C%Geos%Clone(iCLONE)%ETotal
       MDEtot%D(iCLONE) = MDEpot%D(iCLONE)+MDKin%D(iCLONE)
       MDTemp%D(iCLONE)= (Two/Three)*MDKin%D(iCLONE)/DBLE(C%Geos%Clone(iCLONE)%NAtms)*HartreesToKelvin
!
       IF(.TRUE.) THEN
          CALL OpenASCII("EnergiesMD.dat",99)
          WRITE(99,'(F10.4,1x,F18.12,1x,F18.12,1x,F18.12)') MDTime%D(iCLONE),MDKin%D(iCLONE),MDEpot%D(iCLONE),MDEtot%D(iCLONE)
          CLOSE(99)
!
          CALL OpenASCII("PositionMD.dat",99)
          Pos(1:3) = C%Geos%Clone(iCLONE)%AbCarts%D(1:3,1)-C%Geos%Clone(iCLONE)%AbCarts%D(1:3,2)
          Dist =  SQRT(Pos(1)**2+Pos(2)**2+Pos(3)**2)
          WRITE(99,'(F10.4,1x,F18.12)') MDTime%D(iCLONE),Dist
          CLOSE(99) 
!
          WRITE(*,*) "Time = ",MDTime%D(iCLONE)," Temperature = ",MDTemp%D(iCLONE)
       ENDIF
!
    ENDDO
  END SUBROUTINE MDVerlet_NVE
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
    INTEGER                        :: iGEO,iCLONE,iATS
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: File,Line
!
    DO iCLONE=1,C%Geos%Clones
!      Open File
       File = TRIM(C%Nams%SCF_NAME)//"_Clone#"//TRIM(IntToChar(iCLONE))//".MDO"
       CALL OpenASCII(File,Out)
!      Add Header
       IF(iGEO==1) THEN
          Line = "##################################################"
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
          IF(C%Dyns%Velcty_Scaling) THEN
          ELSE
             Line = "# Velocity     : Scaling is Off"
             WRITE(Out,97) Line
          ENDIF
          IF(C%Dyns%Const_Temp) THEN
          ELSE
             Line = "# Temperature  : Scaling is Off"
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
          Line = "##################################################"
          WRITE(Out,97) Line
          WRITE(Out,*)
       ENDIF
!      Add MD Timestep Configuration
       Line = "----------------------------------------------"//TRIM(IntToChar(iGEO))// &
              "----------------------------------------------"
       WRITE(Out,97) Line
       WRITE(Out,98) "MD Time        = ",MDTime%D(iCLONE)
       WRITE(Out,98) "MD Kinetic     = ",MDKin%D(iCLONE)
       WRITE(Out,98) "MD Potential   = ",MDEPot%D(iCLONE)
       WRITE(Out,98) "MD Total       = ",MDEtot%D(iCLONE)
       WRITE(Out,98) "MD Temperature = ",MDTemp%D(iCLONE)
       WRITE(Out,96) "MD L. Momentum = ",MDLinP%D(1:3,iCLONE)
       DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
          WRITE(Out,99) TRIM(C%Geos%Clone(iCLONE)%AtNam%C(iATS)),  &
                        C%Geos%Clone(iCLONE)%AbCarts%D(1:3,iATS),  &
                        C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS)
       ENDDO
       CLOSE(Out)
    ENDDO
97  FORMAT(a128)
98  FORMAT(a18,F16.10)
96  FORMAT(a18,3(F16.10,1x))
99  FORMAT(a3,6(1x,F14.8))
  END SUBROUTINE OutputMD
!--------------------------------------------------------------
! Set an initial Temperature
!--------------------------------------------------------------
  SUBROUTINE SetTempMaxBoltDist(C,Temp)
    TYPE(Controls)        :: C
    REAL(DOUBLE)          :: Temp,Mass,TVel,VX,VY,VZ
    INTEGER               :: iCLONE,iATS,I,J,Jmax
!
    Jmax = 20
    DO iCLONE=1,C%Geos%Clones
       DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
          Mass  =  C%Geos%Clone(iCLONE)%AtMss%D(iATS) 
          TVel  = SQRT(Three*Temp*KelvinToHartrees/Mass)
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
       ENDDO
    ENDDO
    CALL ResetMomentum(C,Zero,Zero,Zero)
    CALL RescaleTemp(C,Temp)
!
  END SUBROUTINE SetTempMaxBoltDist
!--------------------------------------------------------------
! Reset the Velocities so that PX=PX0,PY=PY0,PZ=PZ0
!--------------------------------------------------------------
  SUBROUTINE ResetMomentum(C,PX0,PY0,PZ0)
    TYPE(Controls)        :: C
    REAL(DOUBLE)          :: Mass,PX,PY,PZ,PX0,PY0,PZ0
    INTEGER               :: iCLONE,iATS
!
    DO iCLONE=1,C%Geos%Clones
!      Calculate Linear Momentum
       PX = Zero
       PY = Zero
       PZ = Zero
       DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
          Mass = C%Geos%Clone(iCLONE)%AtMss%D(iATS)
          PX  = PX + Mass*C%Geos%Clone(iCLONE)%Velocity%D(1,iATS)
          PY  = PY + Mass*C%Geos%Clone(iCLONE)%Velocity%D(2,iATS)
          PZ  = PZ + Mass*C%Geos%Clone(iCLONE)%Velocity%D(3,iATS)
       ENDDO
!      Reset Linear Momentum
       DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
          Mass = C%Geos%Clone(iCLONE)%AtMss%D(iATS)
          C%Geos%Clone(iCLONE)%Velocity%D(1,iATS) = C%Geos%Clone(iCLONE)%Velocity%D(1,iATS)-(PX-PX0)/(C%Geos%Clone(iCLONE)%NAtms*Mass)
          C%Geos%Clone(iCLONE)%Velocity%D(2,iATS) = C%Geos%Clone(iCLONE)%Velocity%D(2,iATS)-(PY-PY0)/(C%Geos%Clone(iCLONE)%NAtms*Mass)
          C%Geos%Clone(iCLONE)%Velocity%D(3,iATS) = C%Geos%Clone(iCLONE)%Velocity%D(3,iATS)-(PZ-PZ0)/(C%Geos%Clone(iCLONE)%NAtms*Mass)
       ENDDO
    ENDDO

  END SUBROUTINE ResetMomentum
!--------------------------------------------------------------
! Rescale the Velocities to Temp
!--------------------------------------------------------------
  SUBROUTINE RescaleTemp(C,Temp)
    TYPE(Controls)        :: C
    REAL(DOUBLE)          :: Temp,Mass,Temp0,SUMV,VX,VY,VZ,Scale
    INTEGER               :: iCLONE,iATS,I,J,Jmax
!  
    DO iCLONE=1,C%Geos%Clones
!      Determine the Temp
       SUMV = Zero
       DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
          Mass = C%Geos%Clone(iCLONE)%AtMss%D(iATS)
          VX = C%Geos%Clone(iCLONE)%Velocity%D(1,iATS) 
          VY = C%Geos%Clone(iCLONE)%Velocity%D(2,iATS) 
          VZ = C%Geos%Clone(iCLONE)%Velocity%D(3,iATS) 
          SUMV = SUMV+Half*Mass*(VX**2+VY**2+VZ**2)
       ENDDO
       Temp0 =  (Two/Three)*SUMV/DBLE(C%Geos%Clone(iCLONE)%NAtms)*HartreesToKelvin
       Scale = SQRT(Temp/Temp0)
       DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
          C%Geos%Clone(iCLONE)%Velocity%D(1,iATS) = Scale*C%Geos%Clone(iCLONE)%Velocity%D(1,iATS)
          C%Geos%Clone(iCLONE)%Velocity%D(2,iATS) = Scale*C%Geos%Clone(iCLONE)%Velocity%D(2,iATS)
          C%Geos%Clone(iCLONE)%Velocity%D(3,iATS) = Scale*C%Geos%Clone(iCLONE)%Velocity%D(3,iATS)
       ENDDO
!
    ENDDO
  END SUBROUTINE RescaleTemp
!--------------------------------------------------------------
! Rename the Last Density Matrix, Remove old ones
!--------------------------------------------------------------
  SUBROUTINE RenameDensityMatrix(C,iREMOVE)
    TYPE(Controls)                 :: C
    INTEGER                        :: iCLONE,iREMOVE
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: chSCF,chBAS,chGEO,chCLONE
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: PoldFile,PnewFile
!      
    chSCF = IntToChar(C%Stat%Current%I(1))
    chBAS = IntToChar(C%Stat%Current%I(2))
    chGEO = IntToChar(C%Stat%Current%I(3))
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
    IF(C%Stat%Current%I(3)-iREMOVE > 0) THEN
       DO iCLONE=1,C%Geos%Clones
          chCLONE = IntToChar(iCLONE)
          chGEO   = IntToChar(C%Stat%Current%I(3)-iREMOVE)
          PoldFile = TRIM(C%Nams%M_SCRATCH)//TRIM(C%Nams%SCF_NAME)  &
                                           //'_G#'//TRIM(chGEO)     &
                                           //'_C#'//TRIM(chCLONE)//'.Dsave'
          CALL SYSTEM('/bin/rm -f  '//PoldFile)
       ENDDO
    ENDIF
!
  END SUBROUTINE RenameDensityMatrix
!
!
!
END MODULE MDynamics
