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
  TYPE(DBL_RNK3)      :: Carts0
  CONTAINS
!--------------------------------------------------------------
! Main MD Subroutine
!--------------------------------------------------------------
 SUBROUTINE  MD(C)  
    TYPE(Controls)  :: C
    INTEGER         :: iBAS,iGEO,iCLONE,I
    REAL(DOUBLE)    :: Temp 
!--------------------------------------------------------------
!   Do Molecular Dynamics:Loop over Time Steps
!   Intitialize
    C%Stat%Previous%I=(/0,1,1/)
    iGEO=1
    CALL New(Carts0,(/3,C%Geos%Clone(1)%NAtms,C%Geos%Clones/))
    Carts0%D=Zero
    CALL New(MDTime ,C%Geos%Clones)
    CALL New(MDKin  ,C%Geos%Clones)
    CALL New(MDEpot ,C%Geos%Clones)
    CALL New(MDETot ,C%Geos%Clones)
    CALL New(MDTemp ,C%Geos%Clones)
!   Initial Guess
    IF(C%Opts%Guess==GUESS_EQ_RESTART) THEN
!      Init the Time
       HDFFileID=OpenHDF(C%Nams%RFile)
       DO iCLONE=1,C%Geos%Clones
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
          CALL Get(MDTime%D(iCLONE),"MDTime")
       ENDDO
!      Do The initial SCF
       iBAS=C%Sets%NBSets
       CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos) 
       CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
       CALL SCF(iBAS,iGEO,C)
    ELSE
!      Hack, Give an intial Maxwell Boltzman Temp
       Temp=20.D0
       CALL SetTempMaxBoltDist(C,Temp)
!      Init the Time
       MDTime%D(:) = Zero
       HDFFileID=OpenHDF(C%Nams%HFile)
       DO iCLONE=1,C%Geos%Clones
          HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
          CALL Put(MDTime%D(iCLONE),"MDTime")
       ENDDO
!      Build the SCF 
       DO iBAS=1,C%Sets%NBSets
          CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos) 
          CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
          CALL SCF(iBAS,iGEO,C)
       ENDDO
    ENDIF
!   Do MD    
    iBAS=C%Sets%NBSets
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
       IF(.TRUE.) THEN
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
!      Remove old Stuff from Scratch
       CALL CleanScratch(C,iGEO)
    ENDDO
  END SUBROUTINE MD
!--------------------------------------------------------------
! The Verlet Algorithmn : NVE
!--------------------------------------------------------------
  SUBROUTINE MDVerlet_NVE(C,iGEO)
    TYPE(Controls)            :: C
    INTEGER                   :: iGEO
    INTEGER                   :: iCLONE,iATS
    REAL(DOUBLE)              :: Mass,dT,Time,Dist
    REAL(DOUBLE),DIMENSION(3) :: Pos0,Pos1,Pos2,Vel,Acc
!--------------------------------------------------------------
!   initialize
    dT = C%Dyns%DTime
!   Clone Loop
    DO iCLONE=1,C%Geos%Clones
!      Move The Atoms
       MDKin%D(iCLONE) = Zero
       DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
          Mass      =  C%Geos%Clone(iCLONE)%AtMss%D(iATS)
          Pos1(1:3) =  C%Geos%Clone(iCLONE)%AbCarts%D(1:3,iATS)
          Acc(1:3)  = -C%Geos%Clone(iCLONE)%Gradients%D(1:3,iATS)/Mass
          IF(iGEO==1) THEN
             Vel(1:3)  = C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS)
             Pos0(1:3) = Pos1(1:3)-Vel(1:3)*dT
          ELSE
             Pos0(1:3) = Carts0%D(1:3,iATS,iCLONE)
          ENDIF
!         r(t+dT) = 2*r(t)-r(t-dT) + a(t)*dT*dT
          Pos2(1:3) = Two*Pos1(1:3) - Pos0(1:3) + Acc(1:3)*dT*dT
!         Velocity
          Vel(1:3) =  (Pos2(1:3)-Pos0(1:3))/(Two*dT)
!         Calculate EKin
          MDKin%D(iCLONE) = MDKin%D(iCLONE) + Half*Mass*(Vel(1)**2+Vel(2)**2+Vel(3)**2)
!         Update
          Carts0%D(1:3,iATS,iCLONE)                 = Pos1(1:3)
          C%Geos%Clone(iCLONE)%AbCarts%D(1:3,iATS)  = Pos2(1:3)
          C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS) = Vel(1:3)
       ENDDO
!
       MDEpot%D(iCLONE) = C%Geos%Clone(iCLONE)%ETotal
       MDEtot%D(iCLONE) = MDEpot%D(iCLONE)+MDKin%D(iCLONE)
       MDTemp%D(iCLONE)= (Two/Three)*MDKin%D(iCLONE)/DBLE(C%Geos%Clone(iCLONE)%NAtms)*HartreesToKelvin
!
       CALL OpenASCII("EnergiesMD.dat",99)
       WRITE(99,'(F12.4,F14.8,F14.8,F14.8)') MDTime%D(iCLONE),MDKin%D(iCLONE),MDEpot%D(iCLONE),MDEtot%D(iCLONE)
       CLOSE(99)
       WRITE(*,*) "Time = ",MDTime%D(iCLONE)," Temperature = ",MDTemp%D(iCLONE)
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
          Line = "# MD Time Step = "//TRIM(DblToMedmChar(C%Dyns%DTime))
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
       DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
          WRITE(Out,99) TRIM(C%Geos%Clone(iCLONE)%AtNam%C(iATS)),  &
                        Carts0%D(1:3,iATS,iCLONE),   &
                        C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS)
       ENDDO
!!$       DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
!!$          WRITE(Out,99) TRIM(C%Geos%Clone(iCLONE)%AtNam%C(iATS)),  &
!!$                        Carts0%D(1:3,iATS,iCLONE)/AngstromsToAU,   &
!!$                        C%Geos%Clone(iCLONE)%Velocity%D(1:3,iATS)/AngstromsToAU
!!$       ENDDO
       CLOSE(Out)
    ENDDO
97  FORMAT(a128)
98  FORMAT(a18,F16.10)
99  FORMAT(a3,6(1x,F14.8))
  END SUBROUTINE OutputMD
!--------------------------------------------------------------
! Set an initial Temperature
!--------------------------------------------------------------
  SUBROUTINE SetTempMaxBoltDist(C,Temp)
    TYPE(Controls)        :: C
    REAL(DOUBLE)          :: Temp,Mass,TVel,VX,VY,VZ,SX,SY,SZ
    INTEGER               :: iCLONE,iATS,I,J,Jmax
!
    Jmax = 20
    DO iCLONE=1,C%Geos%Clones
       SX = Zero
       SY = Zero
       SZ = Zero
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
          SX = SX + VX/DBLE(Jmax)
          SY = SX + VY/DBLE(Jmax)         
          SZ = SX + VZ/DBLE(Jmax)
       ENDDO
!      Reset Linear Momentum
       DO iATS=1,C%Geos%Clone(iCLONE)%NAtms
          C%Geos%Clone(iCLONE)%Velocity%D(1,iATS) = C%Geos%Clone(iCLONE)%Velocity%D(1,iATS)-VX/C%Geos%Clone(iCLONE)%NAtms
          C%Geos%Clone(iCLONE)%Velocity%D(2,iATS) = C%Geos%Clone(iCLONE)%Velocity%D(2,iATS)-VY/C%Geos%Clone(iCLONE)%NAtms
          C%Geos%Clone(iCLONE)%Velocity%D(3,iATS) = C%Geos%Clone(iCLONE)%Velocity%D(3,iATS)-VZ/C%Geos%Clone(iCLONE)%NAtms
       ENDDO
    ENDDO
    CALL RescaleTemp(C,Temp)
!
  END SUBROUTINE SetTempMaxBoltDist
!--------------------------------------------------------------
! Set an initial Temperature
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
END MODULE MDynamics
