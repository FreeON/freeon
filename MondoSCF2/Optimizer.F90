MODULE Optimizer
  USE SCFs
  USE InOut
  USE MatFunk
  USE Numerics
  USE AtomPairs
  USE ControlStructures  
  USE InCoords
  USE GDIISMod
  USE GeomOptKeys     
  USE PunchHDF
  USE GlobalScalars
  IMPLICIT NONE
CONTAINS
  !=====================================================================================
  !
  !=====================================================================================  
  SUBROUTINE Descender(C)
    TYPE(Controls) :: C
    !----------------------------------------------------------------------------------!
    IF(C%Opts%Coordinates==GRAD_CART_OPT) THEN
      IF(C%Opts%DoGDIIS)THEN
         ! Follow extrapolated Cartesian gradient down hill
         CALL GDicer(C)
      ELSE
         ! Follow Cartesian gradient down hill
         CALL SteepD(C)
      ENDIF
    ELSE
      CALL IntOpt(C)
    ENDIF 
  END SUBROUTINE Descender
  !=====================================================================================
  !
  !=====================================================================================  
  SUBROUTINE SteepD(C)
    TYPE(Controls) :: C
    INTEGER        :: iBAS,iGEO,iCLONE
    REAL(DOUBLE),DIMENSION(C%Geos%Clones,C%Opts%NSteps) :: Energy
    INTEGER        :: NatmsLoc,IStart
    TYPE(DBL_RNK2) :: OldXYZ
    !-------------------------------------------------------------------------!
    ! initial geometry
    iGEO=C%Stat%Previous%I(3)
    NatmsLoc=C%Geos%Clone(1)%Natms
    ! Build the guess 
    DO iBAS=1,C%Sets%NBSets
       CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)    
       CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
       CALL SCF(iBAS,iGEO,C)
    ENDDO
    ! Print the starting coordinates and energy
    DO iCLONE=1,C%Geos%Clones
       CALL PPrint(C%Geos%Clone(iCLONE),C%Nams%GFile,Geo,C%Opts%GeomPrint)
    ENDDO
    ! Follow the gradient down hill
    iBAS=C%Sets%NBSets
    IStart=iGEO
    DO iGEO=IStart,C%Opts%NSteps
       ! Compute new gradients
       CALL Force(iBAS,iGEO,C%Nams,C%Opts,C%Stat,C%Geos,C%Sets,C%MPIs)
         DO iCLONE=1,C%Geos%Clones
           C%Geos%Clone(iCLONE)%Displ%D=C%Geos%Clone(iCLONE)%AbCarts%D
         ENDDO
       IF(SteepStep(iBAS,iGEO,Energy(:,iGEO),C))THEN
          DO iCLONE=1,C%Geos%Clones
             IF(Energy(iCLONE,iGEO+1)<Energy(iCLONE,iGEO))THEN
                ! On convergence only write the lowest energy geometry
                CALL PPrint(C%Geos%Clone(iCLONE),C%Nams%GFile,Geo,C%Opts%GeomPrint)
             ENDIF  
          ENDDO
          EXIT
       ELSE
          DO iCLONE=1,C%Geos%Clones
             CALL PPrint(C%Geos%Clone(iCLONE),C%Nams%GFile,Geo,C%Opts%GeomPrint)
          ENDDO
       ENDIF
         CALL New(OldXYZ,(/3,NatmsLoc/))
         DO iCLONE=1,C%Geos%Clones
           OldXYZ%D=C%Geos%Clone(iCLONE)%Displ%D
           C%Geos%Clone(iCLONE)%Displ%D=C%Geos%Clone(iCLONE)%AbCarts%D
           C%Geos%Clone(iCLONE)%AbCarts%D=OldXYZ%D
         ENDDO
         CALL Delete(OldXYZ)
       CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)    
         DO iCLONE=1,C%Geos%Clones
           C%Geos%Clone(iCLONE)%AbCarts%D=C%Geos%Clone(iCLONE)%Displ%D
         ENDDO
       CALL GeomArchive(iBAS,iGEO+1,C%Nams,C%Sets,C%Geos)    
    ENDDO
  END SUBROUTINE SteepD
  !=====================================================================================
  !
  !=====================================================================================  
  SUBROUTINE GDicer(C)
    TYPE(Controls)         :: C
    TYPE(DBL_VECT)         :: GradMax,GradRMS,Grad
    REAL(DOUBLE)           :: DIISErr,GRMSQ,GMAXQ
    ! Initial step 
    REAL(DOUBLE),PARAMETER :: StepLength=2D0 
    INTEGER                :: iBAS,iGEO,iCLONE,AccL,IStart
    INTEGER                :: Relaxations=3   ! This should be an input variable at some point
    LOGICAL                :: Converged,Steep
    CHARACTER(LEN=DCL)     :: Mssg
    TYPE(DBL_RNK2)         :: OldXYZ,NewXYZ
    INTEGER                :: NatmsLoc
    !----------------------------------------------------------------------------------!
    ! Start with the first geometry
    iGEO=C%Stat%Previous%I(3)
    NatmsLoc=C%Geos%Clone(1)%Natms
    ! Build the guess 
    DO iBAS=1,C%Sets%NBSets
       CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)    
       CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
       CALL SCF(iBAS,iGEO,C)
    ENDDO    
    ! Space for accumulating convergence statistics
    CALL New(GradMax,C%Opts%NSteps)
    CALL New(GradRMS,C%Opts%NSteps)
    ! Use just the last basis 
    iBAS=C%Sets%NBSets
    AccL=C%Opts%AccuracyLevels(iBAS)
    ! Use simple Cartesian GDIIS to go down hill
    IStart=IGeo
    DO iGEO=IStart,C%Opts%NSteps
       ! Compute new gradients
       CALL Force(iBAS,iGEO,C%Nams,C%Opts,C%Stat,C%Geos,C%Sets,C%MPIs)
       ! Convergence statistics for the gradient
       GradMax%D(iGEO)=Zero
       GradRMS%D(iGEO)=Zero
       DO iCLONE=1,C%Geos%Clones
          GradMAX%D(iGEO)=MAX(GradMax%D(iGEO),C%Geos%Clone(iCLONE)%GradMax)
          GradRMS%D(iGEO)=MAX(GradRMS%D(iGEO),C%Geos%Clone(iCLONE)%GradRMS)
       ENDDO
       Mssg='Gmax='//TRIM(DblToShrtChar(GradMAX%D(iGEO)))//', Grms='//TRIM(DblToShrtChar(GradRMS%D(iGEO)))
       ! Go downhill by following the gradient or with GDIIS
       IF(iGEO>Relaxations)THEN
          ! GDIIS extrapolation 
          ! Put the geometries to HDF
            DO iCLONE=1,C%Geos%Clones
              C%Geos%Clone(iCLONE)%Displ%D=C%Geos%Clone(iCLONE)%AbCarts%D
            ENDDO
            CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)    
          CALL ForceDStep(iGEO,C%Nams,C%Geos,DIISErr)
            CALL New(OldXYZ,(/3,NatmsLoc/))
            DO iCLONE=1,C%Geos%Clones
              OldXYZ%D=C%Geos%Clone(iCLONE)%Displ%D
              C%Geos%Clone(iCLONE)%Displ%D=C%Geos%Clone(iCLONE)%AbCarts%D
              C%Geos%Clone(iCLONE)%AbCarts%D=OldXYZ%D
            ENDDO
            CALL Delete(OldXYZ)
            CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)    
            DO iCLONE=1,C%Geos%Clones
              C%Geos%Clone(iCLONE)%AbCarts%D=C%Geos%Clone(iCLONE)%Displ%D
            ENDDO
            CALL GeomArchive(iBAS,iGEO+1,C%Nams,C%Sets,C%Geos)    
          Mssg=TRIM(Mssg)//', Ediis='//TRIM(DblToShrtChar(DIISErr))
          ! Check for absolute convergence
          Converged=.FALSE.
          IF(GradMAX%D(iGEO)<GTol(AccL).AND.GradRMS%D(iGEO)<GTol(AccL))THEN
             Converged=.TRUE.
             Mssg=ProcessName('GDicer',' Converged!')//TRIM(Mssg)     
          ELSE
             Mssg=ProcessName('GDicer',' GDIIS # '//TRIM(IntToChar(iGEO)))//TRIM(Mssg)          
          ENDIF
          CALL OpenASCII(C%Nams%OFile,Out)
          WRITE(Out,*)TRIM(Mssg)
          CLOSE(Out,STATUS='KEEP')
          IF(Converged)EXIT
       ELSE
          ! Take a few small steps along the gradient to start
          DO iCLONE=1,C%Geos%Clones
             C%Geos%Clone(iCLONE)%Displ%D=C%Geos%Clone(iCLONE)%AbCarts%D &
                                      -StepLength*C%Geos%Clone(iCLONE)%Vects%D
          ENDDO
            CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)    
          DO iCLONE=1,C%Geos%Clones
             C%Geos%Clone(iCLONE)%AbCarts%D=C%Geos%Clone(iCLONE)%Displ%D
          ENDDO
            CALL GeomArchive(iBAS,iGEO+1,C%Nams,C%Sets,C%Geos)    
          ! Put the geometries to HDF
          ! And here is some 
          Mssg=TRIM(Mssg)//', Step='//TRIM(DblToShrtChar(StepLength))
          Mssg=ProcessName('GDicer',' Relax # '//TRIM(IntToChar(iGEO)))//TRIM(Mssg)
          CALL OpenASCII(C%Nams%OFile,Out)
          WRITE(Out,*)TRIM(Mssg)
          WRITE(*,*)TRIM(Mssg)
          CLOSE(Out,STATUS='KEEP')
       ENDIF
       ! Compute a new energy
       CALL SCF(iBAS,iGEO+1,C)
    ENDDO
    ! Clean up
    CALL Delete(GradMax)
    CALL Delete(GradRMS)
  END SUBROUTINE GDicer
  !=====================================================================================
  ! EXTRAPOLATED GEOMETRY STEP USEING GDIIS WITH CARTESIAN FORCE ERRORS
  !=====================================================================================  
  SUBROUTINE ForceDStep(cGEO,N,G,DIISError)
    TYPE(FileNames)        :: N
    TYPE(Geometries)       :: G
    TYPE(BasisSets)        :: B
    TYPE(DBL_RNK2)         :: A,AInvM,Carts
    TYPE(DBL_VECT)         :: V,DIISCo,GradI,GradJ
    REAL(DOUBLE)           :: DIISError,CoNo,Ratio,AMax
    INTEGER                :: cGEO,iCLONE,iGEO,mGEO,iDIIS,jDIIS,nDIIS,iATS,I,J,K,Info
    INTEGER,     PARAMETER :: MaxCoef=8
    REAL(DOUBLE),PARAMETER :: MaxCoNo=1D5
    !----------------------------------------------------------------------------------!
    ! Initial dimension of the A matrix
    nDIIS=MIN(cGEO+1,MaxCoef+1)
    ! DIIS extrapolation for each clone
    HDFFileID=OpenHDF(N%HFile)
    DO iCLONE=1,G%Clones
       HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
       CALL New(GradI,G%Clone(iCLONE)%NAtms*3)
       CALL New(GradJ,G%Clone(iCLONE)%NAtms*3)
1      DIISError=Zero
       ! How far back we will go 
       mGEO=MAX(1,cGEO-nDIIS+2)
       ! Allocate work space
       CALL SetDSYEVWork(nDIIS**2)
       CALL New(A,(/nDIIS,nDIIS/))
       CALL New(DIISCo,nDIIS)
       CALL New(AInvM,(/nDIIS,nDIIS/))
       CALL New(V,nDIIS)
       DO I=1,nDIIS-1
          iDIIS=mGEO+I-1
          CALL Get(GradI,'grade',Tag_O=IntToChar(iDIIS))
          jDIIS=iDIIS
          DO J=I,nDIIS-1
             CALL Get(GradJ,'grade',Tag_O=IntToChar(jDIIS))
             ! The A matrix
             A%D(I,J)=DOT_PRODUCT(GradI%D,GradJ%D)
             A%D(J,I)=A%D(I,J)
             jDIIS=jDIIS+1
          ENDDO
       ENDDO
       A%D(nDIIS,1:nDIIS)=One
       A%D(1:nDIIS,nDIIS)=One
       A%D(nDIIS,nDIIS)=Zero
       ! Set up the eigensystem
       V%D=Zero
       V%D(nDIIS)=One
       AInvM%D=Zero
       ! Solve the inverse least squares problem
       CALL FunkOnSqMat(nDIIS,Inverse,A%D,AInvM%D,CoNo_O=CoNo,Unit_O=6)
       CALL DGEMV('N',nDIIS,nDIIS,One,AInvM%D,nDIIS,V%D,1,Zero,DIISCo%D,1)
       ! Here is the DIIS error (resdiual)       
       DIISError=MAX(DIISError,ABS(DIISCo%D(nDIIS)))
       IF(CoNo>MaxCoNo)THEN
          ! Clean up ...
          CALL Delete(V)
          CALL Delete(A)
          CALL Delete(AInvM)
          CALL Delete(DIISCo)
          CALL UnSetDSYEVWork()
          ! ... and start over with a smaller subspace
          nDIIS=nDIIS-1
          ! Sometimes ya just gota have a goto
          GOTO 1
       ELSEIF(nDIIS>2)THEN
          ! Extrapolate
          iDIIS=1
          G%Clone(iCLONE)%AbCarts%D=Zero
          CALL New(Carts,(/3,G%Clone(iCLONE)%NAtms/))
          ! DIIS extrapolation of the position ...
          DO iGEO=mGEO,cGEO
             ! ... using the unwrapped Cartesian coordinates
             CALL Get(Carts,'abcartesians',Tag_O=IntToChar(iGEO))
             G%Clone(iCLONE)%AbCarts%D=G%Clone(iCLONE)%AbCarts%D+DIISCo%D(iDIIS)*Carts%D
             iDIIS=iDIIS+1
          ENDDO
          CALL Delete(Carts)
       ELSE
          ! Stalled, freshen things up with a steepest descent move
          G%Clone(iCLONE)%AbCarts%D=G%Clone(iCLONE)%AbCarts%D-5D-1*G%Clone(iCLONE)%Vects%D          
       ENDIF
       ! Clean up a bit 
       CALL Delete(GradI)
       CALL Delete(GradJ)
       CALL Delete(V)
       CALL Delete(A)
       CALL Delete(AInvM)
       CALL Delete(DIISCo)
       CALL UnSetDSYEVWork()
    ENDDO
  END SUBROUTINE ForceDStep
  !=====================================================================================
  !
  !=====================================================================================  
  FUNCTION SteepStep(cBAS,cGEO,Energies,C) RESULT(Converged)
    TYPE(Controls)                           :: C
    TYPE(DBL_RNK2), DIMENSION(C%Geos%Clones) :: Carts
    REAL(DOUBLE),   DIMENSION(C%Geos%Clones) :: Energies
    REAL(DOUBLE)                             :: StepLength,RelErrE,MAXGrad,RMSGrad, &
         ETest,GTest
    INTEGER                                  :: cBAS,cGEO,iSTEP,iCLONE,iATS,AL,K
    INTEGER, PARAMETER                       :: MaxSTEP=8
    LOGICAL                                  :: Converged,ECnvrgd,XCnvrgd,GCnvrgd
    CHARACTER(LEN=3)                         :: chGEO
    !----------------------------------------------------------------------------------!
    chGEO=IntToChar(cGEO)
    AL=C%Opts%AccuracyLevels(cBAS)
    ! Store the current minimum energies ...
    Energies(:)=C%Geos%Clone(:)%ETotal
    ! Temp Carts variable to allow backtracking of wrapped coordinates
    RMSGrad=Zero
    MAXGrad=Zero
    DO iCLONE=1,C%Geos%Clones
       CALL New(Carts(iCLONE),(/3,C%Geos%Clone(iCLONE)%NAtms/))
       Carts(iCLONE)%D=C%Geos%Clone(iCLONE)%AbCarts%D
       MAXGrad=MAX(MAXGrad,C%Geos%Clone(iCLONE)%GradMax)
       RMSGrad=MAX(RMSGrad,C%Geos%Clone(iCLONE)%GradRMS)
    ENDDO
    ! Take some steps, more conservative if we are doing NEB ...
    IF(C%Opts%Grad==GRAD_TS_SEARCH_NEB)THEN
       StepLength=75D-2       
       ! Take a step, any step
       DO iCLONE=1,C%Geos%Clones
          C%Geos%Clone(iCLONE)%AbCarts%D=Carts(iCLONE)%D-StepLength*C%Geos%Clone(iCLONE)%Vects%D
       ENDDO
       ! Archive geometries
       CALL GeomArchive(cBAS,cGEO+1,C%Nams,C%Sets,C%Geos)    
       ! Evaluate energies at the new geometry
       CALL SCF(cBAS,cGEO+1,C)          
       ! Relative change in the total Energy
       RelErrE=-1D10
       DO iCLONE=1,C%Geos%Clones
          RelErrE=MAX(RelErrE,(Energies(iCLONE)-C%Geos%Clone(iCLONE)%ETotal) &
               /C%Geos%Clone(iCLONE)%ETotal)
       ENDDO
       ! Gradients only convergence criteria
       GCnvrgd=RMSGrad<GTol(AL).AND.MaxGrad<GTol(AL)       
       IF(GCnvrgd)THEN
          ! Cool, we are done
          Converged=.TRUE.   
          Mssg=ProcessName('SteepStep','Converged #'//TRIM(chGEO))
       ELSE
          Converged=.FALSE.  
          Mssg=ProcessName('SteepStep','Descent #'//TRIM(chGEO))
       ENDIF
    ELSE
       ! Take some steps 
       StepLength=2D0
       DO iSTEP=1,MaxSTEP
          StepLength=StepLength/Two
          ! Step the absolute positions
          DO iCLONE=1,C%Geos%Clones
             C%Geos%Clone(iCLONE)%AbCarts%D=Carts(iCLONE)%D-StepLength*C%Geos%Clone(iCLONE)%Vects%D
          ENDDO
          ! Archive geometries
          CALL GeomArchive(cBAS,cGEO+1,C%Nams,C%Sets,C%Geos)    
          ! Evaluate energies at the new geometry
          CALL SCF(cBAS,cGEO+1,C)          
          ! Relative change in the total Energy
          RelErrE=-1D10
          DO iCLONE=1,C%Geos%Clones
             RelErrE=MAX(RelErrE,(Energies(iCLONE)-C%Geos%Clone(iCLONE)%ETotal) &
                  /C%Geos%Clone(iCLONE)%ETotal)
          ENDDO
          ! Check for going downhill, convergence or stall  
          IF(MaxStep==1)THEN
             ECnvrgd=ABS(RelErrE)<1D1*ETol(AL)
          ELSE
             ECnvrgd=.FALSE.
          ENDIF
          GCnvrgd=RMSGrad<GTol(AL).AND.MaxGrad<GTol(AL)
          !       IF(ECnvrgd.OR.(GCnvrgd.AND.XCnvrgd))THEN
          IF(ECnvrgd.OR.GCnvrgd)THEN
             ! Cool, we are done
             Converged=.TRUE.   
             Mssg=ProcessName('SteepStep','Converged #'//TRIM(chGEO))
             EXIT
          ELSEIF(RelErrE<Zero)THEN
             ! Went down hill but not converged
             Converged=.FALSE.  
             Mssg=ProcessName('SteepStep','Descent #'//TRIM(chGEO))
             EXIT
          ELSEIF(iSTEP==MaxSTEP)THEN
             ! Probably need to readjust thresholds/accuracy goals 
             CALL MondoHalt(DRIV_ERROR,' Reached max resolution in energy = '   &
                  //TRIM(DblToShrtChar(RelErrE))//' in SteepStep.')          
          ELSE
             ! Need to shorten the step length.  Try again ...
             Mssg=ProcessName('SteepStep','BkTrack #'//TRIM(chGEO))          
             Mssg=TRIM(Mssg)//' dE= '//TRIM(DblToShrtChar(RelErrE))  &
                  //', Grms= '//TRIM(DblToShrtChar(RMSGrad))         &
                  //', Gmax= '//TRIM(DblToShrtChar(MAXGrad))         &
                  //', Step= '//TRIM(DblToShrtChar(StepLength))
             !    WRITE(*,*)TRIM(Mssg)             
             CALL OpenASCII(OutFile,Out)
             WRITE(*,*)TRIM(Mssg)             
             WRITE(Out,*)TRIM(Mssg)             
             CLOSE(Out)
          ENDIF
       ENDDO
    ENDIF
    Mssg=TRIM(Mssg)//' dE= '//TRIM(DblToShrtChar(RelErrE)) &
         //', Grms= '//TRIM(DblToShrtChar(RMSGrad))        & 
         //', Gmax= '//TRIM(DblToShrtChar(MAXGrad))        &
         //', Step= '//TRIM(DblToShrtChar(StepLength))
    !    WRITE(*,*)TRIM(Mssg)             
    CALL OpenASCII(OutFile,Out)
    WRITE(*,*)TRIM(Mssg)             
    WRITE(Out,*)TRIM(Mssg)             
    CLOSE(Out)
    ! Clean up
    DO iCLONE=1,C%Geos%Clones
       CALL Delete(Carts(iCLONE))
    ENDDO
    ! Keep geometry
    CALL GeomArchive(cBAS,cGEO,C%Nams,C%Sets,C%Geos)    
    ! Store the current minimum energies ...
    Energies(:)=C%Geos%Clone(:)%ETotal
  END FUNCTION SteepStep
!
!------------------------------------------------------------------
!
   SUBROUTINE IntOpt(C)
     TYPE(Controls)            :: C
     INTEGER                   :: iBAS,iGEO,iCLONE
     INTEGER                   :: AccL 
     INTEGER                   :: FirstGeom,NatmsLoc
     INTEGER                   :: ConvgdAll,MaxSteps,IStart
     TYPE(INT_VECT)            :: Convgd
     !
     ! Set geometry optimization controls
     CALL SetGeOpCtrl(C%GOpt,C%Geos,C%Opts,C%Sets,C%Nams)
     ! initial geometry
     iGEO=C%Stat%Previous%I(3)
     MaxSteps=C%GOpt%GConvCrit%MaxGeOpSteps
     ! Build the guess 
     DO iBAS=1,C%Sets%NBSets-1
       CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)    
       CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
       CALL SCF(iBAS,iGEO,C)
     ENDDO
     iBAS=C%Sets%NBSets
     ! Print the starting coordinates and energy
     DO iCLONE=1,C%Geos%Clones
       CALL PPrint(C%Geos%Clone(iCLONE),C%Nams%GFile, &
         Geo,C%Opts%GeomPrint)
     ENDDO
     !
     CALL New(Convgd,C%Geos%Clones)
     Convgd%I=0
     !
     ! Start optimization                     
     !
     IStart=iGEO
     DO iGEO=IStart,MaxSteps
       CALL PrintClones(IGeo,C%Nams,C%Geos)
       CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)    
       CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
       !
       ! Calculate energy and force for all clones at once.
       !
       CALL SCF(iBAS,iGEO,C)
       IF(iGEO>1) CALL BackTrack(iBAS,iGEO,C)
       CALL Force(iBAS,iGEO,C%Nams,C%Opts,C%Stat, &
                  C%Geos,C%Sets,C%MPIs)
       !
       ! Loop over all clones and modify geometries.
       !
       ConvgdAll=1
       DO iCLONE=1,C%Geos%Clones
         Convgd%I=0
         CALL OptSingleMol(C%GOpt,C%Nams,C%Opts,C%Sets,C%Geos, &
           C%Geos%Clone(iCLONE),Convgd%I,iGEO,iCLONE)
         ConvgdAll=ConvgdAll*Convgd%I(iCLONE)
       ENDDO 
       !
       ! Archive displaced geometries 
       !
       CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)    
       !
       ! Fill in new geometries 
       !
       DO iCLONE=1,C%Geos%Clones
         CALL NewGeomFill(C%Geos%Clone(iCLONE))
       ENDDO
       !
       ! Do GDIIS and print geometries
       !
       DO iCLONE=1,C%Geos%Clones
         CALL MixGeoms(C%Opts,C%Nams,C%GOpt,Convgd%I, &
                       C%Geos%Clone(iCLONE),iCLONE,iGEO)
         CALL PPrint(C%Geos%Clone(iCLONE),C%Nams%GFile,Geo, &
                     C%Opts%GeomPrint)
       ENDDO
       !
       C%Stat%Previous%I(3)=IGeo
       C%Stat%Current%I(3)=IGeo+1
       !
       ! Continue optimization?
       !
       IF(ConvgdAll==1) EXIT
     ENDDO
     CALL Delete(Convgd)
     !
     IGeo=C%Stat%Current%I(3)
     IF(IGeo>=MaxSteps) THEN
       CALL OpenASCII(OutFile,Out)
       IF(ConvgdAll/=1) THEN
         WRITE(Out,700) 
         WRITE(*,700) 
       ELSE
         WRITE(Out,460) iGeo-1
         WRITE(*,460) iGeo-1
       ENDIF          
       CLOSE(Out,STATUS='KEEP')
     ENDIF
     !
     ! Convergence is reached at this point, calculate final energy
     ! and finish optimization.
     !
     CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)    
     CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
     CALL SCF(iBAS,iGEO,C)
     !
     CALL OpenASCII(OutFile,Out)
     WRITE(Out,500)  
     WRITE(*,500)  
     WRITE(Out,600)  
     WRITE(*,600)  
     DO iCLONE=1,C%Geos%Clones
       WRITE(Out,400) iCLONE,C%Geos%Clone(iCLONE)%ETotal
       WRITE(*,400) iCLONE,C%Geos%Clone(iCLONE)%ETotal
     ENDDO
     CLOSE(Out,STATUS='KEEP')
     !
     400  FORMAT(I6,F20.10)
     500  FORMAT('Energies of final structures:')
     600  FORMAT(' Clone #         Energy         ')
     460  FORMAT('Geometry Optimization converged in ',I6,' steps.')
     700  FORMAT('Maximum number of optimization'&
            //' steps exceeded, optimization did not converge.')
   END SUBROUTINE IntOpt
!
!------------------------------------------------------------------
!
   SUBROUTINE ModifyGeom(GOpt,XYZ,AtNum,GradIn,LagrMult,GradMult, &
                  LagrDispl,Convgd,ETot,ELagr,iGEO,iCLONE,SCRPath,Print)
     TYPE(GeomOpt)               :: GOpt
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ,GradIn
     REAL(DOUBLE),DIMENSION(:)   :: AtNum,LagrMult,GradMult,LagrDispl
     REAL(DOUBLE)                :: ETot,ELagr
     INTEGER,DIMENSION(:)        :: Convgd
     INTEGER                     :: NIntC,NCart
     INTEGER                     :: NatmsLoc,iGEO,iCLONE
     TYPE(INTC)                  :: IntCs
     TYPE(DBL_VECT)              :: IntOld
     INTEGER                     :: Refresh
     CHARACTER(LEN=*)            :: SCRPath
     INTEGER                     :: Print
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     !
     ! Do we have to refresh internal coord defs?   
     !
     CALL IntCReDef(GOpt,Refresh)
     !
     ! Get internal coord defs.
     !
     IF(Refresh/=0) THEN
       CALL GetIntCs(XYZ,AtNum, &
         IntCs,NIntC,Refresh,SCRPath,GOpt%CoordCtrl,GOpt%Constr)
       IF(NIntC==0) CALL Halt('Molecule has dissociated,'// &
                     'optimizer did not find any internal coordinates.')
     ENDIF
     IF(NIntC/=0) THEN
       CALL INTCValue(IntCs,XYZ,GOpt%CoordCtrl%LinCrit)
       CALL New(IntOld,NIntC)
       IntOld%D=IntCs%Value
     ENDIF
     !
     ! Print current set of internals for debugging
     !
     IF(Print==DEBUG_GEOP_MAX) CALL PrtIntCoords(IntCs, &
       IntCs%Value,'Internals at step #'//TRIM(IntToChar(iGEO)))
     !
     ! Calculate the B matrix and put it onto the disk.
     !
     CALL RefreshBMatInfo(IntCs,XYZ,GOpt%TrfCtrl%DoClssTrf,Print, &
                     GOpt%CoordCtrl%LinCrit,GOpt%TrfCtrl%ThreeAt, &
                     SCRPath)
     !
     ! Overwrite Cartesian Energy gradients with Lagrangian ones.
     !
     IF(GOpt%Constr%DoLagr) THEN
       CALL LagrGradMult(GradMult,XYZ,IntCs, &
         GOpt%CoordCtrl%LinCrit,GOpt%Constr%NConstr)
       CALL CalcLagrMult(GradMult,LagrMult,LagrDispl,XYZ,IntCs, &
                         GradIn,GOpt%GrdTrf,GOpt%CoordCtrl, &
                         GOpt%TrfCtrl,GOpt%Constr,GOpt%Hessian, &
                         Print,SCRPath)
       CALL LagrGradCart(GradIn,IntCs,LagrDispl,SCRPath,iGEO)
     ENDIF
     CALL LagrEnergy(LagrDispl,IntCs,ELagr)
     !
     ! Calculate simple relaxation step from an inverse Hessian
     !
     CALL RelaxGeom(GOpt,XYZ,AtNum,GradIn,GradMult, &
                    LagrMult,LagrDispl,IntCs,SCRPath,Print) 
     !
     ! Check convergence
     !
     CALL LagrConv(GOpt%GOptStat,LagrDispl,LagrMult,GradMult)
     CALL GeOpConv(GOpt%Constr,GOpt%GOptStat,GOpt%CoordCtrl, &
                   GOpt%GConvCrit,XYZ,ETot,ELagr,IntCs,IntOld, &
                   iCLONE,iGEO)
     !
     IF(GOpt%GOptStat%GeOpConvgd) THEN
       Convgd(iCLONE)=1
     ENDIF
     !
     ! tidy up
     !
     CALL Delete(IntOld)
     CALL Delete(IntCs)
   END SUBROUTINE ModifyGeom
!
!--------------------------------------------------------------------
!
   SUBROUTINE RelaxGeom(GOpt,XYZ,AtNum,GradIn,GradMult, &
                        LagrMult,LagrDispl,IntCs,SCRPath,Print)
     !
     ! Simple Relaxation step
     !
     TYPE(GeomOpt)                  :: GOpt
     REAL(DOUBLE),DIMENSION(:,:)    :: XYZ,GradIn
     REAL(DOUBLE),DIMENSION(:)      :: AtNum,GradMult,LagrMult,LagrDispl
     TYPE(DBL_VECT)                 :: Displ
     TYPE(DBL_VECT)                 :: IntGrad,Grad,CartGrad
     TYPE(INTC)                     :: IntCs
     INTEGER                        :: I,J,NDim
     INTEGER                        :: NatmsLoc,NCart,NIntC,Print
     CHARACTER(LEN=*)               :: SCRPath 
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     NIntC=SIZE(IntCs%Def)
     !
     CALL NewDispl(GOpt,Displ,NCart,NIntC)
     !
     NDim =SIZE(Displ%D)
     IF(NIntC/=NDim.AND.GOpt%TrfCtrl%DoInternals) &
         CALL Halt('Dimensionality error in RelaxGeom')
     !
     CALL New(Grad,NDim)
     CALL New(CartGrad,NCart)
     CALL CartRNK2ToCartRNK1(CartGrad%D,GradIn)
     !
     CALL GetCGradMax(GOpt%Constr,CartGrad%D,NCart,  &
                      IntCs,GOpt%GOptStat%IMaxCGrad, &
                      GOpt%GOptStat%MaxCGrad)
     !
     IF(GOpt%TrfCtrl%DoInternals) THEN
       ! calculate internal coord gradients
       CALL New(IntGrad,NDim)
       CALL CartToInternal(XYZ,IntCs,CartGrad%D,IntGrad%D,&
         GOpt%GrdTrf,GOpt%CoordCtrl,GOpt%TrfCtrl,Print,SCRPath)
       Grad%D=IntGrad%D
       CALL Delete(IntGrad)
     ELSE
       Grad%D=CartGrad%D
     ENDIF
     CALL Delete(CartGrad)
     !
     ! Check for gradient-convergence
     !
     CALL GrdConvrgd(GOpt%GOptStat,IntCs,Grad%D)
     !
     ! Use Hessian matrix to calculate displacement
     !
     SELECT CASE(GOpt%Optimizer)
       CASE(GRAD_STPDESC_OPT) 
         CALL SteepestDesc(GOpt%CoordCtrl,GOpt%Hessian, &
                           Grad,Displ,XYZ)
       CASE(GRAD_DIAGHESS_OPT) 
        !IF(GOpt%Constr%DoLagr.AND.GOpt%Constr%NConstr/=0) THEN
        !  CALL DiagHessLagr(GOpt%CoordCtrl,GOpt%Hessian,Grad%D,Displ%D,&
        !             IntCs,XYZ,SCRPath,LagrMult,GradMult,LagrDispl)
        !ELSE
           CALL DiagHess(GOpt%CoordCtrl,GOpt%Hessian,Grad,Displ, &
                         IntCs,AtNum,XYZ)
        !ENDIF
     END SELECT
     !
     ! Set constraints on the displacements
     !
     IF(.NOT.GOpt%Constr%DoLagr) THEN
       CALL SetConstraint(IntCs,XYZ,Displ,GOpt%CoordCtrl%LinCrit, &
         GOpt%Constr%NConstr,GOpt%TrfCtrl%DoInternals)
     ENDIF
     CALL Delete(Grad)
     !
     ! Add Cartesian or internal coord. displacements 
     ! to current structure
     !
     IF(GOpt%TrfCtrl%DoInternals) THEN 
       CALL INTCValue(IntCs,XYZ,GOpt%CoordCtrl%LinCrit)
       IntCs%Value=IntCs%Value+Displ%D
       CALL InternalToCart(XYZ,IntCs,IntCs%Value,Print, &
                           GOpt%BackTrf,GOpt%TrfCtrl,GOpt%CoordCtrl,&
                           GOpt%Constr,SCRPath)
     ELSE
       CALL CartRNK1ToCartRNK2(Displ%D,XYZ,.TRUE.)
     ENDIF
     !
     CALL Delete(Displ)
   END SUBROUTINE RelaxGeom
!
!-------------------------------------------------------
!
   SUBROUTINE NewStructure(Print,GBackTrf,GTrfCtrl,GCoordCtrl, &
     GConstr,SCRPath,XYZ,Displ,IntCs)
     INTEGER                        :: Print
     TYPE(BackTrf)                  :: GBackTrf
     TYPE(TrfCtrl)                  :: GTrfCtrl
     TYPE(CoordCtrl)                :: GCoordCtrl
     TYPE(Constr)                   :: GConstr
     CHARACTER(LEN=*)               :: SCRPath
     REAL(DOUBLE),DIMENSION(:,:)    :: XYZ   
     TYPE(INTC)                     :: IntCs
     INTEGER                        :: I,J,II,NDim,NIntc
     INTEGER                        :: NatmsLoc,NCart,InitGDIIS
     TYPE(DBL_VECT)                 :: Displ
     LOGICAL                        :: DoInternals
     REAL(DOUBLE)                   :: LinCrit
     !
     ! In the present version there is no line search, only GDIIS
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc   
     NDim =SIZE(Displ%D)
     NIntC =SIZE(IntCs%Def)
     LinCrit=GCoordCtrl%LinCrit
     DoInternals=GTrfCtrl%DoInternals
     IF(NIntC/=NDim.AND.DoInternals) &
         CALL Halt('Dimensionality error in NewStructure.')
     !
     ! Construct new structure either by Cartesian or by internal
     ! displacements
     !
     IF(DoInternals) THEN 
       CALL INTCValue(IntCs,XYZ,LinCrit)
       IntCs%Value=IntCs%Value+Displ%D
       CALL InternalToCart(XYZ,IntCs,IntCs%Value, &
         Print,GBackTrf,GTrfCtrl,GCoordCtrl,GConstr,SCRPath)
     ELSE
       CALL CartRNK1ToCartRNK2(Displ%D,XYZ,.TRUE.)
     ENDIF
     !
   END SUBROUTINE NewStructure
!
!-----------------------------------------------------------------
!
   SUBROUTINE IntCReDef(GOpt,Refresh)
     TYPE(GeomOpt)     :: GOpt
     INTEGER           :: Refresh
     !
     !WARNING! refresh may change the number of internal coordinates!
     !
       Refresh=1
     !IF(GOpt%GOptStat%ActStep==1) THEN
     !  Refresh=1
     !  IF(GOpt%CoordCtrl%RefreshIn==4) Refresh=4
     !ELSE
     !  Refresh=GOpt%CoordCtrl%RefreshIn
     !ENDIF
       GOpt%CoordCtrl%Refresh=Refresh
   END SUBROUTINE IntCReDef
!
!------------------------------------------------------------------
!
   SUBROUTINE NewDispl(GOpt,Displ,NCart,NIntC)
     TYPE(GeomOpt)  :: GOpt 
     TYPE(DBL_VECT) :: Displ
     INTEGER        :: NCart,NIntC
     !
     IF(GOpt%CoordCtrl%CoordType==CoordType_Cartesian) THEN
       CALL New(Displ,NCart)
       Displ%D(:)=Zero
     ELSE
       CALL New(Displ,NIntC)
       Displ%D(:)=Zero
     ENDIF
   END SUBROUTINE NewDispl
!
!-------------------------------------------------------
!
   SUBROUTINE SteepestDesc(CoordC,Hess,Grad,Displ,XYZ)
     TYPE(CoordCtrl)             :: CoordC
     TYPE(Hessian)               :: Hess  
     TYPE(DBL_VECT)              :: Grad,Displ 
     INTEGER                     :: NCart
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     !
     NCart=3*SIZE(XYZ,2)
     IF(CoordC%CoordType==CoordType_Cartesian) THEN
       Displ%D=-Hess%StpDescInvH*Grad%D 
       !!! take translations and rotations off
       CALL TranslsOff(Displ%D,.FALSE.)
       CALL RotationsOff(Displ%D,XYZ,.FALSE.)
     ELSE IF(CoordC%CoordType==CoordType_PrimInt) THEN
       Displ%D=-2.D0*Grad%D 
       !CALL RedundancyOff(Displ%D,XYZ,DoSet_O=.TRUE.)
     ELSE
       Displ%D=-0.5D0*Grad%D 
     ENDIF 
   END SUBROUTINE SteepestDesc
!
!-------------------------------------------------------
!
   SUBROUTINE DiagHess(CoordC,Hess,Grad,Displ,IntCs,AtNum,XYZ)
     TYPE(CoordCtrl)   :: CoordC
     TYPE(Hessian)     :: Hess
     TYPE(DBL_VECT)    :: Grad,Displ 
     TYPE(INTC)        :: IntCs       
     INTEGER           :: I,J,NIntC,NCart,NConstr,NatmsLoc
     REAL(DOUBLE)      :: HStre,HBend,HLinB,HOutP,HTors
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(:)   :: AtNum
     TYPE(DBL_VECT)    :: DHess
     !
     NIntC=SIZE(IntCs%Def)
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     CALL New(DHess,NIntC)
     !
     CALL DiagHess_Vals(DHess%D,AtNum,XYZ,IntCs,CoordC, &
                        Hess,CoordC%NStre,'ThreeVals')
    !CALL DiagHess_Vals(DHess%D,AtNum,XYZ,IntCs,CoordC, &
    !                   Hess,CoordC%NStre,'Lindh')
     !
     IF(CoordC%CoordType==CoordType_Cartesian) THEN
       Displ%D=-1.D0*Grad%D !!!! equivalent with stpdesc
     ELSE IF(CoordC%CoordType==CoordType_PrimInt) THEN
       DO I=1,NIntC
         Displ%D(I)=-DHess%D(I)*Grad%D(I)
       ENDDO
       ! CALL RedundancyOff(Displ%D,XYZ,DoSet_O=.TRUE.)
     ELSE
       CALL Halt('Only Primitiv Internals are available yet.')
     ENDIF 
     ! 
     CALL Delete(DHess)
   END SUBROUTINE DiagHess
!
!-------------------------------------------------------
!
   SUBROUTINE DiagHessRFO(GOpt,Grad,Displ,IntCs,XYZ)
     TYPE(GeomOpt)               :: GOpt
     TYPE(DBL_VECT)              :: Grad,Displ
     TYPE(DBL_RNK2)              :: Hessian
     TYPE(INTC)                  :: IntCs
     INTEGER                     :: I,J,NIntC,NCart,Info
     REAL(DOUBLE)                :: HStre,HBend,HLinB,HOutP,HTors
     REAL(DOUBLE)                :: Sum
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     !
     NIntC=SIZE(IntCs%Def)
     NCart=3*SIZE(XYZ,2)
     CALL New(Hessian,(/NIntC+1,NIntC+1/))
     Hessian%D=Zero
     !
     IF(GOpt%CoordCtrl%CoordType==CoordType_Cartesian) THEN
       Displ%D=-1.D0*Grad%D !!!! equivalent with StpDesc
     ELSE IF(GOpt%CoordCtrl%CoordType==CoordType_PrimInt) THEN
       HStre=GOpt%Hessian%Stre
       HBend=GOpt%Hessian%Bend
       HLinB=GOpt%Hessian%LinB
       HOutP=GOpt%Hessian%OutP
       HTors=GOpt%Hessian%Tors
       DO I=1,NIntC
         IF(IntCs%Def(I)(1:4)=='STRE') THEN
           Hessian%D(I,I)=HStre
         ELSE IF(IntCs%Def(I)(1:4)=='BEND') THEN
           Hessian%D(I,I)=HBend
         ELSE IF(IntCs%Def(I)(1:4)=='LINB') THEN
           Hessian%D(I,I)=HLinB
         ELSE IF(IntCs%Def(I)(1:4)=='OUTP') THEN
           Hessian%D(I,I)=HOutP
         ELSE IF(IntCs%Def(I)(1:4)=='TORS') THEN
           Hessian%D(I,I)=HTors
         ENDIF
       ENDDO
       Hessian%D(1:NIntC,NIntC+1)=Grad%D
       Hessian%D(NIntC+1,1:NIntC)=Grad%D
       !
       ! diagonalize RFO matrix 
       !
       CALL SetDSYEVWork(NIntC+1)
       !
       BLKVECT%D=Hessian%D
       CALL DSYEV('V','U',NIntC+1,BLKVECT%D,BIGBLOK,BLKVALS%D, &
         BLKWORK%D,BLKLWORK,INFO)
       IF(INFO/=SUCCEED) &
       CALL Halt('DSYEV hosed in RotationsOff. INFO='&
                  //TRIM(IntToChar(INFO)))
       !
       ! Choose the eigenvector of the lowest eigenvalue for step,
       ! after RFO scaling
       !
       Sum=One/BLKVECT%D(NIntC+1,1)  
       Displ%D=BLKVECT%D(:,1)*Sum
     CALL UnSetDSYEVWork()
     !
     !project out redundancy
     !         CALL RedundancyOff(Displ%D,XYZ,DoSet_O=.TRUE.)
     ELSE
       CALL Halt('Only Primitiv Internals are available yet.')
     ENDIF 
     CALL Delete(Hessian)
   END SUBROUTINE DiagHessRFO
!
!---------------------------------------------------------------
!
   SUBROUTINE GeOpConv(CtrlConstr,CtrlStat,CtrlCoord,GConvCr, &
                       XYZ,Etot,ELagr,IntCs,IntOld,iCLONE,iGEO)
     TYPE(GOptStat)             :: CtrlStat
     TYPE(Constr)               :: CtrlConstr
     TYPE(CoordCtrl)            :: CtrlCoord  
     TYPE(GConvCrit)            :: GConvCr    
     REAL(DOUBLE),DIMENSION(:,:):: XYZ
     TYPE(INTC)                 :: IntCs
     TYPE(DBL_VECT)             :: IntOld,AuxVect
     INTEGER                    :: iCLONE,iGEO
     !
     REAL(DOUBLE)               :: MaxStreDispl,MaxBendDispl
     REAL(DOUBLE)               :: MaxLinBDispl,MaxOutPDispl
     REAL(DOUBLE)               :: MaxTorsDispl
     REAL(DOUBLE)               :: RMSIntDispl
     !
     REAL(DOUBLE)               :: StreConvCrit,BendConvCrit
     REAL(DOUBLE)               :: LinBConvCrit,OutPConvCrit
     REAL(DOUBLE)               :: TorsConvCrit
     !                      
     INTEGER                    :: I,J,K,L
     INTEGER                    :: NCart,NIntC,NatmsLoc
     REAL(DOUBLE)               :: RMSGrad,MaxGrad,MaxCGrad
     REAL(DOUBLE)               :: RMSGradNoConstr,MaxGradNoConstr
     REAL(DOUBLE)               :: Etot,ELagr,Sum
     !                      
     INTEGER                    :: NStreGeOp,NBendGeOp,NLinBGeOp
     INTEGER                    :: NOutPGeOp,NTorsGeOp
     INTEGER                    :: MaxStre,MaxBend,MaxLinB
     INTEGER                    :: MaxOutP,MaxTors
     !                      
     INTEGER                    :: IMaxGrad,IMaxCGrad,IMaxGradNoConstr
     LOGICAL                    :: GradConv
     !
     NIntC=SIZE(IntCs%Def)
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     !
     RMSGrad=CtrlStat%RMSGrad
     RMSGradNoConstr=CtrlStat%RMSGradNoConstr
     MaxGrad=CtrlStat%MaxGrad
     MaxCGrad=CtrlStat%MaxCGrad
     MaxGradNoConstr=CtrlStat%MaxGradNoConstr
     IMaxGrad=CtrlStat%IMaxGrad
     IMaxCGrad=CtrlStat%IMaxCGrad
     IMaxGradNoConstr=CtrlStat%IMaxGradNoConstr
     !
     NStreGeOp=CtrlCoord%NStre
     NBendGeOp=CtrlCoord%NBend
     NLinBGeOp=CtrlCoord%NLinB
     NOutPGeOp=CtrlCoord%NOutP
     NTorsGeOp=CtrlCoord%NTors
     !
     ! Size of internal coordinate changes
     !
     CALL INTCValue(IntCs,XYZ,CtrlCoord%LinCrit)
     IntOld%D=IntCs%Value-IntOld%D
     CALL MapAngleDispl(IntCs,NIntC,IntOld%D)
     MaxStre=0
     MaxBend=0
     MaxLinB=0
     MaxOutP=0
     MaxTors=0
     MaxStreDispl=-One
     MaxBendDispl=-One
     MaxLinBDispl=-One
     MaxOutPDispl=-One
     MaxTorsDispl=-One
     DO I=1,NIntC 
       IF(.NOT.IntCs%Active(I)) CYCLE
        IF(IntCs%Def(I)(1:4)=='STRE') THEN 
          Sum=ABS(IntOld%D(I))
          IF(MaxStreDispl<Sum) THEN
            MaxStre=I
            MaxStreDispl=Sum
          ENDIF
        ELSE IF(IntCs%Def(I)(1:4)=='BEND') THEN
          Sum=ABS(IntOld%D(I))
          IF(MaxBendDispl<Sum) THEN
            MaxBend=I
            MaxBendDispl=Sum
          ENDIF
        ELSE IF(IntCs%Def(I)(1:4)=='LINB') THEN
          Sum=ABS(IntOld%D(I))
          IF(MaxLinBDispl<Sum) THEN
            MaxLinB=I
            MaxLinBDispl=Sum
          ENDIF
        ELSE IF(IntCs%Def(I)(1:4)=='OUTP') THEN
          Sum=ABS(IntOld%D(I))
          IF(MaxOutPDispl<Sum) THEN
            MaxOutP=I
            MaxOutPDispl=Sum
          ENDIF
        ELSE IF(IntCs%Def(I)(1:4)=='TORS') THEN
          Sum=ABS(IntOld%D(I))
          IF(MaxTorsDispl<Sum) THEN
            MaxTors=I
            MaxTorsDispl=Sum
          ENDIF
        ENDIF
     ENDDO
     RMSIntDispl=SQRT(DOT_PRODUCT(IntOld%D,IntOld%D)/DBLE(NIntC))
     !
     CtrlStat%MaxStreDispl=MaxStreDispl
     CtrlStat%MaxBendDispl=MaxBendDispl
     CtrlStat%MaxLinBDispl=MaxLinBDispl
     CtrlStat%MaxOutPDispl=MaxOutPDispl
     CtrlStat%MaxTorsDispl=MaxTorsDispl
     CtrlStat%RMSIntDispl =RMSIntDispl
     !
     GradConv=MaxCGrad<GConvCr%Grad.AND. &
              ((RMSGrad<GConvCr%Grad.AND.MaxGrad<GConvCr%Grad).OR. &
              (CtrlConstr%NConstr/=0.AND..NOT.CtrlConstr%DoLagr))
     CtrlStat%GeOpConvgd=GradConv.AND. &
                         CtrlStat%MaxLGrad<GConvCr%Grad.AND. &
                         MaxStreDispl<GConvCr%Stre.AND. &
                         MaxBendDispl<GConvCr%Bend.AND. &
                         MaxLinBDispl<GConvCr%LinB.AND. &
                         MaxOutPDispl<GConvCr%OutP.AND. &
                         MaxTorsDispl<GConvCr%Tors
     !
     ! Review iterations
     !
     WRITE(*,399) iCLONE,iGEO,ETot
     WRITE(Out,399) iCLONE,iGEO,ETot
     IF(CtrlConstr%NConstr/=0.AND.CtrlConstr%DoLagr) THEN
       WRITE(*,499) ELagr
       WRITE(Out,499) ELagr
     ENDIF
     !
     MaxStreDispl=MaxStreDispl/AngstromsToAu
     MaxBendDispl=MaxBendDispl*180.D0/PI
     MaxLinBDispl=MaxLinBDispl*180.D0/PI
     MaxOutPDispl=MaxOutPDispl*180.D0/PI
     MaxTorsDispl=MaxTorsDispl*180.D0/PI
     !
     WRITE(*,410) MaxGrad,IntCs%Atoms(IMaxGrad,1:4)
     WRITE(*,140) MaxCGrad,(IMaxCGrad-1)/3+1
     WRITE(*,420) RMSGrad
     WRITE(Out,410) MaxGrad,IntCs%Atoms(IMaxGrad,1:4)
     WRITE(Out,140) MaxCGrad,(IMaxCGrad-1)/3+1
     WRITE(Out,420) RMSGrad
     IF(CtrlConstr%NConstr/=0) THEN
       IF(CtrlConstr%DoLagr) THEN
         WRITE(*,930) CtrlStat%MaxLGrad,CtrlStat%IMaxLGrad
         WRITE(Out,930) CtrlStat%MaxLGrad,CtrlStat%IMaxLGrad
       ELSE IF(IMaxGradNoConstr/=0) THEN
         WRITE(*,510) MaxGradNoConstr, &
                      IntCs%Atoms(IMaxGradNoConstr,1:4)
         WRITE(*,520) RMSGradNoConstr
         WRITE(Out,510) MaxGradNoConstr, &
                        IntCs%Atoms(IMaxGradNoConstr,1:4)
         WRITE(Out,520) RMSGradNoConstr
       ENDIF
     ENDIF
     !
     IF(MaxStre/=0) THEN
       WRITE(*,430) MaxStreDispl,IntCs%Atoms(MaxStre,1:2)
       WRITE(Out,430) MaxStreDispl,IntCs%Atoms(MaxStre,1:2)
     ENDIF
     IF(MaxBend/=0) THEN
       WRITE(*,435) MaxBendDispl,IntCs%Atoms(MaxBend,1:3)
       WRITE(Out,435) MaxBendDispl,IntCs%Atoms(MaxBend,1:3)
     ENDIF
     IF(MaxLinB/=0) THEN
       WRITE(*,436) MaxLinBDispl,IntCs%Atoms(MaxLinB,1:3)
       WRITE(Out,436) MaxLinBDispl,IntCs%Atoms(MaxLinB,1:3)
     ENDIF
     IF(MaxOutP/=0) THEN
       WRITE(*,437) MaxOutPDispl,IntCs%Atoms(MaxOutP,1:4)
       WRITE(Out,437) MaxOutPDispl,IntCs%Atoms(MaxOutP,1:4)
     ENDIF
     IF(MaxTors/=0) THEN
       WRITE(*,438) MaxTorsDispl,IntCs%Atoms(MaxTors,1:4)
       WRITE(Out,438) MaxTorsDispl,IntCs%Atoms(MaxTors,1:4)
     ENDIF
     IF(CtrlConstr%NConstr/=0.AND.CtrlConstr%DoLagr) THEN
       WRITE(*,940) CtrlStat%MaxDMult
       WRITE(Out,940) CtrlStat%MaxDMult
     ENDIF
     !
     WRITE(*,440) RMSIntDispl
     WRITE(Out,440) RMSIntDispl
     !
399 FORMAT('       Clone = ',I6,' GeOp step = ',I6,' Total Energy = ',F20.8)
400 FORMAT('Total Energy at Current Geometry = ',F20.8)
401 FORMAT('                    Total Energy = ',F20.8)
499 FORMAT('               ',6X,'             ',2X,'       Lagrangian = ',F20.8)
410 FORMAT('                        Max Grad = ',F12.6,' between atoms ',4I4)
140 FORMAT('     Max Unconstrained Cart Grad = ',F12.6,'      on atom  ',4I4)
420 FORMAT('                        RMS Grad = ',F12.6)
510 FORMAT('Max Grad on Unconstrained Coords = ',F12.6,' between atoms ',4I4)
520 FORMAT('RMS Grad on Unconstrained Coords = ',F12.6)
930 FORMAT('                  Max LagrM Grad = ',F12.6,' on constr. #  ',I4)
430 FORMAT('                  Max STRE Displ = ',F12.6,' between atoms ',4I4)
435 FORMAT('                  Max BEND Displ = ',F12.6,' between atoms ',4I4)
436 FORMAT('                  Max LINB Displ = ',F12.6,' between atoms ',4I4)
437 FORMAT('                  Max OUTP Displ = ',F12.6,' between atoms ',4I4)
438 FORMAT('                  Max TORS Displ = ',F12.6,' between atoms ',4I4)
440 FORMAT('                       RMS Displ = ',F12.6)
940 FORMAT('                  Max Lagr Displ = ',F12.6)
        !
   END SUBROUTINE GeOpConv
!
!---------------------------------------------------------------
!
   SUBROUTINE SetGeOpCtrl(GOpt,Geos,Opts,Sets,Nams)
     !
     TYPE(GeomOpt)    :: GOpt
     TYPE(Options)    :: Opts
     TYPE(BasisSets)  :: Sets
     TYPE(Geometries) :: Geos
     TYPE(FileNames)  :: Nams
     INTEGER          :: NatmsLoc,NCart
     REAL(DOUBLE)     :: Sum,GCrit
     INTEGER          :: AccL
     !
     AccL    =Opts%AccuracyLevels(Sets%NBSets)
     NatmsLoc=Geos%Clone(1)%Natms
     NCart=3*NatmsLoc
     !
     CALL SetCoordCtrl(GOpt%CoordCtrl)
     CALL   SetHessian(GOpt%Hessian)
     CALL     SetStepS(GOpt%StepSize)
     CALL SetGConvCrit(GOpt%GConvCrit,GOpt%Hessian,AccL,NatmsLoc)
     CALL     SetGDIIS(GOpt%GDIIS)
     CALL    SetGrdTrf(GOpt%GrdTrf,GOpt%GConvCrit)
     CALL   SetBackTrf(GOpt%BackTrf,GOpt%GConvCrit)
     CALL    SetConstr(GOpt%Constr,GOpt%BackTrf)
     CALL   SetTrfCtrl(GOpt%TrfCtrl,GOpt%CoordCtrl,GOpt%Constr)
   END SUBROUTINE SetGeOpCtrl
!
!---------------------------------------------------------------
!
   SUBROUTINE PrintClones(IStep,Nams,Geos)
     TYPE(FileNames)    :: Nams
     TYPE(Geometries)   :: Geos
     INTEGER            :: IStep,NatmsLoc,NClones,iCLONE,I
     CHARACTER(LEN=DCL) :: FileName
     TYPE(DBL_RNK2)     :: XYZ
     TYPE(DBL_VECT)     :: AtNum
     REAL(DOUBLE)       :: Dist,Sum
     ! 
     ! The whole set of clones is going to be printed in a single file.
     ! All clones will be seen by a single view of the resulting file.
     ! 
     NClones=Geos%Clones
     NatmsLoc=0
     DO iCLONE=1,NClones
       NatmsLoc=NatmsLoc+Geos%Clone(iCLONE)%Natms
     ENDDO
     !
     CALL New(XYZ,(/3,NatmsLoc/))
     CALL New(AtNum,NatmsLoc)
     !
     NatmsLoc=0
     Sum=Zero
     Dist=7.D0
     DO iCLONE=1,NClones
       Sum=(iCLONE-1)*Dist
       DO I=1,Geos%Clone(iCLONE)%Natms
         NatmsLoc=NatmsLoc+1
         AtNum%D(NatmsLoc)=Geos%Clone(iCLONE)%AtNum%D(I)
         XYZ%D(1:2,NatmsLoc)=Geos%Clone(iCLONE)%AbCarts%D(1:2,I)
         XYZ%D(3,NatmsLoc)=Geos%Clone(iCLONE)%AbCarts%D(3,I)+Sum
       ENDDO
     ENDDO
     !
     FileName=TRIM(Nams%SCF_NAME)//'.Clones.xyz'
     !
     CALL PrtXYZ(Atnum%D,XYZ%D,FileName,&
                 'Step='//TRIM(IntToChar(IStep)))
     !
     CALL Delete(AtNum)
     CALL Delete(XYZ)
   END SUBROUTINE PrintClones
!
!-------------------------------------------------------------------
!
   SUBROUTINE OptSingleMol(GOpt,Nams,Opts,Sets,Geos, &
                           GMLoc,Convgd,iGEO,iCLONE)
     TYPE(GeomOpt)        :: GOpt
     TYPE(FileNames)      :: Nams
     TYPE(Options)        :: Opts
     TYPE(BasisSets)      :: Sets
     TYPE(Geometries)     :: Geos
     TYPE(CRDS)           :: GMLoc
     INTEGER,DIMENSION(:) :: Convgd
     INTEGER              :: I,iGEO,iCLONE,NatmsNew
     INTEGER              :: InitGDIIS,NConstr,NCart,NatmsLoc
     LOGICAL              :: NoGDIIS,GDIISOn
     CHARACTER(LEN=DCL)   :: SCRPath
     TYPE(DBL_RNK2)       :: XYZNew,GradNew
     TYPE(DBL_VECT)       :: AtNumNew,CartGrad
     !
     SCRPath  =TRIM(Nams%M_SCRATCH)//TRIM(Nams%SCF_NAME)// &
             '.'//TRIM(IntToChar(iCLONE))
     GMLoc%Displ%D=GMLoc%AbCarts%D
     !
     ! Now, cut out that part of the molecule, which is 
     ! not rigidly fixed and pass it in to the optimizer
     !
     NatmsNew=0
     DO I=1,GMLoc%Natms
       IF(GMLoc%CConstrain%I(I)/=2) NatmsNew=NatmsNew+1
     ENDDO
     CALL New(XYZNew,(/3,NatmsNew/))
     CALL New(AtNumNew,NatmsNew)
     CALL New(GradNew,(/3,NatmsNew/))
     ! fill atomic data
     NatmsNew=0
     DO I=1,GMLoc%Natms
       IF(GMLoc%CConstrain%I(I)/=2) THEN
         NatmsNew=NatmsNew+1
         XYZNew%D(1:3,NatmsNew)=GMLoc%Displ%D(1:3,I)
         GradNew%D(1:3,NatmsNew)=GMLoc%Vects%D(1:3,I)
         AtNumNew%D(NatmsNew)=GMLoc%AtNum%D(NatmsNew)
       ENDIF
     ENDDO
     !
     !--------------------------------------------
     CALL OpenASCII(OutFile,Out)
       CALL ModifyGeom(GOpt,XYZNew%D,AtNumNew%D,GradNew%D, &
                       GMLoc%LagrMult%D,GMLoc%GradMult%D, &
                       GMLoc%LagrDispl%D,Convgd,GMLoc%Etotal, &
                       GMLoc%ELagr,IGeo,iCLONE,SCRPath,Opts%PFlags%GeOp)
     CLOSE(Out,STATUS='KEEP')
     !--------------------------------------------
     !
     ! Put back modified geometry into GMLoc array
     ! Also, put back gradients, which may be Lagrangian ones now.
     !
     NatmsNew=0
     DO I=1,GMLoc%Natms
       IF(GMLoc%CConstrain%I(I)/=2) THEN
         NatmsNew=NatmsNew+1
         GMLoc%Displ%D(1:3,I)=XYZNew%D(1:3,NatmsNew)
         GMLoc%Vects%D(1:3,I)=GradNew%D(1:3,NatmsNew)
       ENDIF
     ENDDO
     !
     ! Put modified (Energy->Lagrangian) gradients into HDF 
     !
     HDFFileID=OpenHDF(Nams%HFile)
     HDF_CurrentID=OpenHDFGroup(HDFFileID, &
                   "Clone #"//TRIM(IntToChar(iCLONE)))
       CALL New(CartGrad,3*GMLoc%Natms)
       CALL CartRNK2ToCartRNK1(CartGrad%D,GMLoc%Vects%D)
       CALL Put(CartGrad,'GradE',Tag_O=IntToChar(iGEO))
       CALL Delete(CartGrad)
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(HDFFileID)
     !
     CALL Delete(XYZNew)
     CALL Delete(AtNumNew)
     CALL Delete(GradNew)
   END SUBROUTINE OptSingleMol
!
!-------------------------------------------------------------------
!
   SUBROUTINE SetHessian(Hess)
     TYPE(Hessian) :: Hess
     Hess%Stre = 0.50D0   
     Hess%Bend = 0.20D0
     Hess%LinB = 0.20D0
     Hess%OutP = 0.10D0 
     Hess%Tors = 0.10D0 
     Hess%VDWStre  = 0.50D0 
     Hess%VDWBend  = 0.20D0 
     Hess%VDWLinB  = 0.20D0 
     Hess%VDWOutP  = 0.10D0 
     Hess%VDWTors  = 0.10D0 
   END SUBROUTINE SetHessian
!
!-------------------------------------------------------------------
!
   SUBROUTINE SetGConvCrit(GConv,Hess,AccL,NatmsLoc)
     TYPE(GConvCrit) :: GConv
     TYPE(Hessian)   :: Hess 
     INTEGER         :: AccL,NatmsLoc
     REAL(DOUBLE)    :: GCrit
     !
     GCrit=GTol(AccL)
     !
     GConv%MaxGeOpSteps=MAX(3*NatmsLoc,600)
     GConv%Grad= GCrit
     GConv%Stre=      GCrit/Hess%Stre
     GConv%Bend= 1.D0*GCrit/Hess%Bend
     GConv%LinB= 1.D0*GCrit/Hess%LinB
     GConv%OutP= 1.D0*GCrit/Hess%OutP
     GConv%Tors= 1.D0*GCrit/Hess%Tors
   END SUBROUTINE SetGConvCrit
!
!-------------------------------------------------------------------
!
   SUBROUTINE SetStepS(StepS)
     TYPE(StepSize) :: StepS
     !
     StepS%Stre = 0.1D0*AngstromsToAu
     StepS%Bend = 4.0D0*PI/180.D0
     StepS%LinB = 4.0D0*PI/180.D0
     StepS%OutP = 4.0D0*PI/180.D0
     StepS%Tors = 4.0D0*PI/180.D0
     StepS%Cart = 0.3D0*AngstromsToAu
   END SUBROUTINE SetStepS
!
!-------------------------------------------------------------------
!
   SUBROUTINE SetGDIIS(GD)
     TYPE(GDIIS)  :: GD
     !
     GD%Init    = 4
     GD%MaxMem  = 6
     GD%On=.TRUE.
   END SUBROUTINE SetGDIIS
!
!-------------------------------------------------------------------
!
   SUBROUTINE SetGrdTrf(GT,GConv)
     TYPE(GrdTrf)    :: GT
     TYPE(GConvCrit) :: GConv
     !
     GT%MaxIt_GrdTrf = 10 
    !GT%GrdTrfCrit   = 0.1D0*GConv%Grad
     GT%GrdTrfCrit   = 1.D-7
     GT%MaxGradDiff  = 5.D+2      
   END SUBROUTINE SetGrdTrf
!
!-------------------------------------------------------------------
!
   SUBROUTINE SetBackTrf(BackT,GConv)
     TYPE(BackTrf)   :: BackT
     TYPE(GConvCrit) :: GConv
     !
     BackT%MaxIt_CooTrf = 20
     BackT%CooTrfCrit   = MIN(GConv%Stre/10.D0,1.D-4)
     BackT%RMSCrit      = 0.75D0 
     BackT%MaxCartDiff  = 0.50D0  
     BackT%DistRefresh  = BackT%MaxCartDiff*BackT%RMSCrit
   END SUBROUTINE SetBackTrf
!
!-------------------------------------------------------------------
!
   SUBROUTINE SetConstr(Con,BackT)
     TYPE(Constr)    :: Con
     TYPE(BackTrf)   :: BackT
     !
     Con%ConstrMaxCrit = BackT%CooTrfCrit*1.D-2
     Con%ConstrMax     = Con%ConstrMaxCrit*10.D0
   END SUBROUTINE SetConstr
!
!-------------------------------------------------------------------
!
   SUBROUTINE SetTrfCtrl(TrfC,CoordC,GConstr)
     TYPE(TrfCtrl)   :: TrfC
     TYPE(CoordCtrl) :: CoordC
     TYPE(Constr)    :: GConstr
     !
     IF(.NOT.TrfC%DoClssTrf) THEN
       TrfC%DoTranslOff=.FALSE.
       TrfC%DoRotOff=.FALSE.
     ENDIF
     IF(GConstr%NCartConstr>0) TrfC%DoTranslOff=.FALSE.
     IF(GConstr%NCartConstr>12) TrfC%DoRotOff=.FALSE.
     IF(CoordC%CoordType/=CoordType_Cartesian) THEN
       TrfC%DoInternals=.TRUE.
     ELSE
       TrfC%DoInternals=.FALSE.
     ENDIF
   END SUBROUTINE SetTrfCtrl
!
!-------------------------------------------------------------------
!
   SUBROUTINE SetCoordCtrl(CoordC)
     TYPE(CoordCtrl) :: CoordC
     !
     CoordC%LinCrit =20.D0
     CoordC%OutPCrit=20.D0
   END SUBROUTINE SetCoordCtrl
!
!-------------------------------------------------------------------
!
   SUBROUTINE NewGeomFill(GMLoc)
     TYPE(CRDS) :: GMLoc
     !
     GMLoc%AbCarts%D=GMLoc%Displ%D
     GMLoc%LagrMult%D=GMLoc%LagrDispl%D
   END SUBROUTINE NewGeomFill
!
!-------------------------------------------------------------------
!
   SUBROUTINE MixGeoms(Opts,Nams,GOpt,Convgd,GMLoc,iCLONE,iGEO)
     TYPE(Options)        :: Opts
     TYPE(FileNames)      :: Nams
     TYPE(GeomOpt)        :: GOpt
     TYPE(CRDS)           :: GMLoc
     INTEGER              :: iGEO,iCLONE
     CHARACTER(LEN=DCL)   :: SCRPath
     INTEGER,DIMENSION(:) :: Convgd
     !
     SCRPath  =TRIM(Nams%M_SCRATCH)//TRIM(Nams%SCF_NAME)// &
             '.'//TRIM(IntToChar(iCLONE))
     !
     IF(Convgd(iCLONE)/=1) THEN
       CALL OPENAscii(OutFile,Out)
       IF((.NOT.GOpt%GDIIS%NoGDIIS).AND.GOpt%GDIIS%On) THEN
         CALL GeoDIIS(GMLoc%AbCarts%D,GMLoc%CConstrain%I, &
           GMLoc%LagrMult%D,GOpt%Constr,GOpt%BackTrf,GOpt%TrfCtrl, &
           GOpt%CoordCtrl,GOpt%GDIIS,Nams%HFile, &
           iCLONE,iGEO,Opts%PFlags%GeOp,SCRPath)
       ELSE
         IF(Opts%PFlags%GeOp>=DEBUG_GEOP_MIN) THEN
           WRITE(*,200)
           WRITE(Out,200)
         ENDIF
         200 FORMAT('No Geometric DIIS is being done in this step.')
       ENDIF
       CLOSE(Out,STATUS='KEEP')
     ENDIF
   END SUBROUTINE MixGeoms
!
!-------------------------------------------------------------------
!
   SUBROUTINE GetCGradMax(GConstr,CartGrad,NCart,IntCs, &
                          IMaxCGrad,MaxCGrad)
     TYPE(INTC)                :: IntCs
     TYPE(Constr)              :: GConstr
     REAL(DOUBLE),DIMENSION(:) :: CartGrad
     REAL(DOUBLE)              :: MaxCGrad,Sum
     INTEGER                   :: NCart,NIntC,I,JJ,J,IMaxCGrad
     LOGICAL,DIMENSION(NCart)  :: ConstrGrad
     !
     ConstrGrad=.TRUE.
     IF(GConstr%NCartConstr/=0) THEN
       NIntC=SIZE(IntCs%Def)
       DO I=1,NIntC
         IF(IntCs%Constraint(I)) THEN 
             JJ=IntCs%Atoms(I,1)
           IF(IntCs%Def(I)(1:5)=='CARTX') THEN
             J=3*(JJ-1)+1
             ConstrGrad(J)=.FALSE.
           ELSE IF(IntCs%Def(I)(1:5)=='CARTY') THEN
             J=3*(JJ-1)+2
             ConstrGrad(J)=.FALSE.
           ELSE IF(IntCs%Def(I)(1:5)=='CARTZ') THEN
             J=3*(JJ-1)+3
             ConstrGrad(J)=.FALSE.
           ENDIF
         ENDIF
       ENDDO
     ENDIF
     !
     IMaxCGrad=1   
     MaxCGrad=Zero
     DO I=1,NCart
       IF(ConstrGrad(I)) THEN
         Sum=ABS(CartGrad(I))
         IF(Sum>MaxCGrad) THEN
           IMaxCGrad=I
            MaxCGrad=Sum
         ENDIF
       ENDIF
     ENDDO
   END SUBROUTINE GetCGradMax
!
!-------------------------------------------------------------------
!
   SUBROUTINE BackTrack(iBAS,iGEO,C)
     ! Go over clones and do backtracking whenever necessary
     TYPE(Controls)   :: C
     INTEGER          :: iBAS,iGEO,iCLONE
     INTEGER          :: NatmsLoc,NCart
     INTEGER          :: MaxBStep,IBStep
     CHARACTER(LEN=DCL):: chGEO
     TYPE(CRDS)       :: GMOld
     LOGICAL          :: DoBackTrack
     REAL(DOUBLE)     :: EOld,ENew,MeanDist
     TYPE(DBL_VECT)   :: DistVect1,DistVect2
     !
     IF(.NOT.C%GOpt%GConvCrit%DoBackTr) THEN
       RETURN
     ENDIF
     MaxBStep=10
     DO iBStep=1,MaxBStep+1
       DoBackTrack=.FALSE.
       HDFFileID=OpenHDF(C%Nams%HFile)
       DO iCLONE=1,C%Geos%Clones
         HDF_CurrentID=OpenHDFGroup(HDFFileID, &
                     "Clone #"//TRIM(IntToChar(iCLONE)))
         chGEO=IntToChar(iGEO-1)
         CALL Get(GMOld,chGEO)
         CALL CloseHDFGroup(HDF_CurrentID)
         EOld=GMOld%ETotal
         ENew=C%Geos%Clone(iCLONE)%ETotal
         !
         NatmsLoc=GMOld%Natms
         NCart=3*NatmsLoc
         CALL New(DistVect1,NCart)
         CALL CartRNK2ToCartRNK1(DistVect1%D,GMOld%AbCarts%D)
         CALL New(DistVect2,NCart)
         CALL CartRNK2ToCartRNK1(DistVect2%D,C%Geos%Clone(iCLONE)%AbCarts%D)
         DistVect2%D=DistVect2%D-DistVect1%D
         MeanDist=SQRT(DOT_PRODUCT(DistVect2%D,DistVect2%D))/NatmsLoc
         CALL Delete(DistVect1)
         CALL Delete(DistVect2)
         !
         IF(EOld<ENew) DoBackTrack=.TRUE.
         !
         IF(iBStep>1.OR.DoBackTrack) THEN  
           CALL OPENAscii(OutFile,Out)
             WRITE(*,200) iBStep-1,EOld,ENew,MeanDist
             WRITE(Out,200) iBStep-1,EOld,ENew,MeanDist
           CLOSE(Out,STATUS='KEEP')
           200 FORMAT('Backtr. step= ',I3,' Old Energy= ', &
                       F14.8,' New Energy= ',F14.8,' Dist= ',F14.8)
         ENDIF
         !
         IF(DoBackTrack) THEN  
           ! do bisection
           C%Geos%Clone(iCLONE)%Carts%D= &
             (C%Geos%Clone(iCLONE)%Carts%D+GMOld%Carts%D)*Half
           C%Geos%Clone(iCLONE)%AbCarts%D= &
             (C%Geos%Clone(iCLONE)%AbCarts%D+GMOld%AbCarts%D)*Half
         ENDIF
         CALL CloseHDF(HDFFileID)
       ENDDO
       !
       IF(DoBackTrack) THEN
         IF(iBStep>MaxBStep) THEN
           CALL OPENAscii(OutFile,Out)
             WRITE(*,100) MaxBStep-1
             WRITE(Out,100) MaxBStep-1
           CLOSE(Out,STATUS='KEEP')
           100 FORMAT('Backtracking has not converged in ', &
                       I3,' Steps. Continue with present geometry.')
           EXIT
         ENDIF
         !
         CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)    
         CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
         CALL SCF(iBAS,iGEO,C)
       ELSE
         EXIT
       ENDIF
       !
     ENDDO
     ! In the present version do not do any overwriting of 
     ! the GMOld%Displ, since it refers 
     ! strictly to the simple relaxation step, for GDIIS.
     !
     CALL Delete(GMOld)
   END SUBROUTINE BackTrack
!
!-------------------------------------------------------------------
!
   SUBROUTINE LagrGradCart(GradIn,IntCs,LagrMult,SCRPath,iGEO)
     REAL(DOUBLE),DIMENSION(:,:)    :: GradIn
     TYPE(INTC)                     :: IntCs
     TYPE(BMATR)                    :: B
     REAL(DOUBLE),DIMENSION(:)      :: LagrMult
     INTEGER                        :: I,J,NIntC,NatmsLoc,NCart,iGEO
     CHARACTER(LEN=*)               :: SCRPath
     !
     NIntC=SIZE(IntCs%Def) 
     NatmsLoc=SIZE(GradIn,2)
     NCart=3*NatmsLoc
     CALL ReadBMATR(B,TRIM(SCRPath)//'B')
     !
     CALL GradAddCarts(GradIn,IntCs,LagrMult)
     CALL GradAddInt(GradIn,IntCs,B,LagrMult)
     !
     CALL Delete(B) 
   END SUBROUTINE LagrGradCart
!
!-------------------------------------------------------------------
!
   SUBROUTINE GradAddCarts(GradIn,IntCs,LagrMult)
     TYPE(INTC)                  :: IntCs
     REAL(DOUBLE),DIMENSION(:,:) :: GradIn
     REAL(DOUBLE),DIMENSION(:)   :: LagrMult
     REAL(DOUBLE)                :: Sum
     INTEGER                     :: I,J,NIntC,NConstr,NDim
     !
     ! Warning! At this point IntCs%Value should contain the actual
     ! values of the coordinates
     !
     NIntC=SIZE(IntCs%Def)
     NDim=SIZE(LagrMult)
     !
     NConstr=0
     DO I=1,NIntC
       IF(IntCs%Constraint(I)) THEN
         NConstr=NConstr+1
         Sum=-LagrMult(NConstr)
         IF(NConstr>NDim) CALL Halt('Dimension error #2 in GradAddCarts')
         J=IntCs%Atoms(I,1)
         IF(IntCs%Def(I)(1:5)=='CARTX') THEN
           GradIn(1,J)=GradIn(1,J)+Sum
         ELSE IF(IntCs%Def(I)(1:5)=='CARTY') THEN
           GradIn(2,J)=GradIn(2,J)+Sum
         ELSE IF(IntCs%Def(I)(1:5)=='CARTZ') THEN
           GradIn(3,J)=GradIn(3,J)+Sum
         ENDIF
       ENDIF
     ENDDO
   END SUBROUTINE GradAddCarts
!
!-------------------------------------------------------------------
!
   SUBROUTINE GradAddInt(GradIn,IntCs,B,LagrMult)
     REAL(DOUBLE),DIMENSION(:,:):: GradIn
     REAL(DOUBLE),DIMENSION(:)  :: LagrMult
     REAL(DOUBLE)               :: Sum
     TYPE(INTC)                 :: IntCs
     TYPE(BMATR)                :: B 
     INTEGER                    :: I,J,K,L,JJ,NintC,NConstr,NDim
     !
     NintC=SIZE(IntCs%Def)
     NDim=SIZE(LagrMult)
     NConstr=0
     IF(NintC/=SIZE(B%IB,1)) CALL Halt('Dimension error in GradAddInt')
     DO I=1,NintC
       IF(IntCs%Constraint(I)) THEN
         NConstr=NConstr+1
         Sum=-LagrMult(NConstr)
         IF(NConstr>NDim) &
           CALL Halt('Dimension error #2 in GradAddInt')
         IF(IntCs%Def(I)(1:4)/='CART') THEN
           DO J=1,4
             L=IntCs%Atoms(I,J)
             IF(L==0) EXIT
             JJ=3*(J-1)
             DO K=1,3 
               GradIn(K,L)=GradIn(K,L)+Sum*B%B(I,JJ+K) 
             ENDDO
           ENDDO
         ENDIF
       ENDIF
     ENDDO
   END SUBROUTINE GradAddInt
!
!-------------------------------------------------------------------
!
   SUBROUTINE LagrGradMult(GradMult,XYZ,IntCs,LinCrit,NConstrIn)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(:) :: GradMult
     REAL(DOUBLE)              :: LinCrit,Sum
     TYPE(INTC)                :: IntCs
     INTEGER                   :: I,J,NIntC,NConstr,NConstrIn
     !
     NIntC=SIZE(IntCs%Def)
     NConstr=0
     CALL INTCValue(IntCs,XYZ,LinCrit)
     IF(SIZE(GradMult)/=NConstrIn) &
       CALL Halt('Dimension error #1 in LagrGradMult')
     DO I=1,NIntC
       IF(IntCs%Constraint(I)) THEN
         NConstr=NConstr+1
         IF(NConstr>NConstrIn) &
           CALL Halt('Dimension error #2 in LagrGradMult')
         Sum=IntCs%Value(I)-IntCs%ConstrValue(I)
         GradMult(NConstr)=-Sum
       ENDIF
     ENDDO
   END SUBROUTINE LagrGradMult
!
!----------------------------------------------------------------
!
   SUBROUTINE DiagHessLagr(GCoordCtrl,GHessian, &
        Grad,Displ,IntCs,XYZ,SCRPath,LagrMult,GradMult,LagrDispl)
     REAL(DOUBLE),DIMENSION(:)  :: Grad,Displ
     REAL(DOUBLE),DIMENSION(:)  :: LagrMult,GradMult,LagrDispl
     REAL(DOUBLE),DIMENSION(:,:):: XYZ
     TYPE(CoordCtrl)            :: GCoordCtrl
     TYPE(Hessian)              :: GHessian
     TYPE(INTC)                 :: IntCs
     CHARACTER(LEN=*)           :: SCRPath
     INTEGER                    :: I,J,K,L,NIntC,NCart
     INTEGER                    :: NatmsLoc,NConstr
     TYPE(INT_VECT)             :: IHessL,JHessL
     TYPE(DBL_VECT)             :: AHessL
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     NConstr=SIZE(LagrMult)
     CALL LagrInvHess(IntCs,SCRPath,GHessian, &
                      LagrMult,NCart,NConstr, &
                      IHessL,JHessL,AHessL)
     CALL DiagDispl(IHessL,JHessL,AHessL,SCRPath, &
                    Grad,GradMult,Displ,LagrDispl)
     ! 
     CALL Delete(IHessL)
     CALL Delete(JHessL)
     CALL Delete(AHessL)
   END SUBROUTINE DiagHessLagr
!
!-------------------------------------------------------------------
!
   SUBROUTINE GrdConvrgd(GStat,IntCs,Grad)
     TYPE(GOptStat)            :: GStat
     TYPE(INTC)                :: IntCs
     REAL(DOUBLE),DIMENSION(:) :: Grad
     REAL(DOUBLE)              :: Sum
     INTEGER                   :: I,J,NDim,NIntC,NCart
     !
     NDim=SIZE(Grad)
     NIntC=SIZE(IntCs%Def)
     !
     GStat%IMaxGrad=1
     GStat%MaxGrad=ABS(Grad(1))
     DO I=2,NDim
       Sum=ABS(Grad(I))
       IF(Sum>GStat%MaxGrad) THEN
         GStat%IMaxGrad=I
          GStat%MaxGrad=Sum
       ENDIF
     ENDDO
     GStat%RMSGrad=SQRT(DOT_PRODUCT(Grad,Grad)/DBLE(NDim))
     !
     ! Check for gradient-convergence in the presence of constraints
     !
     GStat%MaxGradNoConstr=Zero
     GStat%RMSGradNoConstr=Zero
     GStat%IMaxGradNoConstr=0
     J=0
     DO I=1,NIntC
       IF(.NOT.IntCs%Constraint(I)) THEN
         J=J+1
         Sum=ABS(Grad(I))
         IF(GStat%MaxGradNoConstr<Sum) THEN
           GStat%IMaxGradNoConstr=I
           GStat%MaxGradNoConstr=Sum
         ENDIF
         GStat%RMSGradNoConstr=GStat%RMSGradNoConstr+Sum*Sum
       ENDIF
     ENDDO
     IF(J/=0) GStat%RMSGradNoConstr=SQRT(GStat%RMSGradNoConstr)/DBLE(J)
     !
     ! Check gradients of translation and rotation
     !
   END SUBROUTINE GrdConvrgd
!
!-------------------------------------------------------------------
!
   SUBROUTINE LagrConv(GStat,LagrDispl,LagrMult,GradMult)
     TYPE(GOptStat)            :: GStat 
     REAL(DOUBLE),DIMENSION(:) :: LagrDispl,LagrMult,GradMult
     REAL(DOUBLE)              :: Sum
     INTEGER                   :: I,J,NLagr
     !
     NLagr=SIZE(LagrMult)
     IF(NLagr==0) THEN
       GStat%IMaxLGrad=0
       GStat%MaxLGrad=Zero
       GStat%MaxDMult=Zero
       RETURN
     ENDIF
     GStat%IMaxLGrad=1
     GStat%MaxLGrad=ABS(GradMult(1))
     GStat%MaxDMult=ABS(LagrDispl(1)-LagrMult(1))
     DO I=2,NLagr
       Sum=ABS(GradMult(I)) 
       IF(GStat%MaxLGrad<Sum) THEN
         GStat%IMaxLGrad=I
         GStat%MaxLGrad=Sum
       ENDIF
       Sum=ABS(LagrDispl(I)-LagrMult(I)) 
       IF(GStat%MaxDMult<Sum) THEN
         GStat%MaxDMult=Sum
       ENDIF
     ENDDO
   END SUBROUTINE LagrConv
!
!-------------------------------------------------------------------
!
   SUBROUTINE CalcLagrMult(GradMult,LagrMult,LagrDispl,XYZ,IntCs, &
                           GradIn,GGrdTrf,GCoordCtrl,GTrfCtrl,GConstr,&
                           GHess,Print,SCRPath)
     !
     REAL(DOUBLE),DIMENSION(:)   :: GradMult,LagrMult,LagrDispl
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ,GradIn
     REAL(DOUBLE)                :: Sum,Hess
     TYPE(INTC)                  :: IntCs
     INTEGER                     :: I,J,NLagr,NConstr,NIntC,NCart
     INTEGER                     :: NatmsLoc 
     INTEGER                     :: Print
     CHARACTER(LEN=*)            :: SCRPath
     TYPE(DBL_VECT)              :: CartGrad,IntGrad
     TYPE(Constr)                :: GConstr
     TYPE(TrfCtrl)               :: GTrfCtrl
     TYPE(CoordCtrl)             :: GCoordCtrl
     TYPE(GrdTrf)                :: GGrdTrf
     TYPE(Hessian)               :: GHess   
     TYPE(INT_VECT)              :: IBc,JBc
     TYPE(DBL_VECT)              :: ABc
     !
     NLagr=SIZE(GradMult)
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     NIntC=SIZE(IntCs%Def)
     !
     ! Calculate constraint B-matrix
     !
   ! CALL BMatrConstr(IBc,JBc,ABc,IntCs,SCRPath,NLagr,NCart)
     !
     ! Calculate new values of Lagrange multipliers 
     ! on internal coord constraints
     !
     CALL New(CartGrad,NCart)
     CALL CartRNK2ToCartRNK1(CartGrad%D,GradIn)
     CALL New(IntGrad,NIntC)
     !
     IF(GConstr%NConstr/=GConstr%NCartConstr) THEN
       CALL CartToInternal(XYZ,IntCs,CartGrad%D,IntGrad%D, &
         GGrdTrf,GCoordCtrl,GTrfCtrl,Print,SCRPath)
     ELSE
       IntGrad%D=Zero
     ENDIF
     !
     NConstr=0
     DO I=1,NIntC
       IF(IntCs%Constraint(I)) THEN
         NConstr=NConstr+1
         IF(NConstr>NLagr) CALL Halt('Dimension error in CalcLagrMult') 
         Sum=Zero
         !Hess=ABS(LagrMult(NConstr))
         Hess=GHess%Stre
         J=IntCs%Atoms(I,1)
         IF(IntCs%Def(I)(1:4)/='CART') THEN
           Sum=IntGrad%D(I)
           IF(IntCs%Def(I)(1:4)=='STRE') Hess=GHess%Stre
           IF(IntCs%Def(I)(1:4)=='BEND') Hess=GHess%Bend
           IF(IntCs%Def(I)(1:4)=='LINB') Hess=GHess%LinB
           IF(IntCs%Def(I)(1:4)=='TORS') Hess=GHess%Tors
           IF(IntCs%Def(I)(1:4)=='OUTP') Hess=GHess%OutP
         ELSE IF(IntCs%Def(I)(1:5)=='CARTX') THEN
           Sum=GradIn(1,J)
         ELSE IF(IntCs%Def(I)(1:5)=='CARTY') THEN
           Sum=GradIn(2,J)
         ELSE IF(IntCs%Def(I)(1:5)=='CARTZ') THEN
           Sum=GradIn(3,J)
         ENDIF
         LagrDispl(NConstr)=Hess*GradMult(NConstr)+Sum
       ENDIF 
     ENDDO
     !
     CALL Delete(CartGrad)
     CALL Delete(IntGrad)
!    CALL Delete(IBc)
!    CALL Delete(JBc)
!    CALL Delete(ABc)
   END SUBROUTINE CalcLagrMult
!
!-------------------------------------------------------------------
!
   SUBROUTINE LagrEnergy(LagrMult,IntCs,ELagr)
     REAL(DOUBLE),DIMENSION(:)       :: LagrMult
     REAL(DOUBLE)                    :: ELagr    
     TYPE(INTC)                      :: IntCs
     INTEGER                         :: I,J,NIntC,NLagr,NConstr
     !
     NIntC=SIZE(IntCs%Def) 
     NLagr=SIZE(LagrMult)
     ELagr=Zero
     IF(NLagr==0) RETURN     
     NConstr=0
     DO I=1,NIntC
       IF(IntCs%Constraint(I)) THEN
         NConstr=NConstr+1
         IF(NConstr>NLagr) CALL Halt('NConstr>NLagr in LagrEnergy')
         ELagr=ELagr-LagrMult(NConstr)*(IntCs%Value(I)-IntCs%ConstrValue(I)) 
       ENDIF
     ENDDO
   END SUBROUTINE LagrEnergy
!
!-------------------------------------------------------------------
!
   SUBROUTINE DiagHess_Vals(DHess,AtNum,XYZ,IntCs,CoordC,Hess,NStre,Char)
     TYPE(Hessian)                 :: Hess
     TYPE(CoordCtrl)               :: CoordC
     REAL(DOUBLE),DIMENSION(:)     :: DHess,AtNum
     REAL(DOUBLE),DIMENSION(:,:)   :: XYZ  
     TYPE(INTC)                    :: IntCs 
     INTEGER                       :: I,J,NIntC,NatmsLoc,NStre
     CHARACTER(LEN=*)              :: Char
     TYPE(INT_VECT)                :: ITop,JTop
     TYPE(DBL_VECT)                :: ATop
     LOGICAL                       :: DoVDW
     !
     NIntC=SIZE(IntCs%Def)
     NatmsLoc=SIZE(XYZ,2)
     IF(NIntC/=SIZE(DHess)) &
       CALL Halt('Dimesion error in DiagHess_Vals.')
     !
     IF(Char=='Lindh') CALL DistMatr(ITop,JTop,ATop,IntCs,NatmsLoc,NStre)
     !
     DO I=1,NIntC
       IF(.NOT.IntCs%Active(I)) THEN
         DHess(I)=Zero
       ELSE 
         DoVDW=(I<CoordC%NCov)
         CALL CalcHess(DHess(I),Char,IntCs%Def(I),Hess,AtNum, &
                       XYZ,IntCs%Atoms(I,1:4),ITop,JTop,ATop,DoVDW)
       ENDIF
     ENDDO
     !
     IF(Char=='Lindh') THEN 
       CALL Delete(ITop)
       CALL Delete(JTop)
       CALL Delete(ATop)
     ENDIF
   END SUBROUTINE DiagHess_Vals
!
!-------------------------------------------------------------------
!
   FUNCTION PeriodicRow(I)
     INTEGER   ::  PeriodicRow,I
     !
     IF(I<=0) THEN
       PeriodicRow=0
     ELSE IF(0<I.AND.I<=2) THEN
       PeriodicRow=1
     ELSE IF(2<I.AND.I<=10) THEN
       PeriodicRow=2
     ELSE IF(10<I.AND.I<=18) THEN
       PeriodicRow=3
     ELSE IF(18<I.AND.I<=36) THEN
       PeriodicRow=4
     ELSE IF(36<I.AND.I<=54) THEN
       PeriodicRow=5
     ELSE IF(54<I.AND.I<=86) THEN
       PeriodicRow=5
     ELSE 
       PeriodicRow=7
     ENDIF
   END FUNCTION PeriodicRow
!
!-------------------------------------------------------------------
!
   SUBROUTINE CalcHess(DHess,Char,Type,Hess,AtNum,XYZ, &
                       Atoms,ITop,JTop,ATop,DoVDW)
     REAL(DOUBLE)                :: DHess
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(:)   :: AtNum
     TYPE(Hessian)               :: Hess
     CHARACTER(LEN=*)            :: Char,Type
     INTEGER,DIMENSION(1:4)      :: Atoms 
     TYPE(INT_VECT)              :: ITop,JTop
     TYPE(DBL_VECT)              :: ATop
     INTEGER                     :: I1Row,I2Row,I3Row,I4Row  
     REAL(DOUBLE)                :: R12,R23,R34
     REAL(DOUBLE)                :: Rho12,Rho23,Rho34
     LOGICAL                     :: DoVDW
     !
     IF(Type(1:4)=='CART') THEN
       DHess=Hess%VDWStre
     ELSE IF(DoVDW) THEN
       IF(Type(1:4)=='STRE') DHess=Hess%VDWStre
       IF(Type(1:4)=='BEND') DHess=Hess%VDWBend
       IF(Type(1:4)=='LINB') DHess=Hess%VDWLinB
       IF(Type(1:4)=='OUTP') DHess=Hess%VDWOutP
       IF(Type(1:4)=='TORS') DHess=Hess%VDWTors
     ELSE
       IF(Char=='ThreeVals') THEN
         IF(Type(1:4)=='STRE') DHess=Hess%Stre
         IF(Type(1:4)=='BEND') DHess=Hess%Bend
         IF(Type(1:4)=='LINB') DHess=Hess%LinB
         IF(Type(1:4)=='OUTP') DHess=Hess%OutP
         IF(Type(1:4)=='TORS') DHess=Hess%Tors
       ELSE IF(Char=='Lindh') THEN
         I1Row=PeriodicRow(INT(AtNum(Atoms(1))))
         I2Row=PeriodicRow(INT(AtNum(Atoms(2))))
         I3Row=PeriodicRow(INT(AtNum(Atoms(3))))
         I4Row=PeriodicRow(INT(AtNum(Atoms(4))))
         IF(Atoms(2)/=0) THEN
           R12=GetR(XYZ,Atoms(1),Atoms(2),ITop,JTop,ATop) 
           Rho12=EXP(Lindh_Alpha(I1Row,I2Row)*(Lindh_R(I1Row,I2Row)**2-R12**2))
         ENDIF
         IF(Atoms(3)/=0) THEN
           R23=GetR(XYZ,Atoms(2),Atoms(3),ITop,JTop,ATop) 
           Rho23=EXP(Lindh_Alpha(I2Row,I3Row)*(Lindh_R(I2Row,I3Row)**2-R23**2))
         ENDIF
         IF(Atoms(4)/=0) THEN
           IF(Type=='OUTP') THEN
             R34=GetR(XYZ,Atoms(2),Atoms(4),ITop,JTop,ATop) 
             Rho34=EXP(Lindh_Alpha(I2Row,I4Row)*(Lindh_R(I2Row,I4Row)**2-R34**2))
           ELSE
             R34=GetR(XYZ,Atoms(3),Atoms(4),ITop,JTop,ATop) 
             Rho34=EXP(Lindh_Alpha(I3Row,I4Row)*(Lindh_R(I3Row,I4Row)**2-R34**2))
           ENDIF
         ENDIF
         DHess=One
         IF(Type(1:4)=='STRE') THEN
           DHess=Lindh_K(1)*Rho12
         ELSE IF(Type(1:4)=='BEND') THEN
           DHess=Lindh_K(2)*Rho12*Rho23
         ELSE IF(Type(1:4)=='LINB') THEN
           DHess=Lindh_K(2)*Rho12*Rho23
         ELSE IF(Type(1:4)=='OUTP') THEN
           DHess=Lindh_K(3)*Rho12*Rho23*Rho34
         ELSE IF(Type(1:4)=='TORS') THEN
           DHess=Lindh_K(3)*Rho12*Rho23*Rho34
         ELSE IF(Type(1:4)=='CART') THEN
           DHess=Lindh_K(1)*Rho12
         ENDIF
       ENDIF
     ENDIF
     DHess=One/DHess
   END SUBROUTINE CalcHess
!
!-------------------------------------------------------------------
!
   FUNCTION GetR(XYZ,I1Row,I2Row,ITop,JTop,ATop)
     REAL(DOUBLE)                  :: GetR
     REAL(DOUBLE),DIMENSION(:,:)   :: XYZ 
     INTEGER                       :: I1Row,I2Row,I,J,III,NDim
     TYPE(INT_VECT)                :: ITop,JTop
     TYPE(DBL_VECT)                :: ATop
     !
     NDim=SIZE(ITop%I)-1
     III=0
     DO J=ITop%I(I1Row),ITop%I(I1Row+1)-1
       IF(JTop%I(J)==I2Row) THEN
         GetR=ATop%D(J)
         III=1
       ENDIF
     ENDDO
     IF(III==0) THEN
       GetR=SQRT((XYZ(1,I1Row)-XYZ(1,I2Row))**2+&
                 (XYZ(2,I1Row)-XYZ(2,I2Row))**2+&
                 (XYZ(3,I1Row)-XYZ(3,I2Row))**2)
     ENDIF 
   END FUNCTION GetR
!
!-------------------------------------------------------------------
!
END MODULE Optimizer
