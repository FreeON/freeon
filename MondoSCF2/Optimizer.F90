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
!---------------------------------------------------------------------
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
       CALL Force(iBAS,iGEO,C%Nams,C%Opts,C%Stat,C%Geos,C%Sets,C%MPIs)
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
       ! Archive displaced geometries and fill in new geoms.
       !
       CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)    
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
   SUBROUTINE ModifyGeom(GOpt,XYZ,AtNum,GradIn,Convgd,ETot, &
                         iGEO,iCLONE,SCRPath,Print)
     TYPE(GeomOpt)               :: GOpt
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ,GradIn
     REAL(DOUBLE),DIMENSION(:)   :: AtNum
     REAL(DOUBLE)                :: ETot 
     INTEGER,DIMENSION(:)        :: Convgd
     INTEGER                     :: NIntC,NCart
     INTEGER                     :: NatmsLoc,iGEO,iCLONE
     TYPE(INTC)                  :: IntCs
     TYPE(DBL_VECT)              :: IntOld,Displ
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
     ENDIF
     IF(NIntC/=0) THEN
       CALL INTCValue(IntCs,XYZ,GOpt%CoordCtrl%LinCrit)
       CALL New(IntOld,NIntC)
       IntOld%D=IntCs%Value
     ENDIF
     !
     ! Print current geometry for debugging
     !
     IF(Print==DEBUG_GEOP_MAX) CALL PrtIntCoords(IntCs, &
       IntCs%Value,'Internals at step #'//TRIM(IntToChar(iGEO)))
     !
     ! Calculate simple relaxation (SR) step from an inverse Hessian
     !
     CALL NewDispl(GOpt,Displ,NCart,NIntC)
     CALL SRStep(GOpt,XYZ,AtNum,GradIn,Displ,IntCs,SCRPath,Print) 
     !
     ! Calculate new geometry 
     !
     CALL NewStructure(Print,GOpt%BackTrf,GOpt%TrfCtrl,GOpt%CoordCtrl, &
       GOpt%Constr,SCRPath,XYZ,Displ,IntCs)
     CALL Delete(Displ)
     !
     ! Check convergence
     !
     CALL GeOpConv(GOpt%Constr,GOpt%GOptStat,GOpt%CoordCtrl, &
                   GOpt%GConvCrit,XYZ,ETot,IntCs,IntOld,iCLONE,iGEO)
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
   SUBROUTINE SRStep(GOpt,XYZ,AtNum,GradIn,Displ,IntCs,SCRPath,Print)
     !
     ! Simple Relaxation step
     !
     TYPE(GeomOpt)                  :: GOpt 
     REAL(DOUBLE),DIMENSION(:,:)    :: XYZ,GradIn
     REAL(DOUBLE),DIMENSION(:)      :: AtNum
     TYPE(DBL_VECT)                 :: Displ
     REAL(DOUBLE)                   :: MaxGrad,RMSGrad
     INTEGER                        :: IMaxGradNoConstr
     REAL(DOUBLE)                   :: MaxGradNoConstr
     REAL(DOUBLE)                   :: RMSGradNoConstr,Sum
     TYPE(DBL_VECT)                 :: IntGrad,Grad,CartGrad
     TYPE(INTC)                     :: IntCs
     INTEGER                        :: I,J,NDim,IMaxGrad
     INTEGER                        :: NatmsLoc,NCart,NIntC,Print
     LOGICAL                        :: DoInternals
     CHARACTER(LEN=*)               :: SCRPath 
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     IF(AllocQ(IntCs%Alloc)) THEN
       NIntC=SIZE(IntCs%Def)
     ELSE
       NIntC=0
     ENDIF
       NDim =SIZE(Displ%D)
     DoInternals=GOpt%TrfCtrl%DoInternals
     IF(NIntC/=NDim.AND.DoInternals) &
         CALL Halt('Dimensionality error in SRStep')
     !
     CALL New(Grad,NDim)
     CALL New(CartGrad,NCart)
     CALL CartRNK2ToCartRNK1(CartGrad%D,GradIn)
     !
     ! If requested, compute internal coord. gradients
     !
     IF(DoInternals) THEN
       NIntC=SIZE(IntCs%Def)
       IF(NIntC/=NDim) CALL Halt('Dimension error in SRStep')
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
     MaxGrad=Zero
     DO I=1,NDim  
       Sum=ABS(Grad%D(I))
       IF(MaxGrad<Sum) THEN
         IMaxGrad=I
          MaxGrad=Sum
       ENDIF
     ENDDO
     RMSGrad=SQRT(DOT_PRODUCT(Grad%D,Grad%D)/DBLE(NDim))
     !
     ! Check for gradient-convergence in the presence of constraints
     !
     MaxGradNoConstr=Zero
     RMSGradNoConstr=Zero
     J=0
     DO I=1,NIntC
       IF(.NOT.IntCs%Constraint(I)) THEN
         J=J+1
         Sum=ABS(Grad%D(I))
         IF(MaxGradNoConstr<Sum) THEN
           IMaxGradNoConstr=I
           MaxGradNoConstr=Sum
         ENDIF
         RMSGradNoConstr=RMSGradNoConstr+Sum*Sum
       ENDIF
     ENDDO
     IF(J/=0) RMSGradNoConstr=SQRT(RMSGradNoConstr)/DBLE(J)
     !
     ! Use Hessian matrix to calculate step
     !
     SELECT CASE(GOpt%Optimizer)
     CASE(GRAD_STPDESC_OPT) 
       CALL SteepestDesc(GOpt%CoordCtrl,GOpt%Hessian, &
                         Grad,Displ,XYZ)
     CASE(GRAD_DIAGHESS_OPT) 
       CALL DiagonalHess(GOpt%CoordCtrl,GOpt%Hessian, &
                         Grad,Displ,IntCs,XYZ)
     ! CALL DiagHessRFO(GOpt,Grad,Displ,IntCs,XYZ)
     ! CALL RedundancyOff(GOpt,Displ%D)
     END SELECT
     !
     ! Set constraints on the displacements
     !
     CALL SetConstraint(IntCs,XYZ,Displ,GOpt%CoordCtrl%LinCrit, &
       GOpt%Constr%NConstr,GOpt%TrfCtrl%DoInternals)
     !
     ! Tidy up
     !
     CALL Delete(Grad)
     GOpt%GOptStat%MaxGrad=MaxGrad
     GOpt%GOptStat%MaxGradNoConstr=MaxGradNoConstr
     GOpt%GOptStat%IMaxGrad=IMaxGrad
     GOpt%GOptStat%RMSGrad=RMSGrad
     GOpt%GOptStat%RMSGradNoConstr=RMSGradNoConstr
     GOpt%GOptStat%IMaxGradNoConstr=IMaxGradNoConstr
   END SUBROUTINE SRStep
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
       Displ%D=-5.D0*Grad%D 
     ENDIF 
   END SUBROUTINE SteepestDesc
!
!-------------------------------------------------------
!
   SUBROUTINE DiagonalHess(CoordC,Hess,Grad,Displ,IntCs,XYZ)
     !
     TYPE(CoordCtrl)   :: CoordC
     TYPE(Hessian)     :: Hess
     TYPE(DBL_VECT)    :: Grad,Displ 
     TYPE(INTC)        :: IntCs       
     INTEGER           :: I,J,NIntC,NCart
     REAL(DOUBLE)      :: HStre,HBend,HLinB,HOutP,HTors
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     !
     NIntC=SIZE(IntCs%Def)
     NCart=3*SIZE(XYZ,2)
     !
     IF(CoordC%CoordType==CoordType_Cartesian) THEN
       Displ%D=-1.D0*Grad%D !!!! equivalent with stpdesc
     ELSE IF(CoordC%CoordType==CoordType_PrimInt) THEN
         HStre=One/Hess%Stre
         HBend=One/Hess%Bend
         HLinB=One/Hess%LinB
         HOutP=One/Hess%OutP
         HTors=One/Hess%Tors
       DO I=1,NIntC
         IF(.NOT.IntCs%Active(I)) THEN
           Displ%D(I)=Zero
         ELSE IF(IntCs%Def(I)(1:4)=='STRE') THEN
           Displ%D(I)=-HStre*Grad%D(I)
         ELSE IF(IntCs%Def(I)(1:4)=='BEND') THEN
           Displ%D(I)=-HBend*Grad%D(I)
         ELSE IF(IntCs%Def(I)(1:4)=='LINB') THEN
           Displ%D(I)=-HLinB*Grad%D(I)
         ELSE IF(IntCs%Def(I)(1:4)=='OUTP') THEN
           Displ%D(I)=-HOutP*Grad%D(I)
         ELSE IF(IntCs%Def(I)(1:4)=='TORS') THEN
           Displ%D(I)=-HTors*Grad%D(I)
         ENDIF
       ENDDO
       !
       !project out redundancy
       !         CALL RedundancyOff(Displ%D,XYZ,DoSet_O=.TRUE.)
     ELSE
       CALL Halt('Only Primitiv Internals are available yet.')
     ENDIF 
   END SUBROUTINE DiagonalHess
!
!-------------------------------------------------------
!
   SUBROUTINE DiagHessRFO(GOpt,Grad,Displ,IntCs,XYZ)
     TYPE(GeomOpt)     :: GOpt
     TYPE(DBL_VECT)    :: Grad,Displ
     TYPE(DBL_RNK2)    :: Hessian
     TYPE(INTC)        :: IntCs
     INTEGER           :: I,J,NIntC,NCart,Info
     REAL(DOUBLE)      :: HStre,HBend,HLinB,HOutP,HTors
     REAL(DOUBLE)      :: Sum
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
                       XYZ,Etot,IntCs,IntOld,iCLONE,iGEO)
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
     REAL(DOUBLE)               :: RMSGrad,MaxGrad
     REAL(DOUBLE)               :: RMSGradNoConstr,MaxGradNoConstr
     REAL(DOUBLE)               :: Etot,Sum
     !                      
     INTEGER                    :: NStreGeOp,NBendGeOp,NLinBGeOp
     INTEGER                    :: NOutPGeOp,NTorsGeOp
     INTEGER                    :: MaxStre,MaxBend,MaxLinB
     INTEGER                    :: MaxOutP,MaxTors
     !                      
     INTEGER                    :: IMaxGrad,IMaxGradNoConstr
     !
     IF(AllocQ(IntCs%Alloc)) THEN
       NIntC=SIZE(IntCs%Def)
     ELSE
       NIntC=0
     ENDIF
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     !
     RMSGrad=CtrlStat%RMSGrad
     RMSGradNoConstr=CtrlStat%RMSGradNoConstr
     MaxGrad=CtrlStat%MaxGrad
     MaxGradNoConstr=CtrlStat%MaxGradNoConstr
     IMaxGrad=CtrlStat%IMaxGrad
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
     IF(CtrlConstr%NConstr/=0) THEN
       CtrlStat%GeOpConvgd=&
                           MaxStreDispl<GConvCr%Stre.AND. &
                           MaxBendDispl<GConvCr%Bend.AND. &
                           MaxLinBDispl<GConvCr%LinB.AND. &
                           MaxOutPDispl<GConvCr%OutP.AND. &
                           MaxTorsDispl<GConvCr%Tors
     ELSE
     ! no constraints
       CtrlStat%GeOpConvgd=RMSGrad<GConvCr%Grad.AND. &
                           MaxGrad<GConvCr%Grad.AND. &
                           MaxStreDispl<GConvCr%Stre.AND. &
                           MaxBendDispl<GConvCr%Bend.AND. &
                           MaxLinBDispl<GConvCr%LinB.AND. &
                           MaxOutPDispl<GConvCr%OutP.AND. &
                           MaxTorsDispl<GConvCr%Tors
     ENDIF
     !
     ! Review iterations
     !
     WRITE(*,399) iCLONE,iGEO,ETot
     WRITE(Out,399) iCLONE,iGEO,ETot
     !
     MaxStreDispl=MaxStreDispl/AngstromsToAu
     MaxBendDispl=MaxBendDispl*180.D0/PI
     MaxLinBDispl=MaxLinBDispl*180.D0/PI
     MaxOutPDispl=MaxOutPDispl*180.D0/PI
     MaxTorsDispl=MaxTorsDispl*180.D0/PI
     !
     WRITE(*,410) MaxGrad,IntCs%Atoms(IMaxGrad,1:4)
     WRITE(*,420) RMSGrad
     WRITE(Out,410) MaxGrad,IntCs%Atoms(IMaxGrad,1:4)
     WRITE(Out,420) RMSGrad
     IF(CtrlConstr%NConstr/=0) THEN
       WRITE(*,510) MaxGradNoConstr,IntCs%Atoms(IMaxGradNoConstr,1:4)
       WRITE(*,520) RMSGradNoConstr
       WRITE(Out,510) MaxGradNoConstr,IntCs%Atoms(IMaxGradNoConstr,1:4)
       WRITE(Out,520) RMSGradNoConstr
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
     !
     WRITE(*,440) RMSIntDispl
     WRITE(Out,440) RMSIntDispl
     !
399 FORMAT('       Clone = ',I6,' GeOp step = ',I6,' Total Energy = ',F20.8)
400 FORMAT('Total Energy at Current Geometry = ',F20.8)
401 FORMAT('                    Total Energy = ',F20.8)
410 FORMAT('                        Max Grad = ',F12.6,' between atoms ',4I4)
420 FORMAT('                        RMS Grad = ',F12.6)
510 FORMAT('Max Grad on Unconstrained Coords = ',F12.6,' between atoms ',4I4)
520 FORMAT('RMS Grad on Unconstrained Coords = ',F12.6)
430 FORMAT('                  Max STRE Displ = ',F12.6,' between atoms ',4I4)
435 FORMAT('                  Max BEND Displ = ',F12.6,' between atoms ',4I4)
436 FORMAT('                  Max LINB Displ = ',F12.6,' between atoms ',4I4)
437 FORMAT('                  Max OUTP Displ = ',F12.6,' between atoms ',4I4)
438 FORMAT('                  Max TORS Displ = ',F12.6,' between atoms ',4I4)
440 FORMAT('                       RMS Displ = ',F12.6)
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
     CALL   SetTrfCtrl(GOpt%TrfCtrl,GOpt%CoordCtrl)
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
     INTEGER              :: iGEO,iCLONE
     INTEGER              :: InitGDIIS,NConstr,NCart,NatmsLoc
     INTEGER              :: Print
     LOGICAL              :: NoGDIIS,GDIISOn
     CHARACTER(LEN=DCL)   :: SCRPath
     TYPE(DBL_RNK2)       :: XYZNew
     !
     SCRPath  =TRIM(Nams%M_SCRATCH)//TRIM(Nams%SCF_NAME)// &
             '.'//TRIM(IntToChar(iCLONE))
     Print    =Opts%PFlags%GeOp
     GMLoc%Displ%D=GMLoc%AbCarts%D
     !
     !--------------------------------------------
     CALL OpenASCII(OutFile,Out)
       CALL ModifyGeom(GOpt,GMLoc%Displ%D,GMLoc%AtNum%D,GMLoc%Vects%D, &
                       Convgd,GMLoc%Etotal,IGeo,iCLONE,SCRPath,Print)
     CLOSE(Out,STATUS='KEEP')
     !--------------------------------------------
     !
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
     GT%GrdTrfCrit   = 0.1D0*GConv%Grad
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
   SUBROUTINE SetTrfCtrl(TrfC,CoordC)
     TYPE(TrfCtrl)   :: TrfC
     TYPE(CoordCtrl) :: CoordC
     !
     IF(.NOT.TrfC%DoClssTrf) THEN
       TrfC%DoTranslOff=.FALSE.
       TrfC%DoRotOff=.FALSE.
     ENDIF
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
   END SUBROUTINE NewGeomFill
!
!-------------------------------------------------------------------
!
   SUBROUTINE MixGeoms(Opts,Nams,GOpt,Convgd,GMLoc,iCLONE,iGEO)
     TYPE(Options)        :: Opts
     TYPE(FileNames)      :: Nams
     TYPE(GeomOpt)        :: GOpt
     TYPE(CRDS)           :: GMLoc
     INTEGER              :: InitGDIIS,NoGDIIS,GDIISOn,iGEO
     INTEGER              :: NConstr,NatmsLoc,NCart,iCLONE
     CHARACTER(LEN=DCL)   :: SCRPath
     INTEGER,DIMENSION(:) :: Convgd
     INTEGER              :: Print
     !
     InitGDIIS=GOpt%GDIIS%Init
     NoGDIIS  =GOpt%GDIIS%NoGDIIS
     GDIISOn  =GOpt%GDIIS%On     
     NConstr  =GOpt%Constr%NConstr
     NatmsLoc =GMLoc%Natms
     NCart    =3*NatmsLoc
     SCRPath  =TRIM(Nams%M_SCRATCH)//TRIM(Nams%SCF_NAME)// &
             '.'//TRIM(IntToChar(iCLONE))
     Print    =Opts%PFlags%GeOp
     !
     IF(Convgd(iCLONE)/=1) THEN
       CALL OPENAscii(OutFile,Out)
       IF((.NOT.NoGDIIS).AND.GDIISOn) THEN
         CALL GeoDIIS(GMLoc%AbCarts%D,GOpt,Nams%HFile,iCLONE, &
           iGEO,Print,SCRPath,InitGDIIS)
       ELSE
         IF(Print>=DEBUG_GEOP_MIN) THEN
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
END MODULE Optimizer
