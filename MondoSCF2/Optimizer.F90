MODULE Optimizer
  USE SCFs
  USE InOut
  USE MatFunk
  USE AtomPairs
  USE Thresholds
  USE ControlStructures  
  IMPLICIT NONE
CONTAINS
  !=====================================================================================
  !
  !=====================================================================================  
  SUBROUTINE Descender(C)
    TYPE(Controls) :: C
    !----------------------------------------------------------------------------------!
    ! Follow Cartesian gradient down hill
!    CALL SteepD(C)
    ! Follow extrapolated Cartesian gradient down hill
    CALL GDicer(C)
  END SUBROUTINE Descender
  !=====================================================================================
  !
  !=====================================================================================  
  SUBROUTINE SteepD(C)
    TYPE(Controls) :: C
    INTEGER        :: iBAS,iGEO
    !----------------------------------------------------------------------------------!
    ! Start with the first geometry
    iGEO=1
    ! Initialize the previous state
    C%Stat%Previous%I=(/0,1,1/)
    ! Archive MPIs and Geos 
    CALL InitClones(C%Nams,C%Geos)
    CALL MPIsArchive(C%Nams,C%Geos,C%MPIs)
    CALL GeomArchive(iGEO,C%Nams,C%Geos)    
    ! Build the guess 
    DO iBAS=1,C%Sets%NBSets
       CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
       CALL SCF(iBAS,iGEO,C)
    ENDDO
    ! Follow the gradient down hill
    iBAS=C%Sets%NBSets
    DO iGEO=1,C%Opts%NSteps
       ! Compute new gradients
       CALL Force(iBAS,iGEO,C%Nams,C%Opts,C%Stat,C%Geos,C%MPIs)
       IF(SteepStep(iBAS,iGEO,C))EXIT
    ENDDO
  END SUBROUTINE SteepD
  !=====================================================================================
  !
  !=====================================================================================  
  SUBROUTINE GDicer(C)
    TYPE(Controls)         :: C
    TYPE(DBL_VECT)         :: GradMax,GradRMS,Grad
    REAL(DOUBLE)           :: DIISErr,GRMSQ,GMAXQ
    REAL(DOUBLE),PARAMETER :: StepLength=3D-1 ! Open issues about this and normalization in GDIIS
    INTEGER                :: iBAS,iGEO,iCLONE,AccL    
    INTEGER                :: Relaxations=3   ! This should be an input variable at some point
    LOGICAL                :: Converged
    CHARACTER(LEN=DCL)     :: Mssg
    !----------------------------------------------------------------------------------!
    ! Start with the first geometry
    iGEO=1
    ! Initialize the previous state
    C%Stat%Previous%I=(/0,1,1/)
    ! Archive MPIs and Geos 
    CALL InitClones(C%Nams,C%Geos)
    CALL MPIsArchive(C%Nams,C%Geos,C%MPIs)
    CALL GeomArchive(iGEO,C%Nams,C%Geos)    
    ! Build the guess 
    DO iBAS=1,C%Sets%NBSets
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
    DO iGEO=1,25  !C%Opts%NSteps
       ! Compute new gradients
       CALL Force(iBAS,iGEO,C%Nams,C%Opts,C%Stat,C%Geos,C%MPIs)       
       ! Convergence statistics for the gradient
       GradMax%D(iGEO)=Zero
       GradRMS%D(iGEO)=Zero
       DO iCLONE=1,C%Geos%Clones
          GradMAX%D(iGEO)=MAX(GradMax%D(iGEO),C%Geos%Clone(iCLONE)%GradMax)
          GradRMS%D(iGEO)=MAX(GradRMS%D(iGEO),C%Geos%Clone(iCLONE)%GradRMS)
       ENDDO
       ! Gradient quotients
       GRMSQ=ABS(GradMax%D(iGEO)-GradMax%D(iGEO-1))/ABS(GradMax%D(iGEO)+1.D-50)    
       GMAXQ=ABS(GradRMS%D(iGEO)-GradRMS%D(iGEO-1))/ABS(GradRMS%D(iGEO)+1.D-50)    
       Mssg='Gmax='//TRIM(DblToShrtChar(GradMAX%D(iGEO)))//', Grms='//TRIM(DblToShrtChar(GradRMS%D(iGEO)))
       ! Go downhill by following the gradient or with GDIIS
       IF(iGEO>Relaxations)THEN
          ! GDIIS extrapolation 
          CALL ForceDStep(iBAS,iGEO,C%Nams,C%Geos,DIISErr)
          Mssg=TRIM(Mssg)//', Ediis='//TRIM(DblToShrtChar(DIISErr))
          ! Check for absolute convergence
          Converged=.FALSE.
          IF(GradMAX%D(iGEO)<GTol(AccL).AND.GradRMS%D(iGEO)<GTol(AccL))THEN
             Converged=.TRUE.
             Mssg=ProcessName('GDicer',' Converged!')//TRIM(Mssg)     
          ELSEIF(GRMSQ<2D-1.AND.GMAXQ<2.D-1.AND.DIISerr<1D-10)THEN  ! Look for non-decreasing errors (stall out)
             ! This part needs more work!!!
!             Converged=.TRUE.
             Mssg=ProcessName('GDicer',' Stalled ')//TRIM(Mssg)     
          ELSE
             Mssg=ProcessName('GDicer',' GDIIS # '//TRIM(IntToChar(iGEO)))//TRIM(Mssg)          
          ENDIF
          CALL OpenASCII(C%Nams%OFile,Out)
          WRITE(Out,*)TRIM(Mssg)
          CLOSE(Out,STATUS='KEEP')
          IF(Converged)EXIT
       ELSE
!          HDFFileID=OpenHDF(C%Nams%HFile)
          ! Take a few small steps along the gradient to start
          DO iCLONE=1,C%Geos%Clones
             C%Geos%Clone(iCLONE)%Carts%D=C%Geos%Clone(iCLONE)%Carts%D &
                  -StepLength*C%Geos%Clone(iCLONE)%Vects%D
             !          CALL WrapAtoms(c%Geos%Clone(iCLONE))
             ! Rescale the gradient by the step for use in GDIIS
!             HDF_CurrentID=OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
!             CALL New(Grad,3*C%Geos%Clone(iCLONE)%NAtms)
!             CALL Get(Grad,'grade',Tag_O=IntToChar(iGEO))
!             Grad%D=Grad%D*StepLength
!             CALL Put(Grad,'grade',Tag_O=IntToChar(iGEO))
!             CALL Delete(Grad)
!             CALL CloseHDFGroup(HDF_CurrentID)
          ENDDO
!          CALL CloseHDF(HDFFileID)
          CALL GeomArchive(iGEO+1,C%Nams,C%Geos)    
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
  SUBROUTINE ForceDStep(cBAS,cGEO,N,G,DIISError)
    TYPE(FileNames)        :: N
    TYPE(Geometries)       :: G
    TYPE(DBL_RNK2)         :: A,AInv,Carts
    TYPE(DBL_VECT)         :: V,DIISCo,GradI,GradJ
    REAL(DOUBLE)           :: DIISError,CoNo
    INTEGER                :: cBAS,cGEO,iCLONE,iGEO,mGEO,iDIIS,jDIIS,nDIIS,iATS,I,J,K
    INTEGER,     PARAMETER :: MaxCoef=5
    REAL(DOUBLE),PARAMETER :: MaxCoNo=1D10
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
       CALL New(AInv,(/nDIIS,nDIIS/))
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
       AInv%D=Zero
       ! Solve the inverse least squares problem
       CALL FunkOnSqMat(nDIIS,Inverse,A%D,AInv%D,CoNo_O=CoNo)
       CALL DGEMV('N',nDIIS,nDIIS,One,AInv%D,nDIIS,V%D,1,Zero,DIISCo%D,1)
       ! Here is the DIIS error (resdiual)       
       DIISError=MAX(DIISError,ABS(DIISCo%D(nDIIS)))
!       WRITE(*,*)' CoNo = ',DblToChar(CoNo)
!       WRITE(*,*)' DIISCo = ',DIISCo%D
       IF(CoNo>MaxCoNo)THEN
          ! Clean up ...
          CALL Delete(V)
          CALL Delete(A)
          CALL Delete(AInv)
          CALL Delete(DIISCo)
          CALL UnSetDSYEVWork()
          ! ... and start over with a smaller subspace
          nDIIS=nDIIS-1
          ! Sometimes ya just gota have a goto
          GOTO 1
       ELSE
          ! Extrapolate
          iDIIS=1
          G%Clone(iCLONE)%Carts%D=Zero
          G%Clone(iCLONE)%AbCarts%D=Zero
          CALL New(Carts,(/3,G%Clone(iCLONE)%NAtms/))
          ! DIIS extrapolation of the position
          DO iGEO=mGEO,cGEO
             ! These are the unwrapped Cartesian coordinates
             !          CALL Get(Carts,'abcarts',Tag_O=IntToChar(iGEO))
             CALL Get(Carts,'cartesians',Tag_O=IntToChar(iGEO))
             !          G%Clone(iCLONE)%AbCarts%D=G%Clone(iCLONE)%AbCarts%D+DIISCo%D(iDIIS)*Carts%D
             G%Clone(iCLONE)%Carts%D=G%Clone(iCLONE)%Carts%D+DIISCo%D(iDIIS)*Carts%D
             iDIIS=iDIIS+1
          ENDDO
          !       CALL WrapAtoms(G%Clone(iCLONE))
          CALL Delete(Carts)
       ENDIF
       ! Clean up a bit ...
       CALL Delete(GradI)
       CALL Delete(GradJ)
    ENDDO
    ! And clean up some more
    CALL Delete(V)
    CALL Delete(A)
    CALL Delete(AInv)
    CALL Delete(DIISCo)
    CALL UnSetDSYEVWork()
    ! Archive geometries
    CALL GeomArchive(cGEO+1,N,G)    
  END SUBROUTINE ForceDStep
  !=====================================================================================
  !
  !=====================================================================================  
  FUNCTION SteepStep(cBAS,cGEO,C) RESULT(Converged)
    TYPE(Controls)                           :: C
    TYPE(DBL_RNK2), DIMENSION(C%Geos%Clones) :: Carts
    REAL(DOUBLE),   DIMENSION(C%Geos%Clones) :: Energies
    REAL(DOUBLE)                             :: StepLength,RelErrE,MAXGrad,RMSGrad, &
         ETest,GTest
    INTEGER                                  :: cBAS,cGEO,iSTEP,iCLONE,iATS,AL,K
    INTEGER, PARAMETER                       :: MaxSTEP=6
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
!       Carts(iCLONE)%D=C%Geos%Clone(iCLONE)%abCarts%D
       Carts(iCLONE)%D=C%Geos%Clone(iCLONE)%Carts%D
       MAXGrad=MAX(MAXGrad,C%Geos%Clone(iCLONE)%GradMax)
       RMSGrad=MAX(RMSGrad,C%Geos%Clone(iCLONE)%GradRMS)
    ENDDO
    WRITE(*,*)'-------------------------------------------------------------- '
    WRITE(*,*)' '
    WRITE(*,*)' '
    WRITE(*,*)cGEO,' CurrentSt = ',RMSGrad,MAXGrad
    WRITE(*,*)' '
    WRITE(*,*)' '
    WRITE(*,*)'-------------------------------------------------------------- '
    ! Take some steps 
    StepLength=One
    DO iSTEP=1,MaxSTEP
       StepLength=StepLength/Two
       ! Step the absolute positions, and wrap to get the Carts array
       DO iCLONE=1,C%Geos%Clones
!          C%Geos%Clone(iCLONE)%abCarts%D=Carts(iCLONE)%D-StepLength*C%Geos%Clone(iCLONE)%Vects%D
          C%Geos%Clone(iCLONE)%Carts%D=Carts(iCLONE)%D-StepLength*C%Geos%Clone(iCLONE)%Vects%D
          !          CALL WrapAtoms(C%Geos%Clone(iCLONE))
       ENDDO
       ! Archive geometries
       CALL GeomArchive(cGEO+1,C%Nams,C%Geos)    
       CALL PPrint(C%Geos%Clone(1),Unit_O=6)
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
!          XCnvrgd=RMSDisp<XTol(AL).AND.MaxDisp<XTol(AL)          
       ELSE
          ECnvrgd=.FALSE.
!          XCnvrgd=.FALSE.
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
          WRITE(*,*)TRIM(Mssg)             
          CALL OpenASCII(OutFile,Out)
          WRITE(Out,*)TRIM(Mssg)             
          CLOSE(Out)
       ENDIF
    ENDDO
    Mssg=TRIM(Mssg)//' dE= '//TRIM(DblToShrtChar(RelErrE)) &
         //', Grms= '//TRIM(DblToShrtChar(RMSGrad))        & 
         //', Gmax= '//TRIM(DblToShrtChar(MAXGrad))        &
         //', Step= '//TRIM(DblToShrtChar(StepLength))
    WRITE(*,*)TRIM(Mssg)             
    CALL OpenASCII(OutFile,Out)
    WRITE(Out,*)TRIM(Mssg)             
    CLOSE(Out)
    ! Clean up
    DO iCLONE=1,C%Geos%Clones
       CALL Delete(Carts(iCLONE))
    ENDDO
    ! Keep geometry
    CALL GeomArchive(cGEO,C%Nams,C%Geos)    
    ! Print something 
  END FUNCTION SteepStep
END MODULE Optimizer
