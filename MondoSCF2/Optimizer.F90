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
  USE HessianMod
  USE SetXYZ
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
    INTEGER        :: NatmsLoc,IStart,gtmp,GBeg,GEnd
    TYPE(DBL_RNK2) :: OldXYZ
    LOGICAL        :: NewFile,ExitQ
    !-------------------------------------------------------------------------!
    IF(C%Opts%Grad==GRAD_TS_SEARCH_NEB)THEN
       GBeg=0
       GEnd=C%Geos%Clones+1
    ELSE	
       GBeg=1
       GEnd=C%Geos%Clones
    ENDIF
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
    IF(C%Opts%GeomPrint=='XSF')CALL XSFPreamble(0,C%Nams%GFile,Geo)
    DO iCLONE=GBeg,GEnd
       CALL PPrint(C%Geos%Clone(iCLONE),TRIM(C%Nams%GFile)//IntToChar(iCLONE),Geo,C%Opts%GeomPrint)
    ENDDO
    CALL MergePrintClones(C%Geos,C%Nams,C%Opts,GBeg,GEnd)
    iBAS=C%Sets%NBSets
    IStart=iGEO
    ! Follow the gradient downhill for NSteps
    DO iGEO=IStart,C%Opts%NSteps
       ! Compute new gradients
       CALL Force(iBAS,iGEO,C%Nams,C%Opts,C%Stat,C%Geos,C%Sets,C%MPIs)
       DO iCLONE=1,C%Geos%Clones
          C%Geos%Clone(iCLONE)%Displ%D=C%Geos%Clone(iCLONE)%AbCarts%D
       ENDDO

       ExitQ=SteepStep(iBAS,iGEO,Energy(:,iGEO),C)

       IF(C%Opts%Grad==GRAD_TS_SEARCH_NEB)THEN
          ! Overwrite the most recent GFile (doing transition state)
          CALL OpenASCII(C%Nams%GFile,Geo,NewFile_O=.TRUE.);CLOSE(Geo)
          IF(C%Opts%GeomPrint=='XSF') &
             CALL XSFPreamble(C%Geos%Clones+2,C%Nams%GFile,Geo)
          DO iCLONE=GBeg,GEnd
             ! If transition states, geometry is the image number
             ! otherwise, it is the step number 
             gtmp=C%Geos%Clone(iCLONE)%Confg
             C%Geos%Clone%Confg=iCLONE+1
            !CALL PPrint(C%Geos%Clone(iCLONE),C%Nams%GFile,Geo,C%Opts%GeomPrint)
             CALL PPrint(C%Geos%Clone(iCLONE),TRIM(C%Nams%GFile)//IntToChar(iCLONE),Geo,C%Opts%GeomPrint)
             C%Geos%Clone%Confg=gtmp
          ENDDO
       ELSE
          ! No transitions states, just good old downhill 
          IF(C%Opts%GeomPrint=='XSF') &
             CALL XSFPreamble(C%Geos%Clone(1)%Confg,C%Nams%GFile,Geo)
          DO iCLONE=1,C%Geos%Clones
             CALL PPrint(C%Geos%Clone(iCLONE),TRIM(C%Nams%GFile)//IntToChar(iCLONE),Geo,C%Opts%GeomPrint)
          ENDDO
       ENDIF
       IF(ExitQ)EXIT
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
       CALL MergePrintClones(C%Geos,C%Nams,C%Opts,GBeg,GEnd)
    ENDDO
  END SUBROUTINE SteepD
  !
  !=====================================================================================
  !
  SUBROUTINE MergePrintClones(Geos,Nams,Opts,GBeg,GEnd)
    TYPE(Geometries) :: Geos
    TYPE(FileNames)  :: Nams
    TYPE(Options)    :: Opts
    INTEGER          :: iCLONE,NatmsMerge,GBeg,GEnd,J,M
    TYPE(CRDS)       :: GMMerge
    REAL(DOUBLE),DIMENSION(3) :: CME1,CME2,TR,DIAG

    NatmsMerge=Geos%Clone(GBeg)%Natms*(GEnd-GBeg+1)
    CME1=Zero
    CME2=Zero
    DO J=1,Geos%Clone(GBeg)%Natms
      CME1=CME1+Geos%Clone(GBeg)%Carts%D(1:3,J)
      CME2=CME2+Geos%Clone(GEnd)%Carts%D(1:3,J)
    ENDDO
    CME1=CME1/DBLE(Geos%Clone(GBeg)%Natms)
    CME2=CME2/DBLE(Geos%Clone(GBeg)%Natms)
    TR=CME2-CME1
    TR=TR/SQRT(DOT_PRODUCT(TR,TR))
  
    ! find bounding box
    Geos%Clone(GBeg)%BndBox%D(1:3,1) = Geos%Clone(GBeg)%Carts%D(1:3,1)
    Geos%Clone(GBeg)%BndBox%D(1:3,2) = Geos%Clone(GBeg)%Carts%D(1:3,1)
    DO J = 2, Geos%Clone(GBeg)%Natms
      Geos%Clone(GBeg)%BndBox%D(1:3,1) = MIN(Geos%Clone(GBeg)%BndBox%D(1:3,1),Geos%Clone(GBeg)%Carts%D(1:3,J))
      Geos%Clone(GBeg)%BndBox%D(1:3,2) = MAX(Geos%Clone(GBeg)%BndBox%D(1:3,2),Geos%Clone(GBeg)%Carts%D(1:3,J))
    ENDDO
    DIAG=Geos%Clone(GBeg)%BndBox%D(1:3,2)-Geos%Clone(GBeg)%BndBox%D(1:3,1)
    TR=TR*SQRT(DOT_PRODUCT(DIAG,DIAG))
 
    GMMerge%Natms=NatmsMerge
    CALL New(GMMerge)
    M=0
    DO iCLONE=GBeg,GEnd
      DO J=1,Geos%Clone(iCLONE)%Natms
        M=M+1
        GMMerge%Carts%D(1:3,M)=Geos%Clone(iCLONE)%Carts%D(1:3,J)+TR*(iCLONE-1)
        GMMerge%AbCarts%D(1:3,M)=Geos%Clone(iCLONE)%AbCarts%D(1:3,J)+TR*(iCLONE-1)
        GMMerge%AtNum%D(M)=Geos%Clone(iCLONE)%AtNum%D(J)
        GMMerge%AtNam%C(M)=Geos%Clone(iCLONE)%AtNam%C(J)
      ENDDO
    ENDDO
    GMMerge%Confg=Geos%Clone(1)%Confg
    CALL PPrint(GMMerge,TRIM(Nams%GFile)//'M',Geo,Opts%GeomPrint)
    CALL Delete(GMMerge)
  END SUBROUTINE MergePrintClones
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
                                      -StepLength*C%Geos%Clone(iCLONE)%Gradients%D
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
    DO iCLONE=1,C%Geos%Clones
      CALL PPrint(C%Geos%Clone(iCLONE),C%Nams%GFile,Geo,C%Opts%GeomPrint)
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
    TYPE(DBL_RNK2)         :: A,AInvM,Carts,XYZAux
    TYPE(DBL_VECT)         :: V,DIISCo,GradI,GradJ
    REAL(DOUBLE)           :: DIISError,CoNo,Ratio,AMax
    INTEGER                :: cGEO,iCLONE,iGEO,mGEO,iDIIS,jDIIS,nDIIS
    INTEGER                :: iATS,I,J,K,Info
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
       CALL New(XYZAux,(/3,G%Clone(iCLONE)%NAtms/))
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
          CALL Get(XYZAux,'gradients',Tag_O=IntToChar(iDIIS))
          CALL CartRNK2ToCartRNK1(GradI%D,XYZAux%D)
          jDIIS=iDIIS
          DO J=I,nDIIS-1
             CALL Get(XYZAux,'gradients',Tag_O=IntToChar(jDIIS))
             CALL CartRNK2ToCartRNK1(GradJ%D,XYZAux%D)
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
       CALL FunkOnSqMat(nDIIS,Inverse,A%D,AInvM%D,CoNo_O=CoNo, &
                        Unit_O=6,PosDefMat_O=.FALSE.)
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
          G%Clone(iCLONE)%AbCarts%D=G%Clone(iCLONE)%AbCarts%D-5D-1*G%Clone(iCLONE)%Gradients%D          
       ENDIF
       ! Clean up a bit 
       CALL Delete(XYZAux)
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
    DO iCLONE=1,C%Geos%Clones
       Energies(iClone)=C%Geos%Clone(iClone)%ETotal
    ENDDO
    ! Temp Carts variable to allow backtracking of wrapped coordinates
    RMSGrad=Zero
    MAXGrad=Zero
    DO iCLONE=1,C%Geos%Clones
       CALL New(Carts(iCLONE),(/3,C%Geos%Clone(iCLONE)%NAtms/))
       Carts(iCLONE)%D=C%Geos%Clone(iCLONE)%AbCarts%D
       MAXGrad=MAX(MAXGrad,C%Geos%Clone(iCLONE)%GradMax)
       RMSGrad=MAX(RMSGrad,C%Geos%Clone(iCLONE)%GradRMS)
    ENDDO
!   Take some steps, more conservative if we are doing NEB ...
    IF(C%Opts%Grad==GRAD_TS_SEARCH_NEB)THEN
      !StepLength=0.5D0
       StepLength = C%Opts%RSL
       ! Take a step, any step
       DO iCLONE=1,C%Geos%Clones
          C%Geos%Clone(iCLONE)%AbCarts%D=Carts(iCLONE)%D-StepLength*C%Geos%Clone(iCLONE)%Gradients%D
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
!       GCnvrgd=RMSGrad<GTol(AL).AND.MaxGrad<GTol(AL)       
       GCnvrgd=MaxGrad<GTol(AL)       
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
             C%Geos%Clone(iCLONE)%AbCarts%D=Carts(iCLONE)%D-StepLength*C%Geos%Clone(iCLONE)%Gradients%D
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
    DO iCLONE=1,C%Geos%Clones
       Energies(iClone)=C%Geos%Clone(iClone)%ETotal
    ENDDO
  END FUNCTION SteepStep
!
!------------------------------------------------------------------
!
   SUBROUTINE IntOpt(C)
     TYPE(Controls)            :: C
     INTEGER                   :: I,iBAS,iGEO,iGEOst,iCLONE
     INTEGER                   :: AccL 
     INTEGER                   :: FirstGeom,NatmsLoc
     INTEGER                   :: ConvgdAll,MaxSteps,IStart
     TYPE(INT_VECT)            :: Convgd
     TYPE(INT_VECT)            :: BPrev,BCur
     !
     ! Set geometry optimization controls
     CALL SetGeOpCtrl(C%GOpt,C%Geos,C%Opts,C%Sets,C%Nams)
     ! initial geometry
     iGEO=C%Stat%Previous%I(3)
     iGEOst=iGEO
     MaxSteps=C%GOpt%GConvCrit%MaxGeOpSteps
     CALL NEW(BPrev,SIZE(C%Stat%Previous%I))
     CALL NEW(BCur ,SIZE(C%Stat%Previous%I))
     ! Build the guess 
     DO iBAS=1,C%Sets%NBSets-1
       CALL ReSetConnect(C%Geos)
       CALL ReDefIntCs(C%Geos,C%Opts,iGEO,iGEOst)
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
       CALL ReSetConnect(C%Geos)
       CALL ReDefIntCs(C%Geos,C%Opts,iGEO,iGEOst)
       CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)    
       CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
       !
       ! Calculate energy and force for all clones at once.
       !
      !IF(C%GOpt%RestartBas) THEN
      !  DO iBAS=1,C%Sets%NBSets-1
      !    CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)    
      !    CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
      !    CALL SCF(iBAS,iGEO,C)
      !  ENDDO
      !  iBAS=C%Sets%NBSets
      !  C%GOpt%RestartBas=.FALSE.
      !ENDIF
       BPrev%I=C%Stat%Previous%I
       BCur%I=C%Stat%Current%I
       CALL SCF(iBAS,iGEO,C)
       IF(iGEO>iGEOst) CALL BackTrack(iBAS,iGEO,C,BPrev%I,BCur%I)
       CALL Force(iBAS,iGEO,C%Nams,C%Opts,C%Stat, &
                  C%Geos,C%Sets,C%MPIs)
     ! CALL PrintClones(IGeo,C%Nams,C%Geos)
       !
       ! Loop over all clones and modify geometries.
       !
       ConvgdAll=1
       DO iCLONE=1,C%Geos%Clones
         Convgd%I=0
         CALL OptSingleMol(C%GOpt,C%Nams,C%Opts, &
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
     ! Convergence is reached at this point, print final energy
     ! and finish optimization.
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
     CALL Delete(BPrev)
     CALL Delete(BCur)
   END SUBROUTINE IntOpt
!
!------------------------------------------------------------------
!
   SUBROUTINE ModifyGeom(GOpt,XYZ,AtNum,IntCs,GradIn, &
                  Bond,AtmB,Convgd,ETot,PBC,iGEO,iCLONE, &
                  SCRPath,PWDPath,DoNEB,Print,HFileIn)
     TYPE(GeomOpt)               :: GOpt
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ,GradIn
     TYPE(PBCInfo)               :: PBC
     REAL(DOUBLE),DIMENSION(:)   :: AtNum
     TYPE(BONDDATA)              :: Bond
     TYPE(ATOMBONDS)             :: AtmB
     REAL(DOUBLE)                :: ETot
     INTEGER,DIMENSION(:)        :: Convgd
     INTEGER                     :: NIntC,NCart
     INTEGER                     :: NatmsLoc,iGEO,iCLONE
     TYPE(INTC)                  :: IntCs
     TYPE(DBL_VECT)              :: IntOld,CartGrad
     INTEGER                     :: Refresh
     CHARACTER(LEN=*)            :: SCRPath,HFileIn,PWDPath
     INTEGER                     :: Print
     LOGICAL                     :: DoNEB,Print2
     TYPE(TOPOLOGY)              :: TOPS
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     Print2= Print>=DEBUG_GEOP_MAX
     !
     ! Do we have to refresh internal coord defs?   
     !
     CALL IntCReDef(GOpt,Refresh,iGEO)
     !
     CALL New(CartGrad,NCart)
     CALL CartRNK2ToCartRNK1(CartGrad%D,GradIn)
     CALL CleanConstrGrad(CartGrad%D,GOpt%ExtIntCs)
     CALL TranslsOff(CartGrad%D,Print2)
     CALL RotationsOff(CartGrad%D,XYZ,Print2)
     !
     CALL GetCGradMax(CartGrad%D,NCart,GOpt%GOptStat%IMaxCGrad,&
                      GOpt%GOptStat%MaxCGrad)
     CALL TurnOnSelect(GOpt%GOptStat%MaxCGrad,GOpt%CoordCtrl%DoSelect)
     !
     ! Get internal coord defs.
     !
     IF(Refresh/=0) THEN
       CALL GetIntCs(XYZ,AtNum,IntCs,NIntC,Refresh,SCRPath, &
                     GOpt%CoordCtrl,GOpt%Constr,GOpt%GDIIS%MaxMem, &
                     GOpt%ExtIntCs,TOPS,Bond,AtmB,PBC, &
                     HFileIn_O=HFileIn,iCLONE_O=iCLONE,iGEO_O=iGEO)
       IF(NIntC==0) CALL Halt('Molecule has dissociated,'// &
                    'optimizer has not found any internal coordinates.')
     ENDIF
     IF(NIntC/=0) THEN
       CALL INTCValue(IntCs,XYZ, &
                      GOpt%CoordCtrl%LinCrit,GOpt%CoordCtrl%TorsLinCrit)
       CALL New(IntOld,NIntC)
       IntOld%D=IntCs%Value%D
     ENDIF
     !
     ! Get B matrices for redundancy projector, etc.
     !
   ! CALL GetMixedBMat(IntCs,XYZ,GOpt%TrfCtrl,GOpt%CoordCtrl,TOPS, &
   !                   GOpt%Constr%NCartConstr,Print,SCRPath,.TRUE.)
     CALL RefreshBMatInfo(IntCs,XYZ,GOpt%TrfCtrl, &
                          GOpt%CoordCtrl,Print,SCRPath,.TRUE.)
   ! CALL BuildUMatr(SCRPath,NCart)
     !
     ! Print current set of internals for debugging
     !
     IF(Print==DEBUG_GEOP_MAX) CALL PrtIntCoords(IntCs, &
       IntCs%Value%D,'Internals at step #'//TRIM(IntToChar(iGEO)))
     !
     CALL GradConv(GOpt%GOptStat,GOpt%GConvCrit,GOpt%Constr)
     !
     ! Calculate simple relaxation step from an inverse Hessian
     !
     IF(.NOT.GOpt%GOptStat%GeOpConvgd) THEN
       CALL RelaxGeom(GOpt,XYZ,AtNum,CartGrad%D,iCLONE,ETot, &
                  IntCs,iGEO,SCRPath,PWDPath,Print,HFileIn,Refresh)
     ELSE
       WRITE(*,200) iCLONE,iGEO
       WRITE(Out,200) iCLONE,iGEO
       200 FORMAT('     Geometry optimization of Clone #',I2,' converged in ',I3,' steps.')
       Convgd(iCLONE)=1
     ENDIF
     CALL GeOpReview(GOpt%Constr,GOpt%GOptStat,GOpt%CoordCtrl, &
                     GOpt%GConvCrit,XYZ,ETot,IntCs,IntOld, &
                     iCLONE,iGEO,DoNEB)
     CALL DoRestartBas(GOpt%RestartBas,GOpt%GOptStat)
     CALL TurnOnGDIIS(GOpt%GOptStat%MaxCGrad,GOpt%GDIIS%On)
     !
     ! tidy up
     !
     CALL Delete(CartGrad)
     CALL Delete(IntOld)
     CALL Delete(TOPS)
   END SUBROUTINE ModifyGeom
!
!--------------------------------------------------------------------
!
   SUBROUTINE GradConv(CtrlStat,GConvCr,CtrlConstr)
     TYPE(GOptStat)             :: CtrlStat
     TYPE(Constr)               :: CtrlConstr
     TYPE(GConvCrit)            :: GConvCr    
     !
     CtrlStat%GeOpConvgd=(CtrlStat%MaxCGrad<GConvCr%Grad) 
   END SUBROUTINE GradConv
!
!--------------------------------------------------------------------
!
   SUBROUTINE RelaxGeom(GOpt,XYZ,AtNum,CartGrad,iCLONE, &
                    ETot,IntCs,IGEO,SCRPath,PWDPath, &
                    Print,HFileIn,Refresh)
     TYPE(GeomOpt)               :: GOpt
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(:)   :: AtNum,CartGrad
     REAL(DOUBLE)                :: ETot     
     TYPE(INTC)                  :: IntCs
     INTEGER                     :: iGEO,iCLONE,Refresh,Print
     CHARACTER(LEN=*)            :: SCRPath,HFileIn,PWDPath 
     TYPE(DBL_RNK2)              :: SRStruct,RefGrad,RefStruct,SRDispl
     !
     ! Relax geometry
     !
     SELECT CASE(GOpt%Optimizer)
     CASE(GRAD_DIAGHESS_OPT)
       CALL RelaxDiagHess(GOpt,SCRPath,CartGrad,IntCs,AtNum, &
                          iGEO,XYZ,Print)
     CASE(GRAD_BiSect_OPT) 
       IF(iGEO<2) THEN
         CALL RelaxDiagHess(GOpt,SCRPath,CartGrad,IntCs,AtNum, &
                            iGEO,XYZ,Print)
       ELSE
         CALL RelaxBiSect(GOpt,SCRPath,PWDPath,HFileIn,CartGrad,IntCs, &
                          AtNum,iGEO,iCLONE,XYZ,Print)
       ENDIF
     END SELECT
   END SUBROUTINE RelaxGeom
!
!-------------------------------------------------------
!
   SUBROUTINE RelaxBiSect(GOpt,SCRPath,PWDPath,HFileIn,CartGrad,IntCs, &
                          AtNum,iGEO,iCLONE,XYZ,Print)
     TYPE(GeomOpt)               :: GOpt
     TYPE(INTC)                  :: IntCs
     REAL(DOUBLE),DIMENSION(:)   :: AtNum,CartGrad
     INTEGER                     :: iGEO,iCLONE,Print
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     CHARACTER(LEN=*)            :: SCRPath,PWDPath,HFileIn
     TYPE(DBL_VECT)              :: Displ,RefPoints,PredVals
     TYPE(DBL_RNK2)              :: SRStruct,RefStruct,RefGrad,SRDispl
     TYPE(DBL_RNK2)              :: IntCValues,IntCGrads,MixMat
     TYPE(INT_VECT)              :: ISpB,JSpB
     TYPE(DBL_VECT)              :: ASpB
     TYPE(DBL_RNK2)              :: FullB
     INTEGER                     :: I,J,NatmsLoc,NDim
     INTEGER                     :: HDFFileID,NCart,NMix
     LOGICAL                     :: DoMix
     !
     NDim=MIN(GOpt%GDIIS%MaxMem,iGEO)
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     DoMix=.FALSE.
     !
     ! calculate internal coord gradients
     !
     CALL CollectPast(SRStruct,RefStruct,RefGrad,SRDispl, &
                      HFileIn,NatmsLoc,NDim,iGEO,iCLONE)
     !
     CALL CollectINTCPast(RefStruct%D,RefGrad%D,IntCValues,IntCGrads, &
                          IntCs,GOpt,SCRPath,Print)
     CALL Delete(SRStruct)
     CALL Delete(RefStruct)
     CALL Delete(RefGrad)
     CALL Delete(SRDispl)
     !
     CALL New(RefPoints,IntCs%N)
     RefPoints%D=IntCValues%D(:,NDim)
     DO I=1,NDim    
       CALL SetBackToRefs(IntCValues%D(:,I),IntCs,RefPoints%D)
     ENDDO
     IF(DoMix) THEN
       CALL ReadBMATR(ISpB,JSpB,ASpB,TRIM(SCRPath)//'B')
       CALL Sp1x1ToFull(ISpB%I,JSpB%I,ASpB%D,IntCs%N,NCart,FullB)
       CALL Delete(ISpB) ; CALL Delete(JSpB) ; CALL Delete(ASpB)
       CALL UnitaryTR(IntCs,FullB%D,IntCValues%D,MixMat,NMix)
       CALL Delete(FullB)
      !CALL UnitaryTR(IntCs,IntCGrads%D,IntCValues%D,MixMat,NMix)
       CALL DisplFit(IntCs,IntCGrads%D,IntCValues%D,GOpt%Hessian, &
                  GOpt%CoordCtrl,PredVals,Displ,PWDPath,SCRPath,NCart, &
                   iGEO,MixMat_O=MixMat%D)
     ELSE
       CALL DisplFit(IntCs,IntCGrads%D,IntCValues%D,GOpt%Hessian, &
                  GOpt%CoordCtrl,PredVals,Displ,PWDPath,SCRPath,NCart, &
                  iGEO)
     ENDIF
     CALL Delete(IntCGrads)
     CALL Delete(IntCValues)
     !
     IF(DoMix) THEN
       CALL InternalToCart(XYZ,IntCs,PredVals%D,RefPoints%D,Print, &
                           GOpt%BackTrf,GOpt%TrfCtrl,GOpt%CoordCtrl,&
                           GOpt%Constr,SCRPath,Mixmat_O=MixMat%D)
     ELSE
       CALL InternalToCart(XYZ,IntCs,PredVals%D,RefPoints%D,Print, &
                           GOpt%BackTrf,GOpt%TrfCtrl,GOpt%CoordCtrl,&
                           GOpt%Constr,SCRPath)
     ENDIF
     ! 
     IF(DoMix) THEN
       CALL Delete(MixMat)
     ENDIF
     CALL Delete(PredVals)
     CALL Delete(RefPoints)
     CALL Delete(Displ)
   END SUBROUTINE RelaxBiSect
!
!-------------------------------------------------------
!
   SUBROUTINE RelaxDiagHess(GOpt,SCRPath,CartGrad,IntCs,AtNum,iGEO, &
                            XYZ,Print)
     TYPE(GeomOpt)               :: GOpt
     TYPE(INTC)                  :: IntCs
     REAL(DOUBLE),DIMENSION(:)   :: AtNum,CartGrad
     INTEGER                     :: iGEO,Print,I,J
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     CHARACTER(LEN=*)            :: SCRPath
     TYPE(DBL_VECT)              :: Grad,Displ,RefPoints
     !
     CALL New(Grad,IntCs%N)
     CALL New(Displ,IntCs%N)
     CALL New(RefPoints,IntCs%N)
     !
     ! calculate internal coord gradients
     !
     CALL CartToInternal(XYZ,IntCs,CartGrad,Grad%D,&
       GOpt%GrdTrf,GOpt%CoordCtrl,GOpt%TrfCtrl,Print,SCRPath)
     !CALL RedundancyOff(Grad%D,SCRPath,Print)
     !
     CALL GrdConvrgd(GOpt%GOptStat,IntCs,Grad%D)
     !
     CALL DiagHess(GOpt%CoordCtrl,GOpt%Hessian,Grad,Displ, &
                   IntCs,AtNum,iGEO,XYZ)
     !
     IntCs%PredVal%D=IntCs%Value%D+Displ%D
     CALL SetConstraints(IntCs,IntCs%PredVal%D)
     Displ%D=IntCs%PredVal%D-IntCs%Value%D
     !
     CALL CutOffDispl(Displ%D,IntCs, &
                      GOpt%CoordCtrl%MaxStre,GOpt%CoordCtrl%MaxAngle)
     !CALL RedundancyOff(Displ%D,SCRPath,Print)  
     ! 
     CALL INTCValue(IntCs,XYZ, &
                    GOpt%CoordCtrl%LinCrit,GOpt%CoordCtrl%TorsLinCrit)
     IntCs%Value%D=IntCs%Value%D+Displ%D
     RefPoints%D=IntCs%Value%D
     CALL InternalToCart(XYZ,IntCs,IntCs%Value%D,RefPoints%D,Print, &
                         GOpt%BackTrf,GOpt%TrfCtrl,GOpt%CoordCtrl,&
                         GOpt%Constr,SCRPath)  
     ! 
     CALL Delete(RefPoints)
     CALL Delete(Displ)
     CALL Delete(Grad)
   END SUBROUTINE RelaxDiagHess
!
!-------------------------------------------------------
!
   SUBROUTINE IntCReDef(GOpt,Refresh,iGEO)
     TYPE(GeomOpt)     :: GOpt
     INTEGER           :: Refresh,iGEO
     !
     !WARNING! refresh may change the number of internal coordinates!
     !
     Refresh=1
     IF(iGeo==1) THEN
       Refresh=1
       IF(GOpt%CoordCtrl%RefreshIn==4) Refresh=4
     ELSE
       Refresh=GOpt%CoordCtrl%RefreshIn
     ENDIF
   END SUBROUTINE IntCReDef
!
!------------------------------------------------------------------
!
   SUBROUTINE SteepestDesc(CoordC,Hess,Grad,Displ,XYZ,IntCs)
     TYPE(CoordCtrl)             :: CoordC
     TYPE(Hessian)               :: Hess  
     TYPE(DBL_VECT)              :: Grad,Displ 
     INTEGER                     :: NCart
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(INTC)                  :: IntCs
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
     !
     IntCs%PredGrad%D=Zero
     IntCs%PredVal%D=IntCs%Value%D+Displ%D
   END SUBROUTINE SteepestDesc
!
!-------------------------------------------------------
!
   SUBROUTINE GeOpReview(CtrlConstr,CtrlStat,CtrlCoord,GConvCr, &
                       XYZ,Etot,IntCs,IntOld,iCLONE,iGEO,DoNEB)
     TYPE(GOptStat)             :: CtrlStat
     TYPE(Constr)               :: CtrlConstr
     TYPE(CoordCtrl)            :: CtrlCoord  
     TYPE(GConvCrit)            :: GConvCr    
     REAL(DOUBLE),DIMENSION(:,:):: XYZ
     TYPE(INTC)                 :: IntCs
     TYPE(DBL_VECT)             :: IntOld,AuxVect
     INTEGER                    :: iCLONE,iGEO
     LOGICAL                    :: DoNEB
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
     REAL(DOUBLE)               :: Etot,Sum
     !                      
     INTEGER                    :: NStreGeOp,NBendGeOp,NLinBGeOp
     INTEGER                    :: NOutPGeOp,NTorsGeOp
     INTEGER                    :: MaxStre,MaxBend,MaxLinB
     INTEGER                    :: MaxOutP,MaxTors
     !                      
     INTEGER                    :: IMaxGrad,IMaxCGrad,IMaxGradNoConstr
     LOGICAL                    :: GradConv
     !
     NIntC=IntCs%N      
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
     IF(CtrlStat%GeOpConvgd) THEN
       WRITE(*,399) iCLONE,iGEO,ETot
       WRITE(Out,399) iCLONE,iGEO,ETot
       WRITE(*,140) MaxCGrad,(IMaxCGrad-1)/3+1
       WRITE(Out,140) MaxCGrad,(IMaxCGrad-1)/3+1
       RETURN
     ENDIF
     !
     NStreGeOp=CtrlCoord%NStre
     NBendGeOp=CtrlCoord%NBend
     NLinBGeOp=CtrlCoord%NLinB
     NOutPGeOp=CtrlCoord%NOutP
     NTorsGeOp=CtrlCoord%NTors
     !
     ! Size of internal coordinate changes
     !
     CALL INTCValue(IntCs,XYZ,CtrlCoord%LinCrit,CtrlCoord%TorsLinCrit)
     IntOld%D=IntCs%Value%D-IntOld%D
     CALL MapAngleDispl(IntCs,IntOld%D)
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
       IF(.NOT.IntCs%Active%L(I)) CYCLE
        IF(IntCs%Def%C(I)(1:4)=='STRE') THEN 
          Sum=ABS(IntOld%D(I))
          IF(MaxStreDispl<Sum) THEN
            MaxStre=I
            MaxStreDispl=Sum
          ENDIF
        ELSE IF(IntCs%Def%C(I)(1:4)=='BEND') THEN
          Sum=ABS(IntOld%D(I))
          IF(MaxBendDispl<Sum) THEN
            MaxBend=I
            MaxBendDispl=Sum
          ENDIF
        ELSE IF(IntCs%Def%C(I)(1:4)=='LINB') THEN
          Sum=ABS(IntOld%D(I))
          IF(MaxLinBDispl<Sum) THEN
            MaxLinB=I
            MaxLinBDispl=Sum
          ENDIF
        ELSE IF(IntCs%Def%C(I)(1:4)=='OUTP') THEN
          Sum=ABS(IntOld%D(I))
          IF(MaxOutPDispl<Sum) THEN
            MaxOutP=I
            MaxOutPDispl=Sum
          ENDIF
        ELSE IF(IntCs%Def%C(I)(1:4)=='TORS') THEN
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
     WRITE(*,410) MaxGrad,IntCs%Atoms%I(IMaxGrad,1:4)
     WRITE(*,140) MaxCGrad,IMaxCGrad
     WRITE(*,420) RMSGrad
     WRITE(Out,410) MaxGrad,IntCs%Atoms%I(IMaxGrad,1:4)
     WRITE(Out,140) MaxCGrad,IMaxCGrad
     WRITE(Out,420) RMSGrad
     IF(CtrlConstr%NConstr/=0) THEN
       WRITE(*,510) MaxGradNoConstr, &
                    IntCs%Atoms%I(IMaxGradNoConstr,1:4)
       WRITE(*,520) RMSGradNoConstr
       WRITE(Out,510) MaxGradNoConstr, &
                      IntCs%Atoms%I(IMaxGradNoConstr,1:4)
       WRITE(Out,520) RMSGradNoConstr
     ENDIF
     !
     IF(MaxStre/=0) THEN
       WRITE(*,430) MaxStreDispl,IntCs%Atoms%I(MaxStre,1:2)
       WRITE(Out,430) MaxStreDispl,IntCs%Atoms%I(MaxStre,1:2)
     ENDIF
     IF(MaxBend/=0) THEN
       WRITE(*,435) MaxBendDispl,IntCs%Atoms%I(MaxBend,1:3)
       WRITE(Out,435) MaxBendDispl,IntCs%Atoms%I(MaxBend,1:3)
     ENDIF
     IF(MaxLinB/=0) THEN
       WRITE(*,436) MaxLinBDispl,IntCs%Atoms%I(MaxLinB,1:3)
       WRITE(Out,436) MaxLinBDispl,IntCs%Atoms%I(MaxLinB,1:3)
     ENDIF
     IF(MaxOutP/=0) THEN
       WRITE(*,437) MaxOutPDispl,IntCs%Atoms%I(MaxOutP,1:4)
       WRITE(Out,437) MaxOutPDispl,IntCs%Atoms%I(MaxOutP,1:4)
     ENDIF
     IF(MaxTors/=0) THEN
       WRITE(*,438) MaxTorsDispl,IntCs%Atoms%I(MaxTors,1:4)
       WRITE(Out,438) MaxTorsDispl,IntCs%Atoms%I(MaxTors,1:4)
     ENDIF
     !
     WRITE(*,440) RMSIntDispl
     WRITE(Out,440) RMSIntDispl
     !
399 FORMAT('       Clone = ',I6,' GeOp step = ',I6,' Total Energy = ',F20.8)
400 FORMAT('Total Energy at Current Geometry = ',F20.8)
401 FORMAT('                    Total Energy = ',F20.8)
410 FORMAT('                   Max intl Grad = ',F12.6,' between atoms ',4I4)
140 FORMAT('     Max Unconstrained Cart Grad = ',F12.6,'      on atom  ',4I4)
420 FORMAT('                   RMS intl Grad = ',F12.6)
510 FORMAT('Max Grad on Unconstrained Coords = ',F12.6,' between atoms ',4I4)
520 FORMAT('RMS Grad on Unconstrained Coords = ',F12.6)
430 FORMAT('                  Max STRE Displ = ',F12.6,' between atoms ',4I4)
435 FORMAT('                  Max BEND Displ = ',F12.6,' between atoms ',4I4)
436 FORMAT('                  Max LINB Displ = ',F12.6,' between atoms ',4I4)
437 FORMAT('                  Max OUTP Displ = ',F12.6,' between atoms ',4I4)
438 FORMAT('                  Max TORS Displ = ',F12.6,' between atoms ',4I4)
440 FORMAT('                       RMS Displ = ',F12.6)
        !
   END SUBROUTINE GeOpReview
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
     GOpt%RestartBas=.FALSE.
     CALL SetCoordCtrl(GOpt%CoordCtrl)
     CALL   SetHessian(GOpt%Hessian)
     CALL SetGConvCrit(GOpt%GConvCrit,GOpt%Hessian,AccL,NatmsLoc)
     CALL     SetGDIIS(GOpt%GDIIS,GOpt%Optimizer)
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
     TYPE(DBL_RNK2)     :: XYZ,Vects
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
     CALL New(Vects,(/3,NatmsLoc/))
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
         Vects%D(3,NatmsLoc)=Geos%Clone(iCLONE)%Gradients%D(3,I)
       ENDDO
     ENDDO
     !
     FileName=TRIM(Nams%SCF_NAME)//'.Clones.xyz'
     !
     CALL PrtXYZ(Atnum%D,XYZ%D,FileName,&
                 'Step='//TRIM(IntToChar(IStep)),Vects_O=Vects%D)
     !
     CALL Delete(AtNum)
     CALL Delete(XYZ)
     CALL Delete(Vects)
   END SUBROUTINE PrintClones
!
!-------------------------------------------------------------------
!
   SUBROUTINE OptSingleMol(GOpt,Nams,Opts, &
                           GMLoc,Convgd,iGEO,iCLONE)
     TYPE(GeomOpt)        :: GOpt
     TYPE(FileNames)      :: Nams
     TYPE(Options)        :: Opts
     TYPE(CRDS)           :: GMLoc
     INTEGER,DIMENSION(:) :: Convgd
     INTEGER              :: I,iGEO,iCLONE,NatmsNew
     INTEGER              :: InitGDIIS,NConstr,NCart,NatmsLoc
     LOGICAL              :: NoGDIIS,GDIISOn,DoNEB
     CHARACTER(LEN=DCL)   :: SCRPath,PWDPath
     TYPE(DBL_RNK2)       :: XYZNew,GradNew
     TYPE(DBL_VECT)       :: AtNumNew,CartGrad
     !
     SCRPath  =TRIM(Nams%M_SCRATCH)//TRIM(Nams%SCF_NAME)// &
             '.'//TRIM(IntToChar(iCLONE))
     PWDPath  =TRIM(Nams%M_PWD)//TRIM(IntToChar(iCLONE))
     GMLoc%Displ%D=GMLoc%AbCarts%D
     DoNEB=(Opts%Grad==GRAD_TS_SEARCH_NEB)
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
         GradNew%D(1:3,NatmsNew)=GMLoc%Gradients%D(1:3,I)
         AtNumNew%D(NatmsNew)=GMLoc%AtNum%D(NatmsNew)
       ENDIF
     ENDDO
     !
     !--------------------------------------------
     CALL OpenASCII(OutFile,Out)
       CALL ModifyGeom(GOpt,XYZNew%D,AtNumNew%D,GMLoc%IntCs,GradNew%D, &
                       GMLoc%Bond,GMLoc%AtmB,Convgd,GMLoc%Etotal, &
                       GMLoc%PBC,iGEO,iCLONE,SCRPath, &
                       PWDPath,DoNEB,Opts%PFlags%GeOp,Nams%HFile)
     CLOSE(Out,STATUS='KEEP')
     !--------------------------------------------
     !
     ! Put back modified geometry into GMLoc array
     !
     NatmsNew=0
     DO I=1,GMLoc%Natms
       IF(GMLoc%CConstrain%I(I)/=2) THEN
         NatmsNew=NatmsNew+1
         GMLoc%Displ%D(1:3,I)=XYZNew%D(1:3,NatmsNew)
         GMLoc%Gradients%D(1:3,I)=GradNew%D(1:3,NatmsNew)
       ENDIF
     ENDDO
     !
     CALL Delete(XYZNew)
     CALL Delete(AtNumNew)
     CALL Delete(GradNew)
   END SUBROUTINE OptSingleMol
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
   ! GCrit=3.D-4
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
   SUBROUTINE SetGDIIS(GD,GOptimizer)
     TYPE(GDIIS)  :: GD
     INTEGER      :: GOptimizer
     !
     GD%Init    = 2
     IF(GOptimizer==GRAD_BiSect_OPT) THEN
       GD%MaxMem  = 7 
     ELSE
       GD%MaxMem  = 7
     ENDIF
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
     BackT%MaxIt_CooTrf = 50
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
    !IF(.NOT.TrfC%DoClssTrf) THEN
    !  TrfC%DoTranslOff=.FALSE.
    !  TrfC%DoRotOff=.FALSE.
    !ENDIF
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
     CoordC%TorsLinCrit =20.D0
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
   SUBROUTINE ReSetConnect(G)
     TYPE(Geometries)   :: G
     INTEGER            :: iCLONE
     !
     DO iCLONE=1,G%Clones
       IF(AllocQ(G%Clone(iCLONE)%Bond%Alloc)) THEN
         CALL Delete(G%Clone(iCLONE)%Bond)
       ELSE
         G%Clone(iCLONE)%Bond%N=0
       ENDIF
       IF(AllocQ(G%Clone(iCLONE)%AtmB%Alloc)) THEN
         CALL Delete(G%Clone(iCLONE)%AtmB)
       ELSE
         G%Clone(iCLONE)%AtmB%N1=0
         G%Clone(iCLONE)%AtmB%N2=0
       ENDIF
     ENDDO
   END SUBROUTINE ReSetConnect
!
!-------------------------------------------------------------------
!
   SUBROUTINE ReDefIntCs(G,O,iGEO,iGEOst)
     TYPE(Geometries)   :: G
     TYPE(Options)      :: O
     INTEGER            :: iCLONE,iGEO,iGEOst
     ! 
     IF(O%Guess==GUESS_EQ_RESTART.AND.iGEO==iGEOst) THEN
       RETURN
     ENDIF
     DO iCLONE=1,G%Clones
       IF(AllocQ(G%Clone(iCLONE)%IntCs%Alloc)) THEN
         CALL Delete(G%Clone(iCLONE)%IntCs)
       ENDIF
       G%Clone(iCLONE)%IntCs%N=IntCPerAtom*G%Clone(iCLONE)%Natms   
       CALL New(G%Clone(iCLONE)%IntCs,G%Clone(iCLONE)%IntCs%N)
       G%Clone(iCLONE)%IntCs%Def%C='BLANK'
     ENDDO
   END SUBROUTINE ReDefIntCs
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
       IF((.NOT.GOpt%GDIIS%NoGDIIS).AND.GOpt%GDIIS%On.AND.&
           iGEO>=GOpt%GDIIS%Init) THEN
         CALL GeoDIIS(GMLoc%AbCarts%D,GOpt%GDIIS,Nams%HFile,iCLONE, &
                      iGEO,Opts%PFlags%GeOp)
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
   SUBROUTINE CleanConstrGrad(CartGrad,IntCs_Extra)
     TYPE(INTC)                :: IntCs_Extra
     REAL(DOUBLE),DIMENSION(:) :: CartGrad
     INTEGER                   :: NIntC,I,JJ,J
     !
     NIntC=IntCs_Extra%N    
     DO I=1,NIntC
       IF(IntCs_Extra%Constraint%L(I)) THEN 
         JJ=IntCs_Extra%Atoms%I(I,1)
         IF(IntCs_Extra%Def%C(I)(1:5)=='CARTX') THEN
           J=3*(JJ-1)+1
           CartGrad(J)=Zero   
         ELSE IF(IntCs_Extra%Def%C(I)(1:5)=='CARTY') THEN
           J=3*(JJ-1)+2
           CartGrad(J)=Zero   
         ELSE IF(IntCs_Extra%Def%C(I)(1:5)=='CARTZ') THEN
           J=3*(JJ-1)+3
           CartGrad(J)=Zero   
         ENDIF
       ENDIF
     ENDDO
   END SUBROUTINE CleanConstrGrad
!
!-------------------------------------------------------------------
!
   SUBROUTINE GetCGradMax(CartGrad,NCart,IMaxCGrad,MaxCGrad)
     REAL(DOUBLE),DIMENSION(:) :: CartGrad
     REAL(DOUBLE),DIMENSION(3) :: Vect    
     REAL(DOUBLE)              :: MaxCGrad,Sum
     INTEGER                   :: Nat,NCart,I,I1,I2,IMaxCGrad
     !
     IMaxCGrad=1   
     MaxCGrad=Zero
     Nat=NCart/3
     DO I=1,Nat   
       I1=3*(I-1)+1
       I2=I1+2
       Vect=CartGrad(I1:I2)
       Sum=SQRT(DOT_PRODUCT(Vect,Vect))
       IF(Sum>MaxCGrad) THEN
         IMaxCGrad=I
          MaxCGrad=Sum
       ENDIF
     ENDDO
   END SUBROUTINE GetCGradMax
!
!-------------------------------------------------------------------
!
   SUBROUTINE BackTrack(iBAS,iGEO,C,BPrev,BCur,DoLineS_O)
     ! Go over clones and do backtracking whenever necessary
     TYPE(Controls)   :: C
     INTEGER          :: iBAS,iGEO,iCLONE
     INTEGER,DIMENSION(:):: BPrev,BCur
     INTEGER          :: NatmsLoc,NCart
     INTEGER          :: MaxBStep,IBStep
     CHARACTER(LEN=DCL):: chGEO
     TYPE(CRDS)       :: GMOld
     LOGICAL          :: DoBackTrack,DoLineS
     LOGICAL,OPTIONAL :: DoLineS_O
     REAL(DOUBLE)     :: EOld,ENew,MeanDist
     TYPE(DBL_VECT)   :: DistVect1,DistVect2
     !
     IF(.NOT.C%GOpt%GConvCrit%DoBackTr) THEN
       RETURN
     ENDIF
     !
     DoLineS=.FALSE.
     IF(PRESENT(DoLineS_O)) DoLineS=DoLineS_O
     IF(DoLineS) THEN
       CALL Force(iBAS,iGEO,C%Nams,C%Opts,C%Stat, &
                  C%Geos,C%Sets,C%MPIs)
     ENDIF
     !
     MaxBStep=1
     IF(DoLineS) MaxBStep=1
     !
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
           IF(DoLineS) THEN
             CALL LineSearch(C%Geos%Clone(iCLONE),GMOld)
           ELSE
             C%Geos%Clone(iCLONE)%Carts%D= &
               (C%Geos%Clone(iCLONE)%Carts%D+GMOld%Carts%D)*Half
             C%Geos%Clone(iCLONE)%AbCarts%D= &
               (C%Geos%Clone(iCLONE)%AbCarts%D+GMOld%AbCarts%D)*Half
           ENDIF
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
         C%Stat%Previous%I=BPrev
         C%Stat%Current%I=BCur
        !C%Stat%Previous%I(3)=iGEO-1
        !C%Stat%Current%I(3)=iGEO
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
   SUBROUTINE LineSearch(GMNew,GMOld)
     TYPE(CRDS)     :: GMNew,GMOld 
     INTEGER        :: I,J,NCart,II
     TYPE(DBL_VECT) :: DeltaX,Grad,PVect,VAux,XV1,XV2,PGrad
     TYPE(DBL_RNK2) :: Mat,Aux,InvMat
     REAL(DOUBLE)   :: GO,GN,DX2,DE,DX,DX3,DX4,X1,X2,X3,E0,E1,E2,E3,D0
     REAL(DOUBLE)   :: DE2DX21,DE2DX22,Fact
     LOGICAL        :: DoHalve
     !
     NCart=3*GMNew%Natms
     CALL New(DeltaX,NCart)
     CALL New(Grad,NCart)
     CALL New(XV1,NCart)
     CALL New(XV2,NCart)
     CALL New(Mat,(/4,4/))
     CALL New(Aux,(/4,4/))
     CALL New(InvMat,(/4,4/))
     CALL New(PVect,4)
     CALL New(PGrad,4)
     CALL New(VAux,4)
     !
     CALL CartRNK2ToCartRNK1(XV1%D,GMOld%AbCarts%D)
     CALL CartRNK2ToCartRNK1(XV2%D,GMNew%AbCarts%D)
     DeltaX%D=XV2%D-XV1%D
     DX2=DOT_PRODUCT(DeltaX%D,DeltaX%D)
     DX=SQRT(DX2)
     DX3=DX*DX2
     DX4=DX2*DX2
     !
     CALL CartRNK2ToCartRNK1(Grad%D,GMNew%Gradients%D)
     GN=DOT_PRODUCT(Grad%D,DeltaX%D)/DX
     CALL CartRNK2ToCartRNK1(Grad%D,GMOld%Gradients%D)
     GO=DOT_PRODUCT(Grad%D,DeltaX%D)/DX
     DE=GMNew%ETotal-GMOld%ETotal
     !
     Mat%D(1,1:4)=(/DX,DX2,DX3,DX4/)
     Mat%D(2,1:4)=(/One,Zero,Zero,Zero/)
     Mat%D(3,1:4)=(/One,Two*DX,3.D0*DX2,4.D0*DX3/)
     Mat%D(4,1:4)=(/Zero,Two,Zero,Zero/)
     DO II=1,2
       PVect%D(1)=DE
       PVect%D(2)=GO
       PVect%D(3)=GN
       PVect%D(4)=Zero
       !
       CALL DGEMM_TNc(4,4,4,One,Zero,Mat%D,Mat%D,Aux%D)
       CALL DIISInvMat(Aux%D,InvMat%D,'Basic')
       CALL DGEMM_TNc(4,4,1,One,Zero,Mat%D,PVect%D,VAux%D)
       CALL DGEMM_NNc(4,4,1,One,Zero,InvMat%D,VAux%D,PVect%D)
       !
       PGrad%D(1)=PVect%D(1)
       PGrad%D(2)=2.D0*PVect%D(2)
       PGrad%D(3)=3.D0*PVect%D(3)
       PGrad%D(4)=4.D0*PVect%D(4)
       CALL CubicRoots(PGrad%D(1),PGrad%D(2),PGrad%D(3),PGrad%D(4), &
                       X1,X2,X3)
       !
       ! Check second derivatives
       !
       DE2DX21=2.D0*PVect%D(2)
       DE2DX22=2.D0*PVect%D(2)+6.D0*PVect%D(3)*DX+12.D0*PVect%D(4)*DX2
       IF(DE2DX21>-1.D-10.AND.DE2DX22>-1.D-10) THEN
         EXIT
       ELSE IF(II<2) THEN
         Mat%D(4,1:4)=(/Zero,Two,6.D0*DX,12.D0*DX2/)
       ELSE
         X1=-One ; X2=-One ; X3=-One
       ENDIF
     ENDDO
     !
     E1=GMOld%Etotal+PVect%D(1)*X1+PVect%D(2)*X1**2+ &
        PVect%D(3)*X1**3+PVect%D(4)*X1**4
     E2=GMOld%Etotal+PVect%D(1)*X2+PVect%D(2)*X2**2+ &
        PVect%D(3)*X2**3+PVect%D(4)*X2**4
     E3=GMOld%Etotal+PVect%D(1)*X3+PVect%D(2)*X3**2+ &
        PVect%D(3)*X3**3+PVect%D(4)*X3**4
     !
     DoHalve=.TRUE.
     E0=E1
     IF(X1>Zero.AND.X1<DX) THEN
       DoHalve=.FALSE. 
       D0=X1
       E0=E1
     ENDIF
     IF((X2>Zero.AND.X2<DX).AND.E2<E0) THEN
       DoHalve=.FALSE. 
       D0=X2
       E0=E2
     ENDIF
     IF((X3>Zero.AND.X3<DX).AND.E3<E0) THEN
       DoHalve=.FALSE. 
       D0=X3
       E0=E3
     ENDIF
     !
     IF(DoHalve) THEN
       GMNew%AbCarts%D=Half*(GMNew%AbCarts%D+GMOld%AbCarts%D)
       GMNew%Carts%D=Half*(GMNew%Carts%D+GMOld%Carts%D)
     ELSE
       Fact=D0/DX
       GMNew%AbCarts%D=GMOld%AbCarts%D+ &
                       Fact*(GMNew%AbCarts%D-GMOld%AbCarts%D)
       GMNew%Carts%D  =GMOld%Carts%D+   &
                       Fact*(GMNew%Carts%D-GMOld%Carts%D)
      !GMNew%ETotal=E0
      !GMNew%Gradients%D=GMOld%Gradients%D+ &
      !                Fact*(GMNew%Gradients%D-GMOld%Gradients%D)
     ENDIF
     !
     CALL Delete(Grad)
     CALL Delete(DeltaX)
     CALL Delete(Mat)
     CALL Delete(Aux)
     CALL Delete(InvMat)
     CALL Delete(PGrad)
     CALL Delete(PVect)
     CALL Delete(VAux)
     CALL Delete(XV1)
     CALL Delete(XV2)
   END SUBROUTINE LineSearch
!
!-------------------------------------------------------------------
!
   SUBROUTINE TurnOnGDIIS(MaxCGrad,On)
     REAL(DOUBLE)  :: MaxCGrad
     LOGICAL       :: On
     !
     On=.FALSE.
     IF(MaxCGrad<0.001D0) THEN
       On=.TRUE.
     ENDIF
   END SUBROUTINE TurnOnGDIIS
!
!-------------------------------------------------------------------
!
   SUBROUTINE TurnOnSelect(MaxCGrad,DoSelect)
     REAL(DOUBLE)  :: MaxCGrad
     LOGICAL       :: DoSelect
     !
    !DoSelect=.FALSE.
     DoSelect=.TRUE.
     IF(MaxCGrad<0.070D0) THEN
       DoSelect=.TRUE.
     ENDIF
   END SUBROUTINE TurnOnSelect
!
!-------------------------------------------------------------------
!
   SUBROUTINE RescaleGrad(Grad,Print)
     REAL(DOUBLE),DIMENSION(:) :: Grad
     REAL(DOUBLE)              :: MaxGrad,SetMax,MaxGradP,MaxGradN,Fact
     INTEGER                   :: Print
     LOGICAL                   :: Print2
     !
     Print2= Print>=DEBUG_GEOP_MAX
     SetMax=0.050D0
     MaxGradP=MAXVAL(Grad)
     MaxGradN=MAXVAL(-Grad)
     MaxGrad=MAX(MaxGradP,MaxGradN)
     IF(MaxGrad>SetMax) THEN
       Fact=SetMax/MaxGrad
       Grad=Fact*Grad
       IF(Print2) THEN
         WRITE(*,100) MaxGrad,SetMax
         WRITE(Out,100) MaxGrad,SetMax
         100 FORMAT('Internal Gradients have been rescaled from ',F10.5,' to ',F10.5)
       ENDIF
      !CALL TranslsOff(Grad,Print2)
      !CALL RotationsOff(Grad,XYZ,Print2)
     ENDIF
   END SUBROUTINE RescaleGrad
!
!-------------------------------------------------------------------
!
   SUBROUTINE DoRestartBas(RestartBas,GStat)
     LOGICAL        :: RestartBas
     TYPE(GOptStat) :: GStat
     REAL(DOUBLE)   :: AngleCrit
     !
     AngleCrit=3.D0/180.D0*PI
     IF(GStat%MaxStreDispl>0.03.OR. &
        GStat%MaxBendDispl>AngleCrit.OR. &
        GStat%MaxLinBDispl>AngleCrit.OR. &
        GStat%MaxTorsDispl>AngleCrit.OR. &
        GStat%MaxOutPDispl>AngleCrit) THEN
       RestartBas=.TRUE.
     ENDIF  
   END SUBROUTINE DoRestartBas
!
!-------------------------------------------------------------------
!
   SUBROUTINE PrepBiSect(Grad,IntCs,Displ)
     INTEGER                     :: I,NIntC
     REAL(DOUBLE)                :: DStre,DAngle,DTors,DOutP
     REAL(DOUBLE)                :: DPhase,Conv,D,Tol
     TYPE(INTC)                  :: IntCs
     REAL(DOUBLE),DIMENSION(:)   :: Grad
     TYPE(DBL_VECT)              :: Displ
     !
     Conv=PI/180.D0
     NIntC=IntCs%N     
     Tol=1.D-8
     !
     DStre=0.010D0 ! in atomic units 
     DAngle=7.D0*Conv ! in radians      
     DTors=7.D0*Conv
     DOutP=7.D0*Conv
     DO I=1,NIntC
       IF(IntCs%Def%C(I)(1:4)=='STRE') THEN
       ! D=RANDOM_DBL((/-DStre,DStre/))
       ! D=DStre
         D=Grad(I)/0.5D0
       ELSE IF(IntCs%Def%C(I)(1:4)=='TORS') THEN
       ! D=RANDOM_DBL((/-DTors,DTors/))
       ! D=DTors
         D=Grad(I)/0.1D0
       ELSE IF(IntCs%Def%C(I)(1:4)=='OUTP') THEN
       ! D=RANDOM_DBL((/-DOutP,DOutP/))
       ! D=DOutP
         D=Grad(I)/0.2D0
       ELSE IF(IntCs%Def%C(I)(1:4)=='BEND') THEN
       ! D=RANDOM_DBL((/-DAngle,DAngle/))
       ! D=DAngle
         D=Grad(I)/0.2D0
       ELSE IF(IntCs%Def%C(I)(1:4)=='LINB') THEN
       ! D=RANDOM_DBL((/-DAngle,DAngle/))
       ! D=DAngle
         D=Grad(I)/0.2D0
       ELSE
       ! D=RANDOM_DBL((/-DStre,DStre/))
       ! D=DStre
         D=Grad(I)/0.5D0
       ENDIF
       IF(ABS(Grad(I))>Tol) THEN
         Displ%D(I)=SIGN(ABS(D),-Grad(I))
       ELSE
         Displ%D(I)=D
       ENDIF
     ENDDO
   END SUBROUTINE PrepBiSect
!
!-------------------------------------------------------------------
!
END MODULE Optimizer
