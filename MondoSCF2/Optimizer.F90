MODULE Optimizer
  USE SCFs
  USE InOut
  USE MatFunk
  USE Numerics
  USE AtomPairs
  USE ControlStructures  
  USE InCoords
  USE QUICCAMod
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
     TYPE(INTC)                :: IntCES
     !
     iGEO=C%Stat%Previous%I(3)
     iGEOst=iGEO
     ! Set geometry optimization controls
     CALL SetGeOpCtrl(C%GOpt,C%Geos,C%Opts,C%Sets,C%Nams,C%PBCs,iGEO)
     ! initial geometry
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
   ! CALL New(IntCES,C%GOpt%ExtIntCs%N)
   ! CALL SetEq(IntCES,C%GOpt%ExtIntCs,1,C%GOpt%ExtIntCs%N,1)
     C%GOpt%GConvCrit%DoLattStep=.TRUE.
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
     ! CALL AdjustLatt(iGEO,C%GOpt%ExtIntCs,IntCES,C%Geos%Clone(1)%PBC%Dimen)
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
     ! !
     ! ! Project out constraints and put them onto disk
     ! !
     ! CALL ProjectConstr(C%Geos)
     ! CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)    
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
   ! CALL Delete(IntCES)
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
   SUBROUTINE ModifyGeom(GOpt,XYZ,RefXYZ,AtNum,IntCs,GradIn, &
                  Bond,AtmB,Convgd,ETot,PBCDim,iGEO,iCLONE, &
                  SCRPath,PWDPath,DoNEB,Print,HFileIn,PBCFit)
     TYPE(GeomOpt)               :: GOpt
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ,GradIn,RefXYZ
     REAL(DOUBLE),DIMENSION(:)   :: AtNum
     TYPE(BONDDATA)              :: Bond
     TYPE(ATOMBONDS)             :: AtmB
     REAL(DOUBLE)                :: ETot
     INTEGER,DIMENSION(:)        :: Convgd
     TYPE(PBCFits)               :: PBCFit
     INTEGER                     :: NCart,K,I,J
     INTEGER                     :: NatmsLoc,iGEO,iCLONE
     TYPE(INTC)                  :: IntCs,IntCL
     TYPE(DBL_VECT)              :: IntOld,CartGrad,Carts,LatOld
     INTEGER                     :: Refresh
     CHARACTER(LEN=*)            :: SCRPath,HFileIn,PWDPath
     INTEGER                     :: Print,PBCDim
     LOGICAL                     :: DoNEB,Print2
     TYPE(TOPOLOGY)              :: TOPS
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     Print2= Print>=DEBUG_GEOP_MAX
     !
     !----------------------------------------------------------------
     !
     CALL LatticeINTC(IntCL,PBCDim,DoVolume_O=.TRUE.)
     CALL New(Gopt%LattIntC%Grad,IntCL%N)
     CALL New(Gopt%LattIntC%Displ,IntCL%N)
     CALL New(CartGrad,NCart)
     !
     ! Project out constraints, translations and rotations
     ! WARNING! For the case of water with a constrained H and H-O-H
     ! angle, the order of projections was very important: 
     ! first constraints, then the rest.
     !
     CALL CartRNK2ToCartRNK1(CartGrad%D,GradIn)
     CALL GetLattGrads(IntCL,CartGrad%D, &
                       XYZ,Gopt%LattIntC%Grad%D,PBCDim)
   ! CALL SYMMSTRESS(CartGrad%D(NCart-8),XYZ,PBCDim)
     CALL CleanConstrCart(XYZ,PBCDim,CartGrad%D,GOpt,SCRPath)
     CALL New(Carts,NCart)
     CALL CartRNK2ToCartRNK1(Carts%D,XYZ)
     IF(GOpt%TrfCtrl%DoTranslOff) &
       CALL TranslsOff(CartGrad%D(1:NCart-9),Print2)
     IF(GOpt%TrfCtrl%DoRotOff) &
       CALL RotationsOff(CartGrad%D,Carts%D,Print2,PBCDim)
     CALL Delete(Carts)
     !
     CALL GetCGradMax(CartGrad%D,NCart,GOpt%GOptStat%IMaxCGrad,&
                      GOpt%GOptStat%MaxCGrad,GOpt%GOptStat%ILMaxCGrad, &
                      GOpt%GOptStat%LMaxCGrad)
     !
     CALL GradConv(GOpt%GOptStat,GOpt%GConvCrit,GOpt%Constr)
     !
     ! Do we have to refresh internal coord defs?   
     !
     CALL IntCReDef(GOpt,Refresh,iGEO)
     !
     ! Get internal coord defs.
     !
     IF(Refresh/=0) THEN
       CALL GetIntCs(XYZ,AtNum,IntCs,Refresh,SCRPath, &
                     GOpt%CoordCtrl,GOpt%Constr,GOpt%GConvCrit, &
                     GOpt%GDIIS%MaxMem,GOpt%ExtIntCs,TOPS,Bond, &
                     AtmB,PBCDim,iGEO, &
                     HFileIn_O=HFileIn,iCLONE_O=iCLONE)
       IF(IntCs%N==0) CALL Halt('Molecule has dissociated,'// &
                    'optimizer has not found any internal coordinates.')
     ENDIF
     IF(IntCs%N/=0) THEN
       CALL INTCValue(IntCs,XYZ,PBCDim, &
                      GOpt%CoordCtrl%LinCrit,GOpt%CoordCtrl%TorsLinCrit)
       CALL New(IntOld,IntCs%N)
       IntOld%D=IntCs%Value%D
       CALL INTCValue(IntCL,XYZ,PBCDim,GOpt%CoordCtrl%LinCrit, &
                      GOpt%CoordCtrl%TorsLinCrit)
       CALL New(LatOld,IntCL%N)
       IF(IntCL%N>0) LatOld%D=IntCL%Value%D
     ENDIF
     !
     ! Get B matrices for redundancy projector, etc.
     !
   ! CALL GetMixedBMat(IntCs,XYZ,PBE,GOpt%TrfCtrl,GOpt%CoordCtrl,TOPS, &
   !                   GOpt%Constr%NCartConstr,Print,SCRPath,.TRUE.)
     CALL RefreshBMatInfo(IntCs,XYZ,GOpt%TrfCtrl,GOpt%GConvCrit, &
                    GOpt%CoordCtrl%LinCrit,GOpt%CoordCtrl%TorsLinCrit,&
                    PBCDim,Print,SCRPath)
   ! CALL BuildUMatr(SCRPath,NCart)
     !
     ! Print current set of internals for debugging
     !
     IF(Print==DEBUG_GEOP_MAX) THEN
       CALL PrtIntCoords(IntCs,IntCs%Value%D,&
           'Internals at step #'//TRIM(IntToChar(iGEO)),PBCDim_O=PBCDim)
       IF(PBCDim>0) CALL PrtPBCs(XYZ)
     ENDIF
     !
     ! Calculate simple relaxation step 
     !
     IF(.NOT.GOpt%GOptStat%GeOpConvgd) THEN
       CALL RelaxGeom(GOpt,XYZ,RefXYZ,AtNum,CartGrad%D,iCLONE,ETot, &
               IntCs,iGEO,SCRPath,PWDPath,PBCDim,Print,HFileIn, &
               Refresh,PBCFit)
     ELSE
       WRITE(*,200) iCLONE,iGEO
       WRITE(Out,200) iCLONE,iGEO
       200 FORMAT('     Geometry optimization of Clone #',I2,' converged in ',I3,' steps.')
       Convgd(iCLONE)=1
     ENDIF
     CALL GeOpReview(GOpt%Constr,GOpt%GOptStat,GOpt%CoordCtrl, &
                     GOpt%GConvCrit,XYZ,ETot,IntCs,IntCL,IntOld,LatOld, &
                     iCLONE,iGEO,DoNEB,GOpt%LattIntC,PBCDim)
     CALL DoRestartBas(GOpt%RestartBas,GOpt%GOptStat)
     CALL TurnOnGDIIS(GOpt%GOptStat%MaxCGrad,GOpt%GDIIS%On)
     !
     ! tidy up
     !
     CALL Delete(CartGrad)
     CALL Delete(LatOld)
     CALL Delete(IntCL)
     CALL Delete(Gopt%LattIntC%Grad)
     CALL Delete(Gopt%LattIntC%Displ)
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
     INTEGER                    :: J,N
     !
     CtrlStat%GeOpConvgd=(CtrlStat%MaxCGrad<GConvCr%Grad).AND. &
                         (CtrlStat%LMaxCGrad<GConvCr%Grad)
   END SUBROUTINE GradConv
!
!--------------------------------------------------------------------
!
   SUBROUTINE RelaxGeom(GOpt,XYZ,RefXYZ,AtNum,CartGrad,iCLONE, &
                    ETot,IntCs,IGEO,SCRPath,PWDPath, &
                    PBCDim,Print,HFileIn,Refresh,PBCFit)
     TYPE(GeomOpt)               :: GOpt
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ,RefXYZ
     REAL(DOUBLE),DIMENSION(:)   :: AtNum,CartGrad
     REAL(DOUBLE)                :: ETot     
     TYPE(INTC)                  :: IntCs
     TYPE(PBCFits)               :: PBCFit
     INTEGER                     :: iGEO,iCLONE,Refresh,Print,PBCDim
     CHARACTER(LEN=*)            :: SCRPath,HFileIn,PWDPath 
     TYPE(DBL_RNK2)              :: SRStruct,RefGrad,RefStruct,SRDispl
     !
     ! Relax geometry
     !
     SELECT CASE(GOpt%Optimizer)
     CASE(GRAD_DIAGHESS_OPT)
       CALL RelaxDiagHess(GOpt,SCRPath,PWDPath,CartGrad,IntCs,AtNum,&
                          iGEO,XYZ,Print,PBCDim)
     CASE(GRAD_BiSect_OPT) 
       IF(iGEO<2) THEN
         CALL RelaxDiagHess(GOpt,SCRPath,PWDPath,CartGrad,IntCs,AtNum,&
                            iGEO,XYZ,Print,PBCDim)
       ELSE
         CALL RelaxBiSect(GOpt,SCRPath,PWDPath,HFileIn,CartGrad,IntCs, &
                          AtNum,PBCDim,PBCFit,iGEO,iCLONE,XYZ,RefXYZ, &
                          Print)
       ENDIF
     END SELECT
   END SUBROUTINE RelaxGeom
!
!-------------------------------------------------------
!
   SUBROUTINE RelaxBiSect(GOpt,SCRPath,PWDPath,HFileIn,CartGrad,IntCs, &
                          AtNum,PBCDim,PBCFit,iGEO,iCLONE,XYZ,RefXYZ, &
                          Print)
     TYPE(GeomOpt)               :: GOpt
     TYPE(INTC)                  :: IntCs
     REAL(DOUBLE),DIMENSION(:)   :: AtNum,CartGrad
     INTEGER                     :: iGEO,iCLONE,Print,PBCDim
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ,RefXYZ
     TYPE(PBCFits)               :: PBCFit
     CHARACTER(LEN=*)            :: SCRPath,PWDPath,HFileIn
     TYPE(DBL_VECT)              :: Displ,RefPoints,PredVals
     TYPE(DBL_RNK2)              :: SRStruct,RefStruct,RefGrad,SRDispl
     TYPE(DBL_RNK2)              :: IntCValues,IntCGrads,MixMat
     TYPE(INT_VECT)              :: ISpB,JSpB
     TYPE(DBL_VECT)              :: ASpB
     TYPE(DBL_RNK2)              :: FullB
     INTEGER                     :: I,J,NatmsLoc,NDim,NDimAux,NMem
     INTEGER                     :: HDFFileID,NCart,NMix
     LOGICAL                     :: DoMix,Print2
     !
     NDimAux=MIN(GOpt%GDIIS%MaxMem,iGEO-GOpt%GDIIS%iGEOStart+1)
     NDim=MIN(NDimAux,iGEO)
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     DoMix=.FALSE.
     Print2=(Print>=DEBUG_GEOP_MAX)
     !
     ! calculate internal coord gradients
     !
     CALL CollectPast(RefXYZ,SRStruct,RefStruct,RefGrad,SRDispl, &
                      HFileIn,NatmsLoc,NDim,iGEO,iCLONE)
     !
     ! Substitute last set of cart grads, as it may contain
     ! information (projection) that is not yet put onto disk
     !
     DO J=1,NCart
       RefGrad%D(J,NDim)=CartGrad(J)
     ENDDO
     !
     CALL CollectINTCPast(RefStruct%D,RefGrad%D,IntCValues,IntCGrads, &
                          IntCs,GOpt,SCRPath,Print,PBCDim)
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
     !
   ! CALL CleanConstrIntc(Displ%D,XYZ,GOpt%ExtIntCs,SCRPath,&
   !                      GOpt%TrfCtrl,GOpt%CoordCtrl,GOpt%GConvCrit, &
   !                      PBCDim,Print)
   ! CALL RedundancyOff(Displ%D,SCRPath,Print)  
   ! CALL POffHardGc(IntCs,XYZ,PBCDim,Displ%D,SCRPath,Print2)
     PredVals%D=IntCs%Value%D+Displ%D
     !
     IF(DoMix) THEN
       CALL InternalToCart(XYZ,AtNum,IntCs,PredVals%D,RefPoints%D,Print, &
                           GOpt%BackTrf,GOpt%TrfCtrl,GOpt%CoordCtrl, &
                           GOpt%GConvCrit,GOpt%Constr,PBCDim,SCRPath, &
                           PWDPath,GOpt%ExtIntCs, &
                           Mixmat_O=MixMat%D,iGEO_O=iGEO)
     ELSE
       CALL InternalToCart(XYZ,AtNum,IntCs,PredVals%D,RefPoints%D,Print, &
                           GOpt%BackTrf,GOpt%TrfCtrl,GOpt%CoordCtrl,&
                           GOpt%GConvCrit,GOpt%Constr,PBCDim,SCRPath, &
                           PWDPath,GOpt%ExtIntCs,iGEO_O=iGEO)
     ENDIF
     ! 
     IF(DoMix) THEN
       CALL Delete(MixMat)
     ENDIF
     ! 
     ! Now, add fitting along correlation coordinates
     ! 
 !   IF(NDim>2) THEN
 !     CALL ModeFit(IntCs,IntCGrads%D,IntCValues%D,RefPoints%D,PBCDim, &
 !                  GOpt,PredVals%D,PWDPath,SCRPath,XYZ,NCart,iGEO, &
 !                  Print)
 !     CALL InternalToCart(XYZ,AtNum,IntCs,PredVals%D,RefPoints%D,Print, &
 !                         GOpt%BackTrf,GOpt%TrfCtrl,GOpt%CoordCtrl,&
 !                         GOpt%Constr,PBCDim,SCRPath,PWDPath, &
 !                         GOpt%ExtIntCs,iGEO_O=iGEO)
 !   ENDIF
     ! 
     CALL Delete(SRStruct)
     CALL Delete(RefStruct)
     CALL Delete(RefGrad)
     CALL Delete(SRDispl)
     ! 
     CALL Delete(IntCGrads)
     CALL Delete(IntCValues)
     CALL Delete(PredVals)
     CALL Delete(RefPoints)
     CALL Delete(Displ)
   END SUBROUTINE RelaxBiSect
!
!-------------------------------------------------------
!
   SUBROUTINE RelaxDiagHess(GOpt,SCRPath,PWDPath,CartGrad,IntCs,AtNum,&
                            iGEO,XYZ,Print,PBCDim)
     TYPE(GeomOpt)               :: GOpt
     TYPE(INTC)                  :: IntCs
     REAL(DOUBLE),DIMENSION(:)   :: AtNum,CartGrad
     INTEGER                     :: iGEO,Print,I,J,PBCDim
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     CHARACTER(LEN=*)            :: SCRPath,PWDPath
     TYPE(DBL_VECT)              :: Grad,Displ,RefPoints
     LOGICAL                     :: Print2
     !
     Print2=(Print>=DEBUG_GEOP_MAX)
     CALL New(Grad,IntCs%N)
     CALL New(Displ,IntCs%N)
     CALL New(RefPoints,IntCs%N)
     !
     CALL INTCValue(IntCs,XYZ,PBCDim, &
                    GOpt%CoordCtrl%LinCrit,GOpt%CoordCtrl%TorsLinCrit)
     RefPoints%D=IntCs%Value%D
     !
     ! calculate internal coord gradients
     !
     CALL CartToInternal(IntCs,CartGrad,Grad%D,XYZ,PBCDim, &
       GOpt%GrdTrf,GOpt%CoordCtrl,GOpt%TrfCtrl,Print,SCRPath)
!CALL ProjectBCol(SCRPath,IntCs,XYZ,Grad%D,PBCDim,.TRUE.)
     IF(PBCDim>0) THEN
       CALL PrtIntCoords(IntCs,Grad%D,&
         'Internal Coordinate forces',PBCDim_O=PBCDim)
     ENDIF
     CALL RedundancyOff(Grad%D,SCRPath,Print,Messg_O='IntC Grads')
   ! CALL POffHardGc(IntCs,XYZ,PBCDim,Grad%D,SCRPath,Print2)
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
    !CALL CutOffDispl(Displ%D,IntCs, &
    !                 GOpt%CoordCtrl%MaxStre,GOpt%CoordCtrl%MaxAngle)
     CALL CutOffDispl(Displ%D,IntCs,0.10D0,0.10D0)
   ! CALL RedundancyOff(Displ%D,SCRPath,Print)  
   ! CALL POffHardGc(IntCs,XYZ,PBCDim,Displ%D,SCRPath,Print2)
     ! 
   ! CALL CleanConstrIntc(Displ%D,XYZ,GOpt%ExtIntCs,SCRPath,&
   !                      GOpt%TrfCtrl,GOpt%CoordCtrl,GOpt%GConvCrit, &
   !                      PBCDim,Print)
!CALL ProjectBCol(SCRPath,IntCs,XYZ,Displ%D,PBCDim,.TRUE.)
     IntCs%PredVal%D=IntCs%Value%D+Displ%D
     CALL InternalToCart(XYZ,AtNum,IntCs,IntCs%PredVal%D, &
                         RefPoints%D,Print,GOpt%BackTrf,GOpt%TrfCtrl, &
                         GOpt%CoordCtrl,GOpt%GConvCrit,GOpt%Constr, &
                         PBCDim,SCRPath, &
                         PWDPath,GOpt%ExtIntCs,iGEO_O=iGEO)  
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
!-------------------------------------------------------
!
   SUBROUTINE GeOpReview(CtrlConstr,CtrlStat,CtrlCoord,GConvCr, &
                       XYZ,Etot,IntCs,IntCL,IntOld,LatOld,iCLONE,iGEO,DoNEB, &
                       LattIntC,PBCDim)
     TYPE(GOptStat)             :: CtrlStat
     TYPE(Constr)               :: CtrlConstr
     TYPE(CoordCtrl)            :: CtrlCoord  
     TYPE(GConvCrit)            :: GConvCr    
     TYPE(LattInfo)             :: LattIntC
     REAL(DOUBLE),DIMENSION(:,:):: XYZ
     TYPE(INTC)                 :: IntCs,IntCL
     TYPE(DBL_VECT)             :: IntOld,LatOld,AuxVect
     INTEGER                    :: iCLONE,iGEO,PBCDim
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
     REAL(DOUBLE)               :: RMSGrad,MaxGrad,MaxCGrad,LMaxCGrad
     REAL(DOUBLE)               :: RMSGradNoConstr,MaxGradNoConstr
     REAL(DOUBLE)               :: Etot,Sum,Fact
     !                      
     INTEGER                    :: NStreGeOp,NBendGeOp,NLinBGeOp
     INTEGER                    :: NOutPGeOp,NTorsGeOp
     INTEGER                    :: MaxStre,MaxBend,MaxLinB
     INTEGER                    :: MaxOutP,MaxTors
     !                      
     INTEGER                    :: IMaxGrad,IMaxCGrad,IMaxGradNoConstr
     INTEGER                    :: ILMaxCGrad
     LOGICAL                    :: GradConv
     CHARACTER(LEN=4)           :: ALatt
     !
     NIntC=IntCs%N      
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     !
     RMSGrad=CtrlStat%RMSGrad
     RMSGradNoConstr=CtrlStat%RMSGradNoConstr
     MaxGrad=CtrlStat%MaxGrad
     MaxCGrad=CtrlStat%MaxCGrad
     LMaxCGrad=CtrlStat%LMaxCGrad
     MaxGradNoConstr=CtrlStat%MaxGradNoConstr
     IMaxGrad=CtrlStat%IMaxGrad
     IMaxCGrad=CtrlStat%IMaxCGrad
     ILMaxCGrad=CtrlStat%ILMaxCGrad
     IMaxGradNoConstr=CtrlStat%IMaxGradNoConstr
     !
     IF(PBCDim>0) THEN
       CALL INTCValue(IntCL,XYZ,PBCDim,CtrlCoord%LinCrit,CtrlCoord%TorsLinCrit)
       LattIntC%Displ%D=Zero
       DO J=1,IntCL%N
         LattIntC%Displ%D(J)=IntCL%Value%D(J)-LatOld%D(J)
       ENDDO
     ENDIF
     IF(CtrlStat%GeOpConvgd) THEN
       WRITE(*,399) iCLONE,iGEO,ETot
       WRITE(Out,399) iCLONE,iGEO,ETot
       WRITE(*,140) MaxCGrad,(IMaxCGrad-1)/3+1
       WRITE(Out,140) MaxCGrad,(IMaxCGrad-1)/3+1
       IF(PBCDim>0) THEN
         IF(ILMaxCGrad==NatmsLoc-2) THEN
           ALatt='   A' 
         ELSE IF(ILMaxCGrad==NatmsLoc-1) THEN
           ALatt='   B' 
         ELSE
           ALatt='   C' 
         ENDIF
         WRITE(*,145) LMaxCGrad,ALAtt
         WRITE(Out,145) LMaxCGrad,ALAtt
         CALL LattReview(IntCL,LatOld,LattIntC,PBCDim)
       ENDIF
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
     CALL INTCValue(IntCs,XYZ,PBCDim,CtrlCoord%LinCrit,CtrlCoord%TorsLinCrit)
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
     IF(ILMaxCGrad==NatmsLoc-2) THEN
       ALatt='   A' 
     ELSE IF(ILMaxCGrad==NatmsLoc-1) THEN
       ALatt='   B' 
     ELSE
       ALatt='   C' 
     ENDIF
     IF(PBCDim>0) THEN
       WRITE(*,145) LMaxCGrad,ALAtt
       WRITE(Out,145) LMaxCGrad,ALAtt
     ENDIF
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
     IF(PBCDim>0) CALL LattReview(IntCL,LatOld,LattIntC,PBCDim)
     !
399 FORMAT('       Clone = ',I6,' GeOp step = ',I6,' Total Energy = ',F20.8)
400 FORMAT('Total Energy at Current Geometry = ',F20.8)
401 FORMAT('                    Total Energy = ',F20.8)
410 FORMAT('                   Max intl Grad = ',F12.6,' between atoms ',4I4)
140 FORMAT('     Max Unconstrained Cart Grad = ',F12.6,'      on atom  ',4I4)
145 FORMAT('     Max Unconstrained Latt Grad = ',F12.6,'      on vect  ',A4)
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
!---------------------------------------------------------------------
!
   SUBROUTINE LattReview(IntCL,LatOld,LattIntC,PBCDim)
     TYPE(INTC)     :: IntCL
     TYPE(LattInfo) :: LattIntC
     TYPE(DBL_VECT) :: LatOld
     INTEGER        :: I,J,PBCDim,NStre,NBend
     REAL(DOUBLE)   :: Fact
     !
     IF(PBCDim>0) THEN
       WRITE(*,1010) 
       WRITE(Out,1010) 
       IF(PBCDim==1) THEN
         NStre=1
         NBend=0
       ENDIF
       IF(PBCDim==2) THEN
         NStre=2
         NBend=1
       ENDIF
       IF(PBCDim==3) THEN
         NStre=3
         NBend=3
       ENDIF
       Fact=One/AngstromsToAu
       DO I=1,NStre
         WRITE(*,1020) IntCL%Def%C(I)(1:6),LatOld%D(I)*Fact, &
                       LattIntC%Grad%D(I),LattIntC%Displ%D(I)*Fact, &
                       IntCL%Value%D(I)*Fact
         WRITE(Out,1020) IntCL%Def%C(I)(1:6),LatOld%D(I)*Fact, &
                         LattIntC%Grad%D(I),LattIntC%Displ%D(I)*Fact, &
                         IntCL%Value%D(I)*Fact
       ENDDO
       Fact=180.D0/PI
       DO I=NStre+1,NStre+NBend
         WRITE(*,1020) IntCL%Def%C(I)(1:6),LatOld%D(I)*Fact, &
                       LattIntC%Grad%D(I),LattIntC%Displ%D(I)*Fact, &
                       IntCL%Value%D(I)*Fact
         WRITE(Out,1020) IntCL%Def%C(I)(1:6),LatOld%D(I)*Fact, &
                         LattIntC%Grad%D(I),LattIntC%Displ%D(I)*Fact, &
                         IntCL%Value%D(I)*Fact
       ENDDO
       !
       IF(PBCDim==3.AND.IntCL%Def%C(7)(1:8)=='VOLM_L  ') THEN
         Fact=One/(AngstromsToAu**3)
         I=7
         WRITE(*,1020) IntCL%Def%C(I)(1:6),LatOld%D(I)*Fact, &
                       LattIntC%Grad%D(I),LattIntC%Displ%D(I)*Fact, &
                       IntCL%Value%D(I)*Fact
         WRITE(Out,1020) IntCL%Def%C(I)(1:6),LatOld%D(I)*Fact, &
                         LattIntC%Grad%D(I),LattIntC%Displ%D(I)*Fact, &
                         IntCL%Value%D(I)*Fact
       ENDIF
       !
       IF(PBCDim==2.AND.IntCL%Def%C(4)(1:8)=='AREA_L  ') THEN
         Fact=One/(AngstromsToAu**2)
         I=4
         WRITE(*,1020) IntCL%Def%C(I)(1:6),LatOld%D(I)*Fact, &
                       LattIntC%Grad%D(I),LattIntC%Displ%D(I)*Fact, &
                       IntCL%Value%D(I)*Fact
         WRITE(Out,1020) IntCL%Def%C(I)(1:6),LatOld%D(I)*Fact, &
                         LattIntC%Grad%D(I),LattIntC%Displ%D(I)*Fact, &
                         IntCL%Value%D(I)*Fact
       ENDIF
     ENDIF
     !
     1010 FORMAT('Lattice Parameter        Old Value      Gradient    Displacement   New Value')
     1020 FORMAT(10X,A6,4X,4F14.6)
     !
   END SUBROUTINE LattReview
!
!---------------------------------------------------------------
!
   SUBROUTINE SetGeOpCtrl(GOpt,Geos,Opts,Sets,Nams,PBCs,iGEO)
     !
     TYPE(GeomOpt)    :: GOpt
     TYPE(Options)    :: Opts
     TYPE(BasisSets)  :: Sets
     TYPE(Geometries) :: Geos
     TYPE(FileNames)  :: Nams
     TYPE(Periodics)  :: PBCs
     INTEGER          :: NatmsLoc,NCart
     REAL(DOUBLE)     :: Sum,GCrit
     INTEGER          :: AccL,iGEO
     !
     AccL    =Opts%AccuracyLevels(Sets%NBSets)
     NatmsLoc=Geos%Clone(1)%Natms
     NCart=3*NatmsLoc
     !
     GOpt%RestartBas=.FALSE.
     CALL SetCoordCtrl(GOpt%CoordCtrl)
     CALL   SetHessian(GOpt%Hessian)
     CALL SetGConvCrit(GOpt%GConvCrit,GOpt%Hessian,AccL,NatmsLoc)
     CALL     SetGDIIS(GOpt%GDIIS,GOpt%Optimizer,iGEO)
     CALL    SetGrdTrf(GOpt%GrdTrf,GOpt%GConvCrit)
     CALL   SetBackTrf(GOpt%BackTrf,GOpt%GConvCrit)
     CALL    SetConstr(GOpt%Constr,GOpt%BackTrf)
     CALL   SetTrfCtrl(GOpt%TrfCtrl,GOpt%CoordCtrl,GOpt%Constr,PBCs%Dimen)
   END SUBROUTINE SetGeOpCtrl
!
!---------------------------------------------------------------
!
   SUBROUTINE OptSingleMol(GOpt,Nams,Opts, &
                           GMLoc,Convgd,iGEO,iCLONE)
     TYPE(GeomOpt)        :: GOpt
     TYPE(FileNames)      :: Nams
     TYPE(Options)        :: Opts
     TYPE(CRDS)           :: GMLoc
     INTEGER,DIMENSION(:) :: Convgd
     INTEGER              :: I,J,K,iGEO,iCLONE,NatmsNew
     INTEGER              :: InitGDIIS,NConstr,NCart,NatmsLoc
     LOGICAL              :: NoGDIIS,GDIISOn,DoNEB
     CHARACTER(LEN=DCL)   :: SCRPath,PWDPath
     TYPE(DBL_RNK2)       :: XYZNew,GradNew,RefXYZ1,RefXYZ
     TYPE(DBL_VECT)       :: AtNumNew,CartGrad
     !
     SCRPath  =TRIM(Nams%M_SCRATCH)//TRIM(Nams%SCF_NAME)// &
             '.'//TRIM(IntToChar(iCLONE))
     PWDPath  =TRIM(Nams%M_PWD)//TRIM(IntToChar(iCLONE))
     GMLoc%Displ%D=GMLoc%AbCarts%D
     DoNEB=(Opts%Grad==GRAD_TS_SEARCH_NEB)
     !
     CALL New(RefXYZ1,(/3,GMLoc%Natms/))
     CALL GetRefXYZ(Nams%HFile,RefXYZ1,iCLONE)
     !
     ! Now, cut out that part of the molecule, which is 
     ! not rigidly fixed and pass it in to the optimizer
     !
     CALL OpenASCII(OutFile,Out)
     NatmsNew=0
     DO I=1,GMLoc%Natms
       IF(GMLoc%CConstrain%I(I)/=2) NatmsNew=NatmsNew+1
     ENDDO
     NatmsNew=NatmsNew+3 ! Added place for PBC data
     CALL New(XYZNew,(/3,NatmsNew/))
     CALL New(AtNumNew,NatmsNew)
     CALL New(GradNew,(/3,NatmsNew/))
     CALL New(RefXYZ,(/3,NatmsNew/))
     ! fill atomic data
     NatmsNew=0
     DO I=1,GMLoc%Natms
       IF(GMLoc%CConstrain%I(I)/=2) THEN
         NatmsNew=NatmsNew+1
         XYZNew%D(1:3,NatmsNew)=GMLoc%Displ%D(1:3,I)
         RefXYZ%D(1:3,NatmsNew)=RefXYZ1%D(1:3,I)
         GradNew%D(1:3,NatmsNew)=GMLoc%Gradients%D(1:3,I)
         AtNumNew%D(NatmsNew)=GMLoc%AtNum%D(NatmsNew)
       ENDIF
     ENDDO
     DO K=1,3
       DO J=1,3
         XYZNew%D(J,NatmsNew+K)=GMLoc%PBC%BoxShape%D(J,K)
         RefXYZ%D(J,NatmsNew+K)=GMLoc%PBC%BoxShape%D(J,K)
         GradNew%D(J,NatmsNew+K)=GMLoc%PBC%LatFrc%D(J,K)
       ENDDO
     ENDDO
     ! ensure proper orientation of (numerical) forces
    !GradNew%D(2:3,NatmsNew+1)=Zero
    !GradNew%D(3,NatmsNew+2)=Zero
     CALL Delete(RefXYZ1)
     AtNumNew%D(NatmsNew+1:NatmsNew+3)=Zero
     CALL ConvertToXYZRef(XYZNew%D,RefXYZ%D,GMLoc%PBC%Dimen)
     !
     !--------------------------------------------
     !
       CALL ModifyGeom(GOpt,XYZNew%D,RefXYZ%D,AtNumNew%D,GMLoc%IntCs, &
                       GradNew%D,GMLoc%Bond,GMLoc%AtmB,Convgd, &
                       GMLoc%Etotal,GMLoc%PBC%Dimen,iGEO,iCLONE, &
                       SCRPath,PWDPath,DoNEB,Opts%PFlags%GeOp, &
                       Nams%HFile,GMLoc%PBCFit)
     CLOSE(Out,STATUS='KEEP')
     !
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
     DO K=1,3
       DO J=1,3
         GMLoc%PBCDispl%BoxShape%D(J,K)=XYZNew%D(J,NatmsNew+K)
         GMLoc%PBC%LatFrc%D(J,K)=GradNew%D(J,NatmsNew+K)
       ENDDO
     ENDDO
     !
     CALL Delete(XYZNew)
     CALL Delete(AtNumNew)
     CALL Delete(GradNew)
     CALL Delete(RefXYZ)
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
   SUBROUTINE SetGDIIS(GD,GOptimizer,iGEO)
     TYPE(GDIIS)  :: GD
     INTEGER      :: GOptimizer,iGEO
     !
     GD%Init    = 3
     IF(GOptimizer==GRAD_BiSect_OPT) THEN
       GD%MaxMem  = 7 
     ELSE
       GD%MaxMem  = 7
     ENDIF
     GD%On=.TRUE.
     GD%iGEOStart=MAX(1,iGEO-GD%MaxMem+1)
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
   ! BackT%CooTrfCrit   = MIN(GConv%Stre/10.D0,1.D-6)
     BackT%CooTrfCrit   = 1.D-5
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
   SUBROUTINE SetTrfCtrl(TrfC,CoordC,GConstr,Dimen)
     TYPE(TrfCtrl)   :: TrfC
     TYPE(CoordCtrl) :: CoordC
     TYPE(Constr)    :: GConstr
     INTEGER         :: Dimen
     !
     TrfC%DoTranslOff=.TRUE.
     TrfC%DoRotOff=.TRUE.
     IF(GConstr%NCartConstr>0) TrfC%DoTranslOff=.FALSE.
     IF(GConstr%NCartConstr>0.AND.Dimen==0) TrfC%DoRotOff=.FALSE.
     !
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
     CoordC%LinCrit = 5.D0 
     CoordC%TorsLinCrit = 5.D0      
     CoordC%OutPCrit=CoordC%LinCrit
   END SUBROUTINE SetCoordCtrl
!
!-------------------------------------------------------------------
!
   SUBROUTINE NewGeomFill(GMLoc)
     TYPE(CRDS)           :: GMLoc
     !
     GMLoc%AbCarts%D=GMLoc%Displ%D
     GMLoc%PBC%BoxShape%D=GMLoc%PBCDispl%BoxShape%D
     CALL PBCInfoFromNewCarts(GMLoc%PBC)
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
       G%Clone(iCLONE)%IntCs%N=0
       CALL New(G%Clone(iCLONE)%IntCs,G%Clone(iCLONE)%IntCs%N)
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
     CHARACTER(LEN=DCL)   :: SCRPath,PWDPath
     INTEGER,DIMENSION(:) :: Convgd
     TYPE(DBL_RNK2)       :: RefXYZ1
     !
     SCRPath  =TRIM(Nams%M_SCRATCH)//TRIM(Nams%SCF_NAME)// &
             '.'//TRIM(IntToChar(iCLONE))
     PWDPath=TRIM(Nams%M_PWD)//TRIM(IntToChar(iCLONE))
     !
     CALL New(RefXYZ1,(/3,GMLoc%Natms/))
     CALL GetRefXYZ(Nams%HFile,RefXYZ1,iCLONE)
     !
     IF(Convgd(iCLONE)/=1) THEN
       CALL OPENAscii(OutFile,Out)
       IF((.NOT.GOpt%GDIIS%NoGDIIS).AND.GOpt%GDIIS%On.AND.&
           iGEO>=GOpt%GDIIS%Init) THEN
         CALL GeoDIIS(GMLoc%AbCarts%D,RefXYZ1%D,GOpt%GDIIS,Nams%HFile,iCLONE, &
                      iGEO,Opts%PFlags%GeOp,PWDPath)
       ELSE
         IF(Opts%PFlags%GeOp>=DEBUG_GEOP_MIN) THEN
           WRITE(*,200)
           WRITE(Out,200)
         ENDIF
         200 FORMAT('No Geometric DIIS is being done in this step.')
       ENDIF
       CLOSE(Out,STATUS='KEEP')
     ENDIF
     CALL Delete(RefXYZ1)
   END SUBROUTINE MixGeoms
!
!---------------------------------------------------------------------
!
   SUBROUTINE GetCGradMax(CartGrad,NCart,IMaxCGrad,MaxCGrad, &
                          ILMaxCGrad,LMaxCGrad)
     REAL(DOUBLE),DIMENSION(:) :: CartGrad
     REAL(DOUBLE),DIMENSION(3) :: Vect    
     REAL(DOUBLE)              :: MaxCGrad,Sum,LMaxCGrad
     INTEGER                   :: Nat,NCart,I,I1,I2,IMaxCGrad,ILMaxCGrad
     !
     IMaxCGrad=1   
     MaxCGrad=Zero
     Nat=NCart/3
     DO I=1,Nat-3
       I1=3*(I-1)+1
       I2=I1+2
       Vect=CartGrad(I1:I2)
       Sum=SQRT(DOT_PRODUCT(Vect,Vect))
       IF(Sum>MaxCGrad) THEN
         IMaxCGrad=I
          MaxCGrad=Sum
       ENDIF
     ENDDO
     !
     LMaxCGrad=Zero
     DO I=Nat-2,Nat
       I1=3*(I-1)+1
       I2=I1+2
       Vect=CartGrad(I1:I2)
       Sum=SQRT(DOT_PRODUCT(Vect,Vect))
       IF(Sum>LMaxCGrad) THEN
         ILMaxCGrad=I
          LMaxCGrad=Sum
       ENDIF
     ENDDO
   END SUBROUTINE GetCGradMax
!
!-------------------------------------------------------------------
!
   SUBROUTINE BackTrack(iBAS,iGEO,C,BPrev,BCur,DoLineS_O)
     ! Go over clones and do backtracking whenever necessary
     TYPE(Controls)   :: C
     INTEGER          :: I,J,iBAS,iGEO,iCLONE,I1,I2
     INTEGER,DIMENSION(:):: BPrev,BCur
     INTEGER          :: NatmsLoc,NCart,AccL
     INTEGER          :: MaxBStep,IBStep
     CHARACTER(LEN=DCL):: chGEO
     TYPE(CRDS)       :: GMOld
     LOGICAL          :: DoBackTrack,DoLineS
     LOGICAL,OPTIONAL :: DoLineS_O
     REAL(DOUBLE)     :: EOld,ENew,MeanDist,FactN,FactO,Fact,MaxDist
     TYPE(DBL_VECT)   :: DistVect1,DistVect2
     TYPE(DBL_RNK2)   :: RefXYZ1
     TYPE(LOG_VECT)   :: NeedBackTr
     !
     AccL=C%Opts%AccuracyLevels(iBAS)
     IF(.NOT.C%GOpt%GConvCrit%DoLattStep) THEN
       C%GOpt%GConvCrit%DoLattStep=.TRUE.
       RETURN
     ENDIF
     IF(.NOT.C%GOpt%GConvCrit%DoAtomBackTr.AND. &
        .NOT.C%GOpt%GConvCrit%DoLattBackTr) THEN
       RETURN
     ENDIF
     chGEO=IntToChar(iGEO-1)
     CALL New(NeedBackTr,C%Geos%Clones)
     HDFFileID=OpenHDF(C%Nams%HFile)
     DO iCLONE=1,C%Geos%Clones
       HDF_CurrentID=OpenHDFGroup(HDFFileID, &
                   "Clone #"//TRIM(IntToChar(iCLONE)))
       CALL Get(EOld,'gm_etot'      ,Tag_O=ChGEO)
       ENew=C%Geos%Clone(iCLONE)%ETotal
       NeedBackTr%L(iCLONE)=((ENew-EOld)/ABS(EOld)>10.D0*ETol(AccL))
       CALL CloseHDFGroup(HDF_CurrentID)
     ENDDO
     CALL CloseHDF(HDFFileID)
     DoBackTrack=.FALSE.
     DO iCLONE=1,C%Geos%Clones
       DoBackTrack=DoBackTrack.OR.NeedBackTr%L(iCLONE)
     ENDDO
     IF(.NOT.DoBackTrack) THEN
       CALL Delete(NeedBackTr)
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
     MaxBStep=11
     IF(DoLineS) MaxBStep=1
     !
     DO iBStep=1,MaxBStep+1
       DoBackTrack=.FALSE.
       DO iCLONE=1,C%Geos%Clones
         NatmsLoc=C%Geos%Clone(iCLONE)%Natms
         NCart=3*NatmsLoc
         CALL New(RefXYZ1,(/3,NatmsLoc/))
         CALL GetRefXYZ(C%Nams%HFile,RefXYZ1,iCLONE)
         HDFFileID=OpenHDF(C%Nams%HFile)
         HDF_CurrentID=OpenHDFGroup(HDFFileID, &
                     "Clone #"//TRIM(IntToChar(iCLONE)))
         chGEO=IntToChar(iGEO-1)
         CALL Get(GMOld,chGEO)
         CALL ConvertToXYZRef(GMOld%Carts%D,RefXYZ1%D, &
                              GMOld%PBC%Dimen,&
                              BoxShape_O=GMOld%PBC%BoxShape%D)
         CALL ConvertToXYZRef(C%Geos%Clone(iCLONE)%Carts%D,RefXYZ1%D,&
                        C%Geos%Clone(iCLONE)%PBC%Dimen, &
                        BoxShape_O=C%Geos%Clone(iCLONE)%PBC%BoxShape%D)
         ! Convert BoxCarts to same reference system
       ! C%Geos%Clone(iCLONE)%PBC%InvBoxSh%D= &
       !              InverseMatrix(C%Geos%Clone(iCLONE)%PBC%BoxShape%D)
       ! GMOld%PBC%InvBoxSh%D=InverseMatrix(GMOld%PBC%BoxShape%D)
         DO I=1,NatmsLoc
           CALL DGEMM_NNc(3,3,1,One,Zero, &
                C%Geos%Clone(iCLONE)%PBC%InvBoxSh%D, &
                C%Geos%Clone(iCLONE)%Carts%D(1:3,I), &
                C%Geos%Clone(iCLONE)%BoxCarts%D(1:3,I))   
           CALL DGEMM_NNc(3,3,1,One,Zero, &
                GMOld%PBC%InvBoxSh%D, &
                GMOld%Carts%D(1:3,I), &
                GMOld%BoxCarts%D(1:3,I))   
         ENDDO
         !
         CALL Delete(RefXYZ1)
         CALL CloseHDFGroup(HDF_CurrentID)
         CALL CloseHDF(HDFFileID)
         EOld=GMOld%ETotal
         ENew=C%Geos%Clone(iCLONE)%ETotal
         !
         CALL New(DistVect1,NCart)
         CALL CartRNK2ToCartRNK1(DistVect1%D,GMOld%Carts%D)
         CALL New(DistVect2,NCart)
         CALL CartRNK2ToCartRNK1(DistVect2%D,C%Geos%Clone(iCLONE)%Carts%D)
         DistVect2%D=DistVect2%D-DistVect1%D
         MeanDist=SQRT(DOT_PRODUCT(DistVect2%D,DistVect2%D))/DBLE(NatmsLoc)
         MaxDist=Zero
         DO I=1,NatmsLoc
           I1=3*(I-1)+1
           I2=3*(I-1)+3
           Fact=SQRT(DOT_PRODUCT(DistVect2%D(I1:I2),DistVect2%D(I1:I2)))
           MaxDist=MAX(MaxDist,Fact)
         ENDDO 
         CALL Delete(DistVect1)
         CALL Delete(DistVect2)
         !
         DoBackTrack=((ENew-EOld)/ABS(EOld)>10.D0*ETol(AccL))
         !
         IF(iBStep>1.OR.DoBackTrack) THEN  
             CALL OPENAscii(OutFile,Out)
             WRITE(*,200) iBStep-1,EOld,ENew,MeanDist,MaxDist
             WRITE(Out,200) iBStep-1,EOld,ENew,MeanDist,MaxDist
           200 FORMAT('Backtr. step= ',I3,' Old Energy= ', &
                       F20.6,' New Energy= ',F20.6,' Dist= ',F14.8, &
                       'MaxDist= ',F14.8)
             CLOSE(Out,STATUS='KEEP')
         ENDIF
         !
         IF(DoBackTrack) THEN  
           ! do bisection
           IF(DoLineS) THEN
             CALL LineSearch(C%Geos%Clone(iCLONE),GMOld)
           ELSE
             !
             ! Set up new box-coordinates
             !
             FactO=Half
             FactN=Half
             IF(C%GOpt%GConvCrit%DoLattBackTr) THEN
               C%Geos%Clone(iCLONE)%PBC%BoxShape%D= &
                 (FactN*C%Geos%Clone(iCLONE)%PBC%BoxShape%D+ &
                                 FactO*GMOld%PBC%BoxShape%D)
               CALL PBCInfoFromNewCarts(C%Geos%Clone(iCLONE)%PBC)
             ENDIF
             !
             ! Set up new atomic positions (fractionals+Carts)
             !
             IF(C%GOpt%GConvCrit%DoAtomBackTr) THEN
               !
               C%Geos%Clone(iCLONE)%BoxCarts%D= &
                 (FactN*C%Geos%Clone(iCLONE)%BoxCarts%D+&
                  FactO*GMOld%BoxCarts%D)
               !
               DO I=1,NatmsLoc
                 CALL DGEMM_NNc(3,3,1,One,Zero, &
                      C%Geos%Clone(iCLONE)%PBC%BoxShape%D, &
                      C%Geos%Clone(iCLONE)%BoxCarts%D(1:3,I), &
                      C%Geos%Clone(iCLONE)%Carts%D(1:3,I))   
               ENDDO
               !
               C%Geos%Clone(iCLONE)%AbCarts%D= &
               C%Geos%Clone(iCLONE)%Carts%D
             ENDIF
             !
           ENDIF
         ENDIF
       ENDDO
       !
       IF(DoBackTrack) THEN
         IF(iBStep>MaxBStep) THEN
             CALL OPENAscii(OutFile,Out)
             WRITE(*,100) MaxBStep-1
             WRITE(Out,100) MaxBStep-1
           100 FORMAT('Backtracking has not converged in ', &
                       I3,' Steps. Continue with present geometry.')
             CLOSE(Out,STATUS='KEEP')
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
     ENDDO
     ! In the present version do not do any overwriting of 
     ! the GMOld%Displ, since it refers 
     ! strictly to the simple relaxation step, for GDIIS.
     !
     CALL Delete(GMOld)
     CALL Delete(NeedBackTr)
   END SUBROUTINE BackTrack
!
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
     IF(MaxCGrad<1.000D0) THEN
    !IF(MaxCGrad<0.001D0) THEN
       On=.TRUE.
     ENDIF
   END SUBROUTINE TurnOnGDIIS
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
   SUBROUTINE ProjectConstr(G)
     TYPE(Geometries)   :: G
     INTEGER            :: iCLONE
     !
     DO iCLONE=1,G%Clones 
     ENDDO
   END SUBROUTINE ProjectConstr
!
!-------------------------------------------------------------------
!
  SUBROUTINE NHessian(C)
    !-------------------------------------------------------------------
    TYPE(Controls) :: C
    !-------------------------------------------------------------------
    INTEGER        :: iBAS,iGEO,iCLONE,AtA,AtB,I,J,II,JJ,IMin,N
    INTEGER        :: NatmsLoc,IStart,GBeg,GEnd
    TYPE(DBL_RNK2) :: Grad0,Hess,P,NMode
    TYPE(DBL_VECT) :: Freqs,EigenV,Work
    REAL(DOUBLE)   :: SMA,SMB
    INTEGER        :: LWORK,Info,NNeg
    ! FreqToWaveNbr=E*E*NA/4*PI*PI*C*C*A0*A0*A0,
    !            E =4.803242E-10 SQRT(G*CM**3)/SEC
    !            NA=6.022045E+23 /G
    !            A0=5.291771E-09 A
    !            C =2.997925E+10 A/SEC
    REAL(DOUBLE), PARAMETER :: FreqToWaveNbr=2.642461D+07
    REAL(DOUBLE), PARAMETER :: ATOM_DISP=0.001D0
    REAL(DOUBLE), EXTERNAL  :: DDOT
    !-------------------------------------------------------------------
    !
    GBeg=1
    GEnd=C%Geos%Clones
    ! initial geometry
    iGEO=C%Stat%Previous%I(3)
    NatmsLoc=C%Geos%Clone(1)%Natms
    N=3*NatmsLoc
    CALL New(Grad0,(/3,NatmsLoc/))
    CALL DBL_VECT_EQ_DBL_SCLR(N,Grad0%D(1,1),0.0D0)
    CALL New(Hess,(/N,N/))
    CALL DBL_VECT_EQ_DBL_SCLR(N**2,Hess%D(1,1),0.0D0)
    CALL New(P,(/N,N/))
    CALL New(Freqs,N)
    CALL New(NMode,(/N,N/))
    !
    ! Build the guess 
    DO iBAS=1,C%Sets%NBSets
       CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)    
       CALL BSetArchive(iBAS,C%Nams,C%Opts,C%Geos,C%Sets,C%MPIs)
       CALL SCF(iBAS,iGEO,C)
    ENDDO
    iBAS=C%Sets%NBSets
    CALL Force(iBAS,iGEO,C%Nams,C%Opts,C%Stat,C%Geos,C%Sets,C%MPIs)
    CALL DCOPY(N,C%Geos%Clone(1)%Gradients%D(1,1),1,Grad0%D(1,1),1)
    CALL Print_DBL_RNK2(Grad0,'Grad0',Unit_O=6)
    !
    iGEO=iGEO+1
    J=1
    DO AtA=1,NAtmsLoc
       DO I=1,3
          !
          ! Move the atom.
          DO iCLONE=1,C%Geos%Clones
             C%Geos%Clone(iCLONE)%AbCarts%D(I,AtA)=C%Geos%Clone(iCLONE)%AbCarts%D(I,AtA)+ATOM_DISP
             CALL MakeGMPeriodic(C%Geos%Clone(iCLONE))
          ENDDO
          CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)
          !
          ! Optimize and Forces.
          CALL SCF(iBAS,iGEO,C)
          CALL Force(iBAS,iGEO,C%Nams,C%Opts,C%Stat,C%Geos,C%Sets,C%MPIs)
          !CALL Print_DBL_RNK2(C%Geos%Clone(1)%Carts,'Posi',Unit_O=6)
          CALL Print_DBL_RNK2(C%Geos%Clone(1)%Gradients,'Grad',Unit_O=6)
          CALL DCOPY(N,C%Geos%Clone(1)%Gradients%D(1,1),1,Hess%D(1,J),1)
          CALL DAXPY(N,-1D0,Grad0%D(1,1),1,Hess%D(1,J),1)
          !
          ! Move back.
          DO iCLONE=1,C%Geos%Clones
             C%Geos%Clone(iCLONE)%AbCarts%D(I,AtA)=C%Geos%Clone(iCLONE)%AbCarts%D(I,AtA)-ATOM_DISP
             CALL MakeGMPeriodic(C%Geos%Clone(iCLONE))
          ENDDO
          CALL GeomArchive(iBAS,iGEO,C%Nams,C%Sets,C%Geos)
          iGEO=iGEO+1
          J=J+1
       ENDDO
    ENDDO
    CALL DSCAL(N**2,1D0/ATOM_DISP,Hess%D(1,1),1)
    !
    ! Symmetrize the Hessian
    CALL SymmMtrx(Hess,N)
    !
    ! Generate mass weighted hessian.
    CALL WghtMtrx(C%Geos%Clone(1)%AtMss,Hess,NatmsLoc,'CToWC')
    !CALL Print_DBL_RNK2(Hess,'Weighted Hessian',Unit_O=6)
    !
    ! Project out translations and rotations.
    CALL ProjOut(C%Geos%Clone(1)%AbCarts,C%Geos%Clone(1)%AtMss,Hess,NatmsLoc)
    !CALL Print_DBL_RNK2(Hess,'Projected Hessian',Unit_O=6)
    !
    ! Get Normal modes.
    CALL New(EigenV,N)
    LWORK=MAX(1,3*N+10)
    CALL New(WORK,LWORK)
    CALL DCOPY(N**2,Hess%D(1,1),1,NMode%D(1,1),1)
    CALL DSYEV('V','U',N,NMode%D(1,1),N,EigenV%D(1),Work%D(1),LWORK,Info)
    IF(Info/=SUCCEED)CALL Halt('DSYEV flaked in NHessian. INFO='//TRIM(IntToChar(Info)))
    CALL DCOPY(N,EigenV%D(1),1,Freqs%D(1),1)
    CALL Delete(EigenV)
    CALL Delete(Work)
    !CALL Print_DBL_RNK2(NMode,'NMode0',Unit_O=6)
    !
    ! Translational and rotational Sayvetz conditions.
    CALL Sayvetz()
    !
    ! Count number of negative eigenvalues.
    NNeg=0
    DO I=1,N
       IF(Freqs%D(I).LT.0.0D0) NNeg=NNeg+1
    ENDDO
    !
    ! Convert frequencies to wavenumbers.
    !CALL Print_DBL_VECT(Freqs,'Frequencies**2',Unit_O=6)
    DO I=1,N
       Freqs%D(I)=SQRT(ABS(FreqToWaveNbr*Freqs%D(I)))
    ENDDO
    !
    ! Convert normal mode displacements from mass 
    ! weighted cartesian to cartesian.
    DO AtB=1,NatmsLoc
       DO J=1,3
          JJ=3*(AtB-1)+J
          DO AtA=1,NatmsLoc
             SMA=1D0/SQRT(C%Geos%Clone(1)%AtMss%D(AtA))
             DO I=1,3
                II=3*(AtA-1)+I
                NMode%D(II,JJ)=NMode%D(II,JJ)*SMA
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    ! Compute reduced mass.
    !DO I=1,N
    !   write(*,*) 'Reduced Mass=',1D0/DDOT(N,NMode%D(1,I),1,NMode%D(1,I),1)
    !END DO
    !
    ! Compute IR intensities.
    CALL CompIRInt()
    !
    ! Compute Raman intensities.
    CALL CompRAMANInt()
    !
    ! Statistic.
    CALL ThermoStat(C%Geos%Clone(1)%AbCarts,C%Geos%Clone(1)%AtMss,Freqs, &
         &          NatmsLoc,C%Geos%Clone(1)%Multp)
    !
    DO I=1,N
       IF(ABS(Freqs%D(I)).GT.1D-3)write(*,'(A,F10.4)') 'Frequencies=',Freqs%D(I)
    ENDDO
    !CALL Print_DBL_VECT(Freqs,'Frequencies' ,Unit_O=6)
    !CALL Print_DBL_RNK2(NMode,'Normal Modes',Unit_O=6)
    !
    ! Clean up.
    CALL Delete(Freqs)
    CALL Delete(NMode)
    CALL Delete(Grad0)
    CALL Delete(P)
    CALL Delete(Hess)
    !
  END SUBROUTINE NHessian
  !
  SUBROUTINE WghtMtrx(AtMss,A,NAtoms,FromTo)
    !-------------------------------------------------------------------
    TYPE(DBL_VECT)   :: AtMss
    TYPE(DBL_RNK2)   :: A
    INTEGER          :: NAtoms
    CHARACTER(LEN=*) :: FromTo
    !-------------------------------------------------------------------
    INTEGER          :: AtA,AtB,I,II,J,JJ,IFromTo
    REAL(DOUBLE)     :: SMA,SMB
    !-------------------------------------------------------------------
    SELECT CASE(FromTo)
    CASE('CToWC');IFromTo=0
    CASE('WCToC');IFromTo=1
    CASE DEFAULT; STOP'Err:WghtMtrx'
    END SELECT
    DO AtB=1,NAtoms
       IF(IFromTo.EQ.0) THEN
          SMB=1.0D0/SQRT(AtMss%D(AtB))
       ELSE
          SMB=SQRT(AtMss%D(AtB))
       ENDIF
       DO J=1,3
          JJ=3*(AtB-1)+J
          DO AtA=1,NAtoms
             IF(IFromTo.EQ.0) THEN
                SMA=1.0D0/SQRT(AtMss%D(AtA))
             ELSE
                SMA=SQRT(AtMss%D(AtA))
             ENDIF
             DO I=1,3
                II=3*(AtA-1)+I
                A%D(II,JJ)=A%D(II,JJ)*SMA*SMB
             ENDDO
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE WghtMtrx
  !
  SUBROUTINE ProjOut(Carts,AtMss,A,NAtoms)
    !-------------------------------------------------------------------
    TYPE(DBL_RNK2) :: Carts,A
    TYPE(DBL_VECT) :: AtMss
    INTEGER        :: NAtoms
    !-------------------------------------------------------------------
    TYPE(DBL_RNK2) :: WCarts,Theta,ITheta,P,Tmp
    TYPE(DBL_VECT) :: EigenV,Work
    INTEGER        :: LWork,Info,N,I,J,II,JJ,AtA,AtB,IMin,IA,IB,IC,JA,JB,JC
    INTEGER        :: KK,LL
    REAL(DOUBLE)   :: MssTot,SMA,SMB,Dum,Tens(3,3,3),CX,CY,CZ,CMss(3)
    DATA Tens/ 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,-1.0D0, &
         &     0.0D0, 1.0D0, 0.0D0, 0.0D0, 0.0D0, 1.0D0, &
         &     0.0D0, 0.0D0, 0.0D0,-1.0D0, 0.0D0, 0.0D0, &
         &     0.0D0,-1.0D0, 0.0D0, 1.0D0, 0.0D0, 0.0D0, &
         &     0.0D0, 0.0D0, 0.0D0  /
    !-------------------------------------------------------------------
    !
    N=3*NAtoms
    !
    ! Allocate arrays.
    CALL New(Theta ,(/3,3/))
    CALL New(ITheta,(/3,3/))
    CALL New(P     ,(/N,N/))
    CALL New(WCarts,(/3,NAtoms/))
    !
    ! Center of mass.
    MssTot=0.0D0
    CMss(1)=0.0D0
    CMss(2)=0.0D0
    CMss(3)=0.0D0
    DO AtA=1,NAtoms
       MssTot=MssTot+AtMss%D(AtA)
       DO I=1,3
          CMss(I)=CMss(I)+AtMss%D(AtA)*Carts%D(I,AtA)
       ENDDO
    ENDDO
    DO I=1,3
       CMss(I)=CMss(I)/MssTot
    ENDDO
    !
    ! Weight coordinates.
    DO AtA=1,NAtoms
       SMA=SQRT(AtMss%D(AtA))
       DO I=1,3
          WCarts%D(I,AtA)=SMA*(Carts%D(I,AtA)-CMss(I))
       ENDDO
    ENDDO
    !
    ! Compute I0.
    CALL New(EigenV,3)
    CALL New(Tmp,(/3,3/))
    CALL Inertia(Carts,AtMss,NAtoms,Theta%D,EigenV%D,Tmp%D)
    CALL Delete(Tmp)
    !
    ! Inverse I0.
    IF(ABS(EigenV%D(1)*EigenV%D(2)*EigenV%D(3)).LT.1.0D-8) THEN
       ! Linear systems.
       WRITE(*,*) 'Linear molecules are not supported yet.'
       WRITE(*,*) 'STOPPED in ProjOut'
       STOP 
    ELSE
       ! Non-linear systems.
       ITheta%D(:,:)=InverseMatrix(Theta%D(:,:))
       !CALL Print_DBL_RNK2(ITheta,'Inv(I0)',Unit_O=6)    
    ENDIF
    CALL Delete(EigenV)
    !
    ! Compute the projection matrix.
    CALL DBL_VECT_EQ_DBL_SCLR(N**2,P%D(1,1),0.0D0)
    DO J=1,NAtoms
       SMB=SQRT(AtMss%D(J))
       JJ=3*(J-1)
       DO JC=1,3
          LL=JJ+JC
          DO I=J,NAtoms
             SMA=SQRT(AtMss%D(I))
             II=3*(I-1)
             IMin=1
             IF(I.EQ.J)IMin=JC
             DO IC=IMin,3
                KK=II+IC
                ! Rotation.
                Dum=0.0D0
                DO IA=1,3
                DO IB=1,3
                   IF(Tens(IA,IB,IC).EQ.ZERO) CYCLE
                   DO JA=1,3
                   DO JB=1,3
                      IF(Tens(JA,JB,JC).EQ.ZERO) CYCLE
                      Dum=Dum+Tens(IA,IB,IC)*Tens(JA,JB,JC)*ITheta%D(IA,JA) &
                           & *WCarts%D(IB,I)*WCarts%D(JB,J)
                   ENDDO
                   ENDDO
                ENDDO
                ENDDO
                P%D(KK,LL)=Dum
                ! Translation.
                IF(IC.EQ.JC) P%D(KK,LL)=P%D(KK,LL)+SMA*SMB/MssTot
                ! Symmetrize.
                P%D(LL,KK)=P%D(KK,LL)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    ! Compute (1-P).
    CALL DSCAL(N**2,-1.0D0,P%D(1,1),1)
    DO I=1,N
       P%D(I,I)=1.0D0+P%D(I,I)
    ENDDO
    !CALL Print_DBL_RNK2(P,'Projector',Unit_O=6)
    !
    ! Let's projecte out.
    CALL New(Tmp,(/N,N/))
    CALL DGEMM('N','N',N,N,N,1D0,A%D(1,1),N,  P%D(1,1),N,0.0D0,Tmp%D(1,1),N)
    CALL DGEMM('T','N',N,N,N,1D0,P%D(1,1),N,Tmp%D(1,1),N,0.0D0,  A%D(1,1),N)
    CALL Delete(Tmp)
    !
    ! Clean up.
    CALL Delete(P)
    CALL Delete(Theta)
    CALL Delete(ITheta)
    CALL Delete(WCarts)
    !
  END SUBROUTINE ProjOut
  !
  SUBROUTINE SymmMtrx(A,N)
    !-------------------------------------------------------------------
    TYPE(DBL_RNK2) :: A
    INTEGER :: N
    !-------------------------------------------------------------------
    INTEGER :: I,J
    !-------------------------------------------------------------------
    DO J=1,N
       DO I=J+1,N
          A%D(I,J)=(A%D(J,I)+A%D(I,J))*0.5D0
          A%D(J,I)=A%D(I,J)
       ENDDO
    ENDDO
  END SUBROUTINE SymmMtrx
  !
  SUBROUTINE Inertia(Carts,AtMss,NAtoms,TInrt,TEig,TVec)
    !-------------------------------------------------------------------
    TYPE(DBL_RNK2) :: Carts
    TYPE(DBL_VECT) :: AtMss
    INTEGER        :: NAtoms
    REAL(DOUBLE)   :: TInrt(3,3),TVec(3,3),TEig(3)
    !-------------------------------------------------------------------
    INTEGER        :: AtA,I,Info,LWork
    REAL(DOUBLE)   :: MssTot,Mss,XC,YC,ZC,Work(20),CMss(3),CTmp(3,NAtoms)
    !-------------------------------------------------------------------
    !
    CALL DCOPY(3*NAtoms,Carts%D(1,1),1,CTmp(1,1),1)
    !
    ! Center of mass.
    MssTot=0.0D0
    CMss(1)=0.0D0
    CMss(2)=0.0D0
    CMss(3)=0.0D0
    DO AtA=1,NAtoms
       MssTot=MssTot+AtMss%D(AtA)
       DO I=1,3
          CMss(I)=CMss(I)+AtMss%D(AtA)*CTmp(I,AtA)
       ENDDO
    ENDDO
    DO I=1,3
       CMss(I)=CMss(I)/MssTot
    ENDDO
    !
    ! Weight coordinates.
    DO AtA=1,NAtoms
       DO I=1,3
          CTmp(I,AtA)=CTmp(I,AtA)-CMss(I)
       ENDDO
    ENDDO
    !
    CALL DBL_VECT_EQ_DBL_SCLR(9,TVec(1,1),0.0D0)
    DO AtA=1,NAtoms
       Mss=AtMss%D(AtA)
       XC=CTmp(1,AtA)
       YC=CTmp(2,AtA)
       ZC=CTmp(3,AtA)
       TVec(1,1)=TVec(1,1)+Mss*(YC*YC+ZC*ZC)
       TVec(2,2)=TVec(2,2)+Mss*(XC*XC+ZC*ZC)
       TVec(3,3)=TVec(3,3)+Mss*(XC*XC+YC*YC)
       TVec(2,1)=TVec(2,1)-Mss*XC*YC
       TVec(3,1)=TVec(3,1)-Mss*XC*ZC
       TVec(3,2)=TVec(3,2)-Mss*YC*ZC
    ENDDO
    TVec(1,2)=TVec(2,1)
    TVec(1,3)=TVec(3,1)
    TVec(2,3)=TVec(3,2)
    CALL DCOPY(9,TVec(1,1),1,TInrt(1,1),1)
    LWORK=20
    CALL DSYEV('V','U',3,TVec(1,1),3,TEig(1),Work(1),LWORK,Info)
    IF(Info/=SUCCEED)CALL Halt('DSYEV flaked in Inertia. INFO='//TRIM(IntToChar(Info)))
    !
  END SUBROUTINE Inertia
  !
  SUBROUTINE ThermoStat(Carts,AtMss,Freqs,NAtoms,Multp)
    !-------------------------------------------------------------------
    TYPE(DBL_RNK2) :: Carts
    TYPE(DBL_VECT) :: AtMss,Freqs
    INTEGER        :: NAtoms,Multp
    !-------------------------------------------------------------------
    REAL(DOUBLE)   :: MssTot,R,T,KT,P,Sigma,Dum,Fac,U,ExpU,ExpMU
    REAL(DOUBLE)   :: QElec,UElec,HElec,CvElec,CpElec,SElec,GElec,AElec
    REAL(DOUBLE)   :: QTran,UTran,HTran,CvTran,CpTran,STran,GTran,ATran
    REAL(DOUBLE)   :: QRot ,URot ,HRot ,CvRot ,CpRot ,SRot ,GRot ,ARot
    REAL(DOUBLE)   :: QVib ,UVib ,HVib ,CvVib ,CpVib ,SVib ,GVib ,AVib
    REAL(DOUBLE)   :: QTot ,UTot ,HTot ,CvTot ,CpTot ,STot ,GTot ,ATot
    REAL(DOUBLE)   :: TInrt(3,3),TEig(3),TEVec(3,3)
    INTEGER        :: I
    !-------------------------------------------------------------------
    REAL(DOUBLE), PARAMETER :: C_Planck=6.626176D-34
    REAL(DOUBLE), PARAMETER :: C_Boltzmann=1.380662D-23
    REAL(DOUBLE), PARAMETER :: C_Light=2.99792458D+08
    REAL(DOUBLE), PARAMETER :: JouleToCal=4.184D+00
    !
    real(double), parameter :: TOKCAL=627.51D+00
    real(double), parameter :: TOCM=2.194746D+05
    !
    T=298.15D0
    KT=C_Boltzmann*T
    R=C_Avogadro*C_Boltzmann
    P=1.01325D+05
    MssTot=SUM(AtMss%D(1:NAtoms))
    !
    ! Electronic.
    ! Since the excitation to excited states is unknown, and we also
    ! know nothing about spin-orbit splittings, We assume that only
    ! one electronic state contributes to the partition function.
    ! In addition, since we don't know what the spatial degeneracy of
    ! this state is, we take the electronic degeneratcy to be just
    ! the spin multiplicity.
    QElec = DBLE(Multp)
    UElec = ZERO
    HElec = ZERO
    CvElec= ZERO
    CpElec= ZERO
    SElec = R*LOG(QElec)
    GElec = HElec-T*Selec
    AElec = UElec-T*Selec
    !
    ! Translation.
    QTran = (TwoPi*(1D-3*MssTot/C_Avogadro/C_Planck)*(KT/C_Planck))**(3D0/2D0)*KT/P
    UTran = 3D0/2D0*R*T
    HTran = 5D0/2D0*R*T
    CvTran= 3D0/2D0*R
    CpTran= 5D0/2D0*R
    STran = R*(LOG(QTran)+5D0/2D0)
    GTran = HTran-T*STran
    ATran = UTran-T*STran
    !
    ! Rotation.
    CALL Inertia(Carts,AtMss,NAtoms,TInrt,TEig,TEVec)
    ! Check for atoms or linear molecules.
    IF(ABS(TEig(1)*TEig(2)*TEig(3)).LT.1D-3) THEN
       WRITE(*,*) 'Atoms or linear molecules are not supported yet'
       STOP 'Err:ThermoStat'
    ENDIF
    !
    !Should be changed to account the molecular symmetry.
    WRITE(*,*) 'WARNING WARNING WARNING WARNING'
    WRITE(*,*) ' In ThermoStat, Sigma is set to 1'
    WRITE(*,*) ' then change it if needed.'
    WRITE(*,*) 'WARNING WARNING WARNING WARNING'
    Sigma=1.0D+00
    !
    Fac=BohrsToAngstroms**2/C_Avogadro
    TEig(1)=Fac*TEig(1)
    TEig(2)=Fac*TEig(2)
    TEig(3)=Fac*TEig(3)
    !
    Dum  = 8D0*Pi2*(KT/C_Planck)*(1.0D-23/C_Planck)
    QRot = SQRT(Pi*(Dum*TEig(1))*(Dum*TEig(2))*(Dum*TEig(3)))/Sigma
    URot = 3D0/2D0*R*T
    CvRot= 3D0/2D0*R
    CpRot= 3D0/2D0*R
    SRot = R*(LOG(QRot)+3D0/2D0)
    HRot = URot
    GRot = HRot-T*SRot
    ARot = URot-T*SRot
    !
    ! Vibration.
    Dum = ZERO
    DO I=1,3*NAtoms
       IF(Freqs%D(I).LT.1D-2) CYCLE
       Dum=Dum+Freqs%D(I)
    ENDDO
    Dum=Dum*TOKCAL*JouleToCal/(TOCM*TWO)
    !
    QVib = ONE
    UVib = ZERO
    CvVib= ZERO
    SVib = ZERO
    DO I=1,3*NAtoms
       IF(Freqs%D(I).LT.1D-2) CYCLE
       U = (C_Planck*Freqs%D(I)*C_Light*1.0D+02)/KT
       ExpU=1D+35
       IF(U.LT.80D0) ExpU = EXP(U)
       ExpMU= EXP(-U)
       QVib = QVib/(ONE-ExpMU)
       UVib = UVib+R*T*U/(ExpU-ONE)
       CvVib= CvVib+R*U*U*(ExpU/(ExpU-ONE))/(ExpU-ONE)
       SVib = SVib+R*(U/(ExpU-ONE)-LOG(ONE-ExpMU))
    ENDDO
    CpVib= CvVib
    UVib = UVib+Dum*1D3
    HVib = UVib
    GVib = HVib-T*SVib
    AVib = UVib-T*SVib
    !
    ! Sum up the contributions.
    QTot = QElec*QTran*QRot*QVib
    !
    UElec= 1D-3*UElec
    HElec= 1D-3*HElec
    GElec= 1D-3*GElec
    AElec= 1D-3*AElec
    !
    UTran= 1D-3*UTran
    HTran= 1D-3*HTran
    GTran= 1D-3*GTran
    ATran= 1D-3*ATran
    !
    URot = 1D-3*URot
    HRot = 1D-3*HRot
    GRot = 1D-3*GRot
    ARot = 1D-3*ARot
    !
    UVib = 1D-3*UVib
    HVib = 1D-3*HVib
    GVib = 1D-3*GVib
    AVib = 1D-3*AVib
    !
    UTot = UElec+UTran+URot+UVib
    HTot = HElec+HTran+HRot+HVib
    GTot = GElec+GTran+GRot+GVib
    ATot = AElec+ATran+ARot+AVib
    STot = SElec+STran+SRot+SVib
    CvTot= CvElec+CvTran+CvRot+CvVib
    CpTot= CpElec+CpTran+CpRot+CpVib

    WRITE(*,  *)
    WRITE(*,907) 
    WRITE(*,908) 'Elec ',QElec,LOG(QElec)
    WRITE(*,908) 'Trans',QTran,LOG(QTran)
    WRITE(*,908) 'Rot  ',QRot ,LOG(QRot )
    WRITE(*,908) 'Vib  ',QVib ,LOG(QVib )
    WRITE(*,908) 'Tot  ',QTot ,LOG(QTot )
    WRITE(*,  *)
    WRITE(*,909) 
    WRITE(*,910) '    KJ/Mol','    KJ/Mol','    KJ/Mol','    J/MolK','    J/MolK','    J/MolK'
    WRITE(*,911) 'Elec ',UElec,HElec,GElec,CvElec,CpElec,SElec
    WRITE(*,911) 'Trans',UTran,HTran,GTran,CvTran,CpTran,STran
    WRITE(*,911) 'Rot  ',URot ,HRot ,GRot ,CvRot ,CpRot ,SRot
    WRITE(*,911) 'Vib  ',UVib ,HVib ,GVib ,CvVib ,CpVib ,SVib
    WRITE(*,911) 'Tot  ',UTot ,HTot ,GTot ,CvTot ,CpTot ,STot

907 FORMAT(/1X,14X,'Q',15X,'ln(Q)')
908 FORMAT(1X,A6,1P,E15.5,0P,F15.6)
909 FORMAT(/14X,'E',9X,'H',9X,'G',9X,'CV',8X,'CP',8X,'S')
910 FORMAT(1X,6X,6A10)
911 FORMAT(1X,A6,6F10.3)
  END SUBROUTINE ThermoStat
  !
  SUBROUTINE CompIRInt()
    ! Compute IR intensities.
    ! Project the dipole derivative onto each normal mode,
    ! and take the square of the norm of this 3 component vector.
  END SUBROUTINE CompIRInt
  !
  SUBROUTINE CompRAMANInt()
    ! Compute Raman intensities.
    ! The formulae can be found in
    ! Vibrational Intensities in IR and Raman spectroscopy,
    ! W.B.Person and G.Zerbi, Eds., Elsevier NY 1982, Page 23.
  END SUBROUTINE CompRAMANInt
  !
  SUBROUTINE Sayvetz()
  END SUBROUTINE Sayvetz
  !
END MODULE Optimizer
