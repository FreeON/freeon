!    DRIVER ROUTINES FOR DYNAMICS, OPTIMIZATION AND TS METHODS
!    Author:  Matt Challacombe and C. J. Tymczak
!------------------------------------------------------------------------------
MODULE GeomOpt
  USE DerivedTypes
  USE GlobalScalars
#ifdef NAG
  USE F90_UNIX_PROC
#endif
  USE DrvFrcs
  USE PrettyPrint
  USE SCFLocals
  USE Overlay
  USE ParsingKeys
  USE Functionals
  USE DrvSCFs
  USE ProcessControl    
  USE InOut
  USE AtomPairs
  USE LinALg
  USE InCoords
  IMPLICIT NONE 
  CONTAINS
!========================================================================================
!
!========================================================================================
    SUBROUTINE GeOp(Ctrl)
       TYPE(SCFControls)         :: Ctrl
       TYPE(CRDS)                :: GM
       INTEGER                   :: JGeo,KS
       REAL(DOUBLE)              :: RMSGrad,MaxGrad
       REAL(DOUBLE),DIMENSION(3) :: PrevState,CrntState
!-------------------------------------------------------
       CALL SetGlobalCtrlIndecies(Ctrl)           
       IF(CBas==1)THEN
!         Initialize file for geometry updates     
          CALL OpenASCII(GeoFile,Geo,NewFile_O=.TRUE.)
          CLOSE(UNIT=Geo,STATUS='KEEP')
       ENDIF
!  Compute starting energy and forces for the 
!  actual basis set and geometry
#ifdef MMech
       IF(HasQM()) Then
         IF(HasMM()) CALL MM_ENERG(Ctrl)
         CALL OneSCF(Ctrl)
         CALL Forces(Ctrl)
       ENDIF
!
       IF(MMOnly()) Then
         Ctrl%Current=(/0,0,1/)
         CALL SetGlobalCtrlIndecies(Ctrl)           
         CALL MM_ENERG(Ctrl)
       ENDIF
#else
       CALL OneSCF(Ctrl)
       CALL Forces(Ctrl)
#endif
!  Print first configuration
#ifdef MMech
       IF(QMOnly()) THEN
         CALL Get(GM,CurGeom)
         CALL Get(GM%ETotal,'Etot',Tag_O=StatsToChar(Ctrl%Current))
         GM%Confg=CGeo
         CALL PPrint(GM,GeoFile,Geo,'PDB')
         CALL Delete(GM)
       ELSE 
         CALL Get(GM_MM,'GM_MM'//CurGeom)
         CALL Get(GM_MM%ETotal,'Etot',Tag_O=StatsToChar(Ctrl%Current))
         GM_MM%Confg=CGeo
         CALL PPrint(GM_MM,GeoFile,Geo,'PDB')
         CALL Delete(GM_MM)
       ENDIF
#else
       CALL Get(GM,CurGeom)
       CALL Get(GM%ETotal,'Etot',Tag_O=StatsToChar(Ctrl%Current))
       GM%Confg=CGeo
       CALL PPrint(GM,GeoFile,Geo,'PDB')
       CALL Delete(GM)
#endif
!  First Hessian update
!      CtrlVect=SetCtrlVect(Ctrl,'BFGSHess') !!!no CtrlVect definition was given before, in front of Invoke('BFGSHess'
       CALL Invoke('BFGSHess',CtrlVect)
       PrevState=Ctrl%Current
       Ctrl%Previous=Ctrl%Current
       Ctrl%Current=(/0,CBas,CGeo+1/) ! step geometry count
       CrntState=Ctrl%Current
       CALL SetGlobalCtrlIndecies(Ctrl)           
       CALL Put(One,'StepSize')
!  Take a step 
       CtrlVect=SetCtrlVect(Ctrl,'QuNew',NoAdvance_O=.TRUE.)
       CALL Invoke('NewStep',CtrlVect)
!  Loop over geometries
       JGeo=CGeo
       DO WHILE(JGeo<=Ctrl%NGeom)
!  Evaluate energies and forces
#ifdef MMech
       IF(HasQM()) Then
         IF(HasMM()) CALL MM_ENERG(Ctrl)
         CALL OneSCF(Ctrl,Sum_O=.FALSE.)
         Ctrl%Previous=PrevState
         CALL Forces(Ctrl)
       ENDIF
!
       IF(MMOnly()) Then
         Ctrl%Current=(/0,0,JGeo/)
         CALL SetGlobalCtrlIndecies(Ctrl)           
         CALL MM_ENERG(Ctrl)
       ENDIF
#else
       CALL OneSCF(Ctrl,Sum_O=.FALSE.)
!  Set to last succesful (minimizing) state
       Ctrl%Previous=PrevState
!  Evaluate forces 
       CALL Forces(Ctrl)
#endif
          KS=KeepStep(Ctrl)
!  Decide if we will keep this step, and check convergence
          IF(KS==1)THEN
!  Update Hessian
             CALL Invoke('BFGSHess',CtrlVect)
!  Write SCF statistics for minimizing geometry
#ifdef MMech
             IF(HasQM()) CALL SCFSummry(Ctrl)
#else
             CALL SCFSummry(Ctrl)
#endif
!  Update status
             JGeo=JGeo+1
             Ctrl%Previous=Ctrl%Current
             PrevState=Ctrl%Current
             Ctrl%Current=(/0,CBas,JGeo/)
             CrntState=Ctrl%Current 
             CALL SetGlobalCtrlIndecies(Ctrl)           
          ELSEIF(KS==0)THEN
             Ctrl%Current=CrntState
             Ctrl%Previous=PrevState
             CALL SetGlobalCtrlIndecies(Ctrl)    
          ELSEIF(KS==-1)THEN
             Ctrl%Previous=Ctrl%Current
             RETURN
          ENDIF
!  Take another step 
          CtrlVect=SetCtrlVect(Ctrl,'QuNew',NoAdvance_O=.TRUE.)
          CALL Invoke('NewStep',CtrlVect)
       ENDDO
       CLOSE(Geo)
   END SUBROUTINE GeOp
!==============================================================================
!     Simple, very naive backtracking (line search)
!==============================================================================
      FUNCTION KeepStep(Ctrl)
         INTEGER            :: KeepStep
         TYPE(SCFControls)  :: Ctrl
         TYPE(CRDS)         :: GM
         INTEGER            :: AccL
         REAL(DOUBLE)       :: StepSz,OldStep,E0,E1,ET,RelErrE,EUncert, &
                               RMSGrad,MAXGrad,RMSDisp,MaxDisp
         LOGICAL            :: ECnvrgd,GCnvrgd,XCnvrgd
         CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
!---------------------------------------------------------------------------
!        Open the InfFile
!        Get the current step size
         CALL Get(StepSz,'StepSize')
         OldStep=StepSz
!        Get the previous and current total energy
         CALL Get(E0,'Etot',Tag_O=StatsToChar(Ctrl%Previous))
         CALL Get(E1,'Etot',Tag_O=StatsToChar(Ctrl%Current))
         CALL Get(RMSDisp,'RMSDisp',CurGeom)
         CALL Get(MaxDisp,'MaxDisp',CurGeom)
         CALL Get(RMSGrad,'RMSGrad',CurGeom)
         CALL Get(MaxGrad,'MaxGrad',CurGeom)
!         WRITE(*,*)'E'//TRIM(StatsToChar(Ctrl%Previous))//' = '//TRIM(DblToChar(E0))
!         WRITE(*,*)'E'//TRIM(StatsToChar(Ctrl%Current))//' = '//TRIM(DblToChar(E1))
!        Check for decreasing energies and need for backtracking
         RelErrE=(E1-E0)/ABS(E0)
!        Check for convergence
         AccL=Ctrl%AccL(CBas)
         ET=ETol(AccL)
         ECnvrgd=ABS(RelErrE)<1D1*ET.AND.StepSz==One
         GCnvrgd=RMSGrad<GTol(AccL).AND.MaxGrad<GTol(AccL)
         XCnvrgd=RMSDisp<XTol(AccL).AND.MaxDisp<XTol(AccL)
         IF(ECnvrgd.OR.(GCnvrgd.AND.XCnvrgd))THEN
            KeepStep=-1
            Mssg=ProcessName('QuNew','Converged #'//TRIM(CurGeom))
          ELSEIF(RelErrE<Zero)THEN
            StepSz=One 
!           Dont back track
            KeepStep=1
            Mssg=ProcessName('QuNew','Descent   #'//TRIM(CurGeom))
         ELSE
!           Decrease step by half (could get fancy here and try interpolation...)
            Mssg=ProcessName('QuNew','Backtrack #'//TRIM(CurGeom))
            StepSz=StepSz*Half
            IF(StepSz<1.D-2) CALL MondoHalt(DRIV_ERROR,                    &
                                            ' Reached max resolution = '   &
                                           //TRIM(DblToShrtChar(RelErrE))  &
                                           //' in minimization of energy')
!           Mark for backtracking
            KeepStep=0
         ENDIF
         CALL Put(StepSz,'StepSize')
         Mssg=TRIM(Mssg)//' dE= '//TRIM(DblToShrtChar(RelErrE))    &
                        //', Grms= '//TRIM(DblToShrtChar(RMSGrad)) &
                        //', Gmax= '//TRIM(DblToShrtChar(MAXGrad)) &
                        //', Xrms= '//TRIM(DblToShrtChar(RMSDisp)) &
                        //', Xmax= '//TRIM(DblToShrtChar(MAXDisp)) &
                        //', Step= '//TRIM(DblToShrtChar(OldStep))

         WRITE(*,*)TRIM(Mssg)             
         CALL OpenASCII(OutFile,Out)
         WRITE(Out,*)TRIM(Mssg)             
         CLOSE(Out)
         IF(ABS(KeepStep)==1)THEN
#ifdef MMech
           IF(QMOnly()) THEN
             CALL Get(GM,Tag_O=CurGeom)
           ELSE IF(HasMM()) THEN
             CALL Get(GM,Tag_O='GM_MM'//CurGeom)
           ENDIF
#else
           CALL Get(GM,CurGeom)
#endif
            GM%Confg=CGeo
            GM%ETotal=E1
            CALL PPrint(GM,GeoFile,Geo,'PDB')
         ENDIF
      END FUNCTION KeepStep
!========================================================================================
!
  SUBROUTINE CALC_GRAD_ONE_FORCE(Ctrl)
    IMPLICIT NONE
    TYPE(SCFControls)     :: Ctrl
    INTEGER :: I,ISet,IGeo
!
#ifdef MMech
   If(HasQM()) Then
#endif
!  Loop first over basis sets 
     DO ISet=1,Ctrl%NSet
        Ctrl%Current=(/0,ISet,1/)
        CALL SetGlobalCtrlIndecies(Ctrl)           
#ifdef MMech
        IF(ISet==1.AND.HasMM()) CALL MM_ENERG(Ctrl)
#endif
        CALL OneSCF(Ctrl)
     ENDDO
!    Calculate the Force (last basis set)
     CALL Forces(Ctrl)
!
     IF(Ctrl%NGeom>1)THEN
!  Go over additional geometries at last basis set
        DO IGeo=2,Ctrl%NGeom
           Ctrl%Current=(/0,Ctrl%NSet,IGeo/)
           CALL SetGlobalCtrlIndecies(Ctrl)           
#ifdef MMech
           IF(HasMM()) CALL MM_ENERG(Ctrl) !!! temporary; overwrites energy terms calcd at prev geoms
#endif
           CALL OneSCF(Ctrl)
           CALL Forces(Ctrl)
        ENDDO
     ENDIF
#ifdef MMech
   EndIf
#endif
!
#ifdef MMech
   If(MMOnly()) Then
     Ctrl%Current=(/0,0,1/)
     CALL SetGlobalCtrlIndecies(Ctrl)           
     CALL MM_ENERG(Ctrl)
   ENDIF
#endif
  END SUBROUTINE CALC_GRAD_ONE_FORCE 
!--------------------------------------------------------
  SUBROUTINE CALC_GRAD_QNEW_OPT(Ctrl)
     IMPLICIT NONE
     TYPE(SCFControls)     :: Ctrl
     INTEGER :: I,ISet,IGeo
!
        CALL SetGlobalCtrlIndecies(Ctrl)           
!
#ifdef MMech
     IF(MMOnly()) THEN
        Ctrl%Current=(/0,0,CGeo/)
        CALL SetGlobalCtrlIndecies(Ctrl)           
        CALL GeOp(Ctrl)
     ELSE
#endif
     DO ISet=1,Ctrl%NSet
!       Optimize geometry for each basis set
        Ctrl%Current=(/0,ISet,CGeo/)
        CALL SetGlobalCtrlIndecies(Ctrl)           
        CALL GeOp(Ctrl)
     ENDDO
#ifdef MMech
     ENDIF
#endif
  END SUBROUTINE CALC_GRAD_QNEW_OPT
!------------------------------------------------------------
  SUBROUTINE CALC_GRAD_QNEW_ONE_OPT(Ctrl)
     IMPLICIT NONE
     TYPE(SCFControls)     :: Ctrl
     INTEGER :: I,ISet,IGeo
!
!optimize only in last basis set and just one step
!
#ifdef MMech
   If(HasQM()) Then
#endif
!  Loop first over basis sets 
     DO ISet=1,Ctrl%NSet-1
        Ctrl%Current=(/0,ISet,1/)
        CALL SetGlobalCtrlIndecies(Ctrl)           
#ifdef MMech
        IF(ISet==1.AND.HasMM()) CALL MM_ENERG(Ctrl)
#endif
        CALL OneSCF(Ctrl)
     ENDDO
     Ctrl%Current=(/0,Ctrl%NSet,1/)
     CALL SetGlobalCtrlIndecies(Ctrl)           
     CALL GeOp(Ctrl)
!
!     IF(Ctrl%NGeom>1)THEN
!!  Go over additional geometries at last basis set
!        DO IGeo=2,Ctrl%NGeom
!           Ctrl%Current=(/0,Ctrl%NSet,IGeo/)
!           CALL SetGlobalCtrlIndecies(Ctrl)           
!           IF(HasMM()) CALL MM_ENERG(Ctrl) !!! temporary; overwrites energy terms calcd at prev geoms
!           CALL OneSCF(Ctrl)
!           Ctrl%Current=(/0,Ctrl%NSet,IGeo/)
!           CALL SetGlobalCtrlIndecies(Ctrl)           
!           CALL GeOp(Ctrl)
!        ENDDO
!     ENDIF
#ifdef MMech
   EndIf
#endif
!
#ifdef MMech
   If(MMOnly()) Then
     Ctrl%Current=(/0,0,1/)
     CALL SetGlobalCtrlIndecies(Ctrl)           
     CALL MM_ENERG(Ctrl)
     CALL GeOp(Ctrl)
   ENDIF
#endif
!    Optimize geometry only in last basis set
  END SUBROUTINE CALC_GRAD_QNEW_ONE_OPT
!--------------------------------------------------------
  SUBROUTINE CALC_GRAD_NO_GRAD(Ctrl)
     IMPLICIT NONE
     TYPE(SCFControls)     :: Ctrl
     INTEGER :: I,ISet,IGeo
!
#ifdef MMech
   If(HasQM()) Then
#endif
!  Loop first over basis sets 
     DO ISet=1,Ctrl%NSet
        Ctrl%Current=(/0,ISet,1/)
        CALL SetGlobalCtrlIndecies(Ctrl)           
#ifdef MMech
        IF(ISet==1.AND.HasMM()) CALL MM_ENERG(Ctrl)
#endif
        CALL OneSCF(Ctrl)
     ENDDO
     IF(Ctrl%NGeom>1)THEN
!  Go over additional geometries at last basis set
        DO IGeo=2,Ctrl%NGeom
           Ctrl%Current=(/0,Ctrl%NSet,IGeo/)
           CALL SetGlobalCtrlIndecies(Ctrl)           
#ifdef MMech
           IF(HasMM()) CALL MM_ENERG(Ctrl) !!! temporary; overwrites energy terms calcd at prev geoms, also for MM multiple geomeptries are not available yet
#endif
           CALL OneSCF(Ctrl)
        ENDDO
     ENDIF
#ifdef MMech
   EndIf
#endif
!
#ifdef MMech
   If(MMOnly()) Then
     Ctrl%Current=(/0,0,1/)
     CALL SetGlobalCtrlIndecies(Ctrl)           
     CALL MM_ENERG(Ctrl)
   ENDIF
#endif
  END SUBROUTINE CALC_GRAD_NO_GRAD
!-------------------------------------------------------------
  SUBROUTINE CALC_SinglePoint(Ctrl)
!
! Energy only calculation, at current geometry and basis set,
! as defined in Ctrl vector.
! Ctrl array is supposed to be set correctly at this point.
!
       IMPLICIT NONE
       TYPE(SCFControls)     :: Ctrl
       INTEGER :: I,ISet,IGeo
!
#ifdef MMech
       IF(HasQM()) THEN
!!! temporary; may overwrite MM energy terms calcd at prev geoms 
         IF(HasMM()) CALL MM_ENERG(Ctrl) 
         CALL OneSCF(Ctrl)
           IF(Ctrl%Grad/=GRAD_NO_GRAD) THEN
             CALL Forces(Ctrl)
           ENDIF
       ENDIF
#else
       CALL OneSCF(Ctrl)
         IF(Ctrl%Grad/=GRAD_NO_GRAD) THEN
           CALL Forces(Ctrl)
         ENDIF
#endif
!
#ifdef MMech
       IF(MMOnly()) Then
         CALL MM_ENERG(Ctrl)
       ENDIF
#endif
!
       END SUBROUTINE CALC_SinglePoint 
!-------------------------------------------------------------
  SUBROUTINE OptimizeNew(Ctrl)
     IMPLICIT NONE
     TYPE(SCFControls)     :: Ctrl
     INTEGER :: I,ISet,IGeo
!
! Geometry optimization, starting at CGeo.
! It is assumed, that CGeo is set by Ctrl at this point.
!
! Initialize file for geometry updates     
!
#ifdef MMech
     IF(CBas==1.OR.MMOnly())THEN
       CALL OpenASCII(GeoFile,Geo,NewFile_O=.TRUE.)
       CLOSE(UNIT=Geo,STATUS='KEEP')
     ENDIF
#else
     IF(CBas==1)THEN
       CALL OpenASCII(GeoFile,Geo,NewFile_O=.TRUE.)
       CLOSE(UNIT=Geo,STATUS='KEEP')
     ENDIF
#endif
!
#ifdef MMech
     IF(MMOnly()) THEN
         CALL GeOpNew(Ctrl)
     ELSE
!
!    Optimize geometry for each basis set
!
       DO ISet=1,Ctrl%NSet
         Ctrl%Current=(/0,ISet,CGeo/)
         CALL SetGlobalCtrlIndecies(Ctrl)
         CALL GeOpNew(Ctrl)
       ENDDO
     ENDIF
#else
!
!    Optimize geometry for each basis set
!
       DO ISet=1,Ctrl%NSet
         Ctrl%Current=(/0,ISet,CGeo/)
         CALL SetGlobalCtrlIndecies(Ctrl)
         CALL GeOpNew(Ctrl)
       ENDDO
#endif
  END SUBROUTINE OptimizeNew
!------------------------------------------------------------
    SUBROUTINE GeOpNew(Ctrl)
       TYPE(SCFControls)         :: Ctrl
       TYPE(CRDS)                :: GMLoc
       INTEGER                   :: JGeo,KS,BlkGeomSize,NCoordsEff,NCoordsEff0
       INTEGER                   :: Refresh,I,J,K,L,II,MaxGeOpSteps
       INTEGER                   :: NCart,NIntC
       REAL(DOUBLE)              :: RMSGrad,MaxGrad
       REAL(DOUBLE)              :: TrixThresh,AInvDistanceThresh
       REAL(DOUBLE)              :: BondConvCrit,AngleConvCrit,GradConvCrit
       REAL(DOUBLE)              :: MaxBondDispl,MaxAngleDispl,RMSIntDispl
       REAL(DOUBLE)              :: Etot
       REAL(DOUBLE),DIMENSION(3) :: PrevState,CrntState
       TYPE(INTC)                :: IntCs
       TYPE(DBL_VECT)            :: IntOld,Displ
       INTEGER                   :: GDIISMemory,FirstGeom
       INTEGER                   :: AccL
       LOGICAL                   :: GCnvrgd     
       LOGICAL                   :: LineSearchOn
       CHARACTER(LEN=DEFAULT_CHR_LEN) :: GMTag
!-------------------------------------------------------
!
       CALL SetGlobalCtrlIndecies(Ctrl)           
!
! Get some thresholds
!
         AInvDistanceThresh=1.D3
         GMTag=''
#ifdef MMech
       IF(HasQM()) THEN
         CALL Get(TrixThresh,'trixneglect',Tag_O=CurBase)
       ELSE
         GMTag='GM_MM'
         CALL Get(TrixThresh,'trixneglect',Tag_O=IntToChar(1))
       ENDIF
#else
         CALL Get(TrixThresh,'trixneglect',Tag_O=CurBase)
#endif
!
! Get GM data base
!
       CALL Get(GMLoc,TRIM(GMTag)//CurGeom)
       GMLoc%Confg=CGeo
!
! Set Convergence Criteria
!
       AccL=Ctrl%AccL(CBas)
       LineSearchOn=.FALSE.
       BondConvCrit=0.0005D0 ! in Angstroems
       AngleConvCrit=0.5D0 ! in degrees   
!
! Define maximum number of optimization steps
!
       NCart=3*GMLoc%Natms
       MaxGeOpSteps=MAX(NCart,200)
!
! Set GDIISmemory, it tells the number of 'remembered' steps in GDIIS
!
       CALL OpenASCII(OutFile,Out)
       IF(Ctrl%Geop%CoordType==CoordType_Cartesian) THEN
         WRITE(*,*) '*****************************************'
         WRITE(*,*) 'Geometry Optimization in Cartesian Coords'
         WRITE(*,*) '*****************************************'
         WRITE(Out,*) '*****************************************'
         WRITE(Out,*) 'Geometry Optimization in Cartesian Coords'
         WRITE(Out,*) '*****************************************'
       ELSE IF(Ctrl%GeOp%Coordtype==CoordType_PrimInt) THEN
         WRITE(*,*) '****************************************'
         WRITE(*,*) 'Geometry Optimization in Internal Coords'
         WRITE(*,*) '****************************************'
         WRITE(Out,*) '****************************************'
         WRITE(Out,*) 'Geometry Optimization in Internal Coords'
         WRITE(Out,*) '****************************************'
       ENDIF
       CLOSE(Out)
!
! Define parameters of block operations for the optimizer,
! this will be the size of the vectors and matrices in the optimizer
!
!!!!       BlkGeomSize=3
!!!!       CALL Put(BlkGeomSize,'BlkGeomSize')
!
! Parametrize block matrix algebra for different cases
!
!!!!!  IF(Ctrl%Geop%CoordType==CoordType_Cartesian) THEN
!!!!!    NCoordsEff0=3*GMLoc%Natms
!!!!!  ELSE IF(Ctrl%GeOp%Coordtype==CoordType_PrimInt) THEN
!!!!!    NCoordsEff0=MAX(NIntC,3*GMLoc%Natms)
!!!!!  ENDIF
!!!!!
!!!!! General block parametrization, once NCoordsEff0 is set.
!!!!!
!!!!         CALL Put(NCoordsEff0,'NCoordsEff0') 
!!!!         NCoordsEff=(NCoordsEff0+BlkGeomSize-1)/BlkGeomSize 
!!!!         CALL Put(NCoordsEff,'NCoordsEff') 
!!!!         Natoms=NCoordsEff !!! number of blocks in a vector
!!!!         NBasF=BlkGeomSize*Natoms
!!!!         MaxAtms=Natoms+1
!!!!         CALL New(BSiz,Natoms)
!!!!         CALL New(OffS,Natoms)
!!!!         BSiz%I=BlkGeomSize
!!!!         OffS%I(1)=1
!!!!         DO I=2,Natoms
!!!!           OffS%I(I)=OffS%I(I-1)+BlkGeomSize
!!!!         ENDDO
!!!!         MaxBlkSize=MAXVAL(BSiz%I)
!!!!         PrintFlags%Mat=DEBUG_MATRICES
!
!
! Start iteration over geometries
!
       FirstGeom=Ctrl%Current(3)
       II=FirstGeom
300    CONTINUE
!
! Check for refresh of int. coord. definitions
!
       IF(II==FirstGeom) THEN !!!! refresh does not work in HDF
         CALL IntCReDef(Ctrl,Refresh)
         IF(Refresh/=0) THEN
           CALL GetIntCs(GMLoc%Carts%D,GMLoc%Natms,InfFile, &
                         IntCs,NIntC,Refresh)
         ENDIF
       ENDIF
       CALL INTCValue(IntCs,GMLoc%Carts%D)
       CALL New(IntOld,NIntC)
       IntOld%D=IntCs%Value
!      CALL PrtIntCoords(IntCs,IntCs%Value, &
!                        'Internals at step '//TRIM(IntToChar(II)))
!
! Get Energy and Cartesian gradient at current geometry
!
       CALL CALC_SinglePoint(Ctrl)
!
! Calculate simple relaxation (SR) step from an inverse Hessian
!
       CALL NewDispl(Ctrl,Displ,NCart,NIntC)
       CALL SRStep(Ctrl,GMLoc,Displ,MaxGrad,RMSGrad,IntCs) 
!
! Calculate new geometry 
!
       CALL NewStructure(Ctrl,GMLoc,Displ,LineSearchOn,IntCs)
       CALL Delete(Displ)
!
! Check convergence
!
       CALL INTCValue(IntCs,GMLoc%Carts%D)
       IntOld%D=IntCs%Value-IntOld%D
       MaxBondDispl=Zero
       MaxAngleDispl=Zero
       DO I=1,NIntC 
          IF(IntCs%Def(I)(1:4)=='STRE') THEN 
            MaxBondDispl=MAX(MaxBondDispl,ABS(IntOld%D(I))) 
          ELSE
            MaxAngleDispl=MAX(MaxAngleDispl,ABS(IntOld%D(I))) 
          ENDIF
       ENDDO
       MaxBondDispl=MaxBondDispl/AngstromsToAU
       MaxAngleDispl=MaxAngleDispl*180.D0/PI
       RMSIntDispl=SQRT(DOT_PRODUCT(IntOld%D,IntOld%D)/DBLE(NIntC))
         CALL Delete(IntOld)
!
       GCnvrgd=RMSGrad<GTol(AccL).AND.MaxGrad<GTol(AccL)
!
! Review iterations
!
       CALL OpenASCII(OutFile,Out)
!
         WRITE(*,399) II-FirstGeom+1
       IF(LineSearchOn) THEN
         CALL Get(Etot,'ETot',StatsToChar(Ctrl%Current))
         WRITE(*,400) Etot
         WRITE(Out,400) Etot
       ELSE
         CALL Get(Etot,'ETot',StatsToChar(Ctrl%Previous))
         WRITE(*,401) Etot
         WRITE(Out,401) Etot
       ENDIF
!
       WRITE(*,410) MaxGrad
       WRITE(*,420) RMSGrad
       WRITE(*,430) MaxBondDispl
       WRITE(*,435) MaxAngleDispl
       WRITE(*,440) RMSIntDispl
       WRITE(Out,410) MaxGrad
       WRITE(Out,420) RMSGrad
       WRITE(Out,430) MaxBondDispl
       WRITE(Out,435) MaxAngleDispl
       WRITE(Out,440) RMSIntDispl
!
       CLOSE(Out,STATUS='KEEP')
!
399 FORMAT('GeOp step= ',I6)
400 FORMAT(' Total Energy at Current Geometry  = ',F20.8)
401 FORMAT('Total Energy at Previous Geometry  = ',F20.8)
410 FORMAT('           ','             MaxGrad= ',F12.6)
420 FORMAT('           ','             RMSGrad= ',F12.6)
430 FORMAT('           ','      Max Bond Displ= ',F12.6)
435 FORMAT('           ','     Max Angle Displ= ',F12.6)
440 FORMAT('           ','           RMS Displ= ',F12.6)
!
! Check convergence
!!!!!!!build in AccL dependent criteria for bonds and angles
!
       IF(II-FirstGeom+1<=MaxGeOpSteps) THEN
         IF(.NOT.GCnvrgd) THEN
!        IF(MaxBondDispl>BondConvCrit .OR. &
!           MaxAngleDispl>AngleConvCrit .OR. &
!           .NOT.GCnvrgd) THEN
           II=II+1
           Ctrl%Current(3)=II
           CALL SetGlobalCtrlIndecies(Ctrl)           
             CALL Put(GMLoc,TRIM(GMTag)//CurGeom)
!#ifdef MMech
!             IF(HasQM()) CALL Halt('put new GM coordinates separately, to be done')
!#else
!             CALL Halt('put new GM coordinates separately, to be done')
!#endif
           GO TO 300
         ENDIF
       ELSE
       ENDIF
!
! Convergence is reached at this point, calculate final energy
! and finish optimization
!
       CALL CALC_SinglePoint(Ctrl)
       CALL Get(Etot,'ETot',StatsToChar(Ctrl%Current))
!
       CALL OpenAscii(OutFile,Out)
!
       WRITE(*,460) II
       WRITE(*,470) Etot
       WRITE(Out,460) II
       WRITE(Out,470) Etot
460  FORMAT('Geometry Optimization converged in ',I6,' steps.')
470  FORMAT('Total Energy at optimum structure= ',F20.8)
!
       CLOSE(Out,STATUS='KEEP')
!
! Tidy up
!
       CALL Delete(IntCs)
       CALL Delete(GMLoc)
!!!!!!      CALL Delete(BSiz)
!!!!!!       CALL Delete(OffS)
!
   END SUBROUTINE GeOpNew
!
!--------------------------------------------------------------------
!
      SUBROUTINE SRStep(Ctrl,GMLoc,Displ,MaxGrad,RMSGrad,IntCs)
!
! Simple Relaxation step
!
      TYPE(SCFControls)              :: Ctrl 
      TYPE(CRDS)                     :: GMLoc
      TYPE(DBL_VECT)                 :: Displ
      REAL(DOUBLE)                   :: MaxGrad,RMSGrad
      TYPE(DBL_VECT)                 :: CartGrad,IntGrad,Grad
      TYPE(INTC)                     :: IntCs
      INTEGER                        :: I,J,NDim
      INTEGER                        :: NatmsLoc,NCart,NIntC
      LOGICAL                        :: DoInternals
!
      NatmsLoc=GMLoc%Natms
      NCart=3*NatmsLoc
      NDim =SIZE(Displ%D)
      DoInternals=.FALSE.
      IF(NCart/=NDim) DoInternals=.TRUE.
!
      CALL New(Grad,NDim)
!
! First, get Cartesian gradient
!
      CALL New(CartGrad,NCart)
      CALL Get(CartGrad,'GradE',Tag_O=CurGeom) 
!
! If requested, compute internal coord. gradients
!
      IF(DoInternals) THEN
        NIntC=SIZE(IntCs%Def)
        IF(NIntC/=NDim) CALL Halt('Dimension error in SRStep')
        CALL New(IntGrad,NDim)
        CALL CartToInternal(GMLoc%Carts%D,IntCs,CartGrad%D,IntGrad%D)
        Grad%D=IntGrad%D
      ELSE
        Grad%D=CartGrad%D
      ENDIF
!
! Steepest descent (may be either internal or Cartesian)
! or other inverse Hessian guess
!
      CALL SteepestDesc(Ctrl,Grad,Displ)
!
! Set constraints on the displacements
!
      CALL SetConstraint(IntCs,GMLoc%Carts%D,Displ)
!
! Check for convergence
!
      MaxGrad=Zero
      DO I=1,NDim ; MaxGrad=MAX(MaxGrad,ABS(Grad%D(I))) ; ENDDO
      RMSGrad=SQRT(DOT_PRODUCT(Grad%D,Grad%D)/DBLE(NDim))
!
! Tidy up
!
      IF(DoInternals) CALL Delete(IntGrad)
      CALL Delete(CartGrad)
      CALL Delete(Grad)
!
      END SUBROUTINE SRStep
!-------------------------------------------------------
      SUBROUTINE NewStructure(Ctrl,GMLoc,Displ,LineSearchOn,IntCs)
!
! Line search. Can be carried out either in internals
! or in Cartesian.
! Here we optimize a single parameter 'Fact',
! to get optimum for the energy at Fact*Displ displacement.
!
      TYPE(SCFControls)              :: Ctrl
      TYPE(CRDS)                     :: GMLoc
      TYPE(INTC)                     :: IntCs
      INTEGER                        :: I,J,II,MaxSteps,NDim,NInTc
      INTEGER                        :: NatmsLoc,NCart
      REAL(DOUBLE)                   :: EStart,Fact
      TYPE(DBL_VECT)                 :: Displ
      TYPE(DBL_RNK2)                 :: ActXYZ
      TYPE(DBL_VECT)                 :: ActDispl
      TYPE(DBL_VECT)                 :: Energy
      CHARACTER(LEN=DEFAULT_CHR_LEN) :: GMTag
      LOGICAL                        :: IntSearch,LineSearchOn
      INTEGER                        :: OldGrad,GDIISMemory
!
! MaxSteps : Maximum number of Line Search Steps
! Displ    : Initial displacement vector from SR step, 
!            may be Cartesian or internal
! ActDispl : Size of actual displacement at the II-th line search step
!            maximum value is ActDispl=MaxFact*Displ
!            MaxFact must be small enough for internal coordinates
!            in order coordinate transformation can converge
!
      MaxSteps=15
!
      GDIISMemory=500
!
      NatmsLoc=GMLoc%Natms
      NCart=3*NatmsLoc   
      NDim =SIZE(Displ%D)
      IntSearch=.FALSE.
      IF(NCart/=NDim) THEN
        IntSearch=.TRUE.
        NIntC =SIZE(IntCs%Def)
        IF(NIntC/=NDim) CALL Halt('Dimensionality error in LineSearch')
      ENDIF
!
! Geometry Tag
!
    GMTag=''
#ifdef MMech
    IF(HasMM()) GMTag='GM_MM'
#endif
!
! Allocate arrays
!
      CALL New(ActDispl,NDim)
      CALL New(ActXYZ,(/3,NatmsLoc/))
      CALL New(Energy,MaxSteps)
!
! Starting energy      
!
      CALL Get(EStart,'Etot',StatsToChar(Current))
!
! Increment Current and prepare for energy only calculation
!
      OldGrad=Ctrl%Grad
      Ctrl%Grad=GRAD_NO_GRAD
      Ctrl%Previous  =Ctrl%Current
      Ctrl%Current(3)=Ctrl%Current(3)+1
      CALL SetGlobalCtrlIndecies(Ctrl)
!
! Initialization
!
      Fact=One
      ActXYZ%D=GMLoc%Carts%D
!
      II=0
100   CONTINUE
      II=II+1
!
      ActDispl%D=Fact*Displ%D
!
! Construct new structure either by Cartesian or by internal
! displacements
!
      IF(IntSearch) THEN 
        CALL INTCValue(IntCs,ActXYZ%D)
        CALL InternalToCart(ActXYZ%D,IntCs,ActDispl%D)
      ELSE
        CALL CartRNK1ToCartRNK2(ActDispl%D,ActXYZ%D,.TRUE.)
      ENDIF
!
! Calculate energy at new structure; CurGeom has a new value now,
! and does not change until new call for LineSearch.
!
      GMLoc%Carts%D=ActXYZ%D
      CALL Put(GMLoc,TRIM(GMTag)//CurGeom)
      IF(.NOT.LineSearchOn) GO TO 200
! Continue here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL CALC_SinglePoint(Ctrl)
      CALL Get(Energy%D(II),'Etot',StatsToChar(Current))
!write(*,*) II,'  energy= ',Energy%D(II)
!
!
!
200   CONTINUE
!
! Set back Ctrl%Grad to its original value
!
      Ctrl%Grad=OldGrad
!
! Now call GDIIS
! and construct GDIIS step from SR step and previous displacements
! May work in either Cartesian or internal displacemets space.
! Subroutine modifies Cartesian Coordinates in GMLoc, now.
!
!     IF(CGeo>2) CALL GDIIS(CGeo,GDIISMemory,GMLoc)
!
!
! Tidy up
!
      CALL Delete(ActDispl)
      CALL Delete(ActXYZ)
      CALL Delete(Energy)
!     IF(IntSearch) CALL Delete(ActIntC)
!
      END SUBROUTINE NewStructure
!
!-----------------------------------------------------------------
!
      SUBROUTINE IntCReDef(Ctrl,Refresh)
      TYPE(SCFControls) :: Ctrl
      INTEGER           :: Refresh
!
! WARNING! refresh may change the number of internal coordinates!
! This may cause problems in HDF!!!!!!
!
         IF(Ctrl%Current(3)==1) THEN
           Ctrl%GeOp%ReDefIntC=1
         ELSE
           Ctrl%GeOp%ReDefIntC=2
         ENDIF
!
         Refresh=Ctrl%GeOp%ReDefIntC
!
      END SUBROUTINE IntCReDef
!------------------------------------------------------------------
!
       SUBROUTINE NewDispl(Ctrl,Displ,NCart,NIntC)
       TYPE(SCFControls) :: Ctrl 
       TYPE(DBL_VECT) :: Displ
       INTEGER        :: NCart,NIntC
!
         IF(Ctrl%Geop%CoordType==CoordType_Cartesian) THEN
           CALL New(Displ,NCart)
           Displ%D(:)=Zero
         ELSE
           CALL New(Displ,NIntC)
           Displ%D(:)=Zero
         ENDIF
!
       END SUBROUTINE NewDispl
!
!-------------------------------------------------------
!
      SUBROUTINE SteepestDesc(Ctrl,Grad,Displ)
!
      TYPE(SCFControls) :: Ctrl
      TYPE(DBL_VECT)    :: Grad,Displ 
!
        IF(Ctrl%GeOp%CoordType==CoordType_Cartesian) THEN
          Displ%D=-1.D0*Grad%D !!!! it oscillates with 2.D0
        ELSE IF(Ctrl%GeOp%CoordType==CoordType_PrimInt) THEN
          Displ%D=-2.D0*Grad%D 
        ELSE
          Displ%D=-5.D0*Grad%D 
        ENDIF 
!
      END SUBROUTINE SteepestDesc
!-------------------------------------------------------
!
 END MODULE GeomOpt
