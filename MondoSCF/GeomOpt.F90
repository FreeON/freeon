!    DRIVER ROUTINES FOR DYNAMICS, OPTIMIZATION AND TS METHODS
!    Author:  Matt Challacombe and C. J. Tymczak
!------------------------------------------------------------------------------
MODULE GeomOpt
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalObjects
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
!
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
!
!-------------------------------------------------------------
!
       SUBROUTINE CALC_SinglePoint(Ctrl)
!
! Energy and Force calculation, at current geometry and basis set,
! as defined in Ctrl vector.
!
       IMPLICIT NONE
       TYPE(SCFControls)     :: Ctrl
       INTEGER :: I,ISet,IGeo,OldGrad
!
#ifdef MMech
       IF(HasQM()) THEN
!!! temporary; may overwrite MM energy terms calcd at prev geoms 
         IF(HasMM()) CALL MM_ENERG(Ctrl) 
         CALL OneSCF(Ctrl)
         CALL Forces(Ctrl)
       ENDIF
#else
         CALL OneSCF(Ctrl)
         CALL Forces(Ctrl)
#endif
!
#ifdef MMech
       IF(MMOnly()) Then
         CALL MM_ENERG(Ctrl)
       ENDIF
#endif
!
       END SUBROUTINE CALC_SinglePoint 
!
!-------------------------------------------------------------
!
       SUBROUTINE CALC_ForceOnly(Ctrl)
!
! Force only calculation, at current geometry and basis set,
! as defined in Ctrl vector.
! It is supposed that energy calculation has already been done.
! In present version MM force and energy cannot be calculated separately
!
       IMPLICIT NONE
       TYPE(SCFControls)     :: Ctrl
       INTEGER :: I,ISet,IGeo,OldGrad
!
#ifdef MMech
       IF(HasQM()) THEN
!!! temporary; may overwrite MM energy terms calcd at prev geoms 
         IF(HasMM()) CALL MM_ENERG(Ctrl) 
         CALL Forces(Ctrl)
       ENDIF
#else
         CALL Forces(Ctrl)
#endif
!
#ifdef MMech
       IF(MMOnly()) Then
         CALL MM_ENERG(Ctrl)
       ENDIF
#endif
!
       END SUBROUTINE CALC_ForceOnly 
!-------------------------------------------------------------
!
       SUBROUTINE CALC_EnergyOnly(Ctrl)
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
       ENDIF
#else
       CALL OneSCF(Ctrl)
#endif
!
#ifdef MMech
       IF(MMOnly()) Then
         CALL MM_ENERG(Ctrl)
       ENDIF
#endif
!
       END SUBROUTINE CALC_EnergyOnly 
!----------------------------------------------------------------
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
       INTEGER                   :: Refresh,I,J,K,L,IStep,MaxGeOpSteps
       INTEGER                   :: NCart,NIntC
       REAL(DOUBLE)              :: RMSGrad,MaxGrad
       INTEGER                   :: BlkGeomSize,NCoordsEff,NCoordsEff0
       REAL(DOUBLE)              :: TrixThresh,AInvDistanceThresh
       REAL(DOUBLE)              :: Etot
       REAL(DOUBLE)              :: StreConvCrit,AngleConvCrit
       REAL(DOUBLE),DIMENSION(3) :: PrevState,CrntState
       TYPE(INTC)                :: IntCs
       TYPE(DBL_VECT)            :: IntOld,Displ,AuxVect
       INTEGER                   :: FirstGeom
       INTEGER                   :: GDIISMaxMem,InitGDIIS
       INTEGER                   :: AccL,ActStep
       LOGICAL                   :: GeOpConvgd
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
       NCart=3*GMLoc%Natms
!
! Set control parameters of geometry optimization
!
       CALL SetGeOpCtrl(Ctrl,GMLoc%Natms)
!
! Initializations.
! Do not use iterative subspace from previous basis calculations!!!!
! That would make convergence slower!
!
       CALL Put(0,'RefMemory')  
       CALL Put(0,'SRMemory')  
       CALL Put(0,'CartGradMemory')  
       CALL Put(0,'IntGradMemory')  
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
! Start iteration over geometries
!
         FirstGeom=Ctrl%Current(3)
         IStep=FirstGeom-1
300    CONTINUE
         IStep=IStep+1
!        Ctrl%Previous(3)=IStep-1
         Ctrl%Current(3)=IStep
         CALL SetGlobalCtrlIndecies(Ctrl)           
         ActStep=IStep-FirstGeom+1
         GeOpCtrl%ActStep=ActStep
!
! Set dynamic changes in control parameters 
!
       CALL SetGeOpCtrl(Ctrl,GMLoc%Natms)
!
! Check for refresh of int. coord. definitions
!
       IF(IStep==FirstGeom) THEN !!!! refresh does not work in HDF
         CALL IntCReDef(Ctrl,Refresh)
         IF(Refresh/=0) THEN
           CALL GetIntCs(GMLoc%Carts%D,GMLoc%Natms,InfFile, &
                         IntCs,NIntC,Refresh)
!                        IntCs,NIntC,4)
         ENDIF
       ENDIF
       IF(NIntC/=0) THEN
         CALL INTCValue(IntCs,GMLoc%Carts%D)
         CALL New(IntOld,NIntC)
         IntOld%D=IntCs%Value
       ENDIF
!
! Compute Energy at current geometry
!
       CALL CALC_SinglePoint(Ctrl)
!      CALL CALC_EnergyOnly(Ctrl)
!
! Print current geometry for debugging
!
       GMLoc%Confg=CGeo
       CALL Get(GMLoc%ETotal,'Etot',Tag_O=StatsToChar(Ctrl%Current))
       CALL PPrint(GMLoc,GeoFile,Geo,'PDB')
       IF(PrintFlags%GeOp==DEBUG_GEOP) THEN
         CALL PrtXYZ(GMLoc%Carts%D,Title_O='geometry at step= ' &
             //TRIM(IntToChar(IStep)),PrtU_O=11,Convert_O=.TRUE.)
         CALL INTCValue(IntCs,GMLoc%Carts%D)
         CALL PrtIntCoords(IntCs, &
           IntCs%Value,'Internals at step #'//TRIM(IntToChar(ActStep)))
       ENDIF
!
! Compute Force at current geometry
!
!      CALL CALC_ForceOnly(Ctrl)
!
! Put current cartesian coords and gradients as a reference for GDIIS!
!
       CALL PutSRStep(XYZ_O=GMLoc%Carts%D,Tag_O='Ref')
!
! Calculate simple relaxation (SR) step from an inverse Hessian
!
       CALL NewDispl(Ctrl,Displ,NCart,NIntC)
       CALL SRStep(Ctrl,GMLoc,Displ,IntCs) 
!
! Calculate new geometry 
!
       CALL NewStructure(Ctrl,GMLoc,Displ,IntCs)
       CALL Delete(Displ)
!
! Check convergence
!
       CALL GeOpConv(Ctrl,GMLoc,IntCs,IntOld)
!
! Save current structure into HDF
!
       CALL Put(GMLoc,TRIM(GMTag)//NxtGeom)
!
! Tidy up and continue optimization if necessary
!
       IF(NIntC/=0)  CALL Delete(IntOld)
       IF(ActStep<=GeOpCtrl%MaxGeOpSteps) THEN
         IF(.NOT.GeOpCtrl%GeOpConvgd) THEN
           GO TO 300
         ENDIF
       ENDIF
!
! Convergence is reached at this point, calculate final energy
! and finish optimization
!
       Ctrl%Current(3)=IStep+1
       CALL SetGlobalCtrlIndecies(Ctrl)           
       CALL CALC_EnergyOnly(Ctrl)
       CALL Get(Etot,'ETot',StatsToChar(Ctrl%Current))
!
       CALL OpenAscii(OutFile,Out)
!
       IF(ActStep<GeOpCtrl%MaxGeOpSteps) THEN
         WRITE(*,460) ActStep
         WRITE(Out,460) ActStep
         WRITE(*,470) Etot
         WRITE(Out,470) Etot
       ELSE
         CALL Halt('Maximum number of optimization steps exceeded, optimization did not converge.')
       ENDIF
       
460  FORMAT('Geometry Optimization converged in ',I6,' steps.')
470  FORMAT('Total Energy at optimum structure= ',F20.8)
!
       CLOSE(Out,STATUS='KEEP')
!
! Tidy up
!
       IF(NIntC/=0) CALL Delete(IntCs)
       CALL Delete(GMLoc)
!!!!!!      CALL Delete(BSiz)
!!!!!!       CALL Delete(OffS)
!
   END SUBROUTINE GeOpNew
!
!--------------------------------------------------------------------
!
      SUBROUTINE SRStep(Ctrl,GMLoc,Displ,IntCs,Step_O,DoSave_O)
!
! Simple Relaxation step
!
        TYPE(SCFControls)              :: Ctrl 
        TYPE(CRDS)                     :: GMLoc
        TYPE(DBL_VECT)                 :: Displ
        REAL(DOUBLE)                   :: MaxGrad,RMSGrad
        INTEGER                        :: IMaxGradNoConstr
        REAL(DOUBLE)                   :: MaxGradNoConstr
        REAL(DOUBLE)                   :: RMSGradNoConstr,Sum
        TYPE(DBL_VECT)                 :: CartGrad,IntGrad,Grad
        TYPE(INTC)                     :: IntCs
        INTEGER                        :: I,J,NDim,IMaxGrad
        INTEGER                        :: NatmsLoc,NCart,NIntC
        LOGICAL                        :: DoInternals
        LOGICAL,OPTIONAL               :: DoSave_O  
        CHARACTER(Len=*),OPTIONAL      :: Step_O
        CHARACTER(LEN=DEFAULT_CHR_LEN) :: TagStep 
!
        NatmsLoc=GMLoc%Natms
        NCart=3*NatmsLoc
        IF(AllocQ(IntCs%Alloc)) THEN
          NIntC=SIZE(IntCs%Def)
        ELSE
          NIntC=0
        ENDIF
          NDim =SIZE(Displ%D)
        DoInternals=GeOpCtrl%DoInternals
        IF(NIntC/=NDim.AND.DoInternals) &
            CALL Halt('Dimensionality error in LineSearch')
!
        IF(PRESENT(Step_O)) THEN
          TagStep=Step_O 
        ELSE
          TagStep=CurGeom
        ENDIF
!
        CALL New(Grad,NDim)
!
! First, get Cartesian gradient
!
        CALL New(CartGrad,NCart)
        CALL Get(CartGrad,'GradE',Tag_O=TRIM(TagStep)) 
!
! Purify Cartesian gradients from rotations and translations
!
        CALL TranslsOff(CartGrad%D)
        CALL RotationsOff(CartGrad%D,GMLoc%Carts%D)
!
        IF(PRESENT(DoSave_O)) THEN
          IF(DoSave_O) THEN
            CALL PutSRStep(Vect_O=CartGrad,Tag_O='CartGrad')
          ENDIF
        ELSE
          CALL PutSRStep(Vect_O=CartGrad,Tag_O='CartGrad')
        ENDIF
!
! If requested, compute internal coord. gradients
!
        IF(DoInternals) THEN
          NIntC=SIZE(IntCs%Def)
          IF(NIntC/=NDim) CALL Halt('Dimension error in SRStep')
          CALL New(IntGrad,NDim)
          CALL CartToInternal(GMLoc%Carts%D,IntCs,CartGrad%D,IntGrad%D)
          Grad%D=IntGrad%D
!         CALL PutSRStep(Vect_O=IntGrad,Tag_O='IntGrad')
!         CALL Put(IntGrad,Tag_O='IntGrad'//NxtGeom)
        ELSE
          Grad%D=CartGrad%D
        ENDIF
!
! Check for gradient-convergence
!
        MaxGrad=Zero
        DO I=1,NDim  
          IF(MaxGrad<ABS(Grad%D(I))) THEN
            IMaxGrad=I
             MaxGrad=ABS(Grad%D(I))
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
            Sum=Grad%D(I)
            IF(MaxGradNoConstr<ABS(Sum)) THEN
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
        SELECT CASE(Ctrl%Grad)
        CASE(GRAD_STPDESC_OPT) 
          CALL SteepestDesc(Ctrl,Grad,Displ,NCart,GMLoc%Carts%D)
        CASE(GRAD_DIAGHESS_OPT) 
          CALL DiagonalHess(Ctrl,Grad,Displ,IntCs,NCart)
!         CALL DiagHessRFO(Ctrl,Grad,Displ,IntCs,NCart)
        END SELECT
!
! Set constraints on the displacements
!
        CALL SetConstraint(IntCs,GMLoc%Carts%D,Displ)
!
! Tidy up
!
        IF(DoInternals) CALL Delete(IntGrad)
        CALL Delete(CartGrad)
        CALL Delete(Grad)
!
        GeOpCtrl%MaxGrad=MaxGrad
        GeOpCtrl%MaxGradNoConstr=MaxGradNoConstr
        GeOpCtrl%IMaxGrad=IMaxGrad
        GeOpCtrl%RMSGrad=RMSGrad
        GeOpCtrl%RMSGradNoConstr=RMSGradNoConstr
        GeOpCtrl%IMaxGradNoConstr=IMaxGradNoConstr
!
      END SUBROUTINE SRStep
!-------------------------------------------------------
      SUBROUTINE NewStructure(Ctrl,GMLoc,Displ,IntCs,DoSave_O)
!
        TYPE(SCFControls)              :: Ctrl
        TYPE(CRDS)                     :: GMLoc
        TYPE(INTC)                     :: IntCs
        INTEGER                        :: I,J,II,NDim,NInTc
        INTEGER                        :: NatmsLoc,NCart,InitGDIIS
        REAL(DOUBLE)                   :: EStart,Fact
        TYPE(DBL_VECT)                 :: Displ
        CHARACTER(LEN=DEFAULT_CHR_LEN) :: GMTag
        LOGICAL                        :: DoInternals
        LOGICAL,OPTIONAL               :: DoSave_O
        INTEGER                        :: OldGrad
!
! In the present version there is no line search, only GDIIS
!
        InitGDIIS=GeOpCtrl%GDIISInit
!
        NatmsLoc=GMLoc%Natms
        NCart=3*NatmsLoc   
        NDim =SIZE(Displ%D)
        NIntC =SIZE(IntCs%Def)
        DoInternals=GeOpCtrl%DoInternals
        IF(NIntC/=NDim.AND.DoInternals) &
            CALL Halt('Dimensionality error in LineSearch')
!
! Construct new structure either by Cartesian or by internal
! displacements
!
        IF(DoInternals) THEN 
          CALL INTCValue(IntCs,GMLoc%Carts%D)
          CALL InternalToCart(GMLoc%Carts%D,IntCs,Displ%D)
        ELSE
          CALL CartRNK1ToCartRNK2(Displ%D,GMLoc%Carts%D,.TRUE.)
        ENDIF
!
! Save simple relaxation geometry (Cartesians) into HDF
!
        IF(PRESENT(DoSave_O)) THEN
          IF(DoSave_O) CALL PutSRStep(XYZ_O=GMLoc%Carts%D,Tag_O='SR')
        ELSE
          CALL PutSRStep(XYZ_O=GMLoc%Carts%D,Tag_O='SR')
        ENDIF
!
! GDIIS optimization acceleration
!
        IF((.NOT.GeOpCtrl%NoGDIIS).AND.&
          (CGeo>InitGDIIS.AND.GeopCtrl%GDIISOn)) THEN
            CALL GDIIS2(GMLoc)
        ENDIF
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
      SUBROUTINE SteepestDesc(Ctrl,Grad,Displ,NCart,XYZ)
!
      TYPE(SCFControls)           :: Ctrl
      TYPE(DBL_VECT)              :: Grad,Displ 
      INTEGER                     :: NCart
      REAL(DOUBLE),DIMENSION(:,:) :: XYZ
!
        IF(Ctrl%GeOp%CoordType==CoordType_Cartesian) THEN
          Displ%D=-1.D0*Grad%D !!!! it oscillates with 2.D0
          !!! take translations and rotations off
          CALL TranslsOff(Displ%D)
          CALL RotationsOff(Displ%D,XYZ)
        ELSE IF(Ctrl%GeOp%CoordType==CoordType_PrimInt) THEN
          Displ%D=-2.D0*Grad%D 
          CALL RedundancyOffFull(Displ%D,NCart)
        ELSE
          Displ%D=-5.D0*Grad%D 
        ENDIF 
!
      END SUBROUTINE SteepestDesc
!
!-------------------------------------------------------
!
      SUBROUTINE DiagonalHess(Ctrl,Grad,Displ,IntCs,NCart)
!
      TYPE(SCFControls) :: Ctrl
      TYPE(DBL_VECT)    :: Grad,Displ 
      TYPE(INTC)        :: IntCs       
      INTEGER           :: I,J,NIntC,NCart
      REAL(DOUBLE)      :: HStre,HBend,HLinB,HOutP,HTors
!
        NIntC=SIZE(IntCs%Def)
!
        IF(Ctrl%GeOp%CoordType==CoordType_Cartesian) THEN
          CALL OpenAscii(OutFile,Out)
          WRITE(*,*) 'WARNING! DIAGONAL HESSIAN FOR '// &
                       'CARTESIANS IS JUST STEEPEST DESCENT'
          WRITE(Out,*) 'WARNING! DIAGONAL HESSIAN FOR '// &
                       'CARTESIANS IS JUST STEEPEST DESCENT'
          CLOSE(Out,STATUS='KEEP')      
          Displ%D=-1.D0*Grad%D !!!! it oscillates with 2.D0
        ELSE IF(Ctrl%GeOp%CoordType==CoordType_PrimInt) THEN
            HStre=One/GeOpCtrl%StreHessian
            HBend=One/GeOpCtrl%BendHessian
            HLinB=One/GeOpCtrl%LinBHessian
            HOutP=One/GeOpCtrl%OutPHessian
            HTors=One/GeOpCtrl%TorsHessian
          DO I=1,NIntC
            IF(IntCs%Def(I)(1:4)=='STRE') THEN
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
!project out redundancy
          CALL RedundancyOffFull(Displ%D,NCart)
        ELSE
          CALL Halt('Only Primitiv Internals are available yet.')
        ENDIF 
!
      END SUBROUTINE DiagonalHess
!
!-------------------------------------------------------
      SUBROUTINE DiagHessRFO(Ctrl,Grad,Displ,IntCs,NCart)
!
      TYPE(SCFControls) :: Ctrl
      TYPE(DBL_VECT)    :: Grad,Displ
      TYPE(DBL_RNK2)    :: Hessian
      TYPE(INTC)        :: IntCs
      INTEGER           :: I,J,NIntC,NCart,Info
      REAL(DOUBLE)      :: HStre,HBend,HLinB,HOutP,HTors
      REAL(DOUBLE)      :: Sum
!
        NIntC=SIZE(IntCs%Def)
        CALL New(Hessian,(/NIntC+1,NIntC+1/))
        Hessian%D=Zero
!
        IF(Ctrl%GeOp%CoordType==CoordType_Cartesian) THEN
          CALL OpenAscii(OutFile,Out)
          WRITE(*,*) 'WARNING! DIAGONAL HESSIAN FOR '// &
                       'CARTESIANS IS JUST STEEPEST DESCENT'
          WRITE(Out,*) 'WARNING! DIAGONAL HESSIAN FOR '// &
                       'CARTESIANS IS JUST STEEPEST DESCENT'
          CLOSE(Out,STATUS='KEEP')      
          Displ%D=-1.D0*Grad%D !!!! it oscillates with 2.D0
        ELSE IF(Ctrl%GeOp%CoordType==CoordType_PrimInt) THEN
            HStre=GeOpCtrl%StreHessian
            HBend=GeOpCtrl%BendHessian
            HLinB=GeOpCtrl%LinBHessian
            HOutP=GeOpCtrl%OutPHessian
            HTors=GeOpCtrl%TorsHessian
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
!
        CALL UnSetDSYEVWork()
!
!project out redundancy
          CALL RedundancyOffFull(Displ%D,NCart)
        ELSE
          CALL Halt('Only Primitiv Internals are available yet.')
        ENDIF 
!
! Tidy up
!
        CALL Delete(Hessian)
!
      END SUBROUTINE DiagHessRFO
!
!-------------------------------------------------------
!
      SUBROUTINE GeOpConv(Ctrl,GMLoc,IntCs,IntOld)
        TYPE(SCFControls)         :: Ctrl
        TYPE(CRDS)                :: GMLoc
        TYPE(INTC)                :: IntCs
        TYPE(DBL_VECT)            :: IntOld,AuxVect
!
        REAL(DOUBLE)              :: MaxStreDispl,MaxBendDispl
        REAL(DOUBLE)              :: MaxLinBDispl,MaxOutPDispl
        REAL(DOUBLE)              :: MaxTorsDispl
        REAL(DOUBLE)              :: RMSIntDispl
!
        REAL(DOUBLE)              :: StreConvCrit,BendConvCrit
        REAL(DOUBLE)              :: LinBConvCrit,OutPConvCrit
        REAL(DOUBLE)              :: TorsConvCrit
!
        INTEGER                   :: I,J,K,L,ActStep     
        INTEGER                   :: NCart,NIntC,NatmsLoc
        REAL(DOUBLE)              :: RMSGrad,MaxGrad
        REAL(DOUBLE)              :: RMSGradNoConstr,MaxGradNoConstr
        REAL(DOUBLE)              :: Etot,Sum
!
        INTEGER                   :: NStreGeOp,NBendGeOp,NLinBGeOp
        INTEGER                   :: NOutPGeOp,NTorsGeOp
        INTEGER                   :: MaxStre,MaxBend,MaxLinB
        INTEGER                   :: MaxOutP,MaxTors
!
        INTEGER                   :: IMaxGrad,IMaxGradNoConstr
!
        IF(AllocQ(IntCs%Alloc)) THEN
          NIntC=SIZE(IntCs%Def)
        ELSE
          NIntC=0
        ENDIF
!
        NatmsLoc=GMLoc%Natms
        NCart=3*NatmsLoc
!
        RMSGrad=GeOpCtrl%RMSGrad
        RMSGradNoConstr=GeOpCtrl%RMSGradNoConstr
        MaxGrad=GeOpCtrl%MaxGrad
        MaxGradNoConstr=GeOpCtrl%MaxGradNoConstr
        ActStep=GeOpCtrl%ActStep
        IMaxGrad=GeOpCtrl%IMaxGrad
        IMaxGradNoConstr=GeOpCtrl%IMaxGradNoConstr
!
! For ionic crystals, until ionic radii is not coded
!
        IF(GeOpCtrl%CoordType==CoordType_Cartesian.AND.NIntC==0) THEN
          GeOpCtrl%GeOpConvgd=RMSGrad<GeOpCtrl%GradCrit.AND. &
                              MaxGrad<GeOpCtrl%GradCrit       
          CALL OpenASCII(OutFile,Out)
          WRITE(*,399) ActStep       
          WRITE(Out,399) ActStep       
          CALL Get(Etot,'ETot',StatsToChar(Ctrl%Current))
          WRITE(*,401) Etot
          WRITE(Out,401) Etot
          WRITE(*,900) RMSGrad,MaxGrad
          WRITE(Out,900) RMSGrad,MaxGrad
          CLOSE(Out,STATUS='KEEP')
900       FORMAT(' RMSGrad= ',F12.6,' MaxGrad= ',F12.6)
          RETURN
        ENDIF
!
        CALL Get(NStreGeOp,'NStreGeOp')
        CALL Get(NBendGeOp,'NBendGeOp')
        CALL Get(NLinBGeOp,'NLinBGeOp')
        CALL Get(NOutPGeOp,'NOutPGeOp')
        CALL Get(NTorsGeOp,'NTorsGeOp')
!
! Determine internal coordinate changes
!
        CALL INTCValue(IntCs,GMLoc%Carts%D)
        IntOld%D=IntCs%Value-IntOld%D
        CALL MapDisplTors(IntCs,NIntC,IntOld%D)
        MaxStre=0
        MaxBend=0
        MaxLinB=0
        MaxOutP=0
        MaxTors=0
        MaxStreDispl=Zero
        MaxBendDispl=Zero
        MaxLinBDispl=Zero
        MaxOutPDispl=Zero
        MaxTorsDispl=Zero
        DO I=1,NIntC 
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
        GeOpCtrl%MaxStreDispl=MaxStreDispl
        GeOpCtrl%MaxBendDispl=MaxBendDispl
        GeOpCtrl%MaxLinBDispl=MaxLinBDispl
        GeOpCtrl%MaxOutPDispl=MaxOutPDispl
        GeOpCtrl%MaxTorsDispl=MaxTorsDispl
        GeOpCtrl%RMSIntDispl=RMSIntDispl
!
        IF(GeOpCtrl%NConstr/=0) THEN
!  constraints
!          GeOpCtrl%GeOpConvgd=(RMSGradNoConstr<GeOpCtrl%GradCrit.AND. &
!                              MaxGradNoConstr<GeOpCtrl%GradCrit).OR. &
          GeOpCtrl%GeOpConvgd=(&
                              MaxStreDispl<GeOpCtrl%StreConvCrit.AND. &
                              MaxBendDispl<GeOpCtrl%BendConvCrit.AND. &
                              MaxLinBDispl<GeOpCtrl%LinBConvCrit.AND. &
                              MaxOutPDispl<GeOpCtrl%OutPConvCrit.AND. &
                              MaxTorsDispl<GeOpCtrl%TorsConvCrit)
        ELSE
! no constraints
          GeOpCtrl%GeOpConvgd=RMSGrad<GeOpCtrl%GradCrit.AND. &
                              MaxGrad<GeOpCtrl%GradCrit.AND. &
                              MaxStreDispl<GeOpCtrl%StreConvCrit.AND. &
                              MaxBendDispl<GeOpCtrl%BendConvCrit.AND. &
                              MaxLinBDispl<GeOpCtrl%LinBConvCrit.AND. &
                              MaxOutPDispl<GeOpCtrl%OutPConvCrit.AND. &
                              MaxTorsDispl<GeOpCtrl%TorsConvCrit
        ENDIF
!
! Review iterations
!
        CALL OpenASCII(OutFile,Out)
!
        WRITE(*,399) ActStep       
        WRITE(Out,399) ActStep       
        CALL Get(Etot,'ETot',StatsToChar(Ctrl%Current))
        WRITE(*,401) Etot
        WRITE(Out,401) Etot
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
        IF(GeOpCtrl%NConstr/=0) THEN
          WRITE(*,510) MaxGradNoConstr,IntCs%Atoms(IMaxGradNoConstr,1:4)
          WRITE(*,520) RMSGradNoConstr
          WRITE(Out,510) MaxGradNoConstr,IntCs%Atoms(IMaxGradNoConstr,1:4)
          WRITE(Out,520) RMSGradNoConstr
        ENDIF
!
        IF(NStreGeOp/=0) THEN
          WRITE(*,430) MaxStreDispl,IntCs%Atoms(MaxStre,1:2)
          WRITE(Out,430) MaxStreDispl,IntCs%Atoms(MaxStre,1:2)
        ENDIF
        IF(NBendGeOp/=0) THEN
          WRITE(*,435) MaxBendDispl,IntCs%Atoms(MaxBend,1:3)
          WRITE(Out,435) MaxBendDispl,IntCs%Atoms(MaxBend,1:3)
        ENDIF
        IF(NLinBGeOp/=0) THEN
          WRITE(*,436) MaxLinBDispl,IntCs%Atoms(MaxLinB,1:3)
          WRITE(Out,436) MaxLinBDispl,IntCs%Atoms(MaxLinB,1:3)
        ENDIF
        IF(NOutPGeOp/=0) THEN
          WRITE(*,437) MaxOutPDispl,IntCs%Atoms(MaxOutP,1:4)
          WRITE(Out,437) MaxOutPDispl,IntCs%Atoms(MaxOutP,1:4)
        ENDIF
        IF(NTorsGeOp/=0) THEN
          WRITE(*,438) MaxTorsDispl,IntCs%Atoms(MaxTors,1:4)
          WRITE(Out,438) MaxTorsDispl,IntCs%Atoms(MaxTors,1:4)
        ENDIF
!
        WRITE(*,440) RMSIntDispl
        WRITE(Out,440) RMSIntDispl
!
        CLOSE(Out,STATUS='KEEP')
!
399 FORMAT('GeOp step= ',I6)
400 FORMAT(' Total Energy at Current Geometry  = ',F20.8)
401 FORMAT('Total Energy at Previous Geometry  = ',F20.8)
410 FORMAT('           ','               Max Grad = ',F12.6,' between atoms ',4I4)
420 FORMAT('           ','               RMS Grad = ',F12.6)
510 FORMAT('  Max Grad on Unconstrained Coords = ',F12.6,' between atoms ',4I4)
520 FORMAT('  RMS Grad on Unconstrained Coords = ',F12.6)
430 FORMAT('           ','         Max STRE Displ = ',F12.6,' between atoms ',4I4)
435 FORMAT('           ','         Max BEND Displ = ',F12.6,' between atoms ',4I4)
436 FORMAT('           ','         Max LINB Displ = ',F12.6,' between atoms ',4I4)
437 FORMAT('           ','         Max OUTP Displ = ',F12.6,' between atoms ',4I4)
438 FORMAT('           ','         Max TORS Displ = ',F12.6,' between atoms ',4I4)
440 FORMAT('           ','              RMS Displ = ',F12.6,' between atoms ',4I4)
!
      END SUBROUTINE GeOpConv
!---------------------------------------------------------------
!
      SUBROUTINE SetGeOpCtrl(Ctrl,NatmsLoc)
!
        TYPE(SCFControls) :: Ctrl
        INTEGER :: AccL,NatmsLoc
        REAL(DOUBLE) :: Sum
!
        AccL=Ctrl%AccL(CBas)
!
        GeOpCtrl%AINVThrsh=1.D-8
        GeOpCtrl%BMatThrsh=1.D-8
!
        GeOpCtrl%AccL         =   AccL         
        GeOpCtrl%CoordType    =   Ctrl%GeOp%CoordType
        IF(GeOpCtrl%CoordType/=CoordType_Cartesian) THEN
          GeOpCtrl%DoInternals=.TRUE.
        ELSE
          GeOpCtrl%DoInternals=.FALSE.
        ENDIF
!
! Blocked sparse matrix algebra for GeOp
!
        GeOpCtrl%BlkGeomSize  =   3
!
! Hessian data, in AU. Data are based on OPLS force-field data.
!
!       GeOpCtrl%StreHessian   =   0.25D0   
!       GeOpCtrl%BendHessian   =   0.025D0
!       GeOpCtrl%LinBHessian   =   0.025D0
!       GeOpCtrl%OutPHessian   =   0.0025D0 
!       GeOpCtrl%TorsHessian   =   0.0025D0 
        GeOpCtrl%StreHessian   =   0.50D0   
        GeOpCtrl%BendHessian   =   0.20D0
        GeOpCtrl%LinBHessian   =   0.20D0
        GeOpCtrl%OutPHessian   =   0.05D0 
        GeOpCtrl%TorsHessian   =   0.05D0 
!
! Number of optimization steps
!
        GeOpCtrl%MaxGeOpSteps =   MAX(3*NatmsLoc,600)
!
! Convergence Criteria
!
        GeOpCtrl%GradCrit = GTol(AccL)
        GeOpCtrl%StreConvCrit= GeOpCtrl%GradCrit/GeOpCtrl%StreHessian
        GeOpCtrl%BendConvCrit= GeOpCtrl%GradCrit/GeOpCtrl%BendHessian
        GeOpCtrl%LinBConvCrit= GeOpCtrl%GradCrit/GeOpCtrl%LinBHessian
        GeOpCtrl%OutPConvCrit= GeOpCtrl%GradCrit/GeOpCtrl%OutPHessian
        GeOpCtrl%TorsConvCrit= GeOpCtrl%GradCrit/GeOpCtrl%TorsHessian
!
! Adaptive Reference GDIIS
!
        GeOpCtrl%NoGDIIS           =   Ctrl%GeOp%NoGDIIS
        GeOpCtrl%GDIISMetricOn     =   .FALSE.
        IF(GeOpCtrl%ActStep==1) THEN
          GeOpCtrl%GDIISMetric     =   1.D0 
        ENDIF
        GeOpCtrl%GDIISInit         =   2
        GeOpCtrl%GDIISMaxMem       =   500
        GeOpCtrl%GDIISBandWidth    =   0.50D0 !(in percent)
        GeOpCtrl%GDIISMinDomCount  =   1
!
! GDIIS restart control
!
        IF(GeOpCtrl%ActStep==1) THEN
          GeOpCtrl%GDIISOn=.TRUE.
        ELSE
!          IF(GeOpCtrl%RMSGrad>GeOpCtrl%OldRMSGrad) THEN
!!           GeOpCtrl%GDIISOn=.FALSE.
!            CALL Put(0,'SRMemory') 
!            CALL Put(0,'RefMemory') 
!          ENDIF
            GeOpCtrl%OldRMSGrad=GeOpCtrl%RMSGrad
        ENDIF
!
! Line Search
!
        GeOpCtrl%LSStepMax=8  
!
! grad. trf.
!
        GeOpCtrl%MaxIt_GrdTrf = 10 
        GeOpCtrl%GrdTrfCrit   = 2*GeOpCtrl%AINVThrsh
        GeOpCtrl%MaxGradDiff  = 5.D+2      
!
! iterative back trf.
!
        GeOpCtrl%MaxIt_CooTrf = 40
        GeOpCtrl%CooTrfCrit = MIN(GeOpCtrl%StreConvCrit, &
                                  GeOpCtrl%BendConvCrit, &
                                  GeOpCtrl%LinBConvCrit, &
                                  GeOpCtrl%OutPConvCrit, &
                                  GeOpCtrl%TorsConvCrit)/100.D0   
        GeOpCtrl%RMSCrit    = 0.90D0 !!! at least 10 percent decrease in RMS at a step
        GeOpCtrl%MaxCartDiff = 0.50D0  
        GeOpCtrl%DistRefresh = GeOpCtrl%MaxCartDiff*0.75D0
!
! Constraints
!
        CALL Get(GeOpCtrl%NConstr,'NConstraints')
        CALL Get(GeOpCtrl%NCartConstr,'NCartConstr')
        GeOpCtrl%ConstrMaxCrit = GeOpCtrl%CooTrfCrit*1.D-2
        GeOpCtrl%ConstrMax=GeOpCtrl%ConstrMaxCrit*10.D0
!
      END SUBROUTINE SetGeOpCtrl
!
!---------------------------------------------------------------
!
      SUBROUTINE VerifyStep(Ctrl,GMLoc,IntCs)
        TYPE(SCFControls)::   Ctrl
        TYPE(CRDS)    :: GMLoc
        TYPE(INTC)    :: IntCs
        INTEGER       :: I,J
        REAL(DOUBLE)  :: EtotCurr,EtotPrev
!
        IF(CGeo==1) RETURN
!
! Get recent total energies
!
        CALL Get(EtotCurr,'ETot',StatsToChar(Current))
        CALL Get(EtotPrev,'ETot',StatsToChar(Previous))
!
! Take care of iterative subspace
!
!!!     CALL GDIISSave(GMLoc%Natms,EtotCurr,EtotPrev)
!
! Carry out LineSearch, if EtotCurr>EtotPrev
!
        IF(EtotCurr>EtotPrev) THEN  
          CALL LineSearch(Ctrl,GMLoc,IntCs,EtotCurr,EtotPrev)
        ENDIF
!
      END SUBROUTINE VerifyStep
!
!---------------------------------------------------------------
!
      SUBROUTINE LineSearch(Ctrl,GMLoc,IntCs,EtotCurr,EtotPrev)
        TYPE(SCFControls)::   Ctrl
        TYPE(CRDS)       ::   GMLoc,GMPrev
        TYPE(INTC)       ::   IntCs
        INTEGER          ::   I,J,K,L,NDim,MaxStep,NStart,NCart,NIntC
        INTEGER          ::   OldGrad,OldGrad2,NatmsLoc
        REAL(DOUBLE)     ::   EtotCurr,EtotPrev,Etot,OldCrit
        REAL(DOUBLE)     ::   Fact,Delta
        TYPE(DBL_VECT)   ::   Displ
        LOGICAL          ::   DoInternals,DoGDIIS
        CHARACTER(LEN=DEFAULT_CHR_LEN) :: GMTag
        CHARACTER(LEN=DEFAULT_CHR_LEN) :: CoordTypeSave
!
! Search for a lower energy point between CGeo-1 and CGeo
! by compressing the range of the GDIIS matrix
!
          GMTag=''
#ifdef MMech
        IF(HasMM()) GMTag='GM_MM'
#endif
!
        CALL Get(GMPrev,TRIM(GMTag)//PrvGeom)
!
        NatmsLoc=GMLoc%Natms
        NCart=3*NatmsLoc
        NIntC=SIZE(IntCs%Def)
        MaxStep=GeOpCtrl%LSStepMax
!
        DoInternals=GeOpCtrl%DoInternals
!
! Calculate steepest descent step
!
        OldGrad=Ctrl%Grad
        Ctrl%Grad=GRAD_STPDESC_OPT
        CoordTypeSave=GeOpCtrl%CoordType
        GeOpCtrl%CoordType=CoordType_Cartesian
!
        CALL NewDispl(Ctrl,Displ,NCart,NIntC)
        CALL SRStep(Ctrl,GMLoc,Displ,IntCs,&
                    Step_O=PrvGeom,DoSave_O=.FALSE.)
!
! Do line search
!
          IF(OldGrad==GRAD_STPDESC_OPT) THEN
            Fact=0.5D0
          ELSE
            Fact=One
          ENDIF
!
        DO I=1,MaxStep
!
          Displ%D=Fact*Displ%D
          DoGDIIS=GeOpCtrl%GDIISOn
          GeOpCtrl%GDIISOn=.FALSE.
          OldCrit=GeOpCtrl%CooTrfCrit
          GeOpCtrl%CooTrfCrit=1.D-8
            CALL NewStructure(Ctrl,GMLoc,Displ,IntCs,DoSave_O=.FALSE.)
          GeOpCtrl%CooTrfCrit=OldCrit
          GeOpCtrl%GDIISOn=DoGDIIS
!
          CALL Put(GMLoc,Tag_O=TRIM(GMTag)//CurGeom)
          OldGrad2=Ctrl%Grad
          Ctrl%Grad=GRAD_NO_GRAD
            CALL CALC_EnergyOnly(Ctrl)
          Ctrl%Grad=OldGrad2    
          CALL Get(Etot,'ETot',StatsToChar(Current))
!
          CALL OpenAscii(OutFile,Out) 
            WRITE(Out,100) I,ETot  
            IF(PrintFlags%Geop==DEBUG_GEOP) THEN
              WRITE(*,100) I,ETot  
            ENDIF
100  FORMAT('LineSearch Step= ',I2,' Energy= ',F12.8)
          Close(Out,STATUS='KEEP')
!
          IF(Etot<EtotPrev) THEN
            EXIT
          ELSE
            Fact=Fact/Two
          ENDIF
        ENDDO
!
        CALL OpenAscii(OutFile,Out) 
          WRITE(Out,200) I            
          WRITE(*,200) I            
        Close(Out,STATUS='KEEP')
200     FORMAT('WARNING! Line Search did not converge in ',I2,' steps')
!
! Set back control
!
        Ctrl%Grad=OldGrad
        GeOpCtrl%CoordType=CoordTypeSave
!
        CALL Delete(Displ)
!
      END SUBROUTINE LineSearch
!
!---------------------------------------------------------------
!
END MODULE GeomOpt
 
