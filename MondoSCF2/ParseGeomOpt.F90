MODULE ParseGeomOpt
   USE DerivedTypes
   USE GlobalScalars
   USE GeomOptKeys
   USE Macros
   USE ControlStructures
   USE OptionKeys
   IMPLICIT NONE
   !
   CONTAINS
   !
   SUBROUTINE LoadGeomOpt(N,GOpt,PBCDim)
     TYPE(FileNames)                :: N
     TYPE(GeomOpt)                  :: GOpt  
     CHARACTER(LEN=DCL)             :: Max_Steps = 'Max_Steps'
     INTEGER                        :: PBCDim
     !
     CALL OpenASCII(N%IFile,Inp)
     ! maximum number of optimization steps
     IF(.NOT.OptIntQ(Inp,Max_Steps,GOpt%GConvCrit%MaxGeOpSteps))THEN
        GOpt%GConvCrit%MaxGeOpSteps=500  
     ENDIF
     !
     ! Optimizer type
     !
     IF(OptKeyQ(Inp,GRADIENTS,OPT_StpDesc))THEN
       GOpt%Optimizer=GRAD_StpDesc_OPT
     ELSEIF(OptKeyQ(Inp,GRADIENTS,OPT_DiagHess))THEN
       GOpt%Optimizer=GRAD_DiagHess_OPT
     ELSEIF(OptKeyQ(Inp,GRADIENTS,OPT_BiSect))THEN
       GOpt%Optimizer=GRAD_BiSect_OPT
     ELSE
       GOpt%Optimizer=GRAD_DiagHess_OPT !default
     ENDIF
     !
     ! Parse for GDIIS 
     !
     GOpt%GDIIS%NoGDIIS=.FALSE.
     IF(OptKeyQ(Inp,GRADIENTS,OPT_NoGDIIS)) THEN
       GOpt%GDIIS%NoGDIIS=.TRUE.
     ENDIF
     !
     ! Parse for minimization type: gradient or gradient norm
     !
     GOpt%DoGradNorm=.FALSE.
     IF(OptKeyQ(Inp,GRADIENTS,OPT_GradNorm)) THEN
       GOpt%DoGradNorm=.TRUE.
     ENDIF
     !
     ! Parse for printing pictures of fits
     !
     GOpt%Pictures=.FALSE.
     IF(OptKeyQ(Inp,GRADIENTS,OPT_Pictures)) THEN
       GOpt%Pictures=.TRUE.
     ENDIF
     !
     ! Parse for addinging explicit lattice coordinates to the optimization
     !
     GOpt%GConvCrit%ExplLatt=.FALSE.
     IF(OptKeyQ(Inp,GRADIENTS,OPT_ExplLatt)) THEN
       GOpt%GConvCrit%ExplLatt=.TRUE.
     ENDIF
     !
     ! Parse for alternating lattice and atomic positions relaxation 
     !
     GOpt%GConvCrit%Alternate=.FALSE.
     GOpt%GConvCrit%FixedAtomsFirst=.FALSE.
     IF(OptKeyQ(Inp,GRADIENTS,OPT_Alternate)) THEN
       GOpt%GConvCrit%Alternate=.TRUE.
       IF(PBCDim==0) GOpt%GConvCrit%Alternate=.FALSE.
       IF(OptKeyQ(Inp,GRADIENTS,OPT_FixedAtomsFirst)) THEN
         GOpt%GConvCrit%FixedAtomsFirst=.TRUE.
       ENDIF
     ENDIF
     !
     ! Parse for energy-back-tracking
     !
     GOpt%GConvCrit%DoAtomBackTr=.TRUE.
     IF(OptKeyQ(Inp,GRADIENTS,OPT_DoAtomBackTr)) THEN
       GOpt%GConvCrit%DoAtomBacktr=.TRUE.
     ENDIF
     GOpt%GConvCrit%DoLattBackTr=.TRUE.
     IF(OptKeyQ(Inp,GRADIENTS,OPT_DoLattBackTr)) THEN
       GOpt%GConvCrit%DoLattBacktr=.TRUE.
     ENDIF
     GOpt%GConvCrit%NoBackTr=.FALSE.
     IF(OptKeyQ(Inp,GRADIENTS,OPT_NoBackTr)) THEN
       GOpt%GConvCrit%NoBackTr=.TRUE.
       GOpt%GConvCrit%DoLattBacktr=.FALSE.
       GOpt%GConvCrit%DoAtomBacktr=.FALSE.
     ENDIF
     !
     ! Parse for projecting out rotations and translations
     ! from geometry displacements.
     !
     GOpt%TrfCtrl%DoRotOff=.TRUE.
     GOpt%TrfCtrl%DoTranslOff=.TRUE.
     IF(OptKeyQ(Inp,GRADIENTS,OPT_NoRotOff)) THEN
        GOpt%TrfCtrl%DoRotOff=.FALSE.
     ENDIF
     IF(OptKeyQ(Inp,GRADIENTS,OPT_NoTranslOff)) THEN
        GOpt%TrfCtrl%DoTranslOff=.FALSE.
     ENDIF
     !
     ! Parse for VDW Radius factor
     !
     IF(.NOT.OptDblQ(Inp,VDWFACT,GOpt%CoordCtrl%VDWFact)) THEN
       GOpt%CoordCtrl%VDWFact=0.8D0 !default value
     ELSE
     ! GOpt%CoordCtrl%VDWFact=GOpt%CoordCtrl%VDWFact*One !default value
     ! GOpt%CoordCtrl%VDWFact=GOpt%CoordCtrl%VDWFact*0.8 !default value
     ENDIF
     !
     ! Parse for MaxAtoms and MaxLattice for Alternating optimization
     !
     IF(.NOT.OptIntQ(Inp,MaxAtoms,GOpt%GConvCrit%MaxAtoms)) THEN
       GOpt%GConvCrit%MaxAtoms=10000 !default value
     ENDIF
     IF(.NOT.OptIntQ(Inp,MaxLattice,GOpt%GConvCrit%MaxLattice)) THEN
       GOpt%GConvCrit%MaxLattice=1 !default value
     ENDIF
     !
     !
     ! Parse for Maximum angle and maximum bondlength displacements
     !
     IF(.NOT.OptDblQ(Inp,MaxAngle,GOpt%CoordCtrl%MaxAngle)) THEN
      !GOpt%CoordCtrl%MaxAngle=Two*PI !default value
       GOpt%CoordCtrl%MaxAngle=10.D0*PI/180.D0  !default value
      !IF(PBCDim==3) THEN
      !  GOpt%CoordCtrl%MaxAngle=5.D0*PI/180.D0
      !ELSE
      !  GOpt%CoordCtrl%MaxAngle=15.D0*PI/180.D0
      !ENDIF
     ELSE
       GOpt%CoordCtrl%MaxAngle=GOpt%CoordCtrl%MaxAngle*PI/180.D0
     ENDIF
     IF(.NOT.OptDblQ(Inp,MaxStre,GOpt%CoordCtrl%MaxStre)) THEN
       GOpt%CoordCtrl%MaxStre=100.000D0 !default value
      !GOpt%CoordCtrl%MaxStre=0.15D0*AngstromsToAu
     ELSE
       GOpt%CoordCtrl%MaxStre=GOpt%CoordCtrl%MaxStre*AngstromsToAu
     ENDIF
     !
     ! Parse for Internal coordinates refresh type
     !
     IF(.NOT.OptIntQ(Inp,INTCREFRESH,GOpt%CoordCtrl%RefreshIn)) THEN
       GOpt%CoordCtrl%RefreshIn=1 !default value
     ENDIF
     !
     ! Parse for coordtype
     !
     IF(OptKeyQ(Inp,GRADIENTS,CoordType_PrimInt)) THEN
        GOpt%CoordCtrl%CoordType=CoordType_PrimInt
     ELSE IF(OptKeyQ(Inp,GRADIENTS,CoordType_Cartesian)) THEN
        GOpt%CoordCtrl%CoordType=CoordType_Cartesian
     ELSE
        GOpt%CoordCtrl%CoordType=CoordType_PrimInt
     ENDIF
     !
     ! Parse for Steepest descent Inverse Hessian
     !
     IF(.NOT.OptDblQ(Inp,STPDESCINVH,GOpt%Hessian%StpDescInvH)) THEN
       GOpt%Hessian%StpDescInvH=0.2 !default value
       IF(GOpt%CoordCtrl%CoordType==CoordType_Cartesian) &
          GOpt%Hessian%StpDescInvH=0.1D0
     ENDIF
     !
     ! Parse for Cholesky-fact.
     !
     GOpt%TrfCtrl%DoNewChol=.FALSE.
     IF(OptKeyQ(Inp,GRADIENTS,CoordType_DoNewChol)) THEN
        GOpt%TrfCtrl%DoNewChol=.TRUE.
     ENDIF
     !
     ! Do classical coordinate transformation? (No cartesian internals)
     !
     GOpt%TrfCtrl%DoClssTrf=.TRUE.
     !GOpt%TrfCtrl%DoClssTrf=.FALSE.
     IF(OptKeyQ(Inp,GRADIENTS,CoordType_DoClssTrf)) THEN
       GOpt%TrfCtrl%DoClssTrf=.TRUE.
     ENDIF
     !
     ! Fix MM coordinates?
     !
     GOpt%Constr%DoFixMM=.FALSE.
     IF(OptKeyQ(Inp,GRADIENTS,CoordType_DoFixMM)) THEN
       GOpt%Constr%DoFixMM=.TRUE.
     ENDIF
     !
     CLOSE(Inp,STATUS='KEEP')
     !
   END SUBROUTINE LoadGeomOpt
  
END MODULE ParseGeomOpt
