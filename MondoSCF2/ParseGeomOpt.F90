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
   SUBROUTINE LoadGeomOpt(N,GOpt)
     TYPE(FileNames)                :: N
     TYPE(GeomOpt)                  :: GOpt  
     CHARACTER(LEN=DCL)             :: Max_Steps = 'Max_Steps'
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
     ! Parse for energy-back-tracking
     !
     GOpt%GConvCrit%DoBackTr=.True.
     IF(OptKeyQ(Inp,GRADIENTS,OPT_NoBackTr)) THEN
       GOpt%GConvCrit%DoBacktr=.FALSE.
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
       GOpt%CoordCtrl%VDWFact=One !default value
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
     !IF(OptKeyQ(Inp,COORDTYPE,CoordType_DoClssTrf)) THEN
     !  GOpt%TrfCtrl%DoClssTrf=.TRUE.
     !ENDIF
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
