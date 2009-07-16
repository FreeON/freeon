!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
MODULE ParseGeomOpt
   USE DerivedTypes
   USE GlobalScalars
   USE GeomOptKeys
   USE Macros
   USE ControlStructures
   USE OptionKeys
   USE MondoLogger

   IMPLICIT NONE

   CONTAINS

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
     ! Parse for printing back-transformation of fits
     !
     GOpt%TrfCtrl%PrtbackTr=.FALSE.
     IF(OptKeyQ(Inp,GRADIENTS,OPT_PrtbackTr)) THEN
       GOpt%TrfCtrl%PrtbackTr=.TRUE.
     ENDIF
     !
     ! Parse for printing back-transformation of fits
     !
     GOpt%TrfCtrl%NoBTRep=.FALSE.
     IF(OptKeyQ(Inp,GRADIENTS,OPT_NoBTRep)) THEN
       GOpt%TrfCtrl%NoBTRep=.TRUE.
     ENDIF
     !
     ! Parse for addinging explicit lattice coordinates to the optimization
     !
     GOpt%GConvCrit%ExplLatt=.FALSE.
     IF(OptKeyQ(Inp,GRADIENTS,OPT_ExplLatt)) THEN
       GOpt%GConvCrit%ExplLatt=.TRUE.
     ENDIF
     !
     ! Parse for VDW-bonding related options
     !
     GOpt%GConvCrit%NonCovBend=.FALSE.
     IF(OptKeyQ(Inp,GRADIENTS,OPT_NonCovBend)) THEN
       GOpt%GConvCrit%NonCovBend=.TRUE.
     ENDIF
     GOpt%GConvCrit%NonCovTors=.FALSE.
     IF(OptKeyQ(Inp,GRADIENTS,OPT_NonCovTors)) THEN
       GOpt%GConvCrit%NonCovTors=.TRUE.
     ENDIF
     GOpt%GConvCrit%HBondOnly=.FALSE.
     IF(OptKeyQ(Inp,GRADIENTS,OPT_HBondOnly)) THEN
       GOpt%GConvCrit%HBondOnly=.TRUE.
     ENDIF
     GOpt%GConvCrit%NoFragmConnect=.FALSE.
     IF(OptKeyQ(Inp,GRADIENTS,OPT_NoFragmConnect)) THEN
       GOpt%GConvCrit%NoFragmConnect=.TRUE.
     ENDIF
     !
     ! Parse for lattice parameter ratios
     !
     CALL FindLattRattio(OPT_RatioABC,GOpt%Constr%RatioABC)
     CALL FindLattRattio(OPT_RatioAlpBetGam,GOpt%Constr%RatioAlpBetGam)
     !
     ! Parse for alternating lattice and atomic positions relaxation
     !
     GOpt%GConvCrit%Alternate=.FALSE.
     GOpt%GConvCrit%LatticeStart=.FALSE.
     IF(OptKeyQ(Inp,GRADIENTS,OPT_Alternate)) THEN
       GOpt%GConvCrit%Alternate=.TRUE.
       IF(PBCDim==0) GOpt%GConvCrit%Alternate=.FALSE.
       IF(OptKeyQ(Inp,GRADIENTS,OPT_LatticeStart)) THEN
         GOpt%GConvCrit%LatticeStart=.TRUE.
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
       GOpt%CoordCtrl%VDWFact=1.0D0 !default value
     ELSE
     ! GOpt%CoordCtrl%VDWFact=GOpt%CoordCtrl%VDWFact*One !default value
     ! GOpt%CoordCtrl%VDWFact=GOpt%CoordCtrl%VDWFact*0.8 !default value
     ENDIF
     !
     ! Parse for MaxAtoms and MaxLattice for Alternating optimization
     !
     IF(.NOT.OptIntQ(Inp,MaxAtomSteps,GOpt%GConvCrit%MaxAtomSteps)) THEN
       GOpt%GConvCrit%MaxAtomSteps=10000 !default value
     ENDIF
     IF(.NOT.OptIntQ(Inp,MaxLatticeSteps,GOpt%GConvCrit%MaxLatticeSteps)) THEN
       GOpt%GConvCrit%MaxLatticeSteps=1 !default value
     ENDIF
     !
     ! Parse for filtering data points in fitting of QUICCA
     !
     GOpt%CoordCtrl%DoQFilter=.FALSE.
     IF(OptKeyQ(Inp,GRADIENTS,OPT_DoQFilter)) THEN
       GOpt%CoordCtrl%DoQFilter=.TRUE.
     ENDIF
     !
     ! Parse for Maximum angle and maximum bondlength displacements
     !
     IF(.NOT.OptDblQ(Inp,MaxAngle,GOpt%CoordCtrl%MaxAngle)) THEN
       GOpt%CoordCtrl%MaxAngle=5.D0*DegToRad  !default value
      !IF(PBCDim==3) THEN
      !  GOpt%CoordCtrl%MaxAngle=5.D0*DegToRad
      !ELSE
      !  GOpt%CoordCtrl%MaxAngle=15.D0*DegToRad
      !ENDIF
     ELSE
       GOpt%CoordCtrl%MaxAngle=GOpt%CoordCtrl%MaxAngle*DegToRad
     ENDIF
     CALL MondoLog(DEBUG_NONE, "FreeON", "using MaxAngle = "//TRIM(FltToChar(GOpt%CoordCtrl%MaxAngle*RadToDeg))//" degrees")

     IF(.NOT.OptDblQ(Inp,MaxStre,GOpt%CoordCtrl%MaxStre)) THEN
       GOpt%CoordCtrl%MaxStre=0.3D0*AngstromsToAU
     ELSE
       GOpt%CoordCtrl%MaxStre=GOpt%CoordCtrl%MaxStre*AngstromsToAU
     ENDIF
     CALL MondoLog(DEBUG_NONE, "FreeON", "using MaxStre = "//TRIM(FltToChar(GOpt%CoordCtrl%MaxStre*AUToAngstroms))//" Angstrom")

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
   !
   !------------------------------------------------------------------
   !
   SUBROUTINE FindLattRattio(QChar,Ratio)
     CHARACTER(LEN=DCL)             :: Line,Aux
     CHARACTER(LEN=*)               :: QChar
     INTEGER                        :: J,ChLen
     REAL(DOUBLE),DIMENSION(3)      :: Ratio
     !
     Ratio=-One
     IF(FindMixedCaseKey(QChar,Inp)) THEN
       ChLen=LEN(QChar)
       DO J=1,LEN(Line) ; Line(J:J)=' ' ; ENDDO
       Line=OPTIONS_BEGIN
       CALL LowCase(Line)
       CALL AlignLowCase(TRIM(Line),Inp)
       DO
         READ(Inp,DEFAULT_CHR_FMT,END=1) Line
         CALL RemoveComments(Line)
         IF(INDEX(Line,QChar)/=0) THEN
           READ(Line,*) Aux(1:ChLen),(Ratio(J),J=1,3)
           EXIT
         ENDIF
         IF(INDEX(Line,OPTIONS_END)/=0) EXIT
       ENDDO
       1 CONTINUE
     ENDIF
   END SUBROUTINE FindLattRattio
   !
END MODULE ParseGeomOpt
