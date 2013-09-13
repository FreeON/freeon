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

MODULE QUICCAMod

   USE DerivedTypes
   USE GlobalScalars
   USE GlobalObjects
   USE GlobalCharacters
   USE InOut
   USE MemMan
   USE LinAlg
   USE AInv
   USE CholFactor
   USE ControlStructures
   USE PunchHDF
   USE SetXYZ
   USE InCoords
   USE HessianMod
   USE Slatec

IMPLICIT NONE

CONTAINS
!
!--------------------------------------------------------------
!
   SUBROUTINE ModeFit(IntCs,IntCGrads,IntCValues,RefPoints,PBCDim,GOpt, &
                      PredVals,PWDPath,SCRPath,XYZ,NCart,iGEO,Print)
     TYPE(INTC)                  :: IntCs,IntCsT
     REAL(DOUBLE),DIMENSION(:,:) :: IntCGrads,IntCValues,XYZ
     REAL(DOUBLE),DIMENSION(:)   :: PredVals,RefPoints
     CHARACTER(LEN=*)            :: PWDPath,SCRPath
     CHARACTER(LEN=DCL)          :: Path2
     INTEGER                     :: iGEO,NMem,NDim,I,J,Print
     INTEGER                     :: NCart,NatmsLoc,PBCDim
     TYPE(DBL_VECT)              :: LastIntCVals,LastIntcGrads
     TYPE(DBL_RNK2)              :: DXVects,DXValues,DXGrads,DelocVals
     TYPE(DBL_RNK2)              :: RangeT,ABCT,LWeightT
     TYPE(INT_VECT)              :: NDegsT
     REAL(DOUBLE)                :: X
     TYPE(GeomOpt)               :: GOpt
     !
     NMem=SIZE(IntCGrads,2)
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     CALL New(LastIntCVals,IntCs%N)
     CALL New(LastIntCGrads,IntCs%N)
     CALL INTCValue(IntCs,XYZ,PBCDim,GOpt%CoordCtrl%LinCrit, &
                    GOpt%CoordCtrl%TorsLinCrit)
     LastIntCVals%D=IntCs%Value%D
     CALL SetBackToRefs(LastIntCVals%D,IntCs,RefPoints)
     LastIntCGrads%D=Zero
    !LastIntCVals%D=IntCValues(:,NMem)
    !LastIntCGrads%D=IntCGrads(:,NMem)
     !
     CALL New(DelocVals,(/NCart,NMem/))
     CALL CollectDelocs(DelocVals%D,IntCValues,LastIntCVals%D,IntCs, &
                        XYZ,GOpt%TrfCtrl,GOpt%CoordCtrl,GOpt%GConvCrit,&
                        SCRPath,Print,PBCDim)
     CALL CollectDXVects(IntCValues,LastIntCVals%D, &
                         DXVects,Delocs_O=DelocVals%D)
    !CALL CollectDXVects(IntCValues,LastIntCVals%D,DXVects)
     CALL Delete(DelocVals)
     !
     NDim=SIZE(DXVects%D,2)
     CALL CollectDXProjection(XYZ,IntCs,GOpt,SCRPath,Print, &
                              IntCValues,IntCGrads,LastIntCVals%D,&
                             LastIntCGrads%D,DXVects%D,DXGrads,DXValues)
     CALL New(LWeightT,(/NDim,NMem/))
     !
     DO I=1,NMem
       X=Zero
       DO J=1,IntCs%N
         X=X+IntCGrads(J,I)**2
       ENDDO
       DO J=1,NDim ; LWeightT%D(J,I)=X ; ENDDO
     ENDDO
     CALL New(NDegsT,NDim)
     CALL New(IntCsT,NDim)
     IntCsT%Active%L=.TRUE.
     CALL New(RangeT,(/NDim,2/))
     CALL New(ABCT,(/NDim,4/))
     ABCT%D(:,:)=Zero
     Path2=TRIM(PWDPath)//'_M_'//TRIM(IntToChar(iGEO))
     !
     CALL LQFit(DXValues%D,DXGrads%D,LWeightT%D,IntCsT,ABCT%D, &
                RangeT%D,NDegsT%I,Zero,.TRUE.)
     CALL DoPredict(ABCT%D,DXValues%D,DXGrads%D,IntCsT, &
                    NDegsT%I,Path2,RangeT%D,GOpt%DoGradNorm)
     CALL CtlrCartRange(IntCsT,RangeT%D)
    !CALL PrtFitM(DXValues%D,DXGrads%D,ABCT%D,IntCsT,Path2)
     !
     ! Now, compose the new geometry
     !
     DO J=1,IntCs%N ; PredVals(J)=LastIntCVals%D(J) ; ENDDO
     DO I=1,NDim
      !IF(RangeT%D(I,1)<IntCsT%PredVal%D(I).AND. &
      !   RangeT%D(I,2)>IntCsT%PredVal%D(I)) THEN
         DO J=1,IntCs%N
           PredVals(J)=PredVals(J)+IntCsT%PredVal%D(I)*DXVects%D(J,I)
         ENDDO
      !ENDIF
     ENDDO
     !
     CALL Delete(IntCsT)
     CALL Delete(RangeT)
     CALL Delete(ABCT)
     CALL Delete(NDegsT)
     CALL Delete(LWeightT)
     CALL Delete(DXValues)
     CALL Delete(DXGrads)
     CALL Delete(DXVects)
     CALL Delete(LastIntCVals)
     CALL Delete(LastIntCGrads)
   END SUBROUTINE ModeFit
!
!-------------------------------------------------------------------
!
   SUBROUTINE CollectDelocs(DelocVals,IntCValues,LastIntCVals,IntCs, &
                            XYZ,GTrfCtrl,GCoordCtrl,GConvCr,SCRPath,&
                            Print,PBCDim)
     TYPE(INTC)                  :: IntCs
     REAL(DOUBLE),DIMENSION(:,:) :: DelocVals,XYZ,IntCValues
     REAL(DOUBLE),DIMENSION(:)   :: LastIntCVals
     CHARACTER(LEN=*)            :: SCRPath
     TYPE(CoordCtrl)             :: GCoordCtrl
     TYPE(TrfCtrl)               :: GTrfCtrl
     INTEGER                     :: I,J,NMem,NCart,Print,PBCDim
     TYPE(INT_VECT)              :: ISpB,JSpB
     TYPE(DBL_VECT)              :: ASpB,CartVect,IntCDispl
     TYPE(Cholesky)              :: CholData
     TYPE(GConvCrit)             :: GConvCr
     !
     CALL RefreshBMatInfo(IntCs,XYZ,GTrfCtrl,GConvCr, &
                          GCoordCtrl%LinCrit,&
                          GCoordCtrl%TorsLinCrit,PBCDim, &
                          Print,SCRPath)
     CALL GetBMatInfo(SCRPath,ISpB,JSpB,ASpB,CholData)
     !
     NMem=SIZE(DelocVals,2)
     NCart=SIZE(DelocVals,1)
     CALL New(IntCDispl,IntCs%N)
     CALL New(CartVect,NCart)
     !
     DO I=1,NMem
       DO J=1,IntCs%N
         IntCDispl%D(J)=IntCValues(J,I)-LastIntCVals(J)
       ENDDO
       CALL CALC_BxVect(ISpB,JSpB,ASpB,IntCDispl%D,CartVect%D,Trp_O=.TRUE.)
       CALL GcInvIter(CartVect%D,ISpB,JSpB,ASpB,CholData,IntCs%N)
       DO J=1,NCart ; DelocVals(J,I)=CartVect%D(J) ; ENDDO
     ENDDO
     !
     CALL Delete(CartVect)
     CALL Delete(IntCDispl)
     CALL Delete(ISpB)
     CALL Delete(JSpB)
     CALL Delete(ASpB)
     CALL Delete(CholData)
   END SUBROUTINE CollectDelocs
!
!-------------------------------------------------------------------
!
   SUBROUTINE CtlrCartRange(IntCsT,RangeT)
     TYPE(INTC)                  :: IntCsT
     REAL(DOUBLE),DIMENSION(:,:) :: RangeT
     INTEGER                     :: I,J
     REAL(DOUBLE)                :: DX,DXP
     !
     DO I=1,IntCsT%N
       IF(HasLattice(IntCsT%Def%C(I))) CYCLE
       DX=RangeT(I,2)-RangeT(I,1)
       DXP=IntCsT%PredVal%D(I)
       IF(DXP>RangeT(I,2)+DX) THEN
         IntCsT%PredVal%D(I)=RangeT(I,2)+DX
       ELSE IF(DXP<RangeT(I,1)-DX) THEN
         IntCsT%PredVal%D(I)=RangeT(I,1)-DX
       ENDIF
     ENDDO
   END SUBROUTINE CtlrCartRange
!
!-------------------------------------------------------------------
!
   SUBROUTINE CollectDXProjection(XYZ,IntCsX,GOpt,SCRPath,Print, &
                                  RefStruct,RefGrad,LastStruct, &
                                  LastGrad,DXVects,DXGrads,DXValues)
     REAL(DOUBLE),DIMENSION(:,:) :: RefStruct,RefGrad,DXVects,XYZ
     REAL(DOUBLE),DIMENSION(:)   :: LastStruct,LastGrad
     TYPE(DBL_RNK2)              :: DXValues,DXGrads
     INTEGER                     :: Print
     INTEGER                     :: II,I,J,NCart,NMem,NDim
     REAL(DOUBLE)                :: G,X
     TYPE(INT_VECT)              :: ISpB,JSpB
     TYPE(DBL_VECT)              :: ASpB
     TYPE(Cholesky)              :: CholData
     TYPE(INTC)                  :: IntCsX
     TYPE(GeomOpt)               :: GOpt
     CHARACTER(LEN=*)            :: SCRPath
     !
     NCart=SIZE(RefStruct,1)
     NMem=SIZE(RefStruct,2)
     NDim=SIZE(DXVects,2)
     CALL New(DXGrads,(/NDim,NMem/))
     CALL New(DXValues,(/NDim,NMem/))
     !
!!! this routines should be used calculate proper overlaps
  !  CALL RefreshBMatInfo(IntCsX,XYZ,GOpt%TrfCtrl,GOpt%GConvCrit, &
  !                 GOpt%CoordCtrl%LinCrit,GOpt%CoordCtrl%TorsLinCrit, &
  !                       PBCDim,Print,SCRPath,Gi_O=.TRUE.)
  !  CALL GetBMatInfo(SCRPath,ISpB,JSpB,ASpB,CholData)
  !  CALL GiInvIter(ISpB,JSpB,ASpB,CholData,IntA1,IntA2, &
  !                 NCart,IntCsX%N)
     !
     ! Project Grads and distance of most recent structure to each previous
     ! structures onto the selected axis system
     !
     DO II=1,NMem
       DO I=1,NDim
         G=Zero
         X=Zero
         DO J=1,NCart
           G=G+DXVects(J,I)*(RefGrad(J,II)-LastGrad(J))
           X=X+DXVects(J,I)*(RefStruct(J,II)-LastStruct(J))
         ENDDO
         DXGrads%D(I,II)=G
         DXValues%D(I,II)=X
       ENDDO
     ENDDO
   END SUBROUTINE CollectDXProjection
!
!-------------------------------------------------------------------
!
   SUBROUTINE CollectDXVects(RefStruct,LastStruct,DXVects,Delocs_O)
     REAL(DOUBLE),DIMENSION(:,:) :: RefStruct
     REAL(DOUBLE),DIMENSION(:)   :: LastStruct
     TYPE(DBL_RNK2)              :: DXVects
     INTEGER                     :: I,J,NMem,NCart,INFO,NDim,NCart2
     TYPE(DBL_RNK2)              :: ScaledD,Overlap,UMat
     TYPE(DBL_VECT)              :: CartVect,ScaleFact
     REAL(DOUBLE)                :: Fact,Crit
     REAL(DOUBLE),DIMENSION(:,:),OPTIONAL :: Delocs_O
     !
     NMem=SIZE(RefStruct,2)
     NCart=SIZE(RefStruct,1)
     Crit=1.D-7
     CALL New(CartVect,NCart)
     CALL New(ScaledD,(/NCart,NMem/))
     CALL New(Overlap,(/NMem,NMem/))
     CALL New(ScaleFact,NMem)
     !
     DO I=1,NMem
       DO J=1,NCart; CartVect%D(J)=RefStruct(J,I)-LastStruct(J) ;ENDDO
       DO J=1,NCart ; ScaledD%D(J,I)=CartVect%D(J) ; ENDDO
     ENDDO
     !
     IF(PRESENT(Delocs_O)) THEN
       NCart2=SIZE(Delocs_O,1)
       CALL DGEMM_TNc(NMem,NCart2,NMem,One,Zero,Delocs_O,Delocs_O,Overlap%D)
     ELSE
       CALL DGEMM_TNc(NMem,NCart,NMem,One,Zero,ScaledD%D,ScaledD%D,Overlap%D)
     ENDIF
     !
     DO I=1,NMem
       ScaleFact%D(I)=One/(Overlap%D(I,I)+1.D-10)
     ENDDO
     DO I=1,NMem
       DO J=1,NMem
         Overlap%D(I,J)=Overlap%D(I,J)*ScaleFact%D(I)*ScaleFact%D(J)
       ENDDO
     ENDDO
     !
     CALL SetDSYEVWork(NMem)
       BLKVECT%D=Overlap%D
       CALL DSYEV('V','U',NMem,BLKVECT%D,BIGBLOK,BLKVALS%D, &
       BLKWORK%D,BLKLWORK,INFO)
       IF(INFO/=SUCCEED) &
       CALL Halt('DSYEV hosed in CollectDXVects. INFO='&
                  //TRIM(IntToChar(INFO)))
       NDim=0
       DO I=1,NMem
         IF(ABS(BLKVALS%D(I))>Crit) NDim=NDim+1
       ENDDO
       !
       CALL New(UMat,(/NMem,NDim/))
       CALL New(DXVects,(/NCart,NDim/))
       NDim=0
       DO I=1,NMem
         IF(ABS(BLKVALS%D(I))>Crit) THEN
           NDim=NDim+1
           DO J=1,NMem
             UMat%D(J,NDim)=BLKVECT%D(J,I)/SQRT(BLKVALS%D(I))*ScaleFact%D(J)
           ENDDO
         ENDIF
       ENDDO
     !
     ! Now, generate orthogonal directions
     !
     CALL DGEMM_NNc(NCart,NMem,NDim,One,Zero,ScaledD%D,Umat%D,DXVects%D)
     !
     CALL UnSetDSYEVWork()
     CALL Delete(Overlap)
     CALL Delete(UMat)
     CALL Delete(ScaledD)
     CALL Delete(CartVect)
     CALL Delete(ScaleFact)
   END SUBROUTINE CollectDXVects
!
!-------------------------------------------------------------------
!
   SUBROUTINE CurviLinOverlap(ScaledD,XYZ,Overlap,IntCs,SCRPath)
     REAL(DOUBLE),DIMENSION(:,:) :: ScaledD,XYZ,Overlap
     TYPE(INTC)                  :: IntCs
     CHARACTER(LEN=*)            :: SCRPath
     INTEGER                     :: NMem,I,J
     !

   END SUBROUTINE CurviLinOverlap
!
!-------------------------------------------------------------------
!
   SUBROUTINE CollectINTCPast(RefStruct,RefGrad,IntCValues,IntCGrads, &
                              IntCs,GOpt,SCRPath,Print,PBCDim)
     REAL(DOUBLE),DIMENSION(:,:)  :: RefStruct,RefGrad
     TYPE(DBL_RNK2)               :: XYZAux,IntCGrads,IntCValues
     TYPE(DBL_VECT)               :: VectC,VectI,VectCG,VectX
     INTEGER                      :: NMem,NatmsLoc,NCart
     INTEGER                      :: Print,I,J,PBCDim
     TYPE(INTC)                   :: IntCs
     TYPE(GeomOpt)                :: GOpt
     CHARACTER(LEN=*)             :: SCRPath
     LOGICAL                      :: Print2
     !
     Print2=(Print>=DEBUG_GEOP_MAX)
     NMem=SIZE(RefStruct,2)
     IF(NMem/=SIZE(RefGrad,2)) CALL Halt('Dim err in CollectINTCPast')
     NCart=SIZE(RefStruct,1)
     NatmsLoc=NCart/3
     IF(NCart/=3*NatmsLoc) CALL Halt('Dim err 2 in CollectINTCPast')
     !
     CALL New(XYZAux,(/3,NatmsLoc/))
     CALL New(VectC,NCart)
     CALL New(VectCG,NCart)
     CALL New(VectX,NMem)
     !
     CALL New(IntCGrads,(/IntCs%N,NMem/))
     CALL New(IntCValues,(/IntCs%N,NMem/))
     CALL New(VectI,IntCs%N)
     !
     ! Calculate IntC gradients using the latest IntC-s
     !
     DO I=1,NMem
       DO J=1,NCart
         VectC%D(J)=RefStruct(J,I)
         VectCG%D(J)=RefGrad(J,I)
       ENDDO
       CALL CartRNK1ToCartRNK2(VectC%D,XYZAux%D)
       CALL RefreshBMatInfo(IntCs,XYZAux%D,GOpt%TrfCtrl,GOPt%GConvCrit,&
                     GOpt%CoordCtrl%LinCrit,GOpt%CoordCtrl%TorsLinCrit,&
                     PBCDim,Print,SCRPath)
       CALL CartToInternal(IntCs,VectCG%D,VectI%D,XYZAux%D,PBCDim, &
         GOpt%GrdTrf,GOpt%CoordCtrl,GOpt%TrfCtrl,Print,SCRPath)
       CALL RedundancyOff(VectI%D,SCRPath,Print,Messg_O='Q IntC Grads')
     ! CALL POffHardGc(IntCs,XYZAux%D,PBCDim,VectI%D,SCRPath,Print2)
       IntCGrads%D(:,I)=VectI%D
       CALL INTCValue(IntCs,XYZAux%D,PBCDim,GOpt%CoordCtrl%LinCrit, &
                      GOpt%CoordCtrl%TorsLinCrit)
       IntCValues%D(:,I)=IntCs%Value%D
     ENDDO
     CALL GrdConvrgd(GOpt%GOptStat,IntCs,VectI%D)
     !
     CALL Delete(VectX)
     CALL Delete(VectI)
     CALL Delete(VectC)
     CALL Delete(VectCG)
     CALL Delete(XYZAux)
   END SUBROUTINE CollectINTCPast
!
!----------------------------------------------------------------------
!
   SUBROUTINE LocalWeight(LWeight,Weights,IntCs,NCart,SCRPath, &
                          USQ_O)
     REAL(DOUBLE),DIMENSION(:,:) :: LWeight,Weights
     CHARACTER(LEN=*)            :: SCRPath
     TYPE(INT_VECT)              :: IGi,JGi,IGiT1,JGiT1,IGiT2,JGiT2
     INTEGER                     :: NCart,NZ,I,J,K1,K2,L1,L2,NMem
     INTEGER                     :: NPBC
     REAL(DOUBLE)                :: Weight,X1,X2,Y,W
     TYPE(INTC)                  :: IntCs
     TYPE(DBL_VECT)              :: Vect1,AGi,AGiT1,AGiT2
     REAL(DOUBLE),DIMENSION(:,:),OPTIONAL :: USQ_O
     !
     ! Calculate Topology (connectivity) of internal coordinates
     !
     NMem=SIZE(LWeight,2)
     CALL New(Vect1,IntCs%N)
     I=MAX(IntCs%N,NMem)
     !
     IF(PRESENT(USQ_O)) THEN
       K1=SIZE(USQ_O,1)
       CALL New(IGi,K1+1)
       CALL New(JGi,K1*K1)
       CALL New(AGi,K1*K1)
       NZ=0
       IGi%I(1)=1
       DO I=1,K1
         DO J=1,K1
           NZ=NZ+1
           AGi%D(NZ)=USQ_O(I,J)
           JGi%I(NZ)=J
         ENDDO
         IGi%I(I+1)=NZ+1
       ENDDO
       CALL ThreshMatr(IGi,JGi,AGi,1.D-7)
     ELSE
       CALL GetPattern(SCRPath,IntCs%N,NCart,IGi,JGi)
     ENDIF
     !
     ! Calculate RMS gradients for all internal coordinates
     ! touching a specific internal coord at a specific geometry
     !
     DO I=1,NMem
       DO J=1,IntCs%N
         IF(.NOT.IntCs%Active%L(J)) THEN
           LWeight(J,I)=1.D99
           CYCLE
         ENDIF
         NZ=0
         DO K1=IGi%I(J),IGi%I(J+1)-1
           L1=JGi%I(K1)
         ! IF(J==L1) CYCLE
           IF(.NOT.IntCs%Active%L(L1)) CYCLE
           NZ=NZ+1
           Vect1%D(NZ)=Weights(L1,I)
         ENDDO
         CALL Reorder1D(Vect1%D(1:NZ)) !!! increasing order
         Weight=Zero
         DO K1=1,NZ
           Weight=Weight+Vect1%D(K1)
         ENDDO
         IF(NZ/=0) THEN
         ! LWeight(J,I)=Weights(J,I)+Weight/DBLE(NZ)
           LWeight(J,I)=Weight/DBLE(NZ)
         ELSE
           LWeight(J,I)=Weights(J,I)
         ENDIF
       ENDDO
     ENDDO
     !
     CALL Delete(Vect1)
     IF(PRESENT(USQ_O)) CALL Delete(AGi)
     CALL Delete(JGi)
     CALL Delete(IGi)
   END SUBROUTINE LocalWeight
!
!------------------------------------------------------------------
!
   SUBROUTINE GetPattern(SCRPath,NIntC,NCart,IGi,JGi)
     CHARACTER(LEN=*)            :: SCRPath
     TYPE(INT_VECT)              :: IGi1,JGi1
     TYPE(INT_VECT)              :: IGi,JGi
     TYPE(INT_VECT)              :: ISpB,JSpB,ISpBt,JSpBt
     TYPE(DBL_VECT)              :: ASpB,ASpBt
     INTEGER                     :: NZ,NIntC,NCart,I,J,K,L
     !
     IF(NIntC<=0) RETURN
     CALL ReadBMATR(ISpB,JSpB,ASpB,TRIM(SCRPath)//'B')
     IF(NIntC/=SIZE(ISpB%I)-1) THEN
       CALL Halt('Dimensionality Error in GetPattern')
     ENDIF
     !
     ! Clean up connections via lattice parameters
     !
     L=NCart-9
     DO I=1,NIntC
       DO J=ISpB%I(I),ISpB%I(I+1)-1
         K=JSpB%I(J)
         IF(K>L) ASpB%D(J)=Zero
       ENDDO
     ENDDO
     CALL ThreshMatr(ISpB,JSpB,ASpB,1.D-7)
     !
     NZ=ISpB%I(NIntC+1)-1
     CALL New(ISpBt,NCart+1)
     CALL New(JSpBt,NZ)
     !
     CALL TransPose1x1(ISpB%I,JSpB%I,ASpB%D,NIntC,NCart, &
                       ISpBt%I,JSpBt%I,ASpB%D,'full')
     CALL MatMulSymbDriver(ISpB%I,JSpB%I,ISpBt%I,JSpBt%I, &
                           NIntC,NCart,NIntC,IGi,JGi)
   ! CALL MatMulSymbDriver(ISpB%I,JSpB%I,ISpBt%I,JSpBt%I, &
   !                       NIntC,NCart,NIntC,IGi1,JGi1)
     CALL Delete(ISpBt)
     CALL Delete(JSpBt)
     CALL Delete(ISpB)
     CALL Delete(JSpB)
     CALL Delete(ASpB)
     !
   ! CALL MatMulSymbDriver(IGi1%I,JGi1%I,IGi1%I,JGi1%I, &
   !                       NIntC,NIntC,NIntC,IGi,JGi)
   ! CALL Delete(IGi1)
   ! CALL Delete(JGi1)
   END SUBROUTINE GetPattern
!
!---------------------------------------------------------------------
!
   SUBROUTINE DisplFit(IntCs,IntCGrads,IntCValues,GHess,GCoordCtrl, &
                       PredVals,Displ,Path,SCRPath,NCart,iGEO, &
                       DoNorm,PrtFits,MixMat_O)
     TYPE(INTC)                 :: IntCs
     TYPE(DBL_VECT)             :: PredVals,Displ,DisplT
     REAL(DOUBLE),DIMENSION(:,:):: IntCGrads,IntCValues
     INTEGER                    :: I,J,NIntC,NDim,iGEO
     INTEGER                    :: NCart,NT
     CHARACTER(LEN=*)           :: Path,SCRPath
     CHARACTER(LEN=DCL)         :: Path2
     TYPE(Hessian)              :: GHess
     TYPE(CoordCtrl)            :: GCoordCtrl
     TYPE(DBL_RNK2)             :: IntCGradsT,IntCValuesT,FittedHessT
     TYPE(DBL_RNK2)             :: WeightsT,LWeightT,ABCT,ABC1T
     TYPE(DBL_RNK2)             :: RangeT,USQ
     TYPE(INTC)                 :: IntCsT
     TYPE(INT_VECT)             :: NDegsT
     REAL(DOUBLE),DIMENSION(:,:),OPTIONAL:: MixMat_O
     LOGICAL                    :: PrtFits,DoNorm
     !
     NIntC=IntCs%N
     NDim=SIZE(IntCGrads,2)
     IF(PRESENT(MixMat_O)) THEN
       NT=SIZE(MixMat_O,2)
       CALL New(USQ,(/NT,NT/))
       USQ%D=One
      !USQ%D=Zero
      !DO J=1,NT ; USQ%D(J,J)=One ; ENDDO
     ELSE
       NT=NIntC
     ENDIF
     CALL New(PredVals,NT)
     CALL New(Displ,NT)
     Path2=TRIM(Path)//'_'//TRIM(IntToChar(IGEO))
     !
     CALL New(ABCT,(/NT,4/))
     CALL New(ABC1T,(/NT,4/))
     CALL New(DisplT,NT)
     CALL New(FittedHessT,(/NT,NDim/))
     CALL New(WeightsT,(/NT,NDim/))
     CALL New(RangeT,(/NT,2/))
     CALL New(IntCGradsT,(/NT,NDim/))
     CALL New(IntCValuesT,(/NT,NDim/))
     CALL New(LWeightT,(/NT,NDim/))
     CALL New(NDegsT,NT)
     CALL New(IntCsT,NT)
     !
     IntCs%PredGrad%D=Zero
     IF(PRESENT(MixMat_O)) THEN
       ABC1T%D(:,:)=Zero
       ABC1T%D(:,2)=One
       ABCT%D=ABC1T%D
       IntCsT%Active%L=.TRUE.
       CALL DGEMM_TNc(NT,NIntC,NDim,One,Zero,MixMat_O,IntCGrads,IntCGradsT%D)
       CALL DGEMM_TNc(NT,NIntC,NDim,One,Zero,MixMat_O,IntCValues,IntCValuesT%D)
     ELSE
       IntCValuesT%D=IntCValues
       IntCGradsT%D=IntCGrads
       CALL SetEq(IntCs,IntCsT,1,NT,1)
       CALL InitialABC(ABC1T%D,IntCs,GHess)
       ABCT%D=ABC1T%D
     ENDIF
     !
     CALL PrepPrimW(WeightsT%D,IntCGradsT%D,IntCsT)
     CALL CalcHessian(FittedHessT%D,ABC1T%D)
     CALL SecondWeight(WeightsT%D,FittedHessT%D)
     IF(PRESENT(MixMat_O)) THEN
       CALL LocalWeight(LWeightT%D,WeightsT%D,IntCsT,NCart, &
                        SCRPath,USQ_O=USQ%D)
     ELSE
       CALL LocalWeight(LWeightT%D,WeightsT%D,IntCsT,NCart, &
                        SCRPath)
     ENDIF
     CALL LQFit(IntCValuesT%D,IntCGradsT%D,LWeightT%D,IntCsT,ABCT%D, &
                RangeT%D,NDegsT%I,Zero,GCoordCtrl%DoQFilter)
              ! RangeT%D,NDegsT%I,Zero,.FALSE.)
     CALL DoPredict(ABCT%D,IntCValuesT%D,IntCGradsT%D,IntCsT, &
                    NDegsT%I,Path2,RangeT%D,DoNorm)
     CALL CleanRange(DisplT%D,RangeT%D,IntCs%Def%C,IntCsT%PredVal%D, &
                     IntCValuesT%D(:,NDim),NDim)
     IntCsT%PredVal%D=IntCValuesT%D(:,NDim)+DisplT%D
     ! set up constraints
     CALL SetConstraints(IntCsT,IntCsT%PredVal%D)
     DisplT%D=IntCsT%PredVal%D-IntCValuesT%D(:,NDim)
     !
     IF(PRTFits) THEN
       CALL PrtFitM(IntCValuesT%D,IntCGradsT%D,ABCT%D,IntCsT,Path2)
     ENDIF
     PredVals%D=IntCsT%PredVal%D
     Displ%D=DisplT%D
     !
     IF(PRESENT(MixMat_O)) THEN
       CALL Delete(USQ)
     ENDIF
     CALL Delete(NDegsT)
     CALL Delete(FittedHessT)
     CALL Delete(WeightsT)
     CALL Delete(IntCsT)
     CALL Delete(DisplT)
     CALL Delete(RangeT)
     CALL Delete(ABC1T)
     CALL Delete(ABCT)
     CALL Delete(IntCGradsT)
     CALL Delete(IntCValuesT)
     CALL Delete(LWeightT)
   END SUBROUTINE DisplFit
!
!---------------------------------------------------------------------
!
   SUBROUTINE CleanRange(DisplT,RangeT,IntCDef,PredVal,Values,NDim)
     REAL(DOUBLE),DIMENSION(:)     :: DisplT,PredVal,Values
     REAL(DOUBLE),DIMENSION(:,:)   :: RangeT
     CHARACTER(LEN=*),DIMENSION(:) :: IntCDef
     REAL(DOUBLE)                  :: Range,Range1,Range2,Displ,X
     INTEGER                       :: I,NT,NDim
     !
     NT=SIZE(DisplT)
     DO I=1,NT
       DisplT(I)=PredVal(I)-Values(I)
     ! IF(HasLattice(IntCDef(I))) CYCLE
       Range1=RangeT(I,1)
       Range2=RangeT(I,2)
       Range=Range2-Range1
      !IF(NDim==2) THEN
      !  Range=Two*Range
      !ENDIF
       Displ=DisplT(I)
       IF(Range1-PredVal(I)>Range) THEN
         PredVal(I)=Range1-Range
       ELSE IF(PredVal(I)-Range2>Range) THEN
         PredVal(I)=Range2+Range
       ENDIF
       DisplT(I)=PredVal(I)-Values(I)
     ENDDO
   END SUBROUTINE CleanRange
!
!---------------------------------------------------------------------
!
   SUBROUTINE BackToPrims(IntCValues,IntCs,IntCsT,DisplT,UMat,Displ)
     REAL(DOUBLE),DIMENSION(:)   :: Displ,DisplT
     REAL(DOUBLE),DIMENSION(:,:) :: UMat,IntCValues
     TYPE(INTC)                  :: IntCs,IntCsT
     INTEGER                     :: NDim,NIntC,J,NT
     !
     NIntC=SIZE(UMat,1)
     NT=SIZE(UMat,2)
     NDim=SIZE(IntCValues,2)
     CALL DGEMM_NNc(NIntC,NT,1,One,Zero, &
                    UMat,DisplT,Displ)
     IntCs%PredVal%D=IntCValues(:,NDim)+Displ
   ! CALL DGEMM_NNc(NIntC,NT,1,One,Zero, &
   !                UMat,IntCsT%PredVal%D,IntCs%PredVal%D)
   ! Displ=IntCs%PredVal%D-IntCValues(:,NDim)
   END SUBROUTINE BackToPrims
!
!---------------------------------------------------------------------
!
   SUBROUTINE UnitaryTR(IntCs,IntCGrads,IntCValues,UMat,NT)
     REAL(DOUBLE),DIMENSION(:,:) :: IntCGrads,IntCvalues
     TYPE(DBL_RNK2)              :: UMatS,UMatSQ,UMat,USQ,ScaledGrads
     TYPE(INTC)                  :: IntCs
     INTEGER                     :: I,J,NDim,NIntC,NT,INFO
     REAL(DOUBLE)                :: FC
     !
     NIntC=SIZE(IntCGrads,1)
     NDim=SIZE(IntCGrads,2)
     CALL New(ScaledGrads,(/NIntC,NDim/))
     CALL New(UMatS,(/NDim,NDim/))
     CALL New(UMatSQ,(/NDim,NDim/))
     !
     DO I=1,NIntC
       IF(IntCs%Def%C(I)(1:4)=='STRE') THEN
         FC=One/0.5D0
       ELSE IF(IntCs%Def%C(I)(1:4)=='BEND') THEN
         FC=One/0.2D0
       ELSE IF(IntCs%Def%C(I)(1:4)=='LINB') THEN
         FC=One/0.2D0
       ELSE IF(IntCs%Def%C(I)(1:4)=='OUTP') THEN
         FC=One/0.2D0
       ELSE IF(IntCs%Def%C(I)(1:4)=='TORS') THEN
         FC=One/0.1D0
       ELSE
         FC=One
       ENDIF
       FC=SQRT(FC)
       DO J=1,NDim
        !ScaledGrads%D(I,J)=IntCGrads(I,J)
         ScaledGrads%D(I,J)=FC*IntCGrads(I,J)
        !ScaledGrads%D(I,J)=IntCGrads(I,J)*FC*IntCGrads(I,J)
       ENDDO
     ENDDO
     !
   ! CALL DGEMM_TNc(NDim,NIntC,NDim,One,Zero, &
   !                IntCGrads,IntCGrads,UmatS%D)
     CALL DGEMM_TNc(NDim,NIntC,NDim,One,Zero, &
                    ScaledGrads%D,ScaledGrads%D,UmatS%D)
       CALL SetDSYEVWork(NDim)
       BLKVECT%D=UMatS%D
       CALL DSYEV('V','U',NDim,BLKVECT%D,BIGBLOK,BLKVALS%D, &
       BLKWORK%D,BLKLWORK,INFO)
       IF(INFO/=SUCCEED) &
       CALL Halt('DSYEV hosed in UnitaryTR. INFO='&
                  //TRIM(IntToChar(INFO)))
       NT=0
       UMatSQ%D=Zero
       DO J=NDim,1,-1
        !IF(BLKVALS%D(J)/ABS(BLKVALS%D(NDim))>1.D-7) THEN
         IF(BLKVALS%D(J)>1.D-14) THEN
           NT=NT+1
           UMatSQ%D(:,NT)=BLKVECT%D(:,J)/SQRT(BLKVALS%D(J))
         ENDIF
       ENDDO
     CALL UnSetDSYEVWork()
     CALL Delete(UMatS)
     !
     CALL New(USQ,(/NDim,NT/))
     DO I=1,NDim
       DO J=1,NT
         USQ%D(I,J)=UMatSQ%D(I,J)
       ENDDO
     ENDDO
     CALL Delete(UMatSQ)
     !
     CALL New(UMat,(/NIntC,NT/))
     CALL DGEMM_NNc(NIntC,NDim,NT,One,Zero, &
                    ScaledGrads%D,USQ%D,UMat%D)
     CALL Delete(USQ)
     !
   ! CALL PPrint(ScaledGrads%D,'ScaledGrads%D',unit_o=6)
   ! CALL PPrint(UMat%D,'UMat%D',unit_o=6)
   ! CALL New(USQ,(/NT,NT/))
   ! CALL DGEMM_TNc(NT,NIntC,NT,One,Zero,UMat%D,UMat%D,USQ%D)
   ! CALL PPrint(USQ,'Unit Mat?',unit_o=6)
   ! CALL Delete(USQ)
     !
     CALL Delete(ScaledGrads)
   END SUBROUTINE UnitaryTR
!
!---------------------------------------------------------------------
!
   SUBROUTINE InitialABC(ABC,IntCs,GHess)
     REAL(DOUBLE),DIMENSION(:,:) :: ABC
     TYPE(INTC)                  :: IntCs
     TYPE(Hessian)               :: GHess
     INTEGER                     :: I
     REAL(DOUBLE)                :: Stre,Bend,LinB,OutP,Tors
     !
    !Stre=1.0D0
    !Bend=0.1D0
    !LinB=0.1D0
    !OutP=0.1D0
    !Tors=0.01D0
     Stre=GHess%Stre
     Bend=GHess%Bend
     LinB=GHess%LinB
     OutP=GHess%OutP
     Tors=GHess%Tors
     DO I=1,IntCs%N
       ABC(I,:)=Zero
       IF(IntCs%Def%C(I)(1:5)=='STRE ') THEN
         ABC(I,2)=Stre
       ELSE IF(IntCs%Def%C(I)(1:4)=='BEND') THEN
         ABC(I,2)=Bend
       ELSE IF(IntCs%Def%C(I)(1:4)=='LINB') THEN
         ABC(I,2)=LinB
       ELSE IF(IntCs%Def%C(I)(1:4)=='OUTP') THEN
         ABC(I,2)=OutP
       ELSE IF(IntCs%Def%C(I)(1:4)=='TORS') THEN
         ABC(I,2)=Tors
       ELSE IF(IntCs%Def%C(I)(1:4)=='CART') THEN
         ABC(I,2)=Stre
       ELSE
         ABC(I,2)=5.D0 ! for lattice
       ENDIF
     ENDDO
   END SUBROUTINE InitialABC
!
!---------------------------------------------------------------------
!
   SUBROUTINE InvertWeight(Weights,InvWeight)
     REAL(DOUBLE),DIMENSION(:,:) :: Weights,InvWeight
     INTEGER                     :: NDim,NIntC,I,J
     !
     NIntC=SIZE(Weights,1)
     NDim=SIZE(Weights,2)
     !
     DO I=1,NDim
         DO J=1,NIntC
           InvWeight(J,I)=One/(Weights(J,I)+1.D-20)
         ENDDO
     ENDDO
   END SUBROUTINE InvertWeight
!
!---------------------------------------------------------------------
!
   SUBROUTINE PrtFitM(IntCValues,IntCGrads,ABC,IntCs,Path)
     REAL(DOUBLE),DIMENSION(:,:) :: IntCValues,IntCGrads,ABC
     TYPE(INTC)                  :: IntCs
     INTEGER                     :: NIntC,NDim,I,J,I1
     CHARACTER(LEN=DCL)          :: Path2
     CHARACTER(LEN=*)            :: Path
     TYPE(DBL_VECT)              :: VectX,VectY
     !
     NIntC=SIZE(IntCValues,1)
     NDim=SIZE(IntCValues,2)
     CALL New(VectX,NDim)
     CALL New(VectY,NDim)
     DO I=1,NIntC
       DO J=1,NDim
         VectX%D(J)=IntCValues(I,J)
         VectY%D(J)=IntCGrads(I,J)
       ENDDO
       Path2=TRIM(Path)//'_'//TRIM(IntToChar(I))
       I1=1
       CALL PrtFits(VectX%D(I1:NDim),VectY%D(I1:NDim),TRIM(Path2), &
                    ABC(I,1),ABC(I,2),ABC(I,3),ABC(I,4), &
                    IntCs%PredVal%D(I),IntCs%PredGrad%D(I), &
                    IntCs%Def%C(I))
     ENDDO
     CALL Delete(VectX)
     CALL Delete(VectY)
   END SUBROUTINE PrtFitM
!
!---------------------------------------------------------------------
!
   SUBROUTINE CleanDispl(IntCValues,IntCs,Displ,MaxStre,MaxAngle)
     REAL(DOUBLE),DIMENSION(:,:) :: IntCValues
     REAL(DOUBLE),DIMENSION(:)   :: Displ
     TYPE(INTC)                  :: IntCs
     INTEGER                     :: NIntC,NDim,I,J,I1
     REAL(DOUBLE)                :: MaxStre,MaxAngle

     CALL MondoLog(DEBUG_NONE, "CleanDispl", "MaxStre = "//TRIM(DblToChar(MaxStre))//", MaxAngle = "//TRIM(DblToChar(MaxAngle)))

     NIntC=SIZE(IntCValues,1)
     NDim=SIZE(IntCValues,2)
     DO I=1,NIntC
       IF(IntCs%Def%C(I)(1:4)=='LINB'.OR. &
          IntCs%Def%C(I)(1:4)=='OUTP'.OR. &
          IntCs%Def%C(I)(1:4)=='TORS') THEN
         CALL Loose2PIs(IntCs%PredVal%D(I))
         CALL LinBTo180(IntCs%PredVal%D(I))
       ELSE IF(IntCs%Def%C(I)(1:4)=='BEND') THEN
         CALL ChkBendLim(IntCValues(I,NDim),IntCs%PredVal%D(I),5.D0)
         CALL BendTo180(IntCs%PredVal%D(I))
       ENDIF

       Displ(I)=IntCs%PredVal%D(I)-IntCValues(I,NDim)
       CALL MapDAngle(IntCs%Def%C(I),IntCValues(I,NDim),Displ(I))

       CALL CtrlDispl(IntCs%Def%C(I),IntCs%Value%D(I),Displ(I),One,MaxStre,MaxAngle)
       I1=1
       IntCs%PredVal%D(I)=IntCValues(I,NDim)+Displ(I)
       IF(IntCs%Def%C(I)(1:4)=='STRE') THEN
         CALL ChkStreLim(IntCValues(I,NDim),IntCs%PredVal%D(I),0.3D0)
       ENDIF
     ENDDO
   END SUBROUTINE CleanDispl
!
!---------------------------------------------------------------------
!
   SUBROUTINE DoPredict(ABC,IntCValues,IntCGrads,IntCs,NDegs, &
                        Path,Range,DoNorm)
     REAL(DOUBLE),DIMENSION(:,:) :: ABC,IntCValues,IntCGrads,Range
     INTEGER,DIMENSION(:)        :: NDegs
     TYPE(INTC)                  :: IntCs
     INTEGER                     :: NIntC,NDim,I,J
     TYPE(DBL_VECT)              :: VectX,VectY
     REAL(DOUBLE)                :: MaxX,MinX
     CHARACTER(LEN=*)            :: Path
     LOGICAL                     :: DoNorm
     !
     NIntC=SIZE(IntCValues,1)
     NDim=SIZE(IntCValues,2)
     CALL New(VectX,NDim)
     CALL New(VectY,NDim)
     DO I=1,NIntC
       IF(.NOT.IntCs%Active%L(I)) THEN
         IntCs%PredVal%D(I)=Zero
         CYCLE
       ENDIF
       DO J=1,NDim
         VectX%D(J)=IntCValues(I,J)
         VectY%D(J)=IntCGrads(I,J)
       ENDDO
       MaxX=MAXVAL(VectX%D)
       MinX=MINVAL(VectX%D)
       CALL Predict(ABC(I,1),ABC(I,2),ABC(I,3),ABC(I,4), &
                    VectX%D,VectY%D,IntCs%PredVal%D(I), &
                    IntCs%PredGrad%D(I),NDegs(I),DoNorm)
     !
     ENDDO
     CALL Delete(VectX)
     CALL Delete(VectY)
    !CALL PrtPred(IntCs,IntCValues,IntCGrads, &
    !             IntCs%PredVal%D,IntCs%PredGrad%D,Path)
   END SUBROUTINE DoPredict
!
!---------------------------------------------------------------------
!
   SUBROUTINE SecondWeight(Weights,Hessian)
     REAL(DOUBLE),DIMENSION(:,:) :: Weights,Hessian
     REAL(DOUBLE)                :: X
     INTEGER                     :: NIntC,NDim,I,J
     !
     NIntC=SIZE(Weights,1)
     NDim=SIZE(Weights,2)
     !
     DO I=1,NDim
         DO J=1,NIntC
           X=Weights(J,I)/(ABS(Hessian(J,I))+1.D-20)
           Weights(J,I)=X
         ENDDO
     ENDDO
   END SUBROUTINE SecondWeight
!
!---------------------------------------------------------------------
!
   SUBROUTINE CalcHessian(Hessian,ABC)
     REAL(DOUBLE),DIMENSION(:,:) :: Hessian,ABC
     REAL(DOUBLE)                :: X
     INTEGER                     :: NIntC,NDim,I,J
     !
     NIntC=SIZE(Hessian,1)
     NDim=SIZE(Hessian,2)
     !
     DO I=1,NDim
       DO J=1,NIntC
        !X=IntCValues(J,I)
        !Hessian(J,I)=ABC(J,2)+Two*ABC(J,3)*X+Three*ABC(J,4)*X*X
         Hessian(J,I)=ABC(J,2)
       ENDDO
     ENDDO
   END SUBROUTINE CalcHessian
!
!---------------------------------------------------------------------
!
   SUBROUTINE LQFit(IntCValues,IntCGrads,Weights,IntCs, &
                    ABC,Range,NDegs,EPS,DoReord)
     REAL(DOUBLE),DIMENSION(:,:) :: IntCValues,IntCGrads,Weights,ABC
     INTEGER,DIMENSION(:)        :: NDegs
     REAL(DOUBLE)                :: MaxX,MinX,MaxY,MinY,Chi2V,EPS,X
     TYPE(INTC)                  :: IntCs
     INTEGER                     :: NIntC,NDim,NDim2,I,J,NDeg,NDeg1,I1
     INTEGER                     :: IStart
     TYPE(DBL_VECT)              :: RMSErr,VectX,VectY,VectAux
     INTEGER                     :: MaxDeg,NS
     TYPE(DBL_VECT)              :: Work1,Work2,VectFit1
     TYPE(INT_VECT)              :: IWork
     REAL(DOUBLE),DIMENSION(:,:) :: Range
     REAL(DOUBLE)                :: Conv
     LOGICAL                     :: DoReord
     !
     NIntC=SIZE(IntCValues,1)
     IntCs%Predgrad%D=Zero
     NDim=SIZE(IntCValues,2)
     NDeg1=NDim
     CALL New(RMSErr,NDim)
     CALL New(VectX,NDim)
     CALL New(VectY,NDim)
     CALL New(VectFit1,NDim)
     CALL New(VectAux,4)
     !
     MaxDeg=3
     NS=MaxDeg+1
     CALL New(Work1,NDim)
     CALL New(Work2,4*NDim*NS+2*NS*NS)
     CALL New(IWork,NDim)
     !
     DO I=1,NIntC
       IF(.NOT.IntCs%Active%L(I)) CYCLE
       DO J=1,NDim
         VectX%D(J)=IntCValues(I,J)
         VectY%D(J)=IntCGrads(I,J)
         RMSErr%D(J)=Weights(I,J)
       ENDDO
       MinX=MINVAL(VectX%D)
       MaxX=MAXVAL(VectX%D)
       !
       IF(MaxX-MinX<1.D-6) THEN
         ! Do steepest descent with empirical force constants
         ! or  5.0D0 if range is too small for fitting.
         ! Control over stepsize is going to be modified to
         ! standard DiagHess range (0.3 a.u.)
         Range(I,1)=VectX%D(NDim)-0.10D0
         Range(I,2)=VectX%D(NDim)+0.10D0
        !ABC(I,2)= 5.0D0
         ABC(I,1)=VectY%D(NDim)-ABC(I,2)*VectX%D(NDim)
         NDegs(I)=1
         CYCLE
       ENDIF
       !
       ! invert to get weights
       DO J=1,NDim
         RMSErr%D(J)=One/(RMSErr%D(J)+1.D-20)
       ENDDO
       RMSErr%D=RMSErr%D/SUM(RMSErr%D)
       !
    !  IF(DoReOrd.AND.NDim>=MaxMem) THEN
       IF(DoReOrd) THEN
    ! warning! do not do reordering back to intcgrads/values
    ! that can mess up the relationship of Cartesians/INTC-s
    ! as predicted displacement refers to last Cartesian set
         CALL QTest(VectX%D,VectY%D,RMSErr%D,Work2%D,IWork%I,NDim,IStart)
         Range(I,1)=MINVAL(VectX%D(IStart:NDim))
         Range(I,2)=MAXVAL(VectX%D(IStart:NDim))
       ELSE
         Range(I,1)=MinX
         Range(I,2)=MaxX
         IStart=1
       ENDIF
       !
       VectAux%D=ABC(I,:)
       CALL BasicFit(VectX%D(IStart:NDim),VectY%D(IStart:NDim), &
                     RMSErr%D(IStart:NDim),Work1%D,Work2%D,EPS,Chi2V, &
                     VectAux%D,NDeg,NDim-IStart+1)
       ABC(I,:)=VectAux%D
       NDegs(I)=NDeg
     ENDDO
     !
     CALL Delete(VectAux)
     CALL Delete(VectFit1)
     CALL Delete(RMSErr)
     CALL Delete(VectX)
     CALL Delete(VectY)
     CALL Delete(Work1)
     CALL Delete(Work2)
     CALL Delete(IWork)
   END SUBROUTINE LQFit
!
!---------------------------------------------------------------------
!
   SUBROUTINE QTest2(VectX,VectY,RMSErr,Work,IWork, &
                     MaxX,MinX,NDim,IStart)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY,RMSErr,Work
     INTEGER,DIMENSION(:)      :: IWork
     INTEGER                   :: IStart,NDim,NDim2,NDim3,I,J
     REAL(DOUBLE)              :: Q,QTab,Range,R1,R2,XC,MaxX,MinX
     !
     IStart=1
     IF(NDim<3) RETURN
     !
     CALL ReorderN(RMSErr,IWork,NDim)
     DO I=1,NDim
       Work(I)=ABS(VectX(I)-VectX(NDim)) !*RMSErr(I)
     ENDDO
     CALL ReorderN(Work,IWork,NDim)
     NDim2=2*NDim
     NDim3=3*NDim
     ! fill old arrays in work-s by new order
     DO I=1,NDim
       J=IWork(I)
       Work(I)=RMSErr(J)
       Work(NDim+I)=VectX(J)
       Work(NDim2+I)=VectY(J)
     ENDDO
     ! reorder original arrays
     DO J=1,NDim
       RMSErr(J)=Work(J)
       VectX(J)=Work(NDim+J)
       VectY(J)=Work(NDim2+J)
      !Work(J)=LOG(Work(J))
     ENDDO
     !
     ! do selection
     DO I=1,NDim
       Work(I)=ABS(VectX(I)-VectX(NDim)) !*RMSErr(I)
     ENDDO
     Range=MAXVAL(Work(1:NDim))
   ! Range=MaxX-MinX
     QTab=QTest90(MIN(10,NDim))
     DO J=NDim-2,1,-1
       Q=ABS(Work(J)-Work(J+1))/Range
       IF(Q>QTab) THEN
         IStart=J+1
         EXIT
       ENDIF
     ENDDO
   END SUBROUTINE QTest2
!
!---------------------------------------------------------------------
!
   SUBROUTINE QTest(VectX,VectY,RMSErr,Work,IWork,NDim,IStart)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY,RMSErr,Work
     INTEGER,DIMENSION(:)      :: IWork
     INTEGER                   :: IStart,NDim,NDim2,NDim3,I,J,ILeft
     REAL(DOUBLE)              :: Q,QTab,Range,Fluct
     !
     ! filtering based on weights
     IStart=1
     ILeft=3
     IF(NDim<ILeft+1) RETURN
     DO J=1,NDim ; Work(J)=RMSErr(J) ; ENDDO
     Q=MAXVAL(RMSErr)
     DO J=NDim-ILeft+1,NDim ; Work(J)=Q+DBLE(J) ; ENDDO
     CALL ReorderI(Work,IWork,NDim)
     NDim2=2*NDim
     NDim3=3*NDim
     ! fill old arrays in work-s by new order
     DO I=1,NDim
       J=IWork(I)
       Work(I)=RMSErr(J)
       Work(NDim+I)=VectX(J)
       Work(NDim2+I)=VectY(J)
     ENDDO
     ! reorder original arrays
     DO J=1,NDim
       RMSErr(J)=Work(J)
       VectX(J)=Work(NDim+J)
       VectY(J)=Work(NDim2+J)
      !Work(J)=LOG(Work(J))
     ENDDO
     !
   ! IF(VectY(NDim)*VectY(NDim-1)<Zero) THEN
   !   IStart=NDim-1
   !   RETURN
   ! ENDIF
     !
     !QTab=QTest90(MIN(10,NDim-I+1))
    !QTab=One/DBLE(NDim)
    !QTab=QTab-QTab*QTab
    !QTab=0.001D0
     J=NDim-ILeft+1
     Q=SUM(RMSErr(J:NDim))
     DO I=J-1,1,-1
      !IF(Q>0.9999D0) THEN
       IF(Q>0.95D0) THEN
         IStart=I+1
         EXIT
       ENDIF
       Q=Q+RMSErr(I)
     ENDDO
   END SUBROUTINE QTest
!
!---------------------------------------------------------------------
!
   SUBROUTINE PrepPrimW(Weights,IntCGrads,IntCs)
     REAL(DOUBLE),DIMENSION(:,:) :: Weights,IntCGrads
     TYPE(INTC)                  :: IntCs
     REAL(DOUBLE)                :: X
     INTEGER                     :: NIntC,NDim,I,J
     !
     NIntC=SIZE(IntCGrads,1)
     NDim=SIZE(IntCGrads,2)
     !
     DO I=1,NDim
       DO J=1,NIntC
         IF(IntCs%Active%L(J)) THEN
           X=IntCGrads(J,I)
           Weights(J,I)=X*X
         ELSE
           Weights(J,I)=1.D99
         ENDIF
       ENDDO
     ENDDO
   END SUBROUTINE PrepPrimW
!
!---------------------------------------------------------------------
!
   SUBROUTINE BasicFit(VectX,VectY,RMSErr,Work1,Work2,EPS,Chi2V, &
                       Params,NDeg,NDim)
     ! arrays for output
     INTEGER       :: MaxDeg,NDeg,NDim,I,II,J
     REAL(DOUBLE)  :: VectX(:),VectY(:),RMSErr(:)
     REAL(DOUBLE)  :: Params(:),Work1(:),Work2(:)
     REAL(DOUBLE)  :: X0,Y0,EPS,Chi2V,MinY,MaxY,NegW
     INTEGER       :: I1,I2,I3,I4,I5,I6,I7,I8,I9,IStart
     LOGICAL       :: DoQFit
     !
     MaxDeg=3
     MaxDeg=MAX(MIN(NDim-2,MaxDeg),1)
     !
     Params=Zero
     CALL Chi2Fit(VectX,VectY,RMSErr,Work1, &
                  Work2,MaxDeg,NDeg,Params,Chi2V,EPS,DoQFit)
   END SUBROUTINE BasicFit
!
!---------------------------------------------------------------------
!
   SUBROUTINE FilterBow(VectX,VectY,MaxDeg)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY
     REAL(DOUBLE)              :: X,X1,X2,MinX,MaxX,Range,Delta
     INTEGER                   :: MaxDeg,I,NDim,II
     !
     NDim=SIZE(VectX)
     MinX=MINVAL(VectX)
     MaxX=MAXVAL(VectX)
     Delta=0.25D0*(MaxX-MinX)
     X1=MinX+Delta
     X2=X1+Delta+Delta
     II=0
     DO I=1,NDim
       X=VectX(I)
       IF(X>X1.AND.X<X2) THEN
         II=1
         EXIT
       ENDIF
     ENDDO
     IF(II/=1) MaxDeg=1
   END SUBROUTINE FilterBow
!
!---------------------------------------------------------------------
!
   SUBROUTINE FilterIEnd(VectX,VectY,IEnd)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY
     INTEGER                   :: IEnd,I,J,NDim
     !
     NDim=SIZE(VectX)
     IF(VectY(NDim)*VectY(NDim-1)<Zero) THEN
       IEnd=NDim-1
       RETURN
     ENDIF
   END SUBROUTINE FilterIEnd
!
!---------------------------------------------------------------------
!
   SUBROUTINE PrtPred(IntCs,IntCValues,IntCGrads, &
                      FitVal,PredGrad,Path)
     INTEGER                   :: I,NDim,NIntC
     TYPE(INTC)                :: IntCs
     REAL(DOUBLE)                :: Conv
     REAL(DOUBLE),DIMENSION(:,:) :: IntCValues,IntCGrads
     REAL(DOUBLE),DIMENSION(:) :: FitVal,PredGrad
     CHARACTER(LEN=*)          :: Path
     !
     NIntC=IntCs%N
     NDim=SIZE(IntCValues,2)
     OPEN(UNIT=91,FILE=TRIM(Path)//'Pred',STATUS='UNKNOWN')
     DO I=1,NIntC
       IF(HasAngle(IntCs%Def%C(I))) THEN
         Conv=180.D0/PI
       ELSE
         Conv=One/AngstromsToAU
       ENDIF
       WRITE(91,12) I,IntCs%Def%C(I)(1:5),Conv*IntCValues(I,NDim), &
                    IntCGrads(I,NDim),Conv*FitVal(I),PredGrad(I)

     ENDDO
     12 FORMAT(I3,2X,A5,2X,3F12.6,F20.6)
     CLOSE(91)
   END SUBROUTINE PrtPred
!
!---------------------------------------------------------------------
!
   SUBROUTINE PrtFits(VectX,VectY,Path,A,B,C,D,FitVal,PredGrad,Def)
     INTEGER                  :: NDim,J
     REAL(DOUBLE),DIMENSION(:):: VectX,VectY
     REAL(DOUBLE)             :: A,B,C,D,FitVal,PredGrad
     REAL(DOUBLE)             :: Conv,ConvI,ConvI2,ConvI3
     CHARACTER(LEN=*)         :: Path,Def
     !
     NDim=SIZE(VectX)
     IF(Def(1:4)=='STRE') THEN
       Conv=One/AngstromsToAU
       ConvI=AngstromsToAU
     ELSE IF(HasAngle(Def)) THEN
       Conv=180.D0/PI
       ConvI=PI/180.D0
     ELSE
       Conv=One
       ConvI=One
     ENDIF
     ConvI2=ConvI*ConvI
     ConvI3=ConvI2*ConvI
     !
     OPEN(UNIT=91,FILE=TRIM(Path)//'_Data',STATUS='UNKNOWN')
       DO J=1,NDim
         WRITE(91,11) J,Conv*VectX(J),VectY(J),FitVal*Conv,PredGrad
       ENDDO
       11 FORMAT(I3,2X,5F40.20)
     CLOSE(91)
     !
     OPEN(UNIT=91,FILE=TRIM(Path)//'_Params',STATUS='UNKNOWN')
         WRITE(91,12) A,B*ConvI,C*ConvI2,D*ConvI3
     CLOSE(91)
     !
     OPEN(UNIT=91,FILE=TRIM(Path)//'_Pred',STATUS='UNKNOWN')
         WRITE(91,12) FitVal*Conv,PredGrad
     CLOSE(91)
     12 FORMAT(5F40.20)
   END SUBROUTINE PrtFits
!
!---------------------------------------------------------------------
!
   SUBROUTINE ChkBendLim(Val,FitVal,DeltaMax)
     REAL(DOUBLE) :: Val,FitVal,DeltaMax
     IF(FitVal>PI.OR.FitVal<Zero) THEN
       FitVal=Val+SIGN(DeltaMax*PI/180.D0,FitVal-Val)
     ENDIF
   END SUBROUTINE ChkBendLim
!
!---------------------------------------------------------------------
!
   SUBROUTINE ChkStreLim(Val,FitVal,DeltaMax)
     REAL(DOUBLE) :: Val,FitVal,DeltaMax
     IF(FitVal<Zero.OR.FitVal>Val*Two) THEN
       FitVal=Val+DeltaMax
     ENDIF
   END SUBROUTINE ChkStreLim
!
!---------------------------------------------------------------------
!
   SUBROUTINE CubicRoots(A1,B1,C1,D1,X1,X2,X3)
     REAL(DOUBLE)              :: A1,B1,C1,D1,X1,X2,X3
     REAL(DOUBLE)              :: A,B,C
     REAL(DOUBLE)              :: Q,R,S,Q3,R2,SQ
     REAL(DOUBLE)              :: AP,BP,Theta
     INTEGER                   :: I,J
     !
     C=A1/D1
     B=B1/D1
     A=C1/D1
     !
     Q=(A*A-3.D0*B)/9.D0
     R=(2.D0*A**3-9.D0*A*B+27.D0*C)/54.D0
     SQ=Two*SQRT(ABS(Q))
     Q3=Q**3
     R2=R*R
     S=R2-Q3
     IF(S<Zero) THEN
       Theta=ACOS(R/SQRT(ABS(Q3)))
       X1=-SQ*COS(Theta/Three)-A/Three
       X2=-SQ*COS((Theta+TwoPi)/Three)-A/Three
       X3=-SQ*COS((Theta-TwoPi)/Three)-A/Three
     ELSE
       AP=-SIGN(One,R)*(ABS(R)+SQRT(S))**(One/Three)
       IF(ABS(AP)>1.D-6) THEN
         BP=Q/AP
       ELSE
         BP=Zero
       ENDIF
       X1=(AP+BP)-A/Three
       IF(ABS(AP-BP)>1.D-6) THEN
         X2=X1
         X3=X2
       ELSE
         X2=(AP+BP)/Two-A/Three
         X3=X2
       ENDIF
     ENDIF
   END SUBROUTINE CubicRoots
!
!---------------------------------------------------------------------
!
   SUBROUTINE Predict(A,B,C,D,VectX,VectY,FitVal,PredGrad,NDeg,DoNorm)
     REAL(DOUBLE)              :: A,B,C,D,FitVal,PredGrad,PredHess
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY
     REAL(DOUBLE)              :: TwoC,ThreeD
     REAL(DOUBLE)              :: X1,X2,G1,G2,X0,Y0,G0,X3,G3,D0,D1,D2,D3
     REAL(DOUBLE)              :: Det,AP,BP,CP,Delta,H0,Q
     REAL(DOUBLE)              :: LastX,LastY,RelErr,AbsErr,Guess
     INTEGER                   :: I,J,NDim,NDeg,IFlag
     LOGICAL                   :: DoNorm
     !
     NDim=SIZE(VectX)
     LastX=VectX(NDim)
     LastY=VectY(NDim)
     TwoC=Two*C
     ThreeD=Three*D
     !
     IF(NDeg==0) THEN
       FitVal=VectX(NDim)
     ELSE IF(NDeg==1) THEN !!! linear fit
         B=SIGN(ABS(B)+1.D-16,B) !!! to handle perfectly flat surfaces
         IF(DoNorm) THEN
           FitVal=-A/B
         ELSE
           IF(B<Zero) THEN
             FitVal=LastX-LastY/ABS(B)
           ELSE
             FitVal=-A/B
           ENDIF
         ENDIF
       RETURN
     ELSE IF(NDeg==2) THEN !!!! quadratic fit
       Det=B*B-4.D0*A*C
       X0=-B/TwoC
       Y0=A+B*X0+C*X0*X0
       G0=Zero
       IF(Det<Zero) THEN
         FitVal=X0
       ELSE
         Det=SQRT(Det)
         IF(ABS(TwoC)>1.D-6) THEN
           Q=-Half*(B+SIGN(One,B)*Det)
           X1=Q/C
           X2=A/Q
         ELSE
           IF(ABS(B)>1.D-6) THEN
             IF(B<Zero) THEN
               FitVal=LastX-LastY/ABS(B)
             ELSE
               FitVal=-A/B
             ENDIF
             RETURN
           ELSE
             FitVal=(VectX(NDim)-VectX(NDim-1))/Two
             RETURN
           ENDIF
         ENDIF
         G1=B+TwoC*X1
         G2=B+TwoC*X2
         IF(G1>Zero) THEN !!! going for minimum
           FitVal=X1
         ELSE
           FitVal=X2
         ENDIF
       ENDIF
       RETURN
     ELSE  !!!! cubic fit
       CALL CubicRoots(A,B,C,D,X1,X2,X3)
       G1=B+TwoC*X1+ThreeD*X1*X1
       G2=B+TwoC*X2+ThreeD*X2*X2
       G3=B+TwoC*X3+ThreeD*X3*X3
       D1=ABS(X1-LastX)
       D2=ABS(X2-LastX)
       D3=ABS(X3-LastX)
       !
       IF(G1>Zero) THEN
         G0=G1 ; X0=X1 ; D0=D1
       ELSE IF(G2>Zero) THEN
         G0=G2 ; X0=X2 ; D0=D2
       ELSE IF(G3>Zero) THEN
         G0=G3 ; X0=X3 ; D0=D3
       ELSE
         G0=B+TwoC*LastX+ThreeD*LastX*LastX
         FitVal=LastX-LastY/ABS(G0)
         RETURN
       ENDIF
       IF(G2>Zero.AND.D2<D0) THEN
         G0=G2 ; X0=X2 ; D0=D2
       ENDIF
       IF(G3>Zero.AND.D3<D0) THEN
         G0=G3 ; X0=X3 ; D0=D3
       ENDIF
       FitVal=X0
     ENDIF
     !
     PredGrad=A+B*FitVal+C*FitVal*FitVal+D*FitVal*FitVal*FitVal
   END SUBROUTINE Predict
!
!-------------------------------------------------------------------
!
   SUBROUTINE Chi2Fit(VectX,VectY,RMSErr,VectFit,Work, &
                      MaxDeg,NDeg,Coeffs,Chi2V,EPSIn,DoQFit)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY,RMSErr,Work
     REAL(DOUBLE),DIMENSION(:) :: VectFit
     REAL(DOUBLE),DIMENSION(:) :: Coeffs
     REAL(DOUBLE)              :: YBalance
     REAL(DOUBLE)              :: Chi2V,Chi2V1,EPSIn,EPS,EPS1,X0,F
     INTEGER                   :: MaxDeg,NDeg,NDeg1
     INTEGER                   :: NDim,II,I,L,IErr
     LOGICAL                   :: DoQFit
     !
     NDim=SIZE(VectX)
   ! IF(EPSIn<Zero) THEN
   !   EPS=EPSIn
   ! ELSE
   !   YBalance=DOT_PRODUCT(VectY,RMSErr)
   !   IF(ABS(YBalance)<1.D-4) THEN
   !     EPS=Zero
   !   ELSE
   !     EPS=-One
   !   ENDIF
   ! ENDIF
     EPS=-One
     IF(NDim==2) EPS=Zero
     !
     DO II=1,2
       CALL POLFIT(NDim,VectX,VectY,RMSErr,MaxDeg,NDeg,EPS, &
                   VectFit,IERR,Work)
       IF(IErr==1.AND.NDeg/=0) THEN
         EXIT
       ELSE
         EPS=Zero
         MaxDeg=1
       ENDIF
     ENDDO
     X0=Zero
     L=NDeg
     Coeffs=Zero
     CALL PCOEF(L,X0,Coeffs,Work)
     CALL Chi2Value(VectX,VectY,RMSErr,VectFit,Chi2V)
     IF(IErr/=1) CALL Halt('Error in Polinomial fit '// &
       'in Subroutine Chi2Fit, IErr= '//IntToChar(IErr))
   END SUBROUTINE Chi2Fit
!
!---------------------------------------------------------------------
!
   SUBROUTINE Chi2Value(VectX,VectY,RMSErr,VectFit,Chi2V)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY,RMSErr,VectFit
     REAL(DOUBLE)              :: Chi2V,DI,DI2
     INTEGER                   :: NDim,I,J
     !
     NDim=SIZE(VectX)
     Chi2V=Zero
     DO I=1,NDim
       DI=VectY(I)-VectFit(I)
       DI2=DI*DI
       Chi2V=Chi2V+DI2*RMSErr(I)
     ENDDO
   END SUBROUTINE Chi2Value
!
!---------------------------------------------------------------------
!
   SUBROUTINE GrdConvrgd(GStat,IntCs,Grad)
     TYPE(GOptStat)            :: GStat
     TYPE(INTC)                :: IntCs
     REAL(DOUBLE),DIMENSION(:) :: Grad
     REAL(DOUBLE)              :: Sum
     INTEGER                   :: I,J,NDim,NIntC,NCart
     !
     NDim=SIZE(Grad)
     NIntC=IntCs%N
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
       IF(.NOT.IntCs%Constraint%L(I)) THEN
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
!----------------------------------------------------------------------
!
   SUBROUTINE CleanPastGrads(RefStruct,RefGrad,PBCDim,GOpt,SCRPath,Print2)
     REAL(DOUBLE),DIMENSION(:,:) :: RefStruct,RefGrad
     TYPE(DBL_VECT)   :: Carts,CartGrad,AuxGrad
     TYPE(DBL_RNK2)   :: XYZAux
     INTEGER          :: NCart,NatmsLoc,I,J,PBCDim,NDim
     TYPE(GeomOpt)    :: GOpt
     CHARACTER(LEN=*) :: SCRPath
     LOGICAL          :: Print2
     REAL(DOUBLE)     :: Vect9
     !
     NCart=SIZE(RefStruct,1)
     NDim=SIZE(RefStruct,2)
     NatmsLoc=NCart/3
     CALL NEW(Carts,NCart)
     CALL NEW(AuxGrad,NCart-9)
     CALL NEW(CartGrad,NCart)
     CALL New(XYZAux,(/3,NatmsLoc/))
     !
     DO I=1,NDim
       DO J=1,NCart ; Carts%D(J)=RefStruct(J,I) ; ENDDO
       DO J=1,NCart ; CartGrad%D(J)=RefGrad(J,I) ; ENDDO
       CALL CartRNK1ToCartRNK2(Carts%D,XYZAux%D)
       !
       CALL TotalLattGrad(XYZAux%D,PBCDim,CartGrad%D)
       CALL CleanConstrCart(XYZAux%D,PBCDim,CartGrad%D,GOpt,SCRPath)
       IF(GOpt%TrfCtrl%DoTranslOff) THEN
         DO J=1,NCart-9 ; AuxGrad%D(J)=CartGrad%D(J) ; ENDDO
         CALL TranslsOff(AuxGrad%D,Print2)
         DO J=1,NCart-9 ; CartGrad%D(J)=AuxGrad%D(J) ; ENDDO
       ENDIF
       IF(GOpt%TrfCtrl%DoRotOff) THEN
         CALL RotationsOff(CartGrad%D,Carts%D,Print2,PBCDim)
       ENDIF
       DO J=1,NCart ; RefGrad(J,I)=CartGrad%D(J) ; ENDDO
     ENDDO
     !
     CALL Delete(AuxGrad)
     CALL Delete(CartGrad)
     CALL Delete(Carts)
     CALL Delete(XYZAux)
   END SUBROUTINE CleanPastGrads
!
!-------------------------------------------------------------------
!
   SUBROUTINE TotalLattGrad(XYZAux,PBCDim,CartGrad)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZAux
     REAL(DOUBLE),DIMENSION(:)   :: CartGrad
     REAL(DOUBLE)                :: BoxShape(3,3),InvBoxSh(3,3)
     REAL(DOUBLE)                :: Vect1(3),Vect2(3),LattGrad(9)
     INTEGER                     :: PBCDim,I,J,NatmsLoc,NCart
     INTEGER                     :: ALPHA,BETA,IA
     !
     RETURN ! Suppose that TotalGrad is provided by default
     IF(PBCDim==0) RETURN
     NatmsLoc=SIZE(XYZAux,2)
     NCart=3*NatmsLoc
     !
     DO I=1,3
       BoxShape(1:3,I)=XYZAux(1:3,NatmsLoc-3+I)
     ENDDO
     InvBoxSh=InverseBoxShape(BoxShape,PBCDim)
     DO I=1,9 ; LattGrad(I)=CartGrad(NCart-9+I) ; ENDDO
     !
     DO I=1,NatmsLoc-3
       Vect1=XYZAux(1:3,I)
       CALL DGEMM_NNc(3,3,1,One,Zero,InvBoxSh,Vect1,Vect2)
       J=3*(I-1)
       IA=0
       DO BETA=1,3
         DO ALPHA=1,3
           IA=IA+1
           LattGrad(IA)=LattGrad(IA)+CartGrad(J+ALPHA)*Vect2(BETA)
         ENDDO
       ENDDO
     ENDDO
     DO I=1,9 ; CartGrad(NCart-9+I)=LattGrad(I) ; ENDDO
     !
   END SUBROUTINE TotalLattGrad
!
!-------------------------------------------------------------------
!
END MODULE QUICCAMod
