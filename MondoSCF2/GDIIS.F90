MODULE GDIISMod
!
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
!
IMPLICIT NONE
!
CONTAINS
!
!--------------------------------------------------------------
!
   SUBROUTINE GeoDIIS(XYZ,GConstr,GBackTrf,GHess, &
              GGrdTrf,GTrfCtrl,GCoordCtrl,GDIISCtrl, &
              ConvCrit,HFileIn,iCLONE,iGEO,Print,SCRPath, &
              Displ_O,Grad_O,IntGrad_O,E_O,PWD_O,IntCs_O)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     CHARACTER(LEN=*)            :: HFileIn
     TYPE(Constr)                :: GConstr
     TYPE(BackTrf)               :: GBackTrf
     TYPE(GrdTrf)                :: GGrdTrf 
     TYPE(TrfCtrl)               :: GTrfCtrl
     TYPE(CoordCtrl)             :: GCoordCtrl 
     TYPE(GDIIS)                 :: GDIISCtrl  
     TYPE(GConvCrit)             :: ConvCrit  
     INTEGER                     :: iCLONE,iGEO,ICount
     INTEGER                     :: Print
     CHARACTER(Len=*)            :: SCRPath
     CHARACTER(Len=*),OPTIONAL   :: PWD_O
     INTEGER                     :: I,II,J,JJ,K,L,NCart,NatmsLoc
     INTEGER                     :: IGeom,HDFFileID,IStart
     INTEGER                     :: SRMemory,RefMemory
     INTEGER                     :: CartGradMemory,GDIISMemory
     TYPE(DBL_RNK2)              :: SRStruct,RefGrad
     TYPE(DBL_RNK2)              :: RefStruct,SRDispl
     TYPE(DBL_RNK2)              :: Aux
     TYPE(DBL_VECT)              :: Vect
     REAL(DOUBLE),DIMENSION(:),OPTIONAL :: Displ_O,Grad_O,IntGrad_O
     REAL(DOUBLE),OPTIONAL       :: E_O
     TYPE(INTC),OPTIONAL         :: IntCs_O
     TYPE(Hessian)               :: GHess
     !
     GDIISMemory=MIN(GDIISCtrl%MaxMem-1,iGEO)
     IStart=iGEO-GDIISMemory+1
     !
     HDFFileID=OpenHDF(HFileIn)
     HDF_CurrentID= &
       OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
     !
     ! Get GDIIS memory
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     !
     ! Get GDIIS memory of Cartesian coords and grads
     !
     CALL New(SRStruct,(/NCart,GDIISMemory/))
     CALL New(RefStruct,(/NCart,GDIISMemory/))
     CALL New(RefGrad,(/NCart,GDIISMemory/))
     CALL New(SRDispl,(/NCart,GDIISMemory/))
     !
     CALL New(Vect,NCart)
     CALL New(Aux,(/3,NatmsLoc/))
     DO IGeom=IStart,iGEO
       ICount=IGeom-IStart+1
       CALL Get(Aux,'Displ',Tag_O=TRIM(IntToChar(IGeom)))
       CALL CartRNK2ToCartRNK1(Vect%D,Aux%D)
       DO J=1,NCart ; SRStruct%D(J,ICount)=Vect%D(J) ; ENDDO
       CALL Get(Aux,'Abcartesians',Tag_O=TRIM(IntToChar(IGeom)))
       CALL CartRNK2ToCartRNK1(Vect%D,Aux%D)
       DO J=1,NCart ; RefStruct%D(J,ICount)=Vect%D(J) ; ENDDO
       CALL Get(Aux,'gradients',Tag_O=TRIM(IntToChar(IGeom)))
       CALL CartRNK2ToCartRNK1(Vect%D,Aux%D)
       DO J=1,NCart ; RefGrad%D(J,ICount)=Vect%D(J) ; ENDDO
     ENDDO
     CALL Delete(Aux)
     CALL Delete(Vect)
       SRDispl%D=SRStruct%D-RefStruct%D
     !
     ! Calculate new Cartesian coordinates 
     !
     IF(PRESENT(Displ_O).AND.PRESENT(Grad_O)) THEN
       CALL IntCFit(XYZ,Grad_O,IntGrad_O,Displ_O,RefStruct%D,RefGrad%D,&
                    IntCs_O,SCRPath,PWD_O,GGrdTrf,GCoordCtrl,GTrfCtrl,&
                    ConvCrit,GHess,Print,GDIISCtrl%MaxMem,iGEO+1)
     ELSE
       CALL BasicGDIIS(XYZ,GConstr,Print,RefStruct,RefGrad, &
                       SRStruct,SRDispl)
     ENDIF
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(HDFFileID)
     !
     ! Tidy up
     !
     CALL Delete(RefGrad)
     CALL Delete(RefStruct)
     CALL Delete(SRStruct)
     CALL Delete(SRDispl)
   END SUBROUTINE GeoDIIS
! 
!-------------------------------------------------------------------
!
   SUBROUTINE BasicGDIIS(XYZ,CtrlConstr,Print, &
     RefStruct,RefGrad,SRStruct,SRDispl)
     !
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(Constr)                :: CtrlConstr
     INTEGER                     :: I,II,J,K,L,NCart,NatmsLoc
     INTEGER                     :: DimGDIIS
     INTEGER                     :: GDIISMemory
     TYPE(DBL_RNK2)              :: SRDispl,SRStruct
     TYPE(DBL_RNK2)              :: RefGrad,RefStruct
     TYPE(DBL_VECT)              :: Vect,Coeffs
     REAL(DOUBLE)                :: Sum
     INTEGER                     :: Print
     !
     NCart=SIZE(XYZ,2)
     DimGDIIS=SIZE(RefStruct%D,1)
     GDIISMemory=SIZE(RefStruct%D,2)
     !
     ! Calculate GDIIS coeffs!
     !
     CALL New(Coeffs,GDIISMemory)
     !
!    IF(CtrlConstr%NConstr/=0) THEN
!      IF(Print>=DEBUG_GEOP_MIN) THEN
!        WRITE(*,200) 
!        WRITE(Out,200) 
!      ENDIF
!      CALL CalcGDCoeffs(SRDispl%D,Coeffs%D,Print)
!    ELSE
       IF(Print>=DEBUG_GEOP_MIN) THEN
         WRITE(*,300) 
         WRITE(Out,300) 
       ENDIF
       CALL CalcGDCoeffs(RefGrad%D,Coeffs%D,Print)
!    ENDIF
     200 FORMAT("Doing Geometric DIIS based on Cartesian displacements.")
     300 FORMAT("Doing Geometric DIIS based on Cartesian gradients.")
     !
     ! Calculate new geometry
     !
     CALL XYZSum(XYZ,SRStruct%D,Coeffs%D)
     !
     ! Tidy up
     !
     CALL Delete(Coeffs)
     !
   END SUBROUTINE BasicGDIIS
!
!----------------------------------------------------------------------
!
   SUBROUTINE CalcGDCoeffs(ErrorVects,Coeffs,PrintIn)
     REAL(DOUBLE),DIMENSION(:,:) :: ErrorVects
     REAL(DOUBLE),DIMENSION(:)   :: Coeffs
     REAL(DOUBLE)                :: Sum,CondNum
     INTEGER                     :: I,J,K,L,PrintIn
     INTEGER                     :: GDIISMemory,DimGDIIS
     TYPE(DBL_RNK2)              :: AMat,InvA
     TYPE(DBL_VECT)              :: Vect,Scale
     LOGICAL                     :: Print
     !
     Print=PrintIn>=DEBUG_GEOP_MIN
     GDIISMemory=SIZE(ErrorVects,2)
     DimGDIIS=SIZE(ErrorVects,1)
     CALL New(AMat,(/GDIISMemory,GDIISMemory/))
     !
     ! Calc error matrix
     !
     CALL DGEMM_TNc(GDIISMemory,DimGDIIS,GDIISMemory,One,Zero,  &
                    ErrorVects,ErrorVects,AMat%D)
     !
     ! Equilibrate A
     !
     CALL New(Scale,GDIISMemory)
     DO I=1,GDIISMemory
       Scale%D(I)=One/SQRT(AMat%D(I,I))
     ENDDO
     DO I=1,GDIISMemory
       DO J=1,GDIISMemory
         AMat%D(I,J)=Scale%D(I)*Scale%D(J)*AMat%D(I,J)
       ENDDO
     ENDDO
     !
     ! Calculate SVD inverse of 'A'.     
     !
     CALL New(InvA,(/GDIISMemory,GDIISMemory/))
     CALL DIISInvMat(AMat%D,InvA%D,'Basic')
     !CALL DIISInvMat(AMat%D,InvA%D,'StpD')
     !
     ! Scale inverse back to original system
     !
     DO I=1,GDIISMemory
       DO J=1,GDIISMemory
         InvA%D(I,J)=Scale%D(I)*Scale%D(J)*InvA%D(I,J)
       ENDDO
     ENDDO
     CALL Delete(Scale)
     !
     ! Calculate GDIIS coeffs 
     !
     CALL New(Vect,GDIISMemory)
     Vect%D=One
     CALL DGEMM_NNc(GDIISMemory,GDIISMemory,1,One,Zero,  &
                    InvA%D,Vect%D,Coeffs)
     CALL Delete(Vect)
     ! 
     ! Rescale coeffs to get a sum of One.
     !
     Sum=Zero
     DO I=1,GDIISMemory 
       Sum=Sum+Coeffs(I)  
     ENDDO
     Sum=One/Sum
     Coeffs=Sum*Coeffs
     !
     CALL Delete(AMat)
     CALL Delete(InvA)
   END SUBROUTINE CalcGDCoeffs
!
!-------------------------------------------------------------------
!
   SUBROUTINE XYZSum(XYZ,SRStruct,Coeffs)
     REAL(DOUBLE),DIMENSION(:,:)        :: SRStruct,XYZ
     REAL(DOUBLE),DIMENSION(:)          :: Coeffs
     INTEGER                            :: I,J,GDIISMemory,DimGDIIS
     INTEGER                            :: NCart
     TYPE(DBL_VECT)                     :: Vect,Vect2
     !
     GDIISMemory=SIZE(SRStruct,2)
     DimGDIIS=SIZE(SRStruct,1)
     !
     CALL New(Vect,DimGDIIS)
     Vect%D=Zero
     DO I=1,GDIISMemory
       Vect%D=Vect%D+Coeffs(I)*SRStruct(1:DimGDIIS,I) 
     ENDDO
     CALL CartRNK1ToCartRNK2(Vect%D,XYZ)
     !
     CALL Delete(Vect)
   END SUBROUTINE XYZSum
!
!-------------------------------------------------------------------
!
   SUBROUTINE DIISInvMat(AMat,InvA,Char,Cond_O)
     REAL(DOUBLE),DIMENSION(:,:) :: AMat,InvA
     REAL(DOUBLE)                :: CondNum,EigMax,TolAbs,EigAux   
     REAL(DOUBLE),OPTIONAL       :: Cond_O
     INTEGER                     :: NDim,I,INFO
     TYPE(DBL_RNK2)              :: EigVects,Aux,EigVals
     CHARACTER(LEN=*)            :: Char
     !
     CondNum=1.D-7
     IF(PRESENT(Cond_O)) CondNum=Cond_O
     NDim=SIZE(AMat,1)
     CALL New(EigVects,(/NDim,NDim/))
     CALL New(Aux,(/NDim,NDim/))
     CALL New(EigVals,(/NDim,NDim/))
     EigVals%D=Zero
     CALL SetDSYEVWork(NDim)
       BLKVECT%D=AMat
       CALL DSYEV('V','U',NDim,BLKVECT%D,BIGBLOK,BLKVALS%D, &
       BLKWORK%D,BLKLWORK,INFO)
       IF(INFO/=SUCCEED) &
       CALL Halt('DSYEV hosed in DIISInvMat. INFO='&
                  //TRIM(IntToChar(INFO)))
       DO I=1,NDim ; EigVals%D(I,I)=BLKVALS%D(I) ; ENDDO
       EigVects%D=BLKVECT%D 
     CALL UnSetDSYEVWork()
     !
     EigMax=Zero
     DO I=1,NDim ; EigMax=MAX(EigMax,ABS(EigVals%D(I,I))) ; ENDDO
     TolAbs=CondNum*EigMax
     DO I=1,NDim
       IF(Char=='StpD') THEN
         EigVals%D=Zero
         EigVals%D(NDim,NDim)=One
       ELSE 
         IF(ABS(EigVals%D(I,I))>TolAbs) THEN
           EigVals%D(I,I)=One/EigVals%D(I,I) 
         ELSE
           EigVals%D(I,I)=Zero
         ENDIF
       ENDIF
     ENDDO
     CALL DGEMM_NNc(NDim,NDim,NDim,One,Zero,EigVects%D,&
                    EigVals%D,Aux%D)
     CALL DGEMM_NTc(NDim,NDim,NDim,One,Zero,Aux%D, &
                    EigVects%D,InvA)
     !
     CALL Delete(Aux)
     CALL Delete(EigVals)
     CALL Delete(EigVects)
   END SUBROUTINE DIISInvMat
!
!---------------------------------------------------------------------
!
   SUBROUTINE IntCFit(XYZ,Grad,IntGrad,Displ,RefStruct,RefGrad, &
                      IntCs,SCRPath, &
                      PWDPath,GGrdTrf,GCoordCtrl,GTrfCtrl,ConvCrit, &
                      GHess,Print,MaxMem,iGEO)
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ,RefStruct,RefGrad
     REAL(DOUBLE),DIMENSION(:)    :: Grad,Displ,IntGrad
     CHARACTER(LEN=*)             :: SCRPath,PWDPath
     INTEGER                      :: Print,iGEO,MaxMem
     INTEGER                      :: I,J,NMem,NCart,NatmsLoc,NIntC
     TYPE(INTC)                   :: IntCs
     TYPE(DBL_RNK2)               :: XYZAux,IntCGrads,IntCValues
     TYPE(DBL_VECT)               :: VectC,VectI,VectCG,VectX
     TYPE(GrdTrf)                 :: GGrdTrf
     TYPE(CoordCtrl)              :: GCoordCtrl
     TYPE(TrfCtrl)                :: GTrfCtrl
     TYPE(GConvCrit)              :: ConvCrit
     TYPE(Hessian)                :: GHess
     !
     NMem=SIZE(RefStruct,2)
     IF(NMem/=SIZE(RefGrad,2)) CALL Halt('Dim err in IntCFit')
     NatmsLoc=SIZE(XYZ,2)
     NCart=NatmsLoc*3
     IF(NCart/=SIZE(RefStruct,1)) CALL Halt('Dim err 2 in IntCFit')
     !
     CALL New(XYZAux,(/3,NatmsLoc/))
     CALL New(VectC,NCart)
     CALL New(VectCG,NCart)
     CALL New(VectX,NMem+1)
     !
     NIntC=IntCs%N   
     !
     CALL New(IntCGrads,(/NIntC,NMem+1/))
     CALL New(IntCValues,(/NIntC,NMem+1/))
     CALL New(VectI,NIntC)
     !
     ! Calculate IntC gradients using the latest IntC-s     
     !
     DO I=1,NMem
       DO J=1,NCart
         VectC%D(J)=RefStruct(J,I)
         VectCG%D(J)=RefGrad(J,I)
       ENDDO
       CALL CartRNK1ToCartRNK2(VectC%D,XYZAux%D)
       CALL RefreshBMatInfo(IntCs,XYZAux%D,GTrfCtrl, &
                          GCoordCtrl,Print,SCRPath,.TRUE.)
       CALL CartToInternal(XYZAux%D,IntCs,VectCG%D,VectI%D,&
         GGrdTrf,GCoordCtrl,GTrfCtrl,Print,SCRPath)
       IntCGrads%D(:,I)=VectI%D
       CALL INTCValue(IntCs,XYZAux%D,GCoordCtrl%LinCrit, &
                      GCoordCtrl%TorsLinCrit)
       IntCValues%D(:,I)=IntCs%Value%D
     ENDDO
       CALL RefreshBMatInfo(IntCs,XYZ,GTrfCtrl, &
                          GCoordCtrl,Print,SCRPath,.TRUE.)
       IntCGrads%D(:,NMem+1)=IntGrad
       CALL INTCValue(IntCs,XYZ,GCoordCtrl%LinCrit,GCoordCtrl%TorsLinCrit)
       IntCValues%D(:,NMem+1)=IntCs%Value%D
     !
     ! Calculate new values of internals by fitting
     !
     CALL DisplFit(IntCs,IntCGrads%D,IntCValues%D,GHess,GCoordCtrl,&
                   Displ,PWDPath,SCRPath,NCart,MaxMem,iGEO,.FALSE.)
     !
     CALL Delete(VectX)
     CALL Delete(IntCGrads) 
     CALL Delete(IntCValues) 
     CALL Delete(VectI) 
     CALL Delete(VectC) 
     CALL Delete(VectCG) 
     CALL Delete(XYZAux) 
   END SUBROUTINE IntCFit
!
!---------------------------------------------------------------------
!
   SUBROUTINE LocalWeight(LWeight,Weights,IntCs,NCart,SCRPath,USQ_O)
     REAL(DOUBLE),DIMENSION(:,:) :: LWeight,Weights
     CHARACTER(LEN=*)            :: SCRPath
     TYPE(INT_VECT)              :: IGi,JGi,IGiT1,JGiT1,IGiT2,JGiT2
     INTEGER                     :: NCart,NZ,I,J,K1,K2,L1,L2,NMem
     REAL(DOUBLE)                :: Weight,X1,X2,Y
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
     INTEGER                     :: NZ,NIntC,NCart
     !
     CALL ReadBMATR(ISpB,JSpB,ASpB,TRIM(SCRPath)//'B')
     IF(NIntC/=SIZE(ISpB%I)-1) THEN
       CALL Halt('Dimensionality Error in LocalWeight')
     ENDIF
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
   SUBROUTINE DisplFit(IntCs,IntCGrads,IntCValues,GHess,GCoordCtrl,&
                       Displ,Path,SCRPath,NCart,MaxMem,iGEO,DoMix)
     TYPE(INTC)                 :: IntCs
     REAL(DOUBLE),DIMENSION(:)  :: Displ
     REAL(DOUBLE),DIMENSION(:,:):: IntCGrads,IntCValues
     INTEGER                    :: I,J,NIntC,NDim,iGEO,MaxMem
     INTEGER                    :: NCart,NT
     CHARACTER(LEN=*)           :: Path,SCRPath
     CHARACTER(LEN=DCL)         :: Path2
     TYPE(Hessian)              :: GHess
     TYPE(CoordCtrl)            :: GCoordCtrl
     TYPE(DBL_VECT)             :: DisplT
     TYPE(DBL_RNK2)             :: IntCGradsT,IntCValuesT,FittedHessT
     TYPE(DBL_RNK2)             :: WeightsT,LWeightT,ABCT,ABC1T
     TYPE(DBL_RNK2)             :: VarianceT,RangeT,UMat,USQ
     TYPE(INTC)                 :: IntCsT
     LOGICAL                    :: DoMix
     !
     NIntC=IntCs%N   
     NDim=SIZE(IntCGrads,2)
     IF(DoMix) THEN
       NT=NDim
       CALL New(UMat,(/NIntC,NT/))
       CALL UnitaryTR(IntCGrads,IntCValues,UMat,USQ,NT)
     ELSE
       NT=NIntC
     ENDIF
     Path2=TRIM(Path)//'_'//TRIM(IntToChar(IGEO))
     !
     CALL New(ABCT,(/NT,3/))
     CALL New(ABC1T,(/NT,3/))
     CALL New(DisplT,NT)
     CALL New(FittedHessT,(/NT,NDim/))
     CALL New(WeightsT,(/NT,NDim/))
     CALL New(VarianceT,(/NT,2/))
     CALL New(RangeT,(/NT,2/))
     CALL New(IntCGradsT,(/NT,NDim/))
     CALL New(IntCValuesT,(/NT,NDim/))
     CALL New(LWeightT,(/NT,NDim/))
     CALL New(IntCsT,NT)
     !
     CALL PrepareData(IntCValues,IntCs)
     IntCs%PredGrad%D=Zero
     IF(DoMix) THEN
       ABC1T%D(:,1)=Zero
       ABC1T%D(:,2)=One 
       ABC1T%D(:,3)=Zero
       ABCT%D=ABC1T%D
       IntCsT%Active%L=.TRUE.
       CALL DGEMM_TNc(NT,NIntC,NDim,One,Zero, &
   	     UMat%D,IntCGrads,IntCGradsT%D)
       CALL DGEMM_TNc(NT,NIntC,NDim,One,Zero, &
   	     UMat%D,IntCValues,IntCValuesT%D)
     ELSE
       IntCValuesT%D=IntCValues
       IntCGradsT%D=IntCGrads
       CALL SetEq(IntCs,IntCsT,1,NT,1)
       CALL InitialABC(ABC1T%D,IntCs,GHess)
     ENDIF
     !
     CALL PrepPrimW(WeightsT%D,IntCGradsT%D,IntCsT)
     CALL LQFit(IntCValuesT%D,IntCGradsT%D,WeightsT%D,IntCsT,ABC1T%D, &
                VarianceT%D,RangeT%D,Zero,MaxMem,.FALSE.)
     CALL CalcHessian(FittedHessT%D,IntCValuesT%D,ABC1T%D)
     CALL SecondWeight(WeightsT%D,FittedHessT%D)
     IF(DoMix) THEN
       CALL LocalWeight(LWeightT%D,WeightsT%D,IntCsT,NCart, &
                        SCRPath,USQ_O=USQ%D)
     ELSE
       CALL LocalWeight(LWeightT%D,WeightsT%D,IntCsT,NCart,SCRPath)
     ENDIF
     CALL LQFit(IntCValuesT%D,IntCGradsT%D,LWeightT%D,IntCsT,ABCT%D, &
                VarianceT%D,RangeT%D,Zero,MaxMem,.TRUE.)
     CALL DoPredict(ABCT%D,IntCValuesT%D,IntCGradsT%D,IntCsT, &
                    VarianceT%D,RangeT%D,Path2)
     CALL CleanRange(DisplT%D,RangeT%D,IntCsT%PredVal%D, &
                     IntCValuesT%D(:,NDim),NDim)
     IntCsT%PredVal%D=IntCValuesT%D(:,NDim)+DisplT%D
   ! CALL PrtFitM(IntCValuesT%D,IntCGradsT%D,ABCT%D,IntCsT,Path2)
     IF(DoMix) THEN
       CALL BackToPrims(IntCValues,IntCs,DisplT%D,UMat%D,Displ)
     ELSE
       IntCs%PredVal%D=IntCsT%PredVal%D
       Displ=DisplT%D
     ENDIF
     CALL CleanDispl(IntCValues,IntCs,Displ, &
                     GCoordCtrl%MaxStre,GCoordCtrl%MaxAngle) 
     !
     IF(DoMix) THEN
       CALL Delete(UMat)
       CALL Delete(USQ)
     ENDIF
     CALL Delete(FittedHessT)
     CALL Delete(WeightsT)
     CALL Delete(IntCsT)
     CALL Delete(DisplT)
     CALL Delete(VarianceT)
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
   SUBROUTINE CleanRange(DisplT,RangeT,PredVal,Values,NDim)
     REAL(DOUBLE),DIMENSION(:)   :: DisplT,PredVal,Values
     REAL(DOUBLE),DIMENSION(:,:) :: RangeT
     REAL(DOUBLE)                :: Range,Range1,Range2,Displ
     INTEGER                     :: I,NT,NDim
     !
     NT=SIZE(DisplT)
     DO I=1,NT
       Range1=RangeT(I,1)
       Range2=RangeT(I,2)
       Range=Range2-Range1
     ! IF(NDim==2) Range=Two*Range
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
   SUBROUTINE BackToPrims(IntCValues,IntCs,DisplT,UMat,Displ)
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
   END SUBROUTINE BackToPrims
!
!---------------------------------------------------------------------
!
   SUBROUTINE UnitaryTR(IntCGrads,IntCValues,UMat,USQ,NT)
     REAL(DOUBLE),DIMENSION(:,:) :: IntCGrads,IntCvalues
     TYPE(DBL_RNK2)              :: UMatS,UMatSQ,UMat,UMat2,USQ
     INTEGER                     :: I,J,NDim,NIntC,NT,INFO
     !
     NIntC=SIZE(IntCGrads,1)     
     NDim=SIZE(IntCGrads,2)     
     CALL New(UMatS,(/NDim,NDim/))
     CALL New(UMatSQ,(/NDim,NT/))
     !
     CALL DGEMM_TNc(NDim,NIntC,NDim,One,Zero, &
		     IntCGrads,IntCGrads,UmatS%D)
     CALL SetDSYEVWork(NDim)
       BLKVECT%D=UMatS%D
       CALL DSYEV('V','U',NDim,BLKVECT%D,BIGBLOK,BLKVALS%D, &
       BLKWORK%D,BLKLWORK,INFO)
       IF(INFO/=SUCCEED) &
       CALL Halt('DSYEV hosed in UnitaryTR. INFO='&
                  //TRIM(IntToChar(INFO)))
       I=0
       UMatSQ%D=Zero
       DO J=NDim,1,-1
         IF(BLKVALS%D(J)/BLKVALS%D(NDim)>1.D-7.AND.I<NT) THEN
           I=I+1
           UMatSQ%D(:,I)=BLKVECT%D(:,J)
         ENDIF
         IF(I==NT) EXIT
       ENDDO
     CALL UnSetDSYEVWork()
     !
     CALL DGEMM_NNc(NIntC,NDim,NT,One,Zero, &
		     IntCGrads,UMatSQ%D,UMat%D)
     NT=I
     DO J=1,NT
       UMat%D(:,J)=UMat%D(:,J)/SQRT(DOT_PRODUCT(UMat%D(:,J),UMat%D(:,J)))
     ENDDO
     CALL New(USQ,(/NT,NT/))
     CALL DGEMM_TNc(NT,NIntC,NT,One,Zero,UMat%D,UMat%D,USQ%D)
     !
     CALL New(UMat2,(/NIntC,NT/))
     UMat2%D(:,1:NT)=UMat%D(:,1:NT)
     CALL Delete(UMat)
     CALL New(UMat,(/NIntC,NT/))
     UMat%D=UMat2%D
     CALL Delete(UMat2)
     !
     CALL Delete(UMatS)
     CALL Delete(UMatSQ)
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
     Stre=GHess%Stre
     Bend=GHess%Bend
     LinB=GHess%LinB
     OutP=GHess%OutP
     Tors=GHess%Tors
     DO I=1,IntCs%N
       ABC(I,1)=Zero
       ABC(I,3)=Zero
       IF(IntCs%Def%C(I)(1:4)=='STRE') THEN
         ABC(I,2)=Stre
       ELSE IF(IntCs%Def%C(I)(1:4)=='BEND') THEN
         ABC(I,2)=Bend
       ELSE IF(IntCs%Def%C(I)(1:4)=='LINB') THEN
         ABC(I,2)=LinB
       ELSE IF(IntCs%Def%C(I)(1:4)=='OUTP') THEN
         ABC(I,2)=OutP
       ELSE IF(IntCs%Def%C(I)(1:4)=='TORS') THEN
         ABC(I,2)=Tors
       ELSE
         ABC(I,2)=Stre
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
                    ABC(I,1),ABC(I,2),ABC(I,3), &
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
     !
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
       !
       Displ(I)=IntCs%PredVal%D(I)-IntCValues(I,NDim)
       CALL MapDAngle(IntCs%Def%C(I),IntCValues(I,NDim),Displ(I))
       !
       CALL CtrlDispl(IntCs%Def%C(I),Displ(I),One,MaxStre,MaxAngle)
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
   SUBROUTINE DoPredict(ABC,IntCValues,IntCGrads,IntCs, &
                        Variance,Range,Path)
     REAL(DOUBLE),DIMENSION(:,:) :: ABC,IntCValues,IntCGrads,Variance
     REAL(DOUBLE),DIMENSION(:,:) :: Range
     TYPE(INTC)                  :: IntCs
     INTEGER                     :: NIntC,NDim,I,J
     TYPE(DBL_VECT)              :: VectX,VectY
     REAL(DOUBLE)                :: MaxX,MinX
     CHARACTER(LEN=*)            :: Path
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
       IF(MaxX-MinX<1.D-6) THEN
         IntCs%PredVal%D(I)=Half*(MaxX+MinX)
         CYCLE
       ENDIF
       CALL Predict(ABC(I,1),ABC(I,2),ABC(I,3), &
                    VectX%D,VectY%D,IntCs%PredVal%D(I), &
                    IntCs%PredGrad%D(I),Variance(I,:),Range(I,:))
     !
     ENDDO 
     CALL Delete(VectX)
     CALL Delete(VectY)
   ! CALL PrtPred(IntCs,IntCValues,IntCGrads, &
   !              IntCs%PredVal%D,IntCs%PredGrad%D,Path)
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
         ! X=Weights(J,I)/(Hessian(J,I)**2+1.D-20)      
         ! Weights(J,I)=SQRT(X)
         ENDDO
     ENDDO
   END SUBROUTINE SecondWeight
!
!---------------------------------------------------------------------
!
   SUBROUTINE CalcHessian(Hessian,IntCValues,ABC)
     REAL(DOUBLE),DIMENSION(:,:) :: Hessian,IntCValues,ABC
     REAL(DOUBLE)                :: X
     INTEGER                     :: NIntC,NDim,I,J
     !
     NIntC=SIZE(IntCValues,1)
     NDim=SIZE(IntCValues,2)
     !
     DO I=1,NDim
       DO J=1,NIntC
         X=IntCValues(J,I)
         Hessian(J,I)=ABC(J,2)+Two*ABC(J,3)*X
       ENDDO
     ENDDO
   END SUBROUTINE CalcHessian
!
!---------------------------------------------------------------------
!
   SUBROUTINE LQFit(IntCValues,IntCGrads,Weights,IntCs, &
                    ABC,Variance,Range,EPS,MaxMem,DoReOrd,ABC1_O)
     REAL(DOUBLE),DIMENSION(:,:) :: IntCValues,IntCGrads,Weights,ABC
     REAL(DOUBLE),DIMENSION(:,:),OPTIONAL :: ABC1_O
     REAL(DOUBLE),DIMENSION(:,:) :: Variance
     REAL(DOUBLE)                :: MaxX,MinX,Chi2V,EPS,X
     TYPE(INTC)                  :: IntCs
     INTEGER                     :: NIntC,NDim,NDim2,I,J,NDeg,NDeg1,I1
     INTEGER                     :: IStart
     TYPE(DBL_VECT)              :: RMSErr,VectX,VectY
     INTEGER                     :: MaxDeg,NS,MaxMem
     TYPE(DBL_VECT)              :: Work,VectFit1,VectFitB
     TYPE(INT_VECT)              :: IWork
     REAL(DOUBLE),DIMENSION(:,:) :: Range  
     LOGICAL                     :: DoReOrd
     !
     NIntC=SIZE(IntCValues,1) 
     IntCs%Predgrad%D=Zero
     NDim=SIZE(IntCValues,2) 
     NDeg1=NDim
     CALL New(RMSErr,NDim)
     CALL New(VectX,NDim)
     CALL New(VectY,NDim)
     CALL New(VectFit1,NDim)
     CALL New(VectFitB,NDim)
     !
     MaxDeg=2
     NS=MaxDeg+1
     CALL SetDSYEVWork(NS)
     CALL New(Work,4*NDim*NS+2*NS*NS)
     CALL New(IWork,NDim)
     !
     DO I=1,NIntC
       MaxDeg=2
       IF(.NOT.IntCs%Active%L(I)) CYCLE
       DO J=1,NDim 
         VectX%D(J)=IntCValues(I,J) 
         VectY%D(J)=IntCGrads(I,J) 
         RMSErr%D(J)=Weights(I,J)
       ENDDO
       MinX=MINVAL(VectX%D)
       MaxX=MAXVAL(VectX%D)
       IF(MaxX-MinX<1.D-6) THEN
         Range(I,1)=MinX
         Range(I,2)=MaxX
         CYCLE  
       ENDIF
       !
    !  IF(DoReOrd.AND.NDim>=MaxMem) THEN
    !  IF(DoReOrd) THEN
    !    X=ABS(VectX%D(NDim)-VectX%D(NDim-1))
    !    CALL QTest(VectX%D,VectY%D,RMSErr%D,Work%D,IWork%I,NDim,IStart)
    !    Range(I,2)=MINVAL(VectX%D(IStart:NDim))
    !    Range(I,2)=MAXVAL(VectX%D(IStart:NDim))
    !  ELSE
       ! IStart=MAX(1,NDim-2)
       ! Range(I,1)=MINVAL(VectX%D(IStart:NDim))
       ! Range(I,2)=MAXVAL(VectX%D(IStart:NDim))
         Range(I,1)=MinX
         Range(I,2)=MaxX
         IStart=1
    !  ENDIF
       ! invert to get weights
       DO J=NDim,IStart,-1
         RMSErr%D(J)=One/(RMSErr%D(J)+1.D-20)
       ENDDO
       RMSErr%D=RMSErr%D/SUM(RMSErr%D)
       !
       VectFit1%D=Zero
       IF(PRESENT(ABC1_O)) THEN
         DO J=NDim,IStart,-1
           X=VectX%D(J)
           VectFit1%D(J)=ABC1_O(I,1)+ABC1_O(I,2)*X+ABC1_O(I,3)*X*X
         ENDDO
       ENDIF
       CALL FilterBow(VectX%D(IStart:NDim),VectY%D(IStart:NDim),MaxDeg)
       CALL BasicFit(VectX%D(IStart:NDim),VectY%D(IStart:NDim), &
                     RMSErr%D(IStart:NDim),MaxDeg,Work%D,EPS,Chi2V, &
                     ABC(I,:),Variance(I,:),NDeg,NDim-IStart+1, &
                     VectFitB%D(IStart:NDim),VectFit1%D(IStart:NDim)) 
     ENDDO
     !
     CALL Delete(VectFitB)
     CALL Delete(VectFit1)
     CALL Delete(RMSErr)
     CALL Delete(VectX)
     CALL Delete(VectY)
     CALL Delete(Work)
     CALL Delete(IWork)
     CALL UnSetDSYEVWork()
   END SUBROUTINE LQFit
!
!---------------------------------------------------------------------
!
   SUBROUTINE QTest(VectX,VectY,RMSErr,Work,IWork,NDim,IStart)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY,RMSErr,Work
     INTEGER,DIMENSION(:)      :: IWork
     INTEGER                   :: IStart,NDim,NDim2,NDim3,I,J
     REAL(DOUBLE)              :: Q,QTab,Range
     !
     IStart=1
     IF(NDim<3) RETURN
     CALL ReorderN(RMSErr,IWork,NDim)
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
     Range=Work(1)-Work(NDim)        
     QTab=QTest90(MIN(10,NDim))
     DO I=NDim-1,2,-1
       Q=(Work(I-1)-Work(I))/Range
       IF(ABS(Q)>QTab) THEN
         IStart=I
         EXIT
       ENDIF
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
   SUBROUTINE PrepareData(IntCValues,IntCs)
     REAL(DOUBLE),DIMENSION(:,:) :: IntCValues
     TYPE(INTC)                  :: IntCs
     INTEGER                     :: NIntC,NDim,I,J
     REAL(DOUBLE)                :: Center
     !
     NIntC=SIZE(IntCValues,1)
     NDim=SIZE(IntCValues,2)
     DO I=1,NIntC
       IF(IntCs%Def%C(I)(1:4)=='LINB'.OR. & 
          IntCs%Def%C(I)(1:4)=='OUTP'.OR. & 
          IntCs%Def%C(I)(1:4)=='TORS') THEN
         Center=IntCValues(I,NDim)
         DO J=1,NDim 
           CALL PAngle1(Center,IntCValues(I,J))
         ENDDO 
         CALL PAngle1(Center,IntCs%PredVal%D(I))
       ENDIF
     ENDDO
     !
   END SUBROUTINE PrepareData
!
!------------------------------------------------------------------
!
   SUBROUTINE BasicFit(VectX,VectY,RMSErr,MaxDegIn,Work,EPS,Chi2V, &
                       Params,Variance,NDeg,NDim,VectFitB,VectFit1)  
     ! arrays for output 
     INTEGER       :: MaxDegIn,MaxDeg,NDeg,NDim,I,II,J,NS,NDimNS
     REAL(DOUBLE)  :: VectX(:),VectY(:),RMSErr(:)
     REAL(DOUBLE)  :: Params(:),Variance(:),Work(:)
     REAL(DOUBLE)  :: VectFitB(:),VectFit1(:)
     REAL(DOUBLE)  :: X0,Y0,EPS,Chi2V,Det
     REAL(DOUBLE)  :: LastX,LastX2,LastY,LastY2
     INTEGER       :: I1,I2,I3,I4,I5,I6,I7,I8,I9,IStart
     !
     IF(MaxDeg>2) CALL Halt('Maximum allowed degree of fit is 2, '// &
       ' change definition of Params Array to allow higher order fits.')
     MaxDeg=MAX(MIN(MaxDegIn,NDim-2),1)
     LastX=VectX(NDim)
     LastY=VectY(NDim)
     LastX2=VectX(NDim-1)
     LastY2=VectY(NDim-1)
     !
     IF(NDim>2) THEN
       DO II=1,2
         Params=Zero
         Variance=Zero
         NS=MaxDeg+1
         CALL Chi2Fit(VectX,VectY,RMSErr,Work(1:NDim),VectFitB, &
                      Work(NDim+1:),MaxDeg,NDeg,Params,Chi2V,EPS)
         MaxDeg=NDeg
         IF(MaxDeg==2.AND.ABS(Params(3))>1.D-6) THEN 
           X0=-Params(2)/(Two*Params(3))
           Y0=Params(1)+Params(2)*X0+Params(3)*X0*X0
           Det=Params(2)**2-4.D0*Params(1)*Params(3)
         ELSE
           Det=Params(2)**2
         ENDIF
         IF(Det>Zero) THEN
           EXIT
         ELSE
           MaxDeg=1
         ENDIF
       ENDDO
     ELSE
       MaxDeg=1
       Params(3)=Zero
       Params(2)=(LastY-LastY2)/(LastX-LastX2)
       Params(1)=LastY-Params(2)*LastX
       Chi2V=Zero
       Variance=Zero
       Work(1:2)=VectY(1:2)
     ENDIF
     NDeg=MaxDeg
       NS=MaxDeg+1
       NDimNS=NDim*NS
       I1=1
       I2=NDim   
       I3=I2+NDimNS
       I4=I3+NS*NS 
       I5=I4+NS*NS
       BIGBLOK=NS
       BLKLWORK=MAX(1,3*BIGBLOK)
     ! CALL CalcVariance(Params,Work(I1:I2),VectX,VectY,RMSErr, &
     !                   VectFit1,NDim,MaxDeg,Work(I2+1:I3), &
     !                   Work(I3+1:I4),Work(I4+1:I5),Variance, &
     !                   BLKVECT%D(1:NS,1:NS),BLKVALS%D(1:NS), &
     !                   BLKWORK%D(1:NS))
   END SUBROUTINE BasicFit
!
!---------------------------------------------------------------------
!
   SUBROUTINE CalcVariance(Params,VectFit,VectX,VectY,RMSErr, &
                           VectFit1,NDim,MaxDeg,Design,S,Covar, &
                           Variance,BLKVECT2,BLKVALS2,BLKWORK2)
     INTEGER      :: NDim,MaxDeg,NS
     REAL(DOUBLE) :: VectX(:),VectY(:),RMSErr(:)
     REAL(DOUBLE) :: Params(:),Variance(:)
     REAL(DOUBLE) :: VectFit(:),VectFit1(:)
     REAL(DOUBLE) :: Design(NDim,MaxDeg+1)
     REAL(DOUBLE) :: Covar(MaxDeg+1,MaxDeg+1),S(MaxDeg+1,MaxDeg+1)
     INTEGER      :: I,J,ISign,INFO
     REAL(DOUBLE) :: Vect(3),Jacobian(3),EigMax,TolAbs,X,MeanDev
     REAL(DOUBLE) :: OneO2C2,OneO2C,Det,SDet,OneOSDet
     REAL(DOUBLE) :: BLKVECT2(MaxDeg+1,MaxDeg+1)
     REAL(DOUBLE) :: BLKVALS2(MaxDeg+1),BLKWORK2(MaxDeg+1)
     !
     IF(VectY(NDim)*VectY(NDim-1)<Zero) THEN
       Variance=Zero
       RETURN
     ENDIF
     NS=MaxDeg+1
     IF(NS==3) THEN
       OneO2C2=One/(Two*Params(3)*Params(3)+1.D-20)
       OneO2C=One/(Two*Params(3)+1.D-20)
       Det=Params(2)*Params(2)-4.D0*Params(1)*Params(3)
       IF(Det<Zero) CALL Halt('Det<Zero in CalcVariance.')
       SDet=SQRT(Det)
       OneOSDet=One/(SDet+1.D-10)
     ENDIF
     !
     DO J=1,NDim 
      !X=VectFit(J)**2+1.D-10
       X=ABS(VectFit(J)-VectY(J))+1.D-10
      !X=ABS(0.1D0*VectY(J))+1.D-10
      !X=X*RMSErr(J)
       Design(J,1)=One/X
     ENDDO
     DO J=1,NDim
       DO I=1,NS-1
         Design(J,I+1)=VectX(J)*Design(J,I) 
       ENDDO
     ENDDO    
     !
     CALL DGEMM_TNc(NS,NDim,NS,One,Zero,Design,Design,S)
     !
     BLKVECT2=S
     CALL DSYEV('V','U',NS,BLKVECT2,BIGBLOK, &
                BLKVALS2,BLKWORK2,BLKLWORK,INFO)
     IF(INFO/=SUCCEED) &
     CALL Halt('DSYEV failed in CalcVariance. INFO='&
                //TRIM(IntToChar(INFO)))
     EigMax=Zero
     DO I=1,NS ; EigMax=MAX(EigMax,BLKVALS2(I)) ; ENDDO
     !
   ! IF(EigMax<1.D-5) RETURN
     !
     TolAbs=EigMax*1.D-5 
     Covar=Zero
     DO I=1,NS
       IF(BLKVALS2(I)>TolAbs) THEN
         BLKVALS2(I)=One/BLKVALS2(I) 
       ELSE
         BLKVALS2(I)=Zero
       ENDIF
       Covar(I,I)=BLKVALS2(I)
     ENDDO
     CALL DGEMM_NNc(NS,NS,NS,One,Zero,BLKVECT2,Covar,S)
     CALL DGEMM_NTc(NS,NS,NS,One,Zero,S,BLKVECT2,Covar)
     !
     Variance=Zero
     ISign=1
     DO I=1,2
       IF(NS==2.OR.ABS(Params(3))<1.D-6) THEN
         Jacobian(1)=-One/(Params(2)+1.D-20)
         Jacobian(2)=Params(1)/(Params(2)*Params(2)+1.D-20)
         Jacobian(3)=Zero
       ELSE
         Jacobian(1)=-DBLE(ISign)*OneOSDet
         Jacobian(2)=OneO2C*(-One+DBLE(ISign)*Params(2)*OneOSDet)
         Jacobian(3)=-OneO2C2*(-Params(2)+DBLE(ISign)*SDet) &
                     +OneO2C*(-DBLE(ISign)*OneOSDet)
       ENDIF
       CALL DGEMM_NNc(NS,NS,1,One,Zero,Covar,Jacobian(1:NS),Vect(1:NS))
       Variance(I)=SQRT(DOT_PRODUCT(Jacobian(1:NS),Vect(1:NS)))
       ISign=-1
     ENDDO
   END SUBROUTINE CalcVariance
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
         Conv=One/AngstromsToAu
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
   SUBROUTINE PrtFits(VectX,VectY,Path,A,B,C,FitVal,PredGrad,Def)
     INTEGER                  :: NDim,J
     REAL(DOUBLE),DIMENSION(:):: VectX,VectY
     REAL(DOUBLE)             :: A,B,C,FitVal,PredGrad,Conv,ConvI,ConvI2
     CHARACTER(LEN=*)         :: Path,Def
     !
     NDim=SIZE(VectX)
     IF(Def(1:4)=='STRE') THEN
       Conv=One/AngstromsToAu
       ConvI=AngstromsToAu
     ELSE IF(HasAngle(Def)) THEN
       Conv=180.D0/PI
       ConvI=PI/180.D0
     ELSE
       Conv=One
       ConvI=One
     ENDIF
     ConvI2=ConvI*ConvI
     !
     OPEN(UNIT=91,FILE=TRIM(Path)//'_Data',STATUS='UNKNOWN')
       DO J=1,NDim
         WRITE(91,11) J,Conv*VectX(J),VectY(J),FitVal*Conv,PredGrad
       ENDDO
       11 FORMAT(I3,2X,5F20.8)
     CLOSE(91)
     !
     OPEN(UNIT=91,FILE=TRIM(Path)//'_Params',STATUS='UNKNOWN')
         WRITE(91,12) A,B*ConvI,C*ConvI2
     CLOSE(91)
     !
     OPEN(UNIT=91,FILE=TRIM(Path)//'_Pred',STATUS='UNKNOWN')
         WRITE(91,12) FitVal*Conv,PredGrad
     CLOSE(91)
     12 FORMAT(5F20.8)
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
   SUBROUTINE Predict(A,B,C,VectX,VectY,FitVal,PredGrad,Variance,Range)
     REAL(DOUBLE)              :: A,B,C,FitVal,PredGrad,Variance(2)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY
     REAL(DOUBLE)              :: Range(2)
     REAL(DOUBLE)              :: Det,TwoC,X1,X2,G1,G2,X0,Y0,G0
     REAL(DOUBLE)              :: LastX,LastY
     INTEGER                   :: I,J,NDim
     !
     NDim=SIZE(VectX)
     LastX=VectX(NDim)
     LastY=VectY(NDim)
     TwoC=Two*C
     Det=B*B-4.D0*A*C
     IF(ABS(C)>1.D-6) THEN !!! quadratic fit
       X0=-B/TwoC
       Y0=A+B*X0+C*X0*X0
       G0=Zero
       IF(Det<Zero) THEN
         FitVal=X0+Half*SIGN(Variance(1)+Variance(2),X0-LastX)
       ELSE 
         Det=SQRT(Det)
         X1=(-B+Det)/TwoC
         X2=(-B-Det)/TwoC
         G1=B+TwoC*X1
         G2=B+TwoC*X2
         IF(G1>Zero) THEN !!! going for minimum
           FitVal=X1+SIGN(Variance(1),X1-LastX)
         ELSE
           FitVal=X2+SIGN(Variance(2),X2-LastX)
         ENDIF
       ENDIF
     ELSE                 !!! linear fit 
       IF(ABS(B)>1.D-6) THEN  !!! steep enough slope
         IF(B>Zero) THEN  
           FitVal=(-A)/B
           FitVal=FitVal+SIGN(Variance(2),FitVal-LastX)
         ELSE
           B=-B
           A=LastY-B*LastX
           FitVal=LastX-LastY/B
           I=MAX(1,NDim-2)
           Range(1)=MINVAL(VectX(I:NDim))
           Range(2)=MAXVAL(VectX(I:NDim))
         ! FitVal=(-A)/B
         ! FitVal=FitVal+SIGN(Variance(2),FitVal-LastX)
         ENDIF
       ELSE                   !!! barely any slope
         FitVal=LastX        
       ENDIF
     ENDIF
     PredGrad=A+B*FitVal+C*FitVal*FitVal
   END SUBROUTINE Predict
!
!-------------------------------------------------------------------
!
   SUBROUTINE PredictOld(A,B,C,LastX,LastY,FitVal,PredGrad,Variance)
     REAL(DOUBLE)              :: A,B,C,AA,FitVal,PredGrad,Variance(:)
     REAL(DOUBLE)              :: Det,TwoC,X1,X2,G1,G2,X0,Y0,G0
     REAL(DOUBLE)              :: LastX,LastY
     INTEGER                   :: I,J
     !
     AA=A-PredGrad
     TwoC=Two*C
     Det=B*B-4.D0*AA*C
     IF(ABS(C)>1.D-6) THEN !!! quadratic fit
       X0=-B/TwoC
       Y0=A+B*X0+C*X0*X0
       G0=Zero
      !IF<Zero.OR.ABS(Y0)<ABS(LastY)) THEN
       IF(Det<Zero) THEN
         PredGrad=Zero
         FitVal=X0
       ELSE 
         Det=SQRT(Det)
         X1=(-B-Det)/TwoC
         X2=(-B+Det)/TwoC
         G1=B+TwoC*X1
         G2=B+TwoC*X2
         IF(G1>Zero) THEN !!! going for minimum
           FitVal=X1
         ELSE
           FitVal=X2
         ENDIF
       ENDIF
     ELSE                 !!! linear fit 
       IF(ABS(B)>1.D-6) THEN  !!! steep enough slope
         IF(B>Zero) THEN  
           FitVal=(PredGrad-A)/B
         ELSE
           B=-B
           A=LastY-B*LastX
           FitVal=(PredGrad-A)/B
         ! FitVal=LastX-LastY/B
         ENDIF
       ELSE                   !!! barely any slope
         FitVal=LastX        
       ENDIF
     ENDIF
   END SUBROUTINE PredictOld
!
!---------------------------------------------------------------------
!
   SUBROUTINE Chi2Fit(VectX,VectY,RMSErr,VectFit,VectFitB,Work, &
                      MaxDegIn,NDeg,Coeffs,Chi2V,EPSIn)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY,RMSErr,Work
     REAL(DOUBLE),DIMENSION(:) :: VectFitB,VectFit
     REAL(DOUBLE),DIMENSION(:) :: Coeffs
     REAL(DOUBLE),DIMENSION(3) :: Coeffs1
     REAL(DOUBLE)              :: Chi2V,Chi2V1,EPSIn,EPS,EPS1,X0,F
     INTEGER                   :: MaxDegIn,MaxDeg,NDeg,NDeg1
     INTEGER                   :: NDim,II,I,L,IErr
     !
     NDim=SIZE(VectX)
     IF(MaxDeg>2) CALL Halt('Maximum allowed degree of fit is to big in Chi2Fit')
     EPS=EPSIn
     MaxDeg=MaxDegIn
     DO II=1,2     
       CALL POLFIT(NDim,VectX,VectY,RMSErr,MaxDeg,NDeg,EPS, &
                   VectFit,IERR,Work)
       IF(IErr/=1) THEN
         EPS=Zero
         MaxDeg=1
         CYCLE    
       ELSE
         X0=Zero
         L=NDeg   
         Coeffs=Zero
         CALL PCOEF(L,X0,Coeffs,Work)
         CALL Chi2Value(VectX,VectY,RMSErr,VectFit,Chi2V)
         IF(EPSIn==Zero.AND.II==1.AND.MaxDegIn>1) THEN
           VectFitB=VectFit
           Coeffs1=Coeffs
           Chi2V1=Chi2V
           NDeg1=NDeg
           EPS1=EPS
           EPS=Zero
           MaxDeg=1
           CYCLE
         ENDIF
         EXIT
       ENDIF
     ENDDO
     IF(IErr/=1) CALL Halt('Error in Polinomial fit '// &
       'in Subroutine Chi2Fit, IErr= '//IntToChar(IErr))
     IF(EPSIn==Zero.AND.MaxDegIn>1) THEN
       F=EPS**2/EPS1**2
       II=MIN(NDim-1,30)
       IF(F>FTest95D(II)) THEN
      !IF(Chi2V1<0.90D0*Chi2V) THEN
         VectFit=VectFitB 
         Chi2V=Chi2V1 
         EPS=EPS1
         NDeg=NDeg1
         Coeffs=Coeffs1
       ENDIF
     ENDIF
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
END MODULE GDIISMod
