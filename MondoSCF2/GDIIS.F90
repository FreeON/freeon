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
   SUBROUTINE GeoDIIS(XYZ,GConstr,GBackTrf, &
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
     TYPE(DBL_VECT)              :: Vect,Energy
     REAL(DOUBLE),DIMENSION(:),OPTIONAL :: Displ_O,Grad_O,IntGrad_O
     REAL(DOUBLE),OPTIONAL       :: E_O
     TYPE(INTC),OPTIONAL         :: IntCs_O
     TYPE(INTC)                  :: IntCsOld
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
     CALL New(Energy,GDIISMemory+1)
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
       CALL Get(Energy%D(ICount),'gm_etot',Tag_O=TRIM(IntToChar(IGeom)))
     ENDDO
     CALL Delete(Aux)
     CALL Delete(Vect)
       SRDispl%D=SRStruct%D-RefStruct%D
     !
     ! Calculate new Cartesian coordinates 
     !
     IF(PRESENT(Displ_O).AND.PRESENT(Grad_O)) THEN
       CALL Get(IntCsOld,'IntCs',Tag_O=TRIM(IntToChar(iGEO)))
       CALL MapOldNew(IntCsOld,IntCs_O,NatmsLoc)
       Energy%D(GDIISMemory+1)=E_O
       CALL IntCFit(XYZ,Grad_O,IntGrad_O,Displ_O,RefStruct%D,RefGrad%D,&
                    IntCs_O,Energy%D,SCRPath,PWD_O,GGrdTrf,GCoordCtrl, &
                    GTrfCtrl,ConvCrit,Print,iGEO+1)
       CALL Delete(IntCsOld)
     ELSE
       CALL BasicGDIIS(XYZ,GConstr,Print,RefStruct,RefGrad, &
                       SRStruct,SRDispl)
     ENDIF
     CALL CloseHDFGroup(HDF_CurrentID)
     CALL CloseHDF(HDFFileID)
     !
     ! Tidy up
     !
     CALL Delete(Energy)
     CALL Delete(RefGrad)
     CALL Delete(RefStruct)
     CALL Delete(SRStruct)
     CALL Delete(SRDispl)
   END SUBROUTINE GeoDIIS
! 
!-------------------------------------------------------------------
!
   SUBROUTINE MapOldNew(IntCsOld,IntCs,NAt)
     TYPE(INTC)     :: IntCsOld,IntCs
     INTEGER        :: I,J,K,L,N1,N2,NZ1,NZ2,NAt
     TYPE(INT_VECT) :: I1,J1,A1
     TYPE(INT_VECT) :: I2,J2,A2
     TYPE(INT_VECT) :: I2t,J2t,A2t
     TYPE(INT_VECT) :: IC,JC,CN
     TYPE(INT_RNK2) :: Top1,Top2
     !
     IntCs%PredGrad%D=Zero
     IntCs%PredVal%D=IntCs%Value%D
     N1=IntCsOld%N
     N2=IntCs%N
     CALL New(Top1,(/N1,5/))
     CALL New(Top2,(/N2,5/))
     !
     CALL IntCToTop(IntCsOld,Top1%I)
     CALL IntCToTop(IntCs,Top2%I)
     CALL TopToSp1x1(Top1%I,I1,J1)
     CALL TopToSp1x1(Top2%I,I2,J2)
     NZ1=I1%I(N1+1)-1
     NZ2=I2%I(N2+1)-1
     CALL New(A1,NZ1)
     CALL New(A2,NZ2)
     CALL New(I2t,NAt+1)
     CALL New(J2t,NZ2)
     A1%I=1
     A2%I=1
     CALL ITransPose1x1(I2%I,J2%I,A2%I,N2,NAt,I2t%I,J2t%I,A2%I,'symb')
     CALL IMatMul_1x1(I1%I,J1%I,A1%I,I2t%I,J2t%I,A2%I, &
                      IC,JC,CN,N1,NAt,N2,SymbOnly_O=.FALSE.)
     CALL InvHessEst(IntCs)
     DO I=1,N1
       DO J=IC%I(I),IC%I(I+1)-1 
         K=JC%I(J)
         IF(Top1%I(I,1)==CN%I(J)) THEN
           IF(IntCsOld%Def%C(I)(1:5)==IntCs%Def%C(K)(1:5)) THEN
             IntCs%PredGrad%D(K)=IntCsOld%PredGrad%D(I)
             IntCs%PredVal%D(K)=IntCsOld%PredVal%D(I)
             IntCs%InvHess%D(K)=IntCsOld%InvHess%D(I)
           ENDIF
         ENDIF
       ENDDO
     ENDDO
     !
     CALL Delete(Top1)
     CALL Delete(Top2)
     CALL Delete(I1)
     CALL Delete(I2)
     CALL Delete(I2t)
     CALL Delete(J1)
     CALL Delete(J2)
     CALL Delete(J2t)
     CALL Delete(A1)
     CALL Delete(A2)
     CALL Delete(IC)
     CALL Delete(JC)
     CALL Delete(CN)
   END SUBROUTINE MapOldNew
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
   SUBROUTINE DIISInvMat(AMat,InvA,Char)
     REAL(DOUBLE),DIMENSION(:,:) :: AMat,InvA
     REAL(DOUBLE)                :: CondNum,EigMax,TolAbs,EigAux   
     INTEGER                     :: NDim,I,INFO
     TYPE(DBL_RNK2)              :: EigVects,Aux,EigVals
     CHARACTER(LEN=*)            :: Char
     !
     CondNum=1.D-7
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
                      IntCs,Energy,SCRPath, &
                      PWDPath,GGrdTrf,GCoordCtrl,GTrfCtrl,ConvCrit, &
                      Print,iGEO)
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ,RefStruct,RefGrad
     REAL(DOUBLE),DIMENSION(:)    :: Grad,Displ,IntGrad,Energy
     CHARACTER(LEN=*)             :: SCRPath,PWDPath
     INTEGER                      :: Print,iGEO
     INTEGER                      :: I,J,NMem,NCart,NatmsLoc,NIntC
     TYPE(INTC)                   :: IntCs
     TYPE(DBL_RNK2)               :: XYZAux,IntCGrads,IntCValues,LWeight
     TYPE(DBL_VECT)               :: VectC,VectI,VectCG,RMSErr,VectX
     TYPE(GrdTrf)                 :: GGrdTrf
     TYPE(CoordCtrl)              :: GCoordCtrl
     TYPE(TrfCtrl)                :: GTrfCtrl
     TYPE(GConvCrit)              :: ConvCrit
     REAL(DOUBLE)                 :: EMin
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
     CALL New(LWeight,(/NIntC,NMem+1/))
     CALL New(RMSErr,NMem+1)
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
     ! Calculate energy weighting
     !
     EMin=MINVAL(Energy)
     !
     CALL LocalWeight(LWeight%D,IntCs,NCart, &
                      IntCGrads%D,IntCValues%D,SCRPath)
     !
     ! Calculate new values of internals by fitting
     !
     CALL DisplFit(IntCs,IntCGrads%D,IntCValues%D,RMSErr%D,Displ,&
                   PWDPath,iGEO,Energy,LWeight%D)
     !
     CALL Delete(VectX)
     CALL Delete(RMSErr)
     CALL Delete(LWeight) 
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
   SUBROUTINE LocalWeight(LWeight,IntCs,NCart, &
                          IntCGrads,IntCValues,SCRPath)
     REAL(DOUBLE),DIMENSION(:,:) :: LWeight,IntCGrads,IntCValues
     CHARACTER(LEN=*)            :: SCRPath
     TYPE(INT_VECT)              :: IGi,JGi
     INTEGER                     :: NCart,NZ,I,J,K1,K2,L1,L2,NMem
     REAL(DOUBLE)                :: Weight,X1,X2,Y
     TYPE(INTC)                  :: IntCs
     ! 
     ! Calculate Topology (connectivity) of internal coordinates
     ! 
     NMem=SIZE(LWeight,2)
     !
     CALL GetPattern(SCRPath,IntCs%N,NCart,IGi,JGi)
     !
     ! Calculate RMS gradients for all internal coordinates
     ! touching a specific internal coord at a specific geometry
     ! 
     DO I=1,NMem
       DO J=1,IntCs%N
         NZ=0
         Weight=Zero
         X2=IntCGrads(J,I)
         DO K1=IGi%I(J),IGi%I(J+1)-1
           L1=JGi%I(K1)
           X1=IntCGrads(L1,I)
         ! Weight=Weight+X1*X1
           Weight=Weight+X1*X1*ABS(IntCs%InvHess%D(L1))
           NZ=NZ+1
         ENDDO
         LWeight(J,I)=Weight/DBLE(NZ) 
       ENDDO
     ENDDO
     ! 
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
   SUBROUTINE InvHessEst(IntCs)
     TYPE(INTC)                  :: IntCs
     REAL(DOUBLE)                :: FAngle,FTors,FStre
     INTEGER                     :: I,J
     !
     FStre=2.D0 ! One/0.5D0
     FAngle=5.D0 ! One/0.2D0
     FTors=10.D0 ! One/0.1D0
     !
     DO I=1,IntCs%N
       IF(IntCs%Def%C(I)(1:4)=='BEND'.OR. &
          IntCs%Def%C(I)(1:4)=='LINB'.OR. & 
          IntCs%Def%C(I)(1:4)=='OUTP') THEN 
         Intcs%InvHess%D(I)=FAngle
       ELSE IF(IntCs%Def%C(I)(1:4)=='TORS') THEN
         IntCs%InvHess%D(I)=FTors
       ELSE
         IntCs%InvHess%D(I)=FStre
       ENDIF  
     ENDDO
   END SUBROUTINE InvHessEst
!
!---------------------------------------------------------------------
!
   SUBROUTINE DisplFit(IntCs,IntCGrads,IntCValues,RMSErr,Displ, &
                       Path,iGEO,Energy,LWeight)
     TYPE(INTC)                 :: IntCs
     REAL(DOUBLE),DIMENSION(:)  :: Displ,RMSErr,Energy
     REAL(DOUBLE),DIMENSION(:,:):: IntCGrads,IntCValues,LWeight
     REAL(DOUBLE)               :: A,B,C,Conv,Range
     REAL(DOUBLE)               :: X,Y
     INTEGER                    :: I,J,NIntC,NDim,iGEO
     INTEGER                    :: MaxDim,IStart,IEnd
     TYPE(DBL_VECT)             :: VectX,VectY,Sig
     TYPE(DBL_VECT)             :: VectFit,Work
     CHARACTER(LEN=*)           :: Path
     CHARACTER(LEN=DCL)         :: Path2
     !
     Conv=180.D0/PI
     NIntC=IntCs%N   
     NDim=SIZE(IntCGrads,2)
     CALL New(Sig,NDim)
     CALL New(VectX,NDim)
     CALL New(VectY,NDim)
     CALL New(VectFit,NDim)
     CALL New(Work,6*NDim+3)
     !
     MaxDim=NDim   
     IStart=MAX(NDim-MaxDim+1,1)
     IEnd=NDim
     !
     DO I=1,NIntC
       IF(.NOT.IntCs%Active%L(I)) THEN
         IntCs%PredVal%D(I)=Zero
         CYCLE
       ENDIF
       DO J=1,NDim 
         VectX%D(J)=IntCValues(I,J) 
         VectY%D(J)=IntCGrads(I,J) 
       ENDDO
       !
       ! Reorder angle values for periodic cases
       !
       IF(IntCs%Def%C(I)(1:4)=='LINB'.OR. & 
          IntCs%Def%C(I)(1:4)=='OUTP'.OR. & 
          IntCs%Def%C(I)(1:4)=='TORS') THEN
         CALL PeriodicAngle(VectX%D,IntCs%PredVal%D(I))
       ENDIF
       Range=MAXVAL(VectX%D)-MINVAL(VectX%D)
       !
       ! Errors for fit
       !
       DO J=1,NDim      
         X=1.D-20+LWeight(I,J)
         RMSErr(J)=One/X
        !RMSErr(J)=LWeight(I,J)
       ENDDO
     ! X=SUM(RMSErr)
     ! RMSErr=RMSErr/X
       !
       Path2=TRIM(Path)//'_'//TRIM(IntToChar(IGEO))//'_'// &
            TRIM(IntToChar(I))
       CALL FitXY2(VectX%D(IStart:IEnd),VectY%D(IStart:IEnd), &
                   RMSErr(IStart:IEnd),VectFit%D(IStart:IEnd), &
                   Work%D(IStart:IEnd), &
                   IntCs%PredVal%D(I),IntCs%PredGrad%D(I), &
                   IntCs%Active%L(I),IntCs%InvHess%D(I),  &
                   A,B,C,IntCs%Def%C(I),Path2,MaxDim)
       !
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
       !
       CALL CtrlDispl(IntCs%Def%C(I),Displ(I),One)
       CALL CtrlRange(Displ(I),Range,NDim)
       IntCs%PredVal%D(I)=IntCValues(I,NDim)+Displ(I)
       IF(IntCs%Def%C(I)(1:4)=='STRE') THEN
         CALL ChkStreLim(IntCValues(I,NDim),IntCs%PredVal%D(I),0.3D0)
       ENDIF
       !
       IStart=MAX(IEnd-5,1)
     ! CALL PrtFits(VectX%D(IStart:IEnd),VectY%D(IStart:IEnd), &
     !              Path2,A,B,C,IntCs%PredVal%D(I), &
     !              IntCs%PredGrad%D(I),IntCs%Def%C(I))
     ENDDO
     !
   ! CALL PrtPred(iGEO,NDim,IntCs,IntCValues,IntCGrads, &
   !              IntCs%PredVal%D,IntCs%PredGrad%D,Path2)
     !
     CALL Delete(Work)
     CALL Delete(VectFit)
     CALL Delete(VectY)
     CALL Delete(VectX)
     CALL Delete(Sig)
   END SUBROUTINE DisplFit
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
   SUBROUTINE PrtPred(iGEO,NDim,IntCs,IntCValues,IntCGrads, &
                      FitVal,PredGrad,Path)
     INTEGER                   :: iGEO,I,NDim,NIntC
     TYPE(INTC)                :: IntCs
     REAL(DOUBLE)                :: Conv
     REAL(DOUBLE),DIMENSION(:,:) :: IntCValues,IntCGrads
     REAL(DOUBLE),DIMENSION(:) :: FitVal,PredGrad
     CHARACTER(LEN=*)          :: Path
     !
     NIntC=IntCs%N    
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
   SUBROUTINE PrtFits(VectX,VectY,Path,A,B,C,FitVal,PredGrad,Def,PrtVect_O,PrtBC_O)
     INTEGER                  :: NDim,J
     REAL(DOUBLE),DIMENSION(:):: VectX,VectY
     REAL(DOUBLE)             :: A,B,C,FitVal,PredGrad,Conv,ConvI,ConvI2
     CHARACTER(LEN=*)         :: Path,Def
     LOGICAL,OPTIONAL         :: PrtVect_O,PrtBC_O
     !
     NDim=SIZE(VectX)
     IF(Def(1:4)=='STRE') THEN
       Conv=One/AngstromsToAu
       ConvI=AngstromsToAu
     ELSE IF(HasAngle(Def)) THEN
       Conv=180.D0/PI
       ConvI=PI/180.D0
     ENDIF
     ConvI2=ConvI*ConvI
     !
     IF(PRESENT(PrtVect_O)) THEN
       IF(PrtVect_O) THEN
         OPEN(UNIT=91,FILE=TRIM(Path),STATUS='UNKNOWN')
         DO J=1,NDim
           WRITE(91,12) Conv*VectX(J),VectY(J)
         ENDDO
         12 FORMAT(5F20.8)
         CLOSE(91)
         IF(PRESENT(PrtBC_O)) THEN
           IF(PrtBC_O) THEN
             OPEN(UNIT=91,FILE=TRIM(Path)//'Param',STATUS='UNKNOWN')
               WRITE(91,12) B,C*ConvI
             CLOSE(91)
           ENDIF
         ENDIF
         RETURN
       ENDIF
     ENDIF
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
   END SUBROUTINE PrtFits
!
!---------------------------------------------------------------------
!
   SUBROUTINE fit(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q)
      INTEGER      :: mwt,ndata
      REAL(DOUBLE) :: a,b,chi2,q,siga,sigb,sig(ndata),x(ndata),y(ndata)
      INTEGER      :: i
      REAL(DOUBLE) :: sigdat,ss,st2,sx,sxoss,sy,t,wt
!U    USES gammq
      sx=Zero
      sy=Zero
      st2=Zero
      b=Zero
      if(mwt.ne.0) then
        ss=Zero
        do 11 i=1,ndata
          wt=One/(sig(i)**2)
          ss=ss+wt
          sx=sx+x(i)*wt
          sy=sy+y(i)*wt
11      continue
      else
        do 12 i=1,ndata
          sx=sx+x(i)
          sy=sy+y(i)
12      continue
        ss=float(ndata)
      endif
      sxoss=sx/ss
      if(mwt.ne.0) then
        do 13 i=1,ndata
          t=(x(i)-sxoss)/sig(i)
          st2=st2+t*t
          b=b+t*y(i)/sig(i)
13      continue
      else
        do 14 i=1,ndata
          t=x(i)-sxoss
          st2=st2+t*t
          b=b+t*y(i)
14      continue
      endif
      b=b/st2
      a=(sy-sx*b)/ss
      siga=sqrt((One+sx*sx/(ss*st2))/ss)
      sigb=sqrt(One/st2)
      chi2=Zero
      if(mwt.eq.0) then
        do 15 i=1,ndata
          chi2=chi2+(y(i)-a-b*x(i))**2
15      continue
        q=One
        sigdat=sqrt(chi2/(ndata-2))
        siga=siga*sigdat
        sigb=sigb*sigdat
      else
        do 16 i=1,ndata
          chi2=chi2+((y(i)-a-b*x(i))/sig(i))**2
16      continue
        q=gammq(Half*(ndata-2),half*chi2)
      endif
    END SUBROUTINE Fit
!
!---------------------------------------------------------------------
!
    FUNCTION GammQ(a,x)
      REAL(DOUBLE) :: a,gammq,x
      REAL(DOUBLE) :: gammcf,gamser,gln
!U    USES gcf,gser
      if(x.lt.Zero.or.a.le.Zero) CALL Halt('bad arguments in gammq')
      if(x.lt.a+One)then
        call gser(gamser,a,x,gln)
        gammq=One-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
    END Function GammQ
!
!---------------------------------------------------------------------
!
    SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER      :: ITMAX
      REAL(DOUBLE) :: a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
      INTEGER      :: i
      REAL(DOUBLE) :: an,b,c,d,del,h
!U    USES gammln
      gln=gammln(a)
      b=x+One-a
      c=One/FPMIN
      d=One/b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=One/d
        del=d*c
        h=h*del
        if(abs(del-One).lt.EPS)goto 1
11    continue
      CALL Halt('a too large, ITMAX too small in gcf')
1     gammcf=exp(-x+a*log(x)-gln)*h
    END SUBROUTINE GCF
!
!---------------------------------------------------------------------
!
    SUBROUTINE GSer(gamser,a,x,gln)
      INTEGER      :: ITMAX
      REAL(DOUBLE) :: a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
      INTEGER      :: n
      REAL(DOUBLE) :: ap,del,sum
!U    USES gammln
      gln=gammln(a)
      if(x.le.Zero)then
        if(x.lt.Zero) CALL Halt('x < 0 in gser')
        gamser=Zero
        return
      endif
      ap=a
      sum=One/a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+One
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      CALL Halt('a too large, ITMAX too small in gser')
1     gamser=sum*exp(-x+a*log(x)-gln)
    END SUBROUTINE GSer
!
!---------------------------------------------------------------------
!
    FUNCTION gammln(xx)
      REAL(DOUBLE) :: gammln,xx
      INTEGER      :: j
      REAL(DOUBLE) :: ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,&
        24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2, &
         -.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
    END FUNCTION GammLN
!
!---------------------------------------------------------------------
!
   SUBROUTINE Locate(XX,N,X,J)
     INTEGER      :: J,N
     REAL(DOUBLE) :: X,XX(N) 
     INTEGER      :: JL,JM,JU,I
     !
     JL=0
     JU=N+1
     DO I=1,N
       IF(JU-JL>1) THEN
         JM=(JU+JL)/2
         IF((XX(N)>XX(1)).EQV.(X>XX(JM))) THEN
           JL=JM
         ELSE
           JU=JM
         ENDIF
       ELSE
         EXIT
       ENDIF
     ENDDO
     J=JL
   END SUBROUTINE Locate
!
!---------------------------------------------------------------------
!
   SUBROUTINE RatInt(XA,YA,N,X,Y,DY)
     INTEGER      :: N,NMax
     REAL(DOUBLE) :: DY,X,Y,XA(N),YA(N),Tiny
     PARAMETER(NMax=10,Tiny=1.D-25)
     INTEGER      :: I,M,NS
     REAL(DOUBLE) :: DD,H,HH,T,W,C(NMax),D(NMax),DistMin
     !
     NS=1
     DistMin=1.D-7
     HH=ABS(X-XA(1))
     DO I=1,N
       H=ABS(X-XA(I))
       IF(H<DistMin) THEN
         Y=YA(I)
         DY=Zero 
         RETURN
       ELSE IF(H<HH) THEN
         NS=I
         HH=H
       ENDIF
       C(I)=YA(I)
       D(I)=YA(I)+Tiny
     ENDDO 
     !
     Y=YA(NS)
     NS=NS-1
     DO M=1,N-1
       DO I=1,N-M
         W=C(I+1)-D(I)
         H=XA(I+M)-X
         T=(XA(I)-X)*D(I)/H
         DD=T-C(I+1)
         IF(ABS(DD)<DistMin) CALL Halt('Failure in RatInt')
         DD=W/DD
         D(I)=C(I+1)*DD
         C(I)=T*DD
       ENDDO
       IF(2*NS<N-M) THEN
         DY=D(NS+1)
       ELSE
         DY=D(NS)
         NS=NS-1
       ENDIF
       Y=Y+DY
     ENDDO
     !
   END SUBROUTINE RatInt
!
!---------------------------------------------------------------------
!
   SUBROUTINE MedFit(x,y,ndata,a,b,abdev)
      INTEGER   ::  ndata,NMAX,ndatat
      PARAMETER (NMAX=1000)
      REAL(DOUBLE) :: a,abdev,b,x(ndata),y(ndata)
      REAL(DOUBLE) :: arr(NMAX),xt(NMAX),yt(NMAX),aa,abdevt
      COMMON /arrays/ xt,yt,arr,aa,abdevt,ndatat
      INTEGER      :: j
      REAL(DOUBLE) :: b1,b2,bb,chisq,del,f,f1,f2,sigb
      REAL(DOUBLE) :: sx,sxx,sxy,sy
      !
      sx=Zero
      sy=Zero
      sxy=Zero
      sxx=Zero
      do 11 j=1,ndata
        xt(j)=x(j)
        yt(j)=y(j)
        sx=sx+x(j)
        sy=sy+y(j)
        sxy=sxy+x(j)*y(j)
        sxx=sxx+x(j)**2
11    continue
      ndatat=ndata
      del=ndata*sxx-sx**2
      aa=(sxx*sy-sx*sxy)/del
      bb=(ndata*sxy-sx*sy)/del
      chisq=Zero
      do 12 j=1,ndata
        chisq=chisq+(y(j)-(aa+bb*x(j)))**2
12    continue
      sigb=sqrt(chisq/del)
      b1=bb
      f1=rofunc(b1)
      b2=bb+sign(3.D0*sigb,f1)
      f2=rofunc(b2)
1     if(f1*f2.gt.Zero)then
        bb=2.D0*b2-b1
        b1=b2
        f1=f2
        b2=bb
        f2=rofunc(b2)
        goto 1
      endif
      sigb=0.01D0*sigb
2     if(abs(b2-b1).gt.sigb)then
        bb=0.5D0*(b1+b2)
        if(bb.eq.b1.or.bb.eq.b2)goto 3
        f=rofunc(bb)
        if(f*f1.ge.Zero)then
          f1=f
          b1=bb
        else
          f2=f
          b2=bb
        endif
        goto 2
      endif
3     a=aa
      b=bb
      abdev=abdevt/ndata
      return
   END SUBROUTINE MedFit
!
!---------------------------------------------------------------------
!
   FUNCTION RoFunc(b)
      INTEGER  ::  NMAX
      REAL(DOUBLE) :: RoFunc,b,EPS
      PARAMETER (NMAX=1000,EPS=1.e-7)
      INTEGER      :: j,ndata
      REAL(DOUBLE) :: aa,abdev,d,sum,arr(NMAX),x(NMAX),y(NMAX)
      COMMON /arrays/ x,y,arr,aa,abdev,ndata
      !
      do 11 j=1,ndata
        arr(j)=y(j)-b*x(j)
11    continue
      if (mod(ndata,2).eq.0) then
        j=ndata/2
        aa=0.5D0*(select(j,ndata,arr)+select(j+1,ndata,arr))
      else
        aa=select((ndata+1)/2,ndata,arr)
      endif
      sum=Zero
      abdev=Zero
      do 12 j=1,ndata
        d=y(j)-(b*x(j)+aa)
        abdev=abdev+abs(d)
        if (ABS(y(j))>EPS) d=d/abs(y(j)) 
        if (abs(d).gt.EPS) sum=sum+x(j)*sign(1.0D0,d)
12    continue
      rofunc=sum
      return
   END FUNCTION RoFunc
!
!---------------------------------------------------------------------
!
   FUNCTION select(k,n,arr)
      INTEGER      :: k,n
      REAL(DOUBLE) :: select,arr(n)
      INTEGER      :: i,ir,j,l,mid
      REAL(DOUBLE) :: a,temp
      l=1
      ir=n
1     if(ir-l.le.1)then
        if(ir-l.eq.1)then
          if(arr(ir).lt.arr(l))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
          endif
        endif
        select=arr(k)
        return
      else
        mid=(l+ir)/2
        temp=arr(mid)
        arr(mid)=arr(l+1)
        arr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        if(j.ge.k)ir=j-1
        if(j.le.k)l=i
      endif
      goto 1
   END FUNCTION Select
!
!---------------------------------------------------------------------
!
   SUBROUTINE TestExtraPol(Val,Grad,FitVal)
     REAL(DOUBLE) :: Val,Grad,FitVal
     !
     IF((FitVal-Val)*(-Grad)<Zero) THEN
       FitVal=Val-Grad
     ENDIF
   END SUBROUTINE TestExtraPol
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
   SUBROUTINE DelocEnMass(IntCValues,IntCGrads, &
                          DelocVals,DelocGrads,SCRPath)
     REAL(DOUBLE),DIMENSION(:,:) :: IntCValues,IntCGrads 
     REAL(DOUBLE),DIMENSION(:,:) :: DelocVals,DelocGrads 
     CHARACTER(LEN=*)            :: SCRPath
     TYPE(INT_VECT)              :: ISpB,JSpB 
     TYPE(DBL_VECT)              :: ASpB
     TYPE(DBL_RNK2)              :: UMatr
     INTEGER                     :: I,J,K,NCart,NIntC,NMem
     !
     CALL ReadBMATR(ISpB,JSpB,ASpB,TRIM(SCRPath)//'UMatr',UMatr_O=UMatr)
     NCart=SIZE(Umatr%D,1)
     NIntC=SIZE(IntCValues,1)
     NMem=SIZE(IntCValues,2)
     !
     DO I=1,NMem
       CALL PrimToDeloc(IntCValues(:,I),DelocVals(:,I), &
                        ISpB,JSpB,ASpB,UMatr)
       CALL PrimToDeloc(IntCGrads(:,I),DelocGrads(:,I), &
                        ISpB,JSpB,ASpB,UMatr)
     ENDDO
     !
     CALL Delete(ISpB)
     CALL Delete(JSpB)
     CALL Delete(ASpB)
     CALL Delete(UMatr)
   END SUBROUTINE DelocEnMass
!
!---------------------------------------------------------------------
!
   SUBROUTINE DelocFit(DelocGrads,DelocVals,RMSErr,Displ,SCRPath)
     REAL(DOUBLE),DIMENSION(:,:) :: DelocGrads,DelocVals
     REAL(DOUBLE),DIMENSION(:)   :: Displ,RMSErr
     CHARACTER(LEN=*)            :: SCRPath
     INTEGER                     :: NMem,NIntC,I,J,NDel
     TYPE(DBL_VECT)              :: VectX,VectY,NewDelocVals
     TYPE(INT_VECT)              :: ISpB,JSpB
     TYPE(DBL_VECT)              :: ASpB
     TYPE(DBL_RNK2)              :: UMatr
     REAL(DOUBLE)                :: FitVal,PredGrad,A,B
     !
     NDel=SIZE(DelocGrads,1)
     NMem=SIZE(DelocGrads,2)
     CALL New(VectX,NMem)
     CALL New(VectY,NMem)
     CALL New(NewDelocVals,NDel)
     !
     DO I=1,NDel
       DO J=1,NMem 
         VectX%D(J)=DelocVals(I,J)  
         VectY%D(J)=DelocGrads(I,J)  
       ENDDO
     ! CALL FitXY(VectX%D,VectY%D,RMSErr,FitVal,PredGrad,A,B,1.D-7)
     ! NewDelocVals%D(I)=FitVal
     ENDDO 
     !
     CALL ReadBMATR(ISpB,JSpB,ASpB,TRIM(SCRPath)//'UMatr',UMatr_O=UMatr)
     !
     ! Deloc and PrimInt displacements 
     !
     NewDelocVals%D=NewDelocVals%D-DelocVals(:,NMem)
     CALL PrimToDeloc(Displ,NewDelocVals%D,ISpB,JSpB,ASpB,UMatr, &
                      Char_O='Back')
     !
     CALL Delete(ISpB)
     CALL Delete(JSpB)
     CALL Delete(ASpB)
     CALL Delete(UMatr)
     CALL Delete(NewDelocVals)
     CALL Delete(VectX)
     CALL Delete(VectY)
   END SUBROUTINE DelocFit
!
!---------------------------------------------------------------------
!
   SUBROUTINE FitXY2(VectX,VectY,RMSErr,VectFit,Work, &
                     FitVal,PredGrad,Active, &
                     InvHess,A,B,C,Def,Path,MaxDim)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY,RMSErr,VectFit,Work
     REAL(DOUBLE)              :: FitVal,PredGrad,A,B,C
     REAL(DOUBLE)              :: MaxX,MinX,X1,X2,G1,G2,Det
     REAL(DOUBLE)              :: OldPGrad,OldPVal,InvHess
     CHARACTER(LEN=*)          :: Def,Path
     INTEGER                   :: I,J,NDim,MaxDim,NFit
     REAL(DOUBLE)              :: PredHess
     LOGICAL                   :: Active
     !
     NDim=SIZE(VectX)
     OldPVal=FitVal     !
     OldPGrad=PredGrad
     !
     A=Zero
     B=Zero
     C=Zero
     MaxX=MAXVAL(VectX)
     MinX=MINVAL(VectX)
     IF(ABS(MaxX-MinX)<1.D-6) THEN
       FitVal=SUM(VectX)/DBLE(NDim)
       PredGrad=Zero
       RETURN
     ENDIF
     !
     CALL QuadraticFit(VectX,VectY,RMSErr,VectFit,Work, &
                       Active,OldPGrad,OldPVal,InvHess,A,B,C,Path, &
                       Def,MaxDim,FitVal,PredGrad,PredHess)  
   END SUBROUTINE FitXY2
!
!-------------------------------------------------------------------
!
   SUBROUTINE Predict(A,B,C,NDim,LastX,LastY,Def,MinX,Range, &
                      FitVal,PredGrad,PredHess,InvHess,EPS)
     REAL(DOUBLE)              :: A,B,C,AA,FitVal,PredGrad,PredHess
     REAL(DOUBLE)              :: Det,TwoC,X1,X2,G1,G2,X0,Y0,G0
     REAL(DOUBLE)              :: LastX,LastY,EPS,Range,MinX,InvHess
     INTEGER                   :: I,J,NDim
     CHARACTER(LEN=*)          :: Def
     !
     AA=A-PredGrad
     TwoC=Two*C
     Det=B*B-4.D0*AA*C
     InvHess=One/(B+TwoC*LastX+1.D-10)
     IF(ABS(C)>1.D-6) THEN !!! quadratic fit
       X0=-B/TwoC
       Y0=A+B*X0+C*X0*X0
       G0=Zero
      !IF<Zero.OR.ABS(Y0)<ABS(LastY)) THEN
       IF(Det<Zero) THEN
         PredGrad=Zero
         PredHess=G0
         FitVal=X0
       ELSE 
         Det=SQRT(Det)
         X1=(-B-Det)/TwoC
         X2=(-B+Det)/TwoC
         G1=B+TwoC*X1
         G2=B+TwoC*X2
         IF(G1>Zero) THEN !!! going for minimum
           PredHess=G1
           FitVal=X1
         ELSE
           PredHess=G2
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
         PredHess=B
       ELSE                   !!! barely any slope
         FitVal=LastX        
         PredHess=B
       ENDIF
     ENDIF
   END SUBROUTINE Predict
!
!---------------------------------------------------------------------
!
   SUBROUTINE QuadraticFit(VectX,VectY,RMSErr,VectFit,Work, &
                           Active,OldPGrad,OldPVal,InvHess,A,B,C, &
                           Path,Def,MaxDim,FitVal,PredGrad,PredHess)  
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY,RMSErr,VectFit,Work
     REAL(DOUBLE)              :: A,B,C,Chi2V
     REAL(DOUBLE)              :: FitVal,PredGrad,PredHess,InvHess
     REAL(DOUBLE)              :: LastX,LastY,LastX2,LastY2,LastG
     REAL(DOUBLE)              :: Conv,Det
     REAL(DOUBLE)              :: Range,MinX,MaxX
     REAL(DOUBLE)              :: OldPGrad,OldPVal
     CHARACTER(LEN=*)          :: Path,Def
     INTEGER                   :: II,I,J,NDim,NFit
     INTEGER                   :: MaxDeg,NDeg,IStart,MaxDim
     REAL(DOUBLE)              :: X0,Y0,EPS,CritChi2V,CritChi2V2
     LOGICAL                   :: Active,DoTwoP
     !
     NDim=SIZE(VectX)
     LastX=VectX(NDim)
     LastY=VectY(NDim)
     LastX2=VectX(NDim-1)
     LastY2=VectY(NDim-1)
     MinX=MINVAL(VectX)
     MaxX=MAXVAL(VectX)
     Range=MaxX-MinX
     CritChi2V=0.01D0
     CritChi2V2=0.03D0
     DoTwoP=.FALSE. 
     IF(HasAngle(Def)) THEN
       Conv=180.D0/PI
     ELSE
       Conv=One/AngstromsToAu
     ENDIF
     !
     NFit=NDim
     IF(ABS(LastY)<ABS(LastY2).AND.LastY*LastY2<Zero) THEN
       NFit=2
       DoTwoP=.TRUE.
     ENDIF
     !
     MaxDeg=MIN(2,NFit-2)
     IStart=NDim-NFit+1
     CALL FilterBow(VectX(IStart:NDim),VectY(IStart:NDim), &
                    MinX,MaxX,MaxDeg)
     CALL DetPredGrad(VectX(IStart:NDim),VectY(IStart:NDim), &
                      RMSErr(IStart:NDim),OldPGrad,OldPVal, &
                      Def,DoTwoP,InvHess,PredGrad)
     IF(NFit>2) THEN
       DO II=1,2
         CALL Chi2Fit(VectX(IStart:NDim),VectY(IStart:NDim), &
                      RMSErr(IStart:NDim),VectFit(IStart:NDim),Work, &
                      MaxDeg,NDeg,A,B,C,EPS,Chi2V)
         IF(ABS(C)>1.D-6) THEN 
           X0=-B/(Two*C)
           Y0=A+B*X0+C*X0*X0
           IF(ABS(Y0)<ABS(LastY)) Chi2V=CritChi2V+One
         ENDIF
         Det=B*B-4.D0*(A-PredGrad)*C
         IF(Det>Zero.AND.Chi2V<CritChi2V) THEN
           EXIT
         ELSE
           MaxDeg=1
         ENDIF
       ENDDO
     ELSE
       C=Zero
       B=(LastY-LastY2)/(LastX-LastX2)
       A=LastY-B*LastX
       Chi2V=Zero
       EPS=Zero
     ENDIF
     !
     CALL Predict(A,B,C,NDim,LastX,LastY,Def,MinX,Range, &
                  FitVal,PredGrad,PredHess,InvHess,EPS)
   END SUBROUTINE QuadraticFit
!
!---------------------------------------------------------------------
!
   SUBROUTINE FilterBow(VectX,VectY,MinX,MaxX,MaxDeg)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY
     REAL(DOUBLE)              :: X,X1,X2,MinX,MaxX,Range,Delta
     INTEGER                   :: MaxDeg,I,NDim,II
     !
     NDim=SIZE(VectX)
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
   SUBROUTINE DetPredGrad(VectX,VectY,RMSErr,OldPGrad,OldPVal, &
                          Def,DoTwoP,InvHess,PredGrad)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY,RMSErr
     REAL(DOUBLE)              :: GradMean,PredGrad,InvHess
     REAL(DOUBLE)              :: OldPGrad,OldPVal,X,X1,X2
     REAL(DOUBLE)              :: DGPred,DGAct,Fact
     INTEGER                   :: NDim,I,J
     CHARACTER(LEN=*)          :: Def
     LOGICAL                   :: DoTwoP
     !
     NDim=SIZE(VectY)
     IF(DoTwoP) THEN
       PredGrad=Zero
       RETURN
     ENDIF
     Fact=Zero  
     IF(Def(1:4)=='STRE') THEN
       Fact=0.15D0
     ELSE IF(Def(1:4)=='BEND') THEN
       Fact=0.20D0
     ELSE IF(Def(1:4)=='OUTP'.OR. & 
             Def(1:4)=='LINB') THEN
       Fact=0.00D0
     ELSE IF(Def(1:4)=='TORS') THEN
       Fact=0.25D0
     ENDIF
     GradMean=Zero
     DO I=1,NDim
       GradMean=GradMean+RMSErr(I)*VectY(I)
     ENDDO 
     GradMean=GradMean/SUM(RMSErr)
     PredGrad=-Fact*GradMean
    !PredGrad=MIN(ABS(VectY(NDim)),ABS(VectY(NDim-1)))
    !PredGrad=Fact*SIGN(PredGrad,-VectY(NDim))
   END SUBROUTINE DetPredGrad
!
!---------------------------------------------------------------------
!
   SUBROUTINE Chi2Fit(VectX,VectY,RMSErr,VectFit,Work, &
                      MaxDeg,NDeg,A,B,C,EPS,Chi2V)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY,RMSErr,VectFit,Work
     REAL(DOUBLE),DIMENSION(10):: Coeffs
     REAL(DOUBLE)              :: A,B,C,Chi2V,EPS,X0
     INTEGER                   :: MaxDeg,NDeg,NDim,I,J,K,L,IErr
     !
     EPS= 0.00D0
     NDim=SIZE(VectX)
     CALL POLFIT(NDim,VectX,VectY,RMSErr,MaxDeg,NDeg,EPS, &
                 VectFit,IERR,Work)
     IF(IErr/=1) CALL Halt('Error in Polinomial fit '// &
       'in Subroutine Chi2Fit, IErr= '//IntToChar(IErr))
     !
     X0=Zero
     L=NDeg   
     Coeffs=Zero
     CALL PCOEF (L,X0,Coeffs,Work)
     A=Coeffs(1) ; B=Coeffs(2) ; C=Coeffs(3)
     !
     CALL Chi2Value(VectX,VectY,RMSErr,VectFit,Chi2V)
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
     Chi2V=Chi2V/DBLE(NDim)
   END SUBROUTINE Chi2Value
!
!---------------------------------------------------------------------
!
END MODULE GDIISMod
