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
   USE InCoords
   USE PunchHDF
!
IMPLICIT NONE
!
CONTAINS
!
!--------------------------------------------------------------
!
   SUBROUTINE GeoDIIS(XYZ,CConstr,LagrMult,GConstr,GBackTrf, &
              GGrdTrf,GTrfCtrl,GCoordCtrl,GDIISCtrl, &
              HFileIn,iCLONE,iGEO,Print,SCRPath,Displ_O,Grad_O)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     CHARACTER(LEN=*)            :: HFileIn
     TYPE(Constr)                :: GConstr
     TYPE(BackTrf)               :: GBackTrf
     TYPE(GrdTrf)                :: GGrdTrf 
     TYPE(TrfCtrl)               :: GTrfCtrl
     TYPE(CoordCtrl)             :: GCoordCtrl 
     TYPE(GDIIS)                 :: GDIISCtrl  
     INTEGER,DIMENSION(:)        :: CConstr
     INTEGER                     :: iCLONE,iGEO,ICount,NLagr
     INTEGER                     :: Print
     CHARACTER(Len=*)            :: SCRPath
     INTEGER                     :: I,II,J,JJ,K,L,NCart,NatmsLoc,NDim
     INTEGER                     :: IGeom,HDFFileID,IStart
     INTEGER                     :: SRMemory,RefMemory
     INTEGER                     :: CartGradMemory,GDIISMemory
     TYPE(DBL_RNK2)              :: SRStruct,RefGrad
     TYPE(DBL_RNK2)              :: RefStruct,SRDispl
     TYPE(DBL_RNK2)              :: Aux
     TYPE(DBL_VECT)              :: Vect,VectLagr
     REAL(DOUBLE),DIMENSION(:),OPTIONAL :: Displ_O,Grad_O
     REAL(DOUBLE),DIMENSION(:)   :: LagrMult
     !
     GDIISMemory=MIN(6,iGEO)
     IF(iGEO<GDIISCtrl%Init) RETURN
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
     CALL Get(NLagr,'nlagr',Tag_O=TRIM(IntToChar(1)))
     NDim=NCart+NLagr
     !
     ! Get GDIIS memory of Cartesian coords and grads
     !
     CALL New(SRStruct,(/NDim,GDIISMemory/))
     CALL New(RefStruct,(/NDim,GDIISMemory/))
     CALL New(RefGrad,(/NDim,GDIISMemory/))
     CALL New(SRDispl,(/NDim,GDIISMemory/))
     !
     CALL New(Vect,NCart)
     CALL New(VectLagr,NLagr)
     CALL New(Aux,(/3,NatmsLoc/))
     DO IGeom=IStart,iGEO
       ICount=IGeom-IStart+1
       CALL Get(Aux,'Displ',Tag_O=TRIM(IntToChar(IGeom)))
       CALL CartRNK2ToCartRNK1(Vect%D,Aux%D)
       DO J=1,NCart ; SRStruct%D(J,ICount)=Vect%D(J) ; ENDDO
       CALL Get(Aux,'Abcartesians',Tag_O=TRIM(IntToChar(IGeom)))
       CALL CartRNK2ToCartRNK1(Vect%D,Aux%D)
       DO J=1,NCart ; RefStruct%D(J,ICount)=Vect%D(J) ; ENDDO
       CALL Get(Vect,'grade',Tag_O=TRIM(IntToChar(IGeom)))
       DO J=1,NCart ; RefGrad%D(J,ICount)=Vect%D(J) ; ENDDO
       DO I=1,NatmsLoc
         IF(CConstr(I)==2) THEN
           K=3*(I-1)
           DO J=K+1,K+3; RefGrad%D(J,ICount)=Zero ; ENDDO
         ENDIF
       ENDDO
       IF(NLagr/=0) THEN
         CALL Get(VectLagr,'LagrDispl',Tag_O=TRIM(IntToChar(IGeom)))
         DO J=1,NLagr ; SRStruct%D(NCart+J,ICount)=VectLagr%D(J) ; ENDDO
         CALL Get(VectLagr,'LagrMult',Tag_O=TRIM(IntToChar(IGeom)))
         DO J=1,NLagr ; RefStruct%D(NCart+J,ICount)=VectLagr%D(J) ; ENDDO
         CALL Get(VectLagr,'GradMult',Tag_O=TRIM(IntToChar(IGeom)))
         DO J=1,NLagr ; RefGrad%D(NCart+J,ICount)=VectLagr%D(J) ; ENDDO
       ENDIF
     ENDDO
     CALL Delete(Aux)
     CALL Delete(Vect)
     CALL Delete(VectLagr)
       SRDispl%D=SRStruct%D-RefStruct%D
     !
     ! Calculate new Cartesian coordinates and Lagr. multipliers
     !
     IF(PRESENT(Displ_O).AND.PRESENT(Grad_O)) THEN
       CALL IntCFit(XYZ,Grad_O,RefStruct%D,RefGrad%D,SCRPath,Displ_O, &
                    GGrdTrf,GCoordCtrl,GTrfCtrl,Print)
     ELSE
       CALL BasicGDIIS(XYZ,GConstr,Print,RefStruct,RefGrad, &
                       SRStruct,SRDispl,LagrMult_O=LagrMult)
       !CALL PIntGDIIS(XYZ,Print,GBackTrf,GTrfCtrl, &
       !  GCoordCtrl,GConstr,SCRPath,RefStruct,RefGrad,SRStruct,SRDispl)
       !
       !CALL DelocGDIIS(SRStruct%D,RefStruct%D,RefGrad%D,SRDispl%D,&
       !  XYZ,Print,SCRPath,GTrfCtrl,GCoordCtrl,GBackTrf,GConstr)
       !
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
   SUBROUTINE DelocGDIIS(SRStruct,RefStruct,RefGrad,SRDispl,XYZ,&
     Print,SCRPath,CtrlTrf,CtrlCoord,CtrlBackTrf,CtrlConstr)
     REAL(DOUBLE),DIMENSION(:,:) :: SRDispl,RefGrad
     REAL(DOUBLE),DIMENSION(:,:) :: SRStruct,RefStruct
     TYPE(INTC)                  :: IntCs
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(DBL_RNK2)              :: XYZRot
     INTEGER                     :: NatmsLoc,NCart,NIntC,GDIISMemory
     INTEGER                     :: I,J,ThreeAt(1:3),Print
     TYPE(BMatr)                 :: B3
     TYPE(INT_VECT)              :: ISpB3,JSpB3
     TYPE(DBL_VECT)              :: ASpB3
     LOGICAL                     :: Linearity
     TYPE(DBL_RNK2)              :: DelocDispl
     TYPE(DBL_RNK2)              :: DelocSR,DelocRef
     TYPE(DBL_VECT)              :: NewDelocs,NewPrims,Displ,Coeffs
     TYPE(Cholesky)              :: CholData3
     TYPE(Constr)                :: CtrlConstr
     TYPE(CoordCtrl)             :: CtrlCoord
     TYPE(BackTrf)               :: CtrlBackTrf
     TYPE(TrfCtrl)               :: CtrlTrf
     CHARACTER(LEN=*)            :: SCRPath
     !
     ! Read internal coord definitions from disk
     !
     CALL ReadIntCs(IntCs,TRIM(SCRPath)//'IntCs')
     NIntC=SIZE(IntCs%Def)
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     GDIISMemory=SIZE(SRDispl,2)
     !
     ! Calculate three-atoms reference system
     ! and rotated structure
     !
     CALL New(XYZRot,(/3,NatmsLoc/))
     XYZRot%D=XYZ
     CALL CALC_XYZRot(XYZRot%D,ThreeAt,Linearity,&
                      CtrlTrf%TranslAt1,CtrlTrf%RotAt2ToX, &
                      CtrlTrf%RotAt3ToXY)
     !
     ! Calculate B matrix for reference system 
     !
     CALL BMatrix(XYZRot%D,NIntC,IntCs,B3, &
                  CtrlCoord%LinCrit,CtrlCoord%TorsLinCrit)
     CALL BtoSpB_1x1(B3,ISpB3,JSpB3,ASpB3)
     !
     ! Calc Gc=B3t*B3 and its Cholesky factors
     !
     CALL CholFact(ISpB3,JSpB3,ASpB3,NCart,NIntC, &
                   CholData3,.FALSE.,Shift_O=0.D0)
     !
     ! Transform Cartesian structures into delocalized internals
     !
     CALL New(DelocRef,(/NCart,GDIISMemory/))
     CALL New(DelocSR,(/NCart,GDIISMemory/))
     CALL New(DelocDispl,(/NCart,GDIISMemory/))
     CALL DelocIntValues(SRStruct,DelocSR%D,B3,CholData3,IntCs, &
                         CtrlCoord%LinCrit,CtrlCoord%TorsLinCrit)
     CALL DelocIntValues(RefStruct,DelocRef%D,B3,CholData3,IntCs, &
                         CtrlCoord%LinCrit,CtrlCoord%TorsLinCrit)
     DelocDispl%D=DelocSR%D-DelocRef%D
     !
     ! Calculate new set of delocalized internals
     !
     CALL New(Coeffs,GDIISMemory)
     CALL CalcGDCoeffs(RefGrad,Coeffs%D,Print)
     !
     CALL New(NewDelocs,NCart) 
     !CALL MixDeloc(Coeffs%D,DelocSR%D,NewDelocs%D) !!! continue here
     !CALL DelocStruct(DelocDispl%D,DelocSR%D,NewDelocs%D)
     !
     ! Turn delocalized internals into primitive ints
     !
     CALL New(NewPrims,NIntC) 
     CALL DelocToPrim(NewDelocs%D,NewPrims%D, &
                      ISpB3,JSpB3,ASpB3,CholData3)
     CALL Delete(B3)
     CALL DeleteBMatInfo(ISpB3,JSpB3,ASpB3,CholData3)
     !
     CALL INTCValue(IntCs,XYZ,CtrlCoord%LinCrit,CtrlCoord%TorsLinCrit)
     CALL New(Displ,NIntC) 
     Displ%D=NewPrims%D-IntCs%Value
     !
     CALL RefreshBMatInfo(IntCs,XYZ,CtrlTrf,CtrlCoord, &
                          Print,SCRPath,.FALSE.)
     CALL InternalToCart(XYZ,IntCs,Displ%D,Print, &
       CtrlBackTrf,CtrlTrf,CtrlCoord,CtrlConstr,SCRPath)
     !
     CALL Delete(Coeffs)
     CALL Delete(Displ)
     CALL Delete(NewPrims)
     CALL Delete(NewDelocs)
     CALL Delete(DelocDispl)
     CALL Delete(DelocSR)
     CALL Delete(DelocRef)
     CALL Delete(XYZRot)
     CALL Delete(IntCs)
   END SUBROUTINE DelocGDIIS
! 
!-------------------------------------------------------------------
!
   SUBROUTINE PrimIntDispl(RefStruct,SRStruct,PrimDispl,IntCs,&
                           LinCrit,Torslincrit)
     TYPE(INTC)     :: IntCs
     REAL(DOUBLE),DIMENSION(:,:) :: RefStruct,SRStruct,PrimDispl
     REAL(DOUBLE)   :: LinCrit,TorsLinCrit
     INTEGER        :: I,J,NCart,NIntC,GDIISMemory,NatmsLoc
     TYPE(DBL_RNK2) :: XYZTmp
     TYPE(DBL_VECT) :: RefValue,SRValue
     !
     NCart=SIZE(RefStruct,1)
     NatmsLoc=NCart/3
     IF(3*NatmsLoc/=NCart) CALL Halt('Dimension error in PrimIntDispl')
     GDIISMemory=SIZE(RefStruct,2)
     NIntC=SIZE(IntCs%Def)
     !
     CALL New(XYZTmp,(/3,NatmsLoc/))
     CALL New(RefValue,NIntC)
     CALL New(SRValue,NIntC)
     DO I=1,GDIISMemory
       CALL CartRNK1ToCartRNK2(RefStruct(1:NCart,I),XYZTmp%D)
       CALL INTCValue(IntCs,XYZTmp%D,LinCrit,TorslinCrit)
       RefValue%D=IntCs%Value        
       CALL CartRNK1ToCartRNK2(SRStruct(1:NCart,I),XYZTmp%D)
       CALL INTCValue(IntCs,XYZTmp%D,LinCrit,TorslinCrit)
       SRValue%D=IntCs%Value        
       PrimDispl(1:NIntC,I)=SRValue%D-RefValue%D
     ENDDO  
     CALL Delete(SRValue)
     CALL Delete(RefValue)
     CALL Delete(XYZTmp)
   END SUBROUTINE PrimIntDispl
! 
!-------------------------------------------------------------------
!
   SUBROUTINE DelocIntValues(CartStruct,DStruct,B3,CholData3, &
                            IntCs,LinCrit,TorslinCrit)
     TYPE(BMatr)                 :: B3
     TYPE(INT_VECT)              :: ISpB3,JSpB3
     TYPE(DBL_VECT)              :: ASpB3
     REAL(DOUBLE),DIMENSION(:,:) :: CartStruct,DStruct
     REAL(DOUBLE)                :: LinCrit,TorsLinCrit
     INTEGER                     :: NCart,NIntC,GDIISMemory
     INTEGER                     :: I,J,NatmsLoc
     TYPE(Cholesky)              :: CholData3
     TYPE(DBL_VECT)              :: VectInt,VectCart
     TYPE(DBL_RNK2)              :: ActCarts
     TYPE(INTC)                  :: IntCs
     !
     NCart=SIZE(DStruct,1)
     NatmsLoc=NCart/3
     NIntC=SIZE(B3%IB,1)
     GDIISMemory=SIZE(CartStruct,2)
     CALL New(VectInt,NIntC)
     CALL New(VectCart,NCart)
     CALL New(ActCarts,(/3,NatmsLoc/))
     CALL BtoSpB_1x1(B3,ISpB3,JSpB3,ASpB3)
     !
     ! B3^t * PrimDispl
     !
     DO I=1,GDIISMemory
       VectCart%D=CartStruct(1:NCart,I)
       CALL CartRNK1ToCartRNK2(VectCart%D,ActCarts%D)
       CALL INTCValue(IntCs,ActCarts%D,LinCrit,TorsLinCrit)
       VectInt%D=IntCs%Value
       CALL PrimToDeloc(VectInt%D,VectCart%D, &
                        ISpB3,JSpB3,ASpB3,CholData3)
       DStruct(1:NCart,I)=VectCart%D
     ENDDO
     !
     CALL Delete(ActCarts)
     CALL Delete(VectCart)
     CALL Delete(VectInt)
   END SUBROUTINE DelocIntValues
! 
!-------------------------------------------------------------------
!
   SUBROUTINE DelocStruct(DelocDispl,DelocSR,NewDelocs)
     REAL(DOUBLE),DIMENSION(:,:) :: DelocDispl,DelocSR
     REAL(DOUBLE),DIMENSION(:)   :: NewDelocs 
     REAL(DOUBLE)                :: Sum,DotP 
     INTEGER                     :: I,J,NCart,GDIISMemory
     TYPE(DBL_VECT)              :: Vect
     !
     NCart=SIZE(DelocDispl,1)
     GDIISMemory=SIZE(DelocDispl,2)
     CALL New(Vect,GDIISMemory)
     !
     DO I=1,NCart
       Vect%D=DelocDispl(I,1:GDIISMemory)
       DotP=DOT_PRODUCT(Vect%D,Vect%D)
       IF(DotP<1.D-3) THEN
         NewDelocs(I)=DelocSR(I,GDIISMemory)
         CYCLE
       ENDIF
       !
       Sum=Zero
       DO J=1,GDIISMemory
         Sum=Sum+Vect%D(J)
       ENDDO
       Sum=Sum/DotP**2
       Vect%D=Sum*Vect%D
       !
       ! ensure, that sum of the components of Vect adds up to One.
       !
       Sum=Zero
       DO J=1,GDIISMemory
         Sum=Sum+Vect%D(J)
       ENDDO
       Sum=One/Sum
       Vect%D=Sum*Vect%D
       !
       ! Calculate new value of the delocalized internals
       !
       Sum=Zero
       DO J=1,GDIISMemory
         Sum=Sum+Vect%D(J)*DelocSR(I,J)
       ENDDO
       NewDelocs(I)=Sum
     ENDDO
     !
     CALL Delete(Vect)
   END SUBROUTINE DelocStruct
!
!-------------------------------------------------------------------
!
   SUBROUTINE BasicGDIIS(XYZ,CtrlConstr,Print, &
     RefStruct,RefGrad,SRStruct,SRDispl,LagrMult_O)
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
     REAL(DOUBLE),DIMENSION(:),OPTIONAL :: LagrMult_O
     !
     NCart=SIZE(XYZ,2)
     DimGDIIS=SIZE(RefStruct%D,1)
     GDIISMemory=SIZE(RefStruct%D,2)
     !
     ! Calculate GDIIS coeffs!
     !
     CALL New(Coeffs,GDIISMemory)
     !
!    IF(CtrlConstr%NConstr/=0.AND..NOT.CtrlConstr%DoLagr) THEN
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
     IF(PRESENT(LagrMult_O)) THEN
       CALL XYZSum(XYZ,SRStruct%D,Coeffs%D,LagrMult_O=LagrMult_O)
     ELSE
       CALL XYZSum(XYZ,SRStruct%D,Coeffs%D)
     ENDIF
     !
     ! Tidy up
     !
     CALL Delete(Coeffs)
     !
   END SUBROUTINE BasicGDIIS
!
!----------------------------------------------------------------------
!
   SUBROUTINE PIntGDIIS(XYZ, &
     PrintIn,GBackTrf,GTrfCtrl,GCoordCtrl,GConstr,SCRPath, &
     RefStruct,RefGrad,SRStruct,SRDispl)
     !
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     LOGICAL                     :: Print
     TYPE(BackTrf)               :: GBackTrf
     TYPE(TrfCtrl)               :: GTrfCtrl
     TYPE(CoordCtrl)             :: GCoordCtrl
     TYPE(Constr)                :: GConstr
     CHARACTER(LEN=*)            :: SCRPath
     INTEGER                     :: PrintIn
     INTEGER                     :: I,II,J,K,L,NCart,NatmsLoc,NIntC
     INTEGER                     :: DimGDIIS
     INTEGER                     :: GDIISMemory
     TYPE(DBL_RNK2)              :: AMat,InvA,SRDispl,SRStruct
     TYPE(DBL_RNK2)              :: RefGrad,RefStruct
     TYPE(DBL_RNK2)              :: XYZ2
     TYPE(DBL_RNK2)              :: PrIDispl,PrIRef,PrISR
     TYPE(DBL_VECT)              :: Coeffs,Vect,Scale
     REAL(DOUBLE)                :: Sum
     TYPE(INTC)                  :: IntCs
     TYPE(INT_VECT)              :: Actives,RangeType
     !
     Print=PrintIn>=DEBUG_GEOP_MIN
     IF(Print) THEN
       WRITE(*,200) 
       WRITE(Out,200) 
     ENDIF
     200 FORMAT('Geometric DIIS, coeffs are calculated from displacements of primitive internal displacements.')
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     DimGDIIS=SIZE(RefStruct%D,1)
     GDIISMemory=SIZE(RefStruct%D,2)
     !
     ! Read internal coord definitions from disk
     !
     CALL ReadIntCs(IntCs,TRIM(SCRPath)//'IntCs')
     NIntC=SIZE(IntCs%Def)
     !
     ! Calculate displacements in internal coordinates
     !
     CALL New(PrIDispl,(/NintC,GDIISMemory/))
     CALL New(PrIRef,(/NintC,GDIISMemory/))
     CALL New(PrISR,(/NintC,GDIISMemory/))
     !
     CALL New(XYZ2,(/3,NatmsLoc/))
     CALL New(Vect,DimGDIIS)
     CALL New(Actives,NIntC)
     CALL New(RangeType,NIntC)
     RangeType%I=0
     Actives%I=1
     DO I=1,GDIISMemory
       Vect%D=SRStruct%D(1:DimGDIIS,I)
       CALL CartRNK1ToCartRNK2(Vect%D,XYZ2%D)
       CALL INTCValue(IntCs,XYZ2%D, &
                      GCoordCtrl%LinCrit,GCoordCtrl%TorsLinCrit)
       PrISR%D(1:NIntC,I)=IntCs%Value
       DO J=1,NIntC 
         IF(.NOT.IntCs%Active(J)) Actives%I(J)=0
       ENDDO
       !
       Vect%D=RefStruct%D(1:DimGDIIS,I)
       CALL CartRNK1ToCartRNK2(Vect%D,XYZ2%D)
       CALL INTCValue(IntCs,XYZ2%D, &
                      GCoordCtrl%LinCrit,GCoordCtrl%TorsLinCrit)
       PrIRef%D(1:NIntC,I)=IntCs%Value
       DO J=1,NIntC 
         IF(.NOT.IntCs%Active(J)) Actives%I(J)=0
       ENDDO
       !
       PrIDispl%D(1:NIntC,I)=PrISR%D(1:NIntC,I)-PrIRef%D(1:NIntC,I)
     ENDDO
     CALL Delete(Vect)
     CALL Delete(XYZ2)
     !
     ! Set types of angle ranges
     !
     CALL SetRangeType(RangeType%I,IntCs,PrIRef%D,PrISR%D)
     !
     ! Convert angle-displacements into degrees and check actives.
     !
     Sum=180.D0/PI
     DO I=1,NIntC
       IF(HasAngle(IntCs%Def(I))) THEN
         PrIDispl%D(I,:)=Sum*PrIDispl%D(I,:)
       ENDIF
       IF(Actives%I(I)==0) PrIDispl%D(I,:)=Zero
     ENDDO
     !
     ! Calculate GDIIS coeffs from internal coord displacements!
     !
     CALL New(Coeffs,GDIISMemory)
     !CALL CalcGDCoeffs(PrIDispl%D,Coeffs%D,PrintIn)
     CALL CalcGDCoeffs(RefGrad%D,Coeffs%D,PrintIn)
     !
     ! Calculate new geometry
     !
     !CALL XYZSum(XYZ,SRStruct%D,Coeffs%D)
     CALL IntCSum(XYZ,PrISR%D,RangeType%I, &
                  Coeffs%D,IntCs,PrintIn,GBackTrf, &
                  GTrfCtrl,GCoordCtrl,GConstr,SCRPath,Actives)
     !
     ! Tidy up
     !
     CALL Delete(Actives)
     CALL Delete(RangeType)
     CALL Delete(Coeffs)
     CALL Delete(PrISR)
     CALL Delete(PrIRef)
     CALL Delete(PrIDispl)
     CALL Delete(IntCs)
   END SUBROUTINE PIntGDIIS
!
!-------------------------------------------------------------------
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
   SUBROUTINE XYZSum(XYZ,SRStruct,Coeffs,LagrMult_O)
     REAL(DOUBLE),DIMENSION(:,:)        :: SRStruct,XYZ
     REAL(DOUBLE),DIMENSION(:)          :: Coeffs
     INTEGER                            :: I,J,GDIISMemory,DimGDIIS
     INTEGER                            :: NCart,NLagr
     TYPE(DBL_VECT)                     :: Vect,Vect2
     REAL(DOUBLE),DIMENSION(:),OPTIONAL :: LagrMult_O
     !
     GDIISMemory=SIZE(SRStruct,2)
     DimGDIIS=SIZE(SRStruct,1)
     !
     CALL New(Vect,DimGDIIS)
     Vect%D=Zero
     DO I=1,GDIISMemory
       Vect%D=Vect%D+Coeffs(I)*SRStruct(1:DimGDIIS,I) 
     ENDDO
     IF(PRESENT(LagrMult_O)) THEN
       NLagr=SIZE(LagrMult_O)
       NCart=DimGDIIS-NLagr
       IF(NCart/=3*(NCart/3)) CALL Halt('Dimesion error in XYZSum.')
       CALL New(Vect2,NCart)
       Vect2%D(1:NCart)=Vect%D(1:NCart)
       CALL CartRNK1ToCartRNK2(Vect2%D,XYZ)
       LagrMult_O(1:NLagr)=Vect%D(NCart+1:DimGDIIS)
       CALL Delete(Vect2)
     ELSE
       CALL CartRNK1ToCartRNK2(Vect%D,XYZ)
     ENDIF
     !
     CALL Delete(Vect)
   END SUBROUTINE XYZSum
!
!-------------------------------------------------------------------
!
   SUBROUTINE IntCSum(XYZ,PrISR,RangeType,Coeffs,IntCs,PrintIn, &
       GBackTrf,GTrfCtrl,GCoordCtrl,GConstr,SCRPath,Actives)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ,PrISR
     REAL(DOUBLE),DIMENSION(:)   :: Coeffs
     INTEGER,DIMENSION(:)        :: RangeType
     TYPE(INTC)                  :: IntCs
     LOGICAL                     :: Print
     TYPE(BackTrf)               :: GBackTrf
     TYPE(TrfCtrl)               :: GTrfCtrl
     TYPE(CoordCtrl)             :: GCoordCtrl
     TYPE(Constr)                :: GConstr
     CHARACTER(LEN=*)            :: SCRPath
     TYPE(DBL_VECT)              :: Vect
     INTEGER                     :: I,NIntC,GDIISMemory
     TYPE(INT_VECT)              :: Actives
     INTEGER                     :: PrintIn
     !
     PRint=PrintIn>=DEBUG_GEOP_MIN
     IF(Print) THEN
       WRITE(*,200) 
       WRITE(Out,200) 
     ENDIF
     200 FORMAT('Doing Geometric DIIS based on the linear combination of primitive internal coordinates.')
     NIntC=SIZE(PrISR,1)
     GDIISMemory=SIZE(PrISR,2)
     !
     ! Calculate new value of internal coordinates
     !
     CALL New(Vect,NIntC)
     Vect%D=Zero
     DO I=1,GDIISMemory
       Vect%D=Vect%D+Coeffs(I)*PrISR(1:NIntC,I)
     ENDDO
     !
     ! For inactive internals set the very last value
     !
     DO I=1,NIntC
       IF(Actives%I(I)==0.OR. &
          RangeType(I)==2) Vect%D(I)=PrISR(I,GDIISMemory)
     ENDDO
     CALL RangeBack(IntCs,RangeType,Vect%D)
       CALL INTCValue(IntCs,XYZ, &
                      GCoordCtrl%LinCrit,GCoordCtrl%TorsLinCrit)
       Vect%D=Vect%D-IntCs%Value
       !CALL RedundancyOff(Vect%D,SCRPath,Print)
     IF(Print) CALL PrtIntCoords(IntCs,Vect%D,'predicted change ')
       Vect%D=IntCs%Value+Vect%D
     IF(Print) CALL PrtIntCoords(IntCs,Vect%D,'predicted internals ')
     !
     CALL InternalToCart(XYZ,IntCs,Vect%D, &
       PrintIn,GBackTrf,GTrfCtrl,GCoordCtrl,GConstr,SCRPath)
     CALL Delete(Vect)
   END SUBROUTINE IntCSum
!
!---------------------------------------------------------------------
!
   SUBROUTINE SetRangeType(RangeType,IntCs,PrIRef,PrISR)
     TYPE(INTC)                  :: IntCs
     INTEGER,DIMENSION(:)        :: RangeType
     REAL(DOUBLE),DIMENSION(:,:) :: PrIRef,PrISR
     INTEGER                     :: I,J,K,NIntC,GDIISMemory
     TYPE(DBL_VECT)              :: VectSR,VectRef
     REAL(DOUBLE)                :: MaxAng,MinAng,MaxRange,PIHalf
     REAL(DOUBLE)                :: Sum1,Sum2,TwoPi,Conv
     !
     ! Set RangeType(I) to 0: no 'remapping', default
     !                     1: do 'remapping at 180 degrees of tors/outp 
     !                     2: range too broad, use last geom's value
     Conv=PI/180.D0  
     MaxRange=10.D0*Conv
     PIHalf=PI*Half
     TwoPi=Two*Pi
     NIntC=SIZE(PrISR,1)
     GDIISMemory=SIZE(PrISR,2) 
     !
     CALL New(VectRef,GDIISMemory)
     CALL New(VectSR,GDIISMemory)
     DO I=1,NIntC
       RangeType(I)=0
       IF(HasAngle(IntCs%Def(I))) THEN
         DO J=1,GDIISMemory ; VectRef%D(J)=PrIRef(I,J) ; ENDDO
         DO J=1,GDIISMemory ; VectSR%D(J)=PrISR(I,J) ; ENDDO
         MaxAng=MaxVal(VectRef%D)
         MaxAng=MAX(MaxAng,MaxVal(VectSR%D))
         MinAng=MinVal(VectRef%D)
         MinAng=MIN(MinAng,MinVal(VectSR%D))
         IF(HasTorsOutP(IntCs%Def(I))) THEN
           IF(MaxAng*MinAng<Zero) THEN
             IF(MaxAng>PiHalf.AND.MinAng<-PiHalf) THEN
               IF(TwoPi+MinAng-MaxAng>MaxRange) THEN
                 RangeType(I)=2
               ELSE
                 RangeType(I)=1
                 DO J=1,GDIISMemory
                   Sum1=PrIRef(I,J)
                   Sum2=PrISR(I,J)
                   PrIRef(I,J)=SIGN(PI,Sum1)-Sum1
                   PrISR(I,J)=SIGN(PI,Sum2)-Sum2
                 ENDDO
               ENDIF
             ENDIF
           ENDIF
         ENDIF
         IF(MaxAng*MinAng<Zero) THEN
           IF(MaxAng+MinAng>MaxRange) RangeType(I)=2
         ELSE
           IF(MaxAng-MinAng>MaxRange) RangeType(I)=2
         ENDIF
       ENDIF
     ENDDO 
     !
     CALL Delete(VectRef)
     CALL Delete(VectSR)
   END SUBROUTINE SetRangeType
!
!---------------------------------------------------------------------
!
   SUBROUTINE RangeBack(IntCs,RangeType,Vect)
     TYPE(INTC)                :: IntCs
     INTEGER,DIMENSION(:)      :: RangeType
     REAL(DOUBLE),DIMENSION(:) :: Vect
     REAL(DOUBLE)              :: Sum 
     INTEGER                   :: I,J,NIntC
     !
     NIntC=SIZE(Vect) 
     IF(NIntC/=SIZE(RangeType))CALL Halt('Dimension error in RangeBack')
     DO I=1,NIntC
       IF(RangeType(I)==1) THEN
         Sum=Vect(I)
         Vect(I)=SIGN(PI,Sum)-Sum
       ENDIF  
     ENDDO
   END SUBROUTINE RangeBack
!
!---------------------------------------------------------------------
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
   SUBROUTINE IntCFit(XYZ,Grad,RefStruct,RefGrad,SCRPath,Displ, &
                      GGrdTrf,GCoordCtrl,GTrfCtrl,Print)
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ,RefStruct,RefGrad
     REAL(DOUBLE),DIMENSION(:)    :: Grad,Displ
     CHARACTER(LEN=*)             :: SCRPath
     INTEGER                      :: Print
     INTEGER                      :: I,J,NMem,NCart,NatmsLoc,NIntC
     TYPE(INTC)                   :: IntCs
     TYPE(DBL_RNK2)               :: XYZAux,IntCGrads,IntCValues
     TYPE(DBL_VECT)               :: VectC,VectI,VectCG
     TYPE(GrdTrf)                 :: GGrdTrf
     TYPE(CoordCtrl)              :: GCoordCtrl
     TYPE(TrfCtrl)                :: GTrfCtrl
     !
     NMem=SIZE(RefStruct,2)
     IF(NMem/=SIZE(RefGrad,2)) CALL Halt('Dim err in IntCFit')
     NatmsLoc=SIZE(XYZ,2)
     NCart=NatmsLoc*3
     IF(NCart/=SIZE(RefStruct,1)) CALL Halt('Dim err 2 in IntCFit')
     !
     CALL New(XYZAux,(/3,NatmsLoc/))
     CALL New(VectC,NCart)
     !
     CALL ReadIntCs(IntCs,TRIM(SCRPath)//'IntCs')
     NIntC=SIZE(IntCs%Def)
     !
     CALL New(IntCGrads,(/NIntC,NMem+1/))
     CALL New(IntCValues,(/NIntC,NMem+1/))
     CALL New(VectI,NIntC)
     !
     ! Calculate IntC gradients using the latest IntC-s     
     !
     DO I=1,NMem
       VectC%D=RefStruct(:,I)
       VectCG%D=RefGrad(:,I)
       CALL CartRNK1ToCartRNK2(VectC%D,XYZAux%D)
       CALL CartToInternal(XYZAux%D,IntCs,VectCG%D,VectI%D,&
         GGrdTrf,GCoordCtrl,GTrfCtrl,Print,SCRPath)
       IntCGrads%D(:,I)=VectI%D
       CALL INTCValue(IntCs,XYZAux%D,GCoordCtrl%LinCrit,GCoordCtrl%TorsLinCrit)
       IntCValues%D(:,I)=IntCs%Value
     ENDDO
       CALL CartToInternal(XYZ,IntCs,Grad,VectI%D,&
         GGrdTrf,GCoordCtrl,GTrfCtrl,Print,SCRPath)
       IntCGrads%D(:,NMem+1)=VectI%D
       CALL INTCValue(IntCs,XYZ,GCoordCtrl%LinCrit,GCoordCtrl%TorsLinCrit)
       IntCValues%D(:,NMem+1)=IntCs%Value
     !
     ! Calculate new values of internals by fitting
     !
     CALL DisplFit(IntCs,IntCGrads%D,IntCValues%D,Displ)
     !
     CALL Delete(IntCGrads) 
     CALL Delete(IntCValues) 
     CALL Delete(VectI) 
     CALL Delete(VectC) 
     CALL Delete(XYZAux) 
     CALL Delete(IntCs) 
   END SUBROUTINE IntCFit
!
!---------------------------------------------------------------------
!
   SUBROUTINE DisplFit(IntCs,IntCGrads,IntCValues,Displ)
     TYPE(INTC)                 :: IntCs
     REAL(DOUBLE),DIMENSION(:)  :: Displ
     REAL(DOUBLE),DIMENSION(:,:):: IntCGrads,IntCValues
     INTEGER                    :: I,J,NIntC,NDim
     TYPE(DBL_VECT)             :: FitVal,VectX,VectY
     !
     NIntC=SIZE(IntCs%Def)
     NDim=SIZE(IntCGrads,2)
     CALL New(VectX,NDim)
     CALL New(VectY,NDim)
     CALL New(FitVal,NIntC)
     !
     DO I=1,NIntC
       IF(.NOT.IntCs%Active(I)) THEN
         FitVal%D(I)=Zero
         CYCLE
       ENDIF
       DO J=1,NDim 
         VectX%D=IntCValues(I,J) 
         VectY%D=IntCGrads(I,J) 
       ENDDO
       CALL Reorder(VectY%D,VectX%D) 
       IF(IntCs%Def(I)=='STRE') THEN
         CALL FitInt(VectX%D,VectY%D,FitVal%D(I),'Morse')
       ELSE IF(HasAngle(IntCs%Def(I))) THEN
         CALL FitInt(VectX%D,VectY%D,FitVal%D(I),'Trigon')
       ENDIF
       Displ(I)=FitVal%D(I)-IntCValues(I,NDim)
     ENDDO
     !
     CALL Delete(VectY)
     CALL Delete(VectX)
     CALL Delete(FitVal)
   END SUBROUTINE DisplFit
!
!---------------------------------------------------------------------
!
   SUBROUTINE FitInt(VectX,VectY,FitVal,Char)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY
     REAL(DOUBLE)              :: FitVal
     INTEGER                   :: I,J,NDim
     CHARACTER(LEN=*)          :: Char
     !
     
     !
   END SUBROUTINE FitInt
!
!---------------------------------------------------------------------
!
   FUNCTION MorsePotDer1(Deq,Req,A,R,D2)
     REAL(DOUBLE) :: MorsePotDer1,Deq,Req,A,R,D2
     !
     MorsePotDer1=Two*Deq*A*(One-exp(-A*(R-Req)))
     D2=Two*Deq*A*A*exp(-A*(R-Req))
   END FUNCTION MorsePotDer1
!
!---------------------------------------------------------------------
!
   SUBROUTINE LocateEndPts(VectX,VectY,XMin,YMin,XMax,YMax)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY
     REAL(DOUBLE)              :: XMin,YMin,XMax,YMax
     REAL(DOUBLE)              :: Tangent,H,D,DY,X,Y,Step,DMax
     INTEGER                   :: I,J,NDim
     !
     ! An ordered set of Y-s must be passed in
     !
     DMax=0.3D0
     NDim=SIZE(VectX)
     IF(VectY(1)<Zero.AND.VectY(NDim)>Zero) THEN ! interpolation
       CALL Locate(VectY,NDim,Zero,J)
       XMin=VectX(J)
       YMin=VectY(J)
       XMax=VectX(J+1)
       YMax=VectY(J+1)
     ELSE ! extrapolation
       IF(ABS(VectY(1))<ABS(VectY(NDim))) THEN
         XMin=VectX(1)
         YMin=VectY(1)
         H=VectX(2)-VectX(1)
         D=VectY(2)-VectY(1)
       ELSE
         XMin=VectX(NDim)
         YMin=VectY(NDim)
         H=VectX(NDim)-VectX(NDim-1)
         D=VectY(NDim)-VectY(NDim-1)
       ENDIF
       ! do linear extrapolation for XMax
       IF(ABS(H)>1.D-7) THEN
         Tangent=D/H
       ELSE
         Tangent=One
       ENDIF
       Step=YMin/Tangent
       X=XMin
       Y=YMin
       DO I=1,10
         X=X+Step
         CALL RatInt(VectX,VectY,NDim,X,Y,DY)
         IF(Y*YMin<Zero) EXIT
       ENDDO
       IF(Y*YMin<Zero) THEN
         XMax=X
         YMax=Y
       ELSE
         XMax=XMin+MIN(SIGN(DMax,(X-XMin)),(X-XMin))
         CALL RatInt(VectX,VectY,NDim,XMax,YMax,DY)
       ENDIF
     ENDIF
   END SUBROUTINE LocateEndPts
!
!---------------------------------------------------------------------
!
   SUBROUTINE Interpolate(VectX,VectY,FitVal)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY
     REAL(DOUBLE)              :: FitVal      
     INTEGER                   :: I,J,NDim
     !
     NDim=SIZE(VectX)
     CALL Locate(VectY,NDim,Zero,J)
   END SUBROUTINE Interpolate
!
!---------------------------------------------------------------------
!
   SUBROUTINE Extrapolate(VectX,VectY,FitVal)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY
     REAL(DOUBLE)              :: FitVal      
     INTEGER                   :: I,J,NDim
   END SUBROUTINE Extrapolate
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
   SUBROUTINE Reorder(VectX,VectY) 
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY
     REAL(DOUBLE)              :: X1,X2,Y1,Y2
     INTEGER                   :: I,J,NDim
     !
     NDim=SIZE(VectX)
     DO I=1,NDim-1
       DO J=1,NDim-I
         X1=VectX(J)
         X2=VectX(J+1)
         Y1=VectY(J)
         Y2=VectY(J+1)
         IF(X1>X2) THEN
           VectX(J)=X2
           VectX(J+1)=X1
           VectY(J)=Y2
           VectY(J+1)=Y1
         ENDIF
       ENDDO
     ENDDO
   END SUBROUTINE Reorder
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
END MODULE GDIISMod
