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
   SUBROUTINE GeoDIIS(XYZ,GOpt,HFileIn,iCLONE, &
                      iGEO,Print,SCRPath,InitGDIIS)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     CHARACTER(LEN=*)            :: HFileIn
     TYPE(GeomOpt)               :: GOpt
     INTEGER                     :: iCLONE,InitGDIIS,iGEO,ICount
     LOGICAL                     :: Print
     CHARACTER(Len=*)            :: SCRPath
     INTEGER                     :: I,II,J,K,L,NCart,NatmsLoc
     INTEGER                     :: IGeom,HDFFileID,IStart
     INTEGER                     :: SRMemory,RefMemory
     INTEGER                     :: CartGradMemory,GDIISMemory
     TYPE(DBL_RNK2)              :: SRStruct,RefGrad,RefStruct,SRDispl
     TYPE(DBL_RNK2)              :: Aux
     TYPE(DBL_VECT)              :: Vect
     !
     GDIISMemory=MIN(6,iGEO)
     IF(GDIISMemory<InitGDIIS) RETURN
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
       SRStruct%D(:,ICount)=Vect%D
       CALL Get(Aux,'Abcartesians',Tag_O=TRIM(IntToChar(IGeom)))
       CALL CartRNK2ToCartRNK1(Vect%D,Aux%D)
       RefStruct%D(:,ICount)=Vect%D
       CALL Get(Vect,'grade',Tag_O=TRIM(IntToChar(IGeom)))
       RefGrad%D(:,ICount)=Vect%D
     ENDDO
     CALL Delete(Aux)
     CALL Delete(Vect)
       SRDispl%D=SRStruct%D-RefStruct%D
     !
     ! Calculate new Cartesian coordinates 
     !
     CALL BasicGDIIS(XYZ,GOpt%Constr,Print, &
                     RefStruct,RefGrad,SRStruct,SRDispl)
     !CALL PIntGDIIS(XYZ,Print,GOpt%BackTrf,GOpt%TrfCtrl, &
     !  GOpt%CoordCtrl,GOpt%Constr,SCRPath, &
     !  RefStruct,RefGrad,SRStruct,SRDispl)
     !
     !CALL DelocGDIIS(SRStruct%D,RefStruct%D,RefGrad%D,SRDispl%D,&
     !  XYZ,IntCs,Print,SCRPath,GOpt%TrfCtrl,GOpt%CoordCtrl, &
     !  GOpt%BackTrf,GOpt%Constr)
     !
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
   SUBROUTINE DelocGDIIS(SRStruct,RefStruct,RefGrad,SRDispl,XYZ,IntCs,&
     Print,SCRPath,CtrlTrf,CtrlCoord,CtrlBackTrf,CtrlConstr)
     TYPE(INTC)                  :: IntCs
     REAL(DOUBLE),DIMENSION(:,:) :: SRDispl,RefGrad
     REAL(DOUBLE),DIMENSION(:,:) :: SRStruct,RefStruct
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(DBL_RNK2)              :: XYZRot
     INTEGER                     :: NatmsLoc,NCart,NIntC,GDIISMemory
     INTEGER                     :: I,J,ThreeAt(1:3)
     TYPE(BMatr)                 :: B3
     LOGICAL                     :: Linearity,Print
     REAL(DOUBLE)                :: LinCrit
     TYPE(DBL_RNK2)              :: DelocDispl
     TYPE(DBL_RNK2)              :: DelocSR,DelocRef
     TYPE(DBL_VECT)              :: NewDelocs,NewPrims,Displ
     TYPE(Cholesky)              :: CholData3
     TYPE(Constr)                :: CtrlConstr
     TYPE(CoordCtrl)             :: CtrlCoord
     TYPE(BackTrf)               :: CtrlBackTrf
     TYPE(TrfCtrl)               :: CtrlTrf
     CHARACTER(LEN=*)            :: SCRPath
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     NIntC=SIZE(IntCs%Def)
     GDIISMemory=SIZE(SRDispl,2)
     LinCrit=CtrlCoord%LinCrit
     !
     ! Calculate three-atoms reference system
     ! and rotated structure
     !
     CALL New(XYZRot,(/3,NatmsLoc/))
     XYZRot%D=XYZ
     CALL CALC_XYZRot(XYZRot%D,ThreeAt,Linearity)
     !
     ! Calculate B matrix for reference system 
     !
     CALL BMatrix(XYZRot%D,NIntC,IntCs,B3,LinCrit,.FALSE.,ThreeAt)
     !
     ! Calc Gc=B3t*B3 and its Cholesky factors
     !
     CALL CholFact(B3,NCart,CholData3,.FALSE.,.FALSE.,ThreeAt)
     !
     ! Transform prim.ints into delocalized internals
     !
     CALL New(DelocRef,(/NCart,GDIISMemory/))
     CALL New(DelocSR,(/NCart,GDIISMemory/))
     CALL New(DelocDispl,(/NCart,GDIISMemory/))
     CALL DelocIntValues(SRStruct,DelocSR%D,B3,CholData3,IntCs,LinCrit)
     CALL DelocIntValues(RefStruct,DelocRef%D,B3,CholData3,IntCs,LinCrit)
     DelocDispl%D=DelocSR%D-DelocRef%D
     !
     ! Calculate new set of delocalized internals
     !
     CALL New(NewDelocs,NCart) 
     CALL DelocStruct(DelocDispl%D,DelocSR%D,NewDelocs%D)
     !
     ! Turn delocalized internals into primitive ints
     !
     CALL New(NewPrims,NIntC) 
     CALL DelocToPrim(NewDelocs%D,NewPrims%D,B3,CholData3)
     CALL Delete(B3)
     CALL Delete(CholData3)
     !
     CALL INTCValue(IntCs,XYZ,LinCrit)
     CALL New(Displ,NIntC) 
     Displ%D=NewPrims%D-IntCs%Value
     !
     CALL RefreshBMatInfo(IntCs,XYZ,B3,CholData3, &
       CtrlTrf%DoClssTrf,Print,LinCrit,CtrlTrf%ThreeAt,SCRPath)
     CALL Delete(B3)
     CALL Delete(CholData3)
     CALL InternalToCart(XYZ,IntCs,Displ%D,Print, &
       CtrlBackTrf,CtrlTrf,CtrlCoord,CtrlConstr,SCRPath)
     !
     CALL Delete(Displ)
     CALL Delete(NewPrims)
     CALL Delete(NewDelocs)
     CALL Delete(DelocDispl)
     CALL Delete(DelocSR)
     CALL Delete(DelocRef)
     CALL Delete(XYZRot)
   END SUBROUTINE DelocGDIIS
! 
!-------------------------------------------------------------------
!
   SUBROUTINE PrimIntDispl(RefStruct,SRStruct,PrimDispl,IntCs,LinCrit)
     TYPE(INTC)     :: IntCs
     REAL(DOUBLE),DIMENSION(:,:) :: RefStruct,SRStruct,PrimDispl
     REAL(DOUBLE)   :: LinCrit
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
       CALL INTCValue(IntCs,XYZTmp%D,LinCrit)
       RefValue%D=IntCs%Value        
       CALL CartRNK1ToCartRNK2(SRStruct(1:NCart,I),XYZTmp%D)
       CALL INTCValue(IntCs,XYZTmp%D,LinCrit)
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
   SUBROUTINE DelocIntValues(CartStruct,DelocStruct,B3,CholData3, &
                            IntCs,LinCrit)
     TYPE(BMatr) :: B3
     REAL(DOUBLE),DIMENSION(:,:) :: CartStruct,DelocStruct
     REAL(DOUBLE)                :: LinCrit
     INTEGER                     :: NCart,NIntC,GDIISMemory
     INTEGER                     :: I,J,NatmsLoc
     TYPE(Cholesky)              :: CholData3
     TYPE(DBL_VECT)              :: VectInt,VectCart
     TYPE(DBL_RNK2)              :: ActCarts
     TYPE(INTC)                  :: IntCs
     !
     NCart=SIZE(DelocStruct,1)
     NatmsLoc=NCart/3
     NIntC=SIZE(B3%IB,1)
     GDIISMemory=SIZE(CartStruct,2)
     CALL New(VectInt,NIntC)
     CALL New(VectCart,NCart)
     CALL New(ActCarts,(/3,NatmsLoc/))
     !
     ! B3^t * PrimDispl
     !
     DO I=1,GDIISMemory
       VectCart%D=CartStruct(1:NCart,I)
       CALL CartRNK1ToCartRNK2(VectCart%D,ActCarts%D)
       CALL INTCValue(IntCs,ActCarts%D,LinCrit)
       VectInt%D=IntCs%Value
       CALL PrimToDeloc(VectInt%D,VectCart%D,B3,CholData3)
       DelocStruct(1:NCart,I)=VectCart%D
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
     LOGICAL                     :: Print
     !
     NCart=SIZE(XYZ,2)
     DimGDIIS=SIZE(RefStruct%D,1)
     GDIISMemory=SIZE(RefStruct%D,2)
     !
     ! Calculate GDIIS coeffs!
     !
     CALL New(Coeffs,GDIISMemory)
     !
     IF(CtrlConstr%NConstr/=0) THEN
       IF(Print) THEN
         WRITE(*,200) 
         WRITE(Out,200) 
       ENDIF
       CALL CalcGDCoeffs(SRDispl%D,Coeffs%D)
     ELSE
       IF(Print) THEN
         WRITE(*,300) 
         WRITE(Out,300) 
       ENDIF
       CALL CalcGDCoeffs(RefGrad%D,Coeffs%D)
     ENDIF
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
   SUBROUTINE PIntGDIIS(XYZ, &
     Print,GBackTrf,GTrfCtrl,GCoordCtrl,GConstr,SCRPath, &
     RefStruct,RefGrad,SRStruct,SRDispl)
     !
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     LOGICAL                     :: Print
     TYPE(BackTrf)               :: GBackTrf
     TYPE(TrfCtrl)               :: GTrfCtrl
     TYPE(CoordCtrl)             :: GCoordCtrl
     TYPE(Constr)                :: GConstr
     CHARACTER(LEN=*)            :: SCRPath
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
     TYPE(INT_VECT)              :: Actives
     !
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
     Actives%I=1
     DO I=1,GDIISMemory
       Vect%D=SRStruct%D(1:DimGDIIS,I)
       CALL CartRNK1ToCartRNK2(Vect%D,XYZ2%D)
       CALL INTCValue(IntCs,XYZ2%D,GCoordCtrl%LinCrit)
       PrISR%D(1:NIntC,I)=IntCs%Value
       DO J=1,NIntC 
         IF(.NOT.IntCs%Active(J)) Actives%I(J)=0
       ENDDO
       !
       Vect%D=RefStruct%D(1:DimGDIIS,I)
       CALL CartRNK1ToCartRNK2(Vect%D,XYZ2%D)
       CALL INTCValue(IntCs,XYZ2%D,GCoordCtrl%LinCrit)
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
     ! Convert angles into degrees and check actives.
     !
     Sum=180.D0/PI
     DO I=1,NIntC
       IF(IntCs%Def(I)(1:4)/='STRE') THEN
         PrIDispl%D(I,:)=Sum*PrIDispl%D(I,:)
       ENDIF
       IF(Actives%I(I)==0) PrIDispl%D(I,:)=Zero
     ENDDO
     !
     ! Calculate GDIIS coeffs from internal coord displacements!
     !
     CALL New(Coeffs,GDIISMemory)
     CALL CalcGDCoeffs(PrIDispl%D,Coeffs%D)
     !
     ! Calculate new geometry
     !
     !CALL XYZSum(XYZ,SRStruct%D,Coeffs%D)
     CALL IntCSum(XYZ,PrISR%D,Coeffs%D,IntCs,Print,GBackTrf, &
       GTrfCtrl,GCoordCtrl,GConstr,SCRPath,Actives)
     !
     ! Tidy up
     !
     CALL Delete(Actives)
     CALL Delete(Coeffs)
     CALL Delete(PrISR)
     CALL Delete(PrIRef)
     CALL Delete(PrIDispl)
     CALL Delete(IntCs)
   END SUBROUTINE PIntGDIIS
!
!-------------------------------------------------------------------
!
   SUBROUTINE CalcGDCoeffs(ErrorVects,Coeffs)
     REAL(DOUBLE),DIMENSION(:,:) :: ErrorVects
     REAL(DOUBLE),DIMENSION(:)   :: Coeffs
     REAL(DOUBLE)                :: Sum   
     INTEGER                     :: I,J,K,L
     INTEGER                     :: GDIISMemory,DimGDIIS
     TYPE(DBL_RNK2)              :: AMat,InvA
     TYPE(DBL_VECT)              :: Vect,Scale
     !
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
     CALL SetDSYEVWork(GDIISMemory**2)
     CALL FunkOnSqMat(GDIISMemory,Inverse,AMat%D,InvA%D,Unit_O=6)
     CALL UnSetDSYEVWork()
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
     REAL(DOUBLE),DIMENSION(:,:) :: SRStruct,XYZ
     REAL(DOUBLE),DIMENSION(:)   :: Coeffs
     INTEGER                     :: I,J,GDIISMemory,DimGDIIS
     TYPE(DBL_VECT)              :: Vect
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
   SUBROUTINE IntCSum(XYZ,PrISR,Coeffs,IntCs,Print,GBackTrf, &
       GTrfCtrl,GCoordCtrl,GConstr,SCRPath,Actives)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ,PrISR
     REAL(DOUBLE),DIMENSION(:)   :: Coeffs
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
     !
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
       IF(Actives%I(I)==0) Vect%D(I)=PrISR(I,GDIISMemory)
     ENDDO
     IF(Print) CALL PrtIntCoords(IntCs,Vect%D,'predicted internals ')
     !
     CALL InternalToCart(XYZ,IntCs,Vect%D, &
       Print,GBackTrf,GTrfCtrl,GCoordCtrl,GConstr,SCRPath)
     CALL Delete(Vect)
   END SUBROUTINE IntCSum
!
!---------------------------------------------------------------------
!
END MODULE GDIISMod
