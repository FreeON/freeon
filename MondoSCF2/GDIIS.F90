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
   USE SetXYZ  
!
IMPLICIT NONE
!
CONTAINS
!
!--------------------------------------------------------------
!
   SUBROUTINE GeoDIIS(XYZ,GConstr,GBackTrf, &
              GGrdTrf,GTrfCtrl,GCoordCtrl,GDIISCtrl, &
              HFileIn,iCLONE,iGEO,Print,SCRPath, &
              Displ_O,Grad_O,IntGrad_O,PWD_O)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     CHARACTER(LEN=*)            :: HFileIn
     TYPE(Constr)                :: GConstr
     TYPE(BackTrf)               :: GBackTrf
     TYPE(GrdTrf)                :: GGrdTrf 
     TYPE(TrfCtrl)               :: GTrfCtrl
     TYPE(CoordCtrl)             :: GCoordCtrl 
     TYPE(GDIIS)                 :: GDIISCtrl  
     INTEGER                     :: iCLONE,iGEO,ICount
     INTEGER                     :: Print
     CHARACTER(Len=*)            :: SCRPath
     CHARACTER(Len=*),OPTIONAL   :: PWD_O
     INTEGER                     :: I,II,J,JJ,K,L,NCart,NatmsLoc,NDim
     INTEGER                     :: IGeom,HDFFileID,IStart
     INTEGER                     :: SRMemory,RefMemory
     INTEGER                     :: CartGradMemory,GDIISMemory
     TYPE(DBL_RNK2)              :: SRStruct,RefGrad
     TYPE(DBL_RNK2)              :: RefStruct,SRDispl
     TYPE(DBL_RNK2)              :: Aux
     TYPE(DBL_VECT)              :: Vect
     REAL(DOUBLE),DIMENSION(:),OPTIONAL :: Displ_O,Grad_O,IntGrad_O
     !
     IF(PRESENT(Displ_O)) THEN
       GDIISMemory=MIN(10,iGEO)
     ELSE
       GDIISMemory=MIN(6,iGEO)
     ENDIF
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
     NDim=NCart
     !
     ! Get GDIIS memory of Cartesian coords and grads
     !
     CALL New(SRStruct,(/NDim,GDIISMemory/))
     CALL New(RefStruct,(/NDim,GDIISMemory/))
     CALL New(RefGrad,(/NDim,GDIISMemory/))
     CALL New(SRDispl,(/NDim,GDIISMemory/))
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
       CALL IntCFit(XYZ,Grad_O,IntGrad_O,Displ_O, &
                    RefStruct%D,RefGrad%D,SCRPath,PWD_O, &
                    GGrdTrf,GCoordCtrl,GTrfCtrl,Print,iGEO+1)
     ELSE
       CALL BasicGDIIS(XYZ,GConstr,Print,RefStruct,RefGrad, &
                       SRStruct,SRDispl)
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
   SUBROUTINE IntCFit(XYZ,Grad,IntGrad,Displ, &
                      RefStruct,RefGrad,SCRPath,PWDPath, &
                      GGrdTrf,GCoordCtrl,GTrfCtrl,Print,iGEO)
     REAL(DOUBLE),DIMENSION(:,:)  :: XYZ,RefStruct,RefGrad
     REAL(DOUBLE),DIMENSION(:)    :: Grad,Displ,IntGrad
     CHARACTER(LEN=*)             :: SCRPath,PWDPath
     INTEGER                      :: Print,iGEO
     INTEGER                      :: I,J,NMem,NCart,NatmsLoc,NIntC
     TYPE(INTC)                   :: IntCs
     TYPE(DBL_RNK2)               :: XYZAux,IntCGrads,IntCValues
     TYPE(DBL_RNK2)               :: DelocVals,DelocGrads
     TYPE(DBL_VECT)               :: VectC,VectI,VectCG,RMSErr,VectX
     TYPE(GrdTrf)                 :: GGrdTrf
     TYPE(CoordCtrl)              :: GCoordCtrl
     TYPE(TrfCtrl)                :: GTrfCtrl
     !
write(*,*) 'chk 1'
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
     CALL New(DelocVals,(/NCart-6,NMem+1/))
     CALL New(DelocGrads,(/NCart-6,NMem+1/))
     !
     CALL ReadIntCs(IntCs,TRIM(SCRPath)//'IntCs')
     NIntC=SIZE(IntCs%Def)
     !
     CALL New(IntCGrads,(/NIntC,NMem+1/))
     CALL New(IntCValues,(/NIntC,NMem+1/))
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
       RMSErr%D(I)=SQRT(DOT_PRODUCT(VectCG%D,VectCG%D)/DBLE(NCart))
       CALL CartRNK1ToCartRNK2(VectC%D,XYZAux%D)
       CALL RefreshBMatInfo(IntCs,XYZAux%D,GTrfCtrl, &
                          GCoordCtrl,Print,SCRPath,.TRUE.)
       CALL CartToInternal(XYZAux%D,IntCs,VectCG%D,VectI%D,&
         GGrdTrf,GCoordCtrl,GTrfCtrl,Print,SCRPath)
       IntCGrads%D(:,I)=VectI%D
       CALL INTCValue(IntCs,XYZAux%D,GCoordCtrl%LinCrit, &
                      GCoordCtrl%TorsLinCrit)
       IntCValues%D(:,I)=IntCs%Value
     ENDDO
       RMSErr%D(NMem+1)=SQRT(DOT_PRODUCT(Grad,Grad)/DBLE(NCart))
       CALL RefreshBMatInfo(IntCs,XYZ,GTrfCtrl, &
                          GCoordCtrl,Print,SCRPath,.TRUE.)
     ! CALL CartToInternal(XYZ,IntCs,Grad,VectI%D,&
     !   GGrdTrf,GCoordCtrl,GTrfCtrl,Print,SCRPath)
       IntCGrads%D(:,NMem+1)=IntGrad
       CALL INTCValue(IntCs,XYZ,GCoordCtrl%LinCrit,GCoordCtrl%TorsLinCrit)
       IntCValues%D(:,NMem+1)=IntCs%Value
     !
     ! Turn primitive internals into delocalized internals
     !
   ! CALL DelocEnMass(IntCValues%D,IntCGrads%D, &
   !                  DelocVals%D,DelocGrads%D,SCRPath)
     !
     ! Calculate new values of internals by fitting
     !
    !CALL DelocFit(DelocGrads%D,DelocVals%D,RMSErr%D,Displ,SCRPath)
     CALL DisplFit(IntCs,IntCGrads%D,IntCValues%D,RMSErr%D,Displ,&
                   PWDPath,iGEO)
     !
     CALL Delete(DelocVals)
     CALL Delete(DelocGrads)
     CALL Delete(VectX)
     CALL Delete(RMSErr)
     CALL Delete(IntCGrads) 
     CALL Delete(IntCValues) 
     CALL Delete(VectI) 
     CALL Delete(VectC) 
     CALL Delete(VectCG) 
     CALL Delete(XYZAux) 
     CALL Delete(IntCs) 
   END SUBROUTINE IntCFit
!
!---------------------------------------------------------------------
!
   SUBROUTINE DisplFit(IntCs,IntCGrads,IntCValues, &
                       RMSErr,Displ,Path,iGEO)
     TYPE(INTC)                 :: IntCs
     REAL(DOUBLE),DIMENSION(:)  :: Displ,RMSErr
     REAL(DOUBLE),DIMENSION(:,:):: IntCGrads,IntCValues
     REAL(DOUBLE)               :: Center,A,B
     INTEGER                    :: I,J,NIntC,NDim,iGEO
     TYPE(DBL_VECT)             :: FitVal,VectX,VectY,Sig,PredGrad
     CHARACTER(LEN=*)           :: Path
     !
     NIntC=SIZE(IntCs%Def)
     NDim=SIZE(IntCGrads,2)
     CALL New(Sig,NDim)
     CALL New(VectX,NDim)
     CALL New(VectY,NDim)
     CALL New(FitVal,NIntC)
     CALL New(PredGrad,NIntC)
     !
     DO I=1,NIntC
       IF(.NOT.IntCs%Active(I)) THEN
         FitVal%D(I)=Zero
         CYCLE
       ENDIF
       DO J=1,NDim 
         VectX%D(J)=IntCValues(I,J) 
         VectY%D(J)=IntCGrads(I,J) 
       ENDDO
       !
       ! Reorder angle values for periodic cases
       !
       IF(IntCs%Def(I)(1:4)=='LINB'.OR. & 
          IntCs%Def(I)(1:4)=='OUTP'.OR. & 
          IntCs%Def(I)(1:4)=='TORS') THEN
         CALL PeriodicAngle(VectX%D)
       ENDIF
       !
       CALL FitXY(VectX%D,VectY%D,RMSErr,FitVal%D(I),PredGrad%D(I),A,B)
       CALL PrtFits(I,IntCs,iGEO,NDim,VectX%D,VectY%D,Path, &
                    A,B,FitVal%D(I),PredGrad%D(I),IntCs%Def(I))
       !
       IF(IntCs%Def(I)(1:4)=='LINB'.OR. & 
          IntCs%Def(I)(1:4)=='OUTP'.OR. & 
          IntCs%Def(I)(1:4)=='TORS') THEN
         CALL Loose2PIs(FitVal%D(I))
         CALL LinBTo180(FitVal%D(I)) 
       ELSE IF(IntCs%Def(I)(1:4)=='BEND') THEN
         CALL ChkBendLim(IntCValues(I,NDim),FitVal%D(I),5.D0)
         CALL BendTo180(FitVal%D(I))
       ELSE IF(IntCs%Def(I)(1:4)=='STRE') THEN
         CALL ChkStreLim(IntCValues(I,NDim),FitVal%D(I),0.3D0)
       ENDIF
       Displ(I)=FitVal%D(I)-IntCValues(I,NDim)
       CALL MapDAngle(IntCs%Def(I),IntCValues(I,NDim),Displ(I))
if(IntCs%Def(I)(1:4)=='STRE') then
write(*,*) 'fit= ',FitVal%D(I),FitVal%D(I)/angstromstoau
write(*,*) iGEO,'displ= ',i,IntCs%Def(I),IntCValues(I,NDim)/angstromstoau,Displ(I)/angstromstoau
write(*,*) igeo,' final line= ',A,B
else if(IntCs%Def(I)(1:4)/='CART') then
write(*,*) 'fit= ',FitVal%D(I),FitVal%D(I)*180.D0/PI
write(*,*) iGEO,'displ= ',i,IntCs%Def(I),IntCValues(I,NDim)*180.D0/PI,Displ(I)*180.D0/PI
write(*,*) igeo,' final line= ',A,B
endif
     ENDDO
!
     CALL PrtPred(iGEO,NDim,IntCs,IntCValues,IntCGrads, &
                  FitVal%D,PredGrad%D,Path)
     !
     CALL Delete(PredGrad)
     CALL Delete(VectY)
     CALL Delete(VectX)
     CALL Delete(FitVal)
     CALL Delete(Sig)
   END SUBROUTINE DisplFit
!
!---------------------------------------------------------------------
!
   SUBROUTINE PrtPred(iGEO,NDim,IntCs,IntCValues,IntCGrads, &
                      FitVal,PredGrad,Path)
     INTEGER                   :: iGEO,I,NDim,NIntC
     TYPE(INTC)                :: IntCs
     REAL(DOUBLE),DIMENSION(:,:) :: IntCValues,IntCGrads
     REAL(DOUBLE),DIMENSION(:) :: FitVal,PredGrad
     CHARACTER(LEN=*)          :: Path
     !
     NIntC=SIZE(IntCs%Def)
     OPEN(UNIT=91,FILE=TRIM(Path)//'Pred_'// &
          TRIM(IntToChar(IGEO)),STATUS='UNKNOWN')
     DO I=1,NIntC
       WRITE(91,12) I,IntCs%Def(I)(1:5),IntCValues(I,NDim), &
                    IntCGrads(I,NDim),FitVal(I),PredGrad(I)
     ENDDO
     12 FORMAT(I3,2X,A5,2X,3F12.6,F20.6)
     CLOSE(91)
   END SUBROUTINE PrtPred
!
!---------------------------------------------------------------------
!
   SUBROUTINE PrtFits(I,IntCs,iGEO,NDim,VectX,VectY,Path, &
                      A,B,FitVal,PredGrad,Def)
     INTEGER                   :: iGEO,NDim,I,J
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY
     REAL(DOUBLE)              :: A,B,FitVal,PredGrad,Conv
     CHARACTER(LEN=*)          :: Path,Def
     TYPE(INTC)                :: IntCs
     !
     IF(Def(1:4)=='STRE') THEN
       Conv=One/AngstromsToAu
     ELSE IF(HasAngle(Def)) THEN
       Conv=180.D0/PI
     ENDIF
     WRITE(*,*) I,IntCs%Def(I),IntCs%Atoms(I,1:4)
     OPEN(UNIT=91,FILE=TRIM(Path)//'Fit_'// &
       TRIM(IntToChar(IGEO))//'_'// &
       TRIM(IntToChar(I)),STATUS='UNKNOWN')
       WRITE(*,*) i,'Data Points= '
       DO J=1,NDim
         WRITE(*,11) iGEO-NDim+J,J,Conv*VectX(J),VectY(J), &
                                   Conv*VectX(J),A+B*VectX(J)
         WRITE(91,11) iGEO-NDim+J,J,Conv*VectX(J),VectY(J), &
                                    Conv*VectX(J),A+B*VectX(J)
       ENDDO
         WRITE(*,11) iGEO+1,iGEO+1,Conv*VectX(NDim),VectY(NDim), &
                                   Conv*FitVal,PredGrad
         WRITE(91,11) iGEO+1,iGEO+1,Conv*VectX(NDim),VectY(NDim), &
                                    Conv*FitVal,PredGrad
       11 FORMAT(I3,2X,I3,2X,5F20.8)
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
write(*,*) 'st2= ',st2
write(*,*) 'ss = ',ss 
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
   SUBROUTINE Extrapolate(VectX,VectY,Sig,FitVal,PredGrad,AFit,BFit)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY,Sig
     REAL(DOUBLE)              :: FitVal,Center,A,B,AbDev,MaxDev
     REAL(DOUBLE)              :: Filter,FitValOld,MinY,MaxY,ValY
     REAL(DOUBLE)              :: SigA,SigB,Chi2,Q,MeanDev,AccGrad
     REAL(DOUBLE)              :: PredGrad,MeanDevOld
     REAL(DOUBLE)              :: AFit,BFit
     INTEGER                   :: I,J,NDim,NFit,K,MWT
     INTEGER                   :: NFit1,NFitStart,IStart
     INTEGER                   :: ICount,ICountOld
     TYPE(DBL_VECT)            :: VectXAct,VectYAct,Devs,SigAux,SigAct
     TYPE(DBL_VECT)            :: VectAuxX,VectAuxY
     !
     AFit=Zero
     BFit=Zero
     FitVal=VectX(NDim)
     PredGrad=VectY(NDim)
     NDim=SIZE(VectX)
     CALL New(Devs,NDim)
     CALL New(VectXAct,NDim)
     CALL New(VectYAct,NDim)
     CALL New(VectAuxX,NDim)
     CALL New(VectAuxY,NDim)
     CALL New(SigAux,NDim)
     CALL New(SigAct,NDim)
     !
     NFit1=MIN(5,NDim)
     NFitStart=3
     MeanDevOld=1.D99
     DO NFit=NFitStart,NFit1    
       IStart=NDim-NFit+1
       !
       MWT=0 
       CALL Fit(VectX(IStart:NDim),VectY(IStart:NDim),NFit, &
                Sig(IStart:NDim),MWT,A,B,SigA,SigB,Chi2,Q)
   !   CALL MedFit(VectX(IStart:NDim),VectY(IStart:NDim), &
   !               NFit,A,B,AbDev)
       write(*,*) 'SigA,SigB,Chi2,Q= ',SigA,SigB,Chi2,Q
       write(*,*) nfit,'line= ',a,b     
       !
       Devs%D=Zero
       DO I=1,NFit
         K=IStart+I-1
         Devs%D(K)=ABS(VectY(K)-(A+B*VectX(K)))
       ENDDO
       write(*,*) 'Devs= ',Devs%D
       MaxDev=MAXVAL(Devs%D)
       MeanDev=SUM(Devs%D)/DBLE(NFit)
       write(*,*) 'MaxDev= ',MaxDev,' MeanDev= ',MeanDev,' MeanDevOld= ',MeanDevOld
       IF(MeanDev<MeanDevOld.OR.NFit==NFitStart) THEN
         IF(ABS(B)>1.D-7) THEN
           IF(B<Zero) THEN
write(*,*) 'original a,b= ',a,b
             A=VectY(NDim)+B*VectX(NDim)
             B=-B
write(*,*) 'modified a,b= ',a,b
           ENDIF
           PredGrad=-VectY(NDim)/Two
           FitVal=(PredGrad-A)/B
           write(*,*) 'extrap 2, predgrad= ',PredGrad 
         ELSE
           FitVal=VectX(NDim)
           PredGrad=VectY(NDim)
           write(*,*) 'flat line no change of value'  
         ENDIF
         MeanDevOld=MeanDev
         AFit=A
         BFit=B
       ENDIF
       write(*,*) ' FitVal = ',FitVal
     ENDDO
     !
     CALL Delete(SigAux)
     CALL Delete(SigAct)
     CALL Delete(VectAuxX)
     CALL Delete(VectAuxY)
     CALL Delete(VectXAct)
     CALL Delete(VectYAct)
     CALL Delete(Devs)
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
write(*,*) 'in TestExtraPol ',Val,Grad,FitVal  
     IF((FitVal-Val)*(-Grad)<Zero) THEN
       FitVal=Val-Grad
write(*,*) 'fitval is changing in TestExtraPol to ',FitVal
     ENDIF
   END SUBROUTINE TestExtraPol
!
!---------------------------------------------------------------------
!
   SUBROUTINE ChkBendLim(Val,FitVal,DeltaMax)
     REAL(DOUBLE) :: Val,FitVal,DeltaMax
     IF(FitVal>PI.OR.FitVal<Zero) THEN
write(*,*) 'fitval is changing in ChkBendLim to ',FitVal
       FitVal=Val+SIGN(DeltaMax*PI/180.D0,FitVal-Val)
     ENDIF
   END SUBROUTINE ChkBendLim
!
!---------------------------------------------------------------------
!
   SUBROUTINE ChkStreLim(Val,FitVal,DeltaMax)
     REAL(DOUBLE) :: Val,FitVal,DeltaMax
     IF(FitVal<Zero.OR.FitVal>Val*Two) THEN
write(*,*) 'fitval is changing in ChkBendLim to ',FitVal
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
write(*,*) i,'IntCValues(:,I)= ',IntCValues(:,I)
       CALL PrimToDeloc(IntCValues(:,I),DelocVals(:,I), &
                        ISpB,JSpB,ASpB,UMatr)
write(*,*) i,'DelocVals(:,I)=  ',DelocVals(:,I)
write(*,*) i,'IntCGrads(:,I)= ',IntCGrads(:,I)
       CALL PrimToDeloc(IntCGrads(:,I),DelocGrads(:,I), &
                        ISpB,JSpB,ASpB,UMatr)
write(*,*) i,'DelocGrads(:,I)= ',DelocGrads(:,I)
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
       CALL FitXY(VectX%D,VectY%D,RMSErr,FitVal,PredGrad,A,B)
       NewDelocVals%D(I)=FitVal
     ENDDO 
     !
     CALL ReadBMATR(ISpB,JSpB,ASpB,TRIM(SCRPath)//'UMatr',UMatr_O=UMatr)
     !
     ! Deloc and PrimInt displacements 
     !
     NewDelocVals%D=NewDelocVals%D-DelocVals(:,NMem)
     CALL PrimToDeloc(Displ,NewDelocVals%D,ISpB,JSpB,ASpB,UMatr, &
                      Char_O='Back')
!if(nmem==3) stop
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
   SUBROUTINE FitXY(VectX,VectY,RMSErr,FitVal,PredGrad,A,B)
     REAL(DOUBLE),DIMENSION(:) :: VectX,VectY,RMSErr
     REAL(DOUBLE)              :: FitVal,G1,G2,X1,X2,A,B,AbsG
     REAL(DOUBLE)              :: PredGrad,MaxX,MinX
     LOGICAL                   :: DoBisect
     INTEGER                   :: I,J,NMem
     !
     A=Zero
     B=Zero
     NMem=SIZE(VectX)
     MaxX=MAXVAL(VectX)
     MinX=MINVAL(VectX)
     IF(ABS(MaxX-MinX)<1.D-6) THEN
       FitVal=(MaxX+MinX)/Two
       PredGrad=Zero
       RETURN
     ENDIF
     DoBisect=(VectY(NMem)*VectY(NMem-1)<Zero)
     IF(DoBisect) THEN
write(*,*) 'doing bisection'
       G1=VectY(NMem) ; G2=VectY(NMem-1)
       X1=VectX(NMem) ; X2=VectX(NMem-1)
       B=(G1-G2)/(X1-X2) ; A=G1-B*X1
       AbsG=MIN(ABS(G1),ABS(G2))
       IF(ABS(B)>1.D-4) THEN
         FitVal=-A/B
         PredGrad=Zero
write(*,*) 'bis 1, predgrad= ',PredGrad 
        !FitVal=(SIGN(AbsG/Two,-G1)-A)/B
       ELSE
         FitVal=(X1+X2)/Two
         PredGrad=G1+G2/Two
write(*,*) 'bis 2, predgrad= ',PredGrad 
       ENDIF
write(*,*) 'fitval= ',FitVal
     ELSE
write(*,*) 'doing extrapolation'
       CALL Extrapolate(VectX,VectY,RMSErr,FitVal,PredGrad,A,B)
     ENDIF
write(*,*) 'predgrad= ',PredGrad
   END SUBROUTINE FitXY
!
!---------------------------------------------------------------------
!
END MODULE GDIISMod
