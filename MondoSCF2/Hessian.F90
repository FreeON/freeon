MODULE HessianMod
USE DerivedTypes
USE GlobalScalars
USE GlobalObjects
USE GlobalCharacters
USE InOut
USE InCoords
USE MemMan
USE SetXYZ
USE ProcessControl
USE PrettyPrint
USE ParsingConstants
USE GeomOptKeys
USE PunchHDF
USE LinAlg
USE AInv
USE CholFactor
USE ControlStructures
IMPLICIT NONE
CONTAINS
!
!-------------------------------------------------------
!
   SUBROUTINE DiagHess(CoordC,Hess,Grad,Displ,IntCs,AtNum,iGEO,XYZ)
     TYPE(CoordCtrl)   :: CoordC
     TYPE(Hessian)     :: Hess
     TYPE(DBL_VECT)    :: Grad,Displ 
     TYPE(INTC)        :: IntCs       
     INTEGER           :: I,J,NIntC,NCart,NConstr,NatmsLoc,iGEO
     REAL(DOUBLE)      :: HStre,HBend,HLinB,HOutP,HTors
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(:)   :: AtNum
     TYPE(DBL_VECT)    :: InvHess
     !
     NIntC=SIZE(IntCs%Def)
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     CALL New(InvHess,NIntC)
     !
    !IF(iGEO==1) THEN
       CALL DiagHess_Vals(InvHess%D,AtNum,XYZ,IntCs,CoordC, &
                          Hess,CoordC%NStre,iGEO,'ThreeVals')
      !CALL DiagHess_Vals(InvHess%D,AtNum,XYZ,IntCs,CoordC, &
      !                   Hess,CoordC%NStre,iGEO,'Lindh')
    !ELSE
    !  CALL DiagHess_Num(InvHess%D
    !ENDIF
     !
     IF(CoordC%CoordType==CoordType_Cartesian) THEN
       Displ%D=-1.D0*Grad%D !!!! equivalent with stpdesc
     ELSE IF(CoordC%CoordType==CoordType_PrimInt) THEN
       DO I=1,NIntC
         Displ%D(I)=-InvHess%D(I)*Grad%D(I)
       ENDDO
     ELSE
       CALL Halt('Only Primitiv Internals are available yet.')
     ENDIF 
     ! 
     CALL Delete(InvHess)
   END SUBROUTINE DiagHess
!
!-------------------------------------------------------
!
   SUBROUTINE DiagHessRFO(GOpt,Grad,Displ,IntCs,XYZ)
     TYPE(GeomOpt)               :: GOpt
     TYPE(DBL_VECT)              :: Grad,Displ
     TYPE(DBL_RNK2)              :: Hessian
     TYPE(INTC)                  :: IntCs
     INTEGER                     :: I,J,NIntC,NCart,Info
     REAL(DOUBLE)                :: HStre,HBend,HLinB,HOutP,HTors
     REAL(DOUBLE)                :: Sum
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     !
     NIntC=SIZE(IntCs%Def)
     NCart=3*SIZE(XYZ,2)
     CALL New(Hessian,(/NIntC+1,NIntC+1/))
     Hessian%D=Zero
     !
     IF(GOpt%CoordCtrl%CoordType==CoordType_Cartesian) THEN
       Displ%D=-1.D0*Grad%D !!!! equivalent with StpDesc
     ELSE IF(GOpt%CoordCtrl%CoordType==CoordType_PrimInt) THEN
       HStre=GOpt%Hessian%Stre
       HBend=GOpt%Hessian%Bend
       HLinB=GOpt%Hessian%LinB
       HOutP=GOpt%Hessian%OutP
       HTors=GOpt%Hessian%Tors
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
     CALL UnSetDSYEVWork()
     !
     !project out redundancy
     !         CALL RedundancyOff(Displ%D,XYZ,DoSet_O=.TRUE.)
     ELSE
       CALL Halt('Only Primitiv Internals are available yet.')
     ENDIF 
     CALL Delete(Hessian)
   END SUBROUTINE DiagHessRFO
!
!---------------------------------------------------------------
!
   SUBROUTINE SetHessian(Hess)
     TYPE(Hessian) :: Hess
     Hess%Stre = 1.20D0   
     Hess%Bend = 0.20D0 
     Hess%LinB = 0.20D0
     Hess%OutP = 0.10D0 
     Hess%Tors = 0.10D0 
     Hess%VDWStre  = 1.20D0 
     Hess%VDWBend  = 0.50D0 
     Hess%VDWLinB  = 0.50D0 
     Hess%VDWOutP  = 0.50D0 
     Hess%VDWTors  = 0.50D0 
   END SUBROUTINE SetHessian
!
!-------------------------------------------------------------------
!
   SUBROUTINE LagrGradCart(GradIn,IntCs,LagrMult,SCRPath,iGEO)
     REAL(DOUBLE),DIMENSION(:,:)    :: GradIn
     TYPE(INTC)                     :: IntCs
     TYPE(BMATR)                    :: B
     REAL(DOUBLE),DIMENSION(:)      :: LagrMult
     INTEGER                        :: I,J,NIntC,NatmsLoc,NCart,iGEO
     CHARACTER(LEN=*)               :: SCRPath
     TYPE(INT_VECT)                 :: ISpB,JSpB
     TYPE(DBL_VECT)                 :: ASpB
     !
     NIntC=SIZE(IntCs%Def) 
     NatmsLoc=SIZE(GradIn,2)
     NCart=3*NatmsLoc
     CALL ReadBMATR(ISpB,JSpB,ASpB,TRIM(SCRPath)//'B')
     !
     CALL GradAddCarts(GradIn,IntCs,LagrMult)
     CALL GradAddInt(GradIn,IntCs,B,LagrMult)
     !
   END SUBROUTINE LagrGradCart
!
!-------------------------------------------------------------------
!
   SUBROUTINE GradAddCarts(GradIn,IntCs,LagrMult)
     TYPE(INTC)                  :: IntCs
     REAL(DOUBLE),DIMENSION(:,:) :: GradIn
     REAL(DOUBLE),DIMENSION(:)   :: LagrMult
     REAL(DOUBLE)                :: Sum
     INTEGER                     :: I,J,NIntC,NConstr,NDim
     !
     ! Warning! At this point IntCs%Value should contain the actual
     ! values of the coordinates
     !
     NIntC=SIZE(IntCs%Def)
     NDim=SIZE(LagrMult)
     !
     NConstr=0
     DO I=1,NIntC
       IF(IntCs%Constraint(I)) THEN
         NConstr=NConstr+1
         Sum=-LagrMult(NConstr)
         IF(NConstr>NDim) CALL Halt('Dimension error #2 in GradAddCarts')
         J=IntCs%Atoms(I,1)
         IF(IntCs%Def(I)(1:5)=='CARTX') THEN
           GradIn(1,J)=GradIn(1,J)+Sum
         ELSE IF(IntCs%Def(I)(1:5)=='CARTY') THEN
           GradIn(2,J)=GradIn(2,J)+Sum
         ELSE IF(IntCs%Def(I)(1:5)=='CARTZ') THEN
           GradIn(3,J)=GradIn(3,J)+Sum
         ENDIF
       ENDIF
     ENDDO
   END SUBROUTINE GradAddCarts
!
!-------------------------------------------------------------------
!
   SUBROUTINE GradAddInt(GradIn,IntCs,B,LagrMult)
     REAL(DOUBLE),DIMENSION(:,:):: GradIn
     REAL(DOUBLE),DIMENSION(:)  :: LagrMult
     REAL(DOUBLE)               :: Sum
     TYPE(INTC)                 :: IntCs
     TYPE(BMATR)                :: B 
     INTEGER                    :: I,J,K,L,JJ,NintC,NConstr,NDim
     !
     NintC=SIZE(IntCs%Def)
     NDim=SIZE(LagrMult)
     NConstr=0
     IF(NintC/=SIZE(B%IB,1)) CALL Halt('Dimension error in GradAddInt')
     DO I=1,NintC
       IF(IntCs%Constraint(I)) THEN
         NConstr=NConstr+1
         Sum=-LagrMult(NConstr)
         IF(NConstr>NDim) &
           CALL Halt('Dimension error #2 in GradAddInt')
         IF(IntCs%Def(I)(1:4)/='CART') THEN
           DO J=1,4
             L=IntCs%Atoms(I,J)
             IF(L==0) EXIT
             JJ=3*(J-1)
             DO K=1,3 
               GradIn(K,L)=GradIn(K,L)+Sum*B%B(I,JJ+K) 
             ENDDO
           ENDDO
         ENDIF
       ENDIF
     ENDDO
   END SUBROUTINE GradAddInt
!
!-------------------------------------------------------------------
!
   SUBROUTINE LagrGradMult(GradMult,XYZ,IntCs, &
                           LinCrit,TorsLinCrit,NConstrIn)
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(:) :: GradMult
     REAL(DOUBLE)              :: LinCrit,Sum,TorsLinCrit
     TYPE(INTC)                :: IntCs
     INTEGER                   :: I,J,NIntC,NConstr,NConstrIn
     !
     NIntC=SIZE(IntCs%Def)
     NConstr=0
     CALL INTCValue(IntCs,XYZ,LinCrit,TorsLinCrit)
     IF(SIZE(GradMult)/=NConstrIn) &
       CALL Halt('Dimension error #1 in LagrGradMult')
     DO I=1,NIntC
       IF(IntCs%Constraint(I)) THEN
         NConstr=NConstr+1
         IF(NConstr>NConstrIn) &
           CALL Halt('Dimension error #2 in LagrGradMult')
         Sum=IntCs%Value(I)-IntCs%ConstrValue(I)
         GradMult(NConstr)=-Sum
       ENDIF
     ENDDO
   END SUBROUTINE LagrGradMult
!
!----------------------------------------------------------------
!
   SUBROUTINE DiagHessLagr(GCoordCtrl,GHessian, &
        Grad,Displ,IntCs,XYZ,SCRPath,LagrMult,GradMult,LagrDispl)
     REAL(DOUBLE),DIMENSION(:)  :: Grad,Displ
     REAL(DOUBLE),DIMENSION(:)  :: LagrMult,GradMult,LagrDispl
     REAL(DOUBLE),DIMENSION(:,:):: XYZ
     TYPE(CoordCtrl)            :: GCoordCtrl
     TYPE(Hessian)              :: GHessian
     TYPE(INTC)                 :: IntCs
     CHARACTER(LEN=*)           :: SCRPath
     INTEGER                    :: I,J,K,L,NIntC,NCart
     INTEGER                    :: NatmsLoc,NConstr
     TYPE(INT_VECT)             :: IHessL,JHessL
     TYPE(DBL_VECT)             :: AHessL
     !
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     NConstr=SIZE(LagrMult)
     CALL LagrInvHess(IntCs,SCRPath,GHessian, &
                      LagrMult,NCart,NConstr, &
                      IHessL,JHessL,AHessL)
     CALL DiagDispl(IHessL,JHessL,AHessL,SCRPath, &
                    Grad,GradMult,Displ,LagrDispl)
     ! 
     CALL Delete(IHessL)
     CALL Delete(JHessL)
     CALL Delete(AHessL)
   END SUBROUTINE DiagHessLagr
!
!-------------------------------------------------------------------
!
   SUBROUTINE LagrConv(GStat,LagrDispl,LagrMult,GradMult)
     TYPE(GOptStat)            :: GStat 
     REAL(DOUBLE),DIMENSION(:) :: LagrDispl,LagrMult,GradMult
     REAL(DOUBLE)              :: Sum
     INTEGER                   :: I,J,NLagr
     !
     NLagr=SIZE(LagrMult)
     IF(NLagr==0) THEN
       GStat%IMaxLGrad=0
       GStat%MaxLGrad=Zero
       GStat%MaxDMult=Zero
       RETURN
     ENDIF
     GStat%IMaxLGrad=1
     GStat%MaxLGrad=ABS(GradMult(1))
     GStat%MaxDMult=ABS(LagrDispl(1)-LagrMult(1))
     DO I=2,NLagr
       Sum=ABS(GradMult(I)) 
       IF(GStat%MaxLGrad<Sum) THEN
         GStat%IMaxLGrad=I
         GStat%MaxLGrad=Sum
       ENDIF
       Sum=ABS(LagrDispl(I)-LagrMult(I)) 
       IF(GStat%MaxDMult<Sum) THEN
         GStat%MaxDMult=Sum
       ENDIF
     ENDDO
   END SUBROUTINE LagrConv
!
!-------------------------------------------------------------------
!
   SUBROUTINE CalcLagrMult(GradMult,LagrMult,LagrDispl,XYZ,IntCs, &
                           GradIn,GGrdTrf,GCoordCtrl,GTrfCtrl,GConstr,&
                           GHess,Print,SCRPath)
     !
     REAL(DOUBLE),DIMENSION(:)   :: GradMult,LagrMult,LagrDispl
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ,GradIn
     REAL(DOUBLE)                :: Sum,Hess
     TYPE(INTC)                  :: IntCs
     INTEGER                     :: I,J,NLagr,NConstr,NIntC,NCart
     INTEGER                     :: NatmsLoc 
     INTEGER                     :: Print
     CHARACTER(LEN=*)            :: SCRPath
     TYPE(DBL_VECT)              :: CartGrad,IntGrad
     TYPE(Constr)                :: GConstr
     TYPE(TrfCtrl)               :: GTrfCtrl
     TYPE(CoordCtrl)             :: GCoordCtrl
     TYPE(GrdTrf)                :: GGrdTrf
     TYPE(Hessian)               :: GHess   
     TYPE(INT_VECT)              :: IBc,JBc
     TYPE(DBL_VECT)              :: ABc
     !
     NLagr=SIZE(GradMult)
     NatmsLoc=SIZE(XYZ,2)
     NCart=3*NatmsLoc
     NIntC=SIZE(IntCs%Def)
     !
     ! Calculate constraint B-matrix
     !
   ! CALL BMatrConstr(IBc,JBc,ABc,XYZ,IntCs,GCoordCtrl, &
   !                  GTrfCtrl,SCRPath,NLagr,NCart)
     !
     ! Calculate new values of Lagrange multipliers 
     ! on internal coord constraints
     !
     CALL New(CartGrad,NCart)
     CALL CartRNK2ToCartRNK1(CartGrad%D,GradIn)
     CALL New(IntGrad,NIntC)
     !
     IF(GConstr%NConstr/=GConstr%NCartConstr) THEN
       CALL CartToInternal(XYZ,IntCs,CartGrad%D,IntGrad%D, &
         GGrdTrf,GCoordCtrl,GTrfCtrl,Print,SCRPath)
     ELSE
       IntGrad%D=Zero
     ENDIF
     !
     NConstr=0
     DO I=1,NIntC
       IF(IntCs%Constraint(I)) THEN
         NConstr=NConstr+1
         IF(NConstr>NLagr) CALL Halt('Dimension error in CalcLagrMult') 
         Sum=Zero
         !Hess=ABS(LagrMult(NConstr))
         Hess=GHess%Stre
         J=IntCs%Atoms(I,1)
         IF(IntCs%Def(I)(1:4)/='CART') THEN
           Sum=IntGrad%D(I)
           IF(IntCs%Def(I)(1:4)=='STRE') Hess=GHess%Stre
           IF(IntCs%Def(I)(1:4)=='BEND') Hess=GHess%Bend
           IF(IntCs%Def(I)(1:4)=='LINB') Hess=GHess%LinB
           IF(IntCs%Def(I)(1:4)=='TORS') Hess=GHess%Tors
           IF(IntCs%Def(I)(1:4)=='OUTP') Hess=GHess%OutP
         ELSE IF(IntCs%Def(I)(1:5)=='CARTX') THEN
           Sum=GradIn(1,J)
         ELSE IF(IntCs%Def(I)(1:5)=='CARTY') THEN
           Sum=GradIn(2,J)
         ELSE IF(IntCs%Def(I)(1:5)=='CARTZ') THEN
           Sum=GradIn(3,J)
         ENDIF
         LagrDispl(NConstr)=Hess*GradMult(NConstr)+Sum
       ENDIF 
     ENDDO
     !
     CALL Delete(CartGrad)
     CALL Delete(IntGrad)
!    CALL Delete(IBc)
!    CALL Delete(JBc)
!    CALL Delete(ABc)
   END SUBROUTINE CalcLagrMult
!
!-------------------------------------------------------------------
!
   SUBROUTINE LagrEnergy(LagrMult,IntCs,ELagr)
     REAL(DOUBLE),DIMENSION(:)       :: LagrMult
     REAL(DOUBLE)                    :: ELagr    
     TYPE(INTC)                      :: IntCs
     INTEGER                         :: I,J,NIntC,NLagr,NConstr
     !
     NIntC=SIZE(IntCs%Def) 
     NLagr=SIZE(LagrMult)
     ELagr=Zero
     IF(NLagr==0) RETURN     
     NConstr=0
     DO I=1,NIntC
       IF(IntCs%Constraint(I)) THEN
         NConstr=NConstr+1
         IF(NConstr>NLagr) CALL Halt('NConstr>NLagr in LagrEnergy')
         ELagr=ELagr-LagrMult(NConstr)*(IntCs%Value(I)-IntCs%ConstrValue(I)) 
       ENDIF
     ENDDO
   END SUBROUTINE LagrEnergy
!
!-------------------------------------------------------------------
!
   SUBROUTINE DiagHess_Vals(DHess,AtNum,XYZ,IntCs,CoordC,Hess, &
                            NStre,iGEO,Char)
     TYPE(Hessian)                 :: Hess
     TYPE(CoordCtrl)               :: CoordC
     REAL(DOUBLE),DIMENSION(:)     :: DHess,AtNum
     REAL(DOUBLE),DIMENSION(:,:)   :: XYZ  
     TYPE(INTC)                    :: IntCs 
     INTEGER                       :: I,J,NIntC,NatmsLoc,NStre,iGEO
     CHARACTER(LEN=*)              :: Char
     TYPE(INT_VECT)                :: ITop,JTop
     TYPE(DBL_VECT)                :: ATop
     LOGICAL                       :: DoVDW
     !
     NIntC=SIZE(IntCs%Def)
     NatmsLoc=SIZE(XYZ,2)
     IF(NIntC/=SIZE(DHess)) &
       CALL Halt('Dimesion error in DiagHess_Vals.')
     !
     IF(Char=='Lindh') CALL DistMatr(ITop,JTop,ATop,IntCs,NatmsLoc,NStre)
     !
     DO I=1,NIntC
       IF(.NOT.IntCs%Active(I)) THEN
         DHess(I)=Zero
       ELSE 
         DoVDW=(I>CoordC%NCov)
         CALL CalcHess(DHess(I),Char,IntCs%Def(I),Hess,AtNum,iGEO, &
                       XYZ,IntCs%Atoms(I,1:4),ITop,JTop,ATop,DoVDW)
       ENDIF
     ENDDO
     !
     IF(Char=='Lindh') THEN 
       CALL Delete(ITop)
       CALL Delete(JTop)
       CALL Delete(ATop)
     ENDIF
   END SUBROUTINE DiagHess_Vals
!
!-------------------------------------------------------------------
!
   FUNCTION PeriodicRow(I)
     INTEGER   ::  PeriodicRow,I
     !
     IF(I<=0) THEN
       PeriodicRow=0
     ELSE IF(0<I.AND.I<=2) THEN
       PeriodicRow=1
     ELSE IF(2<I.AND.I<=10) THEN
       PeriodicRow=2
     ELSE IF(10<I.AND.I<=18) THEN
       PeriodicRow=3
     ELSE IF(18<I.AND.I<=36) THEN
       PeriodicRow=4
     ELSE IF(36<I.AND.I<=54) THEN
       PeriodicRow=5
     ELSE IF(54<I.AND.I<=86) THEN
       PeriodicRow=5
     ELSE 
       PeriodicRow=7
     ENDIF
   END FUNCTION PeriodicRow
!
!-------------------------------------------------------------------
!
   SUBROUTINE CalcHess(DHess,Char,Type,Hess,AtNum,iGEO,XYZ, &
                       Atoms,ITop,JTop,ATop,DoVDW)
     REAL(DOUBLE)                :: DHess
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     REAL(DOUBLE),DIMENSION(:)   :: AtNum
     TYPE(Hessian)               :: Hess
     CHARACTER(LEN=*)            :: Char,Type
     INTEGER,DIMENSION(1:4)      :: Atoms 
     TYPE(INT_VECT)              :: ITop,JTop
     TYPE(DBL_VECT)              :: ATop
     INTEGER                     :: I1Row,I2Row,I3Row,I4Row,iGEO
     REAL(DOUBLE)                :: R12,R23,R34
     REAL(DOUBLE)                :: Rho12,Rho23,Rho34
     LOGICAL                     :: DoVDW
     !
     IF(Type(1:4)=='CART') THEN
       DHess=Hess%VDWStre
     ELSE IF(DoVDW) THEN
       IF(Type(1:4)=='STRE') DHess=Hess%VDWStre
       IF(Type(1:4)=='BEND') DHess=Hess%VDWBend
       IF(Type(1:4)=='LINB') DHess=Hess%VDWLinB
       IF(Type(1:4)=='OUTP') DHess=Hess%VDWOutP
       IF(Type(1:4)=='TORS') DHess=Hess%VDWTors
     ELSE
       IF(Char=='ThreeVals') THEN
         IF(Type(1:4)=='STRE') DHess=Hess%Stre
         IF(Type(1:4)=='BEND') DHess=Hess%Bend
         IF(Type(1:4)=='LINB') DHess=Hess%LinB
         IF(Type(1:4)=='OUTP') DHess=Hess%OutP
         IF(Type(1:4)=='TORS') DHess=Hess%Tors
       ELSE IF(Char=='Lindh') THEN
         I1Row=PeriodicRow(INT(AtNum(Atoms(1))))
         I2Row=PeriodicRow(INT(AtNum(Atoms(2))))
         I3Row=PeriodicRow(INT(AtNum(Atoms(3))))
         I4Row=PeriodicRow(INT(AtNum(Atoms(4))))
         IF(Atoms(2)/=0) THEN
           R12=GetR(XYZ,Atoms(1),Atoms(2),ITop,JTop,ATop) 
           Rho12=EXP(Lindh_Alpha(I1Row,I2Row)*(Lindh_R(I1Row,I2Row)**2-R12**2))
         ENDIF
         IF(Atoms(3)/=0) THEN
           R23=GetR(XYZ,Atoms(2),Atoms(3),ITop,JTop,ATop) 
           Rho23=EXP(Lindh_Alpha(I2Row,I3Row)*(Lindh_R(I2Row,I3Row)**2-R23**2))
         ENDIF
         IF(Atoms(4)/=0) THEN
           IF(Type=='OUTP') THEN
             R34=GetR(XYZ,Atoms(2),Atoms(4),ITop,JTop,ATop) 
             Rho34=EXP(Lindh_Alpha(I2Row,I4Row)*(Lindh_R(I2Row,I4Row)**2-R34**2))
           ELSE
             R34=GetR(XYZ,Atoms(3),Atoms(4),ITop,JTop,ATop) 
             Rho34=EXP(Lindh_Alpha(I3Row,I4Row)*(Lindh_R(I3Row,I4Row)**2-R34**2))
           ENDIF
         ENDIF
         DHess=One
         IF(Type(1:4)=='STRE') THEN
           DHess=Lindh_K(1)*Rho12
         ELSE IF(Type(1:4)=='BEND') THEN
           DHess=Lindh_K(2)*Rho12*Rho23
         ELSE IF(Type(1:4)=='LINB') THEN
           DHess=Lindh_K(2)*Rho12*Rho23
         ELSE IF(Type(1:4)=='OUTP') THEN
           DHess=Lindh_K(3)*Rho12*Rho23*Rho34
         ELSE IF(Type(1:4)=='TORS') THEN
           DHess=Lindh_K(3)*Rho12*Rho23*Rho34
         ELSE IF(Type(1:4)=='CART') THEN
           DHess=Lindh_K(1)*Rho12
         ENDIF
       ENDIF
     ENDIF
     DHess=One/DHess
   END SUBROUTINE CalcHess
!
!-------------------------------------------------------------------
!
   FUNCTION GetR(XYZ,I1Row,I2Row,ITop,JTop,ATop)
     REAL(DOUBLE)                  :: GetR
     REAL(DOUBLE),DIMENSION(:,:)   :: XYZ 
     INTEGER                       :: I1Row,I2Row,I,J,III,NDim
     TYPE(INT_VECT)                :: ITop,JTop
     TYPE(DBL_VECT)                :: ATop
     !
     NDim=SIZE(ITop%I)-1
     III=0
     DO J=ITop%I(I1Row),ITop%I(I1Row+1)-1
       IF(JTop%I(J)==I2Row) THEN
         GetR=ATop%D(J)
         III=1
       ENDIF
     ENDDO
     IF(III==0) THEN
       GetR=SQRT((XYZ(1,I1Row)-XYZ(1,I2Row))**2+&
                 (XYZ(2,I1Row)-XYZ(2,I2Row))**2+&
                 (XYZ(3,I1Row)-XYZ(3,I2Row))**2)
     ENDIF 
   END FUNCTION GetR
!
!-------------------------------------------------------------------
!
   SUBROUTINE LagrInvHess(IntCs,SCRPath,GHessian, &
                          LagrMult,NCart,NConstr, &
                          IHessM,JHessM,AHessM)
     TYPE(BMATR)      :: B
     TYPE(Hessian)    :: GHessian
     TYPE(INTC)       :: IntCs
     REAL(DOUBLE),DIMENSION(:) :: LagrMult
     INTEGER          :: I,J,NIntC,NCart,NConstr
     CHARACTER(LEN=*) :: SCRPath 
     TYPE(INT_VECT)   :: IHessXX,JHessXX,IHessXL,JHessXL
     TYPE(INT_VECT)   :: IHessM,JHessM
     TYPE(DBL_VECT)   :: AHessXX,AHessXL,AHessM
     TYPE(INT_VECT)   :: ISpB,JSpB
     TYPE(DBL_VECT)   :: ASpB
     !
     NIntC=SIZE(IntCs%Def)
     !
     CALL ReadBMATR(ISpB,JSpB,ASpB,TRIM(SCRPath)//'B')
     IF(SIZE(B%IB,1)/=NIntC) &
       CALL Halt('Dimension error in LagrInvHess.')
     CALL GetHessXX(IntCs,B,GHessian,LagrMult, &
                    IHessXX,JHessXX,AHessXX,NCart,SCRPath_O=SCRPath)
     CALL GetHessXL(IntCs,B,IHessXL,JHessXL,AHessXL, &
                    NCart,NConstr,SCRPath_O=SCRPath)
     CALL MergeXLXX(IHessXL,JHessXL,AHessXL, &
                    IHessXX,JHessXX,AHessXX, &
                    IHessM,JHessM,AHessM,SCRPath_O=SCRPath)
     !
     CALL Delete(IHessXX)
     CALL Delete(JHessXX)
     CALL Delete(AHessXX)
     CALL Delete(IHessXL)
     CALL Delete(JHessXL)
     CALL Delete(AHessXL)
     CALL Delete(B)
   END SUBROUTINE LagrInvHess
!
!------------------------------------------------------------------
!
   SUBROUTINE GetHessXX(IntCs,B,GHessian,LagrMult, &
                        IHessXX,JHessXX,AHessXX,NCart,SCRPath_O)
     TYPE(INTC)       :: IntCs
     INTEGER          :: I,J,NIntC,NCart,NZSpB,NConstraint
     REAL(DOUBLE),DIMENSION(:) :: LagrMult
     TYPE(Hessian)    :: GHessian
     TYPE(BMATR)      :: B,BSc
     REAL(DOUBLE)     :: BScale
     TYPE(INT_VECT)   :: ISpB,JSpB,ISpBSc,JSpBSc,IHessXX,JHessXX
     TYPE(INT_VECT)   :: ISpBt,JSpBt
     TYPE(DBL_VECT)   :: ASpB,ASpBSc,AHessXX
     TYPE(DBL_VECT)   :: ASpBt
     CHARACTER(LEN=*) :: SCRPath_O
     !
     CALL Set_BMATR_EQ_BMATR(BSc,B)
     NIntC=SIZE(IntCs%Def)
     !
     NConstraint=0
     DO I=1,NIntC
       IF(IntCs%Def(I)(1:4)=='STRE') THEN
         BScale=GHessian%Stre
       ELSE IF(IntCs%Def(I)(1:4)=='BEND') THEN
         BScale=GHessian%Bend
       ELSE IF(IntCs%Def(I)(1:4)=='TORS') THEN
         BScale=GHessian%Tors
       ELSE IF(IntCs%Def(I)(1:4)=='OUTP') THEN
         BScale=GHessian%OutP
       ELSE IF(IntCs%Def(I)(1:4)=='LINB') THEN
         BScale=GHessian%LinB
       ELSE IF(IntCs%Def(I)(1:4)=='CART') THEN
         BScale=Zero
       ELSE
         BScale=One
       ENDIF
       IF(IntCs%Constraint(I)) THEN
         NConstraint=NConstraint+1
         !IF(IntCs%Def(I)(1:4)/='CART') THEN
         !  BScale=BScale-LagrMult(NConstraint)
         !ENDIF
       ENDIF
       DO J=1,12 ; BSc%B(I,J)=BScale*BSc%B(I,J) ; ENDDO
     ENDDO
     !
     CALL BtoSpB_1x1(B,ISpB,JSpB,ASpB)
     CALL BtoSpB_1x1(BSc,ISpBSc,JSpBSc,ASpBSc)
     !
     CALL New(ISpBt,NCart+1)
     NZSpB=ISpB%I(NIntC+1)-1
     CALL New(JSpBt,NZSpB)
     CALL New(ASpBt,NZSpB)
     CALL TransPose1x1(ISpB%I,JSpB%I,ASpB%D,NIntC,NCart, &
          ISpBt%I,JSpBt%I,ASpBt%D,'full')
     !
     CALL MatMul_1x1(ISpBt%I,JSpBt%I,ASpBt%D, &
                     ISpBSc%I,JSpBSc%I,ASpBSc%D, &
                     IHessXX,JHessXX,AHessXX,NCart,NIntC,NCart)
     !
     CALL Delete(ISpB)
     CALL Delete(JSpB)
     CALL Delete(ASpB)
     CALL Delete(ISpBt)
     CALL Delete(JSpBt)
     CALL Delete(ASpBt)
     CALL Delete(ISpBSc)
     CALL Delete(JSpBSc)
     CALL Delete(ASpBSc)
     CALL Delete(BSc)
   END SUBROUTINE GetHessXX
!
!------------------------------------------------------------------
!
   SUBROUTINE GetHessXL(IntCs,B,IHessXLs,JHessXLs,AHessXLs, &
                        NCart,NConstr,SCRPath_O)
     TYPE(INTC)       :: IntCs
     TYPE(BMATR)      :: B
     TYPE(INT_VECT)   :: IHessXL,JHessXL,IHessXLt,JHessXLt
     TYPE(INT_VECT)   :: IHessXLs,JHessXLs
     TYPE(DBL_VECT)   :: AHessXL,AHessXLt,AHessXLs
     INTEGER          :: NCart,NConstr,NDim,NNew,NZ,I,J,K,JJ,KK,NIntC
     CHARACTER(LEN=*),OPTIONAL :: SCRPath_O
     !
     NDim=NCart+NConstr 
     NIntC=SIZE(IntCs%Def)
     CALL New(IHessXL,NDim+1)
     IHessXL%I(1:NCart+1)=1
     NNew=NCart
     DO I=1,NIntC
       IF(IntCs%Constraint(I)) THEN
         NNew=NNew+1
         IF(NNew-NCart>NConstr) CALL Halt('Dimension error in GetHessXL') 
         NZ=0
         DO J=1,4 
           JJ=B%IB(I,J)
           IF(JJ==0) EXIT
           NZ=NZ+3
         ENDDO
         IHessXL%I(NNew+1)=IHessXL%I(NNew)+NZ
       ENDIF
     ENDDO
     !
     NZ=IHessXL%I(NDim+1)-1
     CALL New(JHessXL,NZ)
     CALL New(AHessXL,NZ)
     !
     NZ=0
     NNew=0
     DO I=1,NIntC
       IF(IntCs%Constraint(I)) THEN
         NNew=NNew+1
         IF(NNew>NConstr)CALL Halt('Dimension error #2 in GetHessXL') 
         DO J=1,4 
           JJ=B%IB(I,J)
           IF(JJ==0) EXIT
           JJ=3*(JJ-1)
           DO K=1,3 
             NZ=NZ+1
             KK=3*(J-1)+K
             JHessXL%I(NZ)=JJ+K
             AHessXL%D(NZ)=-B%B(I,KK)
           ENDDO
         ENDDO
       ENDIF
     ENDDO
     !
     ! Build transpose of XL and merge XL with it.
     !
     CALL New(IHessXLt,NDim+1)
     CALL New(JHessXLt,NZ)
     CALL New(AHessXLt,NZ)
     CALL TransPose1x1(IHessXL%I,JHessXL%I,AHessXL%D, &
       NDim,NDim,IHessXLt%I,JHessXLt%I,AHessXLt%D,'full')
     CALL AddMat_1x1(IHessXL%I,JHessXL%I,AHessXL%D, &
                     IHessXLt%I,JHessXLt%I,AHessXLt%D, &
                     IHessXLs,JHessXLs,AHessXLs,NDim,NDim)
     !IF(PRESENT(SCRPath_O)) THEN
     !  CALL Plot_1x1(IHessXL%I,JHessXL%I,TRIM(SCRPath_O)//'XL',NDim)
     !  CALL Plot_1x1(IHessXLt%I,JHessXLt%I,TRIM(SCRPath_O)//'XLt',NDim)
     !  CALL Plot_1x1(IHessXLs%I,JHessXLs%I,TRIM(SCRPath_O)//'XLs',NDim)
     !ENDIF
     CALL Delete(IHessXLt)
     CALL Delete(JHessXLt)
     CALL Delete(AHessXLt)
     CALL Delete(IHessXL)
     CALL Delete(JHessXL)
     CALL Delete(AHessXL)
   END SUBROUTINE GetHessXL
!
!------------------------------------------------------------------
!
   SUBROUTINE MergeXLXX(IHessXL,JHessXL,AHessXL, &
                        IHessXX,JHessXX,AHessXX, &
                        IHessM,JHessM,AHessM,SCRPath_O)
     TYPE(INT_VECT)  :: IHessXL,JHessXL,IHessXX,JHessXX,IHessP
     TYPE(INT_VECT)  :: IHessM,JHessM
     TYPE(DBL_VECT)  :: AHessXL,AHessXX,AHessM
     INTEGER         :: NDimXL,NZXL,NDimXX,NZXX
     CHARACTER(LEN=*),OPTIONAL :: SCRPath_O
     !
     NDimXL=SIZE(IHessXL%I)-1
     NZXL=IHessXL%I(NDimXL+1)-1 
     NDimXX=SIZE(IHessXX%I)-1
     NZXX=IHessXX%I(NDimXX+1)-1 
     IF(NDimXX>NDimXL) CALL Halt('Dimension error in MergeXLXX')
     CALL New(IHessP,NDimXL+1)  
     IHessP%I(1:NDimXX+1)=IHessXX%I(1:NDimXX+1) 
     IHessP%I((NDimXX+2):(NDimXL+1))=IHessP%I(NDimXX+1)
     !
     CALL AddMat_1x1(IHessP%I,JHessXX%I,AHessXX%D, &
                     IHessXL%I,JHessXL%I,AHessXL%D, &
                     IHessM,JHessM,AHessM,NDimXL,NDimXL)
     !IF(PRESENT(SCRPath_O)) THEN
     !  CALL Plot_1x1(IHessP%I,JHessXX%I,TRIM(SCRPath_O)//'XXP',NDimXL)
     !  CALL Plot_1x1(IHessXL%I,JHessXL%I,TRIM(SCRPath_O)//'XL',NDimXL)
     !  CALL Plot_1x1(IHessM%I,JHessM%I,TRIM(SCRPath_O)//'M',NDimXL)
     !ENDIF
     CALL Delete(IHessP)
   END SUBROUTINE MergeXLXX
!
!------------------------------------------------------------------
!
   SUBROUTINE DiagDispl(IHessL,JHessL,AHessL,SCRPath, &
                        Grad,GradMult,Displ,LagrDispl)
     TYPE(INT_VECT)             :: IHessL,JHessL
     TYPE(DBL_VECT)             :: AHessL,Vect0,Vect1,Vect2
     REAL(DOUBLE),DIMENSION(:)  :: Grad,GradMult,Displ,LagrDispl
     TYPE(DBL_RNK2)             :: FullMat,InvMat
     INTEGER                    :: I,J,NLagr,NIntC,NDim,NCart
     CHARACTER(LEN=*)           :: SCRPath
     TYPE(INT_VECT)             :: ISpB,JSpB
     TYPE(DBL_VECT)             :: ASpB
     !
     NLagr=SIZE(GradMult)
     IF(NLagr/=SIZE(LagrDispl)) &
       CALL Halt('Dimension error #1 in DiagDispl.')
     NIntC=SIZE(Grad)
     IF(NIntC/=SIZE(Displ)) &
       CALL Halt('Dimension error #2 in DiagDispl.')
     NDim=SIZE(IHessL%I)-1
     NCart=NDim-NLagr
     !
     CALL ReadBMATR(ISpB,JSpB,ASpB,TRIM(SCRPath)//'B')
     IF(NIntC/=SIZE(ISpB%I,1)-1) &
       CALL Halt('Dimension error #3 in DiagDispl.')
     !
     CALL New(Vect0,NCart)
     CALL CALC_BxVect(ISpB,JSpB,ASpB,Grad,Vect0%D,Trp_O=.TRUE.)
     !
     CALL Sp1x1ToFull(IHessL%I,JHessL%I,AHessL%D,NDim,NDim,FullMat)
     !
     ! Calc. Abs.Val. Inverse of the Hessian
     !
     CALL New(InvMat,(/NDim,NDim/))
     CALL SetDSYEVWork(NDim**2)
     CALL FunkOnSqMat(NDim,AbsInv,FullMat%D,InvMat%D,Unit_O=6)
     CALL UnSetDSYEVWork()
     !
     CALL New(Vect1,NDim)
     Vect1%D(1:NCart)=-Vect0%D(1:NCart)
     IF(NDim>NCart) Vect1%D(NCart+1:NDim)=-GradMult(1:NLagr)
     !
     CALL New(Vect2,NDim)
       CALL DGEMM_NNc(NDim,NDim,1,One,Zero,InvMat%D,Vect1%D,Vect2%D)
     CALL Delete(Vect1)
     !
     Vect0%D(1:NCart)=Vect2%D(1:NCart)
     CALL CALC_BxVect(ISpB,JSpB,ASpB,Displ,Vect0%D)
     IF(NDim>NCart) LagrDispl(1:NLagr)=Vect2%D(NCart+1:NDim)
     CALL Delete(Vect2)
     !
     CALL Delete(Vect0)
     CALL Delete(InvMat)
     CALL Delete(FullMat)
     CALL Delete(ISpB)
     CALL Delete(JSpB)
     CALL Delete(ASpB)
   END SUBROUTINE DiagDispl
!
!------------------------------------------------------------------
!
   SUBROUTINE BMatrConstr(IBc,JBc,ABc,XYZ,IntCs,GCoordCtrl, &
                          GTrfCtrl,SCRPath,NLagr,NCart)
     TYPE(INT_VECT)        :: IBc,JBc,JBc2
     TYPE(DBL_VECT)        :: ABc,ABc2
     INTEGER               :: I,J,L,LL,JJ,NLagr,NCart,NIntC,NConstr,NZ
     CHARACTER(LEN=*)      :: SCRPath
     TYPE(INTC)            :: IntCs
     TYPE(BMATR)           :: B     
     REAL(DOUBLE),DIMENSION(:,:) :: XYZ
     TYPE(TrfCtrl)         :: GTrfCtrl
     TYPE(CoordCtrl)       :: GCoordCtrl
     !
     NIntC=SIZE(IntCs%Def)
     !!! calling the whole B matrix is not necessary, 
     !!! effort can be saved by calling it only for constraints !!!
     CALL BMatrix(XYZ,NIntC,IntCs,B, &
                  GCoordCtrl%LinCrit,GCoordCtrl%TorsLinCrit)   
     !
     CALL New(IBc,NLagr+1)
     CALL New(JBc2,NLagr*12)
     CALL New(ABc2,NLagr*12)
     IBc%I(1)=1
     NZ=0
     NConstr=0 
     DO I=1,NIntC
       IF(IntCs%Constraint(I)) THEN
         NConstr=NConstr+1
         IF(NConstr>NLagr)CALL Halt('Dimension error #2 in BMatrConstr')
         DO J=1,4
           LL=IntCs%Atoms(I,J)
           IF(LL==0) EXIT
           JJ=(J-1)*3
           LL=(LL-1)*3
           DO L=1,3
             NZ=NZ+1
             JBc2%I(NZ)=LL+L 
             ABc2%D(NZ)=B%B(I,JJ+L)
           ENDDO
         ENDDO
         IBc%I(NConstr+1)=NZ+1
       ENDIF
     ENDDO
     !
     CALL New(JBc,NZ)
     CALL New(ABc,NZ)
     JBc%I(1:NZ)=JBc2%I(1:NZ)
     ABc%D(1:NZ)=ABc2%D(1:NZ)
     CALL Delete(JBc2)
     CALL Delete(ABc2)
     !
     CALL Delete(B) 
   END SUBROUTINE BMatrConstr
!
!------------------------------------------------------------------
!
!  SUBROUTINE DiagHess_Num(InvHess%D,HFileIn,iCLONE)
!    CHARACTER(LEN=*)            :: HFileIn
!    INTEGER                     :: iCLONE
!    !
!    HDFFileID=OpenHDF(HFileIn)
!    HDF_CurrentID= &
!      OpenHDFGroup(HDFFileID,"Clone #"//TRIM(IntToChar(iCLONE)))
!    CALL CloseHDFGroup(HDF_CurrentID)
!    CALL CloseHDF(HDFFileID)
!  END SUBROUTINE DiagHess_Num
!
!------------------------------------------------------------------
!
END MODULE HessianMod

