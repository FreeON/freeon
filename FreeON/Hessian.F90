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
     NIntC=IntCs%N
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
     IntCs%PredGrad%D=Zero
     IntCs%PredVal%D=IntCs%Value%D+Displ%D
     IntCs%InvHess%D=InvHess%D
     CALL Delete(InvHess)
   END SUBROUTINE DiagHess
!
!-------------------------------------------------------
!
   SUBROUTINE SetHessian(Hess)
     TYPE(Hessian) :: Hess
    !Hess%Stre = 0.80D0
    !Hess%Bend = 0.30D0
    !Hess%LinB = 0.30D0
    !Hess%OutP = 0.20D0
    !Hess%Tors = 0.20D0
     Hess%Stre = 0.50D0
     Hess%Bend = 0.20D0
     Hess%LinB = 0.20D0
     Hess%OutP = 0.10D0
     Hess%Tors = 0.10D0
   END SUBROUTINE SetHessian
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
     NIntC=IntCs%N
     NatmsLoc=SIZE(XYZ,2)
     IF(NIntC/=SIZE(DHess)) &
       CALL Halt('Dimesion error in DiagHess_Vals.')
     !
     IF(Char=='Lindh') CALL DistMatr(ITop,JTop,ATop,IntCs,NatmsLoc,NStre)
     !
     DO I=1,NIntC
       IF(.NOT.IntCs%Active%L(I)) THEN
         DHess(I)=Zero
       ELSE
         CALL CalcHess(DHess(I),Char,IntCs%Def%C(I)(1:5),Hess,AtNum, &
                       iGEO,XYZ,IntCs%Atoms%I(I,1:4),ITop,JTop,ATop)
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
                       Atoms,ITop,JTop,ATop)
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
     !
     IF(Type(1:4)=='CART') THEN
       DHess=Hess%Stre
     ELSE
       IF(Char=='ThreeVals') THEN
         IF(Type(1:5)=='STRE ') THEN
           DHess=Hess%Stre
         ELSE IF(Type(1:4)=='BEND') THEN
           DHess=Hess%Bend
         ELSE IF(Type(1:4)=='LINB') THEN
           DHess=Hess%LinB
         ELSE IF(Type(1:4)=='OUTP') THEN
           DHess=Hess%OutP
         ELSE IF(Type(1:4)=='TORS') THEN
           DHess=Hess%Tors
         ELSE
           DHess=0.5D0 ! for inverse hessian of lattice
         ENDIF
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
