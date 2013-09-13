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
MODULE MondoPoles
  USE Derivedtypes
  USE Globals
  USE GlobalScalars
  USE GlobalObjects
  USE PoleGlobals
  USE PoleNodeType
  USE ProcessControl
  USE Indexing
  USE Parse
  USE InOut
  USE Macros
  USE Thresholding
  USE BoundingBox
  USE McMurchie
  USE QCTCIndexing
  IMPLICIT NONE
  LOGICAL MMPRINT
  !=====================================================================================================
  !
  !=====================================================================================================
  INTERFACE HGToSP
     MODULE PROCEDURE HGToSP_PoleNode,HGToSP_Bra
  END INTERFACE
  TYPE(DBL_RNK2) :: CMTmp,SMTmp
CONTAINS
  !====================================================================================
  !     Q->P
  !====================================================================================
  SUBROUTINE XLate(Q,P)
    TYPE(Pole)                :: Q,P
    REAL(DOUBLE),DIMENSION(3) :: QP
    REAL(DOUBLE)              :: PQ2
    INTEGER                   :: LP,LQ,LPQ
    !------------------------------------------------------------------------------------
    LP=P%Ell
    LQ=Q%Ell
    LPQ=MAX(LP,LQ)
    QP=Q%Center-P%Center
    CALL Regular(LPQ,QP(1),QP(2),QP(3))
    CALL XLate77(LP,LQ,P%C(0),P%S(0),Cpq(0),Spq(0),Q%C(0),Q%S(0))
  END SUBROUTINE XLate
  !====================================================================================
  !     Contract P with Q
  !====================================================================================
  SUBROUTINE CTrax(P,Q,SPKetC,SPKetS)
    TYPE(PrimPair)            :: P
    TYPE(Pole)                :: Q
    REAL(DOUBLE),DIMENSION(3) :: PQ
    INTEGER                   :: LP,LQ,LPQ
    REAL(DOUBLE),DIMENSION(0:):: SPKetC,SPKetS
    !------------------------------------------------------------------------------------
    LP=P%Ell
    LQ=Q%Ell
    LPQ=LP+LQ
    PQ=P%P-Q%Center
    CALL IrRegular(LPQ,PQ(1),PQ(2),PQ(3))
    CALL CTraX77(LP,LQ,SPKetC,SPKetS,Cpq,Spq,Q%C,Q%S)
  END SUBROUTINE CTrax
  !====================================================================================
  !
  !====================================================================================
  SUBROUTINE HGToSP_PoleNode(HG,SP)
    TYPE(Herm) :: HG
    TYPE(Pole) :: SP
    REAL(DOUBLE),DIMENSION(10)      :: W
    REAL(DOUBLE),DIMENSION(3)       :: QP
    INTEGER        :: Ell,LenHG,LenSP,I,J,EllP,EllQ,K
    REAL(DOUBLE)   :: PiZ
    !
    TYPE(Pole)                :: Q,P
    !----------------------------------------------------------------------------
    SP%C=0D0
    SP%S=0D0
    DO Ell=0,HG%Ell
       LenHG=LHGTF(Ell)
       SELECT CASE(Ell)
          INCLUDE "HGToSP_PoleNode.Inc"
       CASE DEFAULT
          CALL Halt('Bad logic in HGToSP,time to remake HGToSP.Inc')
       ENDSELECT
       IF(HG%NQ(Ell).NE.0)THEN
          EllQ=Ell
          EllP=MaxPoleEll
          CALL XLateV(MaxCluster,HG%NQ(Ell),EllP,EllQ,XLLen(EllP,EllQ),      &
               XLSgn(EllP,EllQ)%D(1,1),XLIdx(EllP,EllQ)%I(1,1),               &
               FactOlm0(0),FactOlm2(0),SP%Center,HG%Cent(EllQ)%D(1,1),        &
               CMTmp%D(1,0),SMTmp%D(1,0),SP%C(0),SP%S(0))
       ENDIF
    ENDDO
  END SUBROUTINE HGToSP_PoleNode

  SUBROUTINE HGToSP_HGRho(PFFFEll,PFFFLen,Rho,Center,RhoC,RhoS)
    INTEGER                           :: PFFFEll,PFFFLen
    TYPE(HGRho)                       :: Rho
    REAL(DOUBLE),DIMENSION(0:PFFFLen) :: RhoC,RhoS
    TYPE(DBL_RNK2)                    :: PTmp,CTmp,STmp
    INTEGER                           :: Zq,Nq,OffQ,OffR,EllP,EllQ,LenP,LenQ
    INTEGER                           :: I,J,K,IAdd,JAdd,L,U,M,N
    REAL(DOUBLE)                      :: Zeta,PiZ,RTmp
    REAL(DOUBLE),DIMENSION(10)        :: W
    REAL(DOUBLE),DIMENSION(3)         :: Center



    TYPE(Pole)                :: SP,Q,P
    !----------------------------------------------------------------------------
    RhoC=0D0
    RhoS=0D0
    EllP=MaxPFFFEll
    LenP=LSP(EllP)
    CALL New(PTmp,(/3,128/))
    DO zq=1,Rho%NExpt
       NQ=Rho%NQ%I(zq)
       Zeta=Rho%Expt%D(zq)
       PiZ=(Pi/Zeta)**(3D0/2D0)
       OffQ=Rho%OffQ%I(zq)
       OffR=Rho%OffR%I(zq)
       EllQ=Rho%Lndx%I(zq)
       LenQ=LSP(EllQ)
       !
       CALL New(CTmp,(/128,LenQ/),(/1,0/))
       CALL New(STmp,(/128,LenQ/),(/1,0/))
       !
       M=MOD(Nq,128)
       L=1
       N=0
       U=M
       !
       SELECT CASE(EllQ)
          INCLUDE "HGToSP_Density.Inc"
       CASE DEFAULT
          CALL Halt('Bad logic in HGToSP,time to remake HGToSP.Inc')
       ENDSELECT
       !
       CALL XLateV(128,M,EllP,EllQ,XLLen(EllP,EllQ),              &
            XLSgn(EllP,EllQ)%D(1,1),XLIdx(EllP,EllQ)%I(1,1),      &
            FactOlm0(0),FactOlm2(0),Center,PTmp%D(1,1),           &
            CTmp%D(1,0),STmp%D(1,0),RhoC(0),RhoS(0))
       !
       L=M+1
       DO I=M+1,Nq,128
          !
          N=0
          U=L+128-1
          !
          SELECT CASE(EllQ)
             INCLUDE "HGToSP_Density.Inc"
          CASE DEFAULT
             CALL Halt('Bad logic in HGToSP,time to remake HGToSP.Inc')
          ENDSELECT
          !
          CALL XLateV(128,128,EllP,EllQ,XLLen(EllP,EllQ),              &
               XLSgn(EllP,EllQ)%D(1,1),XLIdx(EllP,EllQ)%I(1,1),      &
               FactOlm0(0),FactOlm2(0),Center,PTmp%D(1,1),           &
               CTmp%D(1,0),STmp%D(1,0),RhoC(0),RhoS(0))
          !
          L=U+1
          !
       ENDDO
       !
       CALL Delete(CTmp)
       CALL Delete(STmp)
       !
    ENDDO
    !
    CALL Delete(PTmp)
    !
  END SUBROUTINE HGToSP_HGRho


  SUBROUTINE PrntSH(Ell,C,S)
    INTEGER Ell,L,M,lmdx
    REAL(DOUBLE) C(0:),S(0:)
    WRITE(*,*)'==============================================='
    DO l=0,Ell
       DO m=0,l
          lmdx=l*(l+1)/2+m
          WRITE(*,99)l,m,lmdx,C(lmdx),S(lmdx)
99        FORMAT(3(i3),2(2x,E12.6))
       ENDDO
    ENDDO
  END SUBROUTINE PrntSH


  !====================================================================================
  !
  !====================================================================================
  SUBROUTINE HGToSP_Bra(Zeta,Ell,HGBra,SPBraC,SPBraS)
    INTEGER                           :: Ell,LenHG,LenSP
    REAL(DOUBLE)                      :: Zeta,PiZ
    REAL(DOUBLE), DIMENSION(1:)       :: HGBra
    REAL(DOUBLE), DIMENSION(0:)       :: SPBraC,SPBraS
    !------------------------------------------------------------------------------------
    !        Transform <Bra| coefficients from HG to SP
    PiZ=(Pi/Zeta)**(ThreeHalves)
    LenHG=LHGTF(Ell)
    LenSP=LSP(Ell)
    CALL HGToSP_Gen(Ell,PiZ,HGBra(1:LenHG),SPBraC(0:LenSP),SPBraS(0:LenSP))
  END SUBROUTINE HGToSP_Bra
  !====================================================================================
  !
  !====================================================================================
  SUBROUTINE HGToSP_Direct(L,LenHGTF,LenSP,PiZ,HGCo,C,S)
    INTEGER                    :: L,LenHGTF,LenSP
    REAL(DOUBLE)               :: PiZ
    REAL(DOUBLE),DIMENSION(1:LenHGTF) :: HGCo
    REAL(DOUBLE),DIMENSION(0:LenSP) :: C,S
    REAL(DOUBLE),DIMENSION(20) :: W
    !
    SELECT CASE(L)
    CASE(0)
       C(0)=PiZ*HGCo(1)
       S(0)=0
    CASE(1)
       C(0)=PiZ*HGCo(1)
       S(0)=0
       C(1)=PiZ*HGCo(4)
       S(1)=0
       C(2)=-5.d-1*PiZ*HGCo(3)
       S(2)=-5.d-1*PiZ*HGCo(2)
    CASE(2)
       C(0)=PiZ*HGCo(1)
       S(0)=0
       C(1)=PiZ*HGCo(4)
       S(1)=0
       C(2)=-5.d-1*PiZ*HGCo(3)
       S(2)=-5.d-1*PiZ*HGCo(2)
       C(3)=-5.d-1*PiZ*(HGCo(5)+HGCo(7)-2.d0*HGCo(10))
       S(3)=0
       C(4)=-5.d-1*PiZ*HGCo(9)
       S(4)=-5.d-1*PiZ*HGCo(8)
       C(5)=-2.5d-1*PiZ*(HGCo(5)-HGCo(7))
       S(5)=2.5d-1*PiZ*HGCo(6)
    CASE(3)
       C(0)=PiZ*HGCo(1)
       S(0)=0
       C(1)=PiZ*HGCo(4)
       S(1)=0
       C(2)=-5.d-1*PiZ*HGCo(3)
       S(2)=-5.d-1*PiZ*HGCo(2)
       C(3)=-5.d-1*PiZ*(HGCo(5)+HGCo(7)-2.d0*HGCo(10))
       S(3)=0
       C(4)=-5.d-1*PiZ*HGCo(9)
       S(4)=-5.d-1*PiZ*HGCo(8)
       C(5)=-2.5d-1*PiZ*(HGCo(5)-HGCo(7))
       S(5)=2.5d-1*PiZ*HGCo(6)
       C(6)=-5.d-1*PiZ*(HGCo(15)+HGCo(17)-2.d0*HGCo(20))
       S(6)=0
       C(7)=1.25d-1*PiZ*(HGCo(12)+3.d0*HGCo(14)-4.d0*HGCo(19))
       S(7)=1.25d-1*PiZ*(3.d0*HGCo(11)+HGCo(13)-4.d0*HGCo(18))
       C(8)=-2.5d-1*PiZ*(HGCo(15)-HGCo(17))
       S(8)=2.5d-1*PiZ*HGCo(16)
       C(9)=1.25d-1*PiZ*(HGCo(12)-HGCo(14))
       S(9)=1.25d-1*PiZ*(HGCo(11)-HGCo(13))
    CASE(4)
       C(0)=PiZ*HGCo(1)
       S(0)=0
       C(1)=PiZ*HGCo(4)
       S(1)=0
       C(2)=-5.d-1*PiZ*HGCo(3)
       S(2)=-5.d-1*PiZ*HGCo(2)
       C(3)=-5.d-1*PiZ*(HGCo(5)+HGCo(7)-2.d0*HGCo(10))
       S(3)=0
       C(4)=-5.d-1*PiZ*HGCo(9)
       S(4)=-5.d-1*PiZ*HGCo(8)
       C(5)=-2.5d-1*PiZ*(HGCo(5)-HGCo(7))
       S(5)=2.5d-1*PiZ*HGCo(6)
       C(6)=-5.d-1*PiZ*(HGCo(15)+HGCo(17)-2.d0*HGCo(20))
       S(6)=0
       C(7)=1.25d-1*PiZ*(HGCo(12)+3.d0*HGCo(14)-4.d0*HGCo(19))
       S(7)=1.25d-1*PiZ*(3.d0*HGCo(11)+HGCo(13)-4.d0*HGCo(18))
       C(8)=-2.5d-1*PiZ*(HGCo(15)-HGCo(17))
       S(8)=2.5d-1*PiZ*HGCo(16)
       C(9)=1.25d-1*PiZ*(HGCo(12)-HGCo(14))
       S(9)=1.25d-1*PiZ*(HGCo(11)-HGCo(13))
       W(1)=3.d0*HGCo(21)+HGCo(23)+3.d0*HGCo(25)
       W(2)=-4.d0*HGCo(30)-4.d0*HGCo(32)+8.d0*HGCo(35)
       C(10)=1.25d-1*PiZ*(W(1)+W(2))
       S(10)=0
       C(11)=1.25d-1*PiZ*(HGCo(27)+3.d0*HGCo(29)-4.d0*HGCo(34))
       S(11)=1.25d-1*PiZ*(3.d0*HGCo(26)+HGCo(28)-4.d0*HGCo(33))
       C(12)=2.5d-1*PiZ*(HGCo(21)-HGCo(25)-HGCo(30)+HGCo(32))
       S(12)=-1.25d-1*PiZ*(HGCo(22)+HGCo(24)-2.d0*HGCo(31))
       C(13)=1.25d-1*PiZ*(HGCo(27)-HGCo(29))
       S(13)=1.25d-1*PiZ*(HGCo(26)-HGCo(28))
       C(14)=6.25d-2*PiZ*(HGCo(21)-HGCo(23)+HGCo(25))
       S(14)=-6.25d-2*PiZ*(HGCo(22)-HGCo(24))
    CASE(5)
       C(0)=PiZ*HGCo(1)
       S(0)=0
       C(1)=PiZ*HGCo(4)
       S(1)=0
       C(2)=-5.d-1*PiZ*HGCo(3)
       S(2)=-5.d-1*PiZ*HGCo(2)
       C(3)=-5.d-1*PiZ*(HGCo(5)+HGCo(7)-2.d0*HGCo(10))
       S(3)=0
       C(4)=-5.d-1*PiZ*HGCo(9)
       S(4)=-5.d-1*PiZ*HGCo(8)
       C(5)=-2.5d-1*PiZ*(HGCo(5)-HGCo(7))
       S(5)=2.5d-1*PiZ*HGCo(6)
       C(6)=-5.d-1*PiZ*(HGCo(15)+HGCo(17)-2.d0*HGCo(20))
       S(6)=0
       C(7)=1.25d-1*PiZ*(HGCo(12)+3.d0*HGCo(14)-4.d0*HGCo(19))
       S(7)=1.25d-1*PiZ*(3.d0*HGCo(11)+HGCo(13)-4.d0*HGCo(18))
       C(8)=-2.5d-1*PiZ*(HGCo(15)-HGCo(17))
       S(8)=2.5d-1*PiZ*HGCo(16)
       C(9)=1.25d-1*PiZ*(HGCo(12)-HGCo(14))
       S(9)=1.25d-1*PiZ*(HGCo(11)-HGCo(13))
       W(1)=3.d0*HGCo(21)+HGCo(23)+3.d0*HGCo(25)
       W(2)=-4.d0*HGCo(30)-4.d0*HGCo(32)+8.d0*HGCo(35)
       C(10)=1.25d-1*PiZ*(W(1)+W(2))
       S(10)=0
       C(11)=1.25d-1*PiZ*(HGCo(27)+3.d0*HGCo(29)-4.d0*HGCo(34))
       S(11)=1.25d-1*PiZ*(3.d0*HGCo(26)+HGCo(28)-4.d0*HGCo(33))
       C(12)=2.5d-1*PiZ*(HGCo(21)-HGCo(25)-HGCo(30)+HGCo(32))
       S(12)=-1.25d-1*PiZ*(HGCo(22)+HGCo(24)-2.d0*HGCo(31))
       C(13)=1.25d-1*PiZ*(HGCo(27)-HGCo(29))
       S(13)=1.25d-1*PiZ*(HGCo(26)-HGCo(28))
       C(14)=6.25d-2*PiZ*(HGCo(21)-HGCo(23)+HGCo(25))
       S(14)=-6.25d-2*PiZ*(HGCo(22)-HGCo(24))
       W(1)=3.d0*HGCo(42)+HGCo(44)+3.d0*HGCo(46)
       W(2)=-4.d0*HGCo(51)-4.d0*HGCo(53)+8.d0*HGCo(56)
       C(15)=1.25d-1*PiZ*(W(1)+W(2))
       S(15)=0
       W(1)=HGCo(37)+HGCo(39)+5.d0*HGCo(41)
       W(2)=-2.d0*HGCo(48)-6.d0*HGCo(50)+8.d0*HGCo(55)
       C(16)=-6.25d-2*PiZ*(W(1)+W(2))
       W(1)=5.d0*HGCo(36)+HGCo(38)+HGCo(40)
       W(2)=-6.d0*HGCo(47)-2.d0*HGCo(49)+8.d0*HGCo(54)
       S(16)=-6.25d-2*PiZ*(W(1)+W(2))
       C(17)=2.5d-1*PiZ*(HGCo(42)-HGCo(46)-HGCo(51)+HGCo(53))
       S(17)=-1.25d-1*PiZ*(HGCo(43)+HGCo(45)-2.d0*HGCo(52))
       W(1)=3.d0*HGCo(37)+HGCo(39)-5.d0*HGCo(41)-4.d0*HGCo(48)+4.d0*HGCo(50)
       C(18)=-3.125d-2*PiZ*W(1)
       W(1)=5.d0*HGCo(36)-HGCo(38)
       W(2)=-3.d0*HGCo(40)-4.d0*HGCo(47)+4.d0*HGCo(49)
       S(18)=-3.125d-2*PiZ*(W(1)+W(2))
       C(19)=6.25d-2*PiZ*(HGCo(42)-HGCo(44)+HGCo(46))
       S(19)=-6.25d-2*PiZ*(HGCo(43)-HGCo(45))
       C(20)=-3.125d-2*PiZ*(HGCo(37)-HGCo(39)+HGCo(41))
       S(20)=-3.125d-2*PiZ*(HGCo(36)-HGCo(38)+HGCo(40))
    CASE(6)
       C(0)=PiZ*HGCo(1)
       S(0)=0
       C(1)=PiZ*HGCo(4)
       S(1)=0
       C(2)=-5.d-1*PiZ*HGCo(3)
       S(2)=-5.d-1*PiZ*HGCo(2)
       C(3)=-5.d-1*PiZ*(HGCo(5)+HGCo(7)-2.d0*HGCo(10))
       S(3)=0
       C(4)=-5.d-1*PiZ*HGCo(9)
       S(4)=-5.d-1*PiZ*HGCo(8)
       C(5)=-2.5d-1*PiZ*(HGCo(5)-HGCo(7))
       S(5)=2.5d-1*PiZ*HGCo(6)
       C(6)=-5.d-1*PiZ*(HGCo(15)+HGCo(17)-2.d0*HGCo(20))
       S(6)=0
       C(7)=1.25d-1*PiZ*(HGCo(12)+3.d0*HGCo(14)-4.d0*HGCo(19))
       S(7)=1.25d-1*PiZ*(3.d0*HGCo(11)+HGCo(13)-4.d0*HGCo(18))
       C(8)=-2.5d-1*PiZ*(HGCo(15)-HGCo(17))
       S(8)=2.5d-1*PiZ*HGCo(16)
       C(9)=1.25d-1*PiZ*(HGCo(12)-HGCo(14))
       S(9)=1.25d-1*PiZ*(HGCo(11)-HGCo(13))
       W(1)=3.d0*HGCo(21)+HGCo(23)+3.d0*HGCo(25)
       W(2)=-4.d0*HGCo(30)-4.d0*HGCo(32)+8.d0*HGCo(35)
       C(10)=1.25d-1*PiZ*(W(1)+W(2))
       S(10)=0
       C(11)=1.25d-1*PiZ*(HGCo(27)+3.d0*HGCo(29)-4.d0*HGCo(34))
       S(11)=1.25d-1*PiZ*(3.d0*HGCo(26)+HGCo(28)-4.d0*HGCo(33))
       C(12)=2.5d-1*PiZ*(HGCo(21)-HGCo(25)-HGCo(30)+HGCo(32))
       S(12)=-1.25d-1*PiZ*(HGCo(22)+HGCo(24)-2.d0*HGCo(31))
       C(13)=1.25d-1*PiZ*(HGCo(27)-HGCo(29))
       S(13)=1.25d-1*PiZ*(HGCo(26)-HGCo(28))
       C(14)=6.25d-2*PiZ*(HGCo(21)-HGCo(23)+HGCo(25))
       S(14)=-6.25d-2*PiZ*(HGCo(22)-HGCo(24))
       W(1)=3.d0*HGCo(42)+HGCo(44)+3.d0*HGCo(46)
       W(2)=-4.d0*HGCo(51)-4.d0*HGCo(53)+8.d0*HGCo(56)
       C(15)=1.25d-1*PiZ*(W(1)+W(2))
       S(15)=0
       W(1)=HGCo(37)+HGCo(39)+5.d0*HGCo(41)
       W(2)=-2.d0*HGCo(48)-6.d0*HGCo(50)+8.d0*HGCo(55)
       C(16)=-6.25d-2*PiZ*(W(1)+W(2))
       W(1)=5.d0*HGCo(36)+HGCo(38)+HGCo(40)
       W(2)=-6.d0*HGCo(47)-2.d0*HGCo(49)+8.d0*HGCo(54)
       S(16)=-6.25d-2*PiZ*(W(1)+W(2))
       C(17)=2.5d-1*PiZ*(HGCo(42)-HGCo(46)-HGCo(51)+HGCo(53))
       S(17)=-1.25d-1*PiZ*(HGCo(43)+HGCo(45)-2.d0*HGCo(52))
       W(1)=3.d0*HGCo(37)+HGCo(39)-5.d0*HGCo(41)-4.d0*HGCo(48)+4.d0*HGCo(50)
       C(18)=-3.125d-2*PiZ*W(1)
       W(1)=5.d0*HGCo(36)-HGCo(38)
       W(2)=-3.d0*HGCo(40)-4.d0*HGCo(47)+4.d0*HGCo(49)
       S(18)=-3.125d-2*PiZ*(W(1)+W(2))
       C(19)=6.25d-2*PiZ*(HGCo(42)-HGCo(44)+HGCo(46))
       S(19)=-6.25d-2*PiZ*(HGCo(43)-HGCo(45))
       C(20)=-3.125d-2*PiZ*(HGCo(37)-HGCo(39)+HGCo(41))
       S(20)=-3.125d-2*PiZ*(HGCo(36)-HGCo(38)+HGCo(40))
       W(1)=5.d0*HGCo(57)+HGCo(59)+HGCo(61)+5.d0*HGCo(63)-6.d0*HGCo(70)
       W(2)=-2.d0*HGCo(72)-6.d0*HGCo(74)
       W(3)=8.d0*HGCo(79)+8.d0*HGCo(81)-1.6d1*HGCo(84)
       C(21)=-6.25d-2*PiZ*(W(1)+W(2)+W(3))
       S(21)=0
       W(1)=HGCo(65)+HGCo(67)+5.d0*HGCo(69)
       W(2)=-2.d0*HGCo(76)-6.d0*HGCo(78)+8.d0*HGCo(83)
       C(22)=-6.25d-2*PiZ*(W(1)+W(2))
       W(1)=5.d0*HGCo(64)+HGCo(66)+HGCo(68)
       W(2)=-6.d0*HGCo(75)-2.d0*HGCo(77)+8.d0*HGCo(82)
       S(22)=-6.25d-2*PiZ*(W(1)+W(2))
       W(1)=1.5d1*HGCo(57)+HGCo(59)-HGCo(61)-1.5d1*HGCo(63)
       W(2)=-1.6d1*HGCo(70)+1.6d1*HGCo(74)+1.6d1*HGCo(79)-1.6d1*HGCo(81)
       C(23)=-1.5625d-2*PiZ*(W(1)+W(2))
       W(1)=5.d0*HGCo(58)+3.d0*HGCo(60)+5.d0*HGCo(62)
       W(2)=-8.d0*HGCo(71)-8.d0*HGCo(73)+1.6d1*HGCo(80)
       S(23)=1.5625d-2*PiZ*(W(1)+W(2))
       W(1)=3.d0*HGCo(65)+HGCo(67)-5.d0*HGCo(69)-4.d0*HGCo(76)+4.d0*HGCo(78)
       C(24)=-3.125d-2*PiZ*W(1)
       W(1)=5.d0*HGCo(64)-HGCo(66)
       W(2)=-3.d0*HGCo(68)-4.d0*HGCo(75)+4.d0*HGCo(77)
       S(24)=-3.125d-2*PiZ*(W(1)+W(2))
       W(1)=3.d0*HGCo(57)-HGCo(59)-HGCo(61)
       W(2)=3.d0*HGCo(63)-2.d0*HGCo(70)+2.d0*HGCo(72)-2.d0*HGCo(74)
       C(25)=-3.125d-2*PiZ*(W(1)+W(2))
       S(25)=6.25d-2*PiZ*(HGCo(58)-HGCo(62)-HGCo(71)+HGCo(73))
       C(26)=-3.125d-2*PiZ*(HGCo(65)-HGCo(67)+HGCo(69))
       S(26)=-3.125d-2*PiZ*(HGCo(64)-HGCo(66)+HGCo(68))
       C(27)=-1.5625d-2*PiZ*(HGCo(57)-HGCo(59)+HGCo(61)-HGCo(63))
       S(27)=1.5625d-2*PiZ*(HGCo(58)-HGCo(60)+HGCo(62))
    CASE(7)
       C(0)=PiZ*HGCo(1)
       S(0)=0
       C(1)=PiZ*HGCo(4)
       S(1)=0
       C(2)=-5.d-1*PiZ*HGCo(3)
       S(2)=-5.d-1*PiZ*HGCo(2)
       C(3)=-5.d-1*PiZ*(HGCo(5)+HGCo(7)-2.d0*HGCo(10))
       S(3)=0
       C(4)=-5.d-1*PiZ*HGCo(9)
       S(4)=-5.d-1*PiZ*HGCo(8)
       C(5)=-2.5d-1*PiZ*(HGCo(5)-HGCo(7))
       S(5)=2.5d-1*PiZ*HGCo(6)
       C(6)=-5.d-1*PiZ*(HGCo(15)+HGCo(17)-2.d0*HGCo(20))
       S(6)=0
       C(7)=1.25d-1*PiZ*(HGCo(12)+3.d0*HGCo(14)-4.d0*HGCo(19))
       S(7)=1.25d-1*PiZ*(3.d0*HGCo(11)+HGCo(13)-4.d0*HGCo(18))
       C(8)=-2.5d-1*PiZ*(HGCo(15)-HGCo(17))
       S(8)=2.5d-1*PiZ*HGCo(16)
       C(9)=1.25d-1*PiZ*(HGCo(12)-HGCo(14))
       S(9)=1.25d-1*PiZ*(HGCo(11)-HGCo(13))
       W(1)=3.d0*HGCo(21)+HGCo(23)+3.d0*HGCo(25)
       W(2)=-4.d0*HGCo(30)-4.d0*HGCo(32)+8.d0*HGCo(35)
       C(10)=1.25d-1*PiZ*(W(1)+W(2))
       S(10)=0
       C(11)=1.25d-1*PiZ*(HGCo(27)+3.d0*HGCo(29)-4.d0*HGCo(34))
       S(11)=1.25d-1*PiZ*(3.d0*HGCo(26)+HGCo(28)-4.d0*HGCo(33))
       C(12)=2.5d-1*PiZ*(HGCo(21)-HGCo(25)-HGCo(30)+HGCo(32))
       S(12)=-1.25d-1*PiZ*(HGCo(22)+HGCo(24)-2.d0*HGCo(31))
       C(13)=1.25d-1*PiZ*(HGCo(27)-HGCo(29))
       S(13)=1.25d-1*PiZ*(HGCo(26)-HGCo(28))
       C(14)=6.25d-2*PiZ*(HGCo(21)-HGCo(23)+HGCo(25))
       S(14)=-6.25d-2*PiZ*(HGCo(22)-HGCo(24))
       W(1)=3.d0*HGCo(42)+HGCo(44)+3.d0*HGCo(46)
       W(2)=-4.d0*HGCo(51)-4.d0*HGCo(53)+8.d0*HGCo(56)
       C(15)=1.25d-1*PiZ*(W(1)+W(2))
       S(15)=0
       W(1)=HGCo(37)+HGCo(39)+5.d0*HGCo(41)
       W(2)=-2.d0*HGCo(48)-6.d0*HGCo(50)+8.d0*HGCo(55)
       C(16)=-6.25d-2*PiZ*(W(1)+W(2))
       W(1)=5.d0*HGCo(36)+HGCo(38)+HGCo(40)
       W(2)=-6.d0*HGCo(47)-2.d0*HGCo(49)+8.d0*HGCo(54)
       S(16)=-6.25d-2*PiZ*(W(1)+W(2))
       C(17)=2.5d-1*PiZ*(HGCo(42)-HGCo(46)-HGCo(51)+HGCo(53))
       S(17)=-1.25d-1*PiZ*(HGCo(43)+HGCo(45)-2.d0*HGCo(52))
       W(1)=3.d0*HGCo(37)+HGCo(39)-5.d0*HGCo(41)-4.d0*HGCo(48)+4.d0*HGCo(50)
       C(18)=-3.125d-2*PiZ*W(1)
       W(1)=5.d0*HGCo(36)-HGCo(38)
       W(2)=-3.d0*HGCo(40)-4.d0*HGCo(47)+4.d0*HGCo(49)
       S(18)=-3.125d-2*PiZ*(W(1)+W(2))
       C(19)=6.25d-2*PiZ*(HGCo(42)-HGCo(44)+HGCo(46))
       S(19)=-6.25d-2*PiZ*(HGCo(43)-HGCo(45))
       C(20)=-3.125d-2*PiZ*(HGCo(37)-HGCo(39)+HGCo(41))
       S(20)=-3.125d-2*PiZ*(HGCo(36)-HGCo(38)+HGCo(40))
       W(1)=5.d0*HGCo(57)+HGCo(59)+HGCo(61)+5.d0*HGCo(63)-6.d0*HGCo(70)
       W(2)=-2.d0*HGCo(72)-6.d0*HGCo(74)
       W(3)=8.d0*HGCo(79)+8.d0*HGCo(81)-1.6d1*HGCo(84)
       C(21)=-6.25d-2*PiZ*(W(1)+W(2)+W(3))
       S(21)=0
       W(1)=HGCo(65)+HGCo(67)+5.d0*HGCo(69)
       W(2)=-2.d0*HGCo(76)-6.d0*HGCo(78)+8.d0*HGCo(83)
       C(22)=-6.25d-2*PiZ*(W(1)+W(2))
       W(1)=5.d0*HGCo(64)+HGCo(66)+HGCo(68)
       W(2)=-6.d0*HGCo(75)-2.d0*HGCo(77)+8.d0*HGCo(82)
       S(22)=-6.25d-2*PiZ*(W(1)+W(2))
       W(1)=1.5d1*HGCo(57)+HGCo(59)-HGCo(61)-1.5d1*HGCo(63)
       W(2)=-1.6d1*HGCo(70)+1.6d1*HGCo(74)+1.6d1*HGCo(79)-1.6d1*HGCo(81)
       C(23)=-1.5625d-2*PiZ*(W(1)+W(2))
       W(1)=5.d0*HGCo(58)+3.d0*HGCo(60)+5.d0*HGCo(62)
       W(2)=-8.d0*HGCo(71)-8.d0*HGCo(73)+1.6d1*HGCo(80)
       S(23)=1.5625d-2*PiZ*(W(1)+W(2))
       W(1)=3.d0*HGCo(65)+HGCo(67)-5.d0*HGCo(69)-4.d0*HGCo(76)+4.d0*HGCo(78)
       C(24)=-3.125d-2*PiZ*W(1)
       W(1)=5.d0*HGCo(64)-HGCo(66)
       W(2)=-3.d0*HGCo(68)-4.d0*HGCo(75)+4.d0*HGCo(77)
       S(24)=-3.125d-2*PiZ*(W(1)+W(2))
       W(1)=3.d0*HGCo(57)-HGCo(59)-HGCo(61)
       W(2)=3.d0*HGCo(63)-2.d0*HGCo(70)+2.d0*HGCo(72)-2.d0*HGCo(74)
       C(25)=-3.125d-2*PiZ*(W(1)+W(2))
       S(25)=6.25d-2*PiZ*(HGCo(58)-HGCo(62)-HGCo(71)+HGCo(73))
       C(26)=-3.125d-2*PiZ*(HGCo(65)-HGCo(67)+HGCo(69))
       S(26)=-3.125d-2*PiZ*(HGCo(64)-HGCo(66)+HGCo(68))
       C(27)=-1.5625d-2*PiZ*(HGCo(57)-HGCo(59)+HGCo(61)-HGCo(63))
       S(27)=1.5625d-2*PiZ*(HGCo(58)-HGCo(60)+HGCo(62))
       W(1)=5.d0*HGCo(93)+HGCo(95)+HGCo(97)+5.d0*HGCo(99)-6.d0*HGCo(106)
       W(2)=-2.d0*HGCo(108)-6.d0*HGCo(110)
       W(3)=8.d0*HGCo(115)+8.d0*HGCo(117)-1.6d1*HGCo(120)
       C(28)=-6.25d-2*PiZ*(W(1)+W(2)+W(3))
       S(28)=0
       W(1)=5.d0*HGCo(86)+3.d0*HGCo(88)
       W(2)=5.d0*HGCo(90)+3.5d1*HGCo(92)-8.d0*HGCo(101)
       W(3)=-8.d0*HGCo(103)-4.d1*HGCo(105)
       W(4)=1.6d1*HGCo(112)+4.8d1*HGCo(114)-6.4d1*HGCo(119)
       C(29)=7.8125d-3*PiZ*(W(1)+W(2)+W(3)+W(4))
       W(1)=3.5d1*HGCo(85)+5.d0*HGCo(87)
       W(2)=3.d0*HGCo(89)+5.d0*HGCo(91)-4.d1*HGCo(100)
       W(3)=-8.d0*HGCo(102)-8.d0*HGCo(104)
       W(4)=4.8d1*HGCo(111)+1.6d1*HGCo(113)-6.4d1*HGCo(118)
       S(29)=7.8125d-3*PiZ*(W(1)+W(2)+W(3)+W(4))
       W(1)=1.5d1*HGCo(93)+HGCo(95)-HGCo(97)-1.5d1*HGCo(99)
       W(2)=-1.6d1*HGCo(106)+1.6d1*HGCo(110)+1.6d1*HGCo(115)-1.6d1*HGCo(117)
       C(30)=-1.5625d-2*PiZ*(W(1)+W(2))
       W(1)=5.d0*HGCo(94)+3.d0*HGCo(96)+5.d0*HGCo(98)
       W(2)=-8.d0*HGCo(107)-8.d0*HGCo(109)+1.6d1*HGCo(116)
       S(30)=1.5625d-2*PiZ*(W(1)+W(2))
       W(1)=9.d0*HGCo(86)+3.d0*HGCo(88)+HGCo(90)-2.1d1*HGCo(92)
       W(2)=-1.2d1*HGCo(101)-4.d0*HGCo(103)
       W(3)=2.d1*HGCo(105)+1.6d1*HGCo(112)-1.6d1*HGCo(114)
       C(31)=7.8125d-3*PiZ*(W(1)+W(2)+W(3))
       W(1)=2.1d1*HGCo(85)-HGCo(87)-3.d0*HGCo(89)-9.d0*HGCo(91)
       W(2)=-2.d1*HGCo(100)+4.d0*HGCo(102)
       W(3)=1.2d1*HGCo(104)+1.6d1*HGCo(111)-1.6d1*HGCo(113)
       S(31)=7.8125d-3*PiZ*(W(1)+W(2)+W(3))
       W(1)=3.d0*HGCo(93)-HGCo(95)-HGCo(97)
       W(2)=3.d0*HGCo(99)-2.d0*HGCo(106)+2.d0*HGCo(108)-2.d0*HGCo(110)
       C(32)=-3.125d-2*PiZ*(W(1)+W(2))
       S(32)=6.25d-2*PiZ*(HGCo(94)-HGCo(98)-HGCo(107)+HGCo(109))
       W(1)=5.d0*HGCo(86)-HGCo(88)-3.d0*HGCo(90)
       W(2)=7.d0*HGCo(92)-4.d0*HGCo(101)+4.d0*HGCo(103)-4.d0*HGCo(105)
       C(33)=7.8125d-3*PiZ*(W(1)+W(2))
       W(1)=7.d0*HGCo(85)-3.d0*HGCo(87)-HGCo(89)
       W(2)=5.d0*HGCo(91)-4.d0*HGCo(100)+4.d0*HGCo(102)-4.d0*HGCo(104)
       S(33)=7.8125d-3*PiZ*(W(1)+W(2))
       C(34)=-1.5625d-2*PiZ*(HGCo(93)-HGCo(95)+HGCo(97)-HGCo(99))
       S(34)=1.5625d-2*PiZ*(HGCo(94)-HGCo(96)+HGCo(98))
       C(35)=7.8125d-3*PiZ*(HGCo(86)-HGCo(88)+HGCo(90)-HGCo(92))
       S(35)=7.8125d-3*PiZ*(HGCo(85)-HGCo(87)+HGCo(89)-HGCo(91))
    CASE DEFAULT
       CALL Halt('Bad logic in HGToSP_Direct, time to remake HGToSP.Inc')
    END SELECT
  END SUBROUTINE HGToSP_Direct
  !====================================================================================
  SUBROUTINE HGToSP_Gen(L,PiZ,HGCo,C,S)
    INTEGER                    :: L
    REAL(DOUBLE)               :: PiZ
    REAL(DOUBLE),DIMENSION(1:) :: HGCo
    REAL(DOUBLE),DIMENSION(0:) :: C,S
    REAL(DOUBLE),DIMENSION(20) :: W
    SELECT CASE(L)
    CASE(0)
       C(0)=PiZ*HGCo(1)
       S(0)=0
    CASE(1)
       C(0)=PiZ*HGCo(1)
       S(0)=0
       C(1)=PiZ*HGCo(4)
       S(1)=0
       C(2)=-5.d-1*PiZ*HGCo(3)
       S(2)=-5.d-1*PiZ*HGCo(2)
    CASE(2)
       C(0)=PiZ*HGCo(1)
       S(0)=0
       C(1)=PiZ*HGCo(4)
       S(1)=0
       C(2)=-5.d-1*PiZ*HGCo(3)
       S(2)=-5.d-1*PiZ*HGCo(2)
       C(3)=-5.d-1*PiZ*(HGCo(5)+HGCo(7)-2.d0*HGCo(10))
       S(3)=0
       C(4)=-5.d-1*PiZ*HGCo(9)
       S(4)=-5.d-1*PiZ*HGCo(8)
       C(5)=-2.5d-1*PiZ*(HGCo(5)-HGCo(7))
       S(5)=2.5d-1*PiZ*HGCo(6)
    CASE(3)
       C(0)=PiZ*HGCo(1)
       S(0)=0
       C(1)=PiZ*HGCo(4)
       S(1)=0
       C(2)=-5.d-1*PiZ*HGCo(3)
       S(2)=-5.d-1*PiZ*HGCo(2)
       C(3)=-5.d-1*PiZ*(HGCo(5)+HGCo(7)-2.d0*HGCo(10))
       S(3)=0
       C(4)=-5.d-1*PiZ*HGCo(9)
       S(4)=-5.d-1*PiZ*HGCo(8)
       C(5)=-2.5d-1*PiZ*(HGCo(5)-HGCo(7))
       S(5)=2.5d-1*PiZ*HGCo(6)
       C(6)=-5.d-1*PiZ*(HGCo(15)+HGCo(17)-2.d0*HGCo(20))
       S(6)=0
       C(7)=1.25d-1*PiZ*(HGCo(12)+3.d0*HGCo(14)-4.d0*HGCo(19))
       S(7)=1.25d-1*PiZ*(3.d0*HGCo(11)+HGCo(13)-4.d0*HGCo(18))
       C(8)=-2.5d-1*PiZ*(HGCo(15)-HGCo(17))
       S(8)=2.5d-1*PiZ*HGCo(16)
       C(9)=1.25d-1*PiZ*(HGCo(12)-HGCo(14))
       S(9)=1.25d-1*PiZ*(HGCo(11)-HGCo(13))
    CASE(4)
       C(0)=PiZ*HGCo(1)
       S(0)=0
       C(1)=PiZ*HGCo(4)
       S(1)=0
       C(2)=-5.d-1*PiZ*HGCo(3)
       S(2)=-5.d-1*PiZ*HGCo(2)
       C(3)=-5.d-1*PiZ*(HGCo(5)+HGCo(7)-2.d0*HGCo(10))
       S(3)=0
       C(4)=-5.d-1*PiZ*HGCo(9)
       S(4)=-5.d-1*PiZ*HGCo(8)
       C(5)=-2.5d-1*PiZ*(HGCo(5)-HGCo(7))
       S(5)=2.5d-1*PiZ*HGCo(6)
       C(6)=-5.d-1*PiZ*(HGCo(15)+HGCo(17)-2.d0*HGCo(20))
       S(6)=0
       C(7)=1.25d-1*PiZ*(HGCo(12)+3.d0*HGCo(14)-4.d0*HGCo(19))
       S(7)=1.25d-1*PiZ*(3.d0*HGCo(11)+HGCo(13)-4.d0*HGCo(18))
       C(8)=-2.5d-1*PiZ*(HGCo(15)-HGCo(17))
       S(8)=2.5d-1*PiZ*HGCo(16)
       C(9)=1.25d-1*PiZ*(HGCo(12)-HGCo(14))
       S(9)=1.25d-1*PiZ*(HGCo(11)-HGCo(13))
       W(1)=3.d0*HGCo(21)+HGCo(23)+3.d0*HGCo(25)
       W(2)=-4.d0*HGCo(30)-4.d0*HGCo(32)+8.d0*HGCo(35)
       C(10)=1.25d-1*PiZ*(W(1)+W(2))
       S(10)=0
       C(11)=1.25d-1*PiZ*(HGCo(27)+3.d0*HGCo(29)-4.d0*HGCo(34))
       S(11)=1.25d-1*PiZ*(3.d0*HGCo(26)+HGCo(28)-4.d0*HGCo(33))
       C(12)=2.5d-1*PiZ*(HGCo(21)-HGCo(25)-HGCo(30)+HGCo(32))
       S(12)=-1.25d-1*PiZ*(HGCo(22)+HGCo(24)-2.d0*HGCo(31))
       C(13)=1.25d-1*PiZ*(HGCo(27)-HGCo(29))
       S(13)=1.25d-1*PiZ*(HGCo(26)-HGCo(28))
       C(14)=6.25d-2*PiZ*(HGCo(21)-HGCo(23)+HGCo(25))
       S(14)=-6.25d-2*PiZ*(HGCo(22)-HGCo(24))
    CASE(5)
       C(0)=PiZ*HGCo(1)
       S(0)=0
       C(1)=PiZ*HGCo(4)
       S(1)=0
       C(2)=-5.d-1*PiZ*HGCo(3)
       S(2)=-5.d-1*PiZ*HGCo(2)
       C(3)=-5.d-1*PiZ*(HGCo(5)+HGCo(7)-2.d0*HGCo(10))
       S(3)=0
       C(4)=-5.d-1*PiZ*HGCo(9)
       S(4)=-5.d-1*PiZ*HGCo(8)
       C(5)=-2.5d-1*PiZ*(HGCo(5)-HGCo(7))
       S(5)=2.5d-1*PiZ*HGCo(6)
       C(6)=-5.d-1*PiZ*(HGCo(15)+HGCo(17)-2.d0*HGCo(20))
       S(6)=0
       C(7)=1.25d-1*PiZ*(HGCo(12)+3.d0*HGCo(14)-4.d0*HGCo(19))
       S(7)=1.25d-1*PiZ*(3.d0*HGCo(11)+HGCo(13)-4.d0*HGCo(18))
       C(8)=-2.5d-1*PiZ*(HGCo(15)-HGCo(17))
       S(8)=2.5d-1*PiZ*HGCo(16)
       C(9)=1.25d-1*PiZ*(HGCo(12)-HGCo(14))
       S(9)=1.25d-1*PiZ*(HGCo(11)-HGCo(13))
       W(1)=3.d0*HGCo(21)+HGCo(23)+3.d0*HGCo(25)
       W(2)=-4.d0*HGCo(30)-4.d0*HGCo(32)+8.d0*HGCo(35)
       C(10)=1.25d-1*PiZ*(W(1)+W(2))
       S(10)=0
       C(11)=1.25d-1*PiZ*(HGCo(27)+3.d0*HGCo(29)-4.d0*HGCo(34))
       S(11)=1.25d-1*PiZ*(3.d0*HGCo(26)+HGCo(28)-4.d0*HGCo(33))
       C(12)=2.5d-1*PiZ*(HGCo(21)-HGCo(25)-HGCo(30)+HGCo(32))
       S(12)=-1.25d-1*PiZ*(HGCo(22)+HGCo(24)-2.d0*HGCo(31))
       C(13)=1.25d-1*PiZ*(HGCo(27)-HGCo(29))
       S(13)=1.25d-1*PiZ*(HGCo(26)-HGCo(28))
       C(14)=6.25d-2*PiZ*(HGCo(21)-HGCo(23)+HGCo(25))
       S(14)=-6.25d-2*PiZ*(HGCo(22)-HGCo(24))
       W(1)=3.d0*HGCo(42)+HGCo(44)+3.d0*HGCo(46)
       W(2)=-4.d0*HGCo(51)-4.d0*HGCo(53)+8.d0*HGCo(56)
       C(15)=1.25d-1*PiZ*(W(1)+W(2))
       S(15)=0
       W(1)=HGCo(37)+HGCo(39)+5.d0*HGCo(41)
       W(2)=-2.d0*HGCo(48)-6.d0*HGCo(50)+8.d0*HGCo(55)
       C(16)=-6.25d-2*PiZ*(W(1)+W(2))
       W(1)=5.d0*HGCo(36)+HGCo(38)+HGCo(40)
       W(2)=-6.d0*HGCo(47)-2.d0*HGCo(49)+8.d0*HGCo(54)
       S(16)=-6.25d-2*PiZ*(W(1)+W(2))
       C(17)=2.5d-1*PiZ*(HGCo(42)-HGCo(46)-HGCo(51)+HGCo(53))
       S(17)=-1.25d-1*PiZ*(HGCo(43)+HGCo(45)-2.d0*HGCo(52))
       W(1)=3.d0*HGCo(37)+HGCo(39)-5.d0*HGCo(41)-4.d0*HGCo(48)+4.d0*HGCo(50)
       C(18)=-3.125d-2*PiZ*W(1)
       W(1)=5.d0*HGCo(36)-HGCo(38)
       W(2)=-3.d0*HGCo(40)-4.d0*HGCo(47)+4.d0*HGCo(49)
       S(18)=-3.125d-2*PiZ*(W(1)+W(2))
       C(19)=6.25d-2*PiZ*(HGCo(42)-HGCo(44)+HGCo(46))
       S(19)=-6.25d-2*PiZ*(HGCo(43)-HGCo(45))
       C(20)=-3.125d-2*PiZ*(HGCo(37)-HGCo(39)+HGCo(41))
       S(20)=-3.125d-2*PiZ*(HGCo(36)-HGCo(38)+HGCo(40))
    CASE(6)
       C(0)=PiZ*HGCo(1)
       S(0)=0
       C(1)=PiZ*HGCo(4)
       S(1)=0
       C(2)=-5.d-1*PiZ*HGCo(3)
       S(2)=-5.d-1*PiZ*HGCo(2)
       C(3)=-5.d-1*PiZ*(HGCo(5)+HGCo(7)-2.d0*HGCo(10))
       S(3)=0
       C(4)=-5.d-1*PiZ*HGCo(9)
       S(4)=-5.d-1*PiZ*HGCo(8)
       C(5)=-2.5d-1*PiZ*(HGCo(5)-HGCo(7))
       S(5)=2.5d-1*PiZ*HGCo(6)
       C(6)=-5.d-1*PiZ*(HGCo(15)+HGCo(17)-2.d0*HGCo(20))
       S(6)=0
       C(7)=1.25d-1*PiZ*(HGCo(12)+3.d0*HGCo(14)-4.d0*HGCo(19))
       S(7)=1.25d-1*PiZ*(3.d0*HGCo(11)+HGCo(13)-4.d0*HGCo(18))
       C(8)=-2.5d-1*PiZ*(HGCo(15)-HGCo(17))
       S(8)=2.5d-1*PiZ*HGCo(16)
       C(9)=1.25d-1*PiZ*(HGCo(12)-HGCo(14))
       S(9)=1.25d-1*PiZ*(HGCo(11)-HGCo(13))
       W(1)=3.d0*HGCo(21)+HGCo(23)+3.d0*HGCo(25)
       W(2)=-4.d0*HGCo(30)-4.d0*HGCo(32)+8.d0*HGCo(35)
       C(10)=1.25d-1*PiZ*(W(1)+W(2))
       S(10)=0
       C(11)=1.25d-1*PiZ*(HGCo(27)+3.d0*HGCo(29)-4.d0*HGCo(34))
       S(11)=1.25d-1*PiZ*(3.d0*HGCo(26)+HGCo(28)-4.d0*HGCo(33))
       C(12)=2.5d-1*PiZ*(HGCo(21)-HGCo(25)-HGCo(30)+HGCo(32))
       S(12)=-1.25d-1*PiZ*(HGCo(22)+HGCo(24)-2.d0*HGCo(31))
       C(13)=1.25d-1*PiZ*(HGCo(27)-HGCo(29))
       S(13)=1.25d-1*PiZ*(HGCo(26)-HGCo(28))
       C(14)=6.25d-2*PiZ*(HGCo(21)-HGCo(23)+HGCo(25))
       S(14)=-6.25d-2*PiZ*(HGCo(22)-HGCo(24))
       W(1)=3.d0*HGCo(42)+HGCo(44)+3.d0*HGCo(46)
       W(2)=-4.d0*HGCo(51)-4.d0*HGCo(53)+8.d0*HGCo(56)
       C(15)=1.25d-1*PiZ*(W(1)+W(2))
       S(15)=0
       W(1)=HGCo(37)+HGCo(39)+5.d0*HGCo(41)
       W(2)=-2.d0*HGCo(48)-6.d0*HGCo(50)+8.d0*HGCo(55)
       C(16)=-6.25d-2*PiZ*(W(1)+W(2))
       W(1)=5.d0*HGCo(36)+HGCo(38)+HGCo(40)
       W(2)=-6.d0*HGCo(47)-2.d0*HGCo(49)+8.d0*HGCo(54)
       S(16)=-6.25d-2*PiZ*(W(1)+W(2))
       C(17)=2.5d-1*PiZ*(HGCo(42)-HGCo(46)-HGCo(51)+HGCo(53))
       S(17)=-1.25d-1*PiZ*(HGCo(43)+HGCo(45)-2.d0*HGCo(52))
       W(1)=3.d0*HGCo(37)+HGCo(39)-5.d0*HGCo(41)-4.d0*HGCo(48)+4.d0*HGCo(50)
       C(18)=-3.125d-2*PiZ*W(1)
       W(1)=5.d0*HGCo(36)-HGCo(38)
       W(2)=-3.d0*HGCo(40)-4.d0*HGCo(47)+4.d0*HGCo(49)
       S(18)=-3.125d-2*PiZ*(W(1)+W(2))
       C(19)=6.25d-2*PiZ*(HGCo(42)-HGCo(44)+HGCo(46))
       S(19)=-6.25d-2*PiZ*(HGCo(43)-HGCo(45))
       C(20)=-3.125d-2*PiZ*(HGCo(37)-HGCo(39)+HGCo(41))
       S(20)=-3.125d-2*PiZ*(HGCo(36)-HGCo(38)+HGCo(40))
       W(1)=5.d0*HGCo(57)+HGCo(59)+HGCo(61)+5.d0*HGCo(63)-6.d0*HGCo(70)
       W(2)=-2.d0*HGCo(72)-6.d0*HGCo(74)
       W(3)=8.d0*HGCo(79)+8.d0*HGCo(81)-1.6d1*HGCo(84)
       C(21)=-6.25d-2*PiZ*(W(1)+W(2)+W(3))
       S(21)=0
       W(1)=HGCo(65)+HGCo(67)+5.d0*HGCo(69)
       W(2)=-2.d0*HGCo(76)-6.d0*HGCo(78)+8.d0*HGCo(83)
       C(22)=-6.25d-2*PiZ*(W(1)+W(2))
       W(1)=5.d0*HGCo(64)+HGCo(66)+HGCo(68)
       W(2)=-6.d0*HGCo(75)-2.d0*HGCo(77)+8.d0*HGCo(82)
       S(22)=-6.25d-2*PiZ*(W(1)+W(2))
       W(1)=1.5d1*HGCo(57)+HGCo(59)-HGCo(61)-1.5d1*HGCo(63)
       W(2)=-1.6d1*HGCo(70)+1.6d1*HGCo(74)+1.6d1*HGCo(79)-1.6d1*HGCo(81)
       C(23)=-1.5625d-2*PiZ*(W(1)+W(2))
       W(1)=5.d0*HGCo(58)+3.d0*HGCo(60)+5.d0*HGCo(62)
       W(2)=-8.d0*HGCo(71)-8.d0*HGCo(73)+1.6d1*HGCo(80)
       S(23)=1.5625d-2*PiZ*(W(1)+W(2))
       W(1)=3.d0*HGCo(65)+HGCo(67)-5.d0*HGCo(69)-4.d0*HGCo(76)+4.d0*HGCo(78)
       C(24)=-3.125d-2*PiZ*W(1)
       W(1)=5.d0*HGCo(64)-HGCo(66)
       W(2)=-3.d0*HGCo(68)-4.d0*HGCo(75)+4.d0*HGCo(77)
       S(24)=-3.125d-2*PiZ*(W(1)+W(2))
       W(1)=3.d0*HGCo(57)-HGCo(59)-HGCo(61)
       W(2)=3.d0*HGCo(63)-2.d0*HGCo(70)+2.d0*HGCo(72)-2.d0*HGCo(74)
       C(25)=-3.125d-2*PiZ*(W(1)+W(2))
       S(25)=6.25d-2*PiZ*(HGCo(58)-HGCo(62)-HGCo(71)+HGCo(73))
       C(26)=-3.125d-2*PiZ*(HGCo(65)-HGCo(67)+HGCo(69))
       S(26)=-3.125d-2*PiZ*(HGCo(64)-HGCo(66)+HGCo(68))
       C(27)=-1.5625d-2*PiZ*(HGCo(57)-HGCo(59)+HGCo(61)-HGCo(63))
       S(27)=1.5625d-2*PiZ*(HGCo(58)-HGCo(60)+HGCo(62))
    CASE(7)
       C(0)=PiZ*HGCo(1)
       S(0)=0
       C(1)=PiZ*HGCo(4)
       S(1)=0
       C(2)=-5.d-1*PiZ*HGCo(3)
       S(2)=-5.d-1*PiZ*HGCo(2)
       C(3)=-5.d-1*PiZ*(HGCo(5)+HGCo(7)-2.d0*HGCo(10))
       S(3)=0
       C(4)=-5.d-1*PiZ*HGCo(9)
       S(4)=-5.d-1*PiZ*HGCo(8)
       C(5)=-2.5d-1*PiZ*(HGCo(5)-HGCo(7))
       S(5)=2.5d-1*PiZ*HGCo(6)
       C(6)=-5.d-1*PiZ*(HGCo(15)+HGCo(17)-2.d0*HGCo(20))
       S(6)=0
       C(7)=1.25d-1*PiZ*(HGCo(12)+3.d0*HGCo(14)-4.d0*HGCo(19))
       S(7)=1.25d-1*PiZ*(3.d0*HGCo(11)+HGCo(13)-4.d0*HGCo(18))
       C(8)=-2.5d-1*PiZ*(HGCo(15)-HGCo(17))
       S(8)=2.5d-1*PiZ*HGCo(16)
       C(9)=1.25d-1*PiZ*(HGCo(12)-HGCo(14))
       S(9)=1.25d-1*PiZ*(HGCo(11)-HGCo(13))
       W(1)=3.d0*HGCo(21)+HGCo(23)+3.d0*HGCo(25)
       W(2)=-4.d0*HGCo(30)-4.d0*HGCo(32)+8.d0*HGCo(35)
       C(10)=1.25d-1*PiZ*(W(1)+W(2))
       S(10)=0
       C(11)=1.25d-1*PiZ*(HGCo(27)+3.d0*HGCo(29)-4.d0*HGCo(34))
       S(11)=1.25d-1*PiZ*(3.d0*HGCo(26)+HGCo(28)-4.d0*HGCo(33))
       C(12)=2.5d-1*PiZ*(HGCo(21)-HGCo(25)-HGCo(30)+HGCo(32))
       S(12)=-1.25d-1*PiZ*(HGCo(22)+HGCo(24)-2.d0*HGCo(31))
       C(13)=1.25d-1*PiZ*(HGCo(27)-HGCo(29))
       S(13)=1.25d-1*PiZ*(HGCo(26)-HGCo(28))
       C(14)=6.25d-2*PiZ*(HGCo(21)-HGCo(23)+HGCo(25))
       S(14)=-6.25d-2*PiZ*(HGCo(22)-HGCo(24))
       W(1)=3.d0*HGCo(42)+HGCo(44)+3.d0*HGCo(46)
       W(2)=-4.d0*HGCo(51)-4.d0*HGCo(53)+8.d0*HGCo(56)
       C(15)=1.25d-1*PiZ*(W(1)+W(2))
       S(15)=0
       W(1)=HGCo(37)+HGCo(39)+5.d0*HGCo(41)
       W(2)=-2.d0*HGCo(48)-6.d0*HGCo(50)+8.d0*HGCo(55)
       C(16)=-6.25d-2*PiZ*(W(1)+W(2))
       W(1)=5.d0*HGCo(36)+HGCo(38)+HGCo(40)
       W(2)=-6.d0*HGCo(47)-2.d0*HGCo(49)+8.d0*HGCo(54)
       S(16)=-6.25d-2*PiZ*(W(1)+W(2))
       C(17)=2.5d-1*PiZ*(HGCo(42)-HGCo(46)-HGCo(51)+HGCo(53))
       S(17)=-1.25d-1*PiZ*(HGCo(43)+HGCo(45)-2.d0*HGCo(52))
       W(1)=3.d0*HGCo(37)+HGCo(39)-5.d0*HGCo(41)-4.d0*HGCo(48)+4.d0*HGCo(50)
       C(18)=-3.125d-2*PiZ*W(1)
       W(1)=5.d0*HGCo(36)-HGCo(38)
       W(2)=-3.d0*HGCo(40)-4.d0*HGCo(47)+4.d0*HGCo(49)
       S(18)=-3.125d-2*PiZ*(W(1)+W(2))
       C(19)=6.25d-2*PiZ*(HGCo(42)-HGCo(44)+HGCo(46))
       S(19)=-6.25d-2*PiZ*(HGCo(43)-HGCo(45))
       C(20)=-3.125d-2*PiZ*(HGCo(37)-HGCo(39)+HGCo(41))
       S(20)=-3.125d-2*PiZ*(HGCo(36)-HGCo(38)+HGCo(40))
       W(1)=5.d0*HGCo(57)+HGCo(59)+HGCo(61)+5.d0*HGCo(63)-6.d0*HGCo(70)
       W(2)=-2.d0*HGCo(72)-6.d0*HGCo(74)
       W(3)=8.d0*HGCo(79)+8.d0*HGCo(81)-1.6d1*HGCo(84)
       C(21)=-6.25d-2*PiZ*(W(1)+W(2)+W(3))
       S(21)=0
       W(1)=HGCo(65)+HGCo(67)+5.d0*HGCo(69)
       W(2)=-2.d0*HGCo(76)-6.d0*HGCo(78)+8.d0*HGCo(83)
       C(22)=-6.25d-2*PiZ*(W(1)+W(2))
       W(1)=5.d0*HGCo(64)+HGCo(66)+HGCo(68)
       W(2)=-6.d0*HGCo(75)-2.d0*HGCo(77)+8.d0*HGCo(82)
       S(22)=-6.25d-2*PiZ*(W(1)+W(2))
       W(1)=1.5d1*HGCo(57)+HGCo(59)-HGCo(61)-1.5d1*HGCo(63)
       W(2)=-1.6d1*HGCo(70)+1.6d1*HGCo(74)+1.6d1*HGCo(79)-1.6d1*HGCo(81)
       C(23)=-1.5625d-2*PiZ*(W(1)+W(2))
       W(1)=5.d0*HGCo(58)+3.d0*HGCo(60)+5.d0*HGCo(62)
       W(2)=-8.d0*HGCo(71)-8.d0*HGCo(73)+1.6d1*HGCo(80)
       S(23)=1.5625d-2*PiZ*(W(1)+W(2))
       W(1)=3.d0*HGCo(65)+HGCo(67)-5.d0*HGCo(69)-4.d0*HGCo(76)+4.d0*HGCo(78)
       C(24)=-3.125d-2*PiZ*W(1)
       W(1)=5.d0*HGCo(64)-HGCo(66)
       W(2)=-3.d0*HGCo(68)-4.d0*HGCo(75)+4.d0*HGCo(77)
       S(24)=-3.125d-2*PiZ*(W(1)+W(2))
       W(1)=3.d0*HGCo(57)-HGCo(59)-HGCo(61)
       W(2)=3.d0*HGCo(63)-2.d0*HGCo(70)+2.d0*HGCo(72)-2.d0*HGCo(74)
       C(25)=-3.125d-2*PiZ*(W(1)+W(2))
       S(25)=6.25d-2*PiZ*(HGCo(58)-HGCo(62)-HGCo(71)+HGCo(73))
       C(26)=-3.125d-2*PiZ*(HGCo(65)-HGCo(67)+HGCo(69))
       S(26)=-3.125d-2*PiZ*(HGCo(64)-HGCo(66)+HGCo(68))
       C(27)=-1.5625d-2*PiZ*(HGCo(57)-HGCo(59)+HGCo(61)-HGCo(63))
       S(27)=1.5625d-2*PiZ*(HGCo(58)-HGCo(60)+HGCo(62))
       W(1)=5.d0*HGCo(93)+HGCo(95)+HGCo(97)+5.d0*HGCo(99)-6.d0*HGCo(106)
       W(2)=-2.d0*HGCo(108)-6.d0*HGCo(110)
       W(3)=8.d0*HGCo(115)+8.d0*HGCo(117)-1.6d1*HGCo(120)
       C(28)=-6.25d-2*PiZ*(W(1)+W(2)+W(3))
       S(28)=0
       W(1)=5.d0*HGCo(86)+3.d0*HGCo(88)
       W(2)=5.d0*HGCo(90)+3.5d1*HGCo(92)-8.d0*HGCo(101)
       W(3)=-8.d0*HGCo(103)-4.d1*HGCo(105)
       W(4)=1.6d1*HGCo(112)+4.8d1*HGCo(114)-6.4d1*HGCo(119)
       C(29)=7.8125d-3*PiZ*(W(1)+W(2)+W(3)+W(4))
       W(1)=3.5d1*HGCo(85)+5.d0*HGCo(87)
       W(2)=3.d0*HGCo(89)+5.d0*HGCo(91)-4.d1*HGCo(100)
       W(3)=-8.d0*HGCo(102)-8.d0*HGCo(104)
       W(4)=4.8d1*HGCo(111)+1.6d1*HGCo(113)-6.4d1*HGCo(118)
       S(29)=7.8125d-3*PiZ*(W(1)+W(2)+W(3)+W(4))
       W(1)=1.5d1*HGCo(93)+HGCo(95)-HGCo(97)-1.5d1*HGCo(99)
       W(2)=-1.6d1*HGCo(106)+1.6d1*HGCo(110)+1.6d1*HGCo(115)-1.6d1*HGCo(117)
       C(30)=-1.5625d-2*PiZ*(W(1)+W(2))
       W(1)=5.d0*HGCo(94)+3.d0*HGCo(96)+5.d0*HGCo(98)
       W(2)=-8.d0*HGCo(107)-8.d0*HGCo(109)+1.6d1*HGCo(116)
       S(30)=1.5625d-2*PiZ*(W(1)+W(2))
       W(1)=9.d0*HGCo(86)+3.d0*HGCo(88)+HGCo(90)-2.1d1*HGCo(92)
       W(2)=-1.2d1*HGCo(101)-4.d0*HGCo(103)
       W(3)=2.d1*HGCo(105)+1.6d1*HGCo(112)-1.6d1*HGCo(114)
       C(31)=7.8125d-3*PiZ*(W(1)+W(2)+W(3))
       W(1)=2.1d1*HGCo(85)-HGCo(87)-3.d0*HGCo(89)-9.d0*HGCo(91)
       W(2)=-2.d1*HGCo(100)+4.d0*HGCo(102)
       W(3)=1.2d1*HGCo(104)+1.6d1*HGCo(111)-1.6d1*HGCo(113)
       S(31)=7.8125d-3*PiZ*(W(1)+W(2)+W(3))
       W(1)=3.d0*HGCo(93)-HGCo(95)-HGCo(97)
       W(2)=3.d0*HGCo(99)-2.d0*HGCo(106)+2.d0*HGCo(108)-2.d0*HGCo(110)
       C(32)=-3.125d-2*PiZ*(W(1)+W(2))
       S(32)=6.25d-2*PiZ*(HGCo(94)-HGCo(98)-HGCo(107)+HGCo(109))
       W(1)=5.d0*HGCo(86)-HGCo(88)-3.d0*HGCo(90)
       W(2)=7.d0*HGCo(92)-4.d0*HGCo(101)+4.d0*HGCo(103)-4.d0*HGCo(105)
       C(33)=7.8125d-3*PiZ*(W(1)+W(2))
       W(1)=7.d0*HGCo(85)-3.d0*HGCo(87)-HGCo(89)
       W(2)=5.d0*HGCo(91)-4.d0*HGCo(100)+4.d0*HGCo(102)-4.d0*HGCo(104)
       S(33)=7.8125d-3*PiZ*(W(1)+W(2))
       C(34)=-1.5625d-2*PiZ*(HGCo(93)-HGCo(95)+HGCo(97)-HGCo(99))
       S(34)=1.5625d-2*PiZ*(HGCo(94)-HGCo(96)+HGCo(98))
       C(35)=7.8125d-3*PiZ*(HGCo(86)-HGCo(88)+HGCo(90)-HGCo(92))
       S(35)=7.8125d-3*PiZ*(HGCo(85)-HGCo(87)+HGCo(89)-HGCo(91))
    CASE DEFAULT
       CALL Halt('Bad logic in HGToSP_Gen, time to remake HGToSP.Inc')
    END SELECT
  END SUBROUTINE HGToSP_Gen
  !====================================================================================
  !
  !====================================================================================
  SUBROUTINE MultipoleSetUp
    INTEGER                          :: L,M,LDex,LMDex,   &
         LP,LQ,MP,MQ
    REAL(DOUBLE)                     :: Sgn,DblFact,TwoTimes, &
         DenomP,DenomQ,NumPQ,  &
         DegenP,DegenQ
    !------------------------------------------------------------------------------------
    !         Factorial(M)=M!
    Factorial(0)=One
    DO L=1,FFEll2
       Factorial(L)=Factorial(L-1)*DBLE(L)
    ENDDO
    !         FactOlm2=
    DO L=2,FFEll2
       LDex=L*(L+1)/2
       DO M=0,L-2
          LMDex=LDex+M
          FactOlm2(LMDex)=One/DBLE((L+M)*(L-M))
       ENDDO
    ENDDO
    !         FactMlm2=
    DO L=0,FFEll2
       LDex=L*(L+1)/2
       DO M=0,FFEll2
          LMDex=LDex+M
          FactMlm2(LMDex)=DBLE((L+M-1)*(L-M-1))
       ENDDO
    ENDDO
    !         FactOlm=
    Sgn=One
    DblFact=One
    TwoTimes=One
    DO M=0,FFEll
       FactOlm0(M)=Sgn*DblFact/Factorial(2*M)
       DblFact=DblFact*TwoTimes
       TwoTimes=TwoTimes+Two
       Sgn=-Sgn
    ENDDO
    !         FactMlm=
    Sgn=One
    DblFact=One
    TwoTimes=One
    DO M=0,FFEll2
       FactMlm0(M)=Sgn*DblFact
       DblFact=DblFact*TwoTimes
       TwoTimes=TwoTimes+Two
       Sgn=-Sgn
    ENDDO
    !         FudgeFactorial(LP,LQ)= (LP+LQ)!/(LP! LQ!)
    DO LP=0,SPEll+1
       DO LQ=0,FFELL
          FudgeFactorial(LP,LQ)=Factorial(LP+LQ)/(Factorial(LP)*Factorial(LQ))
       ENDDO
    ENDDO
  END SUBROUTINE MultipoleSetup
  !
  !
  !====================================================================================
  !     Regular Function
  !====================================================================================
  SUBROUTINE Regular(Ell,PQx,PQy,PQz)
    INTEGER                   :: Ell,LenSP
    INTEGER                   :: L,M,M1,M2,MDex,MDex1,LDex,LDex0,LDex1,LDex2,LMDex
    REAL(DOUBLE)              :: PQx2,PQy2,PQxy,PQ,OneOvPQ,CoTan,TwoC,Sq,RS,&
         CoFact,PQToThPlsL,PQx,PQy,PQz
    REAL(DOUBLE) :: RSum
    !------------------------------------------------------------------------------------
    ! Cpq=Zero
    ! Spq=Zero

!!!         RSum=0D0
!!!

    PQx2=PQx*PQx
    PQy2=PQy*PQy
    PQ=SQRT(PQx2+PQy2+PQz*PQz)
    IF(PQ==Zero)THEN
       LenSP=LSP(Ell)
       CALL DBL_VECT_EQ_DBL_SCLR(LenSP+1,Cpq(0),Zero)
       CALL DBL_VECT_EQ_DBL_SCLR(LenSP+1,Spq(0),Zero)
       Cpq(0) = One
       RETURN
    ENDIF
    OneOvPQ=One/PQ
    CoTan=PQz*OneOvPQ
    !        Sine and Cosine by recursion
    Cosine(0)=One
    Sine(  0)=Zero
    PQxy=SQRT(PQx2+PQy2)
    IF(PQxy/=Zero)THEN
       Sine(1)=PQx/PQxy
       Cosine(1)=PQy/PQxy
    ELSE
       Sine(1)  =0.70710678118654752D0
       Cosine(1)=0.70710678118654752D0
    ENDIF
    TwoC=Two*Cosine(1)
    DO M=2,Ell
       M1=M-1
       M2=M-2
       Sine(M)=TwoC*Sine(M1)-Sine(M2)
       Cosine(M)=TwoC*Cosine(M1)-Cosine(M2)
    ENDDO
    !        Associated Legendre Polynomials by recursion
    Sq=SQRT(ABS(One-CoTan**2))
    RS=One
!    WRITE(*,*)' FactO = ',FactOlm0(0)
!    WRITE(*,*)' Cos =  ',Cosine(0)

    DO M=0,Ell
       MDex=LTD(M)+M
       ALegendreP(MDex)=FactOlm0(M)*RS
       RS=RS*Sq
    ENDDO
    DO M=0,Ell-1
       MDex=LTD(M)+M
       MDex1=LTD(M+1)+M
       ALegendreP(MDex1)=CoTan*ALegendreP(MDex)
    ENDDO

    DO L=2,Ell
       CoFact=CoTan*DBLE(2*l-1)
       LDex0=LTD(L)
       LDex1=LTD(L-1)
       LDex2=LTD(L-2)
       DO M=0,L-2
          ALegendreP(LDex0+M)=(CoFact*ALegendreP(LDex1+M)-ALegendreP(LDex2+M))*FactOlm2(LDex0+M)
       ENDDO
    ENDDO
    !        Regular Spharical Harmonics
    PQToThPlsL=One
    DO L=0,Ell
       LDex=LTD(L)
       DO M=0,L
          LMDex=LDex+M
          Spq(LMDex)=PQToThPlsL*ALegendreP(LMDex)*Sine(M)
          Cpq(LMDex)=PQToThPlsL*ALegendreP(LMDex)*Cosine(M)
       ENDDO
       PQToThPlsL=PQToThPlsL*PQ
    ENDDO

!!!101      CONTINUE
!!!         WRITE(*,22)(CPQ(M),M=0,LSP(Ell))
!!! 22      FORMAT('A ',100(D12.6,", "))

!!!         WRITE(*,*)'A ',1,Ell,RSum


  END SUBROUTINE Regular
  !====================================================================================
  !     Irregular Function
  !====================================================================================
  SUBROUTINE IrRegular(Ell,PQx,PQy,PQz,NoR_O,Print_O)
    INTEGER                    :: Ell
    INTEGER                    :: L,M,M1,M2,MDex,MDex1,LDex,LDex0,LDex1,LDex2,LMDex
    REAL(DOUBLE)               :: PQx2,PQy2,PQxy,PQ,OneOvPQ,CoTan,TwoC,Sq,RS,&
         CoFact,PQToThMnsL,PQx,PQy,PQz
    LOGICAL, OPTIONAL :: NoR_O,Print_O
    LOGICAL           :: NoR
    !------------------------------------------------------------------------------------
    IF(PRESENT(NoR_O))THEN
       NoR=NoR_O
    ELSE
       NoR=.FALSE.
    ENDIF

!    IF(PRESENT(Print_O))WRITE(*,*)' R = ',PQx,PQy,PQz

!!!!!!!!!!!! This is fucking dead slow!!!
!!$         Cpq = Zero
!!$         Spq = Zero
    PQx2=PQx*PQx
    PQy2=PQy*PQy
    PQ=SQRT(PQx2+PQy2+PQz*PQz)
    OneOvPQ=One/PQ
    CoTan=PQz*OneOvPQ
    !        Sine and Cosine by recursion
    Cosine(0)=One
    Sine(  0)=Zero
    PQxy=SQRT(PQx2+PQy2)
    IF(PQxy .GT. 1.D-12)THEN
       Sine(1)=PQx/PQxy
       Cosine(1)=PQy/PQxy
    ELSE
       Sine(1)=0.70710678118654752D0
       Cosine(1)=0.70710678118654752D0
    ENDIF

!!$       IF(Present(Print_O))THEN
!!$          WRITE(*,*)' M = ',1,' Cos = ',Cosine(1)
!!$       ENDIF

    TwoC=Two*Cosine(1)
    DO M=2,Ell
       M1=M-1
       M2=M-2
       Sine(M)=TwoC*Sine(M1)-Sine(M2)
       Cosine(M)=TwoC*Cosine(M1)-Cosine(M2)
!!$
!!$
!!$       IF(Present(Print_O))THEN
!!$          WRITE(*,*)' M = ',M,' Cos = ',Cosine(M)
!!$       ENDIF

    ENDDO
    !        Associated Legendre Polynomials by recursion
    Sq=SQRT(ABS(One-CoTan*CoTan))
    RS=One
    DO M=0,Ell
       MDex=LTD(M)+M
       ALegendreP(MDex)=FactMlm0(M)*RS
!!$
!!$       IF(Present(Print_O))THEN
!!$          WRITE(*,*)' M = ',M, 'Fact = ',FactMlm0(M) ,' RS = ',RS,' ALegendre = ',ALegendreP(M)
!!$       ENDIF

       RS=RS*Sq
    ENDDO

    DO M=0,Ell-1
       MDex=LTD(M)+M
       MDex1=LTD(M+1)+M
       ALegendreP(MDex1)=CoTan*DBLE(2*M+1)*ALegendreP(MDex)
    ENDDO
    DO L=2,Ell
       CoFact=CoTan*DBLE(2*l-1)
       LDex0=LTD(L)
       LDex1=LTD(L-1)
       LDex2=LTD(L-2)
       DO M=0,L-2
          ALegendreP(LDex0+M)=CoFact*ALegendreP(LDex1+M)-FactMlm2(LDex0+M)*ALegendreP(LDex2+M)
       ENDDO
    ENDDO
    ! Irregular Spharical Harmonics
    IF(NoR)THEN
       DO L=0,Ell
          LDex=LTD(L)
          DO M=0,L
             LMDex=LDex+M
             Spq(LMDex)=ALegendreP(LMDex)*Sine(M)
             Cpq(LMDex)=ALegendreP(LMDex)*Cosine(M)
          ENDDO
       ENDDO
    ELSE
       PQToThMnsL=OneOvPQ
       DO L=0,Ell
          LDex=LTD(L)
          DO M=0,L
             LMDex=LDex+M
             Spq(LMDex)=PQToThMnsL*ALegendreP(LMDex)*Sine(M)
             Cpq(LMDex)=PQToThMnsL*ALegendreP(LMDex)*Cosine(M)
!!$
!!$       IF(Present(Print_O))THEN
!!$          WRITE(*,33)L,M,PQToThMnsL,ALegendreP(LMDex),Sine(M),PQToThMnsL*ALegendreP(LMDex)*Sine(M)
!!$33        format(2(I3,', ',),4(D12.6,', '))
!!$       ENDIF

          ENDDO
          PQToThMnsL=PQToThMnsL*OneOvPQ
       ENDDO
    ENDIF


  END SUBROUTINE IrRegular
  !====================================================================================
  !     Compute a multipole strength O_L based on Unsolds theorem
  !====================================================================================
  FUNCTION Unsold2(Llow,Lhig,C,S)
    INTEGER                     :: L,Llow,Lhig
    REAL(DOUBLE)                :: Unsold2,U,Exp
    REAL(DOUBLE), DIMENSION(0:) :: C,S
    !
    Unsold2=Zero
    DO L=Llow,Lhig
       IF(L==0) THEN
          U   = Unsold0(L,C,S)
          Unsold2=MAX(Unsold2,U)
       ELSE
          Exp = DBLE(Lhig)/DBLE(L)
          U   = Unsold0(L,C,S)**Exp
          Unsold2=MAX(Unsold2,U)
       ENDIF
    ENDDO
  END FUNCTION Unsold2
  !====================================================================================
  !     Compute a multipole strength O_L based on Unsolds theorem
  !====================================================================================
  FUNCTION Unsold1(Llow,Lhig,C,S)
    INTEGER                     :: Llow,Lhig,LL
    REAL(DOUBLE)                :: Unsold1,U,Exp
    REAL(DOUBLE), DIMENSION(0:) :: C,S
    !
    Unsold1=Zero
    DO LL=Llow,Lhig
       Exp = One/DBLE(MAX(1,LL))
       U   = Unsold0(LL,C,S)**Exp
       Unsold1=MAX(Unsold1,U)
    ENDDO
    !
  END FUNCTION Unsold1
  !====================================================================================
  !     Compute a multipole strength O_L based on Unsolds theorem
  !====================================================================================
  FUNCTION Unsold0(L,C,S)
    INTEGER                     :: I,K,L,M
    REAL(DOUBLE)                :: Unsold0
    REAL(DOUBLE), DIMENSION(0:) :: C,S

    K=LTD(L)
    ! Unsold0=(C(K)**2+S(K)**2)*Factorial(L)**2
    Unsold0=(C(K)*C(K)+S(K)*S(K))*(Factorial(L)*Factorial(L))
    DO M=1,L
       ! Unsold0=Unsold0+Two*(C(K+M)**2+S(K+M)**2)*Factorial(L+M)*Factorial(L-M)
       Unsold0=Unsold0+Two*(C(K+M)*C(K+M)+S(K+M)*S(K+M))*Factorial(L+M)*Factorial(L-M)
    ENDDO
    Unsold0 = SQRT(ABS(Unsold0))

  END FUNCTION Unsold0

  !====================================================================================
  !     Xlate in F90
  !
  !        OO_{l1,m1}[Q] = Sum_{l2,m2} OO_{l1-l2,m1-m2}[Q-P] OO_{l2,m2}[P]
  !
  !        |M1-M2| .LE. |L1-L2|= L3 ==> M2 : MAX(-L2,M1-L3) to MIN(L2,M1+L3)
  !
  !====================================================================================
  SUBROUTINE XLate90(LP,LQ,Cp,Sp,Cpq,Spq,Cq,Sq)
    INTEGER                         :: LP,LQ
    INTEGER                         :: L1,L2,L3,M1,M2,M3,LDX1,LDX2,LDX3,ABSM2,ABSM3
    REAL(DOUBLE)                    :: CN,SN,CMN,SMN
    REAL(DOUBLE),DIMENSION(0:FFLen) :: Cp,Sp,Cq,Sq
    REAL(DOUBLE),DIMENSION(0:FFLen2):: Cpq,Spq
    DO L1 = 0,LP
       DO M1 = 0,L1
          LDX1 = LTD(L1)+M1
          DO L2 = 0,MIN(L1,LQ)
             L3    = L1-L2
             DO M2 = -L2,L2
                ABSM2 = ABS(M2)
                LDX2  = LTD(L2)+ABSM2
                M3    = M1-M2
                ABSM3 = ABS(M3)
                LDX3  = LTD(L3)+ABSM3
                IF(ABSM3 .LE. L3) THEN
                   IF(M2 .LT. 0) THEN
                      CN = (-One)**(ABSM2)
                      SN = -CN
                   ELSE
                      CN = One
                      SN = One
                   ENDIF
                   IF(M3 .LT. 0) THEN
                      CMN = (-One)**(ABSM3)
                      SMN = -CMN
                   ELSE
                      CMN = One
                      SMN = One
                   ENDIF
                   Cp(LDX1) = Cp(LDX1)+CN*CMN*Cq(LDX2)*Cpq(LDX3) &
                        -SN*SMN*Sq(LDX2)*Spq(LDX3)
                   Sp(LDX1) = Sp(LDX1)+CN*SMN*Cq(LDX2)*Spq(LDX3) &
                        +SN*CMN*Sq(LDX2)*Cpq(LDX3)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
  END SUBROUTINE XLate90


#ifdef LDSJFLJSDFLSDJFL
  !====================================================================================
  !
  !====================================================================================
  SUBROUTINE Print_PoleNode(Node,Tag_O)
    TYPE(PoleNode)                 :: Node
    CHARACTER(LEN=*),OPTIONAL      :: Tag_O
    CHARACTER(LEN=DEFAULT_CHR_LEN) :: Line
    INTEGER                        :: L,M,LMDx
    REAL(DOUBLE)                   :: Cpp,Spp
    !
    WRITE(*,*)'======================================================'
    IF(PRESENT(Tag_O))WRITE(*,*)Tag_O
    IF(Node%Leaf)THEN
       Line='Q[#'//TRIM(IntToChar(Node%Box%Number))//',Tier' &
            //TRIM(IntToChar(Node%Box%Tier))//',Leaf] = (' &
            //TRIM(DblToMedmChar(Node%Box%Center(1)))//', '        &
            //TRIM(DblToMedmChar(Node%Box%Center(2)))//', '        &
            //TRIM(DblToMedmChar(Node%Box%Center(3)))//') '
    ELSE
       Line='Q[#'//TRIM(IntToChar(Node%Box%Number))//',Tier' &
            //TRIM(IntToChar(Node%Box%Tier))//'] = (' &
            //TRIM(DblToMedmChar(Node%Box%Center(1)))//', '        &
            //TRIM(DblToMedmChar(Node%Box%Center(2)))//', '        &
            //TRIM(DblToMedmChar(Node%Box%Center(3)))//') '
    ENDIF
    WRITE(*,*)TRIM(Line)
    DO l=0,Node%Ell
       DO m=0,l
          lmdx=LTD(l)+m
          Cpp = Node%C(lmdx)
          Spp = Node%S(lmdx)
          Line='L = '//TRIM(IntToChar(l))//' M = '//TRIM(IntToChar(m)) &
               //' Cq = '//TRIM(DblToMedmChar(Cpp)) &
               //' Sq = '//TRIM(DblToMedmChar(Spp))
          WRITE(*,*)TRIM(Line)
       ENDDO
    ENDDO
  END SUBROUTINE Print_PoleNode
  !========================================================================================
  ! Print Fields
  !========================================================================================
  SUBROUTINE Print_SP(Ell,C,S,Tag_O,Pre_O)
    CHARACTER(LEN=*),OPTIONAL        :: Tag_O,Pre_O
    CHARACTER(LEN=DEFAULT_CHR_LEN)   :: Line
    INTEGER                          :: Ell,L,M,LMDx,IPre
    REAL(DOUBLE)                     :: Cpp,Spp
    REAL(DOUBLE),DIMENSION(0:)       :: C,S
    WRITE(*,*)'======================================================'
    IF(PRESENT(Tag_O)) WRITE(*,*) Tag_O
    IF(PRESENT(Pre_O)) THEN
       IF(Pre_O == 'Short') IPre = 0
       IF(Pre_O == 'Med'  ) IPre = 1
       IF(Pre_O == 'Long' ) IPre = 2
    ELSE
       IPre = 1
    ENDIF
    !
    IF(IPre == 0) THEN
       DO l=0,Ell
          DO m=0,l
             lmdx=LTD(l)+m
             Cpp = C(lmdx)
             Spp = S(lmdx)
             IF(ABS(Cpp) .LT. 1.D-12) Cpp = Zero
             IF(ABS(Spp) .LT. 1.D-12) Spp = Zero
             Line='L = '//TRIM(IntToChar(L))//' M = '//TRIM(IntToChar(M)) &
                  //' Cq = '//TRIM(DblToShrtChar(Cpp)) &
                  //' Sq = '//TRIM(DblToShrtChar(Spp))
             WRITE(*,*) TRIM(Line)
          ENDDO
       ENDDO
    ELSEIF(IPre == 1) THEN
       DO l=0,Ell
          DO m=0,l
             lmdx=LTD(l)+m
             Cpp = C(lmdx)
             Spp = S(lmdx)
             IF(ABS(Cpp) .LT. 1.D-12) Cpp = Zero
             IF(ABS(Spp) .LT. 1.D-12) Spp = Zero
             Line='L = '//TRIM(IntToChar(L))//' M = '//TRIM(IntToChar(M)) &
                  //' Cq = '//TRIM(DblToMedmChar(Cpp)) &
                  //' Sq = '//TRIM(DblToMedmChar(Spp))
             WRITE(*,*)TRIM(Line)
          ENDDO
       ENDDO
    ELSEIF(IPre == 2) THEN
       DO l=0,Ell
          DO m=0,l
             lmdx=LTD(l)+m
             Cpp = C(lmdx)
             Spp = S(lmdx)
             IF(ABS(Cpp) .LT. 1.D-12) Cpp = Zero
             IF(ABS(Spp) .LT. 1.D-12) Spp = Zero
             Line='L = '//TRIM(IntToChar(L))//' M = '//TRIM(IntToChar(M)) &
                  //' Cq = '//TRIM(DblToChar(Cpp)) &
                  //' Sq = '//TRIM(DblToChar(Spp))
             !             IF(.NOT. (Cpp == Zero .AND. Spp == Zero)) THEN
             WRITE(*,*)TRIM(Line)
             !             ENDIF
          ENDDO
       ENDDO
    ENDIF
    !
  END SUBROUTINE Print_SP
#endif
END MODULE MondoPoles
