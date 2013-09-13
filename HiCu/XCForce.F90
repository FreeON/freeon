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
!    COMPUTE DERIVATIVES OF THE EXCHANGE CORRELATION ENERGY $E_{xc}$
!    WRT TO NUCLEAR COORDINATES IN O(N) USING HIERARCHICAL CUBATURE
!    Author: Matt Challacombe
!------------------------------------------------------------------------------
PROGRAM XCForce
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE ProcessControl
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE Thresholding
  USE RhoTree
  USE CubeTree
  USE Functionals
  USE dXCBlok
#ifdef PARALLEL
  USE ParallelHiCu
  USE FastMatrices
#endif
  IMPLICIT NONE
  TYPE(ARGMT)                    :: Args

#ifdef PARALLEL
  TYPE(BCSR)                     :: P
  TYPE(DBL_VECT)                 :: TotXCFrc
  INTEGER                        :: IErr,TotFrcComp
  REAL(DOUBLE)                   :: XCFrcBegTm,XCFrcEndTm,XCFrcTm
  TYPE(DBL_VECT)                 :: TmXCFrcArr
  TYPE(DBL_RNK2)                 :: TmpLatFrc_XC
#else
  TYPE(BCSR)                     :: P
#endif
  TYPE(TIME)                     :: TimeRhoToTree,TimeGridGen
  TYPE(DBL_VECT)                 :: XCFrc,Frc
  TYPE(AtomPair)                 :: Pair
  REAL(DOUBLE)                   :: Electrons,XCFrcChk
  INTEGER                        :: AtA,AtB,MA,NB,MN1,A1,A2,JP,Q,I,J
  CHARACTER(LEN=3)               :: SCFCycle
  CHARACTER(LEN=7),PARAMETER     :: Prog='XCForce'
  CHARACTER(LEN=15),PARAMETER    :: Sub1='XCForce.RhoTree'
  CHARACTER(LEN=15),PARAMETER    :: Sub2='XCForce.GridGen'
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg
  INTEGER                        :: NCA,NCB
  REAL(DOUBLE),DIMENSION(3)      :: A,B,nlm
  REAL(DOUBLE),DIMENSION(6)      :: F_nlm
  TYPE(DBL_RNK2)                 :: LatFrc_XC,LatFrc_XC_S
  TYPE(BBox)                     :: WBox,WBoxTmp
  REAL(DOUBLE)                   :: VolRho,VolExc,DelBox,Exc_old,Etot_old,Etot,dum0,dum1
  LOGICAL                        :: DoingMD
#ifdef PARALLEL
!  REAL(DOUBLE),EXTERNAL    :: MondoTimer
#endif
!---------------------------------------------------------------------------------------
! Macro the start up
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
! Get basis set, geometry, thresholds and model type
  CALL Get(BS,CurBase)
  CALL Get(GM,CurGeom)
  NEl=GM%NElec
! Set local integration thresholds
  CALL Get(DoingMD ,'DoingMD')
  IF(DoingMD) THEN
     CALL SetLocalThresholds(Thresholds%Cube*1.D-2)
  ELSE
     CALL SetLocalThresholds(Thresholds%Cube*1.D-1)
  ENDIF
  CALL SetAACoef()
#ifdef PARALLEL
  CALL ParaInitRho(Args)
  NSDen=Rho%NSDen
!  NSMat=1
!  IF(NSDen.EQ.3) NSMat=2 !<<< SPIN
  write(*,*) '_XCForce_Parallel_: NSDen',NSDen,MyID

  CALL GetBBox()
  CALL SendBBox()
  CALL DistDist()
  CALL ParaRhoToTree()
  CALL ParaGridGen()
#else
! Convert density to a 5-D BinTree
  CALL RhoToTree(Args)

  NSDen=Rho%NSDen
!  NSMat=1
!  IF(NSDen.EQ.3) NSMat=2 !<<< SPIN
  CALL MondoLog(DEBUG_NONE, TRIM(Prog), "NSDen = "//TRIM(IntToChar(NSDen)))

! Generate the grid as a 3-D BinTree
  WBox%BndBox(1:3,1:2) = RhoRoot%Box%BndBox(1:3,1:2)
! Make Box Periodic
  CALL MakeBoxPeriodic(WBox)
  CALL CalCenterAndHalf(WBox)
  CALL GridGen(WBox,VolRho,VolExc)
#endif
! Redue the Energy to reflect a more accurate calculation of Exc
  CALL Get(Exc_old, 'Exc')
  CALL Get(Etot_old,'Etot')
  Etot = Etot_old-Exc_old+Exc
  CALL Put(Exc, 'Exc')
  CALL Put(Etot,'Etot')
! Delete the density
  CALL DeleteRhoTree(RhoRoot)
! More allocations
  CALL NewBraBlok(BS,Gradients_O=.TRUE.)
  CALL New(XCFrc,3*NAtoms)
  CALL New(LatFrc_XC,  (/3,3/))
  CALL New(LatFrc_XC_S,(/3,3/))
#ifdef PARALLEL
  CALL New(P,OnAll_O=.TRUE.)
#endif

  ! Index check...
  CALL MondoLog(DEBUG_NONE, TRIM(Prog), "Index Check: getting P from "//TRIM(TrixFile('D',Args,1)))
  CALL Get(P,TrixFile('D',Args,1),BCast_O=.TRUE.)

!----------------------------------------------------------------------
! Compute the exchange-correlation contribution to the force in O(N)
!
  XCFrc%D      =Zero
  LatFrc_XC%D  =Zero
  LatFrc_XC_S%D=Zero
#ifdef PARALLEL
  XCFrcBegTm = MondoTimer()
#endif
  DO AtA=1,NAtoms
     MA=BSiz%I(AtA)
     A1=3*(AtA-1)+1
     A2=3*AtA
     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1
        AtB=P%ColPt%I(JP)
        IF(SetAtomPair(GM,BS,AtA,AtB,Pair))THEN
           Q=P%BlkPt%I(JP)
           NB=BSiz%I(AtB)
           MN1=MA*NB*P%NSMat-1 !<<< SPIN
           A=Pair%A
           B=Pair%B
           DO NCA=1,CS_OUT%NCells
              Pair%A=A+CS_OUT%CellCarts%D(:,NCA)
              DO NCB=1,CS_OUT%NCells
                 Pair%B=B+CS_OUT%CellCarts%D(:,NCB)
                 Pair%AB2=(Pair%A(1)-Pair%B(1))**2 &
                         +(Pair%A(2)-Pair%B(2))**2 &
                         +(Pair%A(3)-Pair%B(3))**2
                  IF(TestAtomPair(Pair,CubeRoot%Box)) THEN
                    F_nlm = dXC(Pair,P%MTrix%D(Q:Q+MN1),P%NSMat)
                    IF(Pair%SameAtom) THEN
                       XCFrc%D(A1:A2) = XCFrc%D(A1:A2) + Two*F_nlm(4:6)
                    ELSE
                       XCFrc%D(A1:A2) = XCFrc%D(A1:A2) + Four*F_nlm(1:3)
                    ENDIF
                    nlm = AtomToFrac(GM,Pair%A)
                    LatFrc_XC%D =  LatFrc_XC%D + Two*LaticeForce(GM,nlm,F_nlm(1:3))
!
                    nlm = AtomToFrac(GM,Pair%B)
                    LatFrc_XC%D =  LatFrc_XC%D + Two*LaticeForce(GM,nlm,F_nlm(4:6)-F_nlm(1:3))
                 ENDIF
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDDO
!--------------------------------------------------------------------------------
! Calculate the Surface term for the X Force
!--------------------------------------------------------------------------------
! Set the Thresholds
  DelBox = 1.D-4
  CALL SetLocalThresholds(Thresholds%Cube*1.D-4)
! Convert density to a 5-D BinTree
#ifdef PARALLEL
! Convert density to a 5-D BinTree
  CALL Delete(LCoor)
  CALL Delete(RCoor)
  CALL ParaInitRho(Args)
  CALL GetBBox()
  CALL SendBBox()
  CALL DistDist()
  CALL ParaRhoToTree()
  DO I = 1,3
     IF(GM%PBC%AutoW%I(I)==1) THEN
        Exc=Zero
        WBox%BndBox(1:3,1) = LCoor%D(1:3,MyID+1)
        WBox%BndBox(1:3,2) = RCoor%D(1:3,MyID+1)
        IF(WBox%BndBox(I,2)+1.D-8 > GM%PBC%BoxShape%D(I,I))  THEN
           CALL CalCenterAndHalf(WBox)
           WBox%BndBox(I,1) = GM%PBC%BoxShape%D(I,I)-DelBox
           WBox%BndBox(I,2) = GM%PBC%BoxShape%D(I,I)+DelBox
           CALL GridGen(WBox,VolRho,VolExc)
           WBox%BndBox(1:3,1) = LCoor%D(1:3,MyID+1)
           WBox%BndBox(1:3,2) = RCoor%D(1:3,MyID+1)
           LatFrc_XC_S%D(I,I) = LatFrc_XC_S%D(I,I)+Exc/(Two*DelBox)
        ENDIF
     ENDIF
  ENDDO
#else
  CALL RhoToTree(Args)
! Generate the grid as a 3-D BinTree
  DO I = 1,3
     Exc=Zero
     IF(GM%PBC%AutoW%I(I)==1) THEN
        WBox%BndBox(I,1)   = GM%PBC%BoxShape%D(I,I)-DelBox
        WBox%BndBox(I,2)   = GM%PBC%BoxShape%D(I,I)+DelBox
        CALL GridGen(WBox,VolRho,VolExc)
        WBox%BndBox(I,1)   = Zero
        WBox%BndBox(I,2)   = GM%PBC%BoxShape%D(I,I)
        LatFrc_XC_S%D(I,I) = LatFrc_XC_S%D(I,I)+Exc/(Two*DelBox)
     ENDIF
  ENDDO
#endif
! Delete the density
  CALL DeleteRhoTree(RhoRoot)
#ifdef PARALLEL
! Collect the timings
  XCFrcEndTm = MondoTimer()
  XCFrcTm    = XCFrcEndTm-XCFrcBegTm
! Collect the Forces
  TotFrcComp = 3*NAtoms
  CALL New(TotXCFrc,TotFrcComp)
  CALL MPI_Reduce(XCFrc%D(1),TotXCFrc%D(1),TotFrcComp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MONDO_COMM,IErr)
  IF(MyID == ROOT) THEN
    XCFrc%D(1:TotFrcComp) = TotXCFrc%D(1:TotFrcComp)
  ENDIF
  CALL Delete(TotXCFrc)
! Collect the Lattice Forces
  CALL New(TmpLatFrc_XC,(/3,3/))
  CALL DBL_VECT_EQ_DBL_SCLR(9,TmpLatFrc_XC%D(1,1),0.0d0)
  CALL MPI_REDUCE(LatFrc_XC%D(1,1),TmpLatFrc_XC%D(1,1),9,MPI_DOUBLE_PRECISION,MPI_SUM,ROOT,MONDO_COMM,IErr)
  IF(MyID == ROOT) THEN
     LatFrc_XC%D   = TmpLatFrc_XC%D
  ENDIF
  CALL DBL_VECT_EQ_DBL_SCLR(9,TmpLatFrc_XC%D(1,1),0.0d0)
  CALL MPI_REDUCE(LatFrc_XC_S%D(1,1),TmpLatFrc_XC%D(1,1),9,MPI_DOUBLE_PRECISION,MPI_SUM,ROOT,MONDO_COMM,IErr)
  IF(MyID == ROOT) THEN
     LatFrc_XC_S%D = TmpLatFrc_XC%D
  ENDIF
  CALL Delete(TmpLatFrc_XC)
#endif
#ifdef PARALLEL
  IF(MyID == ROOT) THEN
#endif
!    Rescale the Forces if needed.
     IF(P%NSMat.GT.1) CALL DSCAL(3*NAtoms,0.5D0,    XCFrc%D(1  ),1)
     IF(P%NSMat.GT.1) CALL DSCAL(       9,0.5D0,LatFrc_XC%D(1,1),1)
!    Zero the Lower Triange
     DO I=1,3
        DO J=1,I-1
           LatFrc_XC%D(I,J)   = Zero
           LatFrc_XC_S%D(I,J) = Zero
        ENDDO
     ENDDO
!    Sum in contribution to total force
     DO AtA=1,NAtoms
        A1=3*(AtA-1)+1
        A2=3*AtA
        GM%Gradients%D(1:3,AtA) =  GM%Gradients%D(1:3,AtA)+XCFrc%D(A1:A2)
     ENDDO
!    Sum in the J contribution to total lattice force
     LatFrc_XC%D     = LatFrc_XC%D+LatFrc_XC_S%D
     GM%PBC%LatFrc%D = GM%PBC%LatFrc%D+LatFrc_XC%D
#ifdef PARALLEL
  ENDIF
#endif
! Do some printing
  CALL Print_Force(GM,XCFrc,'XC Force')
  CALL Print_Force(GM,XCFrc,'XC Force',Unit_O=6)
  CALL Print_LatForce(GM,LatFrc_XC%D,'XC Lattice Force')
  CALL Print_LatForce(GM,LatFrc_XC%D,'XC Lattice Force',Unit_O=6)
  CALL Print_LatForce(GM,LatFrc_XC_S%D,'XC Dipole Lattice Force: Surface Term')
  CALL Print_LatForce(GM,LatFrc_XC_S%D,'XC Dipole Lattice Force: Surface Term',Unit_O=6)
! Do some checksumming and IO
  CALL PChkSum(XCFrc,    'dXC/dR',Proc_O=Prog)
  CALL PChkSum(LatFrc_XC,'LFrcXC',Proc_O=Prog)
! Save Forces to Disk
  CALL Put(GM,Tag_O=CurGeom)
! Tidy up
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(P)
  CALL Delete(XCFrc)
  CALL Delete(LatFrc_XC)
  CALL Delete(LatFrc_XC_S)
  CALL DeleteBraBlok(Gradients_O=.TRUE.)
#ifdef PARALLEL
  CALL New(TmXCFrcArr,NPrc)
  CALL MPI_Gather(XCFrcTm,1,MPI_DOUBLE_PRECISION,TmXCFrcArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL Delete(TmXCFrcArr)
#endif
! didn't count flops, any accumulation is residual
! from matrix routines
  PerfMon%FLOP=Zero
  CALL ShutDown(Prog)
END PROGRAM XCForce
