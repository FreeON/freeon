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
  TYPE(DBL_VECT)      :: TotXCFrc
  INTEGER :: IErr,TotFrcComp
#else
  TYPE(BCSR)                     :: P
#endif

  TYPE(TIME)                     :: TimeRhoToTree,TimeGridGen
  TYPE(DBL_VECT)                 :: XCFrc,Frc
  TYPE(AtomPair)                 :: Pair
  REAL(DOUBLE)                   :: Electrons,XCFrcChk
  INTEGER                        :: AtA,AtB,MA,NB,MN1,A1,A2,JP,Q
  CHARACTER(LEN=3)               :: SCFCycle
  CHARACTER(LEN=7),PARAMETER     :: Prog='XCForce'
  CHARACTER(LEN=15),PARAMETER    :: Sub1='XCForce.RhoTree' 
  CHARACTER(LEN=15),PARAMETER    :: Sub2='XCForce.GridGen' 
  CHARACTER(LEN=DEFAULT_CHR_LEN) :: Mssg 
#ifdef PERIODIC        
  INTEGER                        :: NCA,NCB
  REAL(DOUBLE),DIMENSION(3)      :: A,B,F_nlm,nlm
  REAL(DOUBLE),DIMENSION(3,3)    :: LatFrc_XC
#endif     
  TYPE(BBox)::WBox
  REAL(DOUBLE)::VolRho,VolExc 

!---------------------------------------------------------------------------------------
! Macro the start up
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
! Get basis set, geometry, thresholds and model type
  CALL Get(BS,CurBase)
  CALL Get(GM,CurGeom)
  CALL Get(ModelChem,'ModelChemistry',CurBase)
  NEl=GM%NElec
! Set local integration thresholds 
  CALL SetLocalThresholds(Thresholds%Cube*1.D-1)
#ifdef PERIODIC
! Get the Outer Cell Set
  CALL Get_CellSet(CS_OUT,'CS_OUT'//CurBase//CurGeom)
  CALL PPrint(CS_OUT,'outer sum',Prog)
#endif 

#ifdef PARALLEL
  CALL ParaInitRho(Args)
  CALL GetBBox()
  CALL SendBBox()
  CALL DistDist()
  CALL ParaRhoToTree()
#else
! Convert density to a 5-D BinTree
  CALL RhoToTree(Args)
#endif

#ifdef PARALLEL
  CALL ParaGridGen()
#else
! Generate the grid as a 3-D BinTree 
  WBox%BndBox(1:3,1:2) = RhoRoot%Box%BndBox(1:3,1:2)
#ifdef PERIODIC
  CALL MakeBoxPeriodic(WBox)
#endif
  CALL CalCenterAndHalf(WBox)
  CALL GridGen(WBox,VolRho,VolExc)

#endif


! Delete the density
  CALL DeleteRhoTree(RhoRoot)

! More allocations 
  CALL NewBraBlok(BS,Gradients_O=.TRUE.)
  CALL New(XCFrc,3*NAtoms)
#ifdef PARALLEL
  CALL New(P,OnAll_O=.TRUE.)
#endif
  CALL Get(P,TrixFile('D',Args,1),BCast_O=.TRUE.)
!----------------------------------------------------------------------
! Compute the exchange-correlation contribution to the force in O(N)
!
  XCFrc%D=Zero
  DO AtA=1,NAtoms
     MA=BSiz%I(AtA)
     A1=3*(AtA-1)+1
     A2=3*AtA
     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1
        AtB=P%ColPt%I(JP)
        IF(SetAtomPair(GM,BS,AtA,AtB,Pair))THEN
           Q=P%BlkPt%I(JP)
           NB=BSiz%I(AtB)
           MN1=MA*NB-1
#ifdef PERIODIC
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
                  ! IF(TestAtomPair(Pair)) THEN
                    F_nlm(1:3)     = dXC(Pair,P%MTrix%D(Q:Q+MN1))
                    XCFrc%D(A1:A2) = XCFrc%D(A1:A2) + F_nlm(1:3)
#ifdef THIS_IS_NOT_WELL_POSED_FOR_DIMEN_ZERO
                    nlm = AtomToFrac(GM,CS_OUT%CellCarts%D(1:3,NCA))+AtomToFrac(GM,CS_OUT%CellCarts%D(1:3,NCB))
                    LatFrc_XC(1,1:3) = LatFrc_XC(1,1:3) + Half*nlm(1)*F_nlm(1:3)
                    LatFrc_XC(2,1:3) = LatFrc_XC(2,1:3) + Half*nlm(2)*F_nlm(1:3)
                    LatFrc_XC(3,1:3) = LatFrc_XC(3,1:3) + Half*nlm(3)*F_nlm(1:3)
#endif
                 ENDIF
              ENDDO
           ENDDO
#else
           XCFrc%D(A1:A2)=XCFrc%D(A1:A2)+dXC(Pair,P%MTrix%D(Q:Q+MN1))
#endif
        ENDIF
     ENDDO
  ENDDO
!--------------------------------------------------------------------------------
! Do some checksumming, resumming and IO 

#ifdef PARALLEL
  TotFrcComp = 3*NAtoms
  CALL New(TotXCFrc,TotFrcComp)
  CALL MPI_Reduce(XCFrc%D(1),TotXCFrc%D(1),TotFrcComp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IErr)
  IF(MyID == 0) THEN
    XCFrc%D(1:TotFrcComp) = TotXCFrc%D(1:TotFrcComp)
  ENDIF
  CALL Delete(TotXCFrc)
#endif
! CALL PPrint(XCFrc,'dXC/dR')

  CALL PChkSum(XCFrc,'dXC/dR',Proc_O=Prog)  
! Sum in contribution to total force
  CALL New(Frc,3*NAtoms)
  CALL Get(Frc,'GradE',Tag_O=CurGeom)
  Frc%D=Frc%D+XCFrc%D
  CALL Put(Frc,'GradE',Tag_O=CurGeom)
!--------------------------------------------------------------------------------
! Tidy up
!--------------------------------------------------------------------------------
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(P)
  CALL Delete(Frc)
  CALL Delete(XCFrc)
  CALL DeleteBraBlok(Gradients_O=.TRUE.)
! didn't count flops, any accumulation is residual
! from matrix routines
  PerfMon%FLOP=Zero 
  CALL ShutDown(Prog)
END PROGRAM XCForce
