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
  IMPLICIT NONE
  TYPE(ARGMT)                    :: Args
  TYPE(BCSR)                     :: P
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
  REAL(DOUBLE),DIMENSION(3)      :: A,B
#endif     
!---------------------------------------------------------------------------------------
! Macro the start up
  CALL StartUp(Args,Prog,Serial_O=.TRUE.)
! Get basis set, geometry, thresholds and model type
  CALL Get(BS,CurBase)
  CALL Get(GM,CurGeom)
  CALL Get(ModelChem,'ModelChemistry',CurBase)
  NEl=GM%NElec
! Set local integration thresholds 
  CALL SetLocalThresholds(Thresholds%Cube*1.D-1)
#ifdef PERIODIC
! Calculate the Number of Cells
  CALL SetCellNumber(GM)
  CALL PPrint(CS_OUT,'outer sum',Prog)
#endif
! Convert density to a 5-D BinTree
  CALL RhoToTree(Args)
! Generate the grid as a 3-D BinTree 
  CALL GridGen()
! Delete the density
  CALL DeleteRhoTree(RhoRoot)
! More allocations 
  CALL NewBraBlok(BS,Gradients_O=.TRUE.)
  CALL New(XCFrc,3*NAtoms)
  CALL Get(P,TrixFile('D',Args,1))
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
                 IF(TestAtomPair(Pair)) THEN
                    XCFrc%D(A1:A2)=XCFrc%D(A1:A2)+dXC(Pair,P%MTrix%D(Q:Q+MN1))
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
  CALL PPrint(XCFrc,'dXC/dR')
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
