!    COMPUTE DERIVATIVES OF THE EXCHANGE CORRELATION ENERGY $E_{xc}$ 
!    WRT TO NUCLEAR COORDINATES IN O(N) USING HIERARCHICAL CUBATURE
!    Author: Matt Challacombe
!------------------------------------------------------------------------------
PROGRAM XCForceDIPMW
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
  USE Functionals
  USE BoundingBox 
  USE DIPMWThresholds
  USE DIPMWTree
  USE KxcGenDIPMW
!
  IMPLICIT NONE
  TYPE(ARGMT)                    :: Args
  TYPE(BCSR)                     :: P
  TYPE(TIME)                     :: TimeRhoToTree,TimeDIPMW,TimeXCForce
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
!---------------------------------------------------------------------------------------
! Macro the start up
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
! Get basis set, geometry, thresholds and model type
  CALL Get(BS,CurBase)
  CALL Get(GM,CurGeom)
  NEl=GM%NElec
! Set local integration thresholds 
! Set local integration thresholds 
  CALL SetLocalThresholdsDIPMW(Thresholds%Cube)
  TauDIPMW = 1.0D-3
  TauRho   = TauDIPMW*1.D-3
  WRITE(*,*) 'TauDIPMW = ',TauDIPMW
  WRITE(*,*) 'TauRho   = ',TauRho
  CALL SetAACoef()
! Convert density to a 5-D BinTree
  CALL RhoToTree(Args)
! Genergate the Wavelet Representation of the XC potential
  CALL DIPMWTreeBuild(20)
! Delete the density
  CALL DeleteRhoTree(RhoRoot)
! More allocations 
  CALL NewBraBlok(BS,Gradients_O=.TRUE.)
  CALL New(XCFrc,3*NAtoms)
  CALL New(LatFrc_XC,  (/3,3/))
  CALL New(LatFrc_XC_S,(/3,3/))
! Get the Density Matrix
!!$  CALL Get(P,TrixFile('D',Args,1),BCast_O=.TRUE.)
!!$!----------------------------------------------------------------------
!!$! Compute the exchange-correlation contribution to the force in O(N)
!!$!
!!$  XCFrc%D      =Zero
!!$  LatFrc_XC%D  =Zero
!!$  LatFrc_XC_S%D=Zero
!!$  DO AtA=1,NAtoms
!!$     MA=BSiz%I(AtA)
!!$     A1=3*(AtA-1)+1
!!$     A2=3*AtA
!!$     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1
!!$        AtB=P%ColPt%I(JP)
!!$        IF(SetAtomPair(GM,BS,AtA,AtB,Pair))THEN
!!$           Q=P%BlkPt%I(JP)
!!$           NB=BSiz%I(AtB)
!!$           MN1=MA*NB*P%NSMat-1 !<<< SPIN
!!$           A=Pair%A
!!$           B=Pair%B
!!$           DO NCA=1,CS_OUT%NCells
!!$              Pair%A=A+CS_OUT%CellCarts%D(:,NCA)
!!$              DO NCB=1,CS_OUT%NCells
!!$                 Pair%B=B+CS_OUT%CellCarts%D(:,NCB)
!!$                 Pair%AB2=(Pair%A(1)-Pair%B(1))**2 &
!!$                         +(Pair%A(2)-Pair%B(2))**2 &
!!$                         +(Pair%A(3)-Pair%B(3))**2
!!$                  IF(TestAtomPair(Pair,CubeRoot%Box)) THEN
!!$!                    F_nlm = dXCDIPMW(Pair,P%MTrix%D(Q:Q+MN1),P%NSMat)
!!$                    IF(Pair%SameAtom) THEN
!!$                       XCFrc%D(A1:A2) = XCFrc%D(A1:A2) + Two*F_nlm(4:6)
!!$                    ELSE
!!$                       XCFrc%D(A1:A2) = XCFrc%D(A1:A2) + Four*F_nlm(1:3)
!!$                    ENDIF
!!$                    nlm = AtomToFrac(GM,Pair%A)
!!$                    LatFrc_XC%D =  LatFrc_XC%D + Two*LaticeForce(GM,nlm,F_nlm(1:3))
!!$!
!!$                    nlm = AtomToFrac(GM,Pair%B)
!!$                    LatFrc_XC%D =  LatFrc_XC%D + Two*LaticeForce(GM,nlm,F_nlm(4:6)-F_nlm(1:3))
!!$                 ENDIF
!!$              ENDDO
!!$           ENDDO
!!$        ENDIF
!!$     ENDDO
!!$  ENDDO
!!$! Zero the Lower Triange
!!$  DO I=1,3
!!$     DO J=1,I-1
!!$        LatFrc_XC%D(I,J)   = Zero
!!$        LatFrc_XC_S%D(I,J) = Zero
!!$     ENDDO
!!$  ENDDO
!!$! Sum in contribution to total force
!!$  DO AtA=1,NAtoms
!!$     A1=3*(AtA-1)+1
!!$     A2=3*AtA
!!$     GM%Gradients%D(1:3,AtA) =  GM%Gradients%D(1:3,AtA)+XCFrc%D(A1:A2)
!!$  ENDDO
! Sum in the J contribution to total lattice force
  LatFrc_XC%D     = LatFrc_XC%D+LatFrc_XC_S%D
  GM%PBC%LatFrc%D = GM%PBC%LatFrc%D+LatFrc_XC%D
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
! didn't count flops, any accumulation is residual
! from matrix routines
  PerfMon%FLOP=Zero 
  CALL ShutDown(Prog)
END PROGRAM XCForceDIPMW







