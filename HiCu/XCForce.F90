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
  TYPE(DBL_RNK2)                 :: LatFrc_XC   
  TYPE(BBox)                     :: WBox,WBoxTmp
  REAL(DOUBLE)                   :: VolRho,VolExc,DelBox,Exc_old,Etot_old,Etot,dum0,dum1
  REAL(DOUBLE),EXTERNAL    :: MondoTimer

!---------------------------------------------------------------------------------------
! Macro the start up
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
! Get basis set, geometry, thresholds and model type
  CALL Get(BS,CurBase)
  CALL Get(GM,CurGeom)
  NEl=GM%NElec
! Set local integration thresholds 
  CALL SetLocalThresholds(Thresholds%Cube*1.D-1)
#ifdef PARALLEL
  CALL ParaInitRho(Args)
  CALL GetBBox()
  CALL SendBBox()
  CALL DistDist()
  CALL ParaRhoToTree()
  CALL ParaGridGen()
#else
! Convert density to a 5-D BinTree
  CALL RhoToTree(Args)
! Generate the grid as a 3-D BinTree 
  WBox%BndBox(1:3,1:2) = RhoRoot%Box%BndBox(1:3,1:2)
! Make Box Periodic
  CALL MakeBoxPeriodic(WBox)
  CALL CalCenterAndHalf(WBox)
  CALL GridGen(WBox,VolRho,VolExc)
! Redue the Energy to reflect a more accurate calculation of Exc
  CALL Get(Exc_old, 'Exc')
  CALL Get(Etot_old,'Etot')
!
  Etot = Etot_old-Exc_old+Exc
  CALL Put(Exc, 'Exc')
  CALL Put(Etot,'Etot')
#endif
! Delete the density
  CALL DeleteRhoTree(RhoRoot)
! More allocations 
  CALL NewBraBlok(BS,Gradients_O=.TRUE.)
  CALL New(XCFrc,3*NAtoms)
  CALL New(LatFrc_XC,(/3,3/))
#ifdef PARALLEL
  CALL New(P,OnAll_O=.TRUE.)
#endif
  CALL Get(P,TrixFile('D',Args,1),BCast_O=.TRUE.)
!----------------------------------------------------------------------
! Compute the exchange-correlation contribution to the force in O(N)
!
  XCFrc%D=Zero
  LatFrc_XC%D=Zero
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
           MN1=MA*NB-1
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
                    F_nlm = dXC(Pair,P%MTrix%D(Q:Q+MN1))
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
  CALL Delete(LCoor)
  CALL Delete(RCoor)
  !
  CALL ParaInitRho(Args)
  CALL GetBBox()
  CALL SendBBox()
  CALL DistDist()
  CALL ParaRhoToTree()
#else
  CALL RhoToTree(Args)
#endif
! Generate the grid as a 3-D BinTree
  DO I = 1,3 
     Exc=Zero
     IF(GM%PBC%AutoW%I(I)==1) THEN
#ifdef PARALLEL
        WBox%BndBox(1:3,1) = LCoor%D(1:3,MyID+1)
        WBox%BndBox(1:3,2) = RCoor%D(1:3,MyID+1)
        CALL CalCenterAndHalf(WBox)
#endif
        WBox%BndBox(I,1) = GM%PBC%BoxShape%D(I,I)-DelBox
        WBox%BndBox(I,2) = GM%PBC%BoxShape%D(I,I)+DelBox
        CALL GridGen(WBox,VolRho,VolExc)
#ifdef PARALLEL
        WBox%BndBox(1:3,1) = LCoor%D(1:3,MyID+1)
        WBox%BndBox(1:3,2) = RCoor%D(1:3,MyID+1)
#else
        WBox%BndBox(I,1) = Zero
        WBox%BndBox(I,2) = GM%PBC%BoxShape%D(I,I)
#endif
        !write(*,'(A,I1,A,I1,A,E26.15,I3)') 'XC_Surface(',i,',',i,')',Exc/(Two*DelBox),MyID
        LatFrc_XC%D(I,I) = LatFrc_XC%D(I,I)+Exc/(Two*DelBox)
     ENDIF
  ENDDO
! Write Out Lattice Force       
!!$  WRITE(*,*) 'LatFrc_XC'
!!$  DO I=1,3
!!$     WRITE(*,*) (LatFrc_XC%D(I,J),J=1,3)
!!$  ENDDO
! Delete the density
  CALL DeleteRhoTree(RhoRoot)

#ifdef PARALLEL
  XCFrcEndTm = MondoTimer()
  XCFrcTm = XCFrcEndTm-XCFrcBegTm
#endif
! Do some checksumming, resumming and IO 
#ifdef PARALLEL
  TotFrcComp = 3*NAtoms
  CALL New(TotXCFrc,TotFrcComp)
  CALL MPI_Reduce(XCFrc%D(1),TotXCFrc%D(1),TotFrcComp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MONDO_COMM,IErr)
  IF(MyID == 0) THEN
    XCFrc%D(1:TotFrcComp) = TotXCFrc%D(1:TotFrcComp)
  ENDIF
  CALL Delete(TotXCFrc)
#endif
  CALL PChkSum(XCFrc,'dXC/dR',Proc_O=Prog)  
! Sum in contribution to total force
  DO AtA=1,NAtoms
     A1=3*(AtA-1)+1
     A2=3*AtA
     GM%Gradients%D(1:3,AtA) =  GM%Gradients%D(1:3,AtA)+XCFrc%D(A1:A2)
  ENDDO
#ifdef PARALLEL
  CALL New(TmpLatFrc_XC,(/3,3/))
  CALL DBL_VECT_EQ_DBL_SCLR(9,TmpLatFrc_XC%D(1,1),0.0d0)
  CALL MPI_REDUCE(LatFrc_XC%D(1,1),TmpLatFrc_XC%D(1,1),9,MPI_DOUBLE_PRECISION, &
       &          MPI_SUM,ROOT,MONDO_COMM,IErr)
  GM%PBC%LatFrc%D = GM%PBC%LatFrc%D+TmpLatFrc_XC%D 
  !if(myid.eq.root) then
  !   do j=1,3
  !      do i=1,3
  !         IF(GM%PBC%AutoW%I(I) == 1 .AND. GM%PBC%AutoW%I(J) == 1) THEN 
  !            write(*,'(A,I1,A,I1,A,E26.15)') 'LatFrc_XC_t(',i,',',j,')=',TmpLatFrc_XC%D(i,j)
  !         ENDIF
  !      enddo
  !   enddo
  !endif
  CALL Delete(TmpLatFrc_XC)
#else
  GM%PBC%LatFrc%D = GM%PBC%LatFrc%D+LatFrc_XC%D
#endif
  CALL Put(GM,Tag_O=CurGeom)
!--------------------------------------------------------------------------------
! Tidy up
!--------------------------------------------------------------------------------
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(P)
  CALL Delete(XCFrc)
  CALL Delete(LatFrc_XC)
  CALL DeleteBraBlok(Gradients_O=.TRUE.)
#ifdef PARALLEL
  CALL New(TmXCFrcArr,NPrc)
  CALL MPI_Gather(XCFrcTm,1,MPI_DOUBLE_PRECISION,TmXCFrcArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
!  IF(MyID == ROOT) THEN
!    CALL PImbalance(TmXCFrcArr,NPrc,Prog_O='XCFrc')
!  ENDIF
  CALL Delete(TmXCFrcArr)
#endif
! didn't count flops, any accumulation is residual
! from matrix routines
  PerfMon%FLOP=Zero 
  CALL ShutDown(Prog)
END PROGRAM XCForce







