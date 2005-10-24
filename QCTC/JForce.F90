!    FAST O(N Lg N) COMPUTATION OF GRADIENTS OF THE COULOMB ENERGY 
!    WRT TO NUCLEAR COORDINATES
!    Author: Matt Challacombe
!==============================================================================
PROGRAM JForce
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE Globals
  USE AtomPairs
  USE BraBloks
  USE PoleTree
  USE PBCFarField
  USE BlokTrPdJ
  USE BlokTrWdS
  USE NuklarE
  USE SetXYZ 
  IMPLICIT NONE

#ifdef PARALLEL
  TYPE(BCSR)                   :: P
  TYPE(DBL_VECT)               :: TotJFrc
  INTEGER                      :: IErr,TotFrcComp,MyAtomNum,TotAtomNum
  REAL(DOUBLE)                 :: JFrcBegTm,JFrcEndTm,JFrcTm
  TYPE(DBL_VECT)               :: TmJFrcArr
  TYPE(DBL_RNK2)               :: TmpLatFrc_J
#else
  TYPE(BCSR)                   :: P
#endif
  TYPE(AtomPair)               :: Pair
  TYPE(DBL_VECT)               :: Frc,JFrc
  INTEGER                      :: AtA,AtB,A1,A2,B1,B2,MA,NB,MN1,MN,JP,Q
  REAL(DOUBLE)                 :: JFrcChk
  CHARACTER(LEN=6),PARAMETER   :: Prog='JForce'
  INTEGER                      :: NC,I,J,K
  REAL(DOUBLE),DIMENSION(3)    :: A,B,nlm
  REAL(DOUBLE),DIMENSION(15)   :: F_nlm
  TYPE(DBL_RNK2)               :: LatFrc_J,LatFrc_J_PFF,LatFrc_J_Dip
  REAL(DOUBLE),DIMENSION(3,3)  :: DivCV
  TYPE(CRDS)                   :: GMLoc
#ifdef PARALLEL
  REAL(DOUBLE),EXTERNAL        :: MondoTimer
#endif
!-------------------------------------------------------------------------------- 
! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GMLoc,Tag_O=CurGeom)
! Allocations 
  CALL New(JFrc,3*GMLoc%Natms)
  JFrc%D(:)=Zero
  CALL NewBraBlok(BS,Gradients_O=.TRUE.)
#ifdef PARALLEL
  CALL New(P,OnAll_O=.TRUE.)
  CALL Get(P,TrixFile('D',Args,1),BCast_O=.TRUE.)
  CALL GetDistrRho('Rho',Args,1)
#else
!  WRITE(*,*) "JForce"
!  WRITE(*,*) "D ",Args%I%I(4)
  CALL Get(P,TrixFile('D',Args,0),BCast_O=.TRUE.)
  CALL Get(Rho,'Rho',Args,1,Bcast_O=.TRUE.)  
#endif
  CALL Get(RhoPoles) 
! Set thresholds local to JForce (for PAC and MAC)
  CALL SetLocalThresholds(Thresholds%TwoE)
! Setup global arrays for computation of multipole tensors
  CALL InitRhoAux
! Setup global arrays for computation of multipole tensors
  CALL MultipoleSetUp()
! Build the global PoleTree representation of the total density
#ifdef PARALLEL
  CALL ParaRhoToPoleTree
#else
  CALL RhoToPoleTree
#endif
! Set the electrostatic background
  CALL PBCFarFieldSetUp(PoleRoot,GMLoc)
! Delete the auxiliary density arrays
  CALL DeleteRhoAux
! Delete the Density
  CALL Delete(Rho)
! Rescale the density matrix for U/G theory.
  IF(P%NSMat.GT.1) CALL DSCAL(P%NNon0,0.5D0,P%MTrix%D(1),1)
!--------------------------------------------------------------------------------
! Compute the Coulomb contribution to the force in O(N Lg N)
!--------------------------------------------------------------------------------
  CALL New(LatFrc_J,(/3,3/))
  LatFrc_J%D = Zero
#ifdef PARALLEL
  MyAtomNum = End%I(MyId)-Beg%I(MyId)+1
  CALL MPI_AllReduce(MyAtomNum,TotAtomNum,1,MPI_INTEGER,MPI_SUM,MONDO_COMM,IErr)
  IF(TotAtomNum /= GMLoc%Natms) THEN
    WRITE(*,*) 'TotAtomNum = ',TotAtomNum
    WRITE(*,*) 'GMLoc%Natms = ',GMLoc%Natms
    STOP 'TotAtomNum not equal to GMLoc%Natms in JForce!'
  ENDIF
  JFrcBegTm = MondoTimer()
  DO AtA=Beg%I(MyId),End%I(MyId)
#else
  DO AtA=1,GMLoc%Natms
#endif
     MA=BSiz%I(AtA)
     A1=3*(AtA-1)+1
     A2=3*AtA
     F_nlm = dNukE(GMLoc,AtA)
     JFrc%D(A1:A2)= Two*F_nlm(1:3) 
!    Store Inner Nuc Lattice Forces
     LatFrc_J%D(1:3,1) = LatFrc_J%D(1:3,1) + F_nlm(7:9)
     LatFrc_J%D(1:3,2) = LatFrc_J%D(1:3,2) + F_nlm(10:12)
     LatFrc_J%D(1:3,3) = LatFrc_J%D(1:3,3) + F_nlm(13:15)
!    Outer Nuc Lattice Forces
     nlm        = AtomToFrac(GMLoc,GMLoc%Carts%D(:,AtA))
     LatFrc_J%D = LatFrc_J%D + Two*LaticeForce(GMLoc,nlm,F_nlm(1:3))
!    Start AtB Loop
     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1 
        AtB=P%ColPt%I(JP)
        B1=3*(AtB-1)+1
        B2=3*AtB
        IF(SetAtomPair(GMLoc,BS,AtA,AtB,Pair)) THEN 
           Q=P%BlkPt%I(JP)
           NB=BSiz%I(AtB)
           MN1=MA*NB-1
           MN=MN1+1
           A=Pair%A
           B=Pair%B
!
           !Quick and dirty!
           SELECT CASE(P%NSMat)
           CASE(1)
              !We don't need to do anything!
           CASE(2)
              !We add up the two density martices!
              CALL DAXPY(MN,1D0,P%MTrix%D(Q+MN),1,P%MTrix%D(Q),1)
           CASE(4)
              !We add up the diagonal density martices!
              CALL DAXPY(MN,1D0,P%MTrix%D(Q+3*MN),1,P%MTrix%D(Q),1)
           CASE DEFAULT;CALL Halt(' JForce: P%NSMat doesn''t have an expected value! ')
           END SELECT
!
           DO NC=1,CS_OUT%NCells 
              Pair%A=A
              Pair%B=B+CS_OUT%CellCarts%D(:,NC)
              Pair%AB2=(Pair%A(1)-Pair%B(1))**2 &
                      +(Pair%A(2)-Pair%B(2))**2 &
                      +(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair))THEN
                 F_nlm = TrPdJ(Pair,P%MTrix%D(Q:Q+MN1),GMLoc)
                 IF(Pair%SameAtom) THEN
                    JFrc%D(A1:A2) = JFrc%D(A1:A2) +  Four*F_nlm(4:6)
                 ELSE
                    JFrc%D(A1:A2) = JFrc%D(A1:A2) + Eight*F_nlm(1:3)
                 ENDIF
!                Store Inner J Lattice Forces
                 LatFrc_J%D(1:3,1) = LatFrc_J%D(1:3,1) + Two*F_nlm(7:9)
                 LatFrc_J%D(1:3,2) = LatFrc_J%D(1:3,2) + Two*F_nlm(10:12)
                 LatFrc_J%D(1:3,3) = LatFrc_J%D(1:3,3) + Two*F_nlm(13:15)
!                Outer Lattice J Forces 
                 nlm        = AtomToFrac(GMLoc,Pair%A) 
                 LatFrc_J%D = LatFrc_J%D+Four*LaticeForce(GMLoc,nlm,F_nlm(1:3))
                 nlm        = AtomToFrac(GMLoc,Pair%B) 
                 LatFrc_J%D = LatFrc_J%D+Four*LaticeForce(GMLoc,nlm,(F_nlm(4:6)-F_nlm(1:3)))
              ENDIF
           ENDDO
!
        ENDIF
     ENDDO
  ENDDO
!
!  A = Zero
!  DO AtA=1,NAtoms
!     A1=3*(AtA-1)+1
!     A2=3*AtA
!     A(1) = A(1) + JFrc%D(A1)
!     A(2) = A(2) + JFrc%D(A1+1)
!     A(3) = A(3) + JFrc%D(A1+2)
!  ENDDO
!  WRITE(*,*) A(1),A(2),A(3)
#ifdef PARALLEL
  IF(MyID == ROOT) THEN
#endif
!    The Dipole Contribution to the Lattice Forces
     CALL New(LatFrc_J_Dip,(/3,3/))
     LatFrc_J_Dip%D = Zero
     IF(GMLoc%PBC%Dimen > 0) THEN
        DivCV   = DivCellVolume(GMLoc%PBC%BoxShape%D,GMLoc%PBC%AutoW%I)
        DO I=1,3
           DO J=1,3
              IF(GMLoc%PBC%AutoW%I(I)==1 .AND. GMLoc%PBC%AutoW%I(J)==1) THEN
                 LatFrc_J_Dip%D(I,J) = LatFrc_J_Dip%D(I,J)  -  E_DP*DivCV(I,J)/GMLoc%PBC%CellVolume
              ENDIF
           ENDDO
        ENDDO
     ENDIF
!    The Farfield Contribution to the Lattice Forces
     CALL New(LatFrc_J_PFF,(/3,3/))
     LatFrc_J_PFF%D = Zero
     IF(GMLoc%PBC%Dimen > 0) THEN
        DO I=1,3
           DO J=1,3
              IF(GMLoc%PBC%AutoW%I(I)==1 .AND. GMLoc%PBC%AutoW%I(J)==1) THEN
                 DO K=0,LSP(MaxEll)
                    LatFrc_J_PFF%D(I,J) = LatFrc_J_PFF%D(I,J) - Two*(RhoC%D(K)*dTenRhoC%D(K,I,J)+RhoS%D(K)*dTenRhoS%D(K,I,J))
                 ENDDO
              ENDIF
           ENDDO
        ENDDO
     ENDIF
#ifdef PARALLEL
  ENDIF
#endif
#ifdef PARALLEL
! Collect the timings
  JFrcEndTm = MondoTimer()
  JFrcTm = JFrcEndTm-JFrcBegTm
! Collect the Forces
  TotFrcComp = 3*GMLoc%Natms
  CALL New(TotJFrc,TotFrcComp)
  CALL MPI_Reduce(JFrc%D(1),TotJFrc%D(1),TotFrcComp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MONDO_COMM,IErr)
  IF(MyID == ROOT) THEN
    JFrc%D(1:TotFrcComp) = TotJFrc%D(1:TotFrcComp)
  ENDIF
  CALL Delete(TotJFrc)
! Collect the Lattice Forces
  CALL New(TmpLatFrc_J,(/3,3/))
  CALL DBL_VECT_EQ_DBL_SCLR(9,TmpLatFrc_J%D(1,1),0.0d0)
  CALL MPI_REDUCE(LatFrc_J%D(1,1),TmpLatFrc_J%D(1,1),9,MPI_DOUBLE_PRECISION,MPI_SUM,ROOT,MONDO_COMM,IErr)
  IF(MyID == ROOT) THEN
     LatFrc_J%D = TmpLatFrc_J%D
  ENDIF
  CALL Delete(TmpLatFrc_J)
#endif
#ifdef PARALLEL
  IF(MyID == ROOT) THEN
#endif
!    Zero the Lower Triange
     DO I=1,3
        DO J=1,I-1
           LatFrc_J%D(I,J)     = Zero
           LatFrc_J_Dip%D(I,J) = Zero
           LatFrc_J_PFF%D(I,J) = Zero
        ENDDO
     ENDDO
!    Do some printing
     CALL Print_LatForce(GMLoc,LatFrc_J%D,'J Lattice Force')
     CALL Print_LatForce(GMLoc,LatFrc_J%D,'J Lattice Force',Unit_O=6)
     CALL Print_LatForce(GMLoc,LatFrc_J_Dip%D,'J Dipole Lattice Force')
     CALL Print_LatForce(GMLoc,LatFrc_J_Dip%D,'J Dipole Lattice Force',Unit_O=6)
     CALL Print_LatForce(GMLoc,LatFrc_J_PFF%D,'J  PFF   Lattice Force')
     CALL Print_LatForce(GMLoc,LatFrc_J_PFF%D,'J  PFF   Lattice Force',Unit_O=6)
!    Sum in the J contribution to total force
     DO AtA=1,NAtoms
        A1=3*(AtA-1)+1
        A2=3*AtA
        GMLoc%Gradients%D(1:3,AtA) = GMLoc%Gradients%D(1:3,AtA)+JFrc%D(A1:A2)
     ENDDO
!    Sum in the J contribution to total lattice force,including Dip and PFF
     LatFrc_J%D         = LatFrc_J%D+LatFrc_J_Dip%D+LatFrc_J_PFF%D
     GMLoc%PBC%LatFrc%D = GMLoc%PBC%LatFrc%D+LatFrc_J%D
!    Tidy Up
     CALL Delete(LatFrc_J_Dip)
     CALL Delete(LatFrc_J_PFF)
#ifdef PARALLEL
  ENDIF
#endif
! Do some printing
  CALL Print_Force(GMLoc,JFrc,'J Force')
  CALL Print_Force(GMLoc,JFrc,'J Force',Unit_O=6)
! Do some checksumming and IO 
  CALL PChkSum(JFrc,    'dJ/dR',Proc_O=Prog)  
  CALL PChkSum(LatFrc_J,'LFrcJ',Proc_O=Prog)  
! Save Forces to Disk
  CALL Put(GMLoc,Tag_O=CurGeom)
! Tidy Up  
  CALL Delete(JFrc)
  CALL Delete(LatFrc_J)
  CALL Delete(P)
  CALL Delete(BS)
  CALL Delete(GMLoc)
  CALL Delete(RhoPoles)
  CALL DeleteBraBlok(Gradients_O=.TRUE.)
#ifdef PARALLEL
  CALL New(TmJFrcArr,NPrc)
  CALL MPI_Gather(JFrcTm,1,MPI_DOUBLE_PRECISION,TmJFrcArr%D(1),1,MPI_DOUBLE_PRECISION,0,MONDO_COMM,IErr)
  CALL Delete(TmJFrcArr)
#endif
! didn't count flops, any accumulation is residual from matrix routines
  PerfMon%FLOP=Zero 
  CALL ShutDown(Prog)
END PROGRAM JForce
