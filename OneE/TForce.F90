!    COMPUTE THE FORCE CORESPONDING TO THE DERIVATIVE OF THE 
!    KINETIC ENERGY MATRIX, TForce=2*Tr{P.dT}
!    Authors: Matt Challacombe and CJ Tymczak
!----------------------------------------------------------------------------------
PROGRAM TForce
  USE DerivedTypes
  USE GlobalScalars
  USE GlobalCharacters
  USE InOut
  USE PrettyPrint
  USE MemMan
  USE Parse
  USE Macros
  USE LinAlg
  USE AtomPairs
  USE BlokTrPdT

#ifdef PARALLEL
  USE MondoMPI
  TYPE(BCSR)         :: P
  TYPE(DBL_VECT)      :: TotTFrc
  INTEGER :: IErr,TotFrcComp
#else
  TYPE(BCSR)          :: P
#endif

#ifdef PERIODIC 
  INTEGER                     :: NC
  REAL(DOUBLE),DIMENSION(3)   :: B,F_nlm,nlm
  REAL(DOUBLE),DIMENSION(3,3) :: LatFrc_T
#endif
  TYPE(AtomPair)      :: Pair
  TYPE(BSET)          :: BS
  TYPE(CRDS)          :: GM
  TYPE(ARGMT)         :: Args
  INTEGER             :: Q,R,AtA,AtB,JP,MB,MA,NB,MN1,A1,A2
  TYPE(HGRho)         :: Rho
  TYPE(DBL_VECT)      :: TFrc,Frc
  REAL(DOUBLE)        :: TFrcChk
 
  CHARACTER(LEN=6),PARAMETER :: Prog='TForce'
!------------------------------------------------------------------------------------- 
! Start up macro
 
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
!----------------------------------------------
! Get basis set and geometry
 
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
#ifdef PARALLEL
  CALL New(P,OnAll_O=.TRUE.)
#endif
  CALL Get(P,TrixFile('D',Args,1),BCast_O=.TRUE.)
  CALL New(TFrc,3*NAtoms)
#ifdef PERIODIC
! Get the Outer Cell Set
  CALL Get_CellSet(CS_OUT,'CS_OUT'//CurBase//CurGeom)
  CALL PPrint(CS_OUT,'outer sum',Prog)
#endif 
!--------------------------------------------------------------------------------
! TForce=2*Tr{P.dT} (Extra 2 to account for symmetry of T in the trace)
!--------------------------------------------------------------------------------
  LatFrc_T = Zero
 
  TFrc%D=Zero
#ifdef PARALLEL
  DO AtA=Beg%I(MyId),End%I(MyId)
#else
  DO AtA=1,NAtoms
#endif
     A1=3*(AtA-1)+1
     A2=3*AtA
     MA=BSiz%I(AtA)
     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1
        AtB=P%ColPt%I(JP)
        IF(SetAtomPair(GM,BS,AtA,AtB,Pair))THEN
           Q=P%BlkPt%I(JP)
           NB=BSiz%I(AtB)
           MN1=MA*NB-1
#ifdef PERIODIC
           B=Pair%B
           DO NC=1,CS_OUT%NCells
              Pair%B=B+CS_OUT%CellCarts%D(:,NC)
              Pair%AB2=(Pair%A(1)-Pair%B(1))**2  &
                      +(Pair%A(2)-Pair%B(2))**2  &
                      +(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair)) THEN
                 F_nlm(1:3)    = Four*TrPdT(BS,Pair,P%MTrix%D(Q:Q+MN1))
                 TFrc%D(A1:A2) = TFrc%D(A1:A2) + F_nlm(1:3)
 
                 nlm = AtomToFrac(GM,CS_OUT%CellCarts%D(1:3,NC))
                 LatFrc_T(1,1:3) = LatFrc_T(1,1:3) + nlm(1)*F_nlm(1:3)
                 LatFrc_T(2,1:3) = LatFrc_T(2,1:3) + nlm(2)*F_nlm(1:3)
                 LatFrc_T(3,1:3) = LatFrc_T(3,1:3) + nlm(3)*F_nlm(1:3)
 
              ENDIF
              
           ENDDO
#else
           TFrc%D(A1:A2)=TFrc%D(A1:A2)  &
                        +Four*TrPdT(BS,Pair,P%MTrix%D(Q:Q+MN1))
#endif
        ENDIF
     ENDDO
  ENDDO
!--------------------------------------------------------------------------------
#ifdef PARALLEL
  TotFrcComp = 3*NAtoms
  CALL New(TotTFrc,TotFrcComp)
  CALL MPI_Reduce(TFrc%D(1),TotTFrc%D(1),TotFrcComp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,IErr)
  IF(MyID == 0) THEN
    TFrc%D(1:TotFrcComp) = TotTFrc%D(1:TotFrcComp)
  ENDIF
  CALL Delete(TotTFrc)
#endif

! Do some checksumming, resumming and IO 
! CALL PPrint(TFrc,'dT/dR')
  CALL PChkSum(TFrc,'dT/dR',Proc_O=Prog)  
! Sum in contribution to total force
  CALL New(Frc,3*NAtoms)
  CALL Get(Frc,'GradE',Tag_O=CurGeom)
  Frc%D=Frc%D+TFrc%D
  CALL Put(Frc,'GradE',Tag_O=CurGeom)
!------------------------------------------------------------------------------
! Tidy up 
  CALL Delete(P)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(TFrc)
  CALL ShutDown(Prog)
END PROGRAM TForce
