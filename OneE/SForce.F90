!    COMPUTE THE FORCE DUE TO CHANGES IN THE DENSITY MATRIX:
!    dP.(2T+J+K)=-2*dS.P.F.P (Early work by McWeeny, need a cite...)
!    Authors: Matt Challacombe and CJ Tymczak
!------------------------------------------------------------------------------
PROGRAM SForce
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
  USE BraBloks
  USE BlokTrWdS

#ifdef PARALLEL
  USE MondoMPI
!  TYPE(BCSR)                 :: T1,F,P
  TYPE(DBCSR)                :: T1,F,P_DBCSR
  TYPE(BCSR)                 :: P
  TYPE(DBL_VECT)             :: TotSFrc
  INTEGER                    :: IErr,TotFrcComp
#else
  TYPE(BCSR)                 :: T1,F,P
#endif

  TYPE(AtomPair)             :: Pair
  TYPE(BSET)                 :: BS
  TYPE(CRDS)                 :: GM
  TYPE(ARGMT)                :: Args
  INTEGER                    :: Q,R,AtA,AtB,NN,JP,MB,MA,NB,MN1,A1,A2
  TYPE(HGRho)                :: Rho
  TYPE(DBL_VECT)             :: SFrc
  REAL(DOUBLE)               :: SFrcChk

  INTEGER                    :: NC,I,J
  REAL(DOUBLE),DIMENSION(3)  :: B,F_nlm,nlm
  TYPE(DBL_RNK2)             :: LatFrc_S
  CHARACTER(LEN=6),PARAMETER :: Prog='SForce'
!---------------------------------------------- 
! Start up macro
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
!----------------------------------------------
! Get basis set and geometry
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
!---------------------------------------------- 
! Allocations 

#ifdef PARALLEL
#else
  CALL New(P,OnAll_O=.TRUE.)
  CALL New(F)
  CALL New(T1)
#endif
!--------------------------------------
! Compute W=P.F.P
#ifdef PARALLEL
  CALL Get(P_DBCSR,TrixFile('D',Args,1))
#else
  CALL Get(P,TrixFile('D',Args,1))
#endif
! Is this a bug if we don't use the extrapolated Fockian?
  CALL Get(F,TrixFile('F',Args,0))
#ifdef PARALLEL
  CALL Multiply(P_DBCSR,F,T1)       ! T1:=P.F
  CALL Multiply(T1,P_DBCSR,F)       ! F:=P.F.P
  CALL Filter(P_DBCSR,F)            ! P=Filter[P.F.P]      
  CALL SetEq(P,P_DBCSR)
  CALL BcastBCSR(P)
  CALL Delete(P_DBCSR)

#else
  CALL Multiply(P,F,T1)       ! T1:=P.F
  CALL Multiply(T1,P,F)       ! F:=P.F.P
  CALL Filter(P,F)            ! P=Filter[P.F.P]      
#endif
  CALL Delete(F)
  CALL Delete(T1)
  CALL New(SFrc,3*NAtoms)
  SFrc%D   = Zero
  CALL NewBraBlok(BS,Gradients_O=.TRUE.)
#ifdef PERIODIC
  CALL New(LatFrc_S,(/3,3/))
  LatFrc_S%D = Zero
#endif
!--------------------------------------------------------------------------------
! SForce=-2*Tr{P.F.P.dS} (Extra 2 to account for symmetry of S in the trace)
!--------------------------------------------------------------------------------
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
           B=Pair%B
           DO NC=1,CS_OUT%NCells
              Pair%B=B+CS_OUT%CellCarts%D(:,NC)
              Pair%AB2=(Pair%A(1)-Pair%B(1))**2  &
                      +(Pair%A(2)-Pair%B(2))**2  &
                      +(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair)) THEN
                 F_nlm(1:3)    = TrWdS(BS,Pair,P%MTrix%D(Q:Q+MN1))
                 IF(.NOT. Pair%SameAtom) THEN
                    SFrc%D(A1:A2) = SFrc%D(A1:A2) - Four*F_nlm(1:3)
                 ENDIF
!                Lattice Forces
                 nlm        = AtomToFrac(GM,Pair%A-Pair%B)
                 LatFrc_S%D = LatFrc_S%D - Two*LaticeForce(GM,nlm,F_nlm)
              ENDIF
           ENDDO

        ENDIF
     ENDDO
  ENDDO
!--------------------------------------------------------------------------------
#ifdef PARALLEL
  TotFrcComp = 3*NAtoms
  CALL New(TotSFrc,TotFrcComp)
  CALL MPI_Reduce(SFrc%D(1),TotSFrc%D(1),TotFrcComp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MONDO_COMM,IErr)
  IF(MyID == 0) THEN
    SFrc%D(1:TotFrcComp) = TotSFrc%D(1:TotFrcComp)
  ENDIF
  CALL Delete(TotSFrc)
#endif
! Do some checksumming and IO 
!  CALL PPrint(SFrc,'dS/dR')
  CALL PChkSum(SFrc,'dS/dR',Proc_O=Prog)  
! Start this off as the first contrib to total gradient 
  CALL Put(SFrc,'GradE',Tag_O=CurGeom)
!
  CALL Put(LatFrc_S,'LatFrc_S',Tag_O=CurGeom)
  CALL Delete(LatFrc_S)
!--------------------------------------------------------------------------------
! Tidy up 
!--------------------------------------------------------------------------------
  CALL Delete(P)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(SFrc)
  CALL DeleteBraBlok(Gradients_O=.TRUE.)
  CALL ShutDown(Prog)
!
END PROGRAM SForce
