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
  TYPE(BCSR)                  :: P
  TYPE(DBL_VECT)              :: TotTFrc
  INTEGER                     :: IErr,TotFrcComp
#else
  TYPE(BCSR)                  :: P
#endif
  INTEGER                     :: NC,I,J
  REAL(DOUBLE),DIMENSION(3)   :: A,B,F_nlm,nlm
  TYPE(DBL_RNK2)              :: LatFrc_T
  TYPE(AtomPair)              :: Pair
  TYPE(BSET)                  :: BS
  TYPE(CRDS)                  :: GM
  TYPE(ARGMT)                 :: Args
  INTEGER                     :: Q,R,AtA,AtB,JP,MB,MA,NB,MN1,A1,A2
  TYPE(HGRho)                 :: Rho
  TYPE(DBL_VECT)              :: TFrc,Frc
  REAL(DOUBLE)                :: TFrcChk
 
  CHARACTER(LEN=6),PARAMETER  :: Prog='TForce'
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
  TFrc%D   = Zero
  CALL New(LatFrc_T,(/3,3/))
  LatFrc_T%D = Zero
!--------------------------------------------------------------------------------
! TForce=2*Tr{P.dT} (Extra 2 to account for symmetry of T in the trace)
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
           A=Pair%A
           B=Pair%B
           DO NC=1,CS_OUT%NCells
              Pair%A=A
              Pair%B=B+CS_OUT%CellCarts%D(:,NC)
              Pair%AB2=(Pair%A(1)-Pair%B(1))**2  &
                      +(Pair%A(2)-Pair%B(2))**2  &
                      +(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair)) THEN
                 F_nlm(1:3)    = TrPdT(BS,Pair,P%MTrix%D(Q:Q+MN1))
                 TFrc%D(A1:A2) = TFrc%D(A1:A2) + Two*F_nlm(1:3)
!                Lattice Forces
                 nlm        = AtomToFrac(GM,Pair%A)
                 LatFrc_T%D = LatFrc_T%D + Two*LaticeForce(GM,nlm,F_nlm)
              ENDIF
              Pair%A=A+CS_OUT%CellCarts%D(:,NC)
              Pair%B=B
              Pair%AB2=(Pair%A(1)-Pair%B(1))**2  &
                      +(Pair%A(2)-Pair%B(2))**2  &
                      +(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair)) THEN
                 F_nlm(1:3)    = TrPdT(BS,Pair,P%MTrix%D(Q:Q+MN1))
                 TFrc%D(A1:A2) = TFrc%D(A1:A2) + Two*F_nlm(1:3)
!                Lattice Forces
                 nlm        = AtomToFrac(GM,Pair%A)
                 LatFrc_T%D = LatFrc_T%D + Two*LaticeForce(GM,nlm,F_nlm)
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDDO
!
!
!
  WRITE(*,*) 'T'
  DO AtA=1,NAtoms
     A1=3*(AtA-1)+1
     A2=3*AtA
     WRITE(*,'(I3,2x,3(D23.16,2X))') AtA,TFrc%D(A1:A2)
  ENDDO
  WRITE(*,*) 
  DO I=1,3
     WRITE(*,*) (LatFrc_T%D(I,J),J=1,3) 
  ENDDO
!--------------------------------------------------------------------------------
#ifdef PARALLEL
  TotFrcComp = 3*NAtoms
  CALL New(TotTFrc,TotFrcComp)
  CALL MPI_Reduce(TFrc%D(1),TotTFrc%D(1),TotFrcComp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MONDO_COMM,IErr)
  IF(MyID == 0) THEN
    TFrc%D(1:TotFrcComp) = TotTFrc%D(1:TotFrcComp)
  ENDIF
  CALL Delete(TotTFrc)
#endif
! Do some checksumming, resumming and IO 
  CALL PChkSum(TFrc,'dT/dR',Proc_O=Prog)  
! Sum in contribution to total force
  DO AtA=1,NAtoms
     A1=3*(AtA-1)+1
     A2=3*AtA
     GM%Gradients%D(1:3,AtA) =  GM%Gradients%D(1:3,AtA)+TFrc%D(A1:A2)
  ENDDO
  GM%PBC%LatFrc%D = GM%PBC%LatFrc%D+LatFrc_T%D
  CALL Put(GM,Tag_O=CurGeom)
!------------------------------------------------------------------------------
! Tidy up 
  CALL Delete(TFrc)
  CALL Delete(LatFrc_T)
  CALL Delete(P)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL ShutDown(Prog)
END PROGRAM TForce
