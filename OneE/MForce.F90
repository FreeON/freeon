PROGRAM MForce
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
  USE BlokTrPdM

#ifdef PARALLEL
  USE MondoMPI
  TYPE(BCSR)                  :: P
  TYPE(DBL_RNK2)              :: MFrcTmp
  INTEGER                     :: IErr
#else
  TYPE(BCSR)                  :: P
#endif
  INTEGER                     :: NC,I,J
  REAL(DOUBLE),DIMENSION(3)   :: A,B,nlm
  TYPE(AtomPair)              :: Pair
  TYPE(BSET)                  :: BS
  TYPE(CRDS)                  :: GM
  TYPE(ARGMT)                 :: Args
  INTEGER                     :: Q,R,AtA,AtB,JP,MB,MA,NB,MN1,A1,B1
  TYPE(HGRho)                 :: Rho
  TYPE(DBL_VECT)              :: COrig
  TYPE(DBL_RNK2)              :: MFrc
  REAL(DOUBLE)                :: F_nlm(18)
  CHARACTER(LEN=*),PARAMETER  :: Prog='MForce'
!------------------------------------------------------------------------------------- 
! Start up macro
!
  CALL StartUp(Args,Prog,Serial_O=.FALSE.)
!-------------------------------------------------------
! Get basis set and geometry
!
  CALL Get(BS,Tag_O=CurBase)
  CALL Get(GM,Tag_O=CurGeom)
!-------------------------------------------------------
! Get Multipole origine. TODO
!
  CALL New(COrig,3)
  COrig%D=0D0
!-------------------------------------------------------
! Get DM
!
#ifdef PARALLEL
  CALL New(P,OnAll_O=.TRUE.)
#endif
  CALL Get(P,TrixFile('D',Args,1),BCast_O=.TRUE.)
  CALL New(MFrc,(/3,3*NAtoms/))
  CALL DBL_VECT_EQ_DBL_SCLR(9*NAtoms,MFrc%D(1,1),0.0d0)
!-------------------------------------------------------
! Let's go.
!
#ifdef PARALLEL
  DO AtA=Beg%I(MyId),End%I(MyId)
#else
  DO AtA=1,NAtoms
#endif
     A1=3*(AtA-1)+1
     MA=BSiz%I(AtA)
     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1
        AtB=P%ColPt%I(JP)
        B1=3*(AtB-1)+1
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
                   &  +(Pair%A(2)-Pair%B(2))**2  &
                   &  +(Pair%A(3)-Pair%B(3))**2
              IF(TestAtomPair(Pair)) THEN
                 F_nlm(1:18)=TrPdM(BS,Pair,P%MTrix%D(Q:Q+MN1),COrig)
                 !Atom A
                 CALL DAXPY(9,1D0,F_nlm( 1),1,MFrc%D(1,A1),1)
                 !Atom B
                 CALL DAXPY(9,1D0,F_nlm(10),1,MFrc%D(1,B1),1)
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDDO
!-------------------------------------------------------
! Reduce to ROOT if PARALLEL and copy to GM.
!
#ifdef PARALLEL
  CALL New(MFrcTmp,(/3,3*NAtoms/))
  CALL MPI_Reduce(MFrc%D(1,1),MFrcTmp%D(1,1),9*NAtoms,MPI_DOUBLE_PRECISION, &
       &          MPI_SUM,0,MONDO_COMM,IErr)
  !CALL DCOPY(9*NAtoms,MFrcTmp%D(1,1),1,GM%DipoleFrc%D(1,1),1)
  CALL DCOPY(9*NAtoms,MFrcTmp%D(1,1),1,MFrc%D(1,1),1)
  CALL Delete(MFrcTmp)
#else
  !CALL DCOPY(9*NAtoms,   MFrc%D(1,1),1,GM%DipoleFrc%D(1,1),1)
#endif
!-------------------------------------------------------
! Put MForce to disk
!
  CALL Put(GM,Tag_O=CurGeom)
!-----------------------------------------------------------
! Printing
!
  CALL Print_DBL_RNK2(MFrc,'dMu/dA',Unit_O=6)
  !CALL PChkSum(MFrc,'MForce',Prog)
  CALL PPrint( MFrc,'MForce')
!-------------------------------------------------------
! Tidy up 
!
  CALL Delete(MFrc)
  CALL Delete(P)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL Delete(COrig)
  CALL ShutDown(Prog)
END PROGRAM MForce
