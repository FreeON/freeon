!    COMPUTE THE FORCE CORESPONDING TO THE DERIVATIVE OF THE 
!    KINETIC ENERGY MATRIX, TForce=2*Tr{P.dT}
!    Authors: Matt Challacombe 
!----------------------------------------------------------------------------------
PROGRAM UForce
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
  USE ECPBlock
#ifdef PARALLEL
  USE MondoMPI
  TYPE(BCSR)                  :: P
  TYPE(DBL_VECT)              :: TotUFrc
  INTEGER                     :: IErr,TotFrcComp
#else
  TYPE(BCSR)                  :: P
#endif
  INTEGER                     :: NC,I,J
  REAL(DOUBLE),DIMENSION(3)   :: B,F
  TYPE(AtomPair)              :: Pair
  TYPE(BSET)                  :: BS
  TYPE(CRDS)                  :: GM
  TYPE(ARGMT)                 :: Args
  INTEGER                     :: Q,R,AtA,AtB,AtC,KC,JP,MB,MA,NB,MN1,A1,A2,C1,C2
  TYPE(HGRho)                 :: Rho
  TYPE(DBL_VECT)              :: UFrc,Frc
  REAL(DOUBLE)                :: Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,UFrcChk,F1a,F2a,F3a,f11,f10
 
  CHARACTER(LEN=6),PARAMETER  :: Prog='UForce'
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
  CALL New(UFrc,3*NAtoms)
  UFrc%D=Zero
  ! First compute 
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

                 ! Go over ECP centers
                 DO AtC=1,NAtoms
                    C1=3*(AtC-1)+1
                    C2=3*AtC
                    KC=GM%AtTyp%I(AtC)
                    IF(BS%NTyp1PF%I(KC)>0)THEN
                       Cx=GM%Carts%D(1,AtC) 
                       Cy=GM%Carts%D(2,AtC) 
                       Cz=GM%Carts%D(3,AtC) 
                       F=Two*TrPdU(BS,Pair,KC,Cx,Cy,Cz,P%MTrix%D(Q:Q+MN1))
                       ! Derivative contribution from bra <A| and ket |B>
                       UFrc%D(A1:A2)=UFrc%D(A1:A2)+F
                       ! then add in the ECP center |C| part
                       UFrc%D(C1:C2)=UFrc%D(C1:C2)-F
                    ENDIF
                 ENDDO
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDDO
!!$!=====================================================================================
!!$     f10=Zero
!!$     DO AtA=1,NAtoms
!!$     A1=3*(AtA-1)+1
!!$     A2=3*AtA
!!$     MA=BSiz%I(AtA)
!!$     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1
!!$        AtB=P%ColPt%I(JP)
!!$           IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
!!$              Q=P%BlkPt%I(JP)
!!$              B=Pair%B
!!$          NB=BSiz%I(AtB)
!!$           MN1=MA*NB-1
!!$              Pair%AB2=(Pair%A(1)-Pair%B(1))**2 &
!!$                   +(Pair%A(2)-Pair%B(2))**2 &
!!$                   +(Pair%A(3)-Pair%B(3))**2
!!$              IF(TestAtomPair(Pair)) THEN
!!$                 ! Go over ECP centers
!!$                 DO AtC=1,NAtoms
!!$                    KC=GM%AtTyp%I(AtC)
!!$                    IF(BS%NTyp1PF%I(KC)>0)THEN
!!$                       Cx=GM%Carts%D(1,AtC)
!!$                       Cy=GM%Carts%D(2,AtC)
!!$                       Cz=GM%Carts%D(3,AtC)
!!$                       Ax = Pair%A(1)
!!$                       Ay = Pair%A(2)
!!$                       Az = Pair%A(3)
!!$                       Bx = Pair%B(1)
!!$                       By = Pair%B(2)
!!$                       Bz = Pair%B(3)
!!$                       f10=f10+DUBlock(BS,Pair,Ax,Ay,Az,Bx,By,Bz,KC,Cx,Cy,Cz,P%MTrix%D(Q:Q+MN1))
!!$                    ENDIF
!!$                 ENDDO
!!$              ENDIF
!!$           ENDIF
!!$        ENDDO
!!$     ENDDO
!!$!=====================================================================================
!!$     GM%Carts%D(1,1)=GM%Carts%D(1,1)+1D-5
!!$     f11=Zero
!!$     DO AtA=1,NAtoms
!!$     A1=3*(AtA-1)+1
!!$     A2=3*AtA
!!$     MA=BSiz%I(AtA)
!!$     DO JP=P%RowPt%I(AtA),P%RowPt%I(AtA+1)-1
!!$        AtB=P%ColPt%I(JP)
!!$           IF(SetAtomPair(GM,BS,AtA,AtB,Pair)) THEN
!!$           Q=P%BlkPt%I(JP)
!!$          NB=BSiz%I(AtB)
!!$           MN1=MA*NB-1
!!$              B=Pair%B
!!$              Pair%AB2=(Pair%A(1)-Pair%B(1))**2 &
!!$                   +(Pair%A(2)-Pair%B(2))**2 &
!!$                   +(Pair%A(3)-Pair%B(3))**2
!!$              IF(TestAtomPair(Pair)) THEN
!!$                 ! Go over ECP centers
!!$                 DO AtC=1,NAtoms
!!$                    KC=GM%AtTyp%I(AtC)
!!$                    IF(BS%NTyp1PF%I(KC)>0)THEN
!!$                       Cx=GM%Carts%D(1,AtC)
!!$                       Cy=GM%Carts%D(2,AtC)
!!$                       Cz=GM%Carts%D(3,AtC)
!!$                       Ax = Pair%A(1)
!!$                       Ay = Pair%A(2)
!!$                       Az = Pair%A(3)
!!$                       Bx = Pair%B(1)
!!$                       By = Pair%B(2)
!!$                       Bz = Pair%B(3)
!!$                       f11=f11+DUBlock(BS,Pair,Ax,Ay,Az,Bx,By,Bz,KC,Cx,Cy,Cz,P%MTrix%D(Q:Q+MN1))
!!$                    ENDIF
!!$                 ENDDO
!!$              ENDIF
!!$           ENDIF
!!$        ENDDO
!!$     ENDDO
UFrc%D=Two*UFrc%D
!--------------------------------------------------------------------------------
#ifdef PARALLEL
  TotFrcComp = 3*NAtoms
  CALL New(TotUFrc,TotFrcComp)
  CALL MPI_Reduce(UFrc%D(1),TotUFrc%D(1),TotFrcComp,MPI_DOUBLE_PRECISION,MPI_SUM,0,MONDO_COMM,IErr)
  IF(MyID == 0) THEN
    UFrc%D(1:TotFrcComp)=TotUFrc%D(1:TotFrcComp)
  ENDIF
  CALL Delete(TotUFrc)
#endif
!  CALL PPrint(UFrc,'dU',Unit_O=6)
! Do some checksumming, resumming and IO 
  CALL PChkSum(UFrc,'dU/dR',Proc_O=Prog)  
! Sum in contribution to total force
  DO AtA=1,NAtoms
     A1=3*(AtA-1)+1
     A2=3*AtA
     GM%Gradients%D(1:3,AtA)=GM%Gradients%D(1:3,AtA)+UFrc%D(A1:A2)
  ENDDO
  CALL Put(GM,Tag_O=CurGeom)
!------------------------------------------------------------------------------
! Tidy up 
  CALL Delete(UFrc)
  CALL Delete(P)
  CALL Delete(BS)
  CALL Delete(GM)
  CALL ShutDown(Prog)
END PROGRAM UForce
